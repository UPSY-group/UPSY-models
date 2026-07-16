module climate_model_basic

  use precisions, only: dp
  use parameters, only: pi, NaN
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use models_basic, only: atype_model, &
    atype_model_context_initialise, atype_model_context_run, atype_model_context_remap
  use climate_model_data, only: atype_climate_model_data
  use mpi_f08, only: MPI_WIN
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry

  implicit none

  private

  public :: atype_climate_model, &
    type_climate_model_context_initialise, type_climate_model_context_run, &
    type_climate_model_context_remap

  type, abstract, extends(atype_climate_model_data) :: atype_climate_model
    !< Stuff that is common to all climate models
    !<
    !< (except for the variables that we want other models to
    !< be able to access, which are already defined in atype_climate_model_data)

      real(dp) :: t_next   !< Time when the climate model should be run next

    contains

      ! These routines all consist of two parts: a 'common' part that is executed for
      ! all models inheriting from atype_climate_model, and a 'specific' part that is
      ! only executed for each specific model class. The specific parts are defined
      ! in the deferred procedures 'allocate_climate_model', 'initialise_climate_model', etc.

      procedure, public :: allocate_climate_model
      procedure, public :: deallocate_model => deallocate_model
      procedure, public :: initialise_model => initialise_model_abs
      procedure, public :: run_model        => run_model_abs
      procedure, public :: remap_model      => remap_model_abs

      procedure(climate_model_allocate_ifc),   deferred :: allocate
      procedure(deallocate_climate_model_ifc), deferred :: deallocate_climate_model
      procedure(initialise_climate_model_ifc), deferred :: initialise_climate_model
      procedure(run_climate_model_ifc),        deferred :: run_climate_model
      procedure(remap_climate_model_ifc),      deferred :: remap_climate_model

      ! Factory functions to create model context objects
      procedure, nopass, public :: ct_initialise
      procedure, nopass, public :: ct_run
      procedure, nopass, public :: ct_remap

  end type atype_climate_model

  ! Context classes for allocate/initialise/run/remap
  ! =================================================

  type, extends(atype_model_context_initialise) :: type_climate_model_context_initialise
    type(type_reference_geometry), pointer :: refgeo_init
    type(type_reference_geometry), pointer :: refgeo_PD
  end type type_climate_model_context_initialise

  type, extends(atype_model_context_run) :: type_climate_model_context_run
    type(type_ice_model), pointer :: ice
  end type type_climate_model_context_run

  type, extends(atype_model_context_remap) :: type_climate_model_context_remap
  end type type_climate_model_context_remap

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine climate_model_allocate_ifc( self, name, region_name, mesh)
      import atype_climate_model, type_mesh
      class(atype_climate_model), intent(inout) :: self
      character(len=*),           intent(in   ) :: name
      character(len=*),           intent(in   ) :: region_name
      type(type_mesh), target,    intent(in   ) :: mesh
    end subroutine climate_model_allocate_ifc

    subroutine deallocate_climate_model_ifc( self)
      import atype_climate_model
      class(atype_climate_model), intent(inout) :: self
    end subroutine deallocate_climate_model_ifc

    subroutine initialise_climate_model_ifc( self, context)
      import atype_climate_model, type_climate_model_context_initialise
      class(atype_climate_model),                          intent(inout) :: self
      type(type_climate_model_context_initialise), target, intent(in   ) :: context
    end subroutine initialise_climate_model_ifc

    subroutine run_climate_model_ifc( self, context)
      import atype_climate_model, type_climate_model_context_run
      class(atype_climate_model),                   intent(inout) :: self
      type(type_climate_model_context_run), target, intent(in   ) :: context
    end subroutine run_climate_model_ifc

    subroutine remap_climate_model_ifc( self, context)
      import atype_climate_model, type_climate_model_context_remap
      class(atype_climate_model),                     intent(inout) :: self
      type(type_climate_model_context_remap), target, intent(in   ) :: context
    end subroutine remap_climate_model_ifc

  end interface

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  interface

    module subroutine deallocate_model( self)
      class(atype_climate_model), intent(inout) :: self
    end subroutine deallocate_model

    module subroutine initialise_model_abs( self, context)
      class(atype_climate_model),                    intent(inout) :: self
      class(atype_model_context_initialise), target, intent(in   ) :: context
    end subroutine initialise_model_abs

    module subroutine run_model_abs( self, context)
      class(atype_climate_model),             intent(inout) :: self
      class(atype_model_context_run), target, intent(in   ) :: context
    end subroutine run_model_abs

    module subroutine remap_model_abs( self, context)
      class(atype_climate_model),               intent(inout) :: self
      class(atype_model_context_remap), target, intent(in   ) :: context
    end subroutine remap_model_abs

    module function ct_initialise( refgeo_init, refgeo_PD) result( context)
      type(type_reference_geometry), target, intent(in) :: refgeo_init, refgeo_PD
      type(type_climate_model_context_initialise)       :: context
    end function ct_initialise

    module function ct_run( time, ice) result( context)
      real(dp),                         intent(in) :: time
      type(type_ice_model),     target, intent(in) :: ice
      type(type_climate_model_context_run)         :: context
    end function ct_run

    module function ct_remap( mesh_new) result( context)
      type(type_mesh), target,    intent(in) :: mesh_new
      type(type_climate_model_context_remap) :: context
    end function ct_remap

  end interface

contains

  subroutine allocate_climate_model( self, name, region_name, mesh)
    !< Allocate stuff that is common to all climate models
    !< (call this from your climate model-specific allocation routine)

    ! In/output variables:
    class(atype_climate_model), intent(inout) :: self
    character(len=*),           intent(in   ) :: name
    character(len=*),           intent(in   ) :: region_name
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_climate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all models
    call self%allocate_model( name, region_name, mesh)

    ! Allocate all the stuff that is specific to climate models

    call self%create_field( self%T2m, self%wT2m, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m', &
      long_name = 'Monthly mean 2-m air temperature', &
      units     = 'K', &
      remap_method = 'reallocate')

    call self%create_field( self%Precip, self%wPrecip, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Precip', &
      long_name = 'Monthly total precipitation', &
      units     = 'm.w.e.', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_climate_model

end module climate_model_basic
