module calving_model_basic

  use precisions, only: dp
  use parameters, only: pi, NaN
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use models_basic, only: atype_model, &
    atype_model_context_initialise, atype_model_context_run, atype_model_context_remap
  use calving_model_data, only: atype_calving_model_data
  use mpi_f08, only: MPI_WIN
  use ice_model_types, only: type_ice_model
  use climate_model_types, only: type_climate_model
  use grid_types, only: type_grid
  use reference_geometry_types, only: type_reference_geometry
  use ice_geometry_calculations, only: type_ice_geometry

  implicit none

  private

  public :: atype_calving_model, &
   type_calving_model_context_initialise , type_calving_model_context_run , &
    type_calving_model_context_remap

  type, abstract, extends(atype_calving_model_data) :: atype_calving_model
    !< Stuff that is common to all calving models
    !<
    !< (except for the variables that we want other models to
    !< be able to access, which are already defined in atype_calving_model_data)

    real(dp) :: t_next   !< Time when the calving model should be run next

    contains

      ! These routines all consist of two parts: a 'common' part that is executed for
      ! all models inheriting from atype_calving_model, and a 'specific' part that is
      ! only executed for each specific model class. The specific parts are defined
      ! in the deferred procedures 'allocate_calving_model', 'initialise_calving_model', etc.

      procedure, public :: allocate_calving_model
      procedure, public :: deallocate_model => deallocate_model
      procedure, public :: initialise_model => initialise_model_abs
      procedure, public :: run_model        => run_model_abs
      procedure, public :: remap_model      => remap_model_abs

      procedure(calving_model_allocate_ifc),   deferred :: allocate
      procedure(deallocate_calving_model_ifc), deferred :: deallocate_calving_model
      procedure(initialise_calving_model_ifc), deferred :: initialise_calving_model
      procedure(run_calving_model_ifc),        deferred :: run_calving_model
      procedure(remap_calving_model_ifc),      deferred :: remap_calving_model

      ! Factory functions to create model context objects
      procedure, nopass, public :: ct_initialise
      procedure, nopass, public :: ct_run
      procedure, nopass, public :: ct_remap

  end type atype_calving_model

  ! Context classes for allocate/initialise/run/remap
  ! =================================================

  type, extends(atype_model_context_initialise) :: type_calving_model_context_initialise
  end type type_calving_model_context_initialise

  type, extends(atype_model_context_run) :: type_calving_model_context_run
    type(type_ice_model),        pointer :: ice
    type(type_ice_geometry),     pointer :: icegeom
  end type type_calving_model_context_run

  type, extends(atype_model_context_remap) :: type_calving_model_context_remap
  end type type_calving_model_context_remap

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine calving_model_allocate_ifc( self, name, region_name, mesh)
      import atype_calving_model, type_mesh
      class(atype_calving_model), intent(inout) :: self
      character(len=*),           intent(in   ) :: name
      character(len=*),           intent(in   ) :: region_name
      type(type_mesh), target,    intent(in   ) :: mesh
    end subroutine calving_model_allocate_ifc

    subroutine deallocate_calving_model_ifc( self)
      import atype_calving_model
      class(atype_calving_model), intent(inout) :: self
    end subroutine deallocate_calving_model_ifc

    subroutine initialise_calving_model_ifc( self, context)
      import atype_calving_model, type_calving_model_context_initialise
      class(atype_calving_model),                          intent(inout) :: self
      type(type_calving_model_context_initialise), target, intent(in   ) :: context
    end subroutine initialise_calving_model_ifc

    subroutine run_calving_model_ifc( self, context)
      import atype_calving_model, type_calving_model_context_run
      class(atype_calving_model),                   intent(inout) :: self
      type(type_calving_model_context_run), target, intent(in   ) :: context
    end subroutine run_calving_model_ifc

    subroutine remap_calving_model_ifc( self, context)
      import atype_calving_model, type_calving_model_context_remap
      class(atype_calving_model),                     intent(inout) :: self
      type(type_calving_model_context_remap), target, intent(in   ) :: context
    end subroutine remap_calving_model_ifc

  end interface

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  interface

    module subroutine deallocate_model( self)
      class(atype_calving_model), intent(inout) :: self
    end subroutine deallocate_model

    module subroutine initialise_model_abs( self, context)
      class(atype_calving_model),                    intent(inout) :: self
      class(atype_model_context_initialise), target, intent(in   ) :: context
    end subroutine initialise_model_abs

    module subroutine run_model_abs( self, context)
      class(atype_calving_model),             intent(inout) :: self
      class(atype_model_context_run), target, intent(in   ) :: context
    end subroutine run_model_abs

    module subroutine remap_model_abs( self, context)
      class(atype_calving_model),               intent(inout) :: self
      class(atype_model_context_remap), target, intent(in   ) :: context
    end subroutine remap_model_abs

    module function ct_initialise( ) result( context)
      type(type_calving_model_context_initialise)       :: context
    end function ct_initialise

    module function ct_run( time, ice, icegeom) result( context)
      real(dp),                         intent(in) :: time
      type(type_ice_model),     target, intent(in) :: ice
      type(type_ice_geometry),  target, intent(in) :: icegeom
      type(type_calving_model_context_run)         :: context
    end function ct_run

    module function ct_remap( mesh_new) result( context)
      type(type_mesh),               target, intent(in) :: mesh_new
      type(type_calving_model_context_remap)            :: context
    end function ct_remap

  end interface

contains

  subroutine allocate_calving_model( self, name, region_name, mesh, nz)
    !< Allocate stuff that is common to all calving models
    !< (call this from your calving model-specific allocation routine)

    ! In/output variables:
    class(atype_calving_model), intent(inout) :: self
    character(len=*),           intent(in   ) :: name
    character(len=*),           intent(in   ) :: region_name
    type(type_mesh), target,    intent(in   ) :: mesh
    integer,                    intent(in   ) :: nz

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_calving_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate stuff that is common to all models
    call self%allocate_model( name, region_name, mesh)

    ! Allocate stuff that is specific to calving models

    call self%create_field( self%Hi_calved, self%wHi_calved, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Hi_calved', &
      long_name = 'Ice thickness after calving', &
      units     = 'm', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_calving_model

end module calving_model_basic
