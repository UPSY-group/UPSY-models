module climate_model_basic

  use precisions, only: dp
  use parameters, only: pi, NaN
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use models_basic, only: atype_model, atype_model_context_allocate, &
    atype_model_context_initialise, atype_model_context_run, atype_model_context_remap
  use climate_model_common, only: type_climate_model_common
  use mpi_f08, only: MPI_WIN
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry

  implicit none

  private

  public :: atype_climate_model, type_climate_model_context_allocate, &
    type_climate_model_context_initialise, type_climate_model_context_run, &
    type_climate_model_context_remap

  type, abstract, extends(type_climate_model_common) :: atype_climate_model

    real(dp) :: t_next   !< Time when the climate model should be run next

    contains

      ! These routines all consist of two parts: a 'common' part that is executed for
      ! all models inheriting from atype_climate_model, and a 'specific' part that is
      ! only executed for each specific model class. The specific parts are defined
      ! in the deferred procedures 'allocate_climate_model', 'initialise_climate_model', etc.

      procedure, public :: allocate_model   => allocate_model_abs
      procedure, public :: deallocate_model => deallocate_model
      procedure, public :: initialise_model => initialise_model_abs
      procedure, public :: run_model        => run_model_abs
      procedure, public :: remap_model      => remap_model_abs

      procedure(allocate_climate_model_ifc),   deferred :: allocate_climate_model
      procedure(deallocate_climate_model_ifc), deferred :: deallocate_climate_model
      procedure(initialise_climate_model_ifc), deferred :: initialise_climate_model
      procedure(run_climate_model_ifc),        deferred :: run_climate_model
      procedure(remap_climate_model_ifc),      deferred :: remap_climate_model

      ! Factory functions to create model context objects

      procedure, nopass, public :: ct_allocate
      procedure, nopass, public :: ct_initialise
      procedure, nopass, public :: ct_run
      procedure, nopass, public :: ct_remap

  end type atype_climate_model

  ! Context classes for allocate/initialise/run/remap
  ! =================================================

  type, extends(atype_model_context_allocate) :: type_climate_model_context_allocate
  end type type_climate_model_context_allocate

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

    subroutine allocate_climate_model_ifc( self, context)
      import atype_climate_model, type_climate_model_context_allocate
      class(atype_climate_model),                        intent(inout) :: self
      type(type_climate_model_context_allocate), target, intent(in   ) :: context
    end subroutine allocate_climate_model_ifc

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

    module subroutine allocate_model_abs( self, context)
      class(atype_climate_model),                  intent(inout) :: self
      class(atype_model_context_allocate), target, intent(in   ) :: context
    end subroutine allocate_model_abs

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

    module function ct_allocate( name, region_name, mesh) result( context)
      character(len=*),              intent(in) :: name
      character(len=*),              intent(in) :: region_name
      type(type_mesh), target,       intent(in) :: mesh
      type(type_climate_model_context_allocate) :: context
    end function ct_allocate

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

end module climate_model_basic
