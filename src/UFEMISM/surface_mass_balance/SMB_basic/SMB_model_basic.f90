module SMB_model_basic

  use precisions, only: dp
  use parameters, only: pi, NaN
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use models_basic, only: atype_model, atype_model_context_allocate, &
    atype_model_context_initialise, atype_model_context_run, atype_model_context_remap
  use SMB_model_common, only: type_SMB_model_common
  use mpi_f08, only: MPI_WIN
  use ice_model_types, only: type_ice_model
  use climate_model_types, only: type_climate_model
  use grid_types, only: type_grid

  implicit none

  private

  public :: atype_SMB_model, type_SMB_model_context_allocate, &
    type_SMB_model_context_initialise, type_SMB_model_context_run, &
    type_SMB_model_context_remap

  type, abstract, extends(type_SMB_model_common) :: atype_SMB_model

    real(dp) :: t_next   !< Time when the SMB model should be run next

    contains

      ! These routines all consist of two parts: a 'common' part that is executed for
      ! all models inheriting from atype_SMB_model, and a 'specific' part that is
      ! only executed for each specific model class. The specific parts are defined
      ! in the deferred procedures 'allocate_SMB_model', 'initialise_SMB_model', etc.

      procedure, public :: allocate_model   => allocate_model_abs
      procedure, public :: deallocate_model => deallocate_model
      procedure, public :: initialise_model => initialise_model_abs
      procedure, public :: run_model        => run_model_abs
      procedure, public :: remap_model      => remap_model_abs

      procedure(allocate_SMB_model_ifc),   deferred :: allocate_SMB_model
      procedure(deallocate_SMB_model_ifc), deferred :: deallocate_SMB_model
      procedure(initialise_SMB_model_ifc), deferred :: initialise_SMB_model
      procedure(run_SMB_model_ifc),        deferred :: run_SMB_model
      procedure(remap_SMB_model_ifc),      deferred :: remap_SMB_model

      ! Factory functions to create model context objects

      procedure, nopass, public :: ct_allocate
      procedure, nopass, public :: ct_initialise
      procedure, nopass, public :: ct_run
      procedure, nopass, public :: ct_remap

  end type atype_SMB_model

  ! Context classes for allocate/initialise/run/remap
  ! =================================================

  type, extends(atype_model_context_allocate) :: type_SMB_model_context_allocate
  end type type_SMB_model_context_allocate

  type, extends(atype_model_context_initialise) :: type_SMB_model_context_initialise
    type(type_ice_model), pointer :: ice
  end type type_SMB_model_context_initialise

  type, extends(atype_model_context_run) :: type_SMB_model_context_run
    type(type_ice_model),     pointer :: ice
    type(type_climate_model), pointer :: climate
    type(type_grid),          pointer :: grid_smooth
  end type type_SMB_model_context_run

  type, extends(atype_model_context_remap) :: type_SMB_model_context_remap
    real(dp) :: time
  end type type_SMB_model_context_remap

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine allocate_SMB_model_ifc( self, context)
      import atype_SMB_model, type_SMB_model_context_allocate
      class(atype_SMB_model),                        intent(inout) :: self
      type(type_SMB_model_context_allocate), target, intent(in   ) :: context
    end subroutine allocate_SMB_model_ifc

    subroutine deallocate_SMB_model_ifc( self)
      import atype_SMB_model
      class(atype_SMB_model), intent(inout) :: self
    end subroutine deallocate_SMB_model_ifc

    subroutine initialise_SMB_model_ifc( self, context)
      import atype_SMB_model, type_SMB_model_context_initialise
      class(atype_SMB_model),                          intent(inout) :: self
      type(type_SMB_model_context_initialise), target, intent(in   ) :: context
    end subroutine initialise_SMB_model_ifc

    subroutine run_SMB_model_ifc( self, context)
      import atype_SMB_model, type_SMB_model_context_run
      class(atype_SMB_model),                   intent(inout) :: self
      type(type_SMB_model_context_run), target, intent(in   ) :: context
    end subroutine run_SMB_model_ifc

    subroutine remap_SMB_model_ifc( self, context)
      import atype_SMB_model, type_SMB_model_context_remap
      class(atype_SMB_model),                     intent(inout) :: self
      type(type_SMB_model_context_remap), target, intent(in   ) :: context
    end subroutine remap_SMB_model_ifc

  end interface

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  interface

    module subroutine allocate_model_abs( self, context)
      class(atype_SMB_model),                      intent(inout) :: self
      class(atype_model_context_allocate), target, intent(in   ) :: context
    end subroutine allocate_model_abs

    module subroutine deallocate_model( self)
      class(atype_SMB_model), intent(inout) :: self
    end subroutine deallocate_model

    module subroutine initialise_model_abs( self, context)
      class(atype_SMB_model),                        intent(inout) :: self
      class(atype_model_context_initialise), target, intent(in   ) :: context
    end subroutine initialise_model_abs

    module subroutine run_model_abs( self, context)
      class(atype_SMB_model),                 intent(inout) :: self
      class(atype_model_context_run), target, intent(in   ) :: context
    end subroutine run_model_abs

    module subroutine remap_model_abs( self, context)
      class(atype_SMB_model),                   intent(inout) :: self
      class(atype_model_context_remap), target, intent(in   ) :: context
    end subroutine remap_model_abs

    module function ct_allocate( name, region_name, mesh) result( context)
      character(len=*),           intent(in) :: name
      character(len=*),           intent(in) :: region_name
      type(type_mesh), target,    intent(in) :: mesh
      type(type_SMB_model_context_allocate) :: context
    end function ct_allocate

    module function ct_initialise( ice) result( context)
      type(type_ice_model), target, intent(in) :: ice
      type(type_SMB_model_context_initialise)  :: context
    end function ct_initialise

    module function ct_run( time, ice, climate, grid_smooth) result( context)
      real(dp),                          intent(in) :: time
      type(type_ice_model),     pointer, intent(in) :: ice
      type(type_climate_model), pointer, intent(in) :: climate
      type(type_grid),          pointer, intent(in) :: grid_smooth
      type(type_SMB_model_context_run)              :: context
    end function ct_run

    module function ct_remap( mesh_new, time) result( context)
      type(type_mesh), target, intent(in) :: mesh_new
      real(dp),                intent(in) :: time
      type(type_SMB_model_context_remap) :: context
    end function ct_remap

  end interface

end module SMB_model_basic
