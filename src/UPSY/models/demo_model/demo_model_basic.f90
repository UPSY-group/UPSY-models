module demo_model_basic

  use precisions, only: dp
  use parameters, only: pi, NaN
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use models_basic, only: atype_model, atype_model_context_allocate, &
    atype_model_context_initialise, atype_model_context_run, atype_model_context_remap
  use demo_model_state, only: type_demo_model_state
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_demo_model, type_demo_model_context_allocate, &
    type_demo_model_context_initialise, type_demo_model_context_run, &
    type_demo_model_context_remap

  type, abstract, extends(atype_model) :: atype_demo_model

      type(type_demo_model_state) :: s

    contains

      ! These routines all consist of two parts: a 'common' part that is executed for
      ! all models inheriting from atype_demo_model, and a 'specific' part that is
      ! only executed for each specific model class. The specific parts are defined
      ! in the deferred procedures 'allocate_demo_model', 'initialise_demo_model', etc.

      procedure, public :: allocate_model   => allocate_model_abs
      procedure, public :: deallocate_model => deallocate_model
      procedure, public :: initialise_model => initialise_model_abs
      procedure, public :: run_model        => run_model_abs
      procedure, public :: remap_model      => remap_model_abs

      procedure(allocate_demo_model_ifc),   deferred :: allocate_demo_model
      procedure(deallocate_demo_model_ifc), deferred :: deallocate_demo_model
      procedure(initialise_demo_model_ifc), deferred :: initialise_demo_model
      procedure(run_demo_model_ifc),        deferred :: run_demo_model
      procedure(remap_demo_model_ifc),      deferred :: remap_demo_model

      ! Factory functions to create model context objects

      procedure, nopass, public :: ct_allocate
      procedure, nopass, public :: ct_initialise
      procedure, nopass, public :: ct_run
      procedure, nopass, public :: ct_remap

  end type atype_demo_model

  ! Context classes for allocate/initialise/run/remap
  ! =================================================

  type, extends(atype_model_context_allocate) :: type_demo_model_context_allocate
    integer  :: nz
  end type type_demo_model_context_allocate

  type, extends(atype_model_context_initialise) :: type_demo_model_context_initialise
    real(dp) :: H0
    real(dp) :: till_friction_angle_uniform
    real(dp) :: beta_sq_uniform
  end type type_demo_model_context_initialise

  type, extends(atype_model_context_run) :: type_demo_model_context_run
    real(dp) :: H_new
    real(dp) :: dH
  end type type_demo_model_context_run

  type, extends(atype_model_context_remap) :: type_demo_model_context_remap
  end type type_demo_model_context_remap

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine allocate_demo_model_ifc( self, context)
      import atype_demo_model, type_demo_model_context_allocate
      class(atype_demo_model),                        intent(inout) :: self
      type(type_demo_model_context_allocate), target, intent(in   ) :: context
    end subroutine allocate_demo_model_ifc

    subroutine deallocate_demo_model_ifc( self)
      import atype_demo_model
      class(atype_demo_model), intent(inout) :: self
    end subroutine deallocate_demo_model_ifc

    subroutine initialise_demo_model_ifc( self, context)
      import atype_demo_model, type_demo_model_context_initialise
      class(atype_demo_model),                          intent(inout) :: self
      type(type_demo_model_context_initialise), target, intent(in   ) :: context
    end subroutine initialise_demo_model_ifc

    subroutine run_demo_model_ifc( self, context)
      import atype_demo_model, type_demo_model_context_run
      class(atype_demo_model),                   intent(inout) :: self
      type(type_demo_model_context_run), target, intent(in   ) :: context
    end subroutine run_demo_model_ifc

    subroutine remap_demo_model_ifc( self, context)
      import atype_demo_model, type_demo_model_context_remap
      class(atype_demo_model),                     intent(inout) :: self
      type(type_demo_model_context_remap), target, intent(in   ) :: context
    end subroutine remap_demo_model_ifc

  end interface

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  interface

    module subroutine allocate_model_abs( self, context)
      class(atype_demo_model),                     intent(inout) :: self
      class(atype_model_context_allocate), target, intent(in   ) :: context
    end subroutine allocate_model_abs

    module subroutine deallocate_model( self)
      class(atype_demo_model), intent(inout) :: self
    end subroutine deallocate_model

    module subroutine initialise_model_abs( self, context)
      class(atype_demo_model),                       intent(inout) :: self
      class(atype_model_context_initialise), target, intent(in   ) :: context
    end subroutine initialise_model_abs

    module subroutine run_model_abs( self, context)
      class(atype_demo_model),                intent(inout) :: self
      class(atype_model_context_run), target, intent(in   ) :: context
    end subroutine run_model_abs

    module subroutine remap_model_abs( self, context)
      class(atype_demo_model),                  intent(inout) :: self
      class(atype_model_context_remap), target, intent(in   ) :: context
    end subroutine remap_model_abs

    module function ct_allocate( name, region_name, mesh, nz) result( context)
      character(len=*),           intent(in) :: name
      character(len=*),           intent(in) :: region_name
      type(type_mesh), target,    intent(in) :: mesh
      integer,                    intent(in) :: nz
      type(type_demo_model_context_allocate) :: context
    end function ct_allocate

    module function ct_initialise( H0, till_friction_angle_uniform, beta_sq_uniform) result( context)
      real(dp),                     intent(in) :: H0
      real(dp),                     intent(in) :: till_friction_angle_uniform
      real(dp),                     intent(in) :: beta_sq_uniform
      type(type_demo_model_context_initialise) :: context
    end function ct_initialise

    module function ct_run( H_new, dH) result( context)
      real(dp),              intent(in) :: H_new
      real(dp),              intent(in) :: dH
      type(type_demo_model_context_run) :: context
    end function ct_run

    module function ct_remap( mesh_new) result( context)
      type(type_mesh), target, intent(in) :: mesh_new
      type(type_demo_model_context_remap) :: context
    end function ct_remap

  end interface

end module demo_model_basic
