module demo_model_b

  use precisions, only: dp
  use parameters, only: pi, NaN
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use demo_model_state, only: type_demo_model_state
  use demo_model_basic, only: atype_demo_model, type_demo_model_context_allocate, &
    type_demo_model_context_initialise, type_demo_model_context_run, &
    type_demo_model_context_remap
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_demo_model_b

  type, extends(atype_demo_model) :: type_demo_model_b

      ! Some additional ice-model-esque data fields specific to demo_model_b
      real(dp), dimension(:  ), contiguous, pointer :: beta_sq => null()
      type(MPI_WIN) :: wbeta_sq

    contains

      procedure, public :: allocate_demo_model   => allocate_demo_model_b_abs
      procedure, public :: deallocate_demo_model => deallocate_demo_model_b
      procedure, public :: initialise_demo_model => initialise_demo_model_b_abs
      procedure, public :: run_demo_model        => run_demo_model_b_abs
      procedure, public :: remap_demo_model      => remap_demo_model_b_abs

  end type type_demo_model_b

contains

  subroutine allocate_demo_model_b_abs( self, context)

    ! In/output variables:
    class(type_demo_model_b),                       intent(inout) :: self
    type(type_demo_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_demo_model_b_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call allocate_demo_model_b( self%mesh, self)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_demo_model_b_abs

  subroutine allocate_demo_model_b( mesh, self)

    ! In/output variables:
    type(type_mesh),         intent(in   ) :: mesh
    type(type_demo_model_b), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_demo_model_b'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%create_field( self%beta_sq, self%wbeta_sq, &
      mesh, Arakawa_grid%a(), &
      name      = 'beta_sq', &
      long_name = 'bed friction coefficient', &
      units     = 'Pa m^−1/m yr^1/m', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_demo_model_b

  subroutine deallocate_demo_model_b( self)

    ! In/output variables:
    class(type_demo_model_b), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'deallocate_demo_model_b'

    ! Add routine to call stack
    call init_routine( routine_name)

    nullify( self%beta_sq)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_demo_model_b

  subroutine initialise_demo_model_b_abs( self, context)

    ! In/output variables:
    class(type_demo_model_b),                         intent(inout) :: self
    type(type_demo_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_demo_model_b_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call initialise_demo_model_b( self%mesh, self, context%beta_sq_uniform)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_demo_model_b_abs

  subroutine initialise_demo_model_b( mesh, self, beta_sq_uniform)

    ! In/output variables:
    type(type_mesh),         intent(in   ) :: mesh
    type(type_demo_model_b), intent(inout) :: self
    real(dp),                intent(in   ) :: beta_sq_uniform

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_demo_model_b'

    ! Add routine to call stack
    call init_routine( routine_name)

    self%beta_sq( mesh%vi1: mesh%vi2) = beta_sq_uniform

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_demo_model_b

  subroutine run_demo_model_b_abs( self, context)

    ! In/output variables:
    class(type_demo_model_b),                  intent(inout) :: self
    type(type_demo_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_demo_model_b_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call run_demo_model_b( self%mesh, self%s, context%H_new)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_demo_model_b_abs

  subroutine run_demo_model_b( mesh, demo, H_new)

    ! In/output variables:
    type(type_mesh),             intent(in   ) :: mesh
    type(type_demo_model_state), intent(inout) :: demo
    real(dp),                    intent(in   ) :: H_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_demo_model_b'

    ! Add routine to call stack
    call init_routine( routine_name)

    demo%H( mesh%vi1: mesh%vi2) = H_new

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_demo_model_b

  subroutine remap_demo_model_b_abs( self, context)

    ! In/output variables:
    class(type_demo_model_b),                    intent(inout) :: self
    type(type_demo_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_demo_model_b_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call remap_demo_model_b( self, context%mesh_new)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_demo_model_b_abs

  subroutine remap_demo_model_b( self, mesh_new)

    ! In/output variables:
    type(type_demo_model_b), intent(inout) :: self
    type(type_mesh),         intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_demo_model_b'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%remap_field( mesh_new, 'beta_sq', self%beta_sq)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_demo_model_b

end module demo_model_b
