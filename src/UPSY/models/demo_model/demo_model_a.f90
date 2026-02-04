module demo_model_a

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

  public :: type_demo_model_a

  type, extends(atype_demo_model) :: type_demo_model_a

      ! Some additional ice-model-esque data fields specific to demo_model_a
      real(dp), dimension(:  ), contiguous, pointer :: till_friction_angle => null()
      type(MPI_WIN) :: wtill_friction_angle

    contains

      procedure, public :: allocate_demo_model   => allocate_demo_model_a_abs
      procedure, public :: deallocate_demo_model => deallocate_demo_model_a
      procedure, public :: initialise_demo_model => initialise_demo_model_a_abs
      procedure, public :: run_demo_model        => run_demo_model_a_abs
      procedure, public :: remap_demo_model      => remap_demo_model_a_abs

  end type type_demo_model_a

contains

  subroutine allocate_demo_model_a_abs( self, context)

    ! In/output variables:
    class(type_demo_model_a),                       intent(inout) :: self
    type(type_demo_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_demo_model_a_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call allocate_demo_model_a( self%mesh, self)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_demo_model_a_abs

  subroutine allocate_demo_model_a( mesh, self)

    ! In/output variables:
    type(type_mesh),         intent(in   ) :: mesh
    type(type_demo_model_a), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_demo_model_a'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%create_field( self%till_friction_angle, self%wtill_friction_angle, &
      mesh, Arakawa_grid%a(), &
      name      = 'till_friction_angle', &
      long_name = 'till friction angle', &
      units     = 'degrees', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_demo_model_a

  subroutine deallocate_demo_model_a( self)

    ! In/output variables:
    class(type_demo_model_a), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'deallocate_demo_model_a'

    ! Add routine to call stack
    call init_routine( routine_name)

    nullify( self%till_friction_angle)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_demo_model_a

  subroutine initialise_demo_model_a_abs( self, context)

    ! In/output variables:
    class(type_demo_model_a),                         intent(inout) :: self
    type(type_demo_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_demo_model_a_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call initialise_demo_model_a( self%mesh, self, context%till_friction_angle_uniform)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_demo_model_a_abs

  subroutine initialise_demo_model_a( mesh, self, till_friction_angle_uniform)

    ! In/output variables:
    type(type_mesh),         intent(in   ) :: mesh
    type(type_demo_model_a), intent(inout) :: self
    real(dp),                intent(in   ) :: till_friction_angle_uniform

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_demo_model_a'

    ! Add routine to call stack
    call init_routine( routine_name)

    self%till_friction_angle( mesh%vi1: mesh%vi2) = till_friction_angle_uniform

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_demo_model_a

  subroutine run_demo_model_a_abs( self, context)

    ! In/output variables:
    class(type_demo_model_a),                  intent(inout) :: self
    type(type_demo_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_demo_model_a_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call run_demo_model_a( self%mesh, self%s, context%dH)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_demo_model_a_abs

  subroutine run_demo_model_a( mesh, demo, dH)

    ! In/output variables:
    type(type_mesh),             intent(in   ) :: mesh
    type(type_demo_model_state), intent(inout) :: demo
    real(dp),                    intent(in   ) :: dH

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_demo_model_a'

    ! Add routine to call stack
    call init_routine( routine_name)

    demo%H( mesh%vi1: mesh%vi2) = demo%H( mesh%vi1: mesh%vi2) + dH

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_demo_model_a

  subroutine remap_demo_model_a_abs( self, context)

    ! In/output variables:
    class(type_demo_model_a),                    intent(inout) :: self
    type(type_demo_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_demo_model_a_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call remap_demo_model_a( self, context%mesh_new)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_demo_model_a_abs

  subroutine remap_demo_model_a( self, mesh_new)

    ! In/output variables:
    type(type_demo_model_a), intent(inout) :: self
    type(type_mesh),         intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_demo_model_a'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%remap_field( mesh_new, 'till_friction_angle', self%till_friction_angle)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_demo_model_a

end module demo_model_a
