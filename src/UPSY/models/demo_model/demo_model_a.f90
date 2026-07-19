module demo_model_a

  use precisions, only: dp
  use parameters, only: pi, NaN
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use demo_model_basic, only: atype_demo_model, type_demo_model_context_remap
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_demo_model_a

  type, extends(atype_demo_model) :: type_demo_model_a

      ! Some additional ice-model-esque data fields specific to demo_model_a
      real(dp), dimension(:  ), contiguous, pointer :: till_friction_angle => null()
      type(MPI_WIN) :: wtill_friction_angle

    contains

      procedure, public :: allocate   => demo_model_a_allocate
      procedure, public :: deallocate => demo_model_a_deallocate
      procedure, public :: initialise => demo_model_a_initialise
      procedure, public :: run        => demo_model_a_run
      procedure, public :: remap_demo_model      => remap_demo_model_a_abs

  end type type_demo_model_a

contains

  subroutine demo_model_a_allocate( self, region_name, mesh, nz)

    ! In/output variables:
    class(type_demo_model_a), intent(inout) :: self
    character(len=*),         intent(in   ) :: region_name
    type(type_mesh), target,  intent(in   ) :: mesh
    integer,                  intent(in   ) :: nz

    ! Local variables:
    character(len=*), parameter :: routine_name = 'demo_model_a_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all demo models
    call self%allocate_demo_model( 'demo_a', region_name, mesh, nz)

    ! Allocate all the stuff that is specific to demo model a

    call self%create_field( self%till_friction_angle, self%wtill_friction_angle, &
      mesh, Arakawa_grid%a(), &
      name      = 'till_friction_angle', &
      long_name = 'till friction angle', &
      units     = 'degrees', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_a_allocate

  subroutine demo_model_a_deallocate( self)

    ! In/output variables:
    class(type_demo_model_a), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'demo_model_a_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is common to all demo models
    call self%deallocate_demo_model()

    ! Deallocate all the stuff that is specific to demo model a

    nullify( self%till_friction_angle)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_a_deallocate

  subroutine demo_model_a_initialise( self, H0, till_friction_angle_uniform, beta_sq_uniform)

    ! In/output variables:
    class(type_demo_model_a), intent(inout) :: self
    real(dp),                 intent(in   ) :: H0
    real(dp),                 intent(in   ) :: till_friction_angle_uniform
    real(dp),                 intent(in   ) :: beta_sq_uniform

    ! Local variables:
    character(len=*), parameter :: routine_name = 'demo_model_a_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise all the stuff that is common to all demo models
    call self%initialise_demo_model( H0)

    ! Initialise all the stuff that is specific to demo model a

    self%till_friction_angle( self%mesh%vi1: self%mesh%vi2) = till_friction_angle_uniform

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_a_initialise

  subroutine demo_model_a_run( self, H_new, dH)

    ! In/output variables:
    class(type_demo_model_a), intent(inout) :: self
    real(dp),                 intent(in   ) :: H_new
    real(dp),                 intent(in   ) :: dH

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'demo_model_a_run'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is common to all demo models
    call self%run_demo_model()

    ! Run all the stuff that is specific to demo model a

    self%H( self%mesh%vi1: self%mesh%vi2) = self%H( self%mesh%vi1: self%mesh%vi2) + dH

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_a_run

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
