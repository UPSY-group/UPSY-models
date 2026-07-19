module demo_model_b

  use precisions, only: dp
  use parameters, only: pi, NaN
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use demo_model_basic, only: atype_demo_model
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_demo_model_b

  type, extends(atype_demo_model) :: type_demo_model_b

      ! Some additional ice-model-esque data fields specific to demo_model_b
      real(dp), dimension(:), contiguous, pointer :: beta_sq => null()
      type(MPI_WIN) :: wbeta_sq

    contains

      procedure, public :: allocate   => demo_model_b_allocate
      procedure, public :: deallocate => demo_model_b_deallocate
      procedure, public :: initialise => demo_model_b_initialise
      procedure, public :: run        => demo_model_b_run
      procedure, public :: remap      => demo_model_b_remap

  end type type_demo_model_b

contains

  subroutine demo_model_b_allocate( self, region_name, mesh, nz)

    ! In/output variables:
    class(type_demo_model_b), intent(inout) :: self
    character(len=*),         intent(in   ) :: region_name
    type(type_mesh), target,  intent(in   ) :: mesh
    integer,                  intent(in   ) :: nz

    ! Local variables:
    character(len=*), parameter :: routine_name = 'demo_model_b_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all demo models
    call self%allocate_demo_model( 'demo_b', region_name, mesh, nz)

    ! Allocate all the stuff that is specific to demo model b

    call self%create_field( self%beta_sq, self%wbeta_sq, &
      mesh, Arakawa_grid%a(), &
      name      = 'beta_sq', &
      long_name = 'bed friction coefficient', &
      units     = 'Pa m^−1/m yr^1/m', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_b_allocate

  subroutine demo_model_b_deallocate( self)

    ! In/output variables:
    class(type_demo_model_b), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'demo_model_b_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is common to all demo models
    call self%deallocate_demo_model()

    ! Deallocate all the stuff that is specific to demo model b

    nullify( self%beta_sq)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_b_deallocate

  subroutine demo_model_b_initialise( self, H0, till_friction_angle_uniform, beta_sq_uniform)

    ! In/output variables:
    class(type_demo_model_b), intent(inout) :: self
    real(dp),                 intent(in   ) :: H0
    real(dp),                 intent(in   ) :: till_friction_angle_uniform
    real(dp),                 intent(in   ) :: beta_sq_uniform

    ! Local variables:
    character(len=*), parameter :: routine_name = 'demo_model_b_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise all the stuff that is common to all demo models
    call self%initialise_demo_model( H0)

    ! Initialise all the stuff that is specific to demo model a

    self%beta_sq( self%mesh%vi1: self%mesh%vi2) = beta_sq_uniform

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_b_initialise

  subroutine demo_model_b_run( self, H_new, dH)

    ! In/output variables:
    class(type_demo_model_b), intent(inout) :: self
    real(dp),                 intent(in   ) :: H_new
    real(dp),                 intent(in   ) :: dH

    ! Local variables:
    character(len=*), parameter :: routine_name = 'demo_model_b_run'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is common to all demo models
    call self%run_demo_model()

    ! Run all the stuff that is specific to demo model b

    self%H( self%mesh%vi1: self%mesh%vi2) = H_new

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_b_run

  subroutine demo_model_b_remap( self, mesh_new)

    ! In/output variables:
    class(type_demo_model_b), intent(inout) :: self
    type(type_mesh), target,  intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'demo_model_b_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap all the stuff that is common to all demo models
    call self%remap_demo_model( mesh_new)

    ! Remap all the stuff that is specific to demo model b

    call self%remap_field( mesh_new, 'beta_sq', self%beta_sq)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_b_remap

end module demo_model_b
