module demo_model_c

  use precisions, only: dp
  use parameters, only: pi, NaN
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use demo_model_basic, only: atype_demo_model
  use mpi_f08, only: MPI_WIN
  use demo_model_b, only: type_demo_model_b

  implicit none

  private

  public :: type_demo_model_c

  type, extends(atype_demo_model) :: type_demo_model_c

      class(type_demo_model_b), allocatable :: b_in_c

      ! Some additional ice-model-esque data fields specific to demo_model_c
      real(dp), dimension(:), contiguous, pointer :: test_field_c => null()
      type(MPI_WIN) :: wtest_field_c

    contains

      ! Overriding deferred procedures for model memory management and operation
      procedure, public :: allocate_demo_model   => demo_model_c_allocate
      procedure, public :: deallocate_demo_model => demo_model_c_deallocate
      procedure, public :: initialise_demo_model => demo_model_c_initialise
      procedure, public :: run_demo_model        => demo_model_c_run
      procedure, public :: remap_demo_model      => demo_model_c_remap

      procedure, public :: get_demo_model_name

  end type type_demo_model_c

contains

  subroutine demo_model_c_allocate( self, nz)

    ! In/output variables:
    class(type_demo_model_c), intent(inout) :: self
    integer,                  intent(in   ) :: nz

    ! Local variables:
    character(len=*), parameter :: routine_name = 'demo_model_c_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is specific to demo_model_c
    allocate( self%b_in_c)
    call self%b_in_c%allocate( self%region_name(), self%mesh, nz)

    call self%create_field( self%test_field_c, self%wtest_field_c, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'test_field_c', &
      long_name = 'test field for demo model c', &
      units     = '', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_c_allocate

  subroutine demo_model_c_deallocate( self)

    ! In/output variables:
    class(type_demo_model_c), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'demo_model_c_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to demo_model_c
    call self%b_in_c%deallocate()
    nullify( self%test_field_c)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_c_deallocate

  subroutine demo_model_c_initialise( self, H0, till_friction_angle_uniform, beta_sq_uniform)

    ! In/output variables:
    class(type_demo_model_c), intent(inout) :: self
    real(dp),                 intent(in   ) :: H0
    real(dp),                 intent(in   ) :: till_friction_angle_uniform
    real(dp),                 intent(in   ) :: beta_sq_uniform

    ! Local variables:
    character(len=*), parameter :: routine_name = 'demo_model_c_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise all the stuff that is specific to demo_model_c
    call self%b_in_c%initialise( H0, till_friction_angle_uniform, beta_sq_uniform)
    self%test_field_c( self%mesh%vi1: self%mesh%vi2) = 13.37_dp

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_c_initialise

  subroutine demo_model_c_run( self, H_new, dH)

    ! In/output variables:
    class(type_demo_model_c), intent(inout) :: self
    real(dp),                 intent(in   ) :: H_new
    real(dp),                 intent(in   ) :: dH

    ! Local variables:
    character(len=*), parameter :: routine_name = 'demo_model_c_run'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is specific to demo model b
    call self%b_in_c%run( H_new, dH)
    self%test_field_c( self%mesh%vi1: self%mesh%vi2) = 42._dp

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_c_run

  subroutine demo_model_c_remap( self, mesh_new)

    ! In/output variables:
    class(type_demo_model_c), intent(inout) :: self
    type(type_mesh), target,  intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'demo_model_c_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap all the stuff that is specific to demo model b
    call self%b_in_c%remap( mesh_new)
    call self%remap_field( mesh_new, 'test_field_c', self%test_field_c)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine demo_model_c_remap

  function get_demo_model_name( self) result( demo_model_name)
    class(type_demo_model_c), intent(in) :: self
    character(len=:), allocatable :: demo_model_name
    demo_model_name = 'c'
  end function get_demo_model_name

end module demo_model_c
