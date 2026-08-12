module ice_velocity_model_dummy

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use ice_velocity_model_basic, only: atype_ice_velocity_model
  use mesh_types, only: type_mesh

  implicit none

  private

  public :: type_ice_velocity_model_dummy

  ! Temporary concrete extension, to be used until all velocity models have been converted
  ! to extensions of the abstract base type
  type, extends(atype_ice_velocity_model) :: type_ice_velocity_model_dummy

    contains

      procedure, public :: allocate_ice_velocity_model   => ice_velocity_model_dummy_allocate
      procedure, public :: deallocate_ice_velocity_model => ice_velocity_model_dummy_deallocate
      procedure, public :: initialise_ice_velocity_model => ice_velocity_model_dummy_initialise
      procedure, public :: run_ice_velocity_model        => ice_velocity_model_dummy_run
      procedure, public :: remap_ice_velocity_model      => ice_velocity_model_dummy_remap

      procedure, public :: get_ice_velocity_model_name

  end type type_ice_velocity_model_dummy

contains

  subroutine ice_velocity_model_dummy_allocate( self)

    ! In/output variables:
    class(type_ice_velocity_model_dummy), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'ice_velocity_model_dummy_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is specific to the dummy ice_velocity model


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ice_velocity_model_dummy_allocate

  subroutine ice_velocity_model_dummy_deallocate( self)

    ! In/output variables:
    class(type_ice_velocity_model_dummy), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'ice_velocity_model_dummy_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to ice_velocity model dummy


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ice_velocity_model_dummy_deallocate

  subroutine ice_velocity_model_dummy_initialise( self)

    ! In/output variables:
    class(type_ice_velocity_model_dummy), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'ice_velocity_model_dummy_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise all the stuff that is specific to ice_velocity model dummy


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ice_velocity_model_dummy_initialise

  subroutine ice_velocity_model_dummy_run( self)

    ! In/output variables:
    class(type_ice_velocity_model_dummy), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_ice_velocity_model_dummy'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is specific to ice_velocity model dummy


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ice_velocity_model_dummy_run

  subroutine ice_velocity_model_dummy_remap( self, mesh_new)

    ! In/output variables:
    class(type_ice_velocity_model_dummy), intent(inout) :: self
    type(type_mesh), target,              intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'ice_velocity_model_dummy_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap all the stuff that is specific to ice_velocity model dummy


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ice_velocity_model_dummy_remap

  function get_ice_velocity_model_name( self) result( ice_velocity_model_name)
    class(type_ice_velocity_model_dummy), intent(in) :: self
    character(len=:), allocatable :: ice_velocity_model_name
    ice_velocity_model_name = 'dummy'
  end function get_ice_velocity_model_name

end module ice_velocity_model_dummy