module momentum_balance_solver_dummy

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_data, only: atype_ice_velocity_model_data

  implicit none

  private

  public :: type_momentum_balance_solver_dummy

  ! Temporary concrete extension, to be used until all momentum balance solvers have been converted
  ! to extensions of the abstract base type
  type, extends(atype_momentum_balance_solver) :: type_momentum_balance_solver_dummy

    contains

      procedure, public :: allocate_momentum_balance_solver   => momentum_balance_solver_dummy_allocate
      procedure, public :: deallocate_momentum_balance_solver => momentum_balance_solver_dummy_deallocate
      procedure, public :: initialise_momentum_balance_solver => momentum_balance_solver_dummy_initialise
      procedure, public :: run_momentum_balance_solver        => momentum_balance_solver_dummy_run
      procedure, public :: remap_momentum_balance_solver      => momentum_balance_solver_dummy_remap

      procedure, public :: get_momentum_balance_solver_name

  end type type_momentum_balance_solver_dummy

contains

  subroutine momentum_balance_solver_dummy_allocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_dummy), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_dummy_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is specific to the dummy momentum balance solver


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_dummy_allocate

  subroutine momentum_balance_solver_dummy_deallocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_dummy), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_dummy_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to momentum balance solver dummy


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_dummy_deallocate

  subroutine momentum_balance_solver_dummy_initialise( self)

    ! In/output variables:
    class(type_momentum_balance_solver_dummy), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_dummy_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise all the stuff that is specific to momentum balance solver dummy


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_dummy_initialise

  subroutine momentum_balance_solver_dummy_run( self, ice, geom, vel)

    ! In/output variables:
    class(type_momentum_balance_solver_dummy), intent(inout) :: self
    class(atype_ice_model_data),               intent(inout) :: ice
    class(atype_ice_geometry_model_data),      intent(in   ) :: geom
    class(atype_ice_velocity_model_data),      intent(inout) :: vel

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_dummy_run'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is specific to momentum balance solver dummy


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_dummy_run

  subroutine momentum_balance_solver_dummy_remap( self, mesh_new)

    ! In/output variables:
    class(type_momentum_balance_solver_dummy), intent(inout) :: self
    type(type_mesh), target,              intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_dummy_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap all the stuff that is specific to momentum balance solver dummy


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_dummy_remap

  function get_momentum_balance_solver_name( self) result( momentum_balance_solver_name)
    class(type_momentum_balance_solver_dummy), intent(in) :: self
    character(len=:), allocatable :: momentum_balance_solver_name
    momentum_balance_solver_name = 'dummy'
  end function get_momentum_balance_solver_name

end module momentum_balance_solver_dummy