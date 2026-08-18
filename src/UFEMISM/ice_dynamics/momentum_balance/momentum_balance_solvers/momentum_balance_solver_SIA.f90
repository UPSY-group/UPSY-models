module momentum_balance_solver_SIA

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_data, only: atype_ice_velocity_model_data
  use parameters, only: grav, ice_density, NaN
  use reallocate_mod, only: reallocate_bounds
  use constitutive_equation, only: calc_ice_rheology_Glen
  use mesh_disc_apply_operators, only: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, map_a_b_3D, ddx_a_a_2D, ddy_a_a_2D
  use mesh_zeta, only: integrate_from_zeta_is_one_to_zeta_is_zetap
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use momentum_balance_solver_plain_SIA, only: type_momentum_balance_solver_plain_SIA

  implicit none

  private

  public :: type_momentum_balance_solver_SIA

  type, extends(atype_momentum_balance_solver) :: type_momentum_balance_solver_SIA

      type(type_momentum_balance_solver_plain_SIA) :: solver

    contains

      ! Procedures for model memory management and operation
      procedure, public :: allocate_momentum_balance_solver   => momentum_balance_solver_SIA_allocate
      procedure, public :: deallocate_momentum_balance_solver => momentum_balance_solver_SIA_deallocate
      procedure, public :: initialise_momentum_balance_solver => momentum_balance_solver_SIA_initialise
      procedure, public :: run_momentum_balance_solver        => momentum_balance_solver_SIA_run
      procedure, public :: remap_momentum_balance_solver      => momentum_balance_solver_SIA_remap

      procedure, public :: get_momentum_balance_solver_name

  end type type_momentum_balance_solver_SIA

contains

  subroutine momentum_balance_solver_SIA_allocate( self, region_name, mesh)

    ! In/output variables:
    class(type_momentum_balance_solver_SIA), intent(inout) :: self
    character(len=*),                                intent(in   ) :: region_name
    type(type_mesh), target,                         intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SIA_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is specific to the SIA momentum balance solver
    call self%solver%allocate( region_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIA_allocate

  subroutine momentum_balance_solver_SIA_deallocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SIA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SIA_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to the SIA momentum balance solver
    call self%solver%deallocate()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIA_deallocate

  subroutine momentum_balance_solver_SIA_initialise( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SIA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SIA_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise all the stuff that is specific to the SIA momentum balance solver
    call self%solver%initialise()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIA_initialise

  subroutine momentum_balance_solver_SIA_run( self, ice, geom, vel)

    ! In/output variables:
    class(type_momentum_balance_solver_SIA), intent(inout) :: self
    class(atype_ice_model_data),                     intent(inout) :: ice
    class(atype_ice_geometry_model_data),            intent(in   ) :: geom
    class(atype_ice_velocity_model_data),            intent(inout) :: vel

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_momentum_balance_solver_SIA'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is specific to the SIA momentum balance solver
    call self%solver%run( ice, geom, vel)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIA_run

  subroutine momentum_balance_solver_SIA_remap( self, mesh_new)

    ! In/output variables:
    class(type_momentum_balance_solver_SIA), intent(inout) :: self
    type(type_mesh), target,                         intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SIA_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap all the stuff that is specific to the SIA momentum balance solver
    call self%solver%remap( mesh_new)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIA_remap

  function get_momentum_balance_solver_name( self) result( model_name)
    class(type_momentum_balance_solver_SIA), intent(in) :: self
    character(len=:), allocatable :: model_name
    model_name = 'SIA'
  end function get_momentum_balance_solver_name

end module momentum_balance_solver_SIA
