module momentum_balance_solver_none

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
  use bed_roughness_model_types, only: type_bed_roughness_model

  implicit none

  private

  public :: type_momentum_balance_solver_none

  type, extends(atype_momentum_balance_solver) :: type_momentum_balance_solver_none

    contains

      ! Procedures for model memory management and operation
      procedure, public :: allocate_momentum_balance_solver   => momentum_balance_solver_none_allocate
      procedure, public :: deallocate_momentum_balance_solver => momentum_balance_solver_none_deallocate
      procedure, public :: initialise_momentum_balance_solver => momentum_balance_solver_none_initialise
      procedure, public :: run_momentum_balance_solver        => momentum_balance_solver_none_run
      procedure, public :: set_velocities_to_solver_results   => momentum_balance_solver_none_set_velocities
      procedure, public :: remap_momentum_balance_solver      => momentum_balance_solver_none_remap

      procedure, public :: get_momentum_balance_solver_name

  end type type_momentum_balance_solver_none

contains

  subroutine momentum_balance_solver_none_allocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_none), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_none_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is specific to the 'none' momentum balance solver

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_none_allocate

  subroutine momentum_balance_solver_none_deallocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_none), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_none_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to the 'none' momentum balance solver

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_none_deallocate

  subroutine momentum_balance_solver_none_initialise( self)

    ! In/output variables:
    class(type_momentum_balance_solver_none), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_none_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise all the stuff that is specific to the 'none' momentum balance solver

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_none_initialise

  subroutine momentum_balance_solver_none_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)

    ! In/output variables:
    class(type_momentum_balance_solver_none), intent(inout) :: self
    class(atype_ice_model_data),              intent(inout) :: ice
    class(atype_ice_geometry_model_data),     intent(in   ) :: geom
    type(type_bed_roughness_model),           intent(in   ) :: bed_roughness
    integer,  dimension(:  ), optional,       intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:  ), optional,       intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(:  ), optional,       intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    integer,  dimension(:,:), optional,       intent(in   ) :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:,:), optional,       intent(in   ) :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
    real(dp), dimension(:,:), optional,       intent(in   ) :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_momentum_balance_solver_none'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is specific to the 'none' momentum balance solver

    self%n_visc_its = 0
    self%n_Axb_its  = 0

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_none_run

  subroutine momentum_balance_solver_none_set_velocities( self, ice, vel)

    ! In/output variables:
    class(type_momentum_balance_solver_none), intent(in   ) :: self
    class(atype_ice_model_data),              intent(inout) :: ice
    class(atype_ice_velocity_model_data),     intent(inout) :: vel

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_none_set_velocities'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_none_set_velocities

  subroutine momentum_balance_solver_none_remap( self, mesh_old, mesh_new)

    ! In/output variables:
    class(type_momentum_balance_solver_none), intent(inout) :: self
    type(type_mesh),                          intent(in   ) :: mesh_old
    type(type_mesh), target,                  intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_none_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap all the stuff that is specific to the 'none' momentum balance solver

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_none_remap

  function get_momentum_balance_solver_name( self) result( model_name)
    class(type_momentum_balance_solver_none), intent(in) :: self
    character(len=:), allocatable :: model_name
    model_name = 'none'
  end function get_momentum_balance_solver_name

end module momentum_balance_solver_none
