module momentum_balance_solver_SSA

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use parameters, only: grav, ice_density, NaN
  use reallocate_mod, only: reallocate_bounds
  use constitutive_equation, only: calc_ice_rheology_Glen
  use mesh_disc_apply_operators, only: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, map_a_b_3D, ddx_a_a_2D, ddy_a_a_2D
  use mesh_zeta, only: integrate_from_zeta_is_one_to_zeta_is_zetap
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use momentum_balance_solver_plain_SSA, only: type_momentum_balance_solver_plain_SSA
  use bed_roughness_model_types, only: type_bed_roughness_model

  implicit none

  private

  public :: type_momentum_balance_solver_SSA

  type, extends(atype_momentum_balance_solver) :: type_momentum_balance_solver_SSA

      type(type_momentum_balance_solver_plain_SSA) :: solver

    contains

      ! Procedures for model memory management and operation
      procedure, public :: allocate_momentum_balance_solver   => momentum_balance_solver_SSA_allocate
      procedure, public :: deallocate_momentum_balance_solver => momentum_balance_solver_SSA_deallocate
      procedure, public :: initialise_momentum_balance_solver => momentum_balance_solver_SSA_initialise
      procedure, public :: run_momentum_balance_solver        => momentum_balance_solver_SSA_run
      procedure, public :: remap_momentum_balance_solver      => momentum_balance_solver_SSA_remap

      procedure, public :: get_momentum_balance_solver_name

  end type type_momentum_balance_solver_SSA

contains

  subroutine momentum_balance_solver_SSA_allocate( self, region_name, mesh)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self
    character(len=*),                        intent(in   ) :: region_name
    type(type_mesh), target,                 intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is specific to the SSA momentum balance solver
    call self%solver%allocate( region_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_allocate

  subroutine momentum_balance_solver_SSA_deallocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to the SSA momentum balance solver
    call self%solver%deallocate()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_deallocate

  subroutine momentum_balance_solver_SSA_initialise( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise all the stuff that is specific to the SSA momentum balance solver
    call self%solver%initialise()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_initialise

  subroutine momentum_balance_solver_SSA_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self
    class(atype_ice_model_data),             intent(inout) :: ice
    class(atype_ice_geometry_model_data),    intent(in   ) :: geom
    type(type_bed_roughness_model),          intent(in   ) :: bed_roughness
    integer,  dimension(:), optional,        intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:), optional,        intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(:), optional,        intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_momentum_balance_solver_SSA'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is specific to the SSA momentum balance solver
    call self%solver%run( ice, geom, bed_roughness, &
      BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)

    self%n_visc_its = self%solver%n_visc_its
    self%n_Axb_its  = self%solver%n_Axb_its

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_run

  subroutine momentum_balance_solver_SSA_remap( self, mesh_old, mesh_new)

    ! In/output variables:
    class(type_momentum_balance_solver_SSA), intent(inout) :: self
    type(type_mesh),                         intent(in   ) :: mesh_old
    type(type_mesh), target,                 intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SSA_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap all the stuff that is specific to the SSA momentum balance solver
    call self%solver%remap( mesh_old, mesh_new)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SSA_remap

  function get_momentum_balance_solver_name( self) result( model_name)
    class(type_momentum_balance_solver_SSA), intent(in) :: self
    character(len=:), allocatable :: model_name
    model_name = 'SSA'
  end function get_momentum_balance_solver_name

end module momentum_balance_solver_SSA
