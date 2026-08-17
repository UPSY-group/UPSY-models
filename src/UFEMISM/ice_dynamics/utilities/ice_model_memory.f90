module ice_model_memory
  !< Routines for administrating the memory for the ice model data.

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_main, only: type_ice_model
  use ice_velocity_model, only: create_ice_velocity_model

  implicit none

  private

  public :: allocate_ice_model

contains

  subroutine allocate_ice_model( mesh, ice, region_name)

    ! In- and output variables
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(  out) :: ice
    character(len=*),     intent(in   ) :: region_name

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_ice_model'

    ! Add routine to path
    call init_routine( routine_name)

    allocate( ice%geom)
    call ice%geom%allocate( region_name, mesh)

    call create_ice_velocity_model( ice%vel, C%choice_stress_balance_approximation)
    call ice%vel%allocate( region_name, mesh)

    ! Geometry changes
    allocate( ice%dHi ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHb ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHs ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHib( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Rates of change
    allocate( ice%dHi_dt         ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHb_dt         ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHs_dt         ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHib_dt        ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHi_dt_raw     ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%dHi_dt_residual( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Target quantities
    allocate( ice%dHi_dt_target   ( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Masks
    allocate( ice%mask_SGD               ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_noice             ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( ice%mask_ROI               ( mesh%vi1:mesh%vi2), source = 0)
    allocate( ice%basin_ID               ( mesh%vi1:mesh%vi2), source = 0)

    ! === Terrain-following coordinate zeta gradients ===
    ! ===================================================

    ! Gradients of the terrain-following (i.e. ice-geometry-dependent) vertical coordinate zeta

    ! On the ak-grid (vertices, vertically regular)
    allocate( ice%dzeta_dt_ak   ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%dzeta_dx_ak   ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%dzeta_dy_ak   ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%dzeta_dz_ak   ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%d2zeta_dx2_ak ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%d2zeta_dxdy_ak( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%d2zeta_dy2_ak ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)

    ! On the bk-grid (triangles, vertically regular)
    allocate( ice%dzeta_dx_bk   ( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( ice%dzeta_dy_bk   ( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( ice%dzeta_dz_bk   ( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( ice%d2zeta_dx2_bk ( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( ice%d2zeta_dxdy_bk( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( ice%d2zeta_dy2_bk ( mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)

    ! On the bks-grid (triangles, vertically staggered)
    allocate( ice%dzeta_dx_bks   ( mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( ice%dzeta_dy_bks   ( mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( ice%dzeta_dz_bks   ( mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( ice%d2zeta_dx2_bks ( mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( ice%d2zeta_dxdy_bks( mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)
    allocate( ice%d2zeta_dy2_bks ( mesh%ti1:mesh%ti2,mesh%nz-1), source = 0._dp)

    ! === Thermodynamics and rheology ===
    ! ===================================

    ! Ice temperatures
    allocate( ice%Ti    ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%Ti_pmp( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%Ti_hom( mesh%vi1:mesh%vi2        ), source = 0._dp)

    ! Physical quantities
    allocate( ice%Cpi( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%Ki ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)

    ! Heating
    allocate( ice%internal_heating  ( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%frictional_heating( mesh%vi1:mesh%vi2        ), source = 0._dp)

    ! Glen's flow law factor
    allocate( ice%A_flow( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)

    ! == Ice flow regime ==
    ! =====================

    allocate( ice%divQ   ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%Qspill ( mesh%vi1:mesh%vi2), source = 0._dp)

    ! == Basal hydrology ==
    ! =====================

    ! Basal hydrology
    allocate( ice%pore_water_pressure  ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%overburden_pressure  ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%effective_pressure   ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%pore_water_likelihood( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%pore_water_fraction  ( mesh%vi1:mesh%vi2), source = 0._dp)

    ! == Basal sliding ==
    ! ===================

    ! Basal friction and shear stress
    allocate( ice%till_yield_stress         ( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%basal_friction_coefficient( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%basal_shear_stress        ( mesh%ti1:mesh%ti2), source = 0._dp)

    ! == Geothermal heat ==
    ! =====================

    allocate( ice%geothermal_heat_flux( mesh%vi1:mesh%vi2), source = 0._dp)

    ! === Ice thickness time stepping ===
    ! ===================================

    ! Predicted model state at next time step
    allocate( ice%Hi_prev( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( ice%Hi_next( mesh%vi1:mesh%vi2), source = 0._dp)

    ! === Ice temperature time stepping ===
    ! =====================================

    ! Predicted model state at next time step
    allocate( ice%Ti_prev( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( ice%Ti_next( mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_ice_model

end module ice_model_memory
