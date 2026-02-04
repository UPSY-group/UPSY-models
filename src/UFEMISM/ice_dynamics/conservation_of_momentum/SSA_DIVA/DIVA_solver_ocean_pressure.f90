module DIVA_solver_ocean_pressure

  ! Routines for calculating ice velocities using the Depth-Integrated Viscosity Approximation (DIVA)

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, &
    MPI_LOR, MPI_LOGICAL, MPI_MIN, MPI_MAX, MPI_SUM
  use mpi_basic, only: par
  use precisions, only: dp
  use parameters, only: grav, ice_density, seawater_density
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph_pair
  use graph_pair_creation, only: create_ice_only_graph_pair, deallocate_graph_pair
  use ice_model_types, only: type_ice_model, type_ice_velocity_solver_DIVA, &
    type_ice_velocity_solver_DIVA_graphs
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist
  use netcdf_io_main
  use sliding_laws, only: calc_basal_friction_coefficient
  use mesh_disc_apply_operators, only: map_a_b_2D, map_a_b_3D, ddx_a_b_2D, ddy_a_b_2D, &
    map_b_a_2D, map_b_a_3D, ddx_b_a_2D, ddy_b_a_2D
  use mesh_zeta, only: integrate_from_zeta_is_one_to_zeta_is_zetap, vertical_average
  use constitutive_equation, only: calc_ice_rheology_Glen, calc_effective_viscosity_Glen_3D_uv_only
  use mpi_distributed_memory, only: gather_to_all
  use petsc_basic, only: solve_matrix_equation_CSR_PETSc
  use reallocate_mod, only: reallocate_bounds, reallocate_clean
  use solve_linearised_SSA_DIVA_ocean_pressure, only: solve_SSA_DIVA_linearised_ocean_pressure
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D, map_from_mesh_to_mesh_with_reallocation_3D
  use bed_roughness_model_types, only: type_bed_roughness_model
  use mpi_f08, only: MPI_WIN
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared
  use mesh_graph_mapping, only: map_mesh_triangles_to_graph, map_mesh_vertices_to_graph, &
    map_graph_to_mesh_triangles, map_graph_to_mesh_vertices
  use CSR_matrix_vector_multiplication, only: multiply_csr_matrix_with_vector_1D_wrapper, &
    multiply_csr_matrix_with_vector_2D_wrapper
  use graph_halo_exchange, only: exchange_halos

  implicit none

  private

  public :: solve_DIVA_ocean_pressure

contains

  ! == Main routines

  subroutine solve_DIVA_ocean_pressure( mesh, ice, bed_roughness, DIVA, n_visc_its, n_Axb_its, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    !< Calculate ice velocities by solving the Depth-Integrated Viscosity Approximation

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(inout) :: ice
    type(type_bed_roughness_model),      intent(in   ) :: bed_roughness
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA
    integer,                             intent(  out) :: n_visc_its               ! Number of non-linear viscosity iterations
    integer,                             intent(  out) :: n_Axb_its                ! Number of iterations in iterative solver for linearised momentum balance
    integer,  dimension(:), optional,    intent(in   ) :: BC_prescr_mask_b         ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:), optional,    intent(in   ) :: BC_prescr_u_b            ! Prescribed velocities in the x-direction
    real(dp), dimension(:), optional,    intent(in   ) :: BC_prescr_v_b            ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'solve_DIVA_ocean_pressure'
    logical                             :: grounded_ice_exists
    integer                             :: ierr
    integer,  dimension(:), allocatable :: BC_prescr_mask_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_u_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_v_b_applied
    integer                             :: viscosity_iteration_i
    logical                             :: has_converged
    real(dp)                            :: L2_uv, L2_uv_prev
    real(dp)                            :: uv_min, uv_max
    real(dp)                            :: visc_it_relax_applied
    real(dp)                            :: Glens_flow_law_epsilon_sq_0_applied
    integer                             :: nit_diverg_consec
    integer                             :: n_Axb_its_visc_it

    real(dp), dimension(:), pointer :: d_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! if there is no grounded ice, no need (in fact, no way) to solve the DIVA
    grounded_ice_exists = any( ice%mask_grounded_ice)
    call MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.not. grounded_ice_exists) then
      DIVA%u_vav_b  = 0._dp
      DIVA%v_vav_b  = 0._dp
      DIVA%u_base_b = 0._dp
      DIVA%v_base_b = 0._dp
      DIVA%u_3D_b   = 0._dp
      DIVA%v_3D_b   = 0._dp
      call finalise_routine( routine_name)
      return
    end if

    ! Handle the optional prescribed u,v boundary conditions
    allocate( BC_prescr_mask_b_applied( mesh%ti1:mesh%ti2))
    allocate( BC_prescr_u_b_applied(    mesh%ti1:mesh%ti2))
    allocate( BC_prescr_v_b_applied(    mesh%ti1:mesh%ti2))
    if (present( BC_prescr_mask_b) .or. present( BC_prescr_u_b) .or. present( BC_prescr_v_b)) then
      ! Safety
      if (.not. (present( BC_prescr_mask_b) .and. present( BC_prescr_u_b) .and. present( BC_prescr_v_b))) then
        call crash('need to provide prescribed u,v fields and mask!')
      end if
      BC_prescr_mask_b_applied = BC_prescr_mask_b
      BC_prescr_u_b_applied    = BC_prescr_u_b
      BC_prescr_v_b_applied    = BC_prescr_v_b
    else
      BC_prescr_mask_b_applied = 0
      BC_prescr_u_b_applied    = 0._dp
      BC_prescr_v_b_applied    = 0._dp
    end if

    call setup_DIVA_solver_on_graphs( mesh, ice, DIVA)

    ! DENK DROM
    call save_graph_pair_as_netcdf( trim( C%output_dir)//'/DIVA_graphs', DIVA%DIVA_graphs%graphs)

    ! Calculate the driving stress and ocean back pressure
    call calc_driving_stress( DIVA%DIVA_graphs)
    call calc_ocean_back_pressure( DIVA%DIVA_graphs)

    ! Adaptive relaxation parameter for the viscosity iteration
    L2_uv                               = 1E9_dp
    nit_diverg_consec                   = 0
    visc_it_relax_applied               = C%visc_it_relax
    Glens_flow_law_epsilon_sq_0_applied = C%Glens_flow_law_epsilon_sq_0

    ! Initialise stability info
    n_visc_its = 0
    n_Axb_its  = 0

    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .false.
    viscosity_iteration: do while (.not. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1

      d_loc => DIVA%DIVA_graphs%u_vav_b( DIVA%DIVA_graphs%graphs%graph_b%ni1: DIVA%DIVA_graphs%graphs%graph_b%ni2)
      call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'u_vav_b0')
      d_loc => DIVA%DIVA_graphs%v_vav_b( DIVA%DIVA_graphs%graphs%graph_b%ni1: DIVA%DIVA_graphs%graphs%graph_b%ni2)
      call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'v_vav_b0')

      ! Calculate the horizontal strain rates for the current velocity solution
      call calc_horizontal_strain_rates( DIVA%DIVA_graphs)

      ! Calculate the vertical shear strain rates
      call calc_vertical_shear_strain_rates( mesh, DIVA%DIVA_graphs)

      ! Calculate the effective viscosity for the current velocity solution
      call calc_effective_viscosity( mesh, ice, DIVA%DIVA_graphs, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the F-integrals (Lipscomb et al. (2019), Eq. 30)
      call calc_F_integrals( mesh, DIVA%DIVA_graphs)

      ! Calculate the "effective" friction coefficient (turning the SSA into the DIVA)
      call calc_effective_basal_friction_coefficient( mesh, ice, bed_roughness, DIVA%DIVA_graphs)

      ! Solve the linearised DIVA to calculate a new velocity solution
      call solve_SSA_DIVA_linearised_ocean_pressure( DIVA%DIVA_graphs, &
        DIVA%PETSc_rtol, DIVA%PETSc_abstol, n_Axb_its_visc_it)

      ! Update stability info
      n_Axb_its = n_Axb_its + n_Axb_its_visc_it

      ! Limit velocities for improved stability
      call apply_velocity_limits( DIVA%DIVA_graphs)

      ! Reduce the change between velocity solutions
      call relax_viscosity_iterations( DIVA%DIVA_graphs, visc_it_relax_applied)

      ! Calculate basal velocities
      call calc_basal_velocities( DIVA%DIVA_graphs)

      ! Calculate basal shear stress
      call calc_basal_shear_stress( DIVA%DIVA_graphs)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      L2_uv_prev = L2_uv
      call calc_L2_norm_uv( DIVA%DIVA_graphs, L2_uv)

      ! if the viscosity iteration diverges, lower the relaxation parameter
      if (L2_uv > L2_uv_prev) then
        nit_diverg_consec = nit_diverg_consec + 1
      else
        nit_diverg_consec = 0
      end if
      if (nit_diverg_consec > 2) then
        nit_diverg_consec = 0
        visc_it_relax_applied               = visc_it_relax_applied               * 0.9_dp
        Glens_flow_law_epsilon_sq_0_applied = Glens_flow_law_epsilon_sq_0_applied * 1.2_dp
      end if
      if (visc_it_relax_applied <= 0.05_dp .or. Glens_flow_law_epsilon_sq_0_applied >= 1E-5_dp) then
        if (visc_it_relax_applied < 0.05_dp) then
          call crash('viscosity iteration still diverges even with very low relaxation factor!')
        elseif (Glens_flow_law_epsilon_sq_0_applied > 1E-5_dp) then
          call crash('viscosity iteration still diverges even with very high effective strain rate regularisation!')
        end if
      end if

      ! DENK DROM
      uv_min = minval( DIVA%DIVA_graphs%u_vav_b)
      uv_max = maxval( DIVA%DIVA_graphs%u_vav_b)
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_doUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_doUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      if (par%primary) WRITE(0,*) '    DIVA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], L2_uv = ', L2_uv

      ! if the viscosity iteration has converged, or has reached the maximum allowed number of iterations, stop it.
      has_converged = .false.
      if (L2_uv < C%visc_it_norm_dUV_tol) then
        has_converged = .true.
      end if

      ! if we've reached the maximum allowed number of iterations without converging, throw a warning
      if (viscosity_iteration_i > C%visc_it_nit) then
        if (par%primary) call warning('viscosity iteration failed to converge within {int_01} iterations!', int_01 = C%visc_it_nit)
        exit viscosity_iteration
      end if

      ! ! DENK DROM
      ! if (viscosity_iteration_i == 50) then
      !   d_loc => DIVA%DIVA_graphs%u_vav_b( DIVA%DIVA_graphs%graphs%graph_b%ni1: DIVA%DIVA_graphs%graphs%graph_b%ni2)
      !   call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'u_vav_b')
      !   d_loc => DIVA%DIVA_graphs%v_vav_b( DIVA%DIVA_graphs%graphs%graph_b%ni1: DIVA%DIVA_graphs%graphs%graph_b%ni2)
      !   call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'v_vav_b')
      !   call crash('whoopsiedaisy')
      ! end if

    end do viscosity_iteration

    ! Calculate 3-D ice velocities
    call calc_3D_velocities( mesh, DIVA%DIVA_graphs)

    ! Stability info
    n_visc_its = viscosity_iteration_i

    ! DENK DROM
    d_loc => DIVA%DIVA_graphs%u_vav_b( DIVA%DIVA_graphs%graphs%graph_b%ni1: DIVA%DIVA_graphs%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'u_vav_b')
    d_loc => DIVA%DIVA_graphs%v_vav_b( DIVA%DIVA_graphs%graphs%graph_b%ni1: DIVA%DIVA_graphs%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'v_vav_b')
    call crash('whoopsiedaisy')

    call finalise_DIVA_solver_on_graphs( mesh, DIVA)

    ! DENK DROM
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), DIVA%u_vav_b, 'u_vav_b')
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), DIVA%v_vav_b, 'v_vav_b')

    call crash('fixme!')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_DIVA_ocean_pressure

  subroutine setup_DIVA_solver_on_graphs( mesh, ice, DIVA)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(inout) :: ice
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_DIVA_solver_on_graphs'

    ! Add routine to path
    call init_routine( routine_name)

    call create_ice_only_graph_pair( mesh, ice, DIVA%DIVA_graphs%graphs)
    call allocate_DIVA_solver_on_graphs( mesh, DIVA%DIVA_graphs)
    call map_DIVA_from_mesh_to_graphs( mesh, ice, DIVA)

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 55)

  end subroutine setup_DIVA_solver_on_graphs

  subroutine finalise_DIVA_solver_on_graphs( mesh, DIVA)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'finalise_DIVA_solver_on_graphs'

    ! Add routine to path
    call init_routine( routine_name)

    call map_DIVA_from_graphs_to_mesh( mesh, DIVA)
    call deallocate_graph_pair( DIVA%DIVA_graphs%graphs)
    call deallocate_DIVA_solver_on_graphs( DIVA%DIVA_graphs)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine finalise_DIVA_solver_on_graphs

  subroutine allocate_DIVA_solver_on_graphs( mesh, DIVA)

    ! In/output variables:
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_DIVA_solver_on_graphs'

    ! Add routine to path
    call init_routine( routine_name)

    ! Ice model forcing data
    call allocate_dist_shared( DIVA%Hi_a                        , DIVA%wHi_a                        , DIVA%graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( DIVA%Hs_a                        , DIVA%wHs_a                        , DIVA%graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( DIVA%Ho_a                        , DIVA%wHo_a                        , DIVA%graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( DIVA%A_flow_3D_a                 , DIVA%wA_flow_3D_a                 , DIVA%graphs%graph_a%pai%n_nih, mesh%nz)
    call allocate_dist_shared( DIVA%basal_friction_coefficient_a, DIVA%wbasal_friction_coefficient_a, DIVA%graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( DIVA%fraction_gr_b               , DIVA%wfraction_gr_b               , DIVA%graphs%graph_b%pai%n_nih)

    ! Solution
    call allocate_dist_shared( DIVA%u_vav_b                     , DIVA%wu_vav_b                     , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%v_vav_b                     , DIVA%wv_vav_b                     , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%u_base_b                    , DIVA%wu_base_b                    , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%v_base_b                    , DIVA%wv_base_b                    , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%u_3D_b                      , DIVA%wu_3D_b                      , DIVA%graphs%graph_b%pai%n_nih, mesh%nz)
    call allocate_dist_shared( DIVA%v_3D_b                      , DIVA%wv_3D_b                      , DIVA%graphs%graph_b%pai%n_nih, mesh%nz)

    ! Intermediate data fields
    call allocate_dist_shared( DIVA%du_dx_a                     , DIVA%wdu_dx_a                     , DIVA%graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( DIVA%du_dy_a                     , DIVA%wdu_dy_a                     , DIVA%graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( DIVA%dv_dx_a                     , DIVA%wdv_dx_a                     , DIVA%graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( DIVA%dv_dy_a                     , DIVA%wdv_dy_a                     , DIVA%graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( DIVA%du_dx_b                     , DIVA%wdu_dx_b                     , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%du_dy_b                     , DIVA%wdu_dy_b                     , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%dv_dx_b                     , DIVA%wdv_dx_b                     , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%dv_dy_b                     , DIVA%wdv_dy_b                     , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%du_dz_3D_a                  , DIVA%wdu_dz_3D_a                  , DIVA%graphs%graph_a%pai%n_nih, mesh%nz)
    call allocate_dist_shared( DIVA%dv_dz_3D_a                  , DIVA%wdv_dz_3D_a                  , DIVA%graphs%graph_a%pai%n_nih, mesh%nz)
    call allocate_dist_shared( DIVA%du_dz_3D_b                  , DIVA%wdu_dz_3D_b                  , DIVA%graphs%graph_b%pai%n_nih, mesh%nz)
    call allocate_dist_shared( DIVA%dv_dz_3D_b                  , DIVA%wdv_dz_3D_b                  , DIVA%graphs%graph_b%pai%n_nih, mesh%nz)
    call allocate_dist_shared( DIVA%eta_3D_a                    , DIVA%weta_3D_a                    , DIVA%graphs%graph_a%pai%n_nih, mesh%nz)
    call allocate_dist_shared( DIVA%eta_3D_b                    , DIVA%weta_3D_b                    , DIVA%graphs%graph_b%pai%n_nih, mesh%nz)
    call allocate_dist_shared( DIVA%eta_vav_a                   , DIVA%weta_vav_a                   , DIVA%graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( DIVA%eta_vav_b                   , DIVA%weta_vav_b                   , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%N_a                         , DIVA%wN_a                         , DIVA%graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( DIVA%N_b                         , DIVA%wN_b                         , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%dN_dx_b                     , DIVA%wdN_dx_b                     , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%dN_dy_b                     , DIVA%wdN_dy_b                     , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%F1_3D_a                     , DIVA%wF1_3D_a                     , DIVA%graphs%graph_a%pai%n_nih, mesh%nz)
    call allocate_dist_shared( DIVA%F2_3D_a                     , DIVA%wF2_3D_a                     , DIVA%graphs%graph_a%pai%n_nih, mesh%nz)
    call allocate_dist_shared( DIVA%F1_3D_b                     , DIVA%wF1_3D_b                     , DIVA%graphs%graph_b%pai%n_nih, mesh%nz)
    call allocate_dist_shared( DIVA%F2_3D_b                     , DIVA%wF2_3D_b                     , DIVA%graphs%graph_b%pai%n_nih, mesh%nz)
    call allocate_dist_shared( DIVA%basal_friction_coefficient_b, DIVA%wbasal_friction_coefficient_b, DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%beta_eff_a                  , DIVA%wbeta_eff_a                  , DIVA%graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( DIVA%beta_eff_b                  , DIVA%wbeta_eff_b                  , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%tau_bx_b                    , DIVA%wtau_bx_b                    , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%tau_by_b                    , DIVA%wtau_by_b                    , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%tau_dx_b                    , DIVA%wtau_dx_b                    , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%tau_dy_b                    , DIVA%wtau_dy_b                    , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%tau_ox_b                    , DIVA%wtau_ox_b                    , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%tau_oy_b                    , DIVA%wtau_oy_b                    , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%u_b_prev                    , DIVA%wu_b_prev                    , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( DIVA%v_b_prev                    , DIVA%wv_b_prev                    , DIVA%graphs%graph_b%pai%n_nih)

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 47)

  end subroutine allocate_DIVA_solver_on_graphs

  subroutine deallocate_DIVA_solver_on_graphs( DIVA)

    ! In/output variables:
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'deallocate_DIVA_solver_on_graphs'

    ! Add routine to path
    call init_routine( routine_name)

    ! Ice model forcing data
    call deallocate_dist_shared( DIVA%Hi_a                        , DIVA%wHi_a                        )
    call deallocate_dist_shared( DIVA%Hs_a                        , DIVA%wHs_a                        )
    call deallocate_dist_shared( DIVA%Ho_a                        , DIVA%wHo_a                        )
    call deallocate_dist_shared( DIVA%A_flow_3D_a                 , DIVA%wA_flow_3D_a                 )
    call deallocate_dist_shared( DIVA%basal_friction_coefficient_a, DIVA%wbasal_friction_coefficient_a)
    call deallocate_dist_shared( DIVA%fraction_gr_b               , DIVA%wfraction_gr_b               )

    ! Solution
    call deallocate_dist_shared( DIVA%u_vav_b                     , DIVA%wu_vav_b                     )
    call deallocate_dist_shared( DIVA%v_vav_b                     , DIVA%wv_vav_b                     )
    call deallocate_dist_shared( DIVA%u_base_b                    , DIVA%wu_base_b                    )
    call deallocate_dist_shared( DIVA%v_base_b                    , DIVA%wv_base_b                    )
    call deallocate_dist_shared( DIVA%u_3D_b                      , DIVA%wu_3D_b                      )
    call deallocate_dist_shared( DIVA%v_3D_b                      , DIVA%wv_3D_b                      )

    ! Intermediate data fields
    call deallocate_dist_shared( DIVA%du_dx_a                     , DIVA%wdu_dx_a                     )
    call deallocate_dist_shared( DIVA%du_dy_a                     , DIVA%wdu_dy_a                     )
    call deallocate_dist_shared( DIVA%dv_dx_a                     , DIVA%wdv_dx_a                     )
    call deallocate_dist_shared( DIVA%dv_dy_a                     , DIVA%wdv_dy_a                     )
    call deallocate_dist_shared( DIVA%du_dx_b                     , DIVA%wdu_dx_b                     )
    call deallocate_dist_shared( DIVA%du_dy_b                     , DIVA%wdu_dy_b                     )
    call deallocate_dist_shared( DIVA%dv_dx_b                     , DIVA%wdv_dx_b                     )
    call deallocate_dist_shared( DIVA%dv_dy_b                     , DIVA%wdv_dy_b                     )
    call deallocate_dist_shared( DIVA%du_dz_3D_a                  , DIVA%wdu_dz_3D_a                  )
    call deallocate_dist_shared( DIVA%dv_dz_3D_a                  , DIVA%wdv_dz_3D_a                  )
    call deallocate_dist_shared( DIVA%du_dz_3D_b                  , DIVA%wdu_dz_3D_b                  )
    call deallocate_dist_shared( DIVA%dv_dz_3D_b                  , DIVA%wdv_dz_3D_b                  )
    call deallocate_dist_shared( DIVA%eta_3D_a                    , DIVA%weta_3D_a                    )
    call deallocate_dist_shared( DIVA%eta_3D_b                    , DIVA%weta_3D_b                    )
    call deallocate_dist_shared( DIVA%eta_vav_a                   , DIVA%weta_vav_a                   )
    call deallocate_dist_shared( DIVA%eta_vav_b                   , DIVA%weta_vav_b                   )
    call deallocate_dist_shared( DIVA%N_a                         , DIVA%wN_a                         )
    call deallocate_dist_shared( DIVA%N_b                         , DIVA%wN_b                         )
    call deallocate_dist_shared( DIVA%dN_dx_b                     , DIVA%wdN_dx_b                     )
    call deallocate_dist_shared( DIVA%dN_dy_b                     , DIVA%wdN_dy_b                     )
    call deallocate_dist_shared( DIVA%F1_3D_a                     , DIVA%wF1_3D_a                     )
    call deallocate_dist_shared( DIVA%F2_3D_a                     , DIVA%wF2_3D_a                     )
    call deallocate_dist_shared( DIVA%F1_3D_b                     , DIVA%wF1_3D_b                     )
    call deallocate_dist_shared( DIVA%F2_3D_b                     , DIVA%wF2_3D_b                     )
    call deallocate_dist_shared( DIVA%basal_friction_coefficient_b, DIVA%wbasal_friction_coefficient_b)
    call deallocate_dist_shared( DIVA%beta_eff_a                  , DIVA%wbeta_eff_a                  )
    call deallocate_dist_shared( DIVA%beta_eff_b                  , DIVA%wbeta_eff_b                  )
    call deallocate_dist_shared( DIVA%tau_bx_b                    , DIVA%wtau_bx_b                    )
    call deallocate_dist_shared( DIVA%tau_by_b                    , DIVA%wtau_by_b                    )
    call deallocate_dist_shared( DIVA%tau_dx_b                    , DIVA%wtau_dx_b                    )
    call deallocate_dist_shared( DIVA%tau_dy_b                    , DIVA%wtau_dy_b                    )
    call deallocate_dist_shared( DIVA%tau_ox_b                    , DIVA%wtau_ox_b                    )
    call deallocate_dist_shared( DIVA%tau_oy_b                    , DIVA%wtau_oy_b                    )
    call deallocate_dist_shared( DIVA%u_b_prev                    , DIVA%wu_b_prev                    )
    call deallocate_dist_shared( DIVA%v_b_prev                    , DIVA%wv_b_prev                    )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_DIVA_solver_on_graphs

  subroutine map_DIVA_from_mesh_to_graphs( mesh, ice, DIVA)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(in   ) :: ice
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_DIVA_from_mesh_to_graphs'

    ! Add routine to path
    call init_routine( routine_name)

    ! Ice model forcing data
    call map_mesh_vertices_to_graph ( mesh, ice%Hi                           , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%Hi_a                        )
    call map_mesh_vertices_to_graph ( mesh, ice%Hs                           , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%Hs_a                        )
    call map_mesh_vertices_to_graph ( mesh, ice%Ho                           , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%Ho_a                        )
    call map_mesh_triangles_to_graph( mesh, ice%fraction_gr_b                , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%fraction_gr_b               )

    ! Solution
    call map_mesh_triangles_to_graph( mesh, DIVA%u_vav_b                     , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%u_vav_b                     )
    call map_mesh_triangles_to_graph( mesh, DIVA%v_vav_b                     , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%v_vav_b                     )
    call map_mesh_triangles_to_graph( mesh, DIVA%u_base_b                    , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%u_base_b                    )
    call map_mesh_triangles_to_graph( mesh, DIVA%v_base_b                    , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%v_base_b                    )
    call map_mesh_triangles_to_graph( mesh, DIVA%u_3D_b                      , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%u_3D_b                      )
    call map_mesh_triangles_to_graph( mesh, DIVA%v_3D_b                      , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%v_3D_b                      )

    ! Intermediate data fields
    call map_mesh_vertices_to_graph ( mesh, DIVA%du_dx_a                     , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%du_dx_a                     )
    call map_mesh_vertices_to_graph ( mesh, DIVA%du_dy_a                     , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%du_dy_a                     )
    call map_mesh_vertices_to_graph ( mesh, DIVA%dv_dx_a                     , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%dv_dx_a                     )
    call map_mesh_vertices_to_graph ( mesh, DIVA%dv_dy_a                     , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%dv_dy_a                     )
    call map_mesh_vertices_to_graph ( mesh, DIVA%du_dz_3D_a                  , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%du_dz_3D_a                  )
    call map_mesh_vertices_to_graph ( mesh, DIVA%dv_dz_3D_a                  , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%dv_dz_3D_a                  )
    call map_mesh_vertices_to_graph ( mesh, DIVA%eta_3D_a                    , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%eta_3D_a                    )
    call map_mesh_triangles_to_graph( mesh, DIVA%eta_3D_b                    , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%eta_3D_b                    )
    call map_mesh_vertices_to_graph ( mesh, DIVA%eta_vav_a                   , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%eta_vav_a                   )
    call map_mesh_vertices_to_graph ( mesh, DIVA%N_a                         , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%N_a                         )
    call map_mesh_triangles_to_graph( mesh, DIVA%N_b                         , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%N_b                         )
    call map_mesh_triangles_to_graph( mesh, DIVA%dN_dx_b                     , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%dN_dx_b                     )
    call map_mesh_triangles_to_graph( mesh, DIVA%dN_dy_b                     , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%dN_dy_b                     )
    call map_mesh_vertices_to_graph ( mesh, DIVA%F1_3D_a                     , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%F1_3D_a                     )
    call map_mesh_vertices_to_graph ( mesh, DIVA%F2_3D_a                     , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%F2_3D_a                     )
    call map_mesh_triangles_to_graph( mesh, DIVA%F1_3D_b                     , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%F1_3D_b                     )
    call map_mesh_triangles_to_graph( mesh, DIVA%F2_3D_b                     , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%F2_3D_b                     )
    call map_mesh_triangles_to_graph( mesh, DIVA%basal_friction_coefficient_b, DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%basal_friction_coefficient_b)
    call map_mesh_vertices_to_graph ( mesh, DIVA%beta_eff_a                  , DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%beta_eff_a                  )
    call map_mesh_triangles_to_graph( mesh, DIVA%beta_eff_b                  , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%beta_eff_b                  )
    call map_mesh_triangles_to_graph( mesh, DIVA%tau_bx_b                    , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%tau_bx_b                    )
    call map_mesh_triangles_to_graph( mesh, DIVA%tau_by_b                    , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%tau_by_b                    )
    call map_mesh_triangles_to_graph( mesh, DIVA%tau_dx_b                    , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%tau_dx_b                    )
    call map_mesh_triangles_to_graph( mesh, DIVA%tau_dy_b                    , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%tau_dy_b                    )
    call map_mesh_triangles_to_graph( mesh, DIVA%u_b_prev                    , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%u_b_prev                    )
    call map_mesh_triangles_to_graph( mesh, DIVA%v_b_prev                    , DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%v_b_prev                    )

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 40)

  end subroutine map_DIVA_from_mesh_to_graphs

  subroutine map_DIVA_from_graphs_to_mesh( mesh, DIVA)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_DIVA_from_graphs_to_mesh'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%u_vav_b , mesh, DIVA%u_vav_b )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%v_vav_b , mesh, DIVA%v_vav_b )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%u_base_b, mesh, DIVA%u_base_b)
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%v_base_b, mesh, DIVA%v_base_b)
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%u_3D_b  , mesh, DIVA%u_3D_b  )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%v_3D_b  , mesh, DIVA%v_3D_b  )

    ! Intermediate data fields
    call map_graph_to_mesh_vertices ( DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%du_dx_a                     , mesh, DIVA%du_dx_a                     )
    call map_graph_to_mesh_vertices ( DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%du_dy_a                     , mesh, DIVA%du_dy_a                     )
    call map_graph_to_mesh_vertices ( DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%dv_dx_a                     , mesh, DIVA%dv_dx_a                     )
    call map_graph_to_mesh_vertices ( DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%dv_dy_a                     , mesh, DIVA%dv_dy_a                     )
    call map_graph_to_mesh_vertices ( DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%du_dz_3D_a                  , mesh, DIVA%du_dz_3D_a                  )
    call map_graph_to_mesh_vertices ( DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%dv_dz_3D_a                  , mesh, DIVA%dv_dz_3D_a                  )
    call map_graph_to_mesh_vertices ( DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%eta_3D_a                    , mesh, DIVA%eta_3D_a                    )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%eta_3D_b                    , mesh, DIVA%eta_3D_b                    )
    call map_graph_to_mesh_vertices ( DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%eta_vav_a                   , mesh, DIVA%eta_vav_a                   )
    call map_graph_to_mesh_vertices ( DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%N_a                         , mesh, DIVA%N_a                         )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%N_b                         , mesh, DIVA%N_b                         )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%dN_dx_b                     , mesh, DIVA%dN_dx_b                     )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%dN_dy_b                     , mesh, DIVA%dN_dy_b                     )
    call map_graph_to_mesh_vertices ( DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%F1_3D_a                     , mesh, DIVA%F1_3D_a                     )
    call map_graph_to_mesh_vertices ( DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%F2_3D_a                     , mesh, DIVA%F2_3D_a                     )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%F1_3D_b                     , mesh, DIVA%F1_3D_b                     )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%F2_3D_b                     , mesh, DIVA%F2_3D_b                     )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%basal_friction_coefficient_b, mesh, DIVA%basal_friction_coefficient_b)
    call map_graph_to_mesh_vertices ( DIVA%DIVA_graphs%graphs%graph_a, DIVA%DIVA_graphs%beta_eff_a                  , mesh, DIVA%beta_eff_a                  )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%beta_eff_b                  , mesh, DIVA%beta_eff_b                  )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%tau_bx_b                    , mesh, DIVA%tau_bx_b                    )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%tau_by_b                    , mesh, DIVA%tau_by_b                    )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%tau_dx_b                    , mesh, DIVA%tau_dx_b                    )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%tau_dy_b                    , mesh, DIVA%tau_dy_b                    )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%u_b_prev                    , mesh, DIVA%u_b_prev                    )
    call map_graph_to_mesh_triangles( DIVA%DIVA_graphs%graphs%graph_b, DIVA%DIVA_graphs%v_b_prev                    , mesh, DIVA%v_b_prev                    )

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 40)

  end subroutine map_DIVA_from_graphs_to_mesh

  ! == Calculate several intermediate terms in the DIVA

  subroutine calc_driving_stress( DIVA)

    ! In/output variables:
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'calc_driving_stress'
    real(dp), dimension(:), pointer :: Hi_b     => null()
    real(dp), dimension(:), pointer :: dHs_dx_b => null()
    real(dp), dimension(:), pointer :: dHs_dy_b => null()
    type(MPI_WIN)                   :: wHi_b, wdHs_dx_b, wdHs_dy_b
    integer                         :: ni

    real(dp), dimension(:), pointer :: d_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( Hi_b    , wHi_b    , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( dHs_dx_b, wdHs_dx_b, DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( dHs_dy_b, wdHs_dy_b, DIVA%graphs%graph_b%pai%n_nih)

    call exchange_halos( DIVA%graphs%graph_a, DIVA%Hi_a)
    call exchange_halos( DIVA%graphs%graph_a, DIVA%Hs_a)

    ! Calculate ice thickness and surface slopes on the b-graph
    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_map_a_b, &
      DIVA%graphs%graph_a%pai, DIVA%Hi_a, DIVA%graphs%graph_b%pai, Hi_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_a%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_ddx_a_b, &
      DIVA%graphs%graph_a%pai, DIVA%Hs_a, DIVA%graphs%graph_b%pai, dHs_dx_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_a%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_ddy_a_b, &
      DIVA%graphs%graph_a%pai, DIVA%Hs_a, DIVA%graphs%graph_b%pai, dHs_dy_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_a%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    ! Calculate the driving stress
    do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2
      DIVA%tau_dx_b( ni) = -ice_density * grav * Hi_b( ni) * dHs_dx_b( ni)
      DIVA%tau_dy_b( ni) = -ice_density * grav * Hi_b( ni) * dHs_dy_b( ni)
    end do

    ! DENK DROM
    d_loc => DIVA%tau_dx_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'tau_dx_b')
    d_loc => DIVA%tau_dy_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'tau_dy_b')

    ! Clean up after yourself
    call deallocate_dist_shared( Hi_b    , wHi_b    )
    call deallocate_dist_shared( dHs_dx_b, wdHs_dx_b)
    call deallocate_dist_shared( dHs_dy_b, wdHs_dy_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_driving_stress

  subroutine calc_ocean_back_pressure( DIVA)

    ! In/output variables:
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'calc_ocean_back_pressure'
    real(dp), dimension(:), pointer :: Hi_b => null()
    real(dp), dimension(:), pointer :: Ho_b => null()
    type(MPI_WIN)                   :: wHi_b, wHo_b
    integer                         :: ni

    real(dp), dimension(:), pointer :: d_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( Hi_b, wHi_b, DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( Ho_b, wHo_b, DIVA%graphs%graph_b%pai%n_nih)

    call exchange_halos( DIVA%graphs%graph_a, DIVA%Hi_a)
    call exchange_halos( DIVA%graphs%graph_a, DIVA%Ho_a)

    ! Calculate ice thickness and adjacent ocean column height on the b-graph
    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_map_a_b, &
      DIVA%graphs%graph_a%pai, DIVA%Hi_a, DIVA%graphs%graph_b%pai, Hi_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_a%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_map_a_b, &
      DIVA%graphs%graph_a%pai, DIVA%Ho_a, DIVA%graphs%graph_b%pai, Ho_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_a%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    ! Calculate the ocean back pressure
    do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2
      DIVA%tau_ox_b( ni) = (&
          0.5_dp * ice_density      * grav * Hi_b( ni)**2 &
        - 0.5_dp * seawater_density * grav * Ho_b( ni)**2) * DIVA%graphs%graph_b%border_nhat( ni,1)
      DIVA%tau_oy_b( ni) = (&
          0.5_dp * ice_density      * grav * Hi_b( ni)**2 &
        - 0.5_dp * seawater_density * grav * Ho_b( ni)**2) * DIVA%graphs%graph_b%border_nhat( ni,2)
    end do

    ! DENK DROM
    d_loc => DIVA%tau_ox_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'tau_ox_b')
    d_loc => DIVA%tau_oy_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'tau_oy_b')

    ! Clean up after yourself
    call deallocate_dist_shared( Hi_b, wHi_b)
    call deallocate_dist_shared( Ho_b, wHo_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ocean_back_pressure

  subroutine calc_horizontal_strain_rates( DIVA)
    !< Calculate the vertically averaged horizontal strain rates

    ! In/output variables:
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_horizontal_strain_rates'

    real(dp), dimension(:), pointer :: d_loc

    ! Add routine to path
    call init_routine( routine_name)

    call exchange_halos( DIVA%graphs%graph_b, DIVA%u_vav_b)
    call exchange_halos( DIVA%graphs%graph_b, DIVA%v_vav_b)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_ddx_b_a, &
      DIVA%graphs%graph_b%pai, DIVA%u_vav_b, DIVA%graphs%graph_a%pai, DIVA%du_dx_a, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_b%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_a%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_ddy_b_a, &
      DIVA%graphs%graph_b%pai, DIVA%u_vav_b, DIVA%graphs%graph_a%pai, DIVA%du_dy_a, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_b%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_a%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_ddx_b_a, &
      DIVA%graphs%graph_b%pai, DIVA%v_vav_b, DIVA%graphs%graph_a%pai, DIVA%dv_dx_a, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_b%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_a%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_ddy_b_a, &
      DIVA%graphs%graph_b%pai, DIVA%v_vav_b, DIVA%graphs%graph_a%pai, DIVA%dv_dy_a, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_b%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_a%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_ddx_b_b, &
      DIVA%graphs%graph_b%pai, DIVA%u_vav_b, DIVA%graphs%graph_b%pai, DIVA%du_dx_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_b%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_ddy_b_b, &
      DIVA%graphs%graph_b%pai, DIVA%u_vav_b, DIVA%graphs%graph_b%pai, DIVA%du_dy_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_b%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_ddx_b_b, &
      DIVA%graphs%graph_b%pai, DIVA%v_vav_b, DIVA%graphs%graph_b%pai, DIVA%dv_dx_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_b%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_ddy_b_b, &
      DIVA%graphs%graph_b%pai, DIVA%v_vav_b, DIVA%graphs%graph_b%pai, DIVA%dv_dy_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_b%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    ! DENK DROM
    d_loc => DIVA%du_dx_a( DIVA%graphs%graph_a%ni1: DIVA%graphs%graph_a%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'du_dx_a')
    d_loc => DIVA%du_dy_a( DIVA%graphs%graph_a%ni1: DIVA%graphs%graph_a%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'du_dy_a')
    d_loc => DIVA%dv_dx_a( DIVA%graphs%graph_a%ni1: DIVA%graphs%graph_a%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'dv_dx_a')
    d_loc => DIVA%dv_dy_a( DIVA%graphs%graph_a%ni1: DIVA%graphs%graph_a%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'dv_dy_a')
    d_loc => DIVA%du_dx_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1d( trim( C%output_dir), d_loc, 'du_dx_b')
    d_loc => DIVA%du_dy_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1d( trim( C%output_dir), d_loc, 'du_dy_b')
    d_loc => DIVA%dv_dx_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1d( trim( C%output_dir), d_loc, 'dv_dx_b')
    d_loc => DIVA%dv_dy_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1d( trim( C%output_dir), d_loc, 'dv_dy_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_horizontal_strain_rates

  subroutine calc_vertical_shear_strain_rates( mesh, DIVA)
    ! Calculate the vertical shear strain rates

    ! In/output variables:
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'calc_vertical_shear_strain_rates'
    integer                               :: ni, k

    real(dp), dimension(:,:), pointer :: d_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate (parameterised) vertical shear strain rates on the b-grid (Lipscomb et al., 2019, Eq. 36)
    do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2
    do k = 1, mesh%nz
      DIVA%du_dz_3D_b( ni,k) = DIVA%tau_bx_b( ni) * mesh%zeta( k) / max( C%visc_eff_min, DIVA%eta_3D_b( ni,k))
      DIVA%dv_dz_3D_b( ni,k) = DIVA%tau_by_b( ni) * mesh%zeta( k) / max( C%visc_eff_min, DIVA%eta_3D_b( ni,k))
    end do
    end do

    call exchange_halos( DIVA%graphs%graph_b, DIVA%du_dz_3D_b)
    call exchange_halos( DIVA%graphs%graph_b, DIVA%dv_dz_3D_b)

    call multiply_CSR_matrix_with_vector_2D_wrapper( DIVA%graphs%M_map_b_a, &
      DIVA%graphs%graph_b%pai, DIVA%du_dz_3D_b, DIVA%graphs%graph_a%pai, DIVA%du_dz_3D_a, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_b%buffer1_gk_nih, &
      buffer_yy_nih = DIVA%graphs%graph_a%buffer2_gk_nih)

    call multiply_CSR_matrix_with_vector_2D_wrapper( DIVA%graphs%M_map_b_a, &
      DIVA%graphs%graph_b%pai, DIVA%dv_dz_3D_b, DIVA%graphs%graph_a%pai, DIVA%dv_dz_3D_a, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_b%buffer1_gk_nih, &
      buffer_yy_nih = DIVA%graphs%graph_a%buffer2_gk_nih)

    ! DENK DROM
    d_loc => DIVA%du_dz_3D_a( DIVA%graphs%graph_a%ni1: DIVA%graphs%graph_a%ni2, 1:mesh%nz)
    call save_variable_as_netcdf_dp_2D( trim( C%output_dir), d_loc, 'du_dz_3D_a')
    d_loc => DIVA%dv_dz_3D_a( DIVA%graphs%graph_a%ni1: DIVA%graphs%graph_a%ni2, 1:mesh%nz)
    call save_variable_as_netcdf_dp_2D( trim( C%output_dir), d_loc, 'dv_dz_3D_a')
    d_loc => DIVA%du_dz_3D_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2, 1:mesh%nz)
    call save_variable_as_netcdf_dp_2d( trim( C%output_dir), d_loc, 'du_dz_3D_b')
    d_loc => DIVA%dv_dz_3D_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2, 1:mesh%nz)
    call save_variable_as_netcdf_dp_2d( trim( C%output_dir), d_loc, 'dv_dz_3D_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_vertical_shear_strain_rates

  subroutine calc_effective_viscosity( mesh, ice, DIVA, Glens_flow_law_epsilon_sq_0_applied)

    ! In/output variables:
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_ice_model),                       intent(inout) :: ice
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA
    real(dp),                                   intent(in   ) :: Glens_flow_law_epsilon_sq_0_applied

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'calc_effective_viscosity'
    integer                           :: ni,k
    real(dp)                          :: A_min, eta_max
    real(dp), dimension( mesh%nz)     :: prof
    real(dp), dimension(:  ), pointer :: Hi_b => null()
    real(dp), dimension(:,:), pointer :: A_flow_3D_b => null()
    type(MPI_WIN)                     :: wHi_b, wA_flow_3D_b

    real(dp), dimension(:), pointer :: d_loc
    real(dp), dimension(:,:), pointer :: d_lock

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( Hi_b       , wHi_b       , DIVA%graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( A_flow_3D_b, wA_flow_3D_b, DIVA%graphs%graph_b%pai%n_nih, mesh%nz)

    ! Calculate maximum allowed effective viscosity, for stability
    A_min = 1E-18_dp
    eta_max = 0.5_dp * A_min**(-1._dp / C%Glens_flow_law_exponent) * (Glens_flow_law_epsilon_sq_0_applied)**((1._dp - C%Glens_flow_law_exponent)/(2._dp*C%Glens_flow_law_exponent))

    ! == On the a-graph ==
    ! ====================

    ! Calculate the effective viscosity eta
    if (C%choice_flow_law == 'Glen') then
      ! Calculate the effective viscosity eta according to Glen's flow law

      ! Calculate flow factors
      call calc_ice_rheology_Glen( mesh, ice)

      call map_mesh_vertices_to_graph( mesh, ice%A_flow, DIVA%graphs%graph_a, DIVA%A_flow_3D_a)

      ! Calculate effective viscosity
      do ni = DIVA%graphs%graph_a%ni1, DIVA%graphs%graph_a%ni2
      do k  = 1, mesh%nz
        DIVA%eta_3D_a( ni,k) = calc_effective_viscosity_Glen_3D_uv_only( &
          Glens_flow_law_epsilon_sq_0_applied, &
          DIVA%du_dx_a( ni), DIVA%du_dy_a( ni), DIVA%du_dz_3D_a( ni,k), &
          DIVA%dv_dx_a( ni), DIVA%dv_dy_a( ni), DIVA%dv_dz_3D_a( ni,k), DIVA%A_flow_3D_a( ni,k))

        ! Safety
        DIVA%eta_3D_a( ni,k) = min( max( DIVA%eta_3D_a( ni,k), C%visc_eff_min), eta_max)
      end do
      end do

    else
      call crash('unknown choice_flow_law "' // TRIM( C%choice_flow_law) // '"!')
    end if

    ! Calculate vertically averaged effective viscosity
    do ni = DIVA%graphs%graph_a%ni1, DIVA%graphs%graph_a%ni2
      prof = DIVA%eta_3D_a( ni,:)
      DIVA%eta_vav_a( ni) = vertical_average( mesh%zeta, prof)
    end do

    ! Calculate the product term N = eta_vav * H
    do ni = DIVA%graphs%graph_a%ni1, DIVA%graphs%graph_a%ni2
      DIVA%N_a( ni) = DIVA%eta_vav_a( ni) * max( 0.1, DIVA%Hi_a( ni))
    end do

    ! DENK DROM
    d_lock => DIVA%eta_3D_a( DIVA%graphs%graph_a%ni1: DIVA%graphs%graph_a%ni2, 1:mesh%nz)
    call save_variable_as_netcdf_dp_2D( trim( C%output_dir), d_lock, 'eta_3D_a')
    d_loc => DIVA%eta_vav_a( DIVA%graphs%graph_a%ni1: DIVA%graphs%graph_a%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'eta_vav_a')
    d_loc => DIVA%N_a( DIVA%graphs%graph_a%ni1: DIVA%graphs%graph_a%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'N_a')

    ! == On the b-graph ==
    ! ====================

    ! Calculate the effective viscosity eta
    if (C%choice_flow_law == 'Glen') then
      ! Calculate the effective viscosity eta according to Glen's flow law

      ! Map flow factor from the a-graph to the b-graph
      call multiply_CSR_matrix_with_vector_2D_wrapper( DIVA%graphs%M_map_a_b, &
        DIVA%graphs%graph_a%pai, DIVA%A_flow_3D_a, DIVA%graphs%graph_b%pai, A_flow_3D_b, &
        xx_is_hybrid = .true., yy_is_hybrid = .true., &
        buffer_xx_nih = DIVA%graphs%graph_a%buffer1_gk_nih, &
        buffer_yy_nih = DIVA%graphs%graph_b%buffer2_gk_nih)

      ! Calculate effective viscosity
      do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2
      do k  = 1, mesh%nz
        DIVA%eta_3D_b( ni,k) = calc_effective_viscosity_Glen_3D_uv_only( &
          Glens_flow_law_epsilon_sq_0_applied, &
          DIVA%du_dx_b( ni), DIVA%du_dy_b( ni), DIVA%du_dz_3D_b( ni,k), &
          DIVA%dv_dx_b( ni), DIVA%dv_dy_b( ni), DIVA%dv_dz_3D_b( ni,k), A_flow_3D_b( ni,k))

        ! Safety
        DIVA%eta_3D_b( ni,k) = min( max( DIVA%eta_3D_b( ni,k), C%visc_eff_min), eta_max)
      end do
      end do

    else
      call crash('unknown choice_flow_law "' // TRIM( C%choice_flow_law) // '"!')
    end if

    ! Calculate vertically averaged effective viscosity
    do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2
      prof = DIVA%eta_3D_b( ni,:)
      DIVA%eta_vav_b( ni) = vertical_average( mesh%zeta, prof)
    end do

    ! Map ice thickness from the a-graph to the b-graph
    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_map_a_b, &
      DIVA%graphs%graph_a%pai, DIVA%Hi_a, DIVA%graphs%graph_b%pai, Hi_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_a%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    ! Calculate the product term N = eta_vav * H
    do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2
      DIVA%N_b( ni) = DIVA%eta_vav_b( ni) * max( 0.1, Hi_b( ni))
    end do

    ! DENK DROM
    d_lock => DIVA%eta_3D_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2, 1:mesh%nz)
    call save_variable_as_netcdf_dp_2d( trim( C%output_dir), d_lock, 'eta_3D_b')
    d_loc => DIVA%eta_vav_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1d( trim( C%output_dir), d_loc, 'eta_vav_b')
    d_loc => DIVA%N_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1d( trim( C%output_dir), d_loc, 'N_b')

    ! == Calculate gradients of N on the b-graph ==
    ! =============================================

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_ddx_a_b, &
      DIVA%graphs%graph_a%pai, DIVA%N_a, DIVA%graphs%graph_b%pai, DIVA%dN_dx_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_a%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_ddy_a_b, &
      DIVA%graphs%graph_a%pai, DIVA%N_a, DIVA%graphs%graph_b%pai, DIVA%dN_dy_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_a%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    ! DENK DROM
    d_loc => DIVA%dN_dx_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'dN_dx_b')
    d_loc => DIVA%dN_dy_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'dN_dy_b')

    ! Clean up after yourself
    call deallocate_dist_shared( Hi_b       , wHi_b       )
    call deallocate_dist_shared( A_flow_3D_b, wA_flow_3D_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_viscosity

  subroutine calc_F_integrals( mesh, DIVA)
    !< Calculate the F-integrals on the a-grid (Lipscomb et al. (2019), Eq. 30)

    ! In/output variables:
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'calc_F_integrals'
    integer                           :: ni,k
    real(dp), dimension( mesh%nz)     :: prof

    real(dp), dimension(:,:), pointer :: d_loc

    ! Add routine to path
    call init_routine( routine_name)

    do ni = DIVA%graphs%graph_a%ni1, DIVA%graphs%graph_a%ni2

      ! F1
      do k = 1, mesh%nz
        prof( k) = (mesh%zeta( k)    / DIVA%eta_3D_a( ni,k))
      end do
      DIVA%F1_3D_a( ni,:) = -max( 0.1_dp, DIVA%Hi_a( ni)) * &
        integrate_from_zeta_is_one_to_zeta_is_zetap( mesh%zeta, prof)

      ! F2
      do k = 1, mesh%nz
        prof( k) = (mesh%zeta( k)**2 / DIVA%eta_3D_a( ni,k))
      end do
      DIVA%F2_3D_a( ni,:) = -max( 0.1_dp, DIVA%Hi_a( ni)) * &
        integrate_from_zeta_is_one_to_zeta_is_zetap( mesh%zeta, prof)

    end do

    ! Map F-integrals from the a-grid to the b-grid

    call exchange_halos( DIVA%graphs%graph_a, DIVA%F1_3D_a)
    call exchange_halos( DIVA%graphs%graph_a, DIVA%F2_3D_a)

    call multiply_CSR_matrix_with_vector_2D_wrapper( DIVA%graphs%M_map_a_b, &
      DIVA%graphs%graph_a%pai, DIVA%F1_3D_a, DIVA%graphs%graph_b%pai, DIVA%F1_3D_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_a%buffer1_gk_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_gk_nih)

    call multiply_CSR_matrix_with_vector_2D_wrapper( DIVA%graphs%M_map_a_b, &
      DIVA%graphs%graph_a%pai, DIVA%F2_3D_a, DIVA%graphs%graph_b%pai, DIVA%F2_3D_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_a%buffer1_gk_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_gk_nih)

      ! DENK DROM
    d_loc => DIVA%F1_3D_a( DIVA%graphs%graph_a%ni1: DIVA%graphs%graph_a%ni2, 1:mesh%nz)
    call save_variable_as_netcdf_dp_2D( trim( C%output_dir), d_loc, 'F1_3D_a')
    d_loc => DIVA%F2_3D_a( DIVA%graphs%graph_a%ni1: DIVA%graphs%graph_a%ni2, 1:mesh%nz)
    call save_variable_as_netcdf_dp_2D( trim( C%output_dir), d_loc, 'F2_3D_a')
    d_loc => DIVA%F1_3D_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2, 1:mesh%nz)
    call save_variable_as_netcdf_dp_2D( trim( C%output_dir), d_loc, 'F1_3D_b')
    d_loc => DIVA%F2_3D_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2, 1:mesh%nz)
    call save_variable_as_netcdf_dp_2D( trim( C%output_dir), d_loc, 'F2_3D_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_F_integrals

  subroutine calc_effective_basal_friction_coefficient( mesh, ice, bed_roughness, DIVA)
    !< Calculate the "effective" friction coefficient (turning the SSA into the DIVA)

    ! In/output variables:
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_ice_model),                       intent(inout) :: ice
    type(type_bed_roughness_model),             intent(in   ) :: bed_roughness
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_effective_basal_friction_coefficient'
    real(dp), dimension(:), pointer        :: u_base_a => null()
    real(dp), dimension(:), pointer        :: v_base_a => null()
    type(MPI_WIN)                          :: wu_base_a, wv_base_a
    real(dp), dimension(mesh%vi1:mesh%vi2) :: u_base_a_m, v_base_a_m
    integer                                :: ni

    real(dp), dimension(:), pointer :: d_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( u_base_a, wu_base_a, DIVA%graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( v_base_a, wv_base_a, DIVA%graphs%graph_a%pai%n_nih)

    ! Map basal velocity to the a-graph

    call exchange_halos( DIVA%graphs%graph_b, DIVA%u_base_b)
    call exchange_halos( DIVA%graphs%graph_b, DIVA%v_base_b)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_map_b_a, &
      DIVA%graphs%graph_b%pai, DIVA%u_base_b, &
      DIVA%graphs%graph_a%pai, u_base_a, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_b%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_a%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_map_b_a, &
      DIVA%graphs%graph_b%pai, DIVA%v_base_b, &
      DIVA%graphs%graph_a%pai, v_base_a, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_b%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_a%buffer2_g_nih)

    ! Calculate the basal friction coefficient beta_b for the current velocity solution
    ! This is where the sliding law is called!
    call map_graph_to_mesh_vertices( DIVA%graphs%graph_a, u_base_a, mesh, u_base_a_m)
    call map_graph_to_mesh_vertices( DIVA%graphs%graph_a, v_base_a, mesh, v_base_a_m)
    call calc_basal_friction_coefficient( mesh, ice, bed_roughness, u_base_a_m, v_base_a_m)
    call map_mesh_vertices_to_graph( mesh, ice%basal_friction_coefficient, DIVA%graphs%graph_a, DIVA%basal_friction_coefficient_a)

    ! DENK DROM
    d_loc => DIVA%basal_friction_coefficient_a( DIVA%graphs%graph_a%ni1: DIVA%graphs%graph_a%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'basal_friction_coefficient_a')

    ! Calculate beta_eff on the a-grid
    if (C%choice_sliding_law == 'no_sliding') then
      ! Exception for the case of no sliding (Lipscomb et al., 2019, Eq. 35)

      do ni = DIVA%graphs%graph_a%ni1, DIVA%graphs%graph_a%ni2
        DIVA%beta_eff_a( ni) = 1._dp / DIVA%F2_3D_a( ni,1)
      end do

    else
      ! Lipscomb et al., 2019, Eq. 33

      do ni = DIVA%graphs%graph_a%ni1, DIVA%graphs%graph_a%ni2
        DIVA%beta_eff_a( ni) = DIVA%basal_friction_coefficient_a( ni) / (1._dp + DIVA%basal_friction_coefficient_a( ni) * DIVA%F2_3D_a( ni,1))
      end do

    end if

    ! DENK DROM
    d_loc => DIVA%beta_eff_a( DIVA%graphs%graph_a%ni1: DIVA%graphs%graph_a%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'beta_eff_a')

    ! Map basal friction coefficient beta_b and effective basal friction coefficient beta_eff to the b-grid

    call exchange_halos( DIVA%graphs%graph_a, DIVA%basal_friction_coefficient_a)
    call exchange_halos( DIVA%graphs%graph_a, DIVA%beta_eff_a)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_map_a_b, &
      DIVA%graphs%graph_a%pai, DIVA%basal_friction_coefficient_a, &
      DIVA%graphs%graph_b%pai, DIVA%basal_friction_coefficient_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_a%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( DIVA%graphs%M_map_a_b, &
      DIVA%graphs%graph_a%pai, DIVA%beta_eff_a, DIVA%graphs%graph_b%pai, DIVA%beta_eff_b, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = DIVA%graphs%graph_a%buffer1_g_nih, &
      buffer_yy_nih = DIVA%graphs%graph_b%buffer2_g_nih)

    ! Apply the sub-grid grounded fraction, and limit the friction coefficient to improve stability
    if (C%do_GL_subgrid_friction) then
      ! On the b-grid
      do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2
        DIVA%beta_eff_b( ni) = DIVA%beta_eff_b( ni) * DIVA%fraction_gr_b( ni)**C%subgrid_friction_exponent_on_B_grid
      end do
    end if

    ! DENK DROM
    d_loc => DIVA%basal_friction_coefficient_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'basal_friction_coefficient_b')
    d_loc => DIVA%beta_eff_b( DIVA%graphs%graph_b%ni1: DIVA%graphs%graph_b%ni2)
    call save_variable_as_netcdf_dp_1D( trim( C%output_dir), d_loc, 'beta_eff_b')

    ! Clean up after yourself
    call deallocate_dist_shared( u_base_a, wu_base_a)
    call deallocate_dist_shared( v_base_a, wv_base_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_basal_friction_coefficient

  subroutine apply_velocity_limits( DIVA)
    !< Limit velocities for improved stability

    ! In/output variables:
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'apply_velocity_limits'
    integer                        :: ni
    real(dp)                       :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2

      ! Calculate absolute speed
      uabs = sqrt( DIVA%u_vav_b( ni)**2 + DIVA%v_vav_b( ni)**2)

      ! Reduce velocities if necessary
      if (uabs > C%vel_max) then
        DIVA%u_vav_b( ni) = DIVA%u_vav_b( ni) * C%vel_max / uabs
        DIVA%v_vav_b( ni) = DIVA%v_vav_b( ni) * C%vel_max / uabs
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_velocity_limits

  subroutine relax_viscosity_iterations( DIVA, visc_it_relax)
    !< Reduce the change between velocity solutions

    ! In/output variables:
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA
    real(dp),                                   intent(in   ) :: visc_it_relax

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'relax_viscosity_iterations'
    integer                        :: ni

    ! Add routine to path
    call init_routine( routine_name)

    do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2
      DIVA%u_vav_b( ni) = (visc_it_relax * DIVA%u_vav_b( ni)) + ((1._dp - visc_it_relax) * DIVA%u_b_prev( ni))
      DIVA%v_vav_b( ni) = (visc_it_relax * DIVA%v_vav_b( ni)) + ((1._dp - visc_it_relax) * DIVA%v_b_prev( ni))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine relax_viscosity_iterations

  subroutine calc_basal_velocities( DIVA)
    !< Calculate basal velocities (Lipscomb et al., 2019, Eq. 32)

    ! In/output variables:
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_basal_velocities'
    integer                        :: ni

    ! Add routine to path
    call init_routine( routine_name)

    if (C%choice_sliding_law == 'no_sliding') then
      ! Exception for the case of no sliding

      DIVA%u_base_b = 0._dp
      DIVA%v_base_b = 0._dp

    else

      ! Calculate basal velocities (Lipscomb et al., 2019, Eq. 32)
      do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2
        DIVA%u_base_b( ni) = DIVA%u_vav_b( ni) / (1._dp + DIVA%basal_friction_coefficient_b( ni) * DIVA%F2_3D_b( ni,1))
        DIVA%v_base_b( ni) = DIVA%v_vav_b( ni) / (1._dp + DIVA%basal_friction_coefficient_b( ni) * DIVA%F2_3D_b( ni,1))
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_basal_velocities

  subroutine calc_basal_shear_stress( DIVA)
    !< Calculate the basal shear stress (Lipscomb et al., 2019, just above Eq. 33)

    ! In/output variables:
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_basal_shear_stress'
    integer                        :: ni

    ! Add routine to path
    call init_routine( routine_name)

    do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2
      ! Lipscomb et al., 2019, just above Eq. 33
      DIVA%tau_bx_b( ni) = DIVA%u_vav_b( ni) * DIVA%beta_eff_b( ni)
      DIVA%tau_by_b( ni) = DIVA%v_vav_b( ni) * DIVA%beta_eff_b( ni)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_basal_shear_stress

  subroutine calc_L2_norm_uv( DIVA, L2_uv)
    !< Calculate the L2-norm of the two consecutive velocity solutions

    ! In/output variables:
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA
    real(dp),                                   intent(  out) :: L2_uv

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_visc_iter_UV_resid'
    integer                        :: ierr
    integer                        :: ni
    real(dp)                       :: res1, res2

    ! Add routine to path
    call init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2

      if (.not. isnan( DIVA%u_vav_b( ni)) .and. .not. isnan( DIVA%v_vav_b( ni))) then

        res1 = res1 + (DIVA%u_vav_b( ni) - DIVA%u_b_prev( ni))**2
        res1 = res1 + (DIVA%v_vav_b( ni) - DIVA%v_b_prev( ni))**2

        res2 = res2 + (DIVA%u_vav_b( ni) + DIVA%u_b_prev( ni))**2
        res2 = res2 + (DIVA%v_vav_b( ni) + DIVA%v_b_prev( ni))**2

      end if

    end do

    ! Combine results from all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate L2-norm
    L2_uv = 2._dp * res1 / max( res2, 1E-8_dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_L2_norm_uv

  subroutine calc_3D_velocities( mesh, DIVA)
    !< Calculate 3D velocities (Lipscomb et al., 2019, Eq. 29)

    ! In/output variables:
    type(type_mesh),                            intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_3D_velocities'
    integer                        :: ni,k

    ! Add routine to path
    call init_routine( routine_name)

    if (C%choice_sliding_law == 'no_sliding') then
      ! Exception for the case of no sliding

      do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2
      do k = 1, mesh%nz
        ! Lipscomb et al., 2019, Eq. 29, and text between Eqs. 33 and 34
        DIVA%u_3D_b( ni,k) = DIVA%tau_bx_b( ni) * DIVA%F1_3D_b( ni,k)
        DIVA%v_3D_b( ni,k) = DIVA%tau_by_b( ni) * DIVA%F1_3D_b( ni,k)
      end do
      end do

    else

      do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2
      do k = 1, mesh%nz
        ! Lipscomb et al., 2019, Eq. 29
        DIVA%u_3D_b( ni,k) = DIVA%u_base_b( ni) * (1._dp + DIVA%basal_friction_coefficient_b( ni) * DIVA%F1_3D_b( ni,k))
        DIVA%v_3D_b( ni,k) = DIVA%v_base_b( ni) * (1._dp + DIVA%basal_friction_coefficient_b( ni) * DIVA%F1_3D_b( ni,k))
      end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_velocities

end module DIVA_solver_ocean_pressure
