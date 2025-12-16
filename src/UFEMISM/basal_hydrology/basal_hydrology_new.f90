MODULE basal_hydrology_new

! ===== Preamble =====
! ====================
! This is preamble is still copied from thermodynamics_3D_heat_equation.f90

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE basal_hydrology_model_types                            , ONLY: type_basal_hydrology_model
  USE climate_model_types                                    , ONLY: type_climate_model
  USE SMB_model_types                                        , ONLY: type_SMB_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  use zeta_gradients, only: calc_zeta_gradients
  USE thermodynamics_utilities                               , ONLY: calc_heat_capacity, calc_thermal_conductivity, calc_pressure_melting_point, &
                                                                     calc_upwind_heat_flux_derivatives, calc_strain_heating, calc_frictional_heating, &
                                                                     replace_Ti_with_robin_solution
  use tridiagonal_solver, only: solve_tridiagonal_matrix_equation
  use netcdf_io_main
  USE mesh_disc_apply_operators                              , ONLY: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, map_b_a_2D, ddx_a_a_2D, ddy_a_a_2D
  use laddie_utilities                                       , ONLY: map_H_a_c
  use mpi_distributed_memory                                 , only: gather_to_all
  use mesh_halo_exchange                                     , only: exchange_halos
  use CSR_matrix_vector_multiplication                       , only: multiply_CSR_matrix_with_vector_1D
  use mesh_utilities                                         , only: find_containing_vertex
  use CSR_matrix_basics                                      , only: finalise_matrix_CSR_dist, add_entry_CSR_dist, add_empty_row_CSR_dist, allocate_matrix_CSR_dist, deallocate_matrix_CSR_dist

  IMPLICIT NONE

CONTAINS

  subroutine basal_hydrology(mesh, ice, basal_hydro, time)
    ! In/output variables:
    type(type_mesh),                  intent(in   ) :: mesh
    type(type_ice_model),             intent(inout) :: ice
    type(type_basal_hydrology_model), intent(inout) :: basal_hydro
    real(dp),                         intent(in)    :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'basal_hydrology'
    integer                        :: vi
    logical, parameter             :: sliding = .true.         ! Whether the ice slides or not when W = 0 and there is ice on land

    real(dp), parameter            :: W_max = 1000000.0_dp     ! Maximum basal water depth (number is placeholder for now, because there is no maximum W)
    real(dp), parameter            :: W_min = 0.0_dp           ! Minimum basal water depth

    real(dp), parameter            :: W_max_til = 2.0_dp       ! Maximum basal water depth in till
    real(dp), parameter            :: W_min_til = 0.0_dp       ! Minimum basal water depth in till

    real(dp), parameter            :: P_min = 0.0_dp           ! Minimum pressure of the ice on the basal water

    real(dp), parameter            :: rho_w = 1000.0_dp   ! Density of water
    real(dp), parameter            :: rho_i = 917.0_dp    ! Density of ice

    ! real(dp), parameter            :: k_coef              ! Coefficient of effective conductivity
    ! real(dp), parameter            :: alpha               ! Exponent used in effective conductivity
    ! real(dp), parameter            :: beta                ! Exponent used in effective conductivity

    real(dp), parameter            :: Cd = 0.0_dp        ! Gradual drain of water in till (m a^-1)

    real(dp), parameter            :: g = 9.81_dp         ! Gravitational acceleration

    ! real(dp), parameter            :: c_1                 ! Scaling coefficient 1 (non-negative opening term)
    ! real(dp), parameter            :: c_2                 ! Scaling coefficient 2 (closing term)

    ! real(dp), parameter            :: W_r                 ! Maximum roughness scale of basal topography

    integer                        :: ti, via, vib, vic   ! Triangle index
    real(dp)                       :: d_ab, d_bc, d_ca    ! Lengths of triangle sides
    real(dp)                       :: d_min               ! Minimum triangle side length
    real(dp)                       :: u_t, v_t, D_t       ! Velocity components on triangle
    real(dp)                       :: ci, vj, ei, u_perp

    real(dp)                       :: dt_crit_CFL, dt_crit_W, dt_crit_P, dt_hydro   ! Critical timesteps for different conditions

    real(dp), parameter            :: phi = 0.01_dp       ! Englacial porosity (This is just a parameter)

    real(dp), parameter            :: correction_factor = 0.9_dp ! To be on the safe side

    logical,  dimension(mesh%nV)   :: mask_floating_ice_tot, mask_icefree_land_tot, mask_icefree_ocean_tot
    real(dp), dimension(mesh%nE)   :: u_c_tot, v_c_tot
    real(dp), dimension(mesh%nV)   :: W_tot

    real(dp)                       :: dt

    ! Add routine to path
    call init_routine( routine_name)

    ! Get the general timestep
    call calc_general_dt(time, basal_hydro, dt)

    ! Initialise basal hydro masks
    call calc_basal_hydro_mask_a_b(mesh, ice, basal_hydro)

    ! Point source
    call point_source( mesh, basal_hydro)

    ! 1) Start with W, W_til and P and make sure they are all within their bounds
    call set_within_bounds(mesh, ice, basal_hydro, W_min, W_max, W_min_til, W_max_til, P_min)

    ! 2) Perform a timestep to get W_til one timestep further, still making sure it is within the bounds
    call calc_W_til_next(mesh, ice, basal_hydro, W_min_til, W_max_til, dt)

    call calc_R(mesh, ice, basal_hydro, .true.)

    call calc_K(mesh, basal_hydro)

    ! 4) Get values on staggered grid 
    call map_all_a_b(mesh, basal_hydro) !calc_D is in here

    call calc_uv(mesh, ice, basal_hydro)

    ! 8) Get the timestep 
    ! Mainly inspired by calc_critical_timestep_SIA subroutine in time_step_criteria.f90
    call get_basal_hydro_timestep(mesh, basal_hydro, dt, dt_hydro)

    if (par%primary) then
      write(*,*) "dt_hydro = ", dt_hydro
    end if

    ! 9) Compute the advective fluxes (Q) on the staggered grid
    call calc_divQ(mesh, ice, basal_hydro)

    ! 11) If icefree set next timestep of P to 0, if floating set to overburden pressure
    ! 11) If W at this timestep is 0 and if icefree and floating are both false, set next timestep of P to 0 (any sliding) or overburden pressure (no sliding)
    ! 11) Otherwise, compute next timestep of P using the equation in the paper
    call calc_P_next(mesh, ice, basal_hydro, P_min, dt_hydro)

    ! 13) If icefree or float, then set next timestep of W to 0.
    ! 13) Otherwise, compute next timestep of W using the equation in the paper
    call calc_W_next(mesh, ice, basal_hydro, W_min, W_max, dt_hydro)

    ! 15) Update time and repeat
    if (par%primary) then
      write(*,*) "Time after basal hydrology step: ", time
    end if

    !call crash('Hello world!')
  
    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine basal_hydrology


  subroutine allocate_basal_hydro( mesh, basal_hydro)
    !< allocate memory for the basal hydrology model variables >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_basal_hydrology_model),   intent(  out) :: basal_hydro

    ! Local variables:
    character(len=1024) :: routine_name = 'allocate_basal_hydro'
    integer             :: vi, ti

    ! Add routine to path
    call init_routine( routine_name)

    allocate(basal_hydro%P_o(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%W(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%W_til(mesh%vi1:mesh%vi2), source =  1.0_dp)
    allocate(basal_hydro%W_til_next(mesh%vi1:mesh%vi2), source =  0.0_dp)
    allocate(basal_hydro%P(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%m(mesh%vi1:mesh%vi2), source = 0.0_dp) ! For now just make this a value
    allocate(basal_hydro%dW_dx_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%W_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%K(mesh%vi1:mesh%vi2), source = 0.001_dp)
    allocate(basal_hydro%D(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%u(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%v(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%dK_dx_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%dD_dx_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%du_dx_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%dv_dx_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%K_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%D_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%u_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%v_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%u_c(mesh%ei1:mesh%ei2), source = 0.0_dp)
    allocate(basal_hydro%v_c(mesh%ei1:mesh%ei2), source = 0.0_dp)
    allocate(basal_hydro%Z(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%C(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%O(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%divQ( mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%old_time, source = 0.0_dp)
    allocate(basal_hydro%mask_a(mesh%vi1:mesh%vi2), source = .false.)
    allocate(basal_hydro%mask_b(mesh%ti1:mesh%ti2), source = .false.)
    allocate(basal_hydro%R(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%dR_dx(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%dR_dy(mesh%vi1:mesh%vi2), source = 0.0_dp)

    do vi = mesh%vi1, mesh%vi2
      ! Initial basal water depth
      basal_hydro%W( vi) = 2.0_dp + sin(mesh%V(vi, 1)*2_dp*pi/80e3_dp)*cos(mesh%V(vi, 2)*2_dp*pi/80e3_dp)

      ! Vortex
      !basal_hydro%u( vi) = 10*sin(mesh%V(vi, 2)*pi/(2_dp*40000_dp))
      !basal_hydro%v( vi) = 10*cos(mesh%V(vi, 1)*pi/(2_dp*40000_dp))

      ! Decreasing flow
      !basal_hydro%u( vi) = 10_dp - mesh%V(vi, 1)/22500_dp

      ! Only diverging flow up and down
      !if (mesh%V(vi, 2) >= 0.0_dp) then
      !  basal_hydro%v( vi) = 10.0_dp
      !else
      !  basal_hydro%v( vi) = -10.0_dp
      !end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_basal_hydro


  subroutine set_within_bounds( mesh, ice, basal_hydro, W_min, W_max, W_min_til, W_max_til, P_min)
    !< Setting P, W, W_til within their bounds >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(in   ) :: ice
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro
    real(dp),                           intent(in   ) :: W_min, W_max, W_min_til, W_max_til, P_min

    ! Local variables:
    character(len=1024) :: routine_name = 'set_within_bounds'
    integer             :: vi
    real(dp), parameter :: g = 9.81_dp
    real(dp), parameter :: rho_i = 917.0_dp


    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      basal_hydro%P_o( vi)     = rho_i * g * ice%Hi( vi) ! Calculate overburden pressure
      basal_hydro%W( vi)       = min( max( basal_hydro%W( vi),       W_min),       W_max)
      basal_hydro%W_til( vi)   = min( max( basal_hydro%W_til( vi),   W_min_til),   W_max_til)
      basal_hydro%P( vi)       = min( max( basal_hydro%P( vi),       P_min),       basal_hydro%P_o( vi))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_within_bounds


  subroutine calc_W_til_next( mesh, ice, basal_hydro, W_min_til, W_max_til, dt)
    !< Calculating W_til for the next timestep >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(in   ) :: ice
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro
    real(dp),                           intent(in   ) :: W_min_til, W_max_til
    real(dp),                           intent(in   ) :: dt

    ! Local variables:
    character(len=1024) :: routine_name = 'calc_W_til_next'
    integer             :: vi
    real(dp), parameter :: Cd = 0.0_dp
    real(dp), parameter :: rho_w = 1000.0_dp


    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      basal_hydro%W_til_next( vi) = basal_hydro%W_til( vi) + dt*(basal_hydro%m( vi)/rho_w - Cd) ! Timestep
      basal_hydro%W_til_next( vi)   = min( max( basal_hydro%W_til_next( vi),   W_min_til),   W_max_til) ! Make sure within bounds
      if (ice%mask_icefree_land( vi) .or. ice%mask_floating_ice( vi) .or. ice%mask_icefree_ocean( vi)) then ! 3) If icefree or floating, set W_til to zero
        basal_hydro%W_til_next( vi) = 0.0_dp
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_W_til_next



  subroutine map_all_a_b( mesh, basal_hydro)
    !< Mapping all the variables that need to go to b grid to the b grid >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro

    ! Local variables:
    character(len=1024) :: routine_name = 'map_all_a_b'

    ! Add routine to path
    call init_routine( routine_name)

    ! 4) Do this for W
    call ddx_a_b_2D(mesh, basal_hydro%W, basal_hydro%dW_dx_b)
    call map_a_b_2D(mesh, basal_hydro%W, basal_hydro%W_b)

    ! 5) Also do this for effective conductivity K
    call ddx_a_b_2D(mesh, basal_hydro%K, basal_hydro%dK_dx_b)
    call map_a_b_2D(mesh, basal_hydro%K, basal_hydro%K_b)

    ! 6) Do this again for velocity components u and v
    call ddx_a_b_2D(mesh, basal_hydro%u, basal_hydro%du_dx_b)
    call ddx_a_b_2D(mesh, basal_hydro%v, basal_hydro%dv_dx_b)
    call map_a_b_2D(mesh, basal_hydro%u, basal_hydro%u_b)
    call map_a_b_2D(mesh, basal_hydro%v, basal_hydro%v_b)

    ! 7) And lastly for the diffusivity D
    call calc_D(mesh, basal_hydro)  ! Calculate D from W and K
    call ddx_a_b_2D(mesh, basal_hydro%D, basal_hydro%dD_dx_b)
    call map_a_b_2D(mesh, basal_hydro%D, basal_hydro%D_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_all_a_b

  

  subroutine get_basal_hydro_timestep( mesh, basal_hydro, dt, dt_hydro)
    !< Get basal hydrology timestep >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_basal_hydrology_model),   intent(in   ) :: basal_hydro
    real(dp),                           intent(in   ) :: dt
    real(dp),                           intent(  out) :: dt_hydro

    ! Local variables:
    character(len=1024) :: routine_name = 'get_basal_hydro_timestep'
    integer             :: ti, via, vib, vic   ! Triangle index and triangle sides
    real(dp)            :: d_ab, d_bc, d_ca    ! Lengths of triangle sides
    real(dp)            :: d_min               ! Minimum triangle side length
    real(dp)            :: u_t, v_t, D_t       ! Velocity components on triangle
    real(dp)            :: dt_crit_CFL, dt_crit_W, dt_crit_P ! Critical timesteps for different conditions
    real(dp), parameter :: phi = 0.01_dp       ! Englacial porosity)
    real(dp), parameter :: correction_factor = 0.9_dp ! To be on the safe side

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise time step with maximum allowed value
    dt_crit_CFL = C%dt_ice_max
    dt_crit_W = C%dt_ice_max
    dt_crit_P = C%dt_ice_max

    do ti = mesh%ti1, mesh%ti2

      ! Calculate shortest triangle side
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)

      d_ab = norm2( mesh%V( vib,:) - mesh%V( via,:))
      d_bc = norm2( mesh%V( vic,:) - mesh%V( vib,:))
      d_ca = norm2( mesh%V( via,:) - mesh%V( vic,:))

      d_min = minval([ d_ab, d_bc, d_ca])

      ! There is one value of u and v per triangle, so we use those
      ! We also assume dx = dy in this case

      ! Velocities on the triangle
      u_t = abs(basal_hydro%u_b( ti))
      v_t = abs(basal_hydro%v_b( ti))

      dt_crit_CFL = min(dt_crit_CFL, d_min/(2*(u_t + v_t)))

      ! Diffusivity on the triangle (Tijn suggested taking the vertices to get a higher value)
      ! Probs have to exchange_halos
      D_t = basal_hydro%D_b( ti) !max(basal_hydro%D( via), basal_hydro%D( vib), basal_hydro%D( vic))

      ! Instead of max of D, we use D on the triangle (I think this makes sense?)
      dt_crit_W = min(dt_crit_W, d_min**2/(8*D_t))

      dt_crit_P = min(dt_crit_P, 2*phi*d_min**2/(8*D_t))

    end do

    ! Timestep we will use here
    dt_hydro = min(correction_factor*min( dt_crit_CFL, dt_crit_W, dt_crit_P), dt)
    if (par%primary) then
      write(*,*) "dt_crit_CFL = ", dt_crit_CFL
      write(*,*) "dt_crit_W   = ", dt_crit_W
      write(*,*) "dt_crit_P   = ", dt_crit_P
    end if

    !call MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_SIA, 1, MPI_doUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    !dt_crit_SIA = min( C%dt_ice_max, dt_crit_SIA)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine get_basal_hydro_timestep



  subroutine calc_divQ( mesh, ice, basal_hydro)
    !< Calculating the divergence of Q >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(in   ) :: ice
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro

    ! Local variables:
    character(len=1024) :: routine_name = 'calc_divQ'
    integer             :: vi, ci, vj, ei
    logical,  dimension(mesh%nV)   :: mask_grounded_ice_tot
    real(dp), dimension(mesh%nE)   :: u_c_tot, v_c_tot
    real(dp), dimension(mesh%nV)   :: W_tot
    real(dp)                       :: u_perp

    ! Add routine to path
    call init_routine( routine_name)

    ! This is very similar to what is done in laddie_thickness.90 in the compute_divQH subroutine

    call gather_to_all(ice%mask_grounded_ice, mask_grounded_ice_tot)
    call calc_M_b_c(mesh, ice, basal_hydro)
    call map_UV_b_c(mesh, basal_hydro)
    call gather_to_all(basal_hydro%u_c, u_c_tot)
    call gather_to_all(basal_hydro%v_c, v_c_tot)
    call gather_to_all(basal_hydro%W, W_tot)

    ! == Loop over vertices ==
    ! =========================

    DO vi = mesh%vi1, mesh%vi2

      ! Initialise divQ with zeros
      basal_hydro%divQ( vi) = 0_dp

      if (.not. ice%mask_grounded_ice( vi)) cycle ! No flux divergence calculation for non-grounded ice

      ! Loop over all connections of vertex vi
      DO ci = 1, mesh%nC( vi)

        ! Connection ci from vertex vi leads through edge ei to vertex vj
        vj = mesh%C(  vi,ci)

        ! Get edge
        ei = mesh%VE( vi,ci)

        ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section
        u_perp = u_c_tot( ei) * mesh%D_x( vi, ci)/mesh%D( vi, ci) + v_c_tot( ei) * mesh%D_y( vi, ci)/mesh%D( vi, ci)

        ! 10) Compute the flux divergence approximations D_ij (this is already here right?)
        ! Calculate upwind momentum divergence
        ! =============================
        ! u_perp > 0: flow is exiting this vertex into vertex vj
        IF (u_perp > 0) THEN
          basal_hydro%divQ( vi) = basal_hydro%divQ( vi) + mesh%Cw( vi, ci) * u_perp * W_tot( vi) / mesh%A( vi)
        ! u_perp < 0: flow is entering this vertex from vertex vj
        ELSE 
          ! Skip connection if neighbour is not grounded. No flux across grounding line
          ! Can be made more flexible when accounting for partial cells (PMP instead of FCMP)
          IF (mask_grounded_ice_tot( vj)) then
            basal_hydro%divQ( vi) = basal_hydro%divQ( vi) + mesh%Cw( vi, ci) * u_perp * W_tot( vj) / mesh%A( vi)
          END IF
        END IF

      END DO ! DO ci = 1, mesh%nC( vi)

    END DO ! DO vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_divQ


  subroutine calc_P_next( mesh, ice, basal_hydro, P_min, dt_hydro)
    !< Calculating P for the next timestep >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(in   ) :: ice
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro
    real(dp),                           intent(in   ) :: P_min
    real(dp),                           intent(in   ) :: dt_hydro

    ! Local variables:
    character(len=1024) :: routine_name = 'calc_P_next'
    integer             :: vi
    real(dp), parameter :: Cd = 0.000_dp
    real(dp), parameter :: rho_w = 1000.0_dp
    real(dp), parameter :: g = 9.81_dp
    real(dp), parameter :: phi = 0.01_dp            ! Englacial porosity
    logical, parameter  :: sliding = .true.         ! Whether the ice slides or not when W = 0 and there is ice on land


    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_icefree_land( vi)) then
        basal_hydro%P( vi) =  0.0_dp
      else if (ice%mask_floating_ice( vi)) then
        basal_hydro%P( vi) = basal_hydro%P_o( vi)
      else if (basal_hydro%W( vi) == 0.0_dp .and. .not. ice%mask_icefree_land( vi) .and. .not. ice%mask_floating_ice( vi) .and. .not. ice%mask_icefree_ocean( vi)) then
        if (sliding) then
          basal_hydro%P( vi) = 0.0_dp
        else
          basal_hydro%P( vi) = basal_hydro%P_o( vi)
        end if
      else
        ! Compute next timestep of P using the equation in the paper
        basal_hydro%Z( vi) = basal_hydro%C( vi) - basal_hydro%O( vi) + basal_hydro%m( vi)/rho_w  - (basal_hydro%W_til_next( vi) - basal_hydro%W_til( vi))/dt_hydro
        basal_hydro%P( vi) = basal_hydro%P( vi) + dt_hydro * ((rho_w * g / phi) * (-basal_hydro%divQ( vi) + basal_hydro%Z( vi)))
      
    ! 12) Make sure P is within its bounds
        basal_hydro%P( vi) = min( max( basal_hydro%P( vi), P_min), basal_hydro%P_o( vi))
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_P_next


  subroutine calc_W_next( mesh, ice, basal_hydro, W_min, W_max, dt_hydro)
    !< Calculating W for the next timestep >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(in   ) :: ice
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro
    real(dp),                           intent(in   ) :: W_min, W_max
    real(dp),                           intent(in   ) :: dt_hydro

    ! Local variables:
    character(len=1024) :: routine_name = 'calc_W_next'
    integer             :: vi
    real(dp), parameter :: rho_w = 1000.0_dp

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_icefree_land( vi) .or. ice%mask_floating_ice( vi) .or. ice%mask_icefree_ocean( vi)) then
        basal_hydro%W( vi) = 0.0_dp
      else
        basal_hydro%W( vi) = basal_hydro%W( vi) + basal_hydro%W_til( vi) - basal_hydro%W_til_next( vi) &
                              + dt_hydro * (-basal_hydro%divQ( vi) + basal_hydro%m( vi)/rho_w) 
      end if
      ! 14) Make sure W is within its bounds (>= 0)
      basal_hydro%W( vi) = min( max( basal_hydro%W( vi), W_min), W_max)
      ! Set W_til to W_til_next for next timestep
      basal_hydro%W_til( vi) = basal_hydro%W_til_next( vi)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_W_next



  subroutine calc_general_dt( time, basal_hydro, dt)
    !< Calculating the general dt since last timestep >!

    ! In/output variables:
    real(dp),                           intent(in   ) :: time
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro
    real(dp),                           intent(  out) :: dt

    ! Local variables:
    character(len=1024) :: routine_name = 'calc_general_dt'


    ! Add routine to path
    call init_routine( routine_name)

    if (time == 0.0_dp) then
      dt = C%dt_ice_max
    else
      dt = time - basal_hydro%old_time
    end if
    basal_hydro%old_time = time

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_general_dt
  

  subroutine calc_D( mesh, basal_hydro)
    !< Calculating the diffusivity D>!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro

    ! Local variables:
    character(len=1024) :: routine_name = 'calc_D'
    integer             :: vi
    real(dp), parameter :: g = 9.81_dp
    real(dp), parameter :: rho_w = 1000.0_dp


    ! Add routine to path
    call init_routine( routine_name)
    !Some are NaNs here. This is due to K being NaN sometimes. But why? Probably because dR_dx and dR_dy are zero.

    do vi = mesh%vi1, mesh%vi2
      if (basal_hydro%mask_a( vi)) then
        basal_hydro%D( vi) = rho_w*g*basal_hydro%K( vi)*basal_hydro%W( vi)
      else
        basal_hydro%D( vi) = 0.0_dp
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_D



  ! Comes from Laddie_velocity
  subroutine map_UV_b_c( mesh, basal_hydro)
    ! Calculate velocities on the c-grid
    !
    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'map_UV_b_c'

    ! Add routine to path
    call init_routine( routine_name)

    call sync

    call multiply_CSR_matrix_with_vector_1D(basal_hydro%M_b_c, &
      mesh%pai_Tri, basal_hydro%u_b, mesh%pai_E, basal_hydro%u_c)
    call multiply_CSR_matrix_with_vector_1D( basal_hydro%M_b_c, &
      mesh%pai_Tri, basal_hydro%v_b, mesh%pai_E, basal_hydro%v_c)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_UV_b_c



  subroutine map_UV_a_c( mesh, basal_hydro)
    ! Calculate velocities from a-grid to the c-grid
    !
    ! Uses a different scheme then the standard mapping operator, as that one is too diffusive

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'map_UV_a_c'

    ! Add routine to path
    call init_routine( routine_name)

    call sync



    call multiply_CSR_matrix_with_vector_1D(basal_hydro%M_a_c, &
      mesh%pai_V, basal_hydro%u, mesh%pai_E, basal_hydro%u_c)
    call multiply_CSR_matrix_with_vector_1D( basal_hydro%M_a_c, &
      mesh%pai_V, basal_hydro%v, mesh%pai_E, basal_hydro%v_c)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_UV_a_c


  ! Comes from Laddie_operators
  subroutine calc_M_b_c( mesh, ice, basal_hydro)
    ! Calculate mapping matrix from b-grid to c-grid

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_ice_model),                   intent(in   )    :: ice
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'calc_M_b_c'
    integer                                               :: ncols, ncols_loc, nrows, nrows_loc, nnz_per_row_est, nnz_est_proc
    integer                                               :: row, ti, n, i, vi, vj, ei, til, tir
    logical, dimension(mesh%nTri)                         :: mask_b_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Get the basal hydrology mask on (a) and b-grid
    call calc_basal_hydro_mask_a_b(mesh, ice, basal_hydro)
    call gather_to_all(basal_hydro%mask_b, mask_b_tot)

    call deallocate_matrix_CSR_dist( basal_hydro%M_b_c)

    ! Matrix size
    ncols           = mesh%nTri        ! from
    ncols_loc       = mesh%nTri_loc
    nrows           = mesh%nE        ! to
    nrows_loc       = mesh%nE_loc
    nnz_per_row_est = 2
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( basal_hydro%M_b_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = mesh%pai_Tri, pai_y = mesh%pai_E)

    ! == Calculate coefficients
    ! =========================

    do row = basal_hydro%M_b_c%i1, basal_hydro%M_b_c%i2

      ! The vertex represented by this matrix row
      ei = mesh%n2ei( row)

      ! Get neighbouring triangles
      til = mesh%ETri( ei, 1)
      tir = mesh%ETri( ei, 2)

      if (til == 0 .and. tir > 0) then
        ! Only triangle on right side exists
        if (mask_b_tot( tir)) then
          ! Within basal_hydro domain, so add
          call add_entry_CSR_dist( basal_hydro%M_b_c, ei, tir, 1._dp)
        else
          ! Outside basal_hydro domain, so omit
          call add_empty_row_CSR_dist( basal_hydro%M_b_c, ei)
        end if
      elseif (tir == 0 .and. til > 0) then
        ! Only triangle on left side exists
        if (mask_b_tot( til)) then
          ! Within basal_hydro domain, so add
          call add_entry_CSR_dist( basal_hydro%M_b_c, ei, til, 1._dp)
        else
          ! Outside basal_hydro domain, so omit
          call add_empty_row_CSR_dist( basal_hydro%M_b_c, ei)
        end if
      elseif (til > 0 .and. tir > 0) then
        ! Both triangles exist
        if (mask_b_tot( til) .or. mask_b_tot( tir)) then
          ! At least one traingle in basal_hydro domain, so add average
          call add_entry_CSR_dist( basal_hydro%M_b_c, ei, til, 0.5_dp)
          call add_entry_CSR_dist( basal_hydro%M_b_c, ei, tir, 0.5_dp)
        else
          ! Both outside basal_hydro domain, so omit
          call add_empty_row_CSR_dist( basal_hydro%M_b_c, ei)
        end if
      else
          call crash('something is seriously wrong with the ETri array of this mesh!')
      end if

    end do

    ! Crop matrix memory
    call finalise_matrix_CSR_dist( basal_hydro%M_b_c)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_M_b_c



  ! Comes from Laddie_operators
  subroutine calc_M_a_c( mesh, ice, basal_hydro)
    ! Calculate mapping matrix from a-grid to c-grid

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_ice_model),                   intent(in   )    :: ice
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'calc_M_a_c'
    integer                                               :: ncols, ncols_loc, nrows, nrows_loc, nnz_per_row_est, nnz_est_proc
    integer                                               :: row, vi, vj, ei
    logical, dimension(mesh%nV)                           :: mask_a_tot
    real(dp), dimension(2)                                :: cM_a_c

    ! Add routine to path
    call init_routine( routine_name)

    ! Get the basal hydrology mask on a (and b) grid
    call calc_basal_hydro_mask_a_b(mesh, ice, basal_hydro)
    call gather_to_all(basal_hydro%mask_a, mask_a_tot)

    call deallocate_matrix_CSR_dist( basal_hydro%M_a_c)

    ! == Initialise the matrix using the native UFEMISM CSR-matrix format
    ! ===================================================================

    ! Matrix size
    ncols           = mesh%nV        ! from
    ncols_loc       = mesh%nV_loc
    nrows           = mesh%nE        ! to
    nrows_loc       = mesh%nE_loc
    nnz_per_row_est = 2
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( basal_hydro%M_a_c, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = mesh%pai_V, pai_y = mesh%pai_E)

    ! == Calculate coefficients
    ! =========================

    do row = basal_hydro%M_a_c%i1, basal_hydro%M_a_c%i2

      ! The vertex represented by this matrix row
      ei = mesh%n2ei( row)

      ! Get neighbouring vertices
      vi = mesh%EV( ei, 1)
      vj = mesh%EV( ei, 2)

      ! Get masked average between the two vertices
      if (mask_a_tot( vi) .and. mask_a_tot( vj)) then
        cM_a_c = [0.5_dp, 0.5_dp]
      elseif (mask_a_tot( vi)) then
        cM_a_c = [1._dp, 0._dp]
      elseif (mask_a_tot( vj)) then
        cM_a_c = [0._dp, 1._dp]
      else
        cM_a_c = 0._dp
      end if

      ! Add weight to matrix
      call add_entry_CSR_dist( basal_hydro%M_a_c, ei, vi, cM_a_c( 1))
      call add_entry_CSR_dist( basal_hydro%M_a_c, ei, vj, cM_a_c( 2))

    end do

    ! Crop matrix memory
    call finalise_matrix_CSR_dist( basal_hydro%M_a_c)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_M_a_c



  subroutine calc_basal_hydro_mask_a_b(mesh, ice, basal_hydro)
    !< Calculating the basal hydrology mask on a and b-grid >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(in   ) :: ice
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro

    ! Local variables:
    character(len=1024)                                   :: routine_name = 'calc_basal_hydro_mask_a_b'
    integer                                               :: vi, i, ti
    logical, dimension(mesh%nV)                           :: mask_a_tot

    ! Add routine to path
    call init_routine( routine_name)
    
    ! Define a-grid mask directly from ice model grounded ice mask
    do vi = mesh%vi1, mesh%vi2
      basal_hydro%mask_a( vi) = .false.
      basal_hydro%mask_a( vi) = ice%mask_grounded_ice( vi)
    end do

    call gather_to_all(basal_hydro%mask_a, mask_a_tot)

    ! Define grounded if any of the three vertices is grounded
    do ti = mesh%ti1, mesh%ti2
      basal_hydro%mask_b( ti) = .false.
      DO i = 1, 3
        vi = mesh%Tri( ti, i)
        IF (mask_a_tot( vi)) THEN
          ! Set true if any of the three vertices is grounded
          basal_hydro%mask_b( ti) = .true.
        END IF
      END DO
    end do
  
    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_basal_hydro_mask_a_b


  subroutine calc_R( mesh, ice, basal_hydro, test)
    ! Calculate subglacial water pressure R

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro
    type(type_ice_model),                   intent(in   )    :: ice
    logical,                                intent(in   )    :: test

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'calc_R'
    integer                                               :: vi
    real(dp), parameter                                   :: g = 9.81_dp
    real(dp), parameter                                   :: rho_w = 1000.0_dp

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      if (test) then                          ! For now we take R without the pressure component
        basal_hydro%R( vi) = ice%Hb( vi)*rho_w*g
      else 
        basal_hydro%R( vi) = ice%Hb( vi)*rho_w*g + basal_hydro%P( vi)
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_R


  subroutine calc_K( mesh, basal_hydro)
    ! Calculate K

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'calc_K'
    integer                                               :: vi
    real(dp), parameter                                   :: k = 0.001_dp
    real(dp), parameter                                   :: alpha = 1.25_dp
    real(dp), parameter                                   :: beta = 1.5_dp

    ! Add routine to path
    call init_routine( routine_name)

    call ddx_a_a_2D(mesh, basal_hydro%R, basal_hydro%dR_dx)
    call ddy_a_a_2D(mesh, basal_hydro%R, basal_hydro%dR_dy)

    ! For some reason the dR_dx and dR_dy values are zero sometimes and if this is precisely 0, this breaks calc_K by dividing by zero.
    do vi = mesh%vi1, mesh%vi2
      if (basal_hydro%mask_a( vi)) then
        basal_hydro%K( vi) = k*basal_hydro%W( vi)**(alpha - 1._dp)*abs(basal_hydro%dR_dx( vi)**2._dp + basal_hydro%dR_dy( vi)**2._dp + 0.00000001_dp)**((beta - 2._dp)/2._dp)
      else
        basal_hydro%K( vi) = 0.0_dp
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_K



  subroutine calc_uv( mesh, ice, basal_hydro)
    ! Calculate u and v from K and R

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_ice_model),                   intent(in   )    :: ice
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'calc_uv'
    integer                                               :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Get the derivatives of R
    call ddx_a_a_2D(mesh, basal_hydro%R, basal_hydro%dR_dx)
    call ddy_a_a_2D(mesh, basal_hydro%R, basal_hydro%dR_dy)

    ! Calculate u and v
    do vi = mesh%vi1, mesh%vi2 ! Not using Ke and Kn as in the paper, but just K on the vertex
      basal_hydro%u( vi) = (- basal_hydro%K( vi) * basal_hydro%dR_dx( vi))!*sec_per_year
      basal_hydro%v( vi) = (- basal_hydro%K( vi) * basal_hydro%dR_dy( vi))!*sec_per_year
    end do

    ! Remap to c-grid velocities
    call calc_M_a_c( mesh, ice, basal_hydro)
    call map_UV_a_c( mesh, basal_hydro)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_uv



  subroutine point_source( mesh, basal_hydro)
    ! Calculate a point source for testing. Allocate W with zeros for best viewable results.

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'point_source'
    integer                                               :: vi, vi_point
    real(dp), dimension(2)                                :: point

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise for 1 point test
    point = [100000.0_dp, 0.0_dp]
    vi_point = 1

    ! Find the vertex closest to the point
    call find_containing_vertex(mesh, point, vi_point)

    ! Set a point source there and leave everything else zero
    if (vi_point >= mesh%vi1 .and. vi_point <= mesh%vi2) then
      basal_hydro%W( vi_point) = 0.1_dp
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine point_source

END MODULE basal_hydrology_new