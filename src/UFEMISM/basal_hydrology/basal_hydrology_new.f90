MODULE basal_hydrology_new

! ===== Preamble =====
! ====================
! This preamble is still copied from thermodynamics_3D_heat_equation.f90

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MIN
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE call_stack_and_comp_time_tracking                      , ONLY: init_routine, finalise_routine
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE basal_hydrology_model_types                            , ONLY: type_basal_hydrology_model
  use netcdf_io_main
  USE mesh_disc_apply_operators                              , ONLY: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, map_b_a_2D, ddx_a_a_2D, ddy_a_a_2D
  use mpi_distributed_memory                                 , only: gather_to_all
  use mesh_halo_exchange                                     , only: exchange_halos
  use CSR_matrix_vector_multiplication                       , only: multiply_CSR_matrix_with_vector_1D_wrapper
  use mesh_utilities                                         , only: find_containing_vertex
  use CSR_matrix_basics                                      , only: finalise_matrix_CSR_dist, add_entry_CSR_dist, add_empty_row_CSR_dist, allocate_matrix_CSR_dist, deallocate_matrix_CSR_dist
  use conservation_of_mass_utilities                         , only: calc_n_interior_neighbours
  use crash_mod                                              , only: crash, warning, happy
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use remapping_main                                         , only: map_from_mesh_to_mesh_with_reallocation_2D
  use checksum_mod                                           , only: checksum

  IMPLICIT NONE

CONTAINS

  ! Run basal hydrology model x times
  subroutine basal_hydrology_leg(mesh, ice, basal_hydro, time)
    !< Main basal hydrology subroutine >!

    ! In/output variables:
    type(type_mesh),                  intent(in   ) :: mesh
    type(type_ice_model),             intent(inout) :: ice
    type(type_basal_hydrology_model), intent(inout) :: basal_hydro
    real(dp),                         intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'basal_hydrology_leg'
    real(dp)                       :: duration
    real(dp)                       :: time_until_convergence
    real(dp)                       :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Start counter on 0
    duration = 0.0_dp

    time_until_convergence = C%basal_hydro_equil_time*sec_per_year  ! Time until convergence 

    ! Loop until convergence time is reached
    time_loop: do while (duration < time_until_convergence)
      call basal_hydrology(mesh, ice, basal_hydro, duration + time*sec_per_year)
      duration = duration + basal_hydro%dt
      !if (par%primary) then
      !  write(*,*) "Duration so far basal hydro: ", duration
      !end if
    end do time_loop

    basal_hydro%old_time = time*sec_per_year + duration

    if (par%primary) then
      !write(*,*) "Done!"
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine basal_hydrology_leg


  ! Run a single leg of basal hydrology model
  subroutine basal_hydrology(mesh, ice, basal_hydro, time)
    ! In/output variables:
    type(type_mesh),                  intent(in   ) :: mesh
    type(type_ice_model),             intent(inout) :: ice
    type(type_basal_hydrology_model), intent(inout) :: basal_hydro
    real(dp),                         intent(in)    :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'basal_hydrology'

    real(dp), parameter            :: W_max = 1000.0_dp     ! Maximum basal water depth (number is placeholder for now, because there is no maximum W)
    real(dp), parameter            :: W_min = 0.0_dp           ! Minimum basal water depth

    real(dp), parameter            :: W_max_til = 2.0_dp       ! Maximum basal water depth in till
    real(dp), parameter            :: W_min_til = 0.0_dp       ! Minimum basal water depth in till

    real(dp), parameter            :: P_min = 0.0_dp           ! Minimum pressure of the ice on the basal water

    real(dp)                       :: dt
    real(dp)                       :: dt_hydro

    ! Add routine to path
    call init_routine( routine_name)

    call convert_ice_to_SI(mesh, ice, basal_hydro)

    ! Initialise basal hydro masks
    call calc_basal_hydro_mask_a_b(mesh, ice, basal_hydro)

    ! 1) Start with W, W_til and P and make sure they are all within their bounds
    call set_within_bounds(mesh, ice, basal_hydro, W_min, W_max, W_min_til, W_max_til, P_min)

    ! Calculate some stuff for timesteps and u and v calculation
    call calc_R(mesh, ice, basal_hydro, .false.)

    call calc_K(mesh, ice, basal_hydro)

    call calc_D(mesh, basal_hydro)

    call calc_uv(mesh, ice, basal_hydro)

    ! 8) Get the timestep 
    ! Mainly inspired by calc_critical_timestep_SIA subroutine in time_step_criteria.f90
    call get_basal_hydro_timestep(mesh, basal_hydro, dt, dt_hydro)

    ! 9) Compute the advective fluxes (Q) on the staggered grid
    call calc_divQ(mesh, ice, basal_hydro)

    ! Compute how much goes in the water layer and how much goes in the till
    call calc_q_til(mesh, ice, basal_hydro, W_max_til)
    
    ! 11) If icefree set next timestep of P to 0, if floating set to overburden pressure
    ! 11) If W at this timestep is 0 and if icefree and floating are both false, set next timestep of P to 0 (any sliding) or overburden pressure (no sliding)
    ! 11) Otherwise, compute next timestep of P using the equation in the paper (Bueler and Van Pelt 2015)
    call calc_P_next(mesh, ice, basal_hydro, P_min)

    ! 13) If icefree or float, then set next timestep of W to 0.
    ! 13) Otherwise, compute next timestep of W using the equation in the paper
    ! Put the computed q to use
    call calc_W_water_W_til_next(mesh, ice, basal_hydro, W_min, W_max, W_min_til, W_max_til)

    ! Calculate output to ice model (effective pressure and yield stress)
    call calc_N_til(mesh, basal_hydro, W_max_til, .true.)
    call calc_yield_stress(mesh, ice, basal_hydro)

    ! Allow boundary conditions to be applied to W
    !call apply_W_thickness_BC_explicit(mesh, ice, basal_hydro)
    !call apply_W_thickness_gl_explicit(mesh, ice, basal_hydro)

    ! 15) Update time and repeat
    if (par%primary) then
      !write(*,*) "Time after basal hydrology step: ", time
    end if
  
    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine basal_hydrology


  subroutine allocate_basal_hydro( mesh, ice, basal_hydro)
    !< allocate memory for the basal hydrology model variables >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(in   ) :: ice
    type(type_basal_hydrology_model),   intent(  out) :: basal_hydro

    ! Local variables:
    character(len=1024) :: routine_name = 'allocate_basal_hydro'
    integer             :: vi, ti
    real(dp), parameter :: g = 9.81_dp
    real(dp), parameter :: rho_i = 917.0_dp
    real(dp), parameter :: rho_w = 1000.0_dp

    ! Add routine to path
    call init_routine( routine_name)

    allocate(basal_hydro%P_o(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%W(mesh%vi1:mesh%vi2), source = 0.01_dp)
    allocate(basal_hydro%W_til(mesh%vi1:mesh%vi2), source =  2.0_dp)
    allocate(basal_hydro%W_til_next(mesh%vi1:mesh%vi2), source =  0.0_dp)
    allocate(basal_hydro%P(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%m(mesh%vi1:mesh%vi2), source = 0.0069_dp*rho_w/sec_per_year) ! basal melt rate kg m^-2 s^-1 (0.0069 m/yr water equivalent)
    allocate(basal_hydro%m_extra(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%dW_dx_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%W_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%K(mesh%vi1:mesh%vi2), source = 0.0_dp)
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
    allocate(basal_hydro%mask_W(mesh%vi1:mesh%vi2), source = .false.)
    allocate(basal_hydro%R(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%dR_dx_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%dR_dy_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%t_next, source = C%start_time_of_run)
    allocate(basal_hydro%dt, source = 0.0_dp)
    allocate(basal_hydro%diff_time, source = 0.0_dp)
    allocate(basal_hydro%W_r(mesh%vi1:mesh%vi2), source = 0.1_dp) !Value used in basal hydrology paper (Bueler and Van Pelt 2015)
    allocate(basal_hydro%Y(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%q_til(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%q_water_layer(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%Cd, source = 0.001_dp/sec_per_year) ! Value for water leaking back from till to water layer (m/s)
    allocate(basal_hydro%ice_u_base(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%ice_v_base(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%ice_w_base(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%N_til(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%tau_c(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%phi(mesh%vi1:mesh%vi2), source = 26.565_dp) !degrees

    do vi = mesh%vi1, mesh%vi2
      ! Initial basal water depth
      !basal_hydro%W( vi) = 2.0_dp + sin(mesh%V(vi, 1)*2_dp*pi/80e3_dp)*cos(mesh%V(vi, 2)*2_dp*pi/80e3_dp)

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

      ! Initialise pressure with overburden pressure (for now at least)
      !basal_hydro%P( vi) = rho_i * g * ice%Hi( vi)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_basal_hydro



  subroutine convert_ice_to_SI(mesh, ice, basal_hydro)
    !< Convert ice model variables to SI units for basal hydrology model >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(in   ) :: ice
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro

    ! Local variables:
    character(len=1024) :: routine_name = 'convert_ice_to_SI'
    integer             :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      ! Convert ice velocities to m/s
      basal_hydro%ice_u_base( vi) = ice%u_base( vi)/sec_per_year  ! Convert to m/s
      basal_hydro%ice_v_base( vi) = ice%v_base( vi)/sec_per_year  ! Convert to m/s
      basal_hydro%ice_w_base( vi) = ice%w_base( vi)/sec_per_year  ! Convert to m/s
    end do
    !call checksum(mesh%pai_V, basal_hydro%ice_u_base, "ice_u_base after conversion to SI")
    !call checksum(mesh%pai_V, basal_hydro%ice_v_base, "ice_v_base after conversion to SI")
    !call checksum(mesh%pai_V, basal_hydro%ice_w_base, "ice_w_base after conversion to SI")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine convert_ice_to_SI


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
    !!call checksum(mesh%pai_V, basal_hydro%P_o, "basal_hydro%P_o")
    !call checksum(mesh%pai_V, basal_hydro%W, "basal_hydro%W")
    !call checksum(mesh%pai_V, basal_hydro%W_til, "basal_hydro%W_til")
    !call checksum(mesh%pai_V, basal_hydro%P, "basal_hydro%P")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_within_bounds



  subroutine get_basal_hydro_timestep( mesh, basal_hydro, dt, dt_hydro)
    !< Get basal hydrology timestep >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro
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
    dt_crit_CFL = C%dt_ice_max*sec_per_year
    dt_crit_W = C%dt_ice_max*sec_per_year
    dt_crit_P = C%dt_ice_max*sec_per_year

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

      ! Diffusivity on the triangle (Tijn suggested taking the vertices to get a higher value (not done yet))
      ! Now added a small value to avoid division by zero, but maybe using if basal_hydro%mask_b(ti) would be better?
      D_t = basal_hydro%D_b( ti) + 0.0000001_dp !max(basal_hydro%D( via), basal_hydro%D( vib), basal_hydro%D( vic))

      ! Instead of max of D, we use D on the triangle
      dt_crit_W = min(dt_crit_W, d_min**2/(8*D_t))

      dt_crit_P = min(dt_crit_P, 2*phi*d_min**2/(8*D_t))

    end do

    ! Timestep we will use here
    dt_hydro = correction_factor*min( dt_crit_CFL, dt_crit_W, dt_crit_P)
    call MPI_ALLREDUCE( MPI_IN_PLACE, dt_hydro, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD)
    basal_hydro%dt = dt_hydro
    
    if (par%i == 1) then
      !write(*,*) "dt_crit_CFL = ", dt_crit_CFL
      !write(*,*) "dt_crit_W   = ", dt_crit_W
      !write(*,*) "dt_crit_P   = ", dt_crit_P
    else
      !write(*,*) "dt_crit_CFL = ", dt_crit_CFL, " 0"
      !write(*,*) "dt_crit_W   = ", dt_crit_W, " 0"
      !write(*,*) "dt_crit_P   = ", dt_crit_P, " 0"
    end if

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
    call gather_to_all(basal_hydro%u_c, u_c_tot)
    call gather_to_all(basal_hydro%v_c, v_c_tot)
    call gather_to_all(basal_hydro%W, W_tot)

    ! == Loop over vertices ==
    ! =========================

    DO vi = mesh%vi1, mesh%vi2

      ! Initialise divQ with zeros
      basal_hydro%divQ( vi) = 0.0_dp

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

    !call checksum(mesh%pai_V, basal_hydro%divQ, "basal_hydro%divQ")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_divQ


  subroutine calc_P_next( mesh, ice, basal_hydro, P_min)
    !< Calculating P for the next timestep >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(in   ) :: ice
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro
    real(dp),                           intent(in   ) :: P_min

    ! Local variables:
    character(len=1024) :: routine_name = 'calc_P_next'
    integer             :: vi
    real(dp), parameter :: Cd = 0.000_dp
    real(dp), parameter :: rho_w = 1000.0_dp
    real(dp), parameter :: g = 9.81_dp
    real(dp), parameter :: phi = 0.01_dp            ! Englacial porosity
    logical, parameter  :: sliding = .false.         ! Whether the ice slides or not when W = 0 and there is ice on land


    ! Add routine to path
    call init_routine( routine_name)

    ! calculate opening O and closing C terms
    call calc_opening_rate(mesh, ice, basal_hydro)
    call calc_closing_rate(mesh, ice, basal_hydro)

    ! Calculate next timestep of P
    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_icefree_land( vi)) then
        basal_hydro%P( vi) =  0.0_dp
      else if (ice%mask_floating_ice( vi)) then
        basal_hydro%P( vi) = basal_hydro%P_o( vi)
      else if (basal_hydro%W( vi) == 0.0_dp .and. basal_hydro%mask_a( vi)) then
        if (sliding) then
          basal_hydro%P( vi) = 0.0_dp
        else
          basal_hydro%P( vi) = basal_hydro%P_o( vi)
        end if
      else
        ! Compute next timestep of P using the equation in the paper
        basal_hydro%Z( vi) = basal_hydro%C( vi) - basal_hydro%O( vi) + basal_hydro%q_water_layer( vi)/basal_hydro%dt
        basal_hydro%P( vi) = basal_hydro%P( vi) + basal_hydro%dt * ((rho_w * g / phi) * (-basal_hydro%divQ( vi) + basal_hydro%Z( vi)))
      
    ! 12) Make sure P is within its bounds
        basal_hydro%P( vi) = min( max( basal_hydro%P( vi), P_min), basal_hydro%P_o( vi))
      end if
    end do
    !call checksum(mesh%pai_V, basal_hydro%P, "basal_hydro%P")
    !call checksum(mesh%pai_V, basal_hydro%Z, "basal_hydro%Z")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_P_next



  subroutine calc_q_til( mesh, ice, basal_hydro, W_max_til)
    !< Calculating the amount of water ending up in the till and the water layer >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(in   ) :: ice
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro
    real(dp),                           intent(in   ) :: W_max_til

    ! Local variables:
    character(len=1024) :: routine_name = 'calc_q_til'
    integer             :: vi
    real(dp), parameter :: rho_w = 1000.0_dp

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      ! Calculate what portion of water goes into till and what portion into water layer
      basal_hydro%q_til( vi) = min(W_max_til - basal_hydro%W_til( vi) + basal_hydro%Cd*basal_hydro%dt, ((basal_hydro%m( vi)/rho_w)*basal_hydro%dt))
      basal_hydro%q_water_layer( vi) = (basal_hydro%Cd + basal_hydro%m( vi)/rho_w)*basal_hydro%dt - basal_hydro%q_til( vi)
    end do
    !call checksum(mesh%pai_V, basal_hydro%q_til, "basal_hydro%q_til")
    !call checksum(mesh%pai_V, basal_hydro%q_water_layer, "basal_hydro%q_water_layer")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_q_til



  subroutine calc_W_water_W_til_next( mesh, ice, basal_hydro, W_min, W_max, W_min_til, W_max_til)
    !< Calculating W in till and in water layer for the next timestep >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_model),               intent(in   ) :: ice
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro
    real(dp),                           intent(in   ) :: W_min, W_max, W_min_til, W_max_til

    ! Local variables:
    character(len=1024) :: routine_name = 'calc_W_water_W_til_next'
    integer             :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      if (.not. basal_hydro%mask_a( vi)) then
        basal_hydro%W( vi) = 0.0_dp
        basal_hydro%W_til( vi) = 0.0_dp
      else
        basal_hydro%W( vi) = basal_hydro%W( vi) + basal_hydro%q_water_layer( vi) &
                               + basal_hydro%dt * (-basal_hydro%divQ( vi))
        basal_hydro%W_til( vi) = basal_hydro%W_til( vi) + basal_hydro%q_til( vi)
      end if
      ! 14) Make sure W is within its bounds (>= 0)
      basal_hydro%W( vi) = min( max( basal_hydro%W( vi), W_min), W_max)
      basal_hydro%W_til( vi) = min( max( basal_hydro%W_til( vi),   W_min_til),   W_max_til)
    end do
    !call checksum(mesh%pai_V, basal_hydro%W, "basal_hydro%W")
    !call checksum(mesh%pai_V, basal_hydro%W_til, "basal_hydro%W_til")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_W_water_W_til_next



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

    if (par%i == 0) then
      !write(*,*) "Time = ", time
    end if
    if (time == 0.0_dp .or. basal_hydro%dt == 0.0_dp) then
      dt = 1_dp!C%dt_ice_max*sec_per_year
    else
      dt = (time - basal_hydro%old_time)*sec_per_year
    end if
    basal_hydro%old_time = time
    ! Get the next timestep time of ice model
    basal_hydro%t_next = time + dt

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
    integer             :: ti
    real(dp), parameter :: g = 9.81_dp
    real(dp), parameter :: rho_w = 1000.0_dp


    ! Add routine to path
    call init_routine( routine_name)

    call map_a_b_2D(mesh, basal_hydro%W, basal_hydro%W_b)

    do ti = mesh%ti1, mesh%ti2
      if (basal_hydro%mask_b( ti)) then
        basal_hydro%D_b( ti) = rho_w*g*basal_hydro%K_b( ti)*basal_hydro%W_b( ti)
      else
        basal_hydro%D_b( ti) = 0.0_dp
      end if
    end do
    !call checksum(mesh%pai_Tri, basal_hydro%D_b, "basal_hydro%D_b")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_D



  ! Comes from Laddie_velocity
  subroutine map_UV_b_c( mesh, basal_hydro)
    ! Calculate velocities on the c-grid

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'map_UV_b_c'

    ! Add routine to path
    call init_routine( routine_name)

    call multiply_CSR_matrix_with_vector_1D_wrapper(basal_hydro%M_b_c, &
      mesh%pai_Tri, basal_hydro%u_b, mesh%pai_E, basal_hydro%u_c)
    call multiply_CSR_matrix_with_vector_1D_wrapper( basal_hydro%M_b_c, &
      mesh%pai_Tri, basal_hydro%v_b, mesh%pai_E, basal_hydro%v_c)
    !call checksum(mesh%pai_E, basal_hydro%u_c, "basal_hydro%u_c")
    !call checksum(mesh%pai_E, basal_hydro%v_c, "basal_hydro%v_c")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_UV_b_c



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
    !call checksum(mesh%pai_V, basal_hydro%mask_a, "basal_hydro%mask_a")
    !call checksum(mesh%pai_Tri, basal_hydro%mask_b, "basal_hydro%mask_b")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_basal_hydro_mask_a_b



  subroutine calc_basal_hydro_mask_W(mesh, basal_hydro)
    !< Calculating the mask for where the water layer is present >!

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_basal_hydrology_model),   intent(inout) :: basal_hydro

    ! Local variables:
    character(len=1024)                                   :: routine_name = 'calc_basal_hydro_mask_W'
    integer                                               :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      basal_hydro%mask_W( vi) = .false.
      if (basal_hydro%W( vi) > 0.0_dp) then
        basal_hydro%mask_W( vi) = .true.
      end if
    end do
    !call checksum(mesh%pai_V, basal_hydro%mask_W, "basal_hydro%mask_W")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_basal_hydro_mask_W


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
      if (test) then                          ! For testing we take R without the pressure component
        basal_hydro%R( vi) = (ice%Hb( vi) + basal_hydro%W( vi))*rho_w*g
      else 
        basal_hydro%R( vi) = (ice%Hb( vi) + basal_hydro%W( vi))*rho_w*g + basal_hydro%P( vi)
      end if
    end do
    !call checksum(mesh%pai_V, basal_hydro%R, "basal_hydro%R")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_R


  subroutine calc_K( mesh, ice, basal_hydro)
    ! Calculate K

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_ice_model),                   intent(in   )    :: ice
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'calc_K'
    integer                                               :: ti
    real(dp), parameter                                   :: k = 0.001_dp
    real(dp), parameter                                   :: alpha = 1.25_dp
    real(dp), parameter                                   :: beta = 1.5_dp

    ! Add routine to path
    call init_routine( routine_name)

    call ddx_a_b_2D(mesh, basal_hydro%R, basal_hydro%dR_dx_b)
    call ddy_a_b_2D(mesh, basal_hydro%R, basal_hydro%dR_dy_b)

    call map_a_b_2D(mesh, basal_hydro%W, basal_hydro%W_b)

    do ti = mesh%ti1, mesh%ti2
      if (basal_hydro%mask_b( ti)) then
        basal_hydro%K_b( ti) = k*basal_hydro%W_b( ti)**(alpha - 1._dp)*abs(basal_hydro%dR_dx_b( ti)**2._dp & 
                               + basal_hydro%dR_dy_b( ti)**2._dp + 0.00000001_dp)**((beta - 2._dp)/2._dp) ! Added small value to avoid dividing by 0
      else
        basal_hydro%K_b( ti) = 0.0_dp
      end if
    end do
    !call checksum(mesh%pai_Tri, basal_hydro%K_b, "basal_hydro%K_b")

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
    integer                                               :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! Get the derivatives of R
    call ddx_a_b_2D(mesh, basal_hydro%R, basal_hydro%dR_dx_b)
    call ddy_a_b_2D(mesh, basal_hydro%R, basal_hydro%dR_dy_b)
    !call checksum(mesh%pai_Tri, basal_hydro%dR_dx_b, "basal_hydro%dR_dx_b")
    !call checksum(mesh%pai_Tri, basal_hydro%dR_dy_b, "basal_hydro%dR_dy_b")

    ! Calculate u and v
    do ti = mesh%ti1, mesh%ti2
      basal_hydro%u_b( ti) = (- basal_hydro%K_b( ti) * basal_hydro%dR_dx_b( ti))
      basal_hydro%v_b( ti) = (- basal_hydro%K_b( ti) * basal_hydro%dR_dy_b( ti))
    end do
    !call checksum(mesh%pai_Tri, basal_hydro%u_b, "basal_hydro%u_b")
    !call checksum(mesh%pai_Tri, basal_hydro%v_b, "basal_hydro%v_b")

    ! Remap to c-grid velocities
    call calc_M_b_c( mesh, ice, basal_hydro)
    call map_UV_b_c( mesh, basal_hydro)
    

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



  subroutine calc_opening_rate( mesh, ice, basal_hydro)
    ! Calculate opening rate

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_ice_model),                   intent(in   )    :: ice
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'calc_opening_rate'
    integer                                               :: vi
    real(dp), parameter                                   :: c1 = 0.5 !Cavitation coefficient (m^-1)

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      ! In the Bueler and Van Pelt 2015 paper Y is defined as W, so we will do that for now too
      basal_hydro%Y( vi) = basal_hydro%W( vi)
      ! Calculate opening rate
      basal_hydro%O( vi) = c1*sqrt(basal_hydro%ice_u_base( vi)**2.0_dp + basal_hydro%ice_v_base( vi)**2.0_dp)&
                            *max((basal_hydro%W_r( vi) - basal_hydro%Y( vi)), 0.0_dp)
    end do

    !call checksum(mesh%pai_V, basal_hydro%O, "basal_hydro%O")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_opening_rate



  subroutine calc_closing_rate( mesh, ice, basal_hydro)
    ! Calculate opening rate

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_ice_model),                   intent(in   )    :: ice
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'calc_closing_rate'
    integer                                               :: vi
    real(dp), parameter                                   :: c2 = 0.04 !Creep closure coefficient 
    real(dp), parameter                                   :: A = 3.1689e-24_dp ! Ice softness parameter (Pa^-3 s^-1)

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      ! In the Bueler and Van Pelt 2015 paper Y is defined as W, so we will do that for now too
      basal_hydro%Y( vi) = basal_hydro%W( vi)
      ! Calculate closing rate
      basal_hydro%C( vi) = c2*A*(basal_hydro%P_o( vi) - basal_hydro%P( vi))**3*basal_hydro%Y( vi)
    end do
    !call checksum(mesh%pai_V, basal_hydro%C, "basal_hydro%C")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_closing_rate



  subroutine apply_W_thickness_BC_explicit( mesh, ice, basal_hydro)
    !< Apply boundary conditions to the water layer thickness on the domain border directly

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    type(type_basal_hydrology_model),       intent(inout) :: basal_hydro

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'apply_W_thickness_BC_explicit'
    integer,  dimension(mesh%vi1:mesh%vi2) :: n_interior_neighbours
    logical,  dimension(mesh%nV)           :: mask_W_tot
    real(dp),  dimension(mesh%nV)          :: W_tot
    integer                                :: vi
    character(len=256)                     :: BC_H
    integer                                :: ci,vj
    real(dp)                               :: W_sum, W_av

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate where water is present
    call calc_basal_hydro_mask_W(mesh, basal_hydro)

    call calc_n_interior_neighbours( mesh, .not. basal_hydro%mask_W, n_interior_neighbours)

    call gather_to_all(basal_hydro%mask_W, mask_W_tot)
    call gather_to_all(basal_hydro%W, W_tot)

    ! == First pass: set values of border vertices to mean of interior neighbours
    !    ...for those border vertices that actually have interior neighbours.
    ! ===========================================================================

    do vi = mesh%vi1, mesh%vi2

      if     (mesh%VBI( vi) == 1 .or. mesh%VBI( vi) == 2) then
        ! Northern domain border
        BC_H = C%BC_H_north
        !BC_H = 'zero'
      elseif (mesh%VBI( vi) == 3 .or. mesh%VBI( vi) == 4) then
        ! Eastern domain border
        BC_H = C%BC_H_east
        !BC_H = 'zero'
      elseif (mesh%VBI( vi) == 5 .or. mesh%VBI( vi) == 6) then
        ! Southern domain border
        BC_H = C%BC_H_south
        !BC_H = 'zero'
      elseif (mesh%VBI( vi) == 7 .or. mesh%VBI( vi) == 8) then
        ! Western domain border
        BC_H = C%BC_H_west
        !BC_H = 'zero'
      else
        ! Free vertex (This is the part that makes sure that not every vortex is meaned)
        cycle
      end if

      select case (BC_H)
      case default
        call crash('unknown BC_H "' // trim( BC_H) // '"')
      case ('zero')
        ! Set water thickness to zero here

        basal_hydro%W( vi) = 0._dp

      case ('infinite')
        ! Set W on this vertex equal to the average value on its neighbours

        if (n_interior_neighbours( vi) > 0) then

          W_sum = 0._dp
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            if (mesh%VBI( vj) == 0 .and. mask_W_tot( vj)) then
              W_sum = W_sum + W_tot( vj)
            end if
          end do
          W_av = W_sum / real( n_interior_neighbours( vi),dp)
          basal_hydro%W( vi) = W_av
        end if

      end select

    end do

    ! == Second pass: set values of border vertices to mean of all neighbours
    !    ...for those border vertices that have no interior neighbours.
    ! =======================================================================


    do vi = mesh%vi1, mesh%vi2

      if     (mesh%VBI( vi) == 1 .or. mesh%VBI( vi) == 2) then
        ! Northern domain border
        BC_H = C%BC_H_north
        !BC_H = 'zero'
      elseif (mesh%VBI( vi) == 3 .or. mesh%VBI( vi) == 4) then
        ! Eastern domain border
        BC_H = C%BC_H_east
        !BC_H = 'zero'
      elseif (mesh%VBI( vi) == 5 .or. mesh%VBI( vi) == 6) then
        ! Southern domain border
        BC_H = C%BC_H_south
        !BC_H = 'zero'
      elseif (mesh%VBI( vi) == 7 .or. mesh%VBI( vi) == 8) then
        ! Western domain border
        BC_H = C%BC_H_west
        !BC_H = 'zero'
      else
        ! Free vertex (This is the part that ensures that not every vortex is meaned)
        cycle
      end if

      select case (BC_H)
      case default
        call crash('unknown BC_H "' // trim( BC_H) // '"')
      case ('zero')
        ! Set water thickness to zero here

        basal_hydro%W( vi) = 0._dp

      case ('infinite')
        ! Set W on this vertex equal to the average value on its neighbours

        if (n_interior_neighbours( vi) == 0) then

          W_sum = 0._dp
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            W_sum = W_sum + W_tot( vj)
          end do
          W_av = W_sum / real( mesh%nC( vi),dp)
          basal_hydro%W( vi) = W_av

        end if

      end select

    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_W_thickness_BC_explicit



  subroutine apply_W_thickness_gl_explicit( mesh, ice, basal_hydro)
    !< Apply boundary conditions to the water layer thickness on the domain border directly

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    type(type_basal_hydrology_model),       intent(inout) :: basal_hydro

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'apply_W_thickness_gl_explicit'
    integer,  dimension(mesh%vi1:mesh%vi2) :: n_interior_neighbours
    logical,  dimension(mesh%nV)           :: mask_W_tot, mask_a_tot, mask_gl_gr_tot
    real(dp),  dimension(mesh%nV)          :: W_tot
    integer                                :: vi
    character(len=256)                     :: BC_H
    integer                                :: ci,vj
    real(dp)                               :: W_sum, W_av

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate where water is present
    call calc_basal_hydro_mask_W(mesh, basal_hydro)

    call gather_to_all(basal_hydro%mask_W, mask_W_tot)
    call gather_to_all(basal_hydro%W, W_tot)
    call gather_to_all(ice%mask_gl_gr, mask_gl_gr_tot)
    call gather_to_all(basal_hydro%mask_a, mask_a_tot)

    ! Slightly adjusted version of calc_n_interior_neighbours
    do vi = mesh%vi1, mesh%vi2

      n_interior_neighbours( vi) = 0

      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        if (mesh%VBI( vj) == 0 .and. mask_W_tot( vj) .and. .not. mask_gl_gr_tot(vj) .and. mask_a_tot( vj)) then
          n_interior_neighbours( vi) = n_interior_neighbours( vi) + 1
        end if
      end do

    end do

    ! == First pass: set values of border vertices to mean of interior neighbours
    !    ...for those border vertices that actually have interior neighbours.
    ! ===========================================================================

    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_gl_gr( vi)) then
        ! For now at least just hard-coded to infinite
        BC_H = 'infinite'

        select case (BC_H)
        case default
          call crash('unknown BC_H "' // trim( BC_H) // '"')
        case ('zero')
          ! Set water thickness to zero here

          basal_hydro%W( vi) = 0._dp

        case ('infinite')
          ! Set W on this vertex equal to the average value on its neighbours

          if (n_interior_neighbours( vi) > 0) then

            W_sum = 0._dp
            do ci = 1, mesh%nC( vi)
              vj = mesh%C( vi,ci)
              if (mesh%VBI( vj) == 0 .and. mask_W_tot( vj) .and. .not. mask_gl_gr_tot(vj) .and. mask_a_tot( vj)) then
                W_sum = W_sum + W_tot( vj)
              end if
            end do
            W_av = W_sum / real( n_interior_neighbours( vi),dp)
            basal_hydro%W( vi) = W_av
          end if

        end select
      else !ice%mask_gl_gr( vi)
        cycle
      end if

      end do

    ! == Second pass: set values of border vertices to mean of all neighbours
    !    ...for those border vertices that have no interior neighbours.
    ! =======================================================================


    do vi = mesh%vi1, mesh%vi2
      if (ice%mask_gl_gr( vi)) then
        BC_H = 'infinite'

        select case (BC_H)
        case default
          call crash('unknown BC_H "' // trim( BC_H) // '"')
        case ('zero')
          ! Set water thickness to zero here

          basal_hydro%W( vi) = 0._dp

        case ('infinite')
          ! Set W on this vertex equal to the average value on its neighbours

          if (n_interior_neighbours( vi) == 0) then

            W_sum = 0._dp
            do ci = 1, mesh%nC( vi)
              vj = mesh%C( vi,ci)
              W_sum = W_sum + W_tot( vj)
            end do
            W_av = W_sum / real( mesh%nC( vi),dp)
            basal_hydro%W( vi) = W_av

          end if

        end select
      else !ice%mask_gl_gr( vi)
        cycle
      end if

    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_W_thickness_gl_explicit



  subroutine calc_N_til(mesh, basal_hydro, W_max_til, testing_water_layer)
    ! Calculate till effective pressure

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro
    real(dp),                               intent(in   )    :: W_max_til
    logical,                                intent(in   )    :: testing_water_layer

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'calc_N_til'
    integer                                               :: vi
    real(dp), parameter                                   :: N0 = 1000_dp !Pa
    real(dp), parameter                                   :: delta = 0.02_dp !Bueler and Van Pelt 2015
    real(dp), parameter                                   :: e0 = 0.69_dp !Bueler and Van Pelt 2015
    real(dp), parameter                                   :: Cc = 0.12_dp !Bueler and Van Pelt 2015
    real(dp)                                              :: s            !W_til/W_max_til
    real(dp), parameter                                   :: g = 9.81_dp
    real(dp), parameter                                   :: rho_w = 1000.0_dp

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      s = basal_hydro%W_til( vi) / W_max_til
      basal_hydro%N_til( vi) = min(basal_hydro%P_o(vi), &
                                   N0*(delta*basal_hydro%P_o( vi)/N0)**(s)*10**(e0/Cc*(1.0_dp - s)))
      if (testing_water_layer) then     !Testing some stuff out when just adding the moving water layer pressure to effective pressure
        basal_hydro%N_til( vi) = basal_hydro%N_til( vi) - rho_w*g*basal_hydro%W( vi)
      end if
    end do
    !call checksum(mesh%pai_V, basal_hydro%N_til, "basal_hydro%N_til")
  
    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_N_til



  subroutine calc_yield_stress(mesh, ice, basal_hydro)
    ! Calculate yield stress of the till

    ! In/output variables:
    type(type_mesh),                        intent(in   )    :: mesh
    type(type_ice_model),                   intent(inout)    :: ice
    type(type_basal_hydrology_model),       intent(inout)    :: basal_hydro

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'calc_yield_stress'
    integer                                               :: vi
    real(dp), parameter                                   :: c0 = 0_dp !Pa
 

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      basal_hydro%tau_c( vi) = c0 + tan(basal_hydro%phi( vi)*pi/180._dp)*basal_hydro%N_til( vi)
      ! This next part is probably overwritten somewhere now?
      ice%till_yield_stress( vi) = basal_hydro%tau_c( vi)
    end do
    !call checksum(mesh%pai_V, basal_hydro%tau_c, "basal_hydro%tau_c")

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_yield_stress



  SUBROUTINE remap_basal_hydro_model_Salle2025( mesh_old, mesh_new, ice, basal_hydro, time)
    ! Remap the BMB model

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_basal_hydrology_model),       INTENT(INOUT) :: basal_hydro
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_basal_hydro_model_Salle2025'
    CHARACTER(LEN=256)                                    :: choice_basal_hydro_model
    integer                                               :: vi
    real(dp), parameter                                   :: rho_w = 1000.0_dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '    Remapping basal hydrology model data to the new mesh...'

    ! Reallocate memory for main variables (Figure out what to calculate again and what to remap)
    call map_from_mesh_to_mesh_with_reallocation_2D(mesh_old, mesh_new, C%output_dir, basal_hydro%W)
    call map_from_mesh_to_mesh_with_reallocation_2D(mesh_old, mesh_new, C%output_dir, basal_hydro%W_til)
    call map_from_mesh_to_mesh_with_reallocation_2D(mesh_old, mesh_new, C%output_dir, basal_hydro%P)

    call reallocate_bounds(basal_hydro%P_o, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%W_til_next, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%m, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%m_extra, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%dW_dx_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds(basal_hydro%W_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds(basal_hydro%K, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%D, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%u, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%v, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%dK_dx_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds(basal_hydro%dD_dx_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds(basal_hydro%du_dx_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds(basal_hydro%dv_dx_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds(basal_hydro%K_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds(basal_hydro%D_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds(basal_hydro%u_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds(basal_hydro%v_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds(basal_hydro%u_c, mesh_new%ei1, mesh_new%ei2)
    call reallocate_bounds(basal_hydro%v_c, mesh_new%ei1, mesh_new%ei2) !Snellius does not like this line?
    call reallocate_bounds(basal_hydro%Z, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%C, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%O, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%divQ, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%mask_a, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%mask_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds(basal_hydro%mask_W, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%R, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%dR_dx_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds(basal_hydro%dR_dy_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds(basal_hydro%W_r, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%Y, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%q_til, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%q_water_layer, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%ice_u_base, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%ice_v_base, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%ice_w_base, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%N_til, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%tau_c, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds(basal_hydro%phi, mesh_new%vi1, mesh_new%vi2)

    ! Initialise m with the same value everywhere
    do vi = mesh_new%vi1, mesh_new%vi2
      basal_hydro%m( vi) = 0.0069_dp*rho_w/sec_per_year
      basal_hydro%phi( vi) = 26.565_dp !degrees
      basal_hydro%W_r( vi) = 0.1_dp
    end do

    !write(*,*) "Doing God's work"

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_basal_hydro_model_Salle2025


END MODULE basal_hydrology_new