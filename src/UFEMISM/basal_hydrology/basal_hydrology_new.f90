MODULE basal_hydrology_new

! ===== Preamble =====
! ====================
! This is preamble is still copied from thermodynamics_3D_heat_equation.f90

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par
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
  USE mesh_disc_apply_operators                              , ONLY: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, map_b_a_2D
  use laddie_utilities                                       , ONLY: map_H_a_c
  use mpi_distributed_memory                                 , only: gather_to_all
  use mesh_halo_exchange                                     , only: exchange_halos
  use CSR_matrix_vector_multiplication                       , only: multiply_CSR_matrix_with_vector_1D

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

    real(dp), parameter            :: Cd = 0.001_dp        ! Gradual drain of water in till (m a^-1)

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

    ! 1) Start with W, W_til and P and make sure they are all within their bounds

    ! All kinds of variables are allocated down below. This is called in UFEMISM_main_model around line 500.

    ! Loop over all vertices and set them within their bounds
    call set_within_bounds(mesh, ice, basal_hydro, W_min, W_max, W_min_til, W_max_til, P_min)

    ! 2) Perform a timestep to get W_til one timestep further, still making sure it is within the bounds

    call calc_W_til_next(mesh, ice, basal_hydro, W_min_til, W_max_til, dt)

    ! 4) Get W values on staggered grid 
    call map_all_a_b(mesh, basal_hydro)

    ! 8) Get the timestep 
    ! Mainly inspired by calc_critical_timestep_SIA subroutine in time_step_criteria.f90

    call get_basal_hydro_timestep(mesh, basal_hydro, dt, dt_hydro)

    write(*,*) "dt_hydro = ", dt_hydro

    ! 9) Compute the advective fluxes (Q) on the staggered grid
    
    call calc_divQ(mesh, ice, basal_hydro)

    ! 11) If icefree set next timestep of P to 0, if floating set to overburden pressure
    ! 11) If W at this timestep is 0 and if icefree and floating are both false, set next timestep of P to 0 (any sliding) or overburden pressure (no sliding)
    ! 11) Otherwise, compute next timestep of P using the equation in the paper
    call calc_P_next(mesh, ice, basal_hydro, P_min, dt_hydro)

    ! 13) If icefree or float, then set next timestep of W to 0.
    ! 13) Otherwise, compute next timestep of W using the equation in the paper
    call calc_W_next(mesh, ice, basal_hydro, W_min, W_max, dt_hydro)

    ! 15) Update time and repeat (not how it is done here?)
    ! time = time + dt_hydro

    write(*,*) "Time after basal hydrology step: ", time

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

    ! Add routine to path
    call init_routine( routine_name)

    allocate(basal_hydro%P_o(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%W(mesh%vi1:mesh%vi2), source = 1.0_dp)
    allocate(basal_hydro%W_til(mesh%vi1:mesh%vi2), source =  1.0_dp)
    allocate(basal_hydro%W_til_next(mesh%vi1:mesh%vi2), source =  0.0_dp)
    allocate(basal_hydro%P(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%m(mesh%vi1:mesh%vi2), source = 0.01_dp) ! For now just make this a value
    allocate(basal_hydro%dW_dx_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%W_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%K(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%D(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%u(mesh%vi1:mesh%vi2), source = 2.0_dp)
    allocate(basal_hydro%v(mesh%vi1:mesh%vi2), source = 1.0_dp)
    allocate(basal_hydro%dK_dx_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%dD_dx_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%du_dx_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%dv_dx_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%K_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%D_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%u_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%v_b(mesh%ti1:mesh%ti2), source = 0.0_dp)
    allocate(basal_hydro%u_c(mesh%ei1:mesh%ei2), source = 1.0_dp)
    allocate(basal_hydro%v_c(mesh%ei1:mesh%ei2), source = 0.0_dp)
    allocate(basal_hydro%Z(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%C(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%O(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%divQ( mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%old_time, source = 0.0_dp)

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
    real(dp), parameter :: Cd = 0.001_dp/sec_per_year
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
    !call map_H_a_c(mesh, basal_hydro%u, basal_hydro%u_c) ! Is there a map function for this not in LADDIE?
    !call map_H_a_c(mesh, basal_hydro%v, basal_hydro%v_c)

    ! 7) And lastly for the diffusivity D
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

      ! Diffusivity on the triangle
      D_t = basal_hydro%D_b( ti)

      ! Instead of max of D, we use D on the triangle (I think this makes sense?)
      dt_crit_W = min(dt_crit_W, d_min**2/(8*D_t))

      dt_crit_P = min(dt_crit_P, 2*phi*d_min**2/(8*D_t))

    end do

    ! Timestep we will use here
    dt_hydro = min(correction_factor*min( dt_crit_CFL, dt_crit_W, dt_crit_P), dt)

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

      ! Initialise

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
          ! Skip connection if neighbour is not grounded or there is no ice. No flux across grounding line
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
    real(dp), parameter :: Cd = 0.001_dp/sec_per_year
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
    logical, parameter  :: sliding = .true.         ! Whether the ice slides or not when W = 0 and there is ice on land


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

  call multiply_CSR_matrix_with_vector_1D(basal_hydro%M_b_c, &
    mesh%pai_Tri, basal_hydro%u, mesh%pai_E, basal_hydro%u_c)
  call multiply_CSR_matrix_with_vector_1D( basal_hydro%M_b_c, &
    mesh%pai_Tri, basal_hydro%v, mesh%pai_E, basal_hydro%v_c)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine map_UV_b_c

END MODULE basal_hydrology_new