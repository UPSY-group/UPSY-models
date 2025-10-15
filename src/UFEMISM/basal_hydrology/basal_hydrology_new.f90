MODULE basal_hydrology_new

! This is for now just a copy of the thermodynamics 3D heat equation file,
! but it should contain the new basal hydrology model in the future.
! should now figure out what I need as variables and such here.

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

  IMPLICIT NONE

CONTAINS

  subroutine basal_hydrology(mesh, ice, basal_hydro, time, dt)
    ! In/output variables:
    type(type_mesh),                  intent(in   ) :: mesh
    type(type_ice_model),             intent(inout) :: ice
    type(type_basal_hydrology_model), intent(inout) :: basal_hydro
    real(dp),                         intent(in)    :: time
    real(dp),                         intent(in)    :: dt

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'basal_hydrology'
    integer                        :: vi
    logical                        :: sliding             ! Whether the ice slides or not when W = 0 and there is ice on land

    real(dp), parameter            :: W_max = 1000000     ! Maximum basal water depth (number is placeholder for now, because there is no maximum W)
    real(dp), parameter            :: W_min = 0           ! Minimum basal water depth

    real(dp), parameter            :: W_max_til = 2       ! Maximum basal water depth in till
    real(dp), parameter            :: W_min_til = 0       ! Minimum basal water depth in till

    real(dp), parameter            :: P_min = 0           ! Minimum pressure of the ice on the basal water

    ! real(dp), parameter            :: rho_w               ! Density of water
    real(dp), parameter            :: rho_i = 917.0_dp    ! Density of ice

    ! real(dp), parameter            :: k_coef              ! Coefficient of effective conductivity
    ! real(dp), parameter            :: alpha               ! Exponent used in effective conductivity
    ! real(dp), parameter            :: beta                ! Exponent used in effective conductivity

    ! real(dp), parameter            :: Cd                  ! Gradual drain of water in till

    real(dp), parameter            :: g = 9.81_dp         ! Gravitational acceleration

    ! real(dp), parameter            :: c_1                 ! Scaling coefficient 1 (non-negative opening term)
    ! real(dp), parameter            :: c_2                 ! Scaling coefficient 2 (closing term)

    ! real(dp), parameter            :: W_r                 ! Maximum roughness scale of basal topography

    ! Add routine to path
    call init_routine( routine_name)

    ! Here we actually do some stuff

    write(*,*) "Mesh vi1: ", mesh%vi1
    write(*,*) "Mesh vi2: ", mesh%vi2
    ! 1) Start with W, W_til and P and make sure they are all within their bounds
    ! Make the arrays for the variables

    ! Should make a function to initialise these arrays and all the other allocatables at once
    allocate(basal_hydro%P_o(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%W(mesh%vi1:mesh%vi2), source = 0.0_dp)
    allocate(basal_hydro%W_til(mesh%vi1:mesh%vi2), source =  0.0_dp)
    allocate(basal_hydro%P(mesh%vi1:mesh%vi2), source = 0.0_dp)

    
    ! Loop over all vertices and set them within their bounds
    do vi = mesh%vi1, mesh%vi2
      basal_hydro%P_o( vi)     = rho_i * g * ice%Hi( vi) ! Calculate overburden pressure
      basal_hydro%W( vi)       = min( max( basal_hydro%W( vi),       W_min),       W_max)
      basal_hydro%W_til( vi)   = min( max( basal_hydro%W_til( vi),   W_min_til),   W_max_til)
      basal_hydro%P( vi)       = min( max( basal_hydro%P( vi),       P_min),       basal_hydro%P_o( vi))
    end do

    call save_variable_as_netcdf_dp_1D(C%output_dir, basal_hydro%W, "W")

    write(*,*) "Overburden pressure for first few vertices: ", basal_hydro%P_o(mesh%vi1:mesh%vi1+5)
    write(*,*) "Water thickness for first few vertices: ", basal_hydro%W(mesh%vi1:mesh%vi1+5)
    write(*,*) "Till water thickness for first few vertices: ", basal_hydro%W_til(mesh%vi1:mesh%vi1+5)
    write(*,*) "Pressure for first few vertices: ", basal_hydro%P(mesh%vi1:mesh%vi1+5)


    ! 2) Perform a timestep to get W_til one timestep further, still making sure it is within the bounds
    ! basal_hydro%W_til = basal_hydro%W_til + dt*(basal_hydro%m/rho_w - Cd) ! Timestep

    ! do vi = mesh%vi1, mesh%vi2 ! Make sure within bounds
    !   basal_hydro%W_til( vi)   = min( max( basal_hydro%W_til( vi),   W_min_til),   W_max_til)
    ! end do
    

    ! 3) If icefree or floating, set W_til to zero (do loops in 2 and 3 can be combined)
    ! do vi = mesh%vi1, mesh%vi2
    !   if (ice%mask_icefree_land( vi) .or. ice%mask_floating_ice( vi)) then
    !     basal_hydro%W_til( vi) = 0.0_dp
    !   end if
    ! end do


    ! 4) Get W values on staggered grid (probably already a function for in this code somewhere)
    ! Does this have to do with something similar like update_laddie_operaters subroutine?
    ! Look in mesh_types.f90 to see the matrices?
    ! Maybe it is just the ti instead of vi? Look at laddie_velocity.f90? laddie_model_types.f90?
    ! laddie_utilities.f90?
    ! basal_hydro%W_stag = mesh%M_ddx_a_b * basal_hydro%W
    ! call ddx_a_b_2D(mesh, basal_hydro%W, basal_hydro%dW_dx_b)

    ! 5) Also do this for effective conductivity K
    ! basal_hydro%K_stag = mesh%M_ddx_a_b * basal_hydro%K


    ! 6) Do this again for velocity components u and v
    ! basal_hydro%u_stag = mesh%M_ddx_a_b * basal_hydro%u
    ! basal_hydro%v_stag = mesh%M_ddy_a_b * basal_hydro%v


    ! 7) And lastly for the diffusivity D
    ! basal_hydro%D_stag = mesh%M_ddx_a_b * basal_hydro%D


    ! 8) Get the timestep (from UFEMISM or use the definition from paper?)
    ! The definitions in the paper use dx and dy, but that does not really exist here.
    !dt = min(dt_CFL, dt_W, dt_P)

    ! 9) Compute the advective fluxes (Q) on the staggered grid (also using a function already written in this code?)
    
    ! This is very similar to what is done in laddie_thickness.90 in the compute_divQH subroutine
    ! Initialise with zeros
    ! basal_hydro%divQH( mesh%vi1:mesh%vi2) = 0.0_dp

    

    ! 10) Compute the flux divergence approximations D_ij


    ! 11) If icefree set next timestep of P to 0, if floating set to overburden pressure
    ! 11) If W at this timestep is 0 and if icefree and floating are both false, set next timestep of P to 0 (any sliding) or overburden pressure (no sliding)
    ! 11) Otherwise, compute next timestep of P using the equation in the paper
    ! do vi = mesh%vi1, mesh%vi2
    !   if (ice%mask_icefree_land( vi)) then
    !     basal_hydro%P( vi) =  0.0_dp
    !   else if (ice%mask_floating_ice( vi)) then
    !     basal_hydro%P( vi) = basal_hydro%P_o( vi)
    !   else if (basal_hydro%W == 0.0_dp .and. .not. ice%mask_icefree_land( vi) .and. .not. ice%mask_floating_ice( vi)) then
    !     if (sliding) then
    !       basal_hydro%P( vi) = 0.0_dp
    !     else
    !       basal_hydro%P( vi) = basal_hydro%P_o( vi)
    !     end if
    !   else
        ! Compute next timestep of P using the equation in the paper
        ! (to be implemented)
      

    ! 12) Make sure P is within its bounds
    !     basal_hydro%P( vi) = min( max( basal_hydro%P( vi), basal_hydro%P_min( vi)), basal_hydro%P_o( vi))
    !   end if
    ! end do



    ! 13) If icefree or float, then set next timestep of W to 0.
    ! 13) Otherwise, compute next timestep of W using the equation in the paper
    ! do vi = mesh%vi1, mesh%vi2
    !   if (ice%mask_icefree_land( vi) .or. ice%mask_floating_ice( vi)) then
    !     basal_hydro%W( vi) = 0.0_dp
    !   else
        ! Compute next timestep of W using the equation in the paper
        ! (to be implemented)

    ! 14) Make sure W is within its bounds (>= 0)
    !   end if
    !   basal_hydro%W( vi) = min( max( basal_hydro%W( vi), basal_hydro%W_min( vi)), basal_hydro%W_max( vi))
    ! end do


    ! 15) Update time and repeat
    ! time = time + dt

    call crash('Hello world!')
  
    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine basal_hydrology

END MODULE basal_hydrology_new