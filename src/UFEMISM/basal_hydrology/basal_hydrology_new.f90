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

  subroutine basal_hydrology(mesh, ice)
    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'basal_hydrology'
    ! Might need some more variables here?

    ! Add routine to path
    call init_routine( routine_name)

    ! Here we actually do some stuff
    ! =========================
    ! =========================
    ! =========================
  
    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine basal_hydrology

END MODULE basal_hydrology_new