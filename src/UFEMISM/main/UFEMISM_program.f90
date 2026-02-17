program UFEMISM_program

  ! ===============================================================================
  ! = The main program of the Utrecht FinitE voluMe Ice Sheet Model (UFEMISM)     =
  ! =                                                                             =
  ! = Main developer:                                                            =
  ! =                                                                             =
  ! = dr. C. J. (Tijn) Berends                                                    =
  ! =   Affiliation: Institute for Marine and Atmospheric Research Utrecht (IMAU) =
  ! =   E-mail     : c.j.berends@uu.nl                                            =
  ! =                                                                             =
  ! ===============================================================================

  use petscksp
  use precisions, only: dp
  use basic_program_info, only: program_name
  use mpi_basic, only: par, initialise_parallelisation
  use parameters, only: initialise_constants
  use call_stack_and_comp_time_tracking, only: initialise_control_and_resource_tracker, &
    reset_resource_tracker
  use basic_model_utilities, only: print_model_start, print_model_end
  use model_configuration, only: C, initialise_model_configuration, initialise_model_configuration_unit_tests
  use netcdf_io_main
  use region_types, only: type_model_region
  use global_forcing_types, only: type_global_forcing
  use UFEMISM_main_model, only: initialise_model_region, run_model_region
  use global_forcings_main, only: initialise_global_forcings, update_global_forcings
  use inversion_utilities, only: MISMIPplus_adapt_flow_factor
  use unit_tests, only: run_all_unit_tests
  use component_tests, only: run_all_component_tests
  use checksum_mod, only: create_checksum_logfile
  use ice_geometry_calculations_mod, only: type_ice_geometry

  implicit none

  character(len=1024)       :: input_argument
  type(type_model_region)   :: NAM, EAS, GRL, ANT          !< The four model regions
  type(type_ice_geometry)   :: geom                       !< The ice geometry 
  type(type_global_forcing) :: forcing                     !< The global forcings
  real(dp)                  :: t_coupling, t_end_models    !< Coupling times
  real(dp)                  :: tstart, tstop, tcomp        !< Computation time tracking
  real(dp)                  :: Hs_prev, Hs_cur             !< Surface elevations for the automated flow factor tuning in MISMIP+
  integer                   :: ierr, perr

  program_name = 'UFEMISM'

  ! Get the input argument (either the path to the config file,
  ! or an instruction to run unit/component tests)
  if (iargc() == 1) then
    call getarg( 1, input_argument)
  else
    stop 'UFEMISM requires a single argument, being the path to the config file, e.g. "mpi_exec  -n 2  UFEMISM_program  config-files/config_test"'
  end if

  ! Initialise MPI parallelisation and PETSc
  call initialise_parallelisation
  call PetscInitialize( PETSC_NULL_CHARACTER, perr)

  ! Initialise constants (pi, NaN, ...)
  call initialise_constants

  ! Start the clock
  tstart = MPI_WTIME()

  ! Print the UFEMISM start message to the terminal
  call print_model_start

  ! Initialise the control and resource tracker
  call initialise_control_and_resource_tracker

  ! Special cases
  if (input_argument == 'unit_tests') then
    call initialise_model_configuration_unit_tests
    call run_all_unit_tests
  elseif (input_argument == 'component_tests') then
    call initialise_model_configuration_unit_tests
    call run_all_component_tests
  else ! An actual model simulation

    ! == Initialisation ==
    ! ====================

    ! Initialise the main model configuration
    call initialise_model_configuration( input_argument)

    ! Create the resource tracking output file
    call create_resource_tracking_file( C%output_dir)
    call create_checksum_logfile( C%output_dir)

    ! Initialise surface elevations for the automated flow factor tuning in MISMIP+
    Hs_cur = 1._dp

    ! Initialise global forcing timeseries
    call initialise_global_forcings(forcing)

    ! Initialise the model regions
    if (C%do_NAM) call initialise_model_region( NAM, 'NAM', forcing, C%start_time_of_run, geom)
    if (C%do_EAS) call initialise_model_region( EAS, 'EAS', forcing, C%start_time_of_run, geom)
    if (C%do_GRL) call initialise_model_region( GRL, 'GRL', forcing, C%start_time_of_run, geom)
    if (C%do_ANT) call initialise_model_region( ANT, 'ANT', forcing, C%start_time_of_run, geom)

    ! == The coupling time loop ==
    ! ============================

    t_coupling = C%start_time_of_run

    do while (t_coupling < C%end_time_of_run)

      ! Run all model regions forward in time for one coupling interval
      t_end_models = MIN( C%end_time_of_run, t_coupling + C%dt_coupling)

      call update_global_forcings(forcing, t_coupling)

      if (C%do_NAM) call run_model_region( NAM, t_end_models, forcing, geom)
      if (C%do_EAS) call run_model_region( EAS, t_end_models, forcing, geom)
      if (C%do_GRL) call run_model_region( GRL, t_end_models, forcing, geom)
      if (C%do_ANT) call run_model_region( ANT, t_end_models, forcing, geom)

      ! Advance coupling time
      t_coupling = t_end_models

      ! MISMIP+ flow factor tuning for GL position
      if (C%refgeo_idealised_MISMIPplus_tune_A) then
        Hs_prev = Hs_cur
        Hs_cur  = MAXVAL( ANT%ice%Hs)
        call MPI_ALLREDUCE( MPI_IN_PLACE, Hs_cur, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
        if (ABS( 1._dp - Hs_cur / Hs_prev) < 5.0E-3_dp) then
          ! The model has converged to a steady state; adapt the flow factor
          call MISMIPplus_adapt_flow_factor( ANT%mesh, ANT%ice)
        end if
      end if

      ! Write to resource tracking file
      call write_to_resource_tracking_file( t_coupling)
      call reset_resource_tracker

    end do

  end if

  ! Stop the clock
  tstop = MPI_WTIME()
  tcomp = tstop - tstart

  ! Print the UFEMISM end message to the terminal
  call print_model_end( tcomp)

  ! Finalise PETSc and MPI parallelisation
  call PetscFinalize( perr)
  call MPI_FINALIZE( ierr)

end program UFEMISM_program
