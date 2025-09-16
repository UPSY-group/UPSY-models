program LADDIE_program
  !
  ! ===============================================================================
  ! = The main program of the one Layer Antarctic model for Dynamical Downscaling =
  ! = of Ice-ocean Exchanges (LADDIE)                                             =
  ! =                                                                             =
  ! = Main developers:                                                            =
  ! =                                                                             =
  ! = dr. E. (Erwin) Lambert                                                      =
  ! =   Affiliation: Royal Dutch Meteorological Institute (KNMI)                  =
  ! =   E-mail     : erwin dot lambert at knmi dot nl                             =
  ! =                                                                             =
  ! = dr. C. J. (Tijn) Berends                                                    =
  ! =   Affiliation: Institute for Marine and Atmospheric Research Utrecht (IMAU) =
  ! =   E-mail     : c.j.berends@uu.nl                                            =
  ! =                                                                             =
  ! ===============================================================================
  !
  ! NOTE: the executable should be run using mpiexec to specify the number of
  !       cores n, and with the path to the config file as the only argument, e.g.:
  !
  !       mpirun -n 2 LADDIE_program  config-files/config_test

! ===== Preamble =====
! ====================

  use petscksp
  use precisions, only: dp
  use basic_program_info, only: program_name
  use mpi_basic, only: par, initialise_parallelisation
  use parameters, only: initialise_constants
  use control_resources_and_error_messaging, only: warning, crash, happy, init_routine, finalise_routine, &
    colour_string, do_colour_strings, initialise_control_and_resource_tracker, reset_resource_tracker, &
    print_model_start, print_model_end
  use model_configuration, only: C, initialise_model_configuration, initialise_model_configuration_unit_tests
  use netcdf_io_main
  use mesh_types, only: type_mesh
  use laddie_model_types, only: type_laddie_model
  use laddie_forcing_types, only: type_laddie_forcing
  use reference_geometry_types, only: type_reference_geometry
  use laddie_forcing_main, only: initialise_forcing
  use LADDIE_main_model, only: run_laddie_model, initialise_laddie_model
  use laddie_unit_tests, only: run_laddie_unit_tests

  implicit none

! ===== Main variables =====
! ==========================

  ! The LADDIE forcing
  type(type_laddie_forcing)              :: forcing

  ! The LADDIE model
  type(type_laddie_model)                :: laddie
  type(type_mesh)                        :: mesh

  ! Computation time tracking
  real(dp)                               :: tstart, tstop, tcomp

  ! Input argument
  character(len=1024)                    :: input_argument

  integer :: ierr, perr
  real(dp), parameter                    :: time = 0.0_dp
  logical, parameter                     :: is_initial = .false.
  logical, parameter                     :: is_standalone = .true.

! ===== START =====
! =================

  program_name = 'LADDIE'

  ! Get the input argument (either the path to the config file,
  ! or an instruction to run unit/component tests)
  if (iargc() == 1) then
    call getarg( 1, input_argument)
  else
    stop 'LADDIE requires a single argument, being the path to the config file, e.g. "mpirun -n 2 LADDIE_program  config-files/config_test"'
  end if

  ! Initialise MPI parallelisation and PETSc
  call initialise_parallelisation( input_argument)
  call PetscInitialize( PETSC_NULL_CHARACTER, perr)

  ! Initialise constants (pi, NaN, ...)
  call initialise_constants

  ! Only the primary process "sees" the input argument; all the others are
  ! initialised by MPI without it. Broadcast it so they know what to do.
  call MPI_BCAST( input_argument, len(input_argument), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

  ! Start the clock
  tstart = MPI_WTIME()

  ! Print the LADDIE start message to the terminal
  call print_model_start

  ! Initialise the control and resource tracker
  call initialise_control_and_resource_tracker

  ! Special cases
  if (input_argument == 'unit_tests') then
    call initialise_model_configuration_unit_tests
    call run_laddie_unit_tests
  else ! An actual model simulation

    ! Initialise the main model configuration
    call initialise_model_configuration

    ! Create the resource tracking output file
    call create_resource_tracking_file( C%output_dir)

    ! == Initialise forcing and mesh ==
    ! ==================================

    call initialise_forcing( mesh, forcing)

    ! == Initialise the model ==
    ! ==========================

    call initialise_laddie_model( mesh, laddie, forcing, is_standalone)

    ! == Run the model ==
    ! ===================

    call run_laddie_model( mesh, laddie, forcing, time, is_initial, is_standalone)

    ! Write to resource tracking file
    call write_to_resource_tracking_file( time)
    call reset_resource_tracker

  end if ! do_unit_test/do_benchmark/run

! ===== FINISH =====
! ==================

  ! Stop the clock
  tstop = MPI_WTIME()
  tcomp = tstop - tstart

  ! Print the LADDIE end message to the terminal
  call print_model_end( tcomp)

  ! Finalise PETSc and MPI parallelisation
  call PetscFinalize( perr)
  call MPI_FINALIZE( ierr)


end program LADDIE_program
