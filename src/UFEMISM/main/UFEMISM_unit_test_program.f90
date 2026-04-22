program UFEMISM_unit_test_program

  use petscksp
  use precisions, only: dp
  use basic_program_info, only: program_name
  use mpi_basic, only: par, initialise_parallelisation
  use parameters, only: initialise_constants
  use call_stack_and_comp_time_tracking, only: initialise_control_and_resource_tracker
  use basic_model_utilities, only: print_model_start, print_model_end
  use model_configuration, only: initialise_model_configuration_unit_tests
  use unit_tests, only: run_all_unit_tests
  use crash_mod, only: crash

  implicit none

  character(len=1024) :: foldername_output
  real(dp)            :: tstart, tstop, tcomp        !< Computation time tracking
  integer             :: ierr, perr

  program_name = 'UFEMISM'

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

  ! Get the input arguments
  if (iargc() == 1) then
    call getarg( 1, foldername_output)                 ! path/to/UPSY-models/automated_testing/UFEMISM/unit_tests/results
    if (par%primary) write(0,*) ''
    if (par%primary) write(0,*) '   Writing unit test results to  ' // trim( foldername_output) // '...'
  else
    call crash('needs input argument foldername_output')
  end if

  call initialise_model_configuration_unit_tests
  call run_all_unit_tests( foldername_output)

  ! Stop the clock
  tstop = MPI_WTIME()
  tcomp = tstop - tstart

  ! Print the UFEMISM end message to the terminal
  call print_model_end( tcomp)

  ! Finalise PETSc and MPI parallelisation
  call PetscFinalize( perr)
  call MPI_FINALIZE( ierr)

end program UFEMISM_unit_test_program
