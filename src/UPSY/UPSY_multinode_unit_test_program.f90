program UPSY_multinode_unit_test_program
  !< A program to run all of UPSY's multi-node unit tests

  use basic_program_info, only: program_name
  use precisions, only: dp
  use petscksp, only: PetscInitialize, PETSC_NULL_CHARACTER, PetscFinalize
  use mpi_basic, only: initialise_parallelisation_multinode_tests
  use parameters, only: initialise_constants
  use control_resources_and_error_messaging, only: initialise_control_and_resource_tracker, &
    print_model_start, print_model_end, crash
  use ut_basic, only: foldername_unit_tests_output, filename_unit_tests_output
  use mpi_f08, only: MPI_WTIME, MPI_FINALIZE

  use ut_mpi_dist_shared_memory, only: unit_tests_mpi_hybrid_distributed_shared_memory_main
  use ut_halo_exchange, only: test_halo_exchange_main
  use ut_halo_exchange_mesh, only: test_mesh_halo_exchange_main
  use ut_mpi_CSR_matrix_algebra, only: test_CSR_matrix_algebra_main

  implicit none

  integer                        :: perr, ierr
  character(len=1024), parameter :: test_name = 'UPSY'
  real(dp)                       :: tstart, tstop, tcomp
  logical                        :: ex

  program_name = 'UPSY_multinode_unit_tests'

  ! Initialise MPI parallelisation and PETSc
  call initialise_parallelisation_multinode_tests
  call PetscInitialize( PETSC_NULL_CHARACTER, perr)

  ! Initialise constants (pi, NaN, ...)
  call initialise_constants

  ! Start the clock
  tstart = MPI_WTIME()

  ! Print the model start message to the terminal
  call print_model_start

  ! Initialise the control and resource tracker
  call initialise_control_and_resource_tracker

  ! Check if a unit tests output file exists (should have been created by UPSY_unit_test_program)
  foldername_unit_tests_output = 'automated_testing/unit_tests/results'
  filename_unit_tests_output = trim( foldername_unit_tests_output) // '/unit_tests_output.txt'
  inquire( file = trim( filename_unit_tests_output), exist = ex)
  if (.not. ex) call crash('Unit tests output file not found - run UPSY_unit_test_program first')

  ! Run all the multi-node unit tests
  call unit_tests_mpi_hybrid_distributed_shared_memory_main( test_name)
  call test_halo_exchange_main( test_name)
  call test_mesh_halo_exchange_main( test_name)
  call test_CSR_matrix_algebra_main( test_name)

  ! Stop the clock
  tstop = MPI_WTIME()
  tcomp = tstop - tstart

  ! Print the UFEMISM end message to the terminal
  call print_model_end( tcomp)

  ! Finalise PETSc and MPI parallelisation
  call PetscFinalize( perr)
  call MPI_FINALIZE( ierr)

end program UPSY_multinode_unit_test_program