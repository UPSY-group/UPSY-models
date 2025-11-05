program UPSY_component_test_program
  !< A program to run all of UPSY's component tests

  use basic_program_info, only: program_name
  use precisions, only: dp
  use mpi_basic, only: par
  use petscksp, only: PetscInitialize, PETSC_NULL_CHARACTER, PetscFinalize
  use mpi_basic, only: initialise_parallelisation
  use parameters, only: initialise_constants
  use control_resources_and_error_messaging, only: initialise_control_and_resource_tracker, routine_path, &
    print_model_start, print_model_end
  use mpi_f08, only: MPI_WTIME, MPI_FINALIZE

  use ct_basic, only: create_component_tests_output_folder
  use ct_create_test_meshes, only: create_all_test_meshes_and_grids
  use ct_discretisation, only: run_all_discretisation_component_tests
  use ct_remapping, only: run_all_remapping_component_tests
  use ct_mesh_focussing, only: run_all_mesh_focussing_component_tests

  implicit none

  integer                                        :: perr, ierr
  character(len=1024), parameter                 :: test_name = 'UPSY'
  character(len=1024)                            :: output_dir
  character(len=1024), dimension(:), allocatable :: test_mesh_filenames
  character(len=1024), dimension(:), allocatable :: test_grid_filenames
  real(dp)                                       :: tstart, tstop, tcomp

  program_name = 'UPSY_component_tests'

  ! Initialise MPI parallelisation and PETSc
  call initialise_parallelisation
  call PetscInitialize( PETSC_NULL_CHARACTER, perr)

  ! Initialise constants (pi, NaN, ...)
  call initialise_constants

  ! Start the clock
  tstart = MPI_WTIME()

  ! Print the model start message to the terminal
  call print_model_start

  ! Initialise the control and resource tracker
  call initialise_control_and_resource_tracker

  output_dir = 'automated_testing/component_tests/results'
  call create_component_tests_output_folder  ( output_dir)
  call create_all_test_meshes_and_grids      ( output_dir, test_mesh_filenames, test_grid_filenames)
  call run_all_mesh_focussing_component_tests( output_dir, test_mesh_filenames)
  call run_all_discretisation_component_tests( output_dir, test_mesh_filenames)
  call run_all_remapping_component_tests     ( output_dir, test_mesh_filenames, test_grid_filenames)

  ! Stop the clock
  tstop = MPI_WTIME()
  tcomp = tstop - tstart

  ! Print the UFEMISM end message to the terminal
  call print_model_end( tcomp)

  ! Finalise PETSc and MPI parallelisation
  call PetscFinalize( perr)
  call MPI_FINALIZE( ierr)

end program UPSY_component_test_program