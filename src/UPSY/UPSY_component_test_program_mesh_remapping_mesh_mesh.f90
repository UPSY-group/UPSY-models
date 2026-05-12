program UPSY_component_test_program_mesh_remapping_mesh_mesh

  use basic_program_info, only: program_name
  use precisions, only: dp
  use mpi_basic, only: par
  use petscksp, only: PetscInitialize, PETSC_NULL_CHARACTER, PetscFinalize
  use mpi_basic, only: initialise_parallelisation
  use parameters, only: initialise_constants
  use call_stack_and_comp_time_tracking, only: initialise_control_and_resource_tracker, routine_path
  use basic_model_utilities, only: print_model_start, print_model_end
  use mpi_f08, only: MPI_WTIME, MPI_FINALIZE
  use crash_mod, only: crash
  use git_commit_hash_and_package_versions, only: print_git_commit_hash_and_package_versions

  use ct_create_test_meshes, only: list_test_meshes_and_grids_in_folder
  use ct_basic, only: create_component_tests_output_folder
  use ct_remapping_mesh_to_mesh, only: run_all_mesh_to_mesh_remapping_tests

  implicit none

  integer                                        :: perr, ierr
  character(len=1024)                            :: foldername_test_meshes_and_grids, foldername_output
  character(len=1024), dimension(:), allocatable :: test_mesh_filenames
  character(len=1024), dimension(:), allocatable :: test_grid_filenames
  real(dp)                                       :: tstart, tstop, tcomp

  program_name = 'UPSY_component_test_mesh_remapping_mesh_mesh'

  ! Initialise MPI parallelisation and PETSc
  call initialise_parallelisation
  call PetscInitialize( PETSC_NULL_CHARACTER, perr)

  call print_git_commit_hash_and_package_versions

  ! Initialise constants (pi, NaN, ...)
  call initialise_constants

  ! Start the clock
  tstart = MPI_WTIME()

  ! Print the model start message to the terminal
  call print_model_start

  ! Initialise the control and resource tracker
  call initialise_control_and_resource_tracker

  ! Get the input arguments
  if (iargc() == 2) then
    call getarg( 1, foldername_test_meshes_and_grids)  ! path/to/UPSY-models/automated_testing/test_meshes_and_grids
    call getarg( 2, foldername_output)                 ! path/to/UPSY-models/automated_testing/component_test_mesh_remapping_mesh_mesh/results
    if (par%primary) write(0,*) ''
    if (par%primary) write(0,*) '   Reading test meshes and grids from ' // trim( foldername_test_meshes_and_grids) // '...'
    if (par%primary) write(0,*) '   Writing component test results to  ' // trim( foldername_output) // '...'
  else
    call crash('needs input arguments foldername_test_meshes_and_grids and foldername_output')
  end if

  call list_test_meshes_and_grids_in_folder( foldername_test_meshes_and_grids, test_grid_filenames, test_mesh_filenames)
  call create_component_tests_output_folder( foldername_output)
  call run_all_mesh_to_mesh_remapping_tests( foldername_output, test_mesh_filenames)

  ! Stop the clock
  tstop = MPI_WTIME()
  tcomp = tstop - tstart

  ! Print the UFEMISM end message to the terminal
  call print_model_end( tcomp)

  ! Finalise PETSc and MPI parallelisation
  call PetscFinalize( perr)
  call MPI_FINALIZE( ierr)

end program UPSY_component_test_program_mesh_remapping_mesh_mesh