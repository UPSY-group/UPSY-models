program UPSY_component_test_program_mesh_creation

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

  use ct_create_test_meshes, only: create_all_test_meshes_and_grids

  implicit none

  integer             :: perr, ierr
  character(len=1024) :: foldername_output
  real(dp)            :: tstart, tstop, tcomp

  program_name = 'UPSY_component_test_mesh_creation'

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

  ! Get the input arguments
  if (iargc() == 1) then
    call getarg( 1, foldername_output)  ! path/to/UPSY-models/automated_testing/UPSY/component_test_mesh_creation/results
    if (par%primary) write(0,*) ''
    if (par%primary) write(0,*) '   Writing component test results to  ' // trim( foldername_output) // '...'
  else
    call crash('needs input argument foldername_output')
  end if

  call create_all_test_meshes_and_grids( foldername_output)

  ! Stop the clock
  tstop = MPI_WTIME()
  tcomp = tstop - tstart

  ! Print the program end message to the terminal
  call print_model_end( tcomp)

  ! Finalise PETSc and MPI parallelisation
  call PetscFinalize( perr)
  call MPI_FINALIZE( ierr)

end program UPSY_component_test_program_mesh_creation