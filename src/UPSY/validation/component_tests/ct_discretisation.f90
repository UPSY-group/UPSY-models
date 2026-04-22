module ct_discretisation

  ! Test everything related to discretisation

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_CHAR
  use precisions, only: dp
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use ct_discretisation_mapping_derivatives, only: run_all_map_deriv_tests
  use ct_discretisation_solve_Laplace_eq, only: run_all_Laplace_eq_solving_tests
  use ct_discretisation_mapping_derivatives_graph, only: run_all_map_deriv_tests_graph

  implicit none

  private

  public :: run_all_discretisation_component_tests

contains

  !> Run all discretisation component tests.
  subroutine run_all_discretisation_component_tests( output_dir, test_mesh_filenames)

    ! In/output variables:
    character(len=*),               intent(in) :: output_dir
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_all_discretisation_component_tests'
    character(len=1024)            :: foldername_discretisation

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '  Running discretisation component tests...'
    if (par%primary) write(0,*) ''

    call run_all_map_deriv_tests         ( output_dir, test_mesh_filenames)
    call run_all_Laplace_eq_solving_tests( output_dir, test_mesh_filenames)
    ! call run_all_map_deriv_tests_graph   ( output_dir, test_mesh_filenames)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_discretisation_component_tests

end module ct_discretisation
