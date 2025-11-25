module ut_mesh_duplicate

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mpi_basic, only: par
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary, crop_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use parameters, only: pi
  use mesh_memory, only: duplicate_mesh_primary
  use mesh_utilities, only: check_if_meshes_are_identical

  use netcdf_io_main

  implicit none

  private

  public :: test_duplicate_mesh

contains

  subroutine test_duplicate_mesh( test_name_parent)
    ! Test the trace_line_grid_a subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_duplicate_mesh'
    character(len=1024), parameter :: test_name_local = 'duplicate_mesh'
    character(len=1024)            :: test_name
    real(dp), parameter            :: xmin = 0._dp
    real(dp), parameter            :: xmax = 1._dp
    real(dp), parameter            :: ymin = 0._dp
    real(dp), parameter            :: ymax = 1._dp
    real(dp)                       :: alpha_min
    real(dp)                       :: res_max
    type(type_mesh)                :: mesh, mesh2
    logical                        :: are_identical

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call allocate_mesh_primary( mesh, 'test_mesh', 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 20._dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call crop_mesh_primary( mesh)

    ! Duplicate the mesh
    call duplicate_mesh_primary( mesh, mesh2)

    ! Verify that it worked
    call check_if_meshes_are_identical( mesh, mesh2, are_identical)
    call unit_test( are_identical, trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_duplicate_mesh

end module ut_mesh_duplicate