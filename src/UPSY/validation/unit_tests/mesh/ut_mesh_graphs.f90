module ut_mesh_graphs

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mpi_basic, only: par
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph
  use parameters, only: pi
  use mesh_memory, only: allocate_mesh_primary, crop_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use create_graphs_from_masked_mesh, only: create_graph_from_masked_mesh_a, create_graph_from_masked_mesh_b
  use ut_mesh_graphs_mapping, only: test_mesh_graph_mapping

  use netcdf_io_main

  implicit none

  private

  public :: test_graphs

contains

  subroutine test_graphs( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graphs'
    character(len=1024), parameter :: test_name_local = 'graphs'
    character(len=1024)            :: test_name
    real(dp), parameter            :: xmin = -1._dp
    real(dp), parameter            :: xmax =  1._dp
    real(dp), parameter            :: ymin = -1._dp
    real(dp), parameter            :: ymax =  1._dp
    real(dp)                       :: alpha_min
    real(dp)                       :: res_max
    type(type_mesh)                :: mesh

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call allocate_mesh_primary( mesh, 'test_mesh', 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 23.2_dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call crop_mesh_primary( mesh)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)

    call test_create_graph_from_masked_mesh( test_name, mesh)
    call test_mesh_graph_mapping( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graphs

  subroutine test_create_graph_from_masked_mesh( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_create_graph_from_masked_mesh'
    character(len=1024), parameter     :: test_name_local = 'create_graph_from_masked_mesh'
    character(len=1024)                :: test_name
    logical, dimension(:), allocatable :: mask_a
    integer                            :: nz
    integer                            :: vi
    type(type_graph)                   :: graph_a, graph_b

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define a masked set of vertices
    allocate( mask_a( mesh%vi1:mesh%vi2), source = .false.)
    do vi = mesh%vi1, mesh%vi2
      if (hypot( mesh%V( vi,1), mesh%V( vi,2)) < mesh%xmax * 0.5_dp) mask_a( vi) = .true.
    end do

    ! Create graphs from the masked vertices and triangles
    nz = 12
    call create_graph_from_masked_mesh_a( mesh, mask_a, nz, graph_a)
    call create_graph_from_masked_mesh_b( mesh, mask_a, nz, graph_b)

    call unit_test( test_graph_is_self_consistent( mesh, graph_a), trim( test_name) // '/a')
    call unit_test( test_graph_is_self_consistent( mesh, graph_b), trim( test_name) // '/b')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_create_graph_from_masked_mesh

end module ut_mesh_graphs