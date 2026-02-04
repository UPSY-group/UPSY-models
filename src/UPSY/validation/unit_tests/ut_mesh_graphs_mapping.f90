module ut_mesh_graphs_mapping

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use mpi_basic, only: par
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph
  use create_graphs_from_masked_mesh, only: create_graph_from_masked_mesh_a, create_graph_from_masked_mesh_b
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD
  use mpi_distributed_memory, only: gather_to_all
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use deallocate_dist_shared_mod, only: deallocate_dist_shared
  use mesh_graph_mapping, only: map_mesh_vertices_to_graph, map_graph_to_mesh_vertices, &
    map_mesh_triangles_to_graph, map_graph_to_mesh_triangles, &
    map_mesh_edges_to_graph, map_graph_to_mesh_edges

  implicit none

  private

  public :: test_mesh_graph_mapping

contains

  subroutine test_mesh_graph_mapping( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_mesh_graph_mapping'
    character(len=1024), parameter     :: test_name_local = 'mapping'
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

    call test_mesh_vertices_to_graph_mapping ( test_name, mesh, graph_a)
    call test_mesh_triangles_to_graph_mapping( test_name, mesh, graph_b)
    call test_mesh_edges_to_graph_mapping    ( test_name, mesh, graph_b)
    call test_graph_to_mesh_vertices_mapping ( test_name, mesh, graph_a)
    call test_graph_to_mesh_triangles_mapping( test_name, mesh, graph_b)
    call test_graph_to_mesh_edges_mapping    ( test_name, mesh, graph_b)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_graph_mapping

  ! Mesh vertices to graph
  ! ======================

  subroutine test_mesh_vertices_to_graph_mapping( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_mesh_vertices_to_graph_mapping'
    character(len=1024), parameter :: test_name_local = 'mesh_vertices_to_graph'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_mesh_vertices_to_graph_mapping_logical_2D( test_name, mesh, graph)
    call test_mesh_vertices_to_graph_mapping_logical_3D( test_name, mesh, graph)
    call test_mesh_vertices_to_graph_mapping_int_2D    ( test_name, mesh, graph)
    call test_mesh_vertices_to_graph_mapping_int_3D    ( test_name, mesh, graph)
    call test_mesh_vertices_to_graph_mapping_dp_2D     ( test_name, mesh, graph)
    call test_mesh_vertices_to_graph_mapping_dp_3D     ( test_name, mesh, graph)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping

  subroutine test_mesh_vertices_to_graph_mapping_logical_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_mesh_vertices_to_graph_mapping_logical_2D'
    character(len=1024), parameter             :: test_name_local = 'logical_2D'
    character(len=1024)                        :: test_name
    logical, dimension(:), allocatable         :: d_mesh_loc_ex
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    logical, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    logical, dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                    :: vi, ni
    logical                                    :: d
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex( mesh%vi1 :mesh%vi2), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_V%n_nih)
    d_mesh_nih_ex( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        d = modulo( ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))),2) == 0
        d_mesh_loc_ex( vi) = d
        d_mesh_nih_ex( vi) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        d = modulo( ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))),2) == 0
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) .eqv. d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) .eqv. d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping_logical_2D

  subroutine test_mesh_vertices_to_graph_mapping_logical_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_mesh_vertices_to_graph_mapping_logical_3D'
    character(len=1024), parameter               :: test_name_local = 'logical_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    logical, dimension(:,:), allocatable         :: d_mesh_loc_ex
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                      :: vi, ni, k
    logical                                      :: d
    logical                                      :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%vi1 :mesh%vi2, 1:nz), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_V%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        do k = 1, nz
          d = modulo( ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))) + k,2) == 0
          d_mesh_loc_ex( vi,k) = d
          d_mesh_nih_ex( vi,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        do k = 1, nz
          d = modulo( ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)) + k),2) == 0
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) .eqv. d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) .eqv. d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping_logical_3D

  subroutine test_mesh_vertices_to_graph_mapping_int_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_mesh_vertices_to_graph_mapping_int_2D'
    character(len=1024), parameter             :: test_name_local = 'int_2D'
    character(len=1024)                        :: test_name
    integer, dimension(:), allocatable         :: d_mesh_loc_ex
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    integer, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    integer, dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                    :: vi, ni
    integer                                    :: d
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex( mesh%vi1 :mesh%vi2), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_V%n_nih)
    d_mesh_nih_ex( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        d = ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2)))
        d_mesh_loc_ex( vi) = d
        d_mesh_nih_ex( vi) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)))
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) == d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) == d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping_int_2D

  subroutine test_mesh_vertices_to_graph_mapping_int_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_mesh_vertices_to_graph_mapping_int_3D'
    character(len=1024), parameter               :: test_name_local = 'int_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    integer, dimension(:,:), allocatable         :: d_mesh_loc_ex
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                      :: vi, ni, k
    integer                                      :: d
    logical                                      :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%vi1 :mesh%vi2, 1:nz), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_V%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))) + k
          d_mesh_loc_ex( vi,k) = d
          d_mesh_nih_ex( vi,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)) + k)
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) == d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) == d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping_int_3D

  subroutine test_mesh_vertices_to_graph_mapping_dp_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_mesh_vertices_to_graph_mapping_dp_2D'
    character(len=1024), parameter              :: test_name_local = 'dp_2D'
    character(len=1024)                         :: test_name
    real(dp), dimension(:), allocatable         :: d_mesh_loc_ex
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                               :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                     :: vi, ni
    real(dp)                                    :: d
    logical                                     :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex( mesh%vi1 :mesh%vi2), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_V%n_nih)
    d_mesh_nih_ex( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        d = 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))
        d_mesh_loc_ex( vi) = d
        d_mesh_nih_ex( vi) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) == d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) == d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping_dp_2D

  subroutine test_mesh_vertices_to_graph_mapping_dp_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_mesh_vertices_to_graph_mapping_dp_3D'
    character(len=1024), parameter                :: test_name_local = 'dp_3D'
    character(len=1024)                           :: test_name
    integer                                       :: nz
    real(dp), dimension(:,:), allocatable         :: d_mesh_loc_ex
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                 :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                       :: vi, ni, k
    real(dp)                                      :: d
    logical                                       :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%vi1 :mesh%vi2, 1:nz), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_V%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        do k = 1, nz
          d = 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2)) + real( k,dp)
          d_mesh_loc_ex( vi,k) = d
          d_mesh_nih_ex( vi,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        do k = 1, nz
          d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)) + real( k,dp)
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) == d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) == d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping_dp_3D

  ! Mesh triangles to graph
  ! =======================

  subroutine test_mesh_triangles_to_graph_mapping( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_mesh_triangles_to_graph_mapping'
    character(len=1024), parameter :: test_name_local = 'mesh_triangles_to_graph'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_mesh_triangles_to_graph_mapping_logical_2D( test_name, mesh, graph)
    call test_mesh_triangles_to_graph_mapping_logical_3D( test_name, mesh, graph)
    call test_mesh_triangles_to_graph_mapping_int_2D    ( test_name, mesh, graph)
    call test_mesh_triangles_to_graph_mapping_int_3D    ( test_name, mesh, graph)
    call test_mesh_triangles_to_graph_mapping_dp_2D     ( test_name, mesh, graph)
    call test_mesh_triangles_to_graph_mapping_dp_3D     ( test_name, mesh, graph)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping

  subroutine test_mesh_triangles_to_graph_mapping_logical_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_mesh_triangles_to_graph_mapping_logical_2D'
    character(len=1024), parameter             :: test_name_local = 'logical_2D'
    character(len=1024)                        :: test_name
    logical, dimension(:), allocatable         :: d_mesh_loc_ex
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    logical, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    logical, dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                    :: ti, ni
    logical                                    :: d
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex( mesh%ti1 :mesh%ti2), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_Tri%n_nih)
    d_mesh_nih_ex( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        d = modulo( ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))),2) == 0
        d_mesh_loc_ex( ti) = d
        d_mesh_nih_ex( ti) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        d = modulo( ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))),2) == 0
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        test_result = test_result .and. (d_graph_nih_ex( ni) .eqv. d_graph_nih_mapped( ni))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        test_result = test_result .and. (d_graph_nih_ex( ni) .eqv. d_graph_nih_mapped( ni))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping_logical_2D

  subroutine test_mesh_triangles_to_graph_mapping_logical_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_mesh_triangles_to_graph_mapping_logical_3D'
    character(len=1024), parameter               :: test_name_local = 'logical_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    logical, dimension(:,:), allocatable         :: d_mesh_loc_ex
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                      :: ti, ni, k
    logical                                      :: d
    logical                                      :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%ti1 :mesh%ti2, 1:nz), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_Tri%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        do k = 1, nz
          d = modulo( ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))) + k,2) == 0
          d_mesh_loc_ex( ti,k) = d
          d_mesh_nih_ex( ti,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        do k = 1, nz
          d = modulo( ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)) + k),2) == 0
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_graph_nih_ex( ni,k) .eqv. d_graph_nih_mapped( ni,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_graph_nih_ex( ni,k) .eqv. d_graph_nih_mapped( ni,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping_logical_3D

  subroutine test_mesh_triangles_to_graph_mapping_int_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_mesh_triangles_to_graph_mapping_int_2D'
    character(len=1024), parameter             :: test_name_local = 'int_2D'
    character(len=1024)                        :: test_name
    integer, dimension(:), allocatable         :: d_mesh_loc_ex
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    integer, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    integer, dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                    :: ti, ni
    integer                                    :: d
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex( mesh%ti1 :mesh%ti2), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_Tri%n_nih)
    d_mesh_nih_ex( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        d = ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2)))
        d_mesh_loc_ex( ti) = d
        d_mesh_nih_ex( ti) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)))
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) == d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) == d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping_int_2D

  subroutine test_mesh_triangles_to_graph_mapping_int_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_mesh_triangles_to_graph_mapping_int_3D'
    character(len=1024), parameter               :: test_name_local = 'int_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    integer, dimension(:,:), allocatable         :: d_mesh_loc_ex
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                      :: ti, ni, k
    integer                                      :: d
    logical                                      :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%ti1 :mesh%ti2, 1:nz), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_Tri%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))) + k
          d_mesh_loc_ex( ti,k) = d
          d_mesh_nih_ex( ti,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)) + k)
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) == d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) == d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping_int_3D

  subroutine test_mesh_triangles_to_graph_mapping_dp_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_mesh_triangles_to_graph_mapping_dp_2D'
    character(len=1024), parameter              :: test_name_local = 'dp_2D'
    character(len=1024)                         :: test_name
    real(dp), dimension(:), allocatable         :: d_mesh_loc_ex
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                               :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                     :: ti, ni
    real(dp)                                    :: d
    logical                                     :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex( mesh%ti1 :mesh%ti2), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_Tri%n_nih)
    d_mesh_nih_ex( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        d = 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))
        d_mesh_loc_ex( ti) = d
        d_mesh_nih_ex( ti) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) == d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) == d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping_dp_2D

  subroutine test_mesh_triangles_to_graph_mapping_dp_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_mesh_triangles_to_graph_mapping_dp_3D'
    character(len=1024), parameter                :: test_name_local = 'dp_3D'
    character(len=1024)                           :: test_name
    integer                                       :: nz
    real(dp), dimension(:,:), allocatable         :: d_mesh_loc_ex
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                 :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                       :: ti, ni, k
    real(dp)                                      :: d
    logical                                       :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%ti1 :mesh%ti2, 1:nz), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_Tri%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        do k = 1, nz
          d = 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2)) + real( k,dp)
          d_mesh_loc_ex( ti,k) = d
          d_mesh_nih_ex( ti,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        do k = 1, nz
          d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)) + real( k,dp)
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) == d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) == d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping_dp_3D

  ! Mesh edges to graph
  ! ===================

  subroutine test_mesh_edges_to_graph_mapping( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_mesh_edges_to_graph_mapping'
    character(len=1024), parameter :: test_name_local = 'mesh_edges_to_graph'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_mesh_edges_to_graph_mapping_logical_2D( test_name, mesh, graph)
    call test_mesh_edges_to_graph_mapping_logical_3D( test_name, mesh, graph)
    call test_mesh_edges_to_graph_mapping_int_2D    ( test_name, mesh, graph)
    call test_mesh_edges_to_graph_mapping_int_3D    ( test_name, mesh, graph)
    call test_mesh_edges_to_graph_mapping_dp_2D     ( test_name, mesh, graph)
    call test_mesh_edges_to_graph_mapping_dp_3D     ( test_name, mesh, graph)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_edges_to_graph_mapping

  subroutine test_mesh_edges_to_graph_mapping_logical_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_mesh_edges_to_graph_mapping_logical_2D'
    character(len=1024), parameter             :: test_name_local = 'logical_2D'
    character(len=1024)                        :: test_name
    logical, dimension(:), allocatable         :: d_mesh_loc_ex
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    logical, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    logical, dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                    :: ei, ni
    logical                                    :: d
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex( mesh%ei1 :mesh%ei2), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_E%n_nih)
    d_mesh_nih_ex( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        d = modulo( ceiling( 1000._dp * hypot( mesh%E( ei,1), mesh%E( ei,2))),2) == 0
        d_mesh_loc_ex( ei) = d
        d_mesh_nih_ex( ei) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        d = modulo( ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))),2) == 0
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_edges_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        test_result = test_result .and. (d_graph_nih_ex( ni) .eqv. d_graph_nih_mapped( ni))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_edges_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        test_result = test_result .and. (d_graph_nih_ex( ni) .eqv. d_graph_nih_mapped( ni))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_edges_to_graph_mapping_logical_2D

  subroutine test_mesh_edges_to_graph_mapping_logical_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_mesh_edges_to_graph_mapping_logical_3D'
    character(len=1024), parameter               :: test_name_local = 'logical_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    logical, dimension(:,:), allocatable         :: d_mesh_loc_ex
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                      :: ei, ni, k
    logical                                      :: d
    logical                                      :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%ei1 :mesh%ei2, 1:nz), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_E%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih, 1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        do k = 1, nz
          d = modulo( ceiling( 1000._dp * hypot( mesh%E( ei,1), mesh%E( ei,2))) + k,2) == 0
          d_mesh_loc_ex( ei,k) = d
          d_mesh_nih_ex( ei,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        do k = 1, nz
          d = modulo( ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)) + k),2) == 0
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_edges_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_graph_nih_ex( ni,k) .eqv. d_graph_nih_mapped( ni,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_edges_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_graph_nih_ex( ni,k) .eqv. d_graph_nih_mapped( ni,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_edges_to_graph_mapping_logical_3D

  subroutine test_mesh_edges_to_graph_mapping_int_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_mesh_edges_to_graph_mapping_int_2D'
    character(len=1024), parameter             :: test_name_local = 'int_2D'
    character(len=1024)                        :: test_name
    integer, dimension(:), allocatable         :: d_mesh_loc_ex
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    integer, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    integer, dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                    :: ei, ni
    integer                                    :: d
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex( mesh%ei1 :mesh%ei2), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_E%n_nih)
    d_mesh_nih_ex( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        d = ceiling( 1000._dp * hypot( mesh%E( ei,1), mesh%E( ei,2)))
        d_mesh_loc_ex( ei) = d
        d_mesh_nih_ex( ei) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)))
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_edges_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) == d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_edges_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) == d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_edges_to_graph_mapping_int_2D

  subroutine test_mesh_edges_to_graph_mapping_int_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_mesh_edges_to_graph_mapping_int_3D'
    character(len=1024), parameter               :: test_name_local = 'int_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    integer, dimension(:,:), allocatable         :: d_mesh_loc_ex
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                      :: ei, ni, k
    integer                                      :: d
    logical                                      :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%ei1 :mesh%ei2, 1:nz), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_E%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih, 1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( mesh%E( ei,1), mesh%E( ei,2))) + k
          d_mesh_loc_ex( ei,k) = d
          d_mesh_nih_ex( ei,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)) + k)
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_edges_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) == d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_edges_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) == d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_edges_to_graph_mapping_int_3D

  subroutine test_mesh_edges_to_graph_mapping_dp_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_mesh_edges_to_graph_mapping_dp_2D'
    character(len=1024), parameter              :: test_name_local = 'dp_2D'
    character(len=1024)                         :: test_name
    real(dp), dimension(:), allocatable         :: d_mesh_loc_ex
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                               :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                     :: ei, ni
    real(dp)                                    :: d
    logical                                     :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex( mesh%ei1 :mesh%ei2), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_E%n_nih)
    d_mesh_nih_ex( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        d = 1000._dp * hypot( mesh%E( ei,1), mesh%E( ei,2))
        d_mesh_loc_ex( ei) = d
        d_mesh_nih_ex( ei) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_edges_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) == d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_edges_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        test_result = test_result .and. d_graph_nih_ex( ni) == d_graph_nih_mapped( ni)
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_edges_to_graph_mapping_dp_2D

  subroutine test_mesh_edges_to_graph_mapping_dp_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_mesh_edges_to_graph_mapping_dp_3D'
    character(len=1024), parameter                :: test_name_local = 'dp_3D'
    character(len=1024)                           :: test_name
    integer                                       :: nz
    real(dp), dimension(:,:), allocatable         :: d_mesh_loc_ex
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                 :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                       :: ei, ni, k
    real(dp)                                      :: d
    logical                                       :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%ei1 :mesh%ei2, 1:nz), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_E%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih, 1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        do k = 1, nz
          d = 1000._dp * hypot( mesh%E( ei,1), mesh%E( ei,2)) + real( k,dp)
          d_mesh_loc_ex( ei,k) = d
          d_mesh_nih_ex( ei,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        do k = 1, nz
          d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)) + real( k,dp)
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map distributed data from the mesh to the graph
    call map_mesh_edges_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) == d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_mesh_edges_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_ex( ni,k) == d_graph_nih_mapped( ni,k)
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_edges_to_graph_mapping_dp_3D

  ! Graph to mesh vertices
  ! ======================

  subroutine test_graph_to_mesh_vertices_mapping( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_to_mesh_vertices_mapping'
    character(len=1024), parameter :: test_name_local = 'graph_to_mesh_vertices'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_graph_to_mesh_vertices_mapping_logical_2D( test_name, mesh, graph)
    call test_graph_to_mesh_vertices_mapping_logical_3D( test_name, mesh, graph)
    call test_graph_to_mesh_vertices_mapping_int_2D    ( test_name, mesh, graph)
    call test_graph_to_mesh_vertices_mapping_int_3D    ( test_name, mesh, graph)
    call test_graph_to_mesh_vertices_mapping_dp_2D     ( test_name, mesh, graph)
    call test_graph_to_mesh_vertices_mapping_dp_3D     ( test_name, mesh, graph)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping

  subroutine test_graph_to_mesh_vertices_mapping_logical_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_graph_to_mesh_vertices_mapping_logical_2D'
    character(len=1024), parameter             :: test_name_local = 'logical_2D'
    character(len=1024)                        :: test_name
    logical, dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    logical, dimension(:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                    :: vi, ni
    logical                                    :: d
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%vi1 :mesh%vi2), source = .false.)
    allocate( d_mesh_loc_mapped( mesh%vi1 :mesh%vi2), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_V%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_V%n_nih)
    d_mesh_nih_ex    ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        d = modulo( ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))),2) == 0
        d_mesh_loc_ex( vi) = d
        d_mesh_nih_ex( vi) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        d = modulo( ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))),2) == 0
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( vi) .eqv. d_mesh_loc_mapped( vi))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( vi) .eqv. d_mesh_nih_mapped( vi))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping_logical_2D

  subroutine test_graph_to_mesh_vertices_mapping_logical_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_graph_to_mesh_vertices_mapping_logical_3D'
    character(len=1024), parameter               :: test_name_local = 'logical_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    logical, dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                      :: vi, ni, k
    logical                                      :: d
    logical                                      :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%vi1 :mesh%vi2, 1:nz), source = .false.)
    allocate( d_mesh_loc_mapped( mesh%vi1 :mesh%vi2, 1:nz), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_V%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_V%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        do k = 1, nz
          d = modulo( ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))) + k,2) == 0
          d_mesh_loc_ex( vi,k) = d
          d_mesh_nih_ex( vi,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        do k = 1, nz
          d = modulo( ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))) + k,2) == 0
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( vi,k) .eqv. d_mesh_loc_mapped( vi,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( vi,k) .eqv. d_mesh_nih_mapped( vi,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping_logical_3D

  subroutine test_graph_to_mesh_vertices_mapping_int_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_graph_to_mesh_vertices_mapping_int_2D'
    character(len=1024), parameter             :: test_name_local = 'int_2D'
    character(len=1024)                        :: test_name
    integer, dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    integer, dimension(:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                    :: vi, ni
    integer                                    :: d
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%vi1 :mesh%vi2), source = 0)
    allocate( d_mesh_loc_mapped( mesh%vi1 :mesh%vi2), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_V%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_V%n_nih)
    d_mesh_nih_ex    ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        d = ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2)))
        d_mesh_loc_ex( vi) = d
        d_mesh_nih_ex( vi) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)))
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( vi) == d_mesh_loc_mapped( vi))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( vi) == d_mesh_nih_mapped( vi))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping_int_2D

  subroutine test_graph_to_mesh_vertices_mapping_int_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_graph_to_mesh_vertices_mapping_int_3D'
    character(len=1024), parameter               :: test_name_local = 'int_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    integer, dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                      :: vi, ni, k
    integer                                      :: d
    logical                                      :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%vi1 :mesh%vi2, 1:nz), source = 0)
    allocate( d_mesh_loc_mapped( mesh%vi1 :mesh%vi2, 1:nz), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_V%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_V%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))) + k
          d_mesh_loc_ex( vi,k) = d
          d_mesh_nih_ex( vi,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))) + k
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( vi,k) == d_mesh_loc_mapped( vi,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( vi,k) == d_mesh_nih_mapped( vi,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping_int_3D

  subroutine test_graph_to_mesh_vertices_mapping_dp_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_graph_to_mesh_vertices_mapping_dp_2D'
    character(len=1024), parameter              :: test_name_local = 'dp_2D'
    character(len=1024)                         :: test_name
    real(dp), dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                               :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                     :: vi, ni
    real(dp)                                    :: d
    logical                                     :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%vi1 :mesh%vi2), source = 0._dp)
    allocate( d_mesh_loc_mapped( mesh%vi1 :mesh%vi2), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_V%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_V%n_nih)
    d_mesh_nih_ex    ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        d = 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))
        d_mesh_loc_ex( vi) = d
        d_mesh_nih_ex( vi) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( vi) == d_mesh_loc_mapped( vi))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( vi) == d_mesh_nih_mapped( vi))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping_dp_2D

  subroutine test_graph_to_mesh_vertices_mapping_dp_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_graph_to_mesh_vertices_mapping_dp_3D'
    character(len=1024), parameter                :: test_name_local = 'dp_3D'
    character(len=1024)                           :: test_name
    integer                                       :: nz
    real(dp), dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                                 :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                       :: vi, ni, k
    real(dp)                                      :: d
    logical                                       :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%vi1 :mesh%vi2, 1:nz), source = 0._dp)
    allocate( d_mesh_loc_mapped( mesh%vi1 :mesh%vi2, 1:nz), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_V%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_V%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        do k = 1, nz
          d = 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))
          d_mesh_loc_ex( vi,k) = d
          d_mesh_nih_ex( vi,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        do k = 1, nz
          d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( vi,k) == d_mesh_loc_mapped( vi,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( vi,k) == d_mesh_nih_mapped( vi,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping_dp_3D

  ! Graph to mesh triangles
  ! ======================-

  subroutine test_graph_to_mesh_triangles_mapping( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_to_mesh_triangles_mapping'
    character(len=1024), parameter :: test_name_local = 'graph_to_mesh_triangles'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_graph_to_mesh_triangles_mapping_logical_2D( test_name, mesh, graph)
    call test_graph_to_mesh_triangles_mapping_logical_3D( test_name, mesh, graph)
    call test_graph_to_mesh_triangles_mapping_int_2D    ( test_name, mesh, graph)
    call test_graph_to_mesh_triangles_mapping_int_3D    ( test_name, mesh, graph)
    call test_graph_to_mesh_triangles_mapping_dp_2D     ( test_name, mesh, graph)
    call test_graph_to_mesh_triangles_mapping_dp_3D     ( test_name, mesh, graph)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping

  subroutine test_graph_to_mesh_triangles_mapping_logical_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_graph_to_mesh_triangles_mapping_logical_2D'
    character(len=1024), parameter             :: test_name_local = 'logical_2D'
    character(len=1024)                        :: test_name
    logical, dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    logical, dimension(:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                    :: ti, ni
    logical                                    :: d
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2), source = .false.)
    allocate( d_mesh_loc_mapped( mesh%ti1 :mesh%ti2), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_Tri%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_Tri%n_nih)
    d_mesh_nih_ex    ( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        d = modulo( ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))),2) == 0
        d_mesh_loc_ex( ti) = d
        d_mesh_nih_ex( ti) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        d = modulo( ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))),2) == 0
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( ti) .eqv. d_mesh_loc_mapped( ti))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( ti) .eqv. d_mesh_nih_mapped( ti))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping_logical_2D

  subroutine test_graph_to_mesh_triangles_mapping_logical_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_graph_to_mesh_triangles_mapping_logical_3D'
    character(len=1024), parameter               :: test_name_local = 'logical_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    logical, dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                      :: ti, ni, k
    logical                                      :: d
    logical                                      :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2, 1:nz), source = .false.)
    allocate( d_mesh_loc_mapped( mesh%ti1 :mesh%ti2, 1:nz), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_Tri%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_Tri%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        do k = 1, nz
          d = modulo( ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))) + k,2) == 0
          d_mesh_loc_ex( ti,k) = d
          d_mesh_nih_ex( ti,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        do k = 1, nz
          d = modulo( ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))) + k,2) == 0
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( ti,k) .eqv. d_mesh_loc_mapped( ti,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( ti,k) .eqv. d_mesh_nih_mapped( ti,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping_logical_3D

  subroutine test_graph_to_mesh_triangles_mapping_int_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_graph_to_mesh_triangles_mapping_int_2D'
    character(len=1024), parameter             :: test_name_local = 'int_2D'
    character(len=1024)                        :: test_name
    integer, dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    integer, dimension(:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                    :: ti, ni
    integer                                    :: d
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2), source = 0)
    allocate( d_mesh_loc_mapped( mesh%ti1 :mesh%ti2), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_Tri%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_Tri%n_nih)
    d_mesh_nih_ex    ( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        d = ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2)))
        d_mesh_loc_ex( ti) = d
        d_mesh_nih_ex( ti) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)))
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( ti) == d_mesh_loc_mapped( ti))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( ti) == d_mesh_nih_mapped( ti))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping_int_2D

  subroutine test_graph_to_mesh_triangles_mapping_int_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_graph_to_mesh_triangles_mapping_int_3D'
    character(len=1024), parameter               :: test_name_local = 'int_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    integer, dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                      :: ti, ni, k
    integer                                      :: d
    logical                                      :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2, 1:nz), source = 0)
    allocate( d_mesh_loc_mapped( mesh%ti1 :mesh%ti2, 1:nz), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_Tri%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_Tri%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))) + k
          d_mesh_loc_ex( ti,k) = d
          d_mesh_nih_ex( ti,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))) + k
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( ti,k) == d_mesh_loc_mapped( ti,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( ti,k) == d_mesh_nih_mapped( ti,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping_int_3D

  subroutine test_graph_to_mesh_triangles_mapping_dp_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_graph_to_mesh_triangles_mapping_dp_2D'
    character(len=1024), parameter              :: test_name_local = 'dp_2D'
    character(len=1024)                         :: test_name
    real(dp), dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                               :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                     :: ti, ni
    real(dp)                                    :: d
    logical                                     :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2), source = 0._dp)
    allocate( d_mesh_loc_mapped( mesh%ti1 :mesh%ti2), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_Tri%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_Tri%n_nih)
    d_mesh_nih_ex    ( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        d = 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))
        d_mesh_loc_ex( ti) = d
        d_mesh_nih_ex( ti) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( ti) == d_mesh_loc_mapped( ti))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( ti) == d_mesh_nih_mapped( ti))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping_dp_2D

  subroutine test_graph_to_mesh_triangles_mapping_dp_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_graph_to_mesh_triangles_mapping_dp_3D'
    character(len=1024), parameter                :: test_name_local = 'dp_3D'
    character(len=1024)                           :: test_name
    integer                                       :: nz
    real(dp), dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                                 :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                       :: ti, ni, k
    real(dp)                                      :: d
    logical                                       :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2, 1:nz), source = 0._dp)
    allocate( d_mesh_loc_mapped( mesh%ti1 :mesh%ti2, 1:nz), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_Tri%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_Tri%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        do k = 1, nz
          d = 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))
          d_mesh_loc_ex( ti,k) = d
          d_mesh_nih_ex( ti,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        do k = 1, nz
          d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( ti,k) == d_mesh_loc_mapped( ti,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( ti,k) == d_mesh_nih_mapped( ti,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping_dp_3D

  ! Graph to mesh edges
  ! ===================

  subroutine test_graph_to_mesh_edges_mapping( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_to_mesh_edges_mapping'
    character(len=1024), parameter :: test_name_local = 'graph_to_mesh_edges'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_graph_to_mesh_edges_mapping_logical_2D( test_name, mesh, graph)
    call test_graph_to_mesh_edges_mapping_logical_3D( test_name, mesh, graph)
    call test_graph_to_mesh_edges_mapping_int_2D    ( test_name, mesh, graph)
    call test_graph_to_mesh_edges_mapping_int_3D    ( test_name, mesh, graph)
    call test_graph_to_mesh_edges_mapping_dp_2D     ( test_name, mesh, graph)
    call test_graph_to_mesh_edges_mapping_dp_3D     ( test_name, mesh, graph)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_edges_mapping

  subroutine test_graph_to_mesh_edges_mapping_logical_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_graph_to_mesh_edges_mapping_logical_2D'
    character(len=1024), parameter             :: test_name_local = 'logical_2D'
    character(len=1024)                        :: test_name
    logical, dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    logical, dimension(:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                    :: ei, ni
    logical                                    :: d
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%ei1 :mesh%ei2), source = .false.)
    allocate( d_mesh_loc_mapped( mesh%ei1 :mesh%ei2), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_E%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_E%n_nih)
    d_mesh_nih_ex    ( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        d = modulo( ceiling( 1000._dp * hypot( mesh%E( ei,1), mesh%E( ei,2))),2) == 0
        d_mesh_loc_ex( ei) = d
        d_mesh_nih_ex( ei) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        d = modulo( ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))),2) == 0
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_edges( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( ei) .eqv. d_mesh_loc_mapped( ei))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_edges( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( ei) .eqv. d_mesh_nih_mapped( ei))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_edges_mapping_logical_2D

  subroutine test_graph_to_mesh_edges_mapping_logical_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_graph_to_mesh_edges_mapping_logical_3D'
    character(len=1024), parameter               :: test_name_local = 'logical_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    logical, dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                      :: ei, ni, k
    logical                                      :: d
    logical                                      :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%ei1 :mesh%ei2, 1:nz), source = .false.)
    allocate( d_mesh_loc_mapped( mesh%ei1 :mesh%ei2, 1:nz), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_E%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_E%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        do k = 1, nz
          d = modulo( ceiling( 1000._dp * hypot( mesh%E( ei,1), mesh%E( ei,2))) + k,2) == 0
          d_mesh_loc_ex( ei,k) = d
          d_mesh_nih_ex( ei,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        do k = 1, nz
          d = modulo( ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))) + k,2) == 0
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_edges( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( ei,k) .eqv. d_mesh_loc_mapped( ei,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_edges( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( ei,k) .eqv. d_mesh_nih_mapped( ei,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_edges_mapping_logical_3D

  subroutine test_graph_to_mesh_edges_mapping_int_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_graph_to_mesh_edges_mapping_int_2D'
    character(len=1024), parameter             :: test_name_local = 'int_2D'
    character(len=1024)                        :: test_name
    integer, dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    integer, dimension(:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                    :: ei, ni
    integer                                    :: d
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%ei1 :mesh%ei2), source = 0)
    allocate( d_mesh_loc_mapped( mesh%ei1 :mesh%ei2), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_E%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_E%n_nih)
    d_mesh_nih_ex    ( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        d = ceiling( 1000._dp * hypot( mesh%E( ei,1), mesh%E( ei,2)))
        d_mesh_loc_ex( ei) = d
        d_mesh_nih_ex( ei) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)))
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_edges( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( ei) == d_mesh_loc_mapped( ei))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_edges( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( ei) == d_mesh_nih_mapped( ei))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_edges_mapping_int_2D

  subroutine test_graph_to_mesh_edges_mapping_int_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_graph_to_mesh_edges_mapping_int_3D'
    character(len=1024), parameter               :: test_name_local = 'int_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    integer, dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                      :: ei, ni, k
    integer                                      :: d
    logical                                      :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%ei1 :mesh%ei2, 1:nz), source = 0)
    allocate( d_mesh_loc_mapped( mesh%ei1 :mesh%ei2, 1:nz), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_E%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_E%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( mesh%E( ei,1), mesh%E( ei,2))) + k
          d_mesh_loc_ex( ei,k) = d
          d_mesh_nih_ex( ei,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))) + k
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_edges( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( ei,k) == d_mesh_loc_mapped( ei,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_edges( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( ei,k) == d_mesh_nih_mapped( ei,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_edges_mapping_int_3D

  subroutine test_graph_to_mesh_edges_mapping_dp_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_graph_to_mesh_edges_mapping_dp_2D'
    character(len=1024), parameter              :: test_name_local = 'dp_2D'
    character(len=1024)                         :: test_name
    real(dp), dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                               :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                     :: ei, ni
    real(dp)                                    :: d
    logical                                     :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%ei1 :mesh%ei2), source = 0._dp)
    allocate( d_mesh_loc_mapped( mesh%ei1 :mesh%ei2), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_E%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_E%n_nih)
    d_mesh_nih_ex    ( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        d = 1000._dp * hypot( mesh%E( ei,1), mesh%E( ei,2))
        d_mesh_loc_ex( ei) = d
        d_mesh_nih_ex( ei) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))
        d_graph_nih_ex( ni) = d
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_edges( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( ei) == d_mesh_loc_mapped( ei))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_edges( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        test_result = test_result .and. (d_mesh_loc_ex( ei) == d_mesh_nih_mapped( ei))
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_edges_mapping_dp_2D

  subroutine test_graph_to_mesh_edges_mapping_dp_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_graph_to_mesh_edges_mapping_dp_3D'
    character(len=1024), parameter                :: test_name_local = 'dp_3D'
    character(len=1024)                           :: test_name
    integer                                       :: nz
    real(dp), dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_ex    => null()
    type(MPI_WIN)                                 :: wd_mesh_nih_ex, wd_mesh_nih_mapped, wd_graph_nih_ex
    integer                                       :: ei, ni, k
    real(dp)                                      :: d
    logical                                       :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%ei1 :mesh%ei2, 1:nz), source = 0._dp)
    allocate( d_mesh_loc_mapped( mesh%ei1 :mesh%ei2, 1:nz), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_E%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_E%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        do k = 1, nz
          d = 1000._dp * hypot( mesh%E( ei,1), mesh%E( ei,2))
          d_mesh_loc_ex( ei,k) = d
          d_mesh_nih_ex( ei,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        do k = 1, nz
          d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))
          d_graph_nih_ex( ni,k) = d
        end do
      end if
    end do

    ! Map data from the graph to the mesh
    call map_graph_to_mesh_edges( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)

    test_result = .true.
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( ei,k) == d_mesh_loc_mapped( ei,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_dist')

    ! Map hybrid distributed/shared data from the mesh to the graph
    call map_graph_to_mesh_edges( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)

    test_result = .true.
    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        do k = 1, nz
          test_result = test_result .and. (d_mesh_loc_ex( ei,k) == d_mesh_nih_mapped( ei,k))
        end do
      end if
    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD)
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_edges_mapping_dp_3D

end module ut_mesh_graphs_mapping