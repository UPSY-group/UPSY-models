module tests_graph

  ! Basic tests for graphs.

  use precisions, only: dp
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mpi_basic, only: par

  implicit none

  private

  public :: test_graph_is_self_consistent

contains

  function test_graph_is_self_consistent( mesh, graph) result( isso)

    ! In/output variables:
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_is_self_consistent'

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.

    isso = isso .and. test_graph_matches_mesh( mesh, graph)
    isso = isso .and. test_graph_connectivity_is_self_consistent( graph)
    isso = isso .and. test_graph_to_mesh_connectivity_is_self_consistent( graph)
    isso = isso .and. test_graph_nodes_are_sorted( graph)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_is_self_consistent

  function test_graph_matches_mesh( mesh, graph) result( isso)

    ! In/output variables:
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_matches_mesh'

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.
    isso = isso .and. test_graph_matches_mesh_vertices ( mesh, graph)
    isso = isso .and. test_graph_matches_mesh_triangles( mesh, graph)
    ! isso = isso .and. test_graph_matches_mesh_edges    ( mesh, graph)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_matches_mesh

  function test_graph_matches_mesh_vertices( mesh, graph) result( isso)

    ! In/output variables:
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_matches_mesh'
    integer                        :: ni, vi

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.
    do ni = 1, graph%n
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        isso = isso .and. norm2( graph%V( ni,:) - mesh%V( vi,:)) < mesh%tol_dist
      end if
    end do

    if (par%primary .and. .not. isso) call warning('graph node coordinates do not match corresponding mesh vertices')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_matches_mesh_vertices

  function test_graph_matches_mesh_triangles( mesh, graph) result( isso)

    ! In/output variables:
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_matches_mesh_triangles'
    integer                        :: ni, ti

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.
    do ni = 1, graph%n
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        isso = isso .and. norm2( graph%V( ni,:) - mesh%TriGC( ti,:)) < mesh%tol_dist
      end if
    end do

    if (par%primary .and. .not. isso) call warning('graph node coordinates do not match corresponding mesh triangles')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_matches_mesh_triangles

  function test_graph_matches_mesh_edges( mesh, graph) result( isso)

    ! In/output variables:
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_matches_mesh_edges'
    integer                        :: ni, ei

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.
    do ni = 1, graph%n
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        isso = isso .and. norm2( graph%V( ni,:) - mesh%E( ei,:)) < mesh%tol_dist
      end if
    end do

    if (par%primary .and. .not. isso) call warning('graph node coordinates do not match corresponding mesh edges')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_matches_mesh_edges

  function test_graph_connectivity_is_self_consistent( graph) result( isso)

    ! In/output variables:
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_connectivity_is_self_consistent'
    integer                        :: ni, ci, nj, cj, nk
    logical                        :: found_reverse

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.

    do ni = 1, graph%n
      do ci = 1, graph%nC( ni)
        nj = graph%C( ni,ci)
        found_reverse = .false.
        do cj = 1, graph%nC( nj)
          nk = graph%C( nj,cj)
          if (nk == ni) then
            found_reverse = .true.
            exit
          end if
        end do
        isso = isso .and. found_reverse
      end do
    end do

    if (par%primary .and. .not. isso) call warning('graph connectivity is not self-consistent')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_connectivity_is_self_consistent

  function test_graph_to_mesh_connectivity_is_self_consistent( graph) result( isso)

    ! In/output variables:
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_to_mesh_connectivity_is_self_consistent'

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.
    isso = isso .and. test_graph_to_mesh_vertex_connectivity_is_self_consistent  ( graph)
    isso = isso .and. test_graph_to_mesh_triangle_connectivity_is_self_consistent( graph)
    isso = isso .and. test_graph_to_mesh_edge_connectivity_is_self_consistent    ( graph)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_to_mesh_connectivity_is_self_consistent

  function test_graph_to_mesh_vertex_connectivity_is_self_consistent( graph) result( isso)

    ! In/output variables:
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_to_mesh_vertex_connectivity_is_self_consistent'
    integer                        :: ni, vi

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.

    do ni = 1, graph%n
      vi = graph%ni2vi( ni)
      if (vi > 0) isso = isso .and. graph%vi2ni( vi) == ni
    end do

    do vi = 1, size( graph%vi2ni)
      ni = graph%vi2ni( vi)
      if (ni > 0) isso = isso .and. graph%ni2vi( ni) == vi
    end do

    if (par%primary .and. .not. isso) call warning('graph to mesh vertex connectivity is not self-consistent')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_to_mesh_vertex_connectivity_is_self_consistent

  function test_graph_to_mesh_triangle_connectivity_is_self_consistent( graph) result( isso)

    ! In/output variables:
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_to_mesh_triangle_connectivity_is_self_consistent'
    integer                        :: ni, ti

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.

    do ni = 1, graph%n
      ti = graph%ni2ti( ni)
      if (ti > 0) isso = isso .and. graph%ti2ni( ti) == ni
    end do

    do ti = 1, size( graph%ti2ni)
      ni = graph%ti2ni( ti)
      if (ni > 0) isso = isso .and. graph%ni2ti( ni) == ti
    end do

    if (par%primary .and. .not. isso) call warning('graph to mesh triangle connectivity is not self-consistent')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_to_mesh_triangle_connectivity_is_self_consistent

  function test_graph_to_mesh_edge_connectivity_is_self_consistent( graph) result( isso)

    ! In/output variables:
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_to_mesh_edge_connectivity_is_self_consistent'
    integer                        :: ni, ei

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.

    do ni = 1, graph%n
      ei = graph%ni2ei( ni)
      if (ei > 0) isso = isso .and. graph%ei2ni( ei) == ni
    end do

    do ei = 1, size( graph%ei2ni)
      ni = graph%ei2ni( ei)
      if (ni > 0) isso = isso .and. graph%ni2ei( ni) == ei
    end do

    if (par%primary .and. .not. isso) call warning('graph to mesh edge connectivity is not self-consistent')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_to_mesh_edge_connectivity_is_self_consistent

  function test_graph_nodes_are_sorted( graph) result( isso)

    ! In/output variables:
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=256), parameter :: routine_name = 'test_graph_nodes_are_sorted'
    integer                       :: ni

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.
    do ni = 2, graph%n
      isso = isso .and. graph%V( ni,1) >= graph%V( ni-1,1)
    end do

    if (par%primary .and. .not. isso) call warning('graph nodes are not sorted in the x-dimension')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_nodes_are_sorted

end module tests_graph
