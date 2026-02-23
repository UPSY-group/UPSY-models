module delete_vertices_local_geometry

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use mesh_types, only: type_mesh
  use assertions_basic, only: assert
  use tests_main
  use mpi_basic, only: par
  use plane_geometry, only: cross2

  implicit none

  private

  public :: type_local_geometry_nCge4, find_local_geometry_nCge4

  type type_local_geometry_nCge4
    integer                            :: vi_kill
    integer                            :: vj_clock, vj_focus, vj_anti
    integer                            :: nvj_opp
    integer, dimension(:), allocatable :: vj_opp
    integer                            :: ti1, ti2, ti3, ti4
    integer                            :: nti_opp
    integer, dimension(:), allocatable :: ti_opp
    integer                            :: tj1, tj2, tj3, tj4
    integer, dimension(:), allocatable :: tj_opp
  end type type_local_geometry_nCge4

contains

  subroutine find_local_geometry_nCge4( mesh, vi, locgeom)

    ! Find the local geometry:
    !
    !                tj_opp()
    !
    !  vj_opp(nvj_opp) ------ vj_opp(1)
    !           /  \ ti_opp() /  \
    !     tj1  /    \        /    \  tj4
    !         / ti1  \      /  ti4 \
    !        /        \    /        \
    !   vj_clock ---- locgeom%vi_kill ---- vj_anti
    !         \          |          /
    !           \  ti2   |   ti3  /
    !             \      |      /
    !        tj2    \    |    /    tj3
    !                 vj_focus

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    integer,                         intent(in   ) :: vi
    type(type_local_geometry_nCge4), intent(  out) :: locgeom

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_nCge4_local_geometry'

    ! Add routine to path
    call init_routine( routine_name)

    call find_local_geometry_nCge4_vertices         ( mesh, vi, locgeom)
    call find_local_geometry_nCge4_inside_triangles ( mesh, locgeom)
    call find_local_geometry_nCge4_outside_triangles( mesh, locgeom)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine find_local_geometry_nCge4

  subroutine find_local_geometry_nCge4_vertices( mesh, vi, locgeom)
    !< Find a valid vj_focus (from where you can "see" all other neighbours of locgeom%vi_kill)
    !< ...which is not guaranteed, as they might span a non-convex polygon

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    integer,                         intent(in   ) :: vi
    type(type_local_geometry_nCge4), intent(inout) :: locgeom

    ! Local variables:
    logical                            :: found_it
    integer, dimension(:), allocatable :: CC
    integer                            :: nit

    locgeom%vi_kill = vi

    CC = mesh%C( locgeom%vi_kill, 1: mesh%nC( locgeom%vi_kill))

    found_it = .false.
    nit = 0
    do while (nit <= mesh%nC( locgeom%vi_kill))
      nit = nit + 1

      CC = [CC( 2: mesh%nC( locgeom%vi_kill)), CC( 1)]

      locgeom%vj_clock = CC( 1)
      locgeom%vj_focus = CC( 2)
      locgeom%vj_anti  = CC( 3)
      locgeom%nvj_opp  = mesh%nC( locgeom%vi_kill) - 3
      locgeom%vj_opp   = CC( 4: mesh%nC( locgeom%vi_kill))

      if (is_valid_order_vj_focus_nCge4( mesh, locgeom)) then
        found_it = .true.
        exit
      end if
    end do

#if (DO_ASSERTIONS)
    call assert( found_it, 'could not find vertices of local geometry')
    call assert( test_eq( size( locgeom%vj_opp), locgeom%nvj_opp), 'invalid nvj_opp')
    call assert( test_ge_le( locgeom%vj_clock, 1, mesh%nV), 'invalid vj_clock')
    call assert( test_ge_le( locgeom%vj_focus, 1, mesh%nV), 'invalid vj_focus')
    call assert( test_ge_le( locgeom%vj_anti , 1, mesh%nV), 'invalid vj_anti')
    call assert( test_ge_le( locgeom%vj_opp  , 1, mesh%nV), 'invalid vj_opp')
#endif

  end subroutine find_local_geometry_nCge4_vertices

  subroutine find_local_geometry_nCge4_inside_triangles( mesh, locgeom)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_local_geometry_nCge4), intent(inout) :: locgeom

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_nCge4_local_geometry_inside_triangles'
    logical                        :: found_it
    integer                        :: iti2, ti, n1, n2, iti1, iti3, iti4

    ! Add routine to path
    call init_routine( routine_name)

    found_it = .false.
    do iti2 = 1, mesh%niTri( locgeom%vi_kill)
      ti = mesh%iTri( locgeom%vi_kill,iti2)
      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1

        if (mesh%Tri( ti,n1) == locgeom%vj_focus .and. mesh%Tri( ti,n2) == locgeom%vi_kill) then

          iti1 = iti2 - 1
          if (iti1 == 0) iti1 = mesh%niTri( locgeom%vi_kill)
          locgeom%ti1 = mesh%iTri( locgeom%vi_kill, iti1)

          locgeom%ti2 = mesh%iTri( locgeom%vi_kill, iti2)

          iti3 = iti2 + 1
          if (iti3 > mesh%niTri( locgeom%vi_kill)) iti3 = 1
          locgeom%ti3 = mesh%iTri( locgeom%vi_kill, iti3)

          iti4 = iti3 + 1
          if (iti4 > mesh%niTri( locgeom%vi_kill)) iti4 = 1
          locgeom%ti4 = mesh%iTri( locgeom%vi_kill, iti4)

          found_it = .true.
          exit

        end if

      end do
      if (found_it) exit
    end do

#if (DO_ASSERTIONS)
    call assert( found_it, 'could not find inside triangles of local geometry')
    call assert( test_ge_le( locgeom%ti1, 1, mesh%nTri), 'invalid ti1')
    call assert( test_ge_le( locgeom%ti2, 1, mesh%nTri), 'invalid ti2')
    call assert( test_ge_le( locgeom%ti3, 1, mesh%nTri), 'invalid ti3')
    call assert( test_ge_le( locgeom%ti4, 1, mesh%nTri), 'invalid ti4')
#endif

    locgeom%nti_opp = mesh%niTri( locgeom%vi_kill) - 4
    if (iti4 > iti1) then
      ! [iti1,iti2,iti3,iti4] form a contiguous block within iTri(locgeom%vi_kill,:)
      locgeom%ti_opp = [mesh%iTri( locgeom%vi_kill, iti4+1:mesh%niTri( locgeom%vi_kill)), mesh%iTri( locgeom%vi_kill, 1:iti1-1)]
    else
      ! The block is non-contiguous
      locgeom%ti_opp = mesh%iTri( locgeom%vi_kill, iti4+1: iti1-1)
    end if

#if (DO_ASSERTIONS)
    call assert( test_eq( size( locgeom%ti_opp), locgeom%nti_opp), 'invalid nti_opp')
    call assert( test_ge_le( locgeom%ti_opp, 1, mesh%nTri), 'invalid ti_opp')
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine find_local_geometry_nCge4_inside_triangles

  subroutine find_local_geometry_nCge4_outside_triangles( mesh, locgeom)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_local_geometry_nCge4), intent(inout) :: locgeom

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_nCge4_local_geometry_outside_triangles'
    logical                        :: found_it
    integer                        :: n, tj, tii, ti

    ! Add routine to path
    call init_routine( routine_name)

    found_it = .false.
    do n = 1, 3
      tj = mesh%TriC( locgeom%ti1,n)
      if (mesh%Tri( locgeom%ti1,n) == locgeom%vi_kill) then
        found_it = .true.
        locgeom%tj1 = tj
        exit
      end if
    end do
#if (DO_ASSERTIONS)
    call assert( found_it, 'could not find tj1')
    call assert( test_ge_le( locgeom%tj1, 1, mesh%nTri), 'invalid tj1')
#endif

    found_it = .false.
    do n = 1, 3
      tj = mesh%TriC( locgeom%ti2,n)
      if (mesh%Tri( locgeom%ti2,n) == locgeom%vi_kill) then
        found_it = .true.
        locgeom%tj2 = tj
        exit
      end if
    end do
#if (DO_ASSERTIONS)
    call assert( found_it, 'could not find tj2')
    call assert( test_ge_le( locgeom%tj2, 1, mesh%nTri), 'invalid tj2')
#endif

    found_it = .false.
    do n = 1, 3
      tj = mesh%TriC( locgeom%ti3,n)
      if (mesh%Tri( locgeom%ti3,n) == locgeom%vi_kill) then
        found_it = .true.
        locgeom%tj3 = tj
        exit
      end if
    end do
#if (DO_ASSERTIONS)
    call assert( found_it, 'could not find tj3')
    call assert( test_ge_le( locgeom%tj3, 1, mesh%nTri), 'invalid tj3')
#endif

    found_it = .false.
    do n = 1, 3
      tj = mesh%TriC( locgeom%ti4,n)
      if (mesh%Tri( locgeom%ti4,n) == locgeom%vi_kill) then
        found_it = .true.
        locgeom%tj4 = tj
        exit
      end if
    end do
#if (DO_ASSERTIONS)
    call assert( found_it, 'could not find tj4')
    call assert( test_ge_le( locgeom%tj4, 1, mesh%nTri), 'invalid tj4')
#endif

    allocate( locgeom%tj_opp( locgeom%nti_opp), source = 0)
    do tii = 1, locgeom%nti_opp
      found_it = .false.
      ti = locgeom%ti_opp( tii)
      do n = 1, 3
        tj = mesh%TriC( ti,n)
        if (mesh%Tri( ti,n) == locgeom%vi_kill) then
          found_it = .true.
          locgeom%tj_opp( tii) = tj
          exit
        end if
      end do
#if (DO_ASSERTIONS)
      call assert( found_it, 'could not find tj_opp')
      call assert( test_ge_le( locgeom%tj_opp( tii), 1, mesh%nTri), 'invalid tj_opp')
#endif
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine find_local_geometry_nCge4_outside_triangles

  function is_valid_order_vj_focus_nCge4( mesh, locgeom) result( isso)
    ! All new triangles resulting from this ordering must have their
    ! vertices ordered counter-clockwise

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_local_geometry_nCge4), intent(in   ) :: locgeom
    logical                                        :: isso

    ! Local variables:
    integer, dimension(3) :: tri
    integer               :: vii

    isso = .true.

    ! ti3
    tri = [locgeom%vj_focus, locgeom%vj_anti, locgeom%vj_opp( 1)]
    isso = isso .and. is_ordered_counterclockwise( mesh, tri)

    ! ti_opp
    do vii = 1, locgeom%nvj_opp-1
      tri = [locgeom%vj_focus, locgeom%vj_opp( vii), locgeom%vj_opp( vii+1)]
      isso = isso .and. is_ordered_counterclockwise( mesh, tri)
    end do

    ! ti2
    tri = [locgeom%vj_focus, locgeom%vj_opp( locgeom%nvj_opp), locgeom%vj_clock]
    isso = isso .and. is_ordered_counterclockwise( mesh, tri)

  end function is_valid_order_vj_focus_nCge4

  pure function is_ordered_counterclockwise( mesh, tri) result( isso)
    !< Check if the vertices in a hypothetical triangle are ordered counter-clockwise

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    integer, dimension(3), intent(in   ) :: tri
    logical                              :: isso

    ! Local arguments:
    real(dp), dimension(2) :: pa, pb, pc

    pa = mesh%V( tri(1),:)
    pb = mesh%V( tri(2),:)
    pc = mesh%V( tri(3),:)

    isso = cross2( pb-pa, pc-pa) > 0._dp

  end function is_ordered_counterclockwise

end module delete_vertices_local_geometry
