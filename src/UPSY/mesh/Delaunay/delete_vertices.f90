module delete_vertices

  ! Delete a vertex from the mesh and update the Delaunay triangulation accordingly.

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use mesh_types, only: type_mesh
  use assertions_basic, only: assert
  use tests_main
  use switch_array_elements, only: switch_rows
  use switch_vertices_triangles, only: switch_vertices, switch_triangles
  use mesh_utilities, only: update_triangle_circumcenter
  use flip_triangles, only: add_triangle_pairs_around_triangle_to_Delaunay_check_stack, &
    flip_triangles_until_Delaunay
  use mpi_basic, only: par
  use delete_vertices_local_geometry, only: type_local_geometry_nCge4, find_local_geometry_nCge4
  use basic_mesh_data_operations, only: replace_vj_in_C_vi_with_vks, remove_vj_in_C_vi, &
    remove_ti_in_iTri_vi, replace_ti_in_iTri_vi_with_tj, add_tjs_in_iTri_vi_after_ti, &
    replace_vi_in_Tri_ti_with_vj, replace_tj_in_TriC_ti_with_tk

  implicit none

  private

  public :: delete_vertex

contains

  subroutine delete_vertex( mesh, vi_kill, vi_new2vi_old, vi_old2vi_new)
    ! Delete vertex vi from the mesh

    ! In/output variables:
    type(type_mesh),                    intent(inout) :: mesh
    integer,                            intent(in   ) :: vi_kill
    integer, dimension(:), allocatable, intent(  out) :: vi_new2vi_old, vi_old2vi_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex'
    integer                        :: ci, vj
    logical                        :: is_double_free

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    is_double_free = .true.
    if (mesh%VBI( vi_kill) > 0) is_double_free = .false.
    do ci = 1, mesh%nC( vi_kill)
      vj = mesh%C( vi_kill,ci)
      if (mesh%VBI( vj) > 0) is_double_free = .false.
    end do
    if (.not. is_double_free) call crash('delete_vertex only works on double-free vertices for now')

    if (mesh%nC( vi_kill) == 3) then
      call crash('delete_vertex_nCeq3 not implemented yet')
      ! call delete_vertex_nCeq3( mesh, vi_kill, vi_new2vi_old, vi_old2vi_new)
    else
      call delete_vertex_nCge4( mesh, vi_kill, vi_new2vi_old, vi_old2vi_new)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex

  ! Delete a vertex which has 4 or more neighbours
  ! ==============================================

  subroutine delete_vertex_nCge4( mesh, vi_kill, vi_new2vi_old, vi_old2vi_new)
    ! Before, local geometry looks like this:
    !
    !                tj_opp()
    !
    !  vj_opp(nvj_opp) ------ vj_opp(1)
    !           /  \ ti_opp() /  \
    !     tj1  /    \        /    \  tj4
    !         / ti1  \      /  ti4 \
    !        /        \    /        \
    !   vj_clock ---- vi_kill ---- vj_anti
    !         \          |          /
    !           \  ti2   |   ti3  /
    !             \      |      /
    !        tj2    \    |    /    tj3
    !                 vj_focus
    !
    ! After, it looks like this:
    !
    !                tj_opp()
    !
    !  vj_opp(nvj_opp) ------ vj_opp(1)
    !           /  \ ti_opp() /  \
    !     tj1  /    |        |    \  tj4
    !         /     \        /     \
    !        /  ti2  |      |  ti3  \
    !   vj_clock     \      /      vj_anti
    !         \       |    |        /
    !           \     \    /      /
    !             \    |  |     /
    !        tj2    \   \ /   /    tj3
    !                 vj_focus

    ! In/output variables:
    type(type_mesh),                    intent(inout) :: mesh
    integer,                            intent(in   ) :: vi_kill
    integer, dimension(:), allocatable, intent(  out) :: vi_new2vi_old, vi_old2vi_new

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'delete_vertex_nCge4'
    integer, dimension(:), allocatable :: ti_new2ti_old, ti_old2ti_new
    integer                            :: vi, ti
    type(type_local_geometry_nCge4)    :: locgeom
    integer                            :: vii, tii

    ! Add routine to path
    call init_routine( routine_name)

    allocate( vi_new2vi_old( mesh%nV))
    allocate( vi_old2vi_new( mesh%nV))
    do vi = 1, mesh%nV
      vi_new2vi_old( vi) = vi
      vi_old2vi_new( vi) = vi
    end do

    allocate( ti_new2ti_old( mesh%nTri))
    allocate( ti_old2ti_new( mesh%nTri))
    do ti = 1, mesh%nTri
      ti_new2ti_old( ti) = ti
      ti_old2ti_new( ti) = ti
    end do

    ! Determine the local geometry
    call find_local_geometry_nCge4( mesh, vi_kill, locgeom)

    ! Adapt all the connectivity data to the new geometry
    call delete_vertex_nCge4_nC_C      ( mesh, locgeom)
    call delete_vertex_nCge4_niTri_iTri( mesh, locgeom)
    call delete_vertex_nCge4_Tri       ( mesh, locgeom)
    call delete_vertex_nCge4_TriC      ( mesh, locgeom)

    ! Delete the now disconnected vertex and triangles from the lists
    call delete_vertex_V  ( mesh, locgeom%vi_kill, vi_new2vi_old, vi_old2vi_new, locgeom)
    call delete_vertex_Tri( mesh, locgeom%ti1    , ti_new2ti_old, ti_old2ti_new, locgeom)
    call delete_vertex_Tri( mesh, locgeom%ti4    , ti_new2ti_old, ti_old2ti_new, locgeom)

#if (DO_ASSERTIONS)
    call assert( test_ge_le( locgeom%vj_clock, 1, mesh%nV  ), 'invalid updated vj_clock')
    call assert( test_ge_le( locgeom%vj_focus, 1, mesh%nV  ), 'invalid updated vj_focus')
    call assert( test_ge_le( locgeom%vj_anti , 1, mesh%nV  ), 'invalid updated vj_anti')
    call assert( test_ge_le( locgeom%vj_opp  , 1, mesh%nV  ), 'invalid updated vj_opp')
    call assert( test_ge_le( locgeom%ti2     , 1, mesh%nTri), 'invalid updated ti2')
    call assert( test_ge_le( locgeom%ti3     , 1, mesh%nTri), 'invalid updated ti3')
    call assert( test_ge_le( locgeom%ti_opp  , 1, mesh%nTri), 'invalid updated ti_opp')
#endif

    ! Calculate circumcenters for the new triangles
    call update_triangle_circumcenter( mesh, locgeom%ti2)
    call update_triangle_circumcenter( mesh, locgeom%ti3)
    do tii = 1, locgeom%nti_opp
      call update_triangle_circumcenter( mesh, locgeom%ti_opp( tii))
    end do

    ! Iteratively flip triangle pairs until the local Delaunay
    ! criterion is satisfied everywhere

    mesh%check_Delaunay_map    = .false.
    mesh%check_Delaunay_stackN = 0

    call add_triangle_pairs_around_triangle_to_Delaunay_check_stack( mesh, locgeom%ti2)
    call add_triangle_pairs_around_triangle_to_Delaunay_check_stack( mesh, locgeom%ti3)
    do tii = 1, locgeom%nti_opp
      call add_triangle_pairs_around_triangle_to_Delaunay_check_stack( mesh, locgeom%ti_opp( tii))
    end do

    call flip_triangles_until_Delaunay( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_nCge4

  subroutine delete_vertex_nCge4_nC_C( mesh, locgeom)

    ! In/output variables:
    type(type_mesh),                 intent(inout) :: mesh
    type(type_local_geometry_nCge4), intent(in   ) :: locgeom

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_nCge4_nC_C'
    integer                        :: vii

    ! Add routine to path
    call init_routine( routine_name)

    ! vj_focus: replace connection to vi_kill by vj_opp(:)
    call replace_vj_in_C_vi_with_vks( mesh, locgeom%vj_focus, locgeom%vi_kill, locgeom%vj_opp)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%C( locgeom%vj_focus,:) == locgeom%vi_kill), 'vi_kill is still listed in C( vj_focus,:)')
    do vii = 1, locgeom%nvj_opp
      call assert( any( mesh%C( locgeom%vj_focus,:) == locgeom%vj_opp( vii)), 'vj_opp is not listed in C( vj_focus,:)')
    end do
#endif

    ! vj_clock: remove connection to vi_kill
    call remove_vj_in_C_vi( mesh, locgeom%vj_clock, locgeom%vi_kill)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%C( locgeom%vj_clock,:) == locgeom%vi_kill), 'vi_kill is still listed in C( vj_clock,:)')
#endif

    ! vj_anti: remove connection to vi_kill
    call remove_vj_in_C_vi( mesh, locgeom%vj_anti, locgeom%vi_kill)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%C( locgeom%vj_anti,:) == locgeom%vi_kill), 'vi_kill is still listed in C( vj_anti,:)')
#endif

    ! vj_opp: replace connection to vi_kill by vj_focus
    do vii = 1, locgeom%nvj_opp
      call replace_vj_in_C_vi_with_vks( mesh, locgeom%vj_opp( vii), locgeom%vi_kill, [locgeom%vj_focus])
#if (DO_ASSERTIONS)
      call assert( .not. any( mesh%C( locgeom%vj_opp( vii),:) == locgeom%vi_kill ), 'vi_kill is still listed in C( vj_opp,:)')
      call assert(       any( mesh%C( locgeom%vj_opp( vii),:) == locgeom%vj_focus), 'vj_focus is not listed in C( vj_opp,:)')
#endif
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_nCge4_nC_C

  subroutine delete_vertex_nCge4_niTri_iTri( mesh, locgeom)

    ! In/output variables:
    type(type_mesh),                 intent(inout) :: mesh
    type(type_local_geometry_nCge4), intent(in   ) :: locgeom

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_nCge4_niTri_iTri'
    integer                        :: tii

    ! Add routine to path
    call init_routine( routine_name)

    ! vj_focus: add ti_opp(:) after ti3
    call add_tjs_in_iTri_vi_after_ti( mesh, locgeom%vj_focus, locgeom%ti3, locgeom%ti_opp)
#if (DO_ASSERTIONS)
    do tii = 1, locgeom%nti_opp
      call assert( any( mesh%iTri( locgeom%vj_focus,:) == locgeom%ti_opp( tii)), 'ti_opp is not listed in iTri( vj_focus,:)')
    end do
#endif

    ! vj_clock: remove ti1 from iTri
    call remove_ti_in_iTri_vi( mesh, locgeom%vj_clock, locgeom%ti1)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%iTri( locgeom%vj_clock,:) == locgeom%ti1), 'ti1 is still listed in iTri( vj_clock,:)')
#endif

    ! vj_anti: remove ti4 from iTri
    call remove_ti_in_iTri_vi( mesh, locgeom%vj_anti, locgeom%ti4)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%iTri( locgeom%vj_anti,:) == locgeom%ti4), 'ti4 is still listed in iTri( vj_anti,:)')
#endif

    ! vj_opp: replace ti1 with ti2, and ti4 with ti3
    call replace_ti_in_iTri_vi_with_tj( mesh, locgeom%vj_opp(locgeom%nvj_opp), locgeom%ti1, locgeom%ti2)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%iTri( locgeom%vj_opp( locgeom%nvj_opp),:) == locgeom%ti1), 'ti1 is still listed in iTri( vj_opp(end),:)')
    call assert(       any( mesh%iTri( locgeom%vj_opp( locgeom%nvj_opp),:) == locgeom%ti2), 'ti2 is not listed in iTri( vj_opp(end),:)')
#endif
    call replace_ti_in_iTri_vi_with_tj( mesh, locgeom%vj_opp( 1             ), locgeom%ti4, locgeom%ti3)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%iTri( locgeom%vj_opp( 1),:) == locgeom%ti4), 'ti1 is still listed in iTri( vj_opp(1),:)')
    call assert(       any( mesh%iTri( locgeom%vj_opp( 1),:) == locgeom%ti3), 'ti2 is not listed in iTri( vj_opp(1),:)')
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_nCge4_niTri_iTri

  subroutine delete_vertex_nCge4_Tri( mesh, locgeom)

    ! In/output variables:
    type(type_mesh),                 intent(inout) :: mesh
    type(type_local_geometry_nCge4), intent(in   ) :: locgeom

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_nCge4_Tri'
    integer                        :: tii

    ! Add routine to path
    call init_routine( routine_name)

    ! ti1: set to -1 (will be deleted)
    mesh%Tri( locgeom%ti1,:) = -1

    ! ti2: replace vi_kill with vj_opp( end)
    call replace_vi_in_Tri_ti_with_vj( mesh, locgeom%ti2, locgeom%vi_kill, locgeom%vj_opp( locgeom%nvj_opp))
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%Tri( locgeom%ti2,:) == locgeom%vi_kill), 'vi_kill is still listed in Tri( ti2,:)')
    call assert(       any( mesh%Tri( locgeom%ti2,:) == locgeom%vj_opp( locgeom%nvj_opp)), 'vj_opp(end) is not listed in Tri( ti2,:)')
#endif

    ! ti3: replace vi_kill with vj_opp( 1)
    call replace_vi_in_Tri_ti_with_vj( mesh, locgeom%ti3, locgeom%vi_kill, locgeom%vj_opp( 1))
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%Tri( locgeom%ti3,:) == locgeom%vi_kill), 'vi_kill is still listed in Tri( ti3,:)')
    call assert(       any( mesh%Tri( locgeom%ti3,:) == locgeom%vj_opp( 1)), 'vj_opp(1) is not listed in Tri( ti3,:)')
#endif

    ! ti4: set to -1 (will be deleted)
    mesh%Tri( locgeom%ti4,:) = -1

    ! ti_opp: replace vi_kill with vj_focus
    do tii = 1, locgeom%nti_opp
      call replace_vi_in_Tri_ti_with_vj( mesh, locgeom%ti_opp( tii), locgeom%vi_kill, locgeom%vj_focus)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%Tri( locgeom%ti_opp( tii),:) == locgeom%vi_kill), 'vi_kill is still listed in Tri( ti_opp,:)')
    call assert(       any( mesh%Tri( locgeom%ti_opp( tii),:) == locgeom%vj_focus), 'vj_focus is not listed in Tri( ti_opp,:)')
#endif
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_nCge4_Tri

  subroutine delete_vertex_nCge4_TriC( mesh, locgeom)

    ! In/output variables:
    type(type_mesh),                 intent(inout) :: mesh
    type(type_local_geometry_nCge4), intent(in   ) :: locgeom

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_nCge4_TriC'

    ! Add routine to path
    call init_routine( routine_name)

    ! ti1: set to -1 (will be deleted)
    mesh%TriC( locgeom%ti1,:) = -1

    ! ti2: replace ti1 in TriC with tj1, and ti3 with ti_opp( end) if the latter exists
    call replace_tj_in_TriC_ti_with_tk( mesh, locgeom%ti2, locgeom%ti1, locgeom%tj1)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%TriC( locgeom%ti2,:) == locgeom%ti1), 'ti1 is still listed in TriC( ti2,:)')
    call assert(       any( mesh%TriC( locgeom%ti2,:) == locgeom%tj1), 'tj1 is not listed in TriC( ti2,:)')
#endif
    if (locgeom%nti_opp >= 1) then
      call replace_tj_in_TriC_ti_with_tk( mesh, locgeom%ti2, locgeom%ti3, locgeom%ti_opp( locgeom%nti_opp))
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%TriC( locgeom%ti2,:) == locgeom%ti3), 'ti3 is still listed in TriC( ti2,:)')
    call assert(       any( mesh%TriC( locgeom%ti2,:) == locgeom%ti_opp( locgeom%nti_opp)), 'ti_opp(end) is not listed in TriC( ti2,:)')
#endif
    end if

    ! ti3: replace ti4 in TriC with tj4, and ti2 with ti_opp( 1) if the latter exists
    call replace_tj_in_TriC_ti_with_tk( mesh, locgeom%ti3, locgeom%ti4, locgeom%tj4)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%TriC( locgeom%ti3,:) == locgeom%ti4), 'ti4 is still listed in TriC( ti3,:)')
    call assert(       any( mesh%TriC( locgeom%ti3,:) == locgeom%tj4), 'tj4 is not listed in TriC( ti3,:)')
#endif
    if (locgeom%nti_opp >= 1) then
      call replace_tj_in_TriC_ti_with_tk( mesh, locgeom%ti3, locgeom%ti2, locgeom%ti_opp( 1))
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%TriC( locgeom%ti3,:) == locgeom%ti2), 'ti2 is still listed in TriC( ti3,:)')
    call assert(       any( mesh%TriC( locgeom%ti3,:) == locgeom%ti_opp( 1)), 'ti_opp(1) is not listed in TriC( ti3,:)')
#endif
    end if

    ! ti4: set to -1 (will be deleted)
    mesh%TriC( locgeom%ti4,:) = -1

    ! if ti_opp( 1) exists: replace ti4 in TriC with ti3
    if (locgeom%nti_opp >= 1) then
      call replace_tj_in_TriC_ti_with_tk( mesh, locgeom%ti_opp( 1), locgeom%ti4, locgeom%ti3)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%TriC( locgeom%ti_opp( 1),:) == locgeom%ti4), 'ti4 is still listed in TriC( ti_opp( 1),:)')
    call assert(       any( mesh%TriC( locgeom%ti_opp( 1),:) == locgeom%ti3), 'ti3 is not listed in TriC( ti_opp( 1),:)')
#endif
    end if

    ! if ti_opp( end) exists: replace ti1 in TriC with ti2
    if (locgeom%nti_opp >= 1) then
      call replace_tj_in_TriC_ti_with_tk( mesh, locgeom%ti_opp( locgeom%nti_opp), locgeom%ti1, locgeom%ti2)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%TriC( locgeom%ti_opp( locgeom%nti_opp),:) == locgeom%ti1), 'ti1 is still listed in TriC( ti_opp( end),:)')
    call assert(       any( mesh%TriC( locgeom%ti_opp( locgeom%nti_opp),:) == locgeom%ti2), 'ti2 is not listed in TriC( ti_opp( end),:)')
#endif
    end if

    ! tj1: replace ti1 in TriC with ti2
    call replace_tj_in_TriC_ti_with_tk( mesh, locgeom%tj1, locgeom%ti1, locgeom%ti2)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%TriC( locgeom%tj1,:) == locgeom%ti1), 'ti1 is still listed in TriC( tj1,:)')
    call assert(       any( mesh%TriC( locgeom%tj1,:) == locgeom%ti2), 'ti2 is not listed in TriC( tj1,:)')
#endif

    ! tj2: nothing changes
    ! tj3: nothing changes

    ! tj4: replace ti4 in TriC with ti3
    call replace_tj_in_TriC_ti_with_tk( mesh, locgeom%tj4, locgeom%ti4, locgeom%ti3)
#if (DO_ASSERTIONS)
    call assert( .not. any( mesh%TriC( locgeom%tj4,:) == locgeom%ti4), 'ti4 is still listed in TriC( tj4,:)')
    call assert(       any( mesh%TriC( locgeom%tj4,:) == locgeom%ti3), 'ti3 is not listed in TriC( tj4,:)')
#endif

    ! tj_opp: nothing changes

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_nCge4_TriC

  subroutine delete_vertex_V( mesh, vi_kill, vi_new2vi_old, vi_old2vi_new, locgeom)

    ! In/output variables:
    type(type_mesh),                    intent(inout) :: mesh
    integer,                            intent(in   ) :: vi_kill
    integer, dimension(:), allocatable, intent(inout) :: vi_new2vi_old, vi_old2vi_new
    type(type_local_geometry_nCge4),    intent(inout) :: locgeom

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_V'
    integer                        :: vii

    ! Add routine to path
    call init_routine( routine_name)

    call switch_vertices( mesh, vi_kill, mesh%nV)

    call switch_rows( vi_new2vi_old, vi_kill, mesh%nV)
    call switch_rows( vi_old2vi_new, vi_kill, mesh%nV)
    vi_new2vi_old = vi_new2vi_old( 1: mesh%nV)
    vi_old2vi_new( vi_kill) = 0

    mesh%nV    = mesh%nV - 1
    mesh%V     = mesh%V(     1:mesh%nV,:)
    mesh%nC    = mesh%nC(    1:mesh%nV  )
    mesh%C     = mesh%C(     1:mesh%nV,:)
    mesh%niTri = mesh%niTri( 1:mesh%nV  )
    mesh%iTri  = mesh%iTri(  1:mesh%nV,:)
    mesh%VBI   = mesh%VBI(   1:mesh%nV  )

    ! Update the local geometry
    locgeom%vi_kill  = -1
    locgeom%vj_clock = vi_old2vi_new( locgeom%vj_clock)
    locgeom%vj_focus = vi_old2vi_new( locgeom%vj_focus)
    locgeom%vj_anti  = vi_old2vi_new( locgeom%vj_anti )
    do vii = 1, locgeom%nvj_opp
      locgeom%vj_opp( vii) = vi_old2vi_new( locgeom%vj_opp( vii))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_V

  subroutine delete_vertex_Tri( mesh, ti_kill, ti_new2ti_old, ti_old2ti_new, locgeom)

    ! In/output variables:
    type(type_mesh),                    intent(inout) :: mesh
    integer,                            intent(in   ) :: ti_kill
    integer, dimension(:), allocatable, intent(inout) :: ti_new2ti_old, ti_old2ti_new
    type(type_local_geometry_nCge4),    intent(inout) :: locgeom

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_Tri'
    logical                        :: found_it
    integer                        :: tii

    ! Add routine to path
    call init_routine( routine_name)

    call switch_triangles( mesh, ti_kill, mesh%nTri)

    call switch_rows( ti_new2ti_old, ti_kill, mesh%nTri)
    call switch_rows( ti_old2ti_new, ti_kill, mesh%nTri)
    ti_new2ti_old = ti_new2ti_old( 1: mesh%nTri)
    ti_old2ti_new( ti_kill) = 0

    mesh%nTri  = mesh%nTri - 1
    mesh%Tri   = mesh%Tri(   1:mesh%nTri,:)
    mesh%Tricc = mesh%Tricc( 1:mesh%nTri,:)
    mesh%TriC  = mesh%TriC(  1:mesh%nTri,:)

    ! Update local geometry
    if (locgeom%ti1 == ti_kill    ) locgeom%ti1 = -1
    if (locgeom%ti1 == mesh%nTri+1) locgeom%ti1 = ti_kill
    if (locgeom%ti2 == ti_kill    ) locgeom%ti2 = -1
    if (locgeom%ti2 == mesh%nTri+1) locgeom%ti2 = ti_kill
    if (locgeom%ti3 == ti_kill    ) locgeom%ti3 = -1
    if (locgeom%ti3 == mesh%nTri+1) locgeom%ti3 = ti_kill
    if (locgeom%ti4 == ti_kill    ) locgeom%ti4 = -1
    if (locgeom%ti4 == mesh%nTri+1) locgeom%ti4 = ti_kill
    do tii = 1, locgeom%nti_opp
      if (locgeom%ti_opp( tii) == ti_kill    ) locgeom%ti_opp( tii) = -1
      if (locgeom%ti_opp( tii) == mesh%nTri+1) locgeom%ti_opp( tii) = ti_kill
      if (locgeom%tj_opp( tii) == ti_kill    ) locgeom%tj_opp( tii) = -1
      if (locgeom%tj_opp( tii) == mesh%nTri+1) locgeom%tj_opp( tii) = ti_kill
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_Tri

end module delete_vertices
