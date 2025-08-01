module flip_triangles

  ! Iteratively flip triangle pairs until the global Delaunay criterion is satisfied

  use tests_main
  use assertions_basic
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use check_Delaunay_criterion, only: are_Delaunay
  use mesh_utilities, only: write_mesh_to_text_file, update_triangle_circumcenter, &
    add_triangle_to_refinement_stack_first, find_triangle_pair_local_geometry

  implicit none

  private

  public :: flip_triangles_until_Delaunay, add_triangle_pairs_around_triangle_to_Delaunay_check_stack, &
    add_triangle_pairs_around_vertex_to_Delaunay_check_stack, initialise_Delaunay_check_stack, &
    add_triangle_pair_to_Delaunay_check_stack

contains

  subroutine flip_triangles_until_Delaunay( mesh)
    ! Iteratively flip triangle pairs until the local Delaunay
    ! criterion is satisfied everywhere

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'flip_triangles_until_Delaunay'
    integer                        :: ti, ni, tj, nj

    ! Add routine to path
    call init_routine( routine_name)

    do while (mesh%check_Delaunay_stackN > 0)

      ! Take the last triangle pair from the Delaunay checking stack
      ti = mesh%check_Delaunay_stack( mesh%check_Delaunay_stackN,1)
      ni = mesh%check_Delaunay_stack( mesh%check_Delaunay_stackN,2)
      tj = mesh%check_Delaunay_stack( mesh%check_Delaunay_stackN,3)
      nj = mesh%check_Delaunay_stack( mesh%check_Delaunay_stackN,4)

      mesh%check_Delaunay_stack( mesh%check_Delaunay_stackN,:) = 0
      mesh%check_Delaunay_stackN = mesh%check_Delaunay_stackN - 1
      mesh%check_Delaunay_map( ti,ni) = .false.
      mesh%check_Delaunay_map( tj,nj) = .false.

#if (DO_ASSERTIONS)
      ! Safety - assert that ti and tj are valid triangles
      call assert( test_ge_le( ti, 1, mesh%nTri), 'invalid value for ti')
      call assert( test_ge_le( tj, 1, mesh%nTri), 'invalid value for tj')
#endif

      ! If triangle pair [ti,tj does not meet the local Delaunay criterion,
      ! flip them, and add mark the new resulting pairs for (re)checking
      if (.not. are_Delaunay( mesh, ti, tj)) then
        call flip_triangle_pair( mesh, ti, tj)
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_triangles_until_Delaunay

  subroutine flip_triangle_pair( mesh, ti, tj)
    ! Flip the triangle pair ti-tj, supposedly because it doesn't meet the
    ! local Delaunay criterion
    !
    ! When going in, the local geometry looks like this:
    !
    !   \ /           \ /           \ /
    ! - -o----------- vic ---------- o- -
    !   / \           / \           / \
    !      \  tib    /   \   tia   /
    !       \       /     \       /
    !        \     /       \     /
    !         \   /   ti    \   /
    !          \ /           \ /
    !      - - via -------- vib - -
    !          / \           / \
    !         /   \         /   \
    !        /     \  tj   /     \
    !       /       \     /       \
    !      /   tjb   \   /   tja   \
    !   \ /           \ /           \ /
    ! - -o----------- vid ---------- o- -
    !   / \           / \           / \
    !
    ! When coming out, it looks like this:
    !
    !   \ /           \ /           \ /
    ! - -o----------- vic ---------- o- -
    !   / \           /|\           / \
    !      \  tib    / | \   tia   /
    !       \       /  |  \       /
    !        \     /   |   \     /
    !         \   / t1 | t2 \   /
    !          \ /     |     \ /
    !      - - via     |     vib - -
    !          / \     |     / \
    !         /   \    |    /   \
    !        /     \   |   /     \
    !       /       \  |  /       \
    !      /   tjb   \ | /   tja   \
    !   \ /           \|/           \ /
    ! - -o----------- vid ---------- o- -
    !   / \           / \           / \

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh
    integer,         intent(in   ) :: ti,tj

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'flip_triangle_pair'
    integer                        :: n
    integer                        :: via, vib, vic, vid
    integer                        :: ci, iti
    integer                        :: tia, tib, tja, tjb, t1, t2
    integer                        :: li_min, li_max
    logical                        :: foundit

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Safety - assert that ti and tj are valid triangles
    call assert( test_ge_le( ti, 1, mesh%nTri), 'invalid value for ti')
    call assert( test_ge_le( tj, 1, mesh%nTri), 'invalid value for tj')
    ! Safety - assert that ti and tj are neighbours
    call assert( test_mesh_triangles_are_neighbours( mesh, ti, tj), 'ti and tj are not meighbours')
#endif

    call remove_triangle_pairs_around_ti_from_Delaunay_check_stack( mesh, ti)
    call remove_triangle_pairs_around_ti_from_Delaunay_check_stack( mesh, tj)

    ! Determine the local geometry around triangle pair [ti,tj]
    call find_triangle_pair_local_geometry( mesh, ti, tj, via, vib, vic, vid, tia, tib, tja, tjb)

    ! == V, Tri

    ! Let triangle t1 be spanned by [via, vid, vic]
    t1 = ti
    mesh%Tri( t1,:) = [via, vid, vic]

    ! Let triangle t2 be spanned by [vib, vic, vid]
    t2 = tj
    mesh%Tri( t2,:) = [vib, vic, vid]

#if (DO_ASSERTIONS)
    ! Safety - check if everything went alright and we didnt create any duplicate triangles
    call assert( test_mesh_triangle_doesnt_have_duplicates( mesh, t1), 'a triangle with the vertices of t1 already exists')
    call assert( test_mesh_triangle_doesnt_have_duplicates( mesh, t2), 'a triangle with the vertices of t2 already exists')
#endif

    ! == nC, C

    ! via: connection to vib is removed
    do ci = 1, mesh%nC( via)
      if (mesh%C( via,ci) == vib) then
        mesh%C( via,:) = [mesh%C( via,1:ci-1), mesh%C( via,ci+1:mesh%nC_mem), 0]
        mesh%nC( via) = mesh%nC( via) - 1
        exit
      end if
    end do
    ! vib: connection to via is removed
    do ci = 1, mesh%nC( vib)
      if (mesh%C( vib,ci) == via) then
        mesh%C( vib,:) = [mesh%C( vib,1:ci-1), mesh%C( vib,ci+1:mesh%nC_mem), 0]
        mesh%nC( vib) = mesh%nC( vib) - 1
        exit
      end if
    end do
    ! vic: a connection to vid is added between via and vib
    mesh%nC( vic) = mesh%nC( vic) + 1
    do ci = 1, mesh%nC( vic)
      if (mesh%C( vic,ci) == via) then
        mesh%C( vic,:) = [mesh%C( vic,1:ci), vid, mesh%C( vic,ci+1:mesh%nC_mem-1)]
        exit
      end if
    end do
    ! vid: a connection to vic is added between vib and via
    mesh%nC( vid) = mesh%nC( vid) + 1
    do ci = 1, mesh%nC( vid)
      if (mesh%C( vid,ci) == vib) then
        mesh%C( vid,:) = [mesh%C( vid,1:ci), vic, mesh%C( vid,ci+1:mesh%nC_mem-1)]
        exit
      end if
    end do

    ! == niTri, iTri

    ! via: tj,ti are replaced by t1
    ! First remove ti
    do iti = 1, mesh%niTri( via)
      if (mesh%iTri( via,iti) == ti) then
        mesh%iTri( via,:) = [mesh%iTri( via,1:iti-1), mesh%iTri( via,iti+1:mesh%nC_mem), 0]
        mesh%niTri( via) = mesh%niTri( via) - 1
        exit
      end if
    end do
    ! then replace tj by t1
    do iti = 1, mesh%niTri( via)
      if (mesh%iTri( via,iti) == tj) then
        mesh%iTri( via,iti) = t1
        exit
      end if
    end do

    ! vib: ti,tj are replaced by t2
    ! First remove tj
    do iti = 1, mesh%niTri( vib)
      if (mesh%iTri( vib,iti) == tj) then
        mesh%iTri( vib,:) = [mesh%iTri( vib,1:iti-1), mesh%iTri( vib,iti+1:mesh%nC_mem), 0]
        mesh%niTri( vib) = mesh%niTri( vib) - 1
        exit
      end if
    end do
    ! then replace ti by t2
    do iti = 1, mesh%niTri( vib)
      if (mesh%iTri( vib,iti) == ti) then
        mesh%iTri( vib,iti) = t2
        exit
      end if
    end do

    ! vic: ti is replaced by t1,t2
    mesh%niTri( vic) = mesh%niTri( vic) + 1
    do iti = 1, mesh%niTri( vic)
      if (mesh%iTri( vic,iti) == ti) then
        mesh%iTri( vic,:) = [mesh%iTri( vic,1:iti-1), t1, t2, mesh%iTri( vic,iti+1:mesh%nC_mem-1)]
        exit
      end if
    end do

    ! vid: tj is replaced by t2,t1
    mesh%niTri( vid) = mesh%niTri( vid) + 1
    do iti = 1, mesh%niTri( vid)
      if (mesh%iTri( vid,iti) == tj) then
        foundit = .true.
        mesh%iTri( vid,:) = [mesh%iTri( vid,1:iti-1), t2, t1, mesh%iTri( vid,iti+1:mesh%nC_mem-1)]
        exit
      end if
    end do

    ! == Border index

    ! No changes.

    ! == TriC

    ! tia: ti is replaced by t2
    if (tia > 0) then
      do n = 1, 3
        if (mesh%TriC( tia,n) == ti) then
          mesh%TriC( tia,n) = t2
          exit
        end if
      end do
    end if
    ! tib: ti is replaced by t1
    if (tib > 0) then
      do n = 1, 3
        if (mesh%TriC( tib,n) == ti) then
          mesh%TriC( tib,n) = t1
          exit
        end if
      end do
    end if
    ! tja: tj is replaced by t2
    if (tja > 0) then
      do n = 1, 3
        if (mesh%TriC( tja,n) == tj) then
          mesh%TriC( tja,n) = t2
          exit
        end if
      end do
    end if
    ! tjb: tj is replaced by t1
    if (tjb > 0) then
      do n = 1, 3
        if (mesh%TriC( tjb,n) == tj) then
          mesh%TriC( tjb,n) = t1
          exit
        end if
      end do
    end if
    ! t1 is adjacent to [t2, tib, tjb]
    mesh%TriC( t1,:) = [t2, tib, tjb]
    ! t2 is adjacent to [t1, tja, tia]
    mesh%TriC( t2,:) = [t1, tja, tia]

    ! == Tricc

    call update_triangle_circumcenter( mesh, t1)
    call update_triangle_circumcenter( mesh, t2)

    ! == Refinement data

    ! Add the two new triangles to the refinement map and stack
    call add_triangle_to_refinement_stack_first( mesh, t1)
    call add_triangle_to_refinement_stack_first( mesh, t2)

    ! Update triangle-line overlap ranges
    li_min = min( mesh%Tri_li( ti,1), mesh%Tri_li( tj,1))
    li_max = max( mesh%Tri_li( ti,2), mesh%Tri_li( tj,2))
    mesh%Tri_li( t1,:) = [li_min, li_max]
    mesh%Tri_li( t2,:) = [li_min, li_max]

    ! Re-mark the changed triangle pairs for Delaunay checking
    if (tia > 0) call add_triangle_pair_to_Delaunay_check_stack( mesh, tia, t2)
    if (tib > 0) call add_triangle_pair_to_Delaunay_check_stack( mesh, tib, t1)
    if (tja > 0) call add_triangle_pair_to_Delaunay_check_stack( mesh, tja, t2)
    if (tjb > 0) call add_triangle_pair_to_Delaunay_check_stack( mesh, tjb, t1)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine flip_triangle_pair

  subroutine remove_triangle_pairs_around_ti_from_Delaunay_check_stack( mesh, ti_remove)

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh
    integer,         intent(in   ) :: ti_remove

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remove_triangle_pairs_around_ti_from_Delaunay_check_stack'
    integer                        :: i, ti, ni, tj, nj, j

    ! Add routine to path
    call init_routine( routine_name)

    i = 1
    do while (i <= mesh%check_Delaunay_stackN)
      if (mesh%check_Delaunay_stack( i,1) == ti_remove .or. &
          mesh%check_Delaunay_stack( i,3) == ti_remove) then

        ti = mesh%check_Delaunay_stack( i,1)
        ni = mesh%check_Delaunay_stack( i,2)
        tj = mesh%check_Delaunay_stack( i,3)
        nj = mesh%check_Delaunay_stack( i,4)

        ! Remove
        do j = 1, 4
          mesh%check_Delaunay_stack( :,j) = &
            [mesh%check_Delaunay_stack( 1  :i-1      ,j), &
             mesh%check_Delaunay_stack( i+1:mesh%nTri,j), &
             0]
        end do
        mesh%check_Delaunay_stackN = mesh%check_Delaunay_stackN - 1

        mesh%check_Delaunay_map( ti,ni) = .false.
        mesh%check_Delaunay_map( tj,nj) = .false.

      else
        i = i+1
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remove_triangle_pairs_around_ti_from_Delaunay_check_stack

  subroutine add_triangle_pairs_around_vertex_to_Delaunay_check_stack( mesh, vi)

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh
    integer,         intent(in   ) :: vi

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_triangle_pairs_around_vertex_to_Delaunay_check_stack'
    integer                        :: iti, ti

    ! Add routine to path
    call init_routine( routine_name)

    do iti = 1, mesh%niTri( vi)
      ti = mesh%iTri( vi,iti)
      call add_triangle_pairs_around_triangle_to_Delaunay_check_stack( mesh, ti)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_triangle_pairs_around_vertex_to_Delaunay_check_stack

  subroutine add_triangle_pairs_around_triangle_to_Delaunay_check_stack( mesh, ti)

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh
    integer,         intent(in   ) :: ti

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_triangle_pairs_around_triangle_to_Delaunay_check_stack'
    integer                        :: n, tj

    ! Add routine to path
    call init_routine( routine_name)

    do n = 1, 3
      tj = mesh%TriC( ti,n)
      if (tj > 0) then
        call add_triangle_pair_to_Delaunay_check_stack( mesh, ti, tj)
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_triangle_pairs_around_triangle_to_Delaunay_check_stack

  subroutine add_triangle_pair_to_Delaunay_check_stack( mesh, ti, tj)

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh
    integer,         intent(in   ) :: ti, tj

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_triangle_pair_to_Delaunay_check_stack'
    logical                        :: found_it
    integer                        :: n, ni, nj

    ! Add routine to path
    call init_routine( routine_name)

    ! Find ni such that TriC( ti,ni) = tj
    found_it = .false.
    ni = 0
    do n = 1, 3
      if (mesh%TriC( ti,n) == tj) then
        found_it = .true.
        ni = n
        exit
      end if
    end do
    if (.not. found_it) call crash('couldnt find connection from ti to tj')

    ! Find nj such that TriC( tj,nj) = ti
    found_it = .false.
    nj = 0
    do n = 1, 3
      if (mesh%TriC( tj,n) == ti) then
        found_it = .true.
        nj = n
        exit
      end if
    end do
    if (.not. found_it) call crash('couldnt find connection from tj to ti')

    ! If this pair has not yet been checked and is not yet marked for checking,
    ! add it to the to-be-checked stack
    if (.not. mesh%check_Delaunay_map( ti,ni)) then
      mesh%check_Delaunay_stackN = mesh%check_Delaunay_stackN + 1;
      mesh%check_Delaunay_stack( mesh%check_Delaunay_stackN,:) = [ti,ni,tj,nj]

      mesh%check_Delaunay_map( ti,ni) = .true.
      mesh%check_Delaunay_map( tj,nj) = .true.
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_triangle_pair_to_Delaunay_check_stack

  subroutine initialise_Delaunay_check_stack( mesh)

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_Delaunay_check_stack'

    ! Add routine to path
    call init_routine( routine_name)

    mesh%check_Delaunay_map    = .false.
    mesh%check_Delaunay_stack  = 0
    mesh%check_Delaunay_stackN = 0

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_Delaunay_check_stack

end module flip_triangles
