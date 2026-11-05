submodule( ugrid_basic) build_from_xy_grid_mod

contains

  subroutine build_from_xy_grid( ugrid, xy_grid, xy_grid_Arakawa_grid)

    ! In/output variables:
    type(type_ugrid),        intent(inout) :: ugrid
    type(type_xy_grid),      intent(in   ) :: xy_grid
    type(type_Arakawa_grid), intent(in   ) :: xy_grid_Arakawa_grid

    ! Local variables:
    character(len=*), parameter :: routine_name = 'build_from_xy_grid'

    ! Add routine to path
    call init_routine( routine_name)

    ! Build from the regular (a-grid) or staggered (b/cx/cy-grid) grid cells
    if     (xy_grid_Arakawa_grid == Arakawa_grid%a()) then
      ugrid%name = 'ugrid_xy_grid_' // trim( xy_grid%name) // '_a'
      call build_from_xy_grid_a( ugrid, xy_grid)
    elseif (xy_grid_Arakawa_grid == Arakawa_grid%b()) then
      ugrid%name = 'ugrid_xy_grid_' // trim( xy_grid%name) // '_b'
      call crash('fixme!')
    elseif (xy_grid_Arakawa_grid == Arakawa_grid%cx()) then
      ugrid%name = 'ugrid_xy_grid_' // trim( xy_grid%name) // '_cx'
      call crash('fixme!')
    elseif (xy_grid_Arakawa_grid == Arakawa_grid%cy()) then
      ugrid%name = 'ugrid_xy_grid_' // trim( xy_grid%name) // '_cy'
      call crash('fixme!')
    else
      call crash('unsupported Arakawa grid')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine build_from_xy_grid

  subroutine build_from_xy_grid_a( ugrid, xy_grid)
    ! Build a ugrid from the regular grid cells of a square x/y-grid

    ! In/output variables:
    type(type_ugrid),   intent(inout) :: ugrid
    type(type_xy_grid), intent(in   ) :: xy_grid

    ! Local variables:
    character(len=*), parameter :: routine_name = 'build_from_xy_grid_a'
    integer                     :: i,j,fi,ii,jj,ni,ni_sw,ni_se,ni_ne,ni_nw,fi_s,fi_n,fi_w,fi_e

    ! Add routine to path
    call init_routine( routine_name)

    ! Number of nodes and faces in the ugrid
    ugrid%nNodes = (xy_grid%nx + 1) * (xy_grid%ny + 1)
    ugrid%nFaces = (xy_grid%nx    ) * (xy_grid%ny    )

    ! x/y-coordinates of nodes and faces
    call allocate_dist_shared( ugrid%node_coordinates, ugrid%wnode_coordinates, [1, ugrid%nNodes], [1, 2])
    if (par%node_primary) then
      do ii = 1, xy_grid%nx+1
      do jj = 1, xy_grid%ny+1
        ni = xy_grid%iijj2ni( ii, jj)
        if (ii == xy_grid%nx+1) then
          ugrid%node_coordinates( ni,1) = xy_grid%x( xy_grid%nx) + xy_grid%dx / 2._dp
        else
          ugrid%node_coordinates( ni,1) = xy_grid%x( ii) - xy_grid%dx / 2._dp
        end if
        if (jj == xy_grid%ny+1) then
          ugrid%node_coordinates( ni,2) = xy_grid%y( xy_grid%ny) + xy_grid%dx / 2._dp
        else
          ugrid%node_coordinates( ni,2) = xy_grid%y( jj) - xy_grid%dx / 2._dp
        end if
      end do
      end do
    end if

    call allocate_dist_shared( ugrid%face_coordinates, ugrid%wface_coordinates, [1, ugrid%nFaces], [1, 2])
    if (par%node_primary) then
      do i = 1, xy_grid%nx
      do j = 1, xy_grid%ny
        fi = xy_grid%ij2fi( i, j)
        ugrid%face_coordinates( fi,1) = xy_grid%x( i)
        ugrid%face_coordinates( fi,2) = xy_grid%y( j)
      end do
      end do
    end if

    ! Face-node connectivity
    ugrid%face_node_connectivity_fill_value = -999
    ugrid%face_node_connectivity_n_max = 4

    call allocate_dist_shared( ugrid%face_node_connectivity_n, ugrid%wface_node_connectivity_n, &
      [1, ugrid%nFaces])
    call allocate_dist_shared( ugrid%face_node_connectivity, ugrid%wface_node_connectivity, &
      [1, ugrid%nFaces], [1, ugrid%face_node_connectivity_n_max])

    if (par%node_primary) then
      ugrid%face_node_connectivity_n = 4
      do i = 1, xy_grid%nx
      do j = 1, xy_grid%ny
        fi    = xy_grid%ij2fi  ( i  , j  )
        ni_sw = xy_grid%iijj2ni( i  , j  )
        ni_se = xy_grid%iijj2ni( i+1, j  )
        ni_nw = xy_grid%iijj2ni( i  , j+1)
        ni_ne = xy_grid%iijj2ni( i+1, j+1)
        ugrid%face_node_connectivity( fi,:) = [ni_sw, ni_se, ni_ne, ni_nw]
      end do
      end do
    end if

    ! Face-face connectivity
    ugrid%face_face_connectivity_fill_value = -999
    ugrid%face_face_connectivity_n_max = 4

    call allocate_dist_shared( ugrid%face_face_connectivity_n, ugrid%wface_face_connectivity_n, &
      [1, ugrid%nFaces])
    call allocate_dist_shared( ugrid%face_face_connectivity, ugrid%wface_face_connectivity, &
      [1, ugrid%nFaces], [1, ugrid%face_face_connectivity_n_max])

    if (par%node_primary) then
      do i = 1, xy_grid%nx
      do j = 1, xy_grid%ny
        fi = xy_grid%ij2fi( i, j)
        ugrid%face_face_connectivity_n( fi) = 0
        ugrid%face_face_connectivity( fi,:) = ugrid%face_face_connectivity_fill_value
        ! Eastern neighbouring face
        if (i < xy_grid%nx) then
          fi_e = xy_grid%ij2fi( i+1, j)
          ugrid%face_face_connectivity_n( fi) = ugrid%face_face_connectivity_n( fi) + 1
          ugrid%face_face_connectivity( fi, ugrid%face_face_connectivity_n( fi)) = fi_e
        end if
        ! Northern neighbouring face
        if (j < xy_grid%ny) then
          fi_n = xy_grid%ij2fi( i, j+1)
          ugrid%face_face_connectivity_n( fi) = ugrid%face_face_connectivity_n( fi) + 1
          ugrid%face_face_connectivity( fi, ugrid%face_face_connectivity_n( fi)) = fi_n
        end if
        ! Western neighbouring face
        if (i > 1) then
          fi_w = xy_grid%ij2fi( i-1, j)
          ugrid%face_face_connectivity_n( fi) = ugrid%face_face_connectivity_n( fi) + 1
          ugrid%face_face_connectivity( fi, ugrid%face_face_connectivity_n( fi)) = fi_w
        end if
        ! Southern neighbouring face
        if (j > 1) then
          fi_s = xy_grid%ij2fi( i, j-1)
          ugrid%face_face_connectivity_n( fi) = ugrid%face_face_connectivity_n( fi) + 1
          ugrid%face_face_connectivity( fi, ugrid%face_face_connectivity_n( fi)) = fi_s
        end if
      end do
      end do
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine build_from_xy_grid_a

end submodule build_from_xy_grid_mod
