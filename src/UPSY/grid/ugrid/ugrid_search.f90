submodule(ugrid_basic) ugrid_search

contains

  subroutine find_nearest_face_center( self, p, face_index)
    ! Use a linear search to find the ugrid face whose center
    ! is nearest to p. Note that this is not necessarily
    ! the face containing p!

    ! In/output variables:
    class(type_ugrid),      intent(in   ) :: self
    real(dp), dimension(2), intent(in   ) :: p
    integer,                intent(inout) :: face_index

    ! Local variables:
    integer  :: ncycle
    integer  :: face_index_prev, vcmin, ci, neighbouring_face_index
    real(dp) :: d, dc, dcmin

    ncycle = 0
    face_index_prev = face_index
    do while (ncycle < self%nFaces)

      d = hypot( &
        self%face_coordinates( face_index,1) - p(1), &
        self%face_coordinates( face_index,2) - p(2))

      dcmin = d*2
      vcmin = 0
      do ci = 1, self%face_face_connectivity_n( face_index)
        neighbouring_face_index = self%face_face_connectivity( face_index, ci)
        if (neighbouring_face_index == face_index_prev) cycle ! This is the neighbour we just came from
        dc = hypot( &
          self%face_coordinates( neighbouring_face_index,1) - p(1), &
          self%face_coordinates( neighbouring_face_index,2) - p(2))
        if (dc < dcmin) then
          dcmin = dc
          vcmin = neighbouring_face_index
        end if
      end do

      if (dcmin < d) then
        face_index_prev = face_index
        face_index = vcmin
      else
        return
      end if

    end do

    ! If we reach this point, we didnt find the containing vertex - should not be possible, so throw an error
    call crash('self%find_nearest_face_center - couldnt find it!')

  end subroutine find_nearest_face_center

  subroutine get_face_as_polygon( self, face_index, poly)
    ! Return the nodes of the polyogn spanning a face of the ugrid

    ! In/output variables:
    class(type_ugrid),                     intent(in   ) :: self
    integer,                               intent(in   ) :: face_index
    real(dp), dimension(:,:), allocatable, intent(inout) :: poly

    ! Local variables:
    integer :: n, ci, node_index

    if (allocated( poly)) deallocate( poly)

    n = self%face_node_connectivity_n( face_index)

    allocate( poly( n,2))

    do ci = 1, n
      node_index = self%face_node_connectivity( face_index, ci)
      poly( ci,:) = self%node_coordinates( node_index,:)
    end do

  end subroutine get_face_as_polygon

end submodule ugrid_search
