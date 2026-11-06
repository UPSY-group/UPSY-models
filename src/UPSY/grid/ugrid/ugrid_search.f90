submodule(ugrid_basic) ugrid_search

contains

  subroutine find_containing_face( self, p, face_index)

    ! In/output variables:
    class(type_ugrid),      intent(in   ) :: self
    real(dp), dimension(2), intent(in   ) :: p
    integer,                intent(inout) :: face_index

    ! Local variables:

    ! Perform a linear search to find a ugrid face that either
    ! contains p, or is at least very close to it.
    call find_approximate_containing_face_linear_search( self, p, face_index)

    ! Then perform a flood-fill search outward from that initial guess
    ! to get the exact answer
    call find_containing_face_floodfill( self, p, face_index)

  end subroutine find_containing_face

  subroutine find_approximate_containing_face_linear_search( self, p, face_index)
    ! Use a linear search to find the "approximate" face of the ugrid
    ! containing the point p. "Approximate", because a linear search only
    ! is guaranteed to return the correct answer when the faces of the!
    ! ugrid are the Voronoi cells of the face coordinates (so for ugrids
    ! based on an UPSY mesh a-grid, or on a square x/y-grid). For e.g.
    ! a ugrid based on an UPSY mesh b-grid (triangles), it might be that
    ! the triangle whose geometric center is closest to p, does not actually
    ! contain p - although it will be very close. So this is still
    ! a good way to get an accurate initial guess.

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
    call crash('self%find_containing_face - couldnt find it!')

  end subroutine find_approximate_containing_face_linear_search

  subroutine find_containing_face_floodfill( self, p, face_index)

    ! In/output variables:
    class(type_ugrid),      intent(in   ) :: self
    real(dp), dimension(2), intent(in   ) :: p
    integer,                intent(inout) :: face_index

    ! Local variables:
    integer, dimension(self%nFaces) :: map, stack
    integer                          :: stackN
    integer                          :: ncycle
    integer                          :: ci, neighbouring_face_index

    ! Initialise the map and stack
    map = 0
    stackN = 1
    stack( stackN) = face_index
    map( face_index) = 1

    ncycle = 0
    do while (stackN > 0)

      ! Safety
      ncycle = ncycle + 1
      if (ncycle > self%nFaces) call crash('find_containing_face_floodfill got stuck')

      ! Take the last face from the stack
      face_index = stack( stackN)
      stackN = stackN - 1

      ! Mark it on the map
      map( face_index) = 2

      if (self%is_in_face( p, face_index)) then
      ! If this is the face that contains p, we've finished
        exit
      else
        ! If not, add its unchecked neighbours to the stack
        do ci = 1, self%face_face_connectivity_n( face_index)
          neighbouring_face_index = self%face_face_connectivity( face_index, ci)
          if (map( neighbouring_face_index) == 0) then
            map( neighbouring_face_index) = 1
            stackN = stackN + 1
            stack( stackN) = neighbouring_face_index
          end if
        end do
      end if

    end do

  end subroutine find_containing_face_floodfill

  function is_in_face( self, p, face_index) result( isso)
    ! Check if the point p lies inside face fi of the ugrid

    ! In/output variables:
    class(type_ugrid),      intent(in   ) :: self
    real(dp), dimension(2), intent(in   ) :: p
    integer,                intent(inout) :: face_index
    logical                               :: isso

    ! Local variables:
    real(dp), dimension(:,:), allocatable :: Vpoly
    integer                               :: n, ci, node_index

    n = self%face_node_connectivity_n( face_index)
    allocate( Vpoly( n,2))
    do ci = 1, n
      node_index = self%face_node_connectivity( face_index, ci)
      Vpoly( ci,:) = self%node_coordinates( node_index,:)
    end do

    isso = is_in_polygon( Vpoly, p)

  end function is_in_face

end submodule ugrid_search
