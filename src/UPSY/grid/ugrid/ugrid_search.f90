submodule(ugrid_basic) ugrid_search

contains

  subroutine find_containing_face( self, p, face_index)

    ! In/output variables:
    class(type_ugrid),      intent(in   ) :: self
    real(dp), dimension(2), intent(in   ) :: p
    integer,                intent(inout) :: face_index

    ! Local variables:
    integer                          :: ncycles
    integer,  dimension(self%nFaces) :: map, stack
    integer                          :: stackN

    ! Safety
    if (p(1) < self%xmin .or. p(1) > self%xmax .or. p(2) < self%ymin .or. p(2) > self%ymax) then
      call crash('type_ugrid%find_containing_face - point outside of ugrid bounding box')
    end if

    ! Get an approximation by finding the face whose center is nearest to p
    call self%find_nearest_face_center( p, face_index)

    ! Search outward from there, floodfill-style
    map = 0
    stackN = 1
    stack( stackN) = face_index
    map( face_index) = 1

    ncycles = 0
    do while (stackN > 0)

      ! Safety
      ncycles = ncycles + 1
      if (ncycles > self%nFaces) call crash('type_ugrid%find_containing_face - floodfill-style outward search got stuck')

      ! Take the last face from the stack
      face_index = stack( stackN)
      stackN = stackN - 1
      map( face_index) = 2

      ! See if p lies inside of it
      if (self%is_in_face( face_index, p)) return

      ! If not, mark it as checked and search outward
      call self%add_nonchecked_neighbours_to_map_and_stack( map, stack, stackN, face_index)

    end do

    ! If we reach this point, we didnt find the containing face - should not be possible, so throw an error
    write(0,*) 'p = [', p(1),',', p(2), ']'
    call crash('type_ugrid%find_containing_face - couldnt find it!')

  end subroutine find_containing_face

  function is_in_face( self, face_index, p) result( isso)

    ! In/output variables:
    class(type_ugrid),      intent(in   ) :: self
    integer,                intent(in   ) :: face_index
    real(dp), dimension(2), intent(in   ) :: p
    logical                               :: isso

    ! Local variables:
    real(dp), dimension(:,:), allocatable :: poly

    call self%get_face_as_polygon( face_index, poly)
    isso = is_in_convex_polygon( poly, p, self%tol_dist)

  end function is_in_face

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

    ! Safety
    if (p(1) < self%xmin .or. p(1) > self%xmax .or. p(2) < self%ymin .or. p(2) > self%ymax) then
      call crash('type_ugrid%find_nearest_face_center - point outside of ugrid bounding box')
    end if

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

    ! If we reach this point, we didnt find the nearest face center - should not be possible, so throw an error
    call crash('type_ugrid%find_nearest_face_center - couldnt find it!')

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

  subroutine add_nonchecked_neighbours_to_map_and_stack( self, map, stack, stackN, face_index)

    ! In/output variables:
    class(type_ugrid),                intent(in   ) :: self
    integer,  dimension(self%nFaces), intent(inout) :: map, stack
    integer,                          intent(inout) :: stackN
    integer,                          intent(in   ) :: face_index

    ! Local variables:
    integer :: ci, neighbour_face_index

    do ci = 1, self%face_face_connectivity_n( face_index)
      neighbour_face_index = self%face_face_connectivity( face_index, ci)
      if (map( neighbour_face_index) == 0) then
        map( neighbour_face_index) = 1
        stackN = stackN + 1
        stack( stackN) = neighbour_face_index
      end if
    end do

  end subroutine add_nonchecked_neighbours_to_map_and_stack

end submodule ugrid_search
