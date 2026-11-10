submodule( remapping_basic) remapping_calc_CSR_operators

use mpi_basic, only: par

contains

  subroutine calc_CSR_remapping_operator( self, ugrid_src, ugrid_dst, choice_remapping_method)

    ! In/output variables:
    class(type_map),             intent(inout) :: self
    type(type_ugrid),            intent(in   ) :: ugrid_src
    type(type_ugrid),            intent(in   ) :: ugrid_dst
    type(type_remapping_method), intent(in   ) :: choice_remapping_method

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_CSR_remapping_operator'

    ! Add routine to path
    call init_routine( routine_name)

    if     (choice_remapping_method == remapping_method%nearest_neighbour()) then
      call calc_CSR_remapping_operator_nearest_neighbour( self, ugrid_src, ugrid_dst)
    elseif (choice_remapping_method == remapping_method%linear()) then
      call calc_CSR_remapping_operator_linear( self, ugrid_src, ugrid_dst)
    elseif (choice_remapping_method == remapping_method%first_order_conservative()) then
      call calc_CSR_remapping_operator_first_order_conservative( self, ugrid_src, ugrid_dst)
    elseif (choice_remapping_method == remapping_method%second_order_conservative()) then
      call calc_CSR_remapping_operator_second_order_conservative( self, ugrid_src, ugrid_dst)
    else
      call crash("unsupported remapping method '" // choice_remapping_method%str() // "'")
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_CSR_remapping_operator

  subroutine calc_CSR_remapping_operator_nearest_neighbour( self, ugrid_src, ugrid_dst)

    ! In/output variables:
    class(type_map),  intent(inout) :: self
    type(type_ugrid), intent(in   ) :: ugrid_src
    type(type_ugrid), intent(in   ) :: ugrid_dst

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_CSR_remapping_operator_nearest_neighbour'
    integer                     :: face_index_src, face_index_dst
    real(dp), dimension(2)      :: p

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_remapping_CSR_matrix( ugrid_src, ugrid_dst, self%M)

    face_index_src = 1
    do face_index_dst = ugrid_dst%pai%i1, ugrid_dst%pai%i2
      p = ugrid_dst%face_coordinates( face_index_dst,:)
      call ugrid_src%find_nearest_face_center( p, face_index_src)
      call self%M%add_entry( face_index_dst, face_index_src, 1._dp)
    end do
    call self%M%finalise()

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_CSR_remapping_operator_nearest_neighbour

  subroutine calc_CSR_remapping_operator_linear( self, ugrid_src, ugrid_dst)

    ! In/output variables:
    class(type_map),  intent(inout) :: self
    type(type_ugrid), intent(in   ) :: ugrid_src
    type(type_ugrid), intent(in   ) :: ugrid_dst

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_CSR_remapping_operator_linear'

    ! Add routine to path
    call init_routine( routine_name)

    call crash('fixme')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_CSR_remapping_operator_linear

  subroutine calc_CSR_remapping_operator_first_order_conservative( self, ugrid_src, ugrid_dst)

    ! In/output variables:
    class(type_map),  intent(inout) :: self
    type(type_ugrid), intent(in   ) :: ugrid_src
    type(type_ugrid), intent(in   ) :: ugrid_dst

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_CSR_remapping_operator_first_order_conservative'
    type(type_CSR_matrix_dp)    :: M_xdy
    type(type_CSR_matrix_dp)    :: M_mxydx
    type(type_CSR_matrix_dp)    :: M_xydy
    integer                     :: row, k1, k2, k
    real(dp)                    :: A_overlap_tot

    ! Add routine to path
    call init_routine( routine_name)

    call calc_line_integral_matrices( ugrid_src, ugrid_dst, M_xdy, M_mxydx, M_xydy)

    self%M = M_xdy

    do row = self%ct%pai_dst%i1, self%ct%pai_dst%i2

      ! Calculate total area of overlap for this dst face
      A_overlap_tot = 0._dp
      k1 = M_xdy%ptr( row)
      k2 = M_xdy%ptr( row+1) - 1
      do k = k1, k2
        A_overlap_tot = A_overlap_tot + M_xdy%val( k)
      end do

      ! Scale LI_xdy to find remapping weights
      do k = k1, k2
        self%M%val( k) = self%M%val( k) / A_overlap_tot
      end do

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_CSR_remapping_operator_first_order_conservative

  subroutine calc_CSR_remapping_operator_second_order_conservative( self, ugrid_src, ugrid_dst)

    ! In/output variables:
    class(type_map),  intent(inout) :: self
    type(type_ugrid), intent(in   ) :: ugrid_src
    type(type_ugrid), intent(in   ) :: ugrid_dst

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_CSR_remapping_operator_second_order_conservative'

    ! Add routine to path
    call init_routine( routine_name)

    call crash('fixme')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_CSR_remapping_operator_second_order_conservative

  subroutine allocate_remapping_CSR_matrix( ugrid_src, ugrid_dst, M)

    ! In/output variables:
    type(type_ugrid),         intent(in   ) :: ugrid_src
    type(type_ugrid),         intent(in   ) :: ugrid_dst
    type(type_CSR_matrix_dp), intent(inout) :: M

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_remapping_CSR_matrix'

    ! Add routine to path
    call init_routine( routine_name)

    call M%allocate( &
      ugrid_dst%pai%n, &
      ugrid_src%pai%n, &
      ugrid_dst%pai%n_loc, &
      ugrid_src%pai%n_loc, &
      1, &
      ugrid_src%pai, &
      ugrid_dst%pai)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_remapping_CSR_matrix

  subroutine calc_line_integral_matrices( ugrid_src, ugrid_dst, M_xdy, M_mxydx, M_xydy)
    ! Calculate the three line integral matrices M_xdy, M_mxydx, M_xydy, where:
    !
    ! M_xdy(i,j) = [the line integral of [x dy] around the intersection of destinatioj
    !               ugrid face i and source ugrid face j]
    !
    ! ... and similarly for M_mxydx and M_xydy

    ! In/output variables:
    type(type_ugrid),         intent(in   ) :: ugrid_src
    type(type_ugrid),         intent(in   ) :: ugrid_dst
    type(type_CSR_matrix_dp), intent(inout) :: M_xdy
    type(type_CSR_matrix_dp), intent(inout) :: M_mxydx
    type(type_CSR_matrix_dp), intent(inout) :: M_xydy

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'calc_line_integral_matrices'
    integer                               :: face_index_src
    integer                               :: face_index_dst
    real(dp), dimension(:,:), allocatable :: poly_src
    real(dp), dimension(:,:), allocatable :: poly_dst
    integer,  dimension(ugrid_src%nFaces) :: map, stack
    integer                               :: stackN
    real(dp), dimension(2)                :: p
    integer                               :: ncycles
    real(dp), dimension(:,:), allocatable :: poly_intersect
    real(dp)                              :: LI_xdy, LI_mxydx, LI_xydy
    real(dp), parameter                   :: tol_dist = 1e-3_dp
    integer                               :: ci, neighbour_face_index_src

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_remapping_CSR_matrix( ugrid_src, ugrid_dst, M_xdy)
    call allocate_remapping_CSR_matrix( ugrid_src, ugrid_dst, M_mxydx)
    call allocate_remapping_CSR_matrix( ugrid_src, ugrid_dst, M_xydy)

    face_index_src = 1
    do face_index_dst = ugrid_dst%pai%i1, ugrid_dst%pai%i2

      ! The polygon describing this destination mesh Voronoi cell
      call ugrid_dst%get_face_as_polygon( face_index_dst, poly_dst)

      ! Initialise map and stack with the source ugrid face
      ! nearest to the destination ugrid face
      p = ugrid_dst%face_coordinates( face_index_dst,:)
      call ugrid_src%find_nearest_face_center( p, face_index_src)

      map = 0
      map( face_index_src) = 1
      stackN = 1
      stack( stackN) = face_index_src

      ncycles = 0
      do while (stackN > 0)

        ! Safety
        ncycles = ncycles + 1
        if (ncycles > ugrid_src%nFaces) call crash('floodfill-style outward search got stuck')

        ! Take the last source vertex from the stack
        face_index_src = stack( stackN)
        stackN = stackN - 1
        map( face_index_src) = 2

        ! Get this source grid face as a polygon
        call ugrid_src%get_face_as_polygon( face_index_src, poly_src)

        ! Calculate intersection
        call calc_convex_polygon_intersection( poly_src, poly_dst, poly_intersect)

        ! If the intersection is empty, move on
        if (.not. allocated( poly_intersect)) cycle

        ! Calculate the three line integrals around the perimeter of the intersection
        call calc_line_integrals_around_polygon( poly_intersect, LI_xdy, LI_mxydx, LI_xydy, tol_dist)

        ! Add line integrals to sparse matrices
        call M_xdy%add_entry  ( face_index_dst, face_index_src, LI_xdy)
        call M_mxydx%add_entry( face_index_dst, face_index_src, LI_mxydx)
        call M_xydy%add_entry ( face_index_dst, face_index_src, LI_xydy)

        ! Add non-checked neighbours to the stack to continue the floodfill-style outward expansion
        do ci = 1, ugrid_src%face_face_connectivity_n( face_index_src)
          neighbour_face_index_src = ugrid_src%face_face_connectivity( face_index_src, ci)
          if (map( neighbour_face_index_src) == 0) then
            map( neighbour_face_index_src) = 1
            stackN = stackN + 1
            stack( stackN) = neighbour_face_index_src
          end if
        end do

      end do

    end do

    call M_xdy%finalise()
    call M_mxydx%finalise()
    call M_xydy%finalise()

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_line_integral_matrices

end submodule remapping_calc_CSR_operators