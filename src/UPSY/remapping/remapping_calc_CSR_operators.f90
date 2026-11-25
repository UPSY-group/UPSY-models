submodule( remapping_basic) remapping_calc_CSR_operators

#include <petsc/finclude/petscsys.h>

contains

  subroutine calc_CSR_remapping_operator_nearest_neighbour( self)

    ! In/output variables:
    class(type_map), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_CSR_remapping_operator_nearest_neighbour'
    type(type_ugrid)            :: ugrid_src
    type(type_ugrid)            :: ugrid_dst
    integer                     :: face_index_src, face_index_dst
    real(dp), dimension(2)      :: p

    ! Add routine to path
    call init_routine( routine_name)

    ! Create ugrids from source and destination grids
    call ugrid_src%build_from_grid( self%ct%grid_src, self%ct%Arakawa_grid_src)
    call ugrid_dst%build_from_grid( self%ct%grid_dst, self%ct%Arakawa_grid_dst)

    call allocate_remapping_CSR_matrix( ugrid_src, ugrid_dst, self%M)

    face_index_src = 1
    do face_index_dst = ugrid_dst%pai%i1, ugrid_dst%pai%i2

      p = ugrid_dst%face_coordinates( face_index_dst,:)

      ! Assign NaN to destination ugrid faces lying outside of the source ugrid bounding box
      if (p(1) < ugrid_src%xmin + ugrid_src%tol_dist .or. &
          p(1) > ugrid_src%xmax - ugrid_src%tol_dist .or. &
          p(2) < ugrid_src%ymin + ugrid_src%tol_dist .or. &
          p(2) > ugrid_src%ymax - ugrid_src%tol_dist) then
        call self%M%add_empty_row( face_index_dst)
        cycle
      end if

      call ugrid_src%find_containing_face( p, face_index_src)
      call self%M%add_entry( face_index_dst, face_index_src, 1._dp)

    end do

    call self%M%finalise()

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_CSR_remapping_operator_nearest_neighbour

  subroutine calc_CSR_remapping_operator_first_order_conservative( self)

    ! In/output variables:
    class(type_map), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_CSR_remapping_operator_first_order_conservative'
    type(type_ugrid)            :: ugrid_src
    type(type_ugrid)            :: ugrid_dst
    type(type_CSR_matrix_dp)    :: M_xdy
    type(type_CSR_matrix_dp)    :: M_mxydx
    type(type_CSR_matrix_dp)    :: M_xydy
    integer                     :: face_index_dst, k1, k2, k
    real(dp)                    :: A_overlap_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Create ugrids from source and destination grids
    call ugrid_src%build_from_grid( self%ct%grid_src, self%ct%Arakawa_grid_src)
    call ugrid_dst%build_from_grid( self%ct%grid_dst, self%ct%Arakawa_grid_dst)

    call calc_line_integral_matrices( ugrid_src, ugrid_dst, M_xdy, M_mxydx, M_xydy)

    self%M = M_xdy

    do face_index_dst = ugrid_dst%pai%i1, ugrid_dst%pai%i2

      ! Calculate total area of overlap for this dst face
      A_overlap_tot = 0._dp
      k1 = M_xdy%ptr( face_index_dst)
      k2 = M_xdy%ptr( face_index_dst+1) - 1
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

  subroutine calc_CSR_remapping_operator_second_order_conservative( self)

    ! In/output variables:
    class(type_map),  intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_CSR_remapping_operator_second_order_conservative'
    type(type_ugrid)            :: ugrid_src
    type(type_ugrid)            :: ugrid_dst
    type(type_CSR_matrix_dp)    :: M0_CSR, M_ddx_CSR, M_ddy_CSR
    type(type_CSR_matrix_dp)    :: M_xdy_CSR, M_mxydx_CSR, M_xydy_CSR
    type(tMat)                  :: w0, w1x, w1y
    type(tMat)                  :: M
    integer                     :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    call calc_2nd_ord_cons_ugrids_and_map_gradient_matrices( self, ugrid_src, ugrid_dst, M0_CSR, M_ddx_CSR, M_ddy_CSR)
    call calc_line_integral_matrices( ugrid_src, ugrid_dst, M_xdy_CSR, M_mxydx_CSR, M_xydy_CSR)
    call calc_2nd_ord_cons_w_matrices( ugrid_src, ugrid_dst, M_xdy_CSR, M_mxydx_CSR, M_xydy_CSR, w0, w1x, w1y)
    call calc_2nd_ord_cons_remapping_matrix( M0_CSR, M_ddx_CSR, M_ddy_CSR, w0, w1x, w1y, M)
    call mat_petsc2CSR( M, self%M)

    ! Clean up after yourself
    PetscCall( MatDestroy( w0, ierr))
    PetscCall( MatDestroy( w1x, ierr))
    PetscCall( MatDestroy( w1y, ierr))
    PetscCall( MatDestroy( M, ierr))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_CSR_remapping_operator_second_order_conservative

  subroutine calc_2nd_ord_cons_ugrids_and_map_gradient_matrices( self, ugrid_src, ugrid_dst, M0_CSR, M_ddx_CSR, M_ddy_CSR)
    ! Calculate the mapping matrix M0 and the gradient matrices M_ddx, M_ddy

    ! In/output variables:
    class(type_map),          intent(inout) :: self
    type(type_ugrid)        , intent(  out) :: ugrid_src
    type(type_ugrid)        , intent(  out) :: ugrid_dst
    type(type_CSR_matrix_dp), intent(  out) :: M0_CSR
    type(type_CSR_matrix_dp), intent(  out) :: M_ddx_CSR
    type(type_CSR_matrix_dp), intent(  out) :: M_ddy_CSR

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_2nd_ord_cons_ugrids_and_map_gradient_matrices'

    ! Add routine to path
    call init_routine( routine_name)

    ! Destination grid is always the same
    call ugrid_dst%build_from_grid( self%ct%grid_dst, self%ct%Arakawa_grid_dst)

    select type (g => self%ct%grid_src)
      class default
        call crash('unsupported grid type')
      class is (type_mesh)

        if     (self%ct%Arakawa_grid_src == Arakawa_grid%a()) then
          call ugrid_src%build_from_grid( self%ct%grid_src, Arakawa_grid%b())
          M0_CSR    = g%M_map_a_b
          M_ddx_CSR = g%M_ddx_a_b
          M_ddy_CSR = g%M_ddy_a_b
        elseif (self%ct%Arakawa_grid_src == Arakawa_grid%b()) then
          call ugrid_src%build_from_grid( self%ct%grid_src, Arakawa_grid%b())
          call M0_CSR%identity( g%nTri, g%nTri_loc, g%pai_Tri, g%pai_Tri)
          M_ddx_CSR = g%M_ddx_b_b
          M_ddy_CSR = g%M_ddy_b_b
        elseif (self%ct%Arakawa_grid_src == Arakawa_grid%c()) then
          call crash('remapping from mesh edges not (yet) supported')
        else
          call crash('unsupported Arakawa grid')
        end if

      class is (type_xy_grid)

        call ugrid_src%build_from_grid( self%ct%grid_src, self%ct%Arakawa_grid_src)
        call M0_CSR%identity( g%n, g%n_loc, g%pai, g%pai)
        call calc_xy_grid_matrix_operators( g, M_ddx_CSR, M_ddy_CSR)

      class is (type_lonlat_grid)
        call crash('fixme')

    end select

    call checksum( M0_CSR   , 'M0_CSR'   )
    call checksum( M_ddx_CSR, 'M_ddx_CSR')
    call checksum( M_ddy_CSR, 'M_ddy_CSR')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_2nd_ord_cons_ugrids_and_map_gradient_matrices

  subroutine calc_2nd_ord_cons_w_matrices( ugrid_src, ugrid_dst, M_xdy_CSR, M_mxydx_CSR, M_xydy_CSR, w0, w1x, w1y)

    ! In/output variables
    type(type_ugrid)        , intent(in   ) :: ugrid_src
    type(type_ugrid)        , intent(in   ) :: ugrid_dst
    type(type_CSR_matrix_dp), intent(in   ) :: M_xdy_CSR
    type(type_CSR_matrix_dp), intent(in   ) :: M_mxydx_CSR
    type(type_CSR_matrix_dp), intent(in   ) :: M_xydy_CSR
    type(tMat),               intent(  out) :: w0, w1x, w1y

    ! Local variables:
    character(len=*), parameter                 :: routine_name = 'calc_2nd_ord_cons_w_matrices'
    type(tMat)                                  :: M_xdy, M_mxydx, M_xydy
    integer                                     :: ierr
    integer                                     :: nnz_per_row_max
    integer                                     :: istart, iend, n, k, face_index
    integer                                     :: ncols
    integer,  dimension(:), allocatable, target :: cols
    real(dp), dimension(:), allocatable, target :: vals, w0_row, w1x_row, w1y_row
    integer,  dimension(:), pointer             :: cols_
    real(dp), dimension(:), pointer             :: vals_, w0_row_, w1x_row_, w1y_row_
    real(dp)                                    :: A_overlap_tot, x_src, y_src

    ! Add routine to path
    call init_routine( routine_name)

    ! Convert matrices from UPSY-CSR to PETSc types
    call mat_CSR2petsc( M_xdy_CSR  , M_xdy  )
    call mat_CSR2petsc( M_mxydx_CSR, M_mxydx)
    call mat_CSR2petsc( M_xydy_CSR , M_xydy )

    ! Calculate w0, w1x, w1y for the mesh-to-grid remapping operator
    PetscCall( MatConvert( M_xdy, MATAIJ, MAT_INITIAL_MATRIX, w0, ierr))
    PetscCall( MatConvert( M_xdy, MATAIJ, MAT_INITIAL_MATRIX, w1x, ierr))
    PetscCall( MatConvert( M_xdy, MATAIJ, MAT_INITIAL_MATRIX, w1y, ierr))

    ! Find maximum number of non-zeros per row
    nnz_per_row_max = M_xdy_CSR%nnz_per_row_max()

    ! allocate memory for a single matrix row
    allocate( cols(    nnz_per_row_max))
    allocate( vals(    nnz_per_row_max))
    allocate( w0_row(  nnz_per_row_max))
    allocate( w1x_row( nnz_per_row_max))
    allocate( w1y_row( nnz_per_row_max))

    cols_    => cols
    vals_    => vals
    w0_row_  => w0_row
    w1x_row_ => w1x_row
    w1y_row_ => w1y_row

    call MatGetOwnershipRange( M_xdy, istart, iend, ierr)

    do n = istart+1, iend ! +1 because PETSc indexes from 0

      ! Calculate area of overlap
      PetscCall( MatGetRow( M_xdy, n-1, ncols, cols_, vals_, ierr))
      A_overlap_tot = sum( vals_( 1:ncols))
      PetscCall( MatRestoreRow( M_xdy, n-1, ncols, cols_, vals_, ierr))

      ! w0
      PetscCall( MatGetRow( M_xdy, n-1, ncols, cols_, vals_, ierr))
      do k = 1, ncols
        w0_row_( k) = vals_( k) / A_overlap_tot
        PetscCall( MatSetValues( w0, 1, [n-1], 1, [cols_( k)], [w0_row_( k)], INSERT_VALUES, ierr))
      end do
      PetscCall( MatRestoreRow( M_xdy, n-1, ncols, cols_, vals_, ierr))

      ! w1x
      PetscCall( MatGetRow( M_mxydx, n-1, ncols, cols_, vals_, ierr))
      do k = 1, ncols
        face_index = cols_( k)+1
        x_src = ugrid_src%face_coordinates( face_index,1)
        w1x_row_( k) = (vals_( k) / A_overlap_tot) - (x_src * w0_row( k))
        PetscCall( MatSetValues( w1x, 1, [n-1], 1, [cols_( k)], [w1x_row_( k)], INSERT_VALUES, ierr))
      end do
      PetscCall( MatRestoreRow( M_mxydx, n-1, ncols, cols_, vals_, ierr))

      ! w1y
      PetscCall( MatGetRow( M_xydy, n-1, ncols, cols_, vals_, ierr))
      do k = 1, ncols
        face_index = cols_( k)+1
        y_src = ugrid_src%face_coordinates( face_index,2)
        w1y_row_( k) = (vals_( k) / A_overlap_tot) - (y_src * w0_row( k))
        PetscCall( MatSetValues( w1y, 1, [n-1], 1, [cols_( k)], [w1y_row_( k)], INSERT_VALUES, ierr))
      end do
      PetscCall( MatRestoreRow( M_xydy, n-1, ncols, cols_, vals_, ierr))

    end do

    PetscCall( MatAssemblyBegin( w0 , MAT_FINAL_ASSEMBLY, ierr))
    PetscCall( MatAssemblyEnd(   w0 , MAT_FINAL_ASSEMBLY, ierr))
    PetscCall( MatAssemblyBegin( w1x, MAT_FINAL_ASSEMBLY, ierr))
    PetscCall( MatAssemblyEnd(   w1x, MAT_FINAL_ASSEMBLY, ierr))
    PetscCall( MatAssemblyBegin( w1y, MAT_FINAL_ASSEMBLY, ierr))
    PetscCall( MatAssemblyEnd(   w1y, MAT_FINAL_ASSEMBLY, ierr))

    ! Clean up after yourself
    PetscCall( MatDestroy( M_xdy  , ierr))
    PetscCall( MatDestroy( M_mxydx, ierr))
    PetscCall( MatDestroy( M_xydy , ierr))

    call checksum( w0 , 'w0' )
    call checksum( w1x, 'w1x')
    call checksum( w1y, 'w1y')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_2nd_ord_cons_w_matrices

  subroutine calc_2nd_ord_cons_remapping_matrix( M0_CSR, M_ddx_CSR, M_ddy_CSR, w0, w1x, w1y, M)
    !< Calculate the remapping matrix M

    ! In/output variables
    type(type_CSR_matrix_dp), intent(in   ) :: M0_CSR, M_ddx_CSR, M_ddy_CSR
    type(tMat),               intent(in   ) :: w0, w1x, w1y
    type(tMat),               intent(  out) :: M

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_2nd_ord_cons_remapping_matrix'
    type(tMat)                  :: M0, M_ddx, M_ddy
    integer                     :: ierr
    type(tMat)                  :: M_cons_1st_order
    type(tMat)                  :: M1, M2

    ! Add routine to path
    call init_routine( routine_name)

    ! Convert matrices from UPSY-CSR to PETSc types
    call mat_CSR2petsc( M0_CSR   , M0   )
    call mat_CSR2petsc( M_ddx_CSR, M_ddx)
    call mat_CSR2petsc( M_ddy_CSR, M_ddy)

    ! 1st-order = w0 * M0
    PetscCall( MatMatMult( w0 , M0, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M_cons_1st_order, ierr))

    ! 2nd-order = 1st-order + w1x * M_ddx + w1y * M_ddy
    PetscCall( MatMatMult( w1x, M_ddx, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M1, ierr))  ! This can be done more efficiently now that the non-zero structure is known...
    PetscCall( MatMatMult( w1y, M_ddy, MAT_INITIAL_MATRIX, PETSC_DEFAULT_REAL, M2, ierr))

    PetscCall( MatConvert( M_cons_1st_order, MATAIJ, MAT_INITIAL_MATRIX, M, ierr))
    PetscCall( MatAXPY( M, 1._dp, M1, DIFFERENT_NONZERO_PATTERN, ierr))
    PetscCall( MatAXPY( M, 1._dp, M2, DIFFERENT_NONZERO_PATTERN, ierr))

    ! Clean up after yourself
    PetscCall( MatDestroy( M0, ierr))
    PetscCall( MatDestroy( M_ddx, ierr))
    PetscCall( MatDestroy( M_ddy, ierr))
    PetscCall( MatDestroy( M_cons_1st_order, ierr))
    PetscCall( MatDestroy( M1, ierr))
    PetscCall( MatDestroy( M2, ierr))

    call checksum( M, 'M')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_2nd_ord_cons_remapping_matrix

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
    ! M_xdy(i,j) = [the line integral of [x dy] around the intersection of destination
    !               ugrid face i and source ugrid face j]
    !
    ! ... and similarly for M_mxydx and M_xydy

    ! In/output variables:
    type(type_ugrid),         intent(inout) :: ugrid_src
    type(type_ugrid),         intent(inout) :: ugrid_dst
    type(type_CSR_matrix_dp), intent(inout) :: M_xdy
    type(type_CSR_matrix_dp), intent(inout) :: M_mxydx
    type(type_CSR_matrix_dp), intent(inout) :: M_xydy

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'calc_line_integral_matrices'
    integer                               :: face_index_src
    integer                               :: face_index_dst
    real(dp), dimension(:,:), allocatable :: poly_src
    real(dp), dimension(:,:), allocatable :: poly_dst
    type(type_floodfill)                  :: floodfill_src
    real(dp), dimension(2)                :: p
    integer                               :: ncycles
    real(dp), dimension(:,:), allocatable :: poly_intersect
    real(dp)                              :: LI_xdy, LI_mxydx, LI_xydy
    logical                               :: found_any_intersection

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_remapping_CSR_matrix( ugrid_src, ugrid_dst, M_xdy)
    call allocate_remapping_CSR_matrix( ugrid_src, ugrid_dst, M_mxydx)
    call allocate_remapping_CSR_matrix( ugrid_src, ugrid_dst, M_xydy)

    call floodfill_src%create( &
      ugrid_src%nFaces, ugrid_src%face_face_connectivity_n, ugrid_src%face_face_connectivity)

    face_index_src = 1
    do face_index_dst = ugrid_dst%pai%i1, ugrid_dst%pai%i2

      p = ugrid_dst%face_coordinates( face_index_dst,:)

      ! Assign NaN to destination ugrid faces lying outside of the source ugrid bounding box
      if (p(1) < ugrid_src%xmin + ugrid_src%tol_dist .or. &
          p(1) > ugrid_src%xmax - ugrid_src%tol_dist .or. &
          p(2) < ugrid_src%ymin + ugrid_src%tol_dist .or. &
          p(2) > ugrid_src%ymax - ugrid_src%tol_dist) then
        call M_xdy%add_empty_row  ( face_index_dst)
        call M_mxydx%add_empty_row( face_index_dst)
        call M_xydy%add_empty_row ( face_index_dst)
        cycle
      end if

      ! The polygon describing this destination mesh Voronoi cell
      call ugrid_dst%get_face_as_polygon( face_index_dst, poly_dst)

      ! Initialise map and stack with the source face
      ! containing the centroid of the destination face
      call ugrid_src%find_containing_face( p, face_index_src)

      call floodfill_src%initialise_floodfill_at_single_node( face_index_src)

      ncycles = 0
      found_any_intersection = .false.
      do while (floodfill_src%stack_front_ii2 >= floodfill_src%stack_front_ii1)

        ! Safety
        ncycles = ncycles + 1
        if (ncycles > ugrid_src%nFaces) call crash('floodfill-style outward search got stuck')

        ! Take the first source vertex from the stack
        call floodfill_src%take_first_element_from_front_stack( face_index_src)

        ! Get this source grid face as a polygon
        call ugrid_src%get_face_as_polygon( face_index_src, poly_src)

        ! Calculate intersection
        call calc_convex_polygon_intersection( poly_src, poly_dst, poly_intersect)

        if (.not. allocated( poly_intersect)) then
          if (ncycles == 1) then
            ! If the first source ugrid face we try doesn't intersect
            ! the destination ugrid face, we must still try its neighbours
            call floodfill_src%add_unchecked_neighbours_to_front_stack( face_index_src)
          end if
          cycle
        else
          found_any_intersection = .true.
        end if

        ! Calculate the three line integrals around the perimeter of the intersection
        call calc_line_integrals_around_polygon( poly_intersect, &
          LI_xdy, LI_mxydx, LI_xydy)

        ! Add line integrals to sparse matrices
        call M_xdy%add_entry  ( face_index_dst, face_index_src, LI_xdy)
        call M_mxydx%add_entry( face_index_dst, face_index_src, LI_mxydx)
        call M_xydy%add_entry ( face_index_dst, face_index_src, LI_xydy)

        ! Add non-checked neighbours to the stack to continue the floodfill-style outward expansion
        call floodfill_src%add_unchecked_neighbours_to_front_stack( face_index_src)

      end do

      ! Safety
      if (.not. found_any_intersection) call crash('could not find any intersecting source ugrid faces ' // &
        'for destination ugrid face {int_01}', int_01 = face_index_dst)

    end do

    call M_xdy%finalise()
    call M_mxydx%finalise()
    call M_xydy%finalise()

    call checksum( M_xdy  , 'M_xdy'  )
    call checksum( M_mxydx, 'M_mxydx')
    call checksum( M_xydy , 'M_xydy' )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_line_integral_matrices

end submodule remapping_calc_CSR_operators