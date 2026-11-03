module xy_grid_basic

  ! Functions for working with simple square x/y-grids

  use precisions, only: dp
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: crash, init_routine, finalise_routine
  use basic_grid_mod, only: atype_grid
  use parameters
  use reallocate_mod, only: reallocate
  use interpolation, only: linint_points
  use projections, only: inverse_oblique_sg_projection
  use mpi_distributed_memory, only: partition_list
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use mpi_f08, only: MPI_ALLREDUCE, MPI_INTEGER, MPI_MIN, MPI_MAX
  use parallel_array_info_type, only: type_par_arr_info

  implicit none

  private

  public :: type_xy_grid, setup_xy_grid, deallocate_xy_grid, calc_secondary_xy_grid_data, &
    calc_xy_grid_matrix_operators

  type, extends(atype_grid) :: type_xy_grid
    ! A regular x/y-grid

    ! Basic properties
    integer                                 :: nx           !           Number of grid cells in the x-direction
    integer                                 :: ny           !           Number of grid cells in the x-direction
    integer                                 :: n            !           Total number of grid cells (= nx * ny)
    real(dp)                                :: dx           ! [m]       Resolution (square grid, so dy = dx)
    real(dp), dimension(:    ), allocatable :: x            ! [m]       x-coordinates
    real(dp), dimension(:    ), allocatable :: y            ! [m]       y-coordinates
    real(dp)                                :: xmin         ! [m]       x and y range of the square covered by the grid
    real(dp)                                :: xmax         ! [m]
    real(dp)                                :: ymin         ! [m]
    real(dp)                                :: ymax         ! [m]

    ! Remapping data
    real(dp)                                :: tol_dist     ! [m]       Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    integer,  dimension(:,:  ), allocatable :: ij2n, n2ij   !           Conversion table for grid-form vs. vector-form data

    ! Lon/lat-coordinates
    real(dp)                                :: lambda_M
    real(dp)                                :: phi_M
    real(dp)                                :: beta_stereo
    real(dp), dimension(:,:  ), allocatable :: lon, lat

    ! Parallelisation
    integer                                 :: n1,n2,n_loc  ! Matrix rows owned by each process
    type(type_par_arr_info)                 :: pai

  contains

    procedure, public :: ij2fi, iijj2ni

  end type type_xy_grid

contains

  subroutine setup_xy_grid( name, xmin, xmax, ymin, ymax, dx, grid, lambda_M, phi_M, beta_stereo)
    !< Set up a square grid that covers the specified domain

    ! In/output variables:
    character(len=*),   intent(in   ) :: name
    real(dp),           intent(in   ) :: xmin, xmax, ymin, ymax        ! [m] Domain
    real(dp),           intent(in   ) :: dx                            ! [m] Resolution
    type(type_xy_grid), intent(  out) :: grid
    real(dp), optional, intent(in   ) :: lambda_M, phi_M, beta_stereo

    ! Local variables:
    character(len=*), parameter :: routine_name = 'setup_xy_grid'
    real(dp)                    :: xmid, ymid
    integer                     :: nx_pos, ny_pos
    integer                     :: i,j,ii,jj

    ! Add routine to path
    call init_routine( routine_name)

    grid%name = name
    grid%dx = dx

    ! Determine how many grid points are needed to fully cover the domain

    xmid = (xmin + xmax) / 2._dp
    nx_pos = 0
    do while (xmid + real( nx_pos,dp) * grid%dx + grid%dx / 2._dp < xmax)
      nx_pos = nx_pos + 1
    end do
    grid%nx = 1 + 2 * nx_pos

    ymid = (ymin + ymax) / 2._dp
    ny_pos = 0
    do while (ymid + real( ny_pos,dp) * grid%dx + grid%dx / 2._dp < ymax)
      ny_pos = ny_pos + 1
    end do
    grid%ny = 1 + 2 * ny_pos

    ! Fill in x and y
    allocate( grid%x( grid%nx))
    allocate( grid%y( grid%ny))

    do i = 1, grid%nx
      ii = i - (nx_pos+1)
      grid%x( i) = xmid + real( ii,dp) * grid%dx
    end do

    do j = 1, grid%ny
      jj = j - (ny_pos+1)
      grid%y( j) = ymid + real( jj,dp) * grid%dx
    end do

    ! Calculate secondary grid geometry data
    call calc_secondary_xy_grid_data( grid, lambda_M, phi_M, beta_stereo)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_xy_grid

  subroutine deallocate_xy_grid( grid)
    !< deallocate memory for a grid object

    ! In/output variables:
    type(type_xy_grid), intent(inout) :: grid

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_xy_grid'

    ! Add routine to path
    call init_routine( routine_name)

    if (allocated( grid%x   )) deallocate( grid%x   )
    if (allocated( grid%y   )) deallocate( grid%y   )
    if (allocated( grid%lon )) deallocate( grid%lon )
    if (allocated( grid%lat )) deallocate( grid%lat )
    if (allocated( grid%ij2n)) deallocate( grid%ij2n)
    if (allocated( grid%n2ij)) deallocate( grid%n2ij)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_xy_grid

  subroutine check_if_grids_are_identical( grid1, grid2, isso)
    !< Check if two grids are identical

    ! In/output variables:
    type(type_xy_grid), intent(in   ) :: grid1, grid2
    logical,         intent(  out) :: isso

    ! Local variables:
    real(dp), parameter :: rtol = 1E-9_dp
    integer             :: i,j

    isso = .true.

    ! Size
    if (grid1%nx /= grid2%nx .or. grid1%ny /= grid2%ny) then
      isso = .false.
      return
    end if

    ! Coordinates
    do i = 1, grid1%nx
      if (abs(grid1%x( i) - grid2%x( i)) > abs(rtol * max(grid1%x(i), grid2%x(i)))) then
        isso = .false.
        return
      end if
    end do
    do j = 1, grid1%ny
      if (abs(grid1%y( j) - grid2%y( j)) > abs(rtol * max(grid1%y(j), grid2%y(j)))) then
        isso = .false.
        return
      end if
    end do

  end subroutine check_if_grids_are_identical

  subroutine calc_secondary_xy_grid_data( grid, lambda_M, phi_M, beta_stereo)
    !< Calculate secondary geometry data for a square grid

    ! In/output variables:
    type(type_xy_grid), intent(inout) :: grid
    real(dp), optional, intent(in   ) :: lambda_M, phi_M, beta_stereo

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_secondary_xy_grid_data'
    integer                     :: i,j
    real(dp), parameter         :: tol = 1E-9_dp

    ! Add routine to path
    call init_routine( routine_name)

    ! Resolution
    grid%dx = abs( grid%x( 2) - grid%x( 1))

    ! Safety
    do i = 1, grid%nx-1
      if (1._dp - abs(grid%x( i+1) - grid%x( i)) / grid%dx > 1E-6_dp) then
        call crash( trim( grid%name) // '" has an irregular x-dimension!')
      end if
    end do
    do j = 1, grid%ny-1
      if (1._dp - abs(grid%y( j+1) - grid%y( j)) / grid%dx > 1E-6_dp) then
        call crash( trim( grid%name) // '" has an irregular y-dimension!')
      end if
    end do

    ! Domain size
    grid%xmin = minval( grid%x)
    grid%xmax = maxval( grid%x)
    grid%ymin = minval( grid%y)
    grid%ymax = maxval( grid%y)

    ! Tolerance; points lying within this distance of each other are treated as identical
    grid%tol_dist = ((grid%xmax - grid%xmin) + (grid%ymax - grid%ymin)) * tol / 2._dp

    ! Conversion tables for grid-form vs. vector-form data
    call calc_field_to_vector_form_translation_tables( grid)

    ! Parallelisation
    call setup_grid_parallelisation( grid)

    ! Lon/lat-coordinates
    if (present( lambda_M) .or. present( phi_M) .or. present( beta_stereo)) then

      ! Safety
      if (.not. (present( lambda_M) .and. present( phi_M) .and. present( beta_stereo))) then
        call crash('need lambda_M, phi_M, and beta_stereo!')
      end if

      ! allocate memory
      allocate( grid%lon( grid%nx, grid%ny))
      allocate( grid%lat( grid%nx, grid%ny))

      ! Calculate lon/lat-coordinates for each grid point using the provided oblique stereographic projection parameters
      do i = 1, grid%nx
      do j = 1, grid%ny
        call inverse_oblique_sg_projection( grid%x( i), grid%y( j), lambda_M, phi_M, beta_stereo, grid%lon( i,j), grid%lat( i,j))
      end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_secondary_xy_grid_data

  subroutine calc_xy_grid_matrix_operators( grid, M_ddx_CSR, M_ddy_CSR)
    !< Calculate matrix operators for partial derivatives on a regular grid
    !< (needed for conservative remapping)

    ! In/output variables:
    type(type_xy_grid),              intent(in   ) :: grid
    type(type_CSR_matrix_dp), intent(  out) :: M_ddx_CSR, M_ddy_CSR

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_xy_grid_matrix_operators'
    integer                     :: ncols, nrows, ncols_loc, nrows_loc, nnz_per_row_est, nnz_est_proc
    integer                     :: row, i, j, col
    real(dp)                    :: valpos, valneg

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the matrices using the native UFEMISM CSR-matrix format
    ! =====================================================================

    ! Matrix size
    ncols           = grid%n      ! from
    ncols_loc       = grid%n_loc
    nrows           = grid%n      ! to
    nrows_loc       = grid%n_loc
    nnz_per_row_est = 2
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call M_ddx_CSR%allocate( nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)
    call M_ddy_CSR%allocate( nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! == Fill matrix coefficients
    ! ===========================

    valpos =  1._dp / (2._dp * grid%dx)
    valneg = -1._dp / (2._dp * grid%dx)

    do row = grid%n1, grid%n2

      ! Grid indices
      i = grid%n2ij( row,1)
      j = grid%n2ij( row,2)

      ! Skip the border
      if (i == 1 .or. i == grid%nx .or. j == 1 .or. j == grid%ny) then
        M_ddx_CSR%ptr( row+1) = M_ddx_CSR%ptr( row)
        M_ddy_CSR%ptr( row+1) = M_ddy_CSR%ptr( row)
        cycle
      end if

    ! == d/dx

      ! Left
      col = grid%ij2n( i-1,j)
      M_ddx_CSR%nnz = M_ddx_CSR%nnz + 1
      M_ddx_CSR%ind( M_ddx_CSR%nnz) = col
      M_ddx_CSR%val( M_ddx_CSR%nnz) = valneg

      ! Left
      col = grid%ij2n( i+1,j)
      M_ddx_CSR%nnz = M_ddx_CSR%nnz + 1
      M_ddx_CSR%ind( M_ddx_CSR%nnz) = col
      M_ddx_CSR%val( M_ddx_CSR%nnz) = valpos

      ! Ptr
      M_ddx_CSR%ptr( row+1) = M_ddx_CSR%ptr( row) + 2

    ! == d/dy

      ! Left
      col = grid%ij2n( i,j-1)
      M_ddy_CSR%nnz = M_ddy_CSR%nnz + 1
      M_ddy_CSR%ind( M_ddy_CSR%nnz) = col
      M_ddy_CSR%val( M_ddy_CSR%nnz) = valneg

      ! Left
      col = grid%ij2n( i,j+1)
      M_ddy_CSR%nnz = M_ddy_CSR%nnz + 1
      M_ddy_CSR%ind( M_ddy_CSR%nnz) = col
      M_ddy_CSR%val( M_ddy_CSR%nnz) = valpos

      ! Ptr
      M_ddy_CSR%ptr( row+1) = M_ddy_CSR%ptr( row) + 2

    end do ! DO row = grid%n1, grid%n2

    call M_ddx_CSR%finalise
    call M_ddy_CSR%finalise

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_xy_grid_matrix_operators

  subroutine calc_field_to_vector_form_translation_tables( grid)
    !< Calculate grid-cell-to-matrix-row translation tables

    ! In/output variables
    type(type_xy_grid), intent(inout) :: grid

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_field_to_vector_form_translation_tables'
    integer                     :: i,j,fi

    ! Add routine to path
    call init_routine( routine_name)

    ! Total number of grid cells
    grid%n = grid%nx * grid%ny

    ! Allocate memory
    if (allocated( grid%ij2n)) deallocate( grid%ij2n)
    if (allocated( grid%n2ij)) deallocate( grid%n2ij)
    allocate( grid%ij2n( grid%nx, grid%ny), source = 0)
    allocate( grid%n2ij( grid%n , 2      ), source = 0)

    ! Fill in tables
    do i = 1, grid%nx
    do j = 1, grid%ny
      fi = grid%ij2fi( i,j)
      grid%ij2n( i,j) = fi
      grid%n2ij( fi,:) = [i,j]
    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_field_to_vector_form_translation_tables

  subroutine setup_grid_parallelisation( grid)

    ! In/output variables
    type(type_xy_grid), intent(inout) :: grid

    ! Local variables:
    character(len=*), parameter :: routine_name = 'setup_grid_parallelisation'
    integer                     :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Parallelisation domains
    call partition_list( grid%n, par%i, par%n, grid%n1, grid%n2)
    grid%n_loc = grid%n2 + 1 - grid%n1

    ! Parallel array info
    grid%pai%n  = grid%n

    grid%pai%i1    = grid%n1
    grid%pai%i2    = grid%n2
    grid%pai%n_loc = grid%n_loc

    call MPI_ALLREDUCE( grid%pai%i1, grid%pai%i1_node, 1, MPI_INTEGER, MPI_MIN, par%mpi_comm_node, ierr)
    call MPI_ALLREDUCE( grid%pai%i2, grid%pai%i2_node, 1, MPI_INTEGER, MPI_MAX, par%mpi_comm_node, ierr)
    grid%pai%n_node = grid%pai%i2_node + 1 - grid%pai%i1_node

    grid%pai%i1_nih = grid%pai%i1_node
    grid%pai%i2_nih = grid%pai%i2_node
    grid%pai%n_nih  = grid%pai%n_node

    grid%pai%i1_hle = 0
    grid%pai%i2_hle = -1
    grid%pai%n_hle  = 0

    grid%pai%i1_hli = 0
    grid%pai%i2_hli = -1
    grid%pai%n_hli  = 0

    grid%pai%i1_hre = 0
    grid%pai%i2_hre = -1
    grid%pai%n_hre  = 0

    grid%pai%i1_hri = 0
    grid%pai%i2_hri = -1
    grid%pai%n_hri  = 0

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_grid_parallelisation

  pure function iijj2ni( self,ii,jj) result( ni)
    class(type_xy_grid), intent(in) :: self
    integer,             intent(in) :: ii,jj
    integer                         :: ni
    ni = (ii-1) * (self%ny+1) + jj
  end function iijj2ni

  pure function ij2fi( self,i,j) result( fi)
    class(type_xy_grid), intent(in) :: self
    integer,             intent(in) :: i,j
    integer                         :: fi
    fi = (i-1) * (self%ny) + j
  end function ij2fi

end module xy_grid_basic
