module smooth_gridded_data

  ! Functions for smoothing data on a square grid

  use precisions, only: dp
  use grid_types, only: type_grid
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary, distribute_gridded_data_from_primary
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash

  implicit none

  private

  public :: smooth_Gaussian_grid, extrapolate_fillvalue_Gaussian_grid

  interface smooth_Gaussian_grid
    procedure smooth_Gaussian_grid_2D
    procedure smooth_Gaussian_grid_3D
  end interface smooth_Gaussian_grid

  interface extrapolate_fillvalue_Gaussian_grid
    procedure extrapolate_fillvalue_Gaussian_grid_2D
    procedure extrapolate_fillvalue_Gaussian_grid_3D
  end interface extrapolate_fillvalue_Gaussian_grid

contains

subroutine smooth_Gaussian_grid_2D( grid, d_grid_vec_partial, r)
  !< Apply a Gaussian smoothing filter with sigma = r to the 2D data field d

  ! In/output variables:
  type(type_grid),        intent(in   ) :: grid
  real(dp), dimension(:), intent(inout) :: d_grid_vec_partial
  real(dp),               intent(in   ) :: r                    !< [m] Smoothing radius

  ! Local variables:
  character(len=256), parameter         :: routine_name = 'smooth_Gaussian_grid_2D'
  integer                               :: n
  real(dp), dimension(:,:), allocatable :: d_grid_tot
  real(dp), dimension(:,:), allocatable :: d_grid_tot_smoothed
  real(dp), dimension(:),   allocatable :: f
  integer                               :: i,j,k,ii,jj

  ! Add routine to path
  call init_routine( routine_name)

  ! Number of cells to extend the data by (3 standard deviations is enough to capture the tails of the normal distribution)
  n = ceiling( r / grid%dx) * 3

  ! Calculate the 1-D smoothing filter
  allocate( f( -n:n))
  f = 0._dp
  do k = -n, n
    f( k) = exp( -0.5_dp * (real( k,dp) * grid%dx/r)**2)
  end do
  f = f / sum(f)

  ! allocate memory
  allocate( d_grid_tot(          grid%nx,grid%ny), source = 0._dp)
  allocate( d_grid_tot_smoothed( grid%nx,grid%ny), source = 0._dp)

  ! Gather data to the primary in grid form
  call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid_tot)

  ! Let the primary do the actual work
  if (par%primary) then

    ! First smooth in the x-direction
    d_grid_tot_smoothed = 0._dp
    do i = 1, grid%nx
    do j = 1, grid%ny
      DO k = -n, n
        ii = max( 1, min( grid%nx, i+k ))
        d_grid_tot_smoothed( i,j) = d_grid_tot_smoothed( i,j) + d_grid_tot( ii,j) * f( k)
      end do
    end do
    end do
    d_grid_tot = d_grid_tot_smoothed

    ! then smooth in the y-direction
    d_grid_tot_smoothed = 0._dp
    do i = 1, grid%nx
    do j = 1, grid%ny
      DO k = -n, n
        jj = max( 1, min( grid%ny, j+k ))
        d_grid_tot_smoothed( i,j) = d_grid_tot_smoothed( i,j) + d_grid_tot( i,jj) * f( k)
      end do
    end do
    end do
    d_grid_tot = d_grid_tot_smoothed

  end if

  ! Add routine to path
  call finalise_routine( routine_name)

  ! Distributed smoothed data back from the primary
  call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid_tot)

end subroutine smooth_Gaussian_grid_2D

subroutine smooth_Gaussian_grid_3D( grid, d_grid_vec_partial, r)
  !< Apply a Gaussian smoothing filter with sigma = r to the 3D data field d

  ! In/output variables:
  type(type_grid),          intent(in   ) :: grid
  real(dp), dimension(:,:), intent(inout) :: d_grid_vec_partial
  real(dp),                 intent(in   ) :: r                    !< [m] Smoothing radius

  ! Local variables:
  character(len=256), parameter           :: routine_name = 'smooth_Gaussian_grid_3D'
  integer                                 :: n
  real(dp), dimension(:,:,:), allocatable :: d_grid_tot
  real(dp), dimension(:,:,:), allocatable :: d_grid_tot_smoothed
  real(dp), dimension(:    ), allocatable :: f
  integer                                 :: i,j,k,ii,jj

  ! Add routine to path
  call init_routine( routine_name)

  ! Number of cells to extend the data by (3 standard deviations is enough to capture the tails of the normal distribution)
  n = ceiling( r / grid%dx) * 3

  ! Calculate the 1-D smoothing filter
  allocate( f( -n:n))
  f = 0._dp
  do k = -n, n
    f( k) = exp( -0.5_dp * (real( k,dp) * grid%dx/r)**2)
  end do
  f = f / sum(f)

  ! allocate memory
  allocate( d_grid_tot(          grid%nx,grid%ny,size( d_grid_vec_partial,2)), source = 0._dp)
  allocate( d_grid_tot_smoothed( grid%nx,grid%ny,size( d_grid_vec_partial,2)), source = 0._dp)

  ! Gather data to the primary in grid form
  call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid_tot)

  ! Let the primary do the actual work
  if (par%primary) then

    ! First smooth in the x-direction
    d_grid_tot_smoothed = 0._dp
    do i = 1, grid%nx
    do j = 1, grid%ny
      do k = -n, n
        ii = max( 1, min( grid%nx, i+k ))
        d_grid_tot_smoothed( i,j,:) = d_grid_tot_smoothed( i,j,:) + d_grid_tot( ii,j,:) * f( k)
      end do
    end do
    end do
    d_grid_tot = d_grid_tot_smoothed

    ! then smooth in the y-direction
    d_grid_tot_smoothed = 0._dp
    do i = 1, grid%nx
    do j = 1, grid%ny
      do k = -n, n
        jj = max( 1, min( grid%ny, j+k ))
        d_grid_tot_smoothed( i,j,:) = d_grid_tot_smoothed( i,j,:) + d_grid_tot( i,jj,:) * f( k)
      end do
    end do
    end do
    d_grid_tot = d_grid_tot_smoothed

  end if

  ! Add routine to path
  call finalise_routine( routine_name)

  ! Distributed smoothed data back from the primary
  call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid_tot)

end subroutine smooth_Gaussian_grid_3D

subroutine extrapolate_fillvalue_Gaussian_grid_2D( grid, d_grid_vec_partial, sigma, fill_value)
  ! Replace fillvalues by Gaussian extrapolation (with radius sigma) of the rest of the data

  ! In/output variables:
  type(type_grid),                      intent(in   ) :: grid
  real(dp), dimension(grid%n1:grid%n2), intent(inout) :: d_grid_vec_partial
  real(dp),                             intent(in   ) :: sigma
  real(dp),                             intent(in   ) :: fill_value

  ! Local variables:
  character(len=*), parameter           :: routine_name = 'extrapolate_fillvalue_Gaussian_grid_2D'
  real(dp), dimension(:,:), allocatable :: d_grid_vec_partial_3D
  integer                               :: n

  ! Add routine to path
  call init_routine( routine_name)

  allocate( d_grid_vec_partial_3D( grid%n1: grid%n2, 1))
  do n = grid%n1, grid%n2
    d_grid_vec_partial_3D( n,1) = d_grid_vec_partial( n)
  end do
  call extrapolate_fillvalue_Gaussian_grid_3D( grid, d_grid_vec_partial_3D, sigma, fill_value)
  do n = grid%n1, grid%n2
    d_grid_vec_partial( n) = d_grid_vec_partial_3D( n,1)
  end do

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine extrapolate_fillvalue_Gaussian_grid_2D

subroutine extrapolate_fillvalue_Gaussian_grid_3D( grid, d_grid_vec_partial, sigma, fill_value)
  ! Replace fillvalues by Gaussian extrapolation (with radius sigma) of the rest of the data

  ! In/output variables:
  type(type_grid),          intent(in   ) :: grid
  real(dp), dimension(:,:), intent(inout) :: d_grid_vec_partial
  real(dp),                 intent(in   ) :: sigma
  real(dp),                 intent(in   ) :: fill_value

  ! Local variables:
  character(len=*), parameter             :: routine_name = 'extrapolate_fillvalue_Gaussian_grid_3D'
  real(dp), dimension(:,:,:), allocatable :: d_grid

  ! Add routine to path
  call init_routine( routine_name)

  ! Not feasibly parallelisable, just let the primary do all the work
  if (par%primary) allocate( d_grid( grid%nx, grid%ny, size( d_grid_vec_partial,2)))
  call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)
  if (par%primary) call extrapolate_fillvalue_Gaussian_grid_3D_primary( grid, d_grid, sigma, fill_value)
  call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid)

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine extrapolate_fillvalue_Gaussian_grid_3D

subroutine extrapolate_fillvalue_Gaussian_grid_3D_primary( grid, d, sigma, fill_value)
  ! Replace fillvalues by Gaussian extrapolation (with radius sigma) of the rest of the data

  ! In/output variables:
  type(type_grid),            intent(in   ) :: grid
  real(dp), dimension(:,:,:), intent(inout) :: d
  real(dp),                   intent(in   ) :: sigma
  real(dp),                   intent(in   ) :: fill_value

  ! Local variables:
  character(len=*), parameter           :: routine_name = 'extrapolate_fillvalue_Gaussian_grid_3D_primary'
  integer                               :: nz
  integer                               :: n_search
  real(dp)                              :: r_search, sigma_scaled
  integer                               :: i,j
  integer                               :: n_fill
  integer                               :: it_floodfill
  integer, dimension(grid%nx,grid%ny)   :: map
  integer, dimension(grid%nx*grid%ny,2) :: stack1, stack2
  integer                               :: stackN1, stackN2
  integer                               :: si,ii,jj

  ! Add routine to path
  call init_routine( routine_name)

  nz = size( d,3)

  sigma_scaled = sigma / grid%dx
  r_search     = sigma_scaled * 3._dp
  n_search     = ceiling( r_search)

  if (.not. par%primary) call crash('should only be called from the primary')
  if (size( d,1) /= grid%nx .or. size( d,2) /= grid%ny) call crash('array size does not match grid')

  ! Quickly check how many cells need filling
  n_fill = 0
  do i = 1, grid%nx
    do j = 1, grid%ny
      if (any( d( i,j,:) == fill_value)) then
        n_fill = n_fill + 1
        ! Safety
        if (.not. all( d( i,j,:) == fill_value)) then
          call crash('not all layers of d have the same to-be-filled cells')
        end if
      end if
    end do
  end do

  if (n_fill == 0) then
    call finalise_routine( routine_name)
    return
  elseif (n_fill == grid%nx * grid%ny) then
    call crash('d consist entirely of fill_values')
  end if

  ! Use a flood-fill style algorithm to perform the extrapolation

  ! Initialise map and stacks
  map = 0
  do i = 1, grid%nx
    do j = 1, grid%ny
      if (d( i,j,1) /= fill_value) map( i,j) = 2
    end do
  end do

  stack1  = 0
  stackN1 = 0
  do i = 1, grid%nx
    do j = 1, grid%ny
      if (map( i,j) == 0) then
        if (is_empty_next_to_filled( map, grid%nx, grid%ny, i, j)) then
          map( i,j) = 1
          stackN1 = stackN1 + 1
          stack1( stackN1,:) = [i,j]
        end if
      end if
    end do
  end do

  it_floodfill = 0
  iterate_floodfill: do while (.true.)

    ! If the stack is empty, we're done
    if (stackN1 == 0) exit iterate_floodfill

    ! DENK DROM
    ! write(0,*) '      iteration ', it_floodfill, ': stackN1 = ', stackN1

    ! Safety
    it_floodfill = it_floodfill + 1
    if (it_floodfill > grid%nx * grid%ny) call crash('flood-fill iteration got stuck!')

    ! Fill all elements in the stack
    do si = 1, stackN1
      i = stack1( si,1)
      j = stack1( si,2)
      call fill_element( d, map, grid%nx, grid%ny, nz, i, j, n_search, r_search, sigma_scaled)
    end do

    ! Mark them as filled, and construct the new stack
    stackN2 = 0
    do si = 1, stackN1
      i = stack1( si,1)
      j = stack1( si,2)
      map( i,j) = 2
      ! Add all non-stacked neighbours of [i,j] to the new stack
      do ii = max( 1, i-1), min( grid%nx, i+1)
      do jj = max( 1, j-1), min( grid%ny, j+1)
        if (map( ii,jj) == 0) then
          map( ii,jj) = 1
          stackN2 = stackN2 + 1
          stack2( stackN2,:) = [ii,jj]
        end if
      end do
      end do
    end do

    ! Cycle stacks
    stackN1 = stackN2
    stack1( 1:stackN1,:) = stack2( 1:stackN1,:)

  end do iterate_floodfill

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine extrapolate_fillvalue_Gaussian_grid_3D_primary

function is_empty_next_to_filled( map,nx,ny,i,j) result( isso)
  integer, dimension(:,:), intent(in) :: map
  integer,                 intent(in) :: nx,ny,i,j
  logical :: isso
  isso = (map( i,j) == 0) .and. (any( map( max( 1,i-1): min( nx,i+1), max( 1,j-1): min( ny,j+1)) == 2))
end function is_empty_next_to_filled

subroutine fill_element( d, map, nx, ny, nz, i, j, n_search, r_search, sigma_scaled)

  ! In/output variables:
  real(dp), dimension(nx,ny,nz), intent(inout) :: d
  integer,  dimension(nx,ny),    intent(in   ) :: map
  integer,                       intent(in   ) :: nx, ny, nz
  integer,                       intent(in   ) :: i,j
  integer,                       intent(in   ) :: n_search
  real(dp),                      intent(in   ) :: r_search
  real(dp),                      intent(in   ) :: sigma_scaled

  ! Local variables:
  integer                 :: ii, jj
  real(dp)                :: dist, w, sum_w
  real(dp), dimension(nz) :: sum_filled_near

  ! Average over the grid cells that lie within r_search of [i,j]
  sum_w           = 0._dp
  sum_filled_near = 0._dp
  do ii = max( 1, i-n_search), min( nx, i+n_search)
  do jj = max( 1, j-n_search), min( ny, j+n_search)
    dist = hypot( real( ii-i,dp), real( jj-j,dp))
    if (dist <= r_search .and. map( ii,jj) == 2) then
      w = exp( -0.5_dp * (dist / sigma_scaled)**2)
      sum_w           = sum_w           + 1
      sum_filled_near = sum_filled_near + d( ii,jj,:)
    end if
  end do
  end do

  d( i,j,:) = sum_filled_near / sum_w

end subroutine fill_element

end module smooth_gridded_data
