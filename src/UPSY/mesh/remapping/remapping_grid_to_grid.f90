module remapping_grid_to_grid

  use mpi_basic, only: par
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use grid_types, only: type_grid
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary, distribute_gridded_data_from_primary
  use line_integrals, only: line_integral_mxydx, line_integral_xydy

  implicit none

  private

  public :: map_from_xy_grid_to_xy_grid_2D

contains

  subroutine map_from_xy_grid_to_xy_grid_2D( grid_src, grid_dst, d_grid_src_vec_partial, d_grid_dst_vec_partial)
    !< Map data from one x/y-grid to another x/y-grid

    ! In/output variables
    type(type_grid),                                      intent(in   ) :: grid_src
    type(type_grid),                                      intent(in   ) :: grid_dst
    real(dp), dimension(grid_src%pai%i1:grid_src%pai%i2), intent(in   ) :: d_grid_src_vec_partial
    real(dp), dimension(grid_dst%pai%i1:grid_dst%pai%i2), intent(  out) :: d_grid_dst_vec_partial

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'map_from_xy_grid_to_xy_grid_2D'
    real(dp), dimension(:,:), allocatable :: d_grid_src
    real(dp), dimension(:,:), allocatable :: d_grid_dst
    real(dp)                              :: Ad
    integer,  dimension(grid_dst%nx)      :: il_src, iu_src
    integer,  dimension(grid_dst%ny)      :: jl_src, ju_src
    integer                               :: i_src, j_src, i_dst, j_dst, i,j
    real(dp), dimension(:,:), allocatable :: ddx_grid_src
    real(dp), dimension(:,:), allocatable :: ddy_grid_src
    real(dp)                              :: Asum, Asd, w0, w1x, w1y
    logical,  dimension(:,:), allocatable :: mask_dst_outside_src
    real(dp)                              :: xomin, xomax, yomin, yomax
    real(dp), dimension(2)                :: sw, se, nw, ne
    integer                               :: igmin, igmax, jgmin, jgmax

    ! Add routine to path
    call init_routine( routine_name)

    ! Let the primary do the work

    ! NOTE: not very elegant or efficient, but right now this code only has one call site,
    ! being the remapping of the GIA-ELRA spatially variable relaxation time from its
    ! NetCDF input grid to the GIA grid. That only happens during model initialisation, and
    ! involves two coarse-resolution grids, so it won't noticeably slow down the model.
    ! And this code will be replaced with the remapping overhaul that is planned for 2027 anyway...

    if (par%primary) then
      allocate( d_grid_src( grid_src%nx, grid_src%ny), source = 0._dp)
      allocate( d_grid_src( grid_src%nx, grid_src%ny), source = 0._dp)
    else
      allocate( d_grid_src( 0,0))
      allocate( d_grid_src( 0,0))
    end if
    call gather_gridded_data_to_primary( grid_src, d_grid_src_vec_partial, d_grid_src)

    if (par%primary) then

      Ad = grid_dst%dx**2

    ! Find overlaps between grids
    ! dst cell [i,j] overlaps with src cells [il_src( i):iu_src( i), jl_src( j):ju_src( j)]

    ! il
    do i_dst = 1, grid_dst%nx

      if (i_dst == 1) then
        il_src( i_dst) = 1
      else
        il_src( i_dst) = il_src( i_dst-1)
      end if

      do while (grid_src%x( il_src( i_dst)) + grid_src%dx / 2._dp < grid_dst%x( i_dst) - grid_dst%dx / 2._dp)
        il_src( i_dst) = il_src( i_dst) + 1
        if (il_src( i_dst) >= grid_src%nx) then
          exit
        end if
      end do
      il_src( i_dst) = max( 1, min( grid_src%nx, il_src( i_dst) ))

    end do

    il_src = il_src - 1
    il_src = max( 1, min( grid_src%nx, il_src))

    ! iu
    do i_dst = grid_dst%nx, 1, -1

      if (i_dst == grid_dst%nx) then
        iu_src( i_dst) = grid_src%nx
      else
        iu_src( i_dst) = il_src( i_dst+1)
      end if

      do while (grid_src%x( iu_src( i_dst)) - grid_src%dx / 2._dp > grid_dst%x( i_dst) + grid_dst%dx / 2._dp)
        iu_src( i_dst) = iu_src( i_dst) - 1
        if (iu_src( i_dst) <= 1) then
          exit
        end if
      end do
      iu_src( i_dst) = max( 1, min( grid_src%nx, iu_src( i_dst) ))

    end do

    iu_src = iu_src + 1
    iu_src = max( 1, min( grid_src%nx, iu_src))

    ! jl
    do j_dst = 1, grid_dst%ny

      if (j_dst == 1) then
        jl_src( j_dst) = 1
      else
        jl_src( j_dst) = jl_src( j_dst-1)
      end if

      do while (grid_src%y( jl_src( j_dst)) + grid_src%dx / 2._dp < grid_dst%y( j_dst) - grid_dst%dx / 2._dp)
        jl_src( j_dst) = jl_src( j_dst) + 1
        if (jl_src( j_dst) >= grid_src%ny) then
          exit
        end if
      end do
      jl_src( j_dst) = max( 1, min( grid_src%ny, jl_src( j_dst) ))

    end do

    jl_src = jl_src - 1
    jl_src = max( 1, min( grid_src%ny, jl_src))

    ! ju
    do j_dst = grid_dst%ny, 1, -1

      if (j_dst == grid_dst%ny) then
        ju_src( j_dst) = grid_src%ny
      else
        ju_src( j_dst) = jl_src( j_dst+1)
      end if

      do while (grid_src%y( ju_src( j_dst)) - grid_src%dx / 2._dp > grid_dst%y( j_dst) + grid_dst%dx / 2._dp)
        ju_src( j_dst) = ju_src( j_dst) - 1
        if (ju_src( j_dst) <= 1) then
          exit
        end if
      end do
      ju_src( j_dst) = max( 1, min( grid_src%ny, ju_src( j_dst) ))

    end do

    ju_src = ju_src + 1
    ju_src = max( 1, min( grid_src%ny, ju_src))

    ! Calculate derivatives of d_src
    allocate( ddx_grid_src( grid_src%nx, grid_src%ny), source = 0._dp)
    allocate( ddy_grid_src( grid_src%nx, grid_src%ny), source = 0._dp)
    do i_src = 2, grid_src%nx-1
    do j_src = 2, grid_src%ny-1
      ddx_grid_src( i_src,j_src) = (d_grid_src( i_src+1,j_src  ) - d_grid_src( i_src-1,j_src  )) / (2._dp * grid_src%dx)
      ddy_grid_src( i_src,j_src) = (d_grid_src( i_src  ,j_src+1) - d_grid_src( i_src  ,j_src-1)) / (2._dp * grid_src%dx)
    end do
    end do

    ! Perform 2nd-order conservative remapping
    allocate( mask_dst_outside_src( grid_dst%nx, grid_dst%ny), source = .false.)
    do i = 1, grid_dst%nx
    do j = 1, grid_dst%ny

      d_grid_dst( i,j) = 0._dp
      Asum             = 0._dp

      ! If this dst cell lies (partly) outside of the src grid, mark it as such;
      ! in that case, use nearest-neighbour extrapolation instead of conservative remapping
      if (grid_dst%x( i) - grid_dst%dx/2._dp < minval( grid_src%x) - grid_src%dx/2._dp .or. &
          grid_dst%x( i) + grid_dst%dx/2._dp > maxval( grid_src%x) + grid_src%dx/2._dp .or. &
          grid_dst%y( j) - grid_dst%dx/2._dp < minval( grid_src%y) - grid_src%dx/2._dp .or. &
          grid_dst%y( j) + grid_dst%dx/2._dp > maxval( grid_src%y) + grid_src%dx/2._dp) then
        mask_dst_outside_src( i,j) = .true.
        cycle
      end if

      do i_src = il_src( i), iu_src( i)
      do j_src = jl_src( j), ju_src( j)

        xomin = max( grid_dst%x( i) - grid_dst%dx/2._dp, grid_src%x( i_src) - grid_src%dx/2._dp)
        xomax = min( grid_dst%x( i) + grid_dst%dx/2._dp, grid_src%x( i_src) + grid_src%dx/2._dp)
        yomin = max( grid_dst%y( j) - grid_dst%dx/2._dp, grid_src%y( j_src) - grid_src%dx/2._dp)
        yomax = min( grid_dst%y( j) + grid_dst%dx/2._dp, grid_src%y( j_src) + grid_src%dx/2._dp)

        if (xomax <= xomin .or. yomax <= yomin) cycle

        Asd  = (xomax - xomin) * (yomax - yomin)
        Asum = Asum + Asd

        sw = [xomin, yomin]
        se = [xomax, yomin]
        nw = [xomin, yomax]
        ne = [xomax, yomax]

        w0  = Asd / Ad

        w1x = 1._dp / Ad * (line_integral_mxydx( sw, se, 1E-9_dp) + &
                            line_integral_mxydx( se, ne, 1E-9_dp) + &
                            line_integral_mxydx( ne, nw, 1E-9_dp) + &
                            line_integral_mxydx( nw, sw, 1E-9_dp)) - w0 * grid_src%x( i_src)

        w1y = 1._dp / Ad * (line_integral_xydy(  sw, se, 1E-9_dp) + &
                            line_integral_xydy(  se, ne, 1E-9_dp) + &
                            line_integral_xydy(  ne, nw, 1E-9_dp) + &
                            line_integral_xydy(  nw, sw, 1E-9_dp)) - w0 * grid_src%y( j_src)

        d_grid_dst( i,j) = d_grid_dst( i,j) + &
          w0  * d_grid_src  ( i_src, j_src) + &
          w1x * ddx_grid_src( i_src, j_src) + &
          w1y * ddy_grid_src( i_src, j_src)

      end do
      end do

      ! Safety
      if (abs( 1._dp - Asum / Ad) > 1E-4_dp) then
        call crash('dst grid cell [{int_01},{int_02}] couldnt be completely filled! Asum = {dp_01}, Ad = {dp_02}', int_01 = j, int_02 = i, dp_01 = Asum, dp_02 = Ad)
      end if

    end do
    end do

    ! Use nearest-neighbour extrapolation for dst cells outside of the src grid
    ! =========================================================================

    ! Find the range of grid cells that were mapped correctly
    igmin = 0
    igmax = 0
    jgmin = 0
    jgmax = 0

    j = int( real( grid_dst%ny,dp)/2._dp)
    do i = 1, grid_dst%nx
      if (.not. mask_dst_outside_src( i,j)) then
        igmin = i
        exit
      end if
    end do
    do i = grid_dst%nx, 1, -1
      if (.not. mask_dst_outside_src( i,j)) then
        igmax = i
        exit
      end if
    end do

    i = int( real( grid_dst%nx,dp)/2._dp)
    do j = 1, grid_dst%ny
      if (.not. mask_dst_outside_src( i,j)) then
        jgmin = j
        exit
      end if
    end do
    do j = grid_dst%ny, 1, -1
      if (.not. mask_dst_outside_src( i,j)) then
        jgmax = j
        exit
      end if
    end do

    ! Corners
    ! Southwest
    d_grid_dst( 1:igmin-1, 1:jgmin-1) = d_grid_dst( igmin, jgmin)
    ! Southeast
    d_grid_dst( igmax+1:grid_dst%nx, 1:jgmin-1) = d_grid_dst( igmax, jgmin)
    ! Northwest
    d_grid_dst( 1:igmin-1, jgmax+1:grid_dst%ny) = d_grid_dst( igmin, jgmax)
    ! Northeast
    d_grid_dst( igmax+1:grid_dst%nx, jgmax+1:grid_dst%ny) = d_grid_dst( igmax, jgmax)

    ! Borders
    do i = max( 1,igmin), min( grid_dst%nx,igmax)
      ! South
      d_grid_dst( i, 1:jgmin-1) = d_grid_dst( i, jgmin)
      ! North
      d_grid_dst( i, jgmax+1:grid_dst%ny) = d_grid_dst( i, jgmax)
    end do
    do j = max( 1,jgmin), min( grid_dst%ny,jgmax)
      ! West
      d_grid_dst( 1:igmin-1, j) = d_grid_dst( igmin, j)
      ! East
      d_grid_dst( igmax+1:grid_dst%nx, j) = d_grid_dst( igmax, j)
    end do

    end if ! if (par%primary) then

    call distribute_gridded_data_from_primary( grid_dst, d_grid_dst_vec_partial, d_grid_dst)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_xy_grid_to_xy_grid_2D

end module remapping_grid_to_grid
