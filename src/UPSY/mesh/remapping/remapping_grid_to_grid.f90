module remapping_grid_to_grid

  use mpi_basic, only: par
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use grid_types, only: type_grid
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use remapping_types, only: type_map

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
    character(len=*), parameter :: routine_name = 'map_from_xy_grid_to_xy_grid_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! ! Allocate shared memory
    ! CALL allocate_shared_dp_2D(  ny_src, nx_src, ddx_src,              wddx_src             )
    ! CALL allocate_shared_dp_2D(  ny_src, nx_src, ddy_src,              wddy_src             )
    ! CALL allocate_shared_int_2D( ny_dst, nx_dst, mask_dst_outside_src, wmask_dst_outside_src)

    ! ! Find grid spacings
    ! dx_src = x_src(2) - x_src(1)
    ! dy_src = y_src(2) - y_src(1)
    ! dx_dst = x_dst(2) - x_dst(1)
    ! dy_dst = y_dst(2) - y_dst(1)
    ! Ad = dx_dst * dy_dst

    ! ! If the grids are equal, the solution is trivial; just copy the data
    ! IF (dx_src == dx_dst .AND. dy_src == dy_dst .AND. nx_src == nx_dst .AND. ny_src == ny_dst) THEN
    !   CALL partition_list( nx_dst, par%i, par%n, i1, i2)
    !   d_dst( :,i1:i2) = d_src( :,i1:i2)
    !   CALL sync
    !   CALL deallocate_shared(  wddx_src             )
    !   CALL deallocate_shared(  wddy_src             )
    !   CALL deallocate_shared(  wmask_dst_outside_src)
    !   CALL finalise_routine( routine_name)
    !   RETURN
    ! END IF

    ! ! Find overlaps between grids
    ! ! dst cell [i,j] overlaps with src cells [il_src( i):iu_src( i), jl_src( j):ju_src( j)]

    ! ! il
    ! DO i_dst = 1, nx_dst

    !   IF (i_dst == 1) THEN
    !     il_src( i_dst) = 1
    !   ELSE
    !     il_src( i_dst) = il_src( i_dst-1)
    !   END IF

    !   DO WHILE (x_src( il_src( i_dst)) + dx_src / 2._dp < x_dst( i_dst) - dx_dst / 2._dp)
    !     il_src( i_dst) = il_src( i_dst) + 1
    !     IF (il_src( i_dst) >= nx_src) THEN
    !       EXIT
    !     END IF
    !   END DO
    !   il_src( i_dst) = MAX( 1, MIN( nx_src, il_src( i_dst) ))

    ! END DO ! DO i_dst = 1, nx_dst

    ! il_src = il_src - 1
    ! il_src = MAX( 1, MIN( nx_src, il_src))

    ! ! iu
    ! DO i_dst = nx_dst, 1, -1

    !   IF (i_dst == nx_dst) THEN
    !     iu_src( i_dst) = nx_src
    !   ELSE
    !     iu_src( i_dst) = il_src( i_dst+1)
    !   END IF

    !   DO WHILE (x_src( iu_src( i_dst)) - dx_src / 2._dp > x_dst( i_dst) + dx_dst / 2._dp)
    !     iu_src( i_dst) = iu_src( i_dst) - 1
    !     IF (iu_src( i_dst) <= 1) THEN
    !       EXIT
    !     END IF
    !   END DO
    !   iu_src( i_dst) = MAX( 1, MIN( nx_src, iu_src( i_dst) ))

    ! END DO ! DO i_dst = nx_dst, 1, -1

    ! iu_src = iu_src + 1
    ! iu_src = MAX( 1, MIN( nx_src, iu_src))

    ! ! jl
    ! DO j_dst = 1, ny_dst

    !   IF (j_dst == 1) THEN
    !     jl_src( j_dst) = 1
    !   ELSE
    !     jl_src( j_dst) = jl_src( j_dst-1)
    !   END IF

    !   DO WHILE (y_src( jl_src( j_dst)) + dx_src / 2._dp < y_dst( j_dst) - dx_dst / 2._dp)
    !     jl_src( j_dst) = jl_src( j_dst) + 1
    !     IF (jl_src( j_dst) >= ny_src) THEN
    !       EXIT
    !     END IF
    !   END DO
    !   jl_src( j_dst) = MAX( 1, MIN( ny_src, jl_src( j_dst) ))

    ! END DO ! DO j_dst = 1, ny_dst

    ! jl_src = jl_src - 1
    ! jl_src = MAX( 1, MIN( ny_src, jl_src))

    ! ! ju
    ! DO j_dst = ny_dst, 1, -1

    !   IF (j_dst == ny_dst) THEN
    !     ju_src( j_dst) = ny_src
    !   ELSE
    !     ju_src( j_dst) = jl_src( j_dst+1)
    !   END IF

    !   DO WHILE (y_src( ju_src( j_dst)) - dx_src / 2._dp > y_dst( j_dst) + dx_dst / 2._dp)
    !     ju_src( j_dst) = ju_src( j_dst) - 1
    !     IF (ju_src( j_dst) <= 1) THEN
    !       EXIT
    !     END IF
    !   END DO
    !   ju_src( j_dst) = MAX( 1, MIN( ny_src, ju_src( j_dst) ))

    ! END DO ! DO j_dst = ny_dst, 1, -1

    ! ju_src = ju_src + 1
    ! ju_src = MAX( 1, MIN( ny_src, ju_src))

    ! ! Get derivatives of d_src
    ! CALL partition_list( nx_src, par%i, par%n, i1, i2)
    ! DO i = MAX(2,i1), MIN(nx_src-1,i2)
    ! DO j = 2, ny_src-1
    !   ddx_src( j,i) = (d_src( j,i+1) - d_src( j,i-1)) / (2._dp * dx_src)
    !   ddy_src( j,i) = (d_src( j+1,i) - d_src( j-1,i)) / (2._dp * dy_src)
    ! END DO
    ! END DO
    ! CALL sync

    ! ! Find parallelisation domains
    ! CALL partition_list( nx_dst, par%i, par%n, i1, i2)
    ! CALL partition_list( ny_dst, par%i, par%n, j1, j2)

    ! DO i = i1, i2
    ! DO j = 1, ny_dst

    !   d_dst( j,i) = 0._dp
    !   Asum        = 0._dp

    !   ! If this dst cell lies (partly) outside of the src grid, mark it as such;
    !   ! in that case, use nearest-neighbour extrapolation instead of conservative remapping
    !   IF (x_dst( i) - dx_dst/2._dp < MINVAL( x_src) - dx_src/2._dp .OR. &
    !       x_dst( i) + dx_dst/2._dp > MAXVAL( x_src) + dx_src/2._dp .OR. &
    !       y_dst( j) - dx_dst/2._dp < MINVAL( y_src) - dx_src/2._dp .OR. &
    !       y_dst( j) + dx_dst/2._dp > MAXVAL( y_src) + dx_src/2._dp) THEN
    !     mask_dst_outside_src( j,i) = 1
    !     CYCLE
    !   ELSE
    !     mask_dst_outside_src( j,i) = 0
    !   END IF

    !   DO i_src = il_src( i), iu_src( i)
    !   DO j_src = jl_src( j), ju_src( j)

    !     xomin = MAX( x_dst( i) - dx_dst/2._dp, x_src( i_src) - dx_src/2._dp)
    !     xomax = MIN( x_dst( i) + dx_dst/2._dp, x_src( i_src) + dx_src/2._dp)
    !     yomin = MAX( y_dst( j) - dy_dst/2._dp, y_src( j_src) - dy_src/2._dp)
    !     yomax = MIN( y_dst( j) + dy_dst/2._dp, y_src( j_src) + dy_src/2._dp)

    !     IF (xomax <= xomin .OR. yomax <= yomin) CYCLE

    !     Asd  = (xomax - xomin) * (yomax - yomin)
    !     Asum = Asum + Asd

    !     w0  = Asd / Ad
    !     w1x = 1._dp / Ad * (line_integral_mxydx( [xomin,yomin], [xomax,yomin], 1E-9_dp) + &
    !                         line_integral_mxydx( [xomax,yomin], [xomax,yomax], 1E-9_dp) + &
    !                         line_integral_mxydx( [xomax,yomax], [xomin,yomax], 1E-9_dp) + &
    !                         line_integral_mxydx( [xomin,yomax], [xomin,yomin], 1E-9_dp)) - w0 * x_src( i_src)
    !     w1y = 1._dp / Ad * (line_integral_xydy(  [xomin,yomin], [xomax,yomin], 1E-9_dp) + &
    !                         line_integral_xydy(  [xomax,yomin], [xomax,yomax], 1E-9_dp) + &
    !                         line_integral_xydy(  [xomax,yomax], [xomin,yomax], 1E-9_dp) + &
    !                         line_integral_xydy(  [xomin,yomax], [xomin,yomin], 1E-9_dp)) - w0 * y_src( j_src)

    !     d_dst( j,i) = d_dst( j,i) + w0  * d_src(   j_src,i_src) + &
    !                                 w1x * ddx_src( j_src,i_src) + &
    !                                 w1y * ddy_src( j_src,i_src)

    !   END DO ! DO j_src = jr_src( j,1), jr_src( j,2)
    !   END DO ! DO i_src = ir_src( i,1), ir_src( i,2)

    !   ! Safety
    !   IF (ABS( 1._dp - Asum / Ad) > 1E-4_dp) THEN
    !     CALL crash('dst grid cell [{int_01},{int_02}] couldnt be completely filled! Asum = {dp_01}, Ad = {dp_02}', int_01 = j, int_02 = i, dp_01 = Asum, dp_02 = Ad)
    !   END IF

    ! END DO ! DO j = 1, ny_dst
    ! END DO ! DO i = i1, i2
    ! CALL sync

    ! ! Use nearest-neighbour extrapolation for dst cells outside of the src grid
    ! ! =========================================================================

    ! ! Find the range of grid cells that were mapped correctly
    ! igmin = 0
    ! igmax = 0
    ! jgmin = 0
    ! jgmax = 0

    ! j = INT( REAL(ny_dst,dp)/2._dp)
    ! DO i = 1, nx_dst
    !   IF (mask_dst_outside_src( j,i) == 0) THEN
    !     igmin = i
    !     EXIT
    !   END IF
    ! END DO
    ! DO i = nx_dst, 1, -1
    !   IF (mask_dst_outside_src( j,i) == 0) THEN
    !     igmax = i
    !     EXIT
    !   END IF
    ! END DO

    ! i = INT( REAL(nx_dst,dp)/2._dp)
    ! DO j = 1, ny_dst
    !   IF (mask_dst_outside_src( j,i) == 0) THEN
    !     jgmin = j
    !     EXIT
    !   END IF
    ! END DO
    ! DO j = ny_dst, 1, -1
    !   IF (mask_dst_outside_src( j,i) == 0) THEN
    !     jgmax = j
    !     EXIT
    !   END IF
    ! END DO

    ! ! Corners
    ! IF (par%master) THEN
    !   ! Southwest
    !   d_dst( 1      :jgmin-1 ,1      :igmin-1) = d_dst( jgmin,igmin)
    !   ! Southeast
    !   d_dst( 1      :jgmin-1 ,igmax+1:nx_dst ) = d_dst( jgmin,igmax)
    !   ! Northwest
    !   d_dst( jgmax+1:ny_dst  ,1      :igmin-1) = d_dst( jgmax,igmin)
    !   ! Northeast
    !   d_dst( jgmax+1:ny_dst  ,igmax+1:nx_dst ) = d_dst( jgmax,igmax)
    ! END IF ! IF (par%master) THEN
    ! CALL sync

    ! ! Borders
    ! DO i = MAX(i1,igmin), MIN(i2,igmax)
    !   ! South
    !   d_dst( 1      :jgmin-1,i) = d_dst( jgmin,i)
    !   ! North
    !   d_dst( jgmax+1:ny_dst ,i) = d_dst( jgmax,i)
    ! END DO
    ! DO j = MAX(j1,jgmin), MIN(j2,jgmax)
    !   ! West
    !   d_dst( j,1      :igmin-1) = d_dst( j,igmin)
    !   ! East
    !   d_dst( j,igmax+1:nx_dst ) = d_dst( j,igmax)
    ! END DO
    ! CALL sync

    ! ! Clean up after yourself
    ! CALL deallocate_shared( wddx_src             )
    ! CALL deallocate_shared( wddy_src             )
    ! CALL deallocate_shared( wmask_dst_outside_src)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_from_xy_grid_to_xy_grid_2D

end module remapping_grid_to_grid
