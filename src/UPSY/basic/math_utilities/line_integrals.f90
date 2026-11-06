module line_integrals

  ! Some simple line integrals (used in conservative remapping)

  use precisions, only: dp

  implicit none

  private

  public :: line_integral_xdy, line_integral_mxydx, line_integral_xydy, &
    calc_line_integrals_around_polygon

contains

  pure function line_integral_xdy(   p, q, tol_dist) result( I_pq)
    ! Calculate the line integral x dy from p to q

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: p, q
    real(dp)              , intent(in) :: tol_dist
    real(dp)                           :: I_pq

    ! Local variables:
    real(dp) :: xp, yp, xq, yq, dx, dy

    xp = p( 1)
    yp = p( 2)
    xq = q( 1)
    yq = q( 2)

    if (abs( yp-yq) < tol_dist) then
      I_pq = 0._dp
      return
    end if

    dx = q( 1) - p( 1)
    dy = q( 2) - p( 2)

    I_pq = xp*dy - yp*dx + (dx / (2._dp*dy)) * (yq**2 - yp**2)

  end function line_integral_xdy

  pure function line_integral_mxydx( p, q, tol_dist) result( I_pq)
    ! Calculate the line integral -xy dx from p to q

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: p, q
    real(dp)              , intent(in) :: tol_dist
    real(dp)                           :: I_pq

    ! Local variables:
    real(dp) :: xp, yp, xq, yq, dx, dy

    xp = p( 1)
    yp = p( 2)
    xq = q( 1)
    yq = q( 2)

    if (abs( xp-xq) < tol_dist) then
      I_pq = 0._dp
      return
    end if

    dx = q( 1) - p( 1)
    dy = q( 2) - p( 2)

    I_pq = (1._dp/2._dp * (xp*dy/dx - yp) * (xq**2-xp**2)) - (1._dp/3._dp * dy/dx * (xq**3-xp**3))

  end function line_integral_mxydx

  pure function line_integral_xydy(  p, q, tol_dist) result( I_pq)
    ! Calculate the line integral xy dy from p to q

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: p, q
    real(dp)              , intent(in) :: tol_dist
    real(dp)                           :: I_pq

    ! Local variables:
    real(dp) :: xp, yp, xq, yq, dx, dy

    xp = p( 1)
    yp = p( 2)
    xq = q( 1)
    yq = q( 2)

    if (abs( yp-yq) < tol_dist) then
      I_pq = 0._dp
      return
    end if

    dx = q( 1 ) - p( 1)
    dy = q( 2 ) - p( 2)

    I_pq = (1._dp/2._dp * (xp - yp*dx/dy) * (yq**2-yp**2)) + (1._dp/3._dp * dx/dy * (yq**3-yp**3))

  end function line_integral_xydy

  subroutine calc_line_integrals_around_polygon( poly, LI_xdy, LI_mxydx, LI_xydy, tol_dist)

    ! In/output variables:
    real(dp), dimension(:,:), intent(in   ) :: poly
    real(dp),                 intent(  out) :: LI_xdy, LI_mxydx, LI_xydy
    real(dp),                 intent(in   ) :: tol_dist

    ! Local variables:
    integer                :: n, i1, i2
    real(dp), dimension(2) :: p,q

    n = size( poly,1)

    LI_xdy   = 0;
    LI_mxydx = 0;
    LI_xydy  = 0;

    do i1 = 1, n
      if (i1 == n) then
        i2 = 1
      else
        i2 = i1 + 1
      end if

      p = poly( i1,:)
      q = poly( i2,:)

      LI_xdy   = LI_xdy   + line_integral_xdy(   p, q, tol_dist)
      LI_mxydx = LI_mxydx + line_integral_mxydx( p, q, tol_dist)
      LI_xydy  = LI_xydy  + line_integral_xydy(  p, q, tol_dist)

    end do

  end subroutine

end module line_integrals
