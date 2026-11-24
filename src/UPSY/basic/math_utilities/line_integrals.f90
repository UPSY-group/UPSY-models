module line_integrals

  ! Some simple line integrals (used in conservative remapping)

  use precisions, only: dp

  implicit none

  private

  public :: line_integral_xdy, line_integral_mxydx, line_integral_xydy, &
    calc_line_integrals_around_polygon

contains

  pure function line_integral_xdy(   p, q) result( I_pq)
    ! Calculate the line integral x dy from p to q

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: p, q
    real(dp)                           :: I_pq

    I_pq = 0.5_dp * (p(1) + q(1)) * (q(2) - p(2))

  end function line_integral_xdy

  pure function line_integral_mxydx( p, q) result( I_pq)
    ! Calculate the line integral -xy dx from p to q

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: p, q
    real(dp)                           :: I_pq

    I_pq = -(q(1) - p(1)) * (2._dp*p(1)*p(2) + 2._dp*q(1)*q(2) + p(1)*q(2) + q(1)*p(2)) / 6._dp

  end function line_integral_mxydx

  pure function line_integral_xydy(  p, q) result( I_pq)
    ! Calculate the line integral xy dy from p to q

    ! In/output variables:
    real(dp), dimension(2), intent(in) :: p, q
    real(dp)                           :: I_pq

    I_pq =  (q(2) - p(2)) * (2._dp*p(1)*p(2) + 2._dp*q(1)*q(2) + p(1)*q(2) + q(1)*p(2)) / 6._dp

  end function line_integral_xydy

  subroutine calc_line_integrals_around_polygon( poly, LI_xdy, LI_mxydx, LI_xydy)

    ! In/output variables:
    real(dp), dimension(:,:), intent(in   ) :: poly
    real(dp),                 intent(  out) :: LI_xdy, LI_mxydx, LI_xydy

    ! Local variables:
    integer                :: n, i1, i2
    real(dp), dimension(2) :: p,q

    n = size( poly,1)

    LI_xdy   = 0
    LI_mxydx = 0
    LI_xydy  = 0

    do i1 = 1, n
      if (i1 == n) then
        i2 = 1
      else
        i2 = i1 + 1
      end if

      p = poly( i1,:)
      q = poly( i2,:)

      LI_xdy   = LI_xdy   + line_integral_xdy(   p, q)
      LI_mxydx = LI_mxydx + line_integral_mxydx( p, q)
      LI_xydy  = LI_xydy  + line_integral_xydy(  p, q)

    end do

  end subroutine

end module line_integrals
