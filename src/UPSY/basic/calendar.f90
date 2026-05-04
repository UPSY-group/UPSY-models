module calendar

  ! Routines to convert times in years to CF times (days since ...)

! ===== Preamble =====
! ====================

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning


  implicit none

  private

  public :: convert_time_to_days

  contains

    subroutine convert_time_to_days( time, days, days_bounds)

      ! In/output variables:
      real(dp),                          intent(in   ) :: time
      real(dp),                          intent(  out) :: days
      real(dp), dimension(2), optional,  intent(  out) :: days_bounds

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'convert_time_to_days'
      real(dp)                       :: frac_year
      integer                        :: full_year
      integer                        :: i
      real(dp), parameter            :: eps = 1.e-8

      ! Add routine to path
      call init_routine( routine_name)

      ! Convert time to full_year integer
      full_year = nint(time)

      ! Check whether input time falls at start/end of year
      if (abs(time-full_year*1._dp) > eps) then
        call crash('Requested time to convert to days is not a full year')
      end if

      if (present(days_bounds)) then
        ! Determine days at 1 July of last year,
        ! Determine bounds to span the last year from 1 Jan to 1 Jan
        call convert_time_to_days_with_bounds( full_year, days, days_bounds)
      else
        ! Determine days at the actual current time at 1 Jan
        call convert_time_to_days_nobounds( full_year, days)
      end if

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine convert_time_to_days

    subroutine convert_time_to_days_nobounds( full_year, days)

      ! In/output variables:
      integer,                 intent(in   ) :: full_year
      real(dp),                intent(  out) :: days

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'convert_time_to_days_nobounds'
      integer                        :: i

      ! Add routine to path
      call init_routine( routine_name)

      ! Initialise to align with compliance checker
      days = -1._dp

      ! Count days for full years from 1850 to full year
      do i = 1850, full_year
        if (is_leap_year(i)) then
          days = days + 366._dp
        else
          days = days + 365._dp
        end if
      end do

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine convert_time_to_days_nobounds

    subroutine convert_time_to_days_with_bounds( full_year, days, days_bounds)

      ! In/output variables:
      integer,                 intent(in   ) :: full_year
      real(dp),                intent(  out) :: days
      real(dp), dimension(2),  intent(  out) :: days_bounds

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'convert_time_to_days_with_bounds'
      integer                        :: i
      real(dp)                       :: days_start, days_end

      ! Add routine to path
      call init_routine( routine_name)

      ! Initialise
      days_start = 0._dp

      ! Count days for full years from 1850 to start of last year
      do i = 1850, full_year- 1._dp
        if (is_leap_year(i)) then
          days_start = days_start + 366._dp
        else
          days_start = days_start + 365._dp
        end if
      end do

      ! Add days of current year to year end
      if (is_leap_year( full_year)) then
        days_end = days_start + 366._dp
      else
        days_end = days_start + 365._dp
      end if

      ! Parse to days_bounds
      days_bounds(1) = days_start
      days_bounds(2) = days_end

      ! Get center value at 1 July of last year
      days = days_end - 184._dp

      ! Finalise routine path
      call finalise_routine( routine_name)

   end subroutine convert_time_to_days_with_bounds

    function is_leap_year(y) result(is_leap)
        implicit none
        integer, intent(in) :: y
        logical :: is_leap
        
        is_leap = (mod(y, 4) == 0) .and. (.not. mod(y, 100) == 0 .or. mod(y, 400) == 0)
    end function is_leap_year

end module calendar
