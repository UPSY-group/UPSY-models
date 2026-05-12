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

    subroutine convert_time_to_days( time, days, days_bounds, calendar, allow_residual)

      ! In/output variables:
      real(dp),                          intent(in   ) :: time
      real(dp),                          intent(  out) :: days
      real(dp), dimension(2), optional,  intent(  out) :: days_bounds
      character(len=*), optional,        intent(in   ) :: calendar
      logical, optional,                 intent(in   ) :: allow_residual

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'convert_time_to_days'
      real(dp)                       :: frac_year
      integer                        :: full_year
      real(dp)                       :: res
      integer                        :: i
      real(dp), parameter            :: eps = 1.e-8
      character(len=1024)            :: calendar_applied
      logical                        :: residual_allowed

      ! Add routine to path
      call init_routine( routine_name)

      ! Determine calendar to use
      if (present(calendar)) then
        calendar_applied = calendar
      else
        calendar_applied = 'standard'
      end if

      ! Determine whether residual days (non-full year) are allowed. Default: not
      if (present(allow_residual)) then
        residual_allowed = allow_residual
      else
        residual_allowed = .false.
      end if

      ! Determine full_year integer and optional residual
      if (residual_allowed) then
        full_year = int(time) ! Rounded down
        res = time - full_year
      else
        full_year = nint(time) ! Closest full year
        res = 0._dp
      end if

      if (present(days_bounds)) then
        ! This is only used for output, which is only allowed at full years,
        ! so that the bounds can be defined properly over the full last year,
        ! and the centered time at mid-year last yera (1 July)

        ! Check whether input time falls at start/end of year
        if (abs(time-full_year*1._dp) > eps) then
          call crash('Requested time to convert to days is not a full year')
        end if

        call convert_time_to_days_with_bounds( full_year, days, days_bounds)
      else
        ! Determine days at the actual current time at 1 Jan
        call convert_time_to_days_nobounds( full_year, res, calendar_applied, days)
      end if

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine convert_time_to_days

    subroutine convert_time_to_days_nobounds( full_year, res, calendar, days)

      ! In/output variables:
      integer,                 intent(in   ) :: full_year
      real(dp),                intent(in   ) :: res
      character(len=*),        intent(in   ) :: calendar
      real(dp),                intent(  out) :: days

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'convert_time_to_days_nobounds'
      integer                        :: i

      ! Add routine to path
      call init_routine( routine_name)

      select case (calendar)
        case default
          call crash ('unknown calendar "' // trim(calendar) // '"')
        case ('standard')
          ! Determine days according to standard/gregorian calendar,
          ! including leap years

          ! Initialise with offset 1 day to align with ISMIP7 compliance checker,
          ! probably to make sure output is 31 Dec for a given year
          days = -1._dp
    
          ! Count days for full years from 1850 to full year
          do i = 1850, full_year
            if (is_leap_year(i)) then
              days = days + 366._dp
            else
              days = days + 365._dp
            end if
          end do

          ! Add approximate residual
          days = days + res * 365.24_dp

        case ('noleap', '365_day')
          ! Determine days according to years of 365 days without leap years
          days = (full_year-1850) * 365._dp

          ! Add residual
          days = days + res * 365._dp

        case ('360_day')
          ! Determine days according to years of 360 days (each month 30 days)
          days = (full_year-1850) * 360._dp

          ! Add residual
          days = days + res * 360._dp

      end select

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine convert_time_to_days_nobounds

    subroutine convert_time_to_days_with_bounds( full_year, days, days_bounds)
      ! This routine is specific for creating output, which we always want on
      ! the standard calendar, hence no other calendars are implemented
      ! This is tailored specifically to ISMIP7 output

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
      do i = 1850, full_year- 1
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
