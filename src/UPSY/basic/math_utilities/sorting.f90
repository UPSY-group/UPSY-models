module sorting

  ! Sorting routines

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine

  implicit none

contains

  subroutine quick_n_dirty_sort( f, ii)
    ! Inefficient but simple sorting algorithm
    ! Sort f ascending and return list of new indices

    ! In/output variables:
    real(dp), dimension(:), intent(inout) :: f
    integer,  dimension(:), intent(inout) :: ii

    ! Local variables:
    character(len=256), parameter :: routine_name = 'quick_n_dirty_sort'
    integer                       :: n, i, j, iia, iib
    real(dp)                      :: a, b

    ! Add routine to path
    call init_routine( routine_name)

    ! Number of elements
    n = size( f,1)

    ! Initialise current indices
    do i = 1, n
      ii( i) = i
    end do

    ! Sort
    do i = 1, n-1
      do j = i+1, n
        if (f( i) > f( j)) then

          a   = f(  i)
          b   = f(  j)
          iia = ii( i)
          iib = ii( j)

          f(  i) = b
          f(  j) = a
          ii( i) = iib
          ii( j) = iia

        end if
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine quick_n_dirty_sort

end module sorting
