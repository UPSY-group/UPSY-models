module switch_array_elements

  ! Switch individual elements or entire rows in 2-D arrays

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: crash

  implicit none

  private

  public :: switch_rows, switch_entries_by_value

  interface switch_rows
    procedure :: switch_rows_int_1D
    procedure :: switch_rows_int_2D
    procedure :: switch_rows_dp_1D
    procedure :: switch_rows_dp_2D
  end interface switch_rows

  interface switch_entries_by_value
    procedure :: switch_entries_by_value_int_2D
    procedure :: switch_entries_by_value_dp_2D
  end interface switch_entries_by_value

contains

  subroutine switch_rows_int_1D( d, i1, i2)
    ! Switch rows i1 and i2 in d

    ! In/output variables:
    integer, dimension(:), intent(inout) :: d
    integer,               intent(in   ) :: i1, i2

    ! Local variables:
    integer :: a, b

    a = d( i1)
    b = d( i2)
    d( i1) = b
    d( i2) = a

  end subroutine switch_rows_int_1D

  subroutine switch_rows_int_2D( d, i1, i2)
    ! Switch rows i1 and i2 in d

    ! In/output variables:
    integer, dimension(:,:), intent(inout) :: d
    integer,                 intent(in   ) :: i1, i2

    ! Local variables:
    integer, dimension(size(d,2)) :: a, b

    a = d( i1,:)
    b = d( i2,:)
    d( i1,:) = b
    d( i2,:) = a

  end subroutine switch_rows_int_2D

  subroutine switch_rows_dp_1D( d, i1, i2)
    ! Switch rows i1 and i2 in d

    ! In/output variables:
    real(dp), dimension(:), intent(inout) :: d
    integer,                intent(in   ) :: i1, i2

    ! Local variables:
    real(dp) :: a, b

    a = d( i1)
    b = d( i2)
    d( i1) = b
    d( i2) = a

  end subroutine switch_rows_dp_1D

  subroutine switch_rows_dp_2D( d, i1, i2)
    ! Switch rows i1 and i2 in d

    ! In/output variables:
    real(dp), dimension(:,:), intent(inout) :: d
    integer,                  intent(in   ) :: i1, i2

    ! Local variables:
    real(dp), dimension(size(d,2)) :: a, b

    a = d( i1,:)
    b = d( i2,:)
    d( i1,:) = b
    d( i2,:) = a

  end subroutine switch_rows_dp_2D

  subroutine switch_entries_by_value_int_2D( d, d1, d2)
    ! Switch entries d1 and d2 in d (i.e. d_new( d==d1) = d2, d_new( d==d2) = d1, elsewhere d_new = d

    ! In/output variables:
    integer, dimension(:,:), intent(inout) :: d
    integer,                 intent(in   ) :: d1, d2

    ! Local variables:
    integer, dimension(:,:), allocatable :: ij_d1, ij_d2
    integer                              :: n_d1, n_d2
    integer                              :: i,j,ii

    ! List all places where d1 and d2 occur in d
    allocate( ij_d1( size(d,1)*size(d,2), 2), source = 0)
    allocate( ij_d2( size(d,1)*size(d,2), 2), source = 0)
    n_d1 = 0
    n_d2 = 0

    do i = 1, size( d,1)
      do j = 1, size( d,2)
        if (d( i,j) == d1) then
          n_d1 = n_d1 + 1
          ij_d1( n_d1,:) = [i,j]
        elseif (d( i,j) == d2) then
          n_d2 = n_d2 + 1
          ij_d2( n_d2,:) = [i,j]
        end if
      end do
    end do

    ! Replace all (old) occurences of d1 with d2
    do ii = 1, n_d1
      i = ij_d1( ii,1)
      j = ij_d1( ii,2)
      d( i,j) = d2
    end do

    ! Replace all (old) occurences of d2 with d1
    do ii = 1, n_d2
      i = ij_d2( ii,1)
      j = ij_d2( ii,2)
      d( i,j) = d1
    end do

  end subroutine switch_entries_by_value_int_2D

  subroutine switch_entries_by_value_dp_2D( d, d1, d2)
    ! Switch entries d1 and d2 in d (i.e. d_new( d==d1) = d2, d_new( d==d2) = d1, elsewhere d_new = d

    ! In/output variables:
    real(dp), dimension(:,:), intent(inout) :: d
    integer,                  intent(in   ) :: d1, d2

    ! Local variables:
    integer, dimension(:,:), allocatable :: ij_d1, ij_d2
    integer                              :: n_d1, n_d2
    integer                              :: i,j,ii

    ! List all places where d1 and d2 occur in d
    allocate( ij_d1( size(d,1)*size(d,2), 2), source = 0)
    allocate( ij_d2( size(d,1)*size(d,2), 2), source = 0)
    n_d1 = 0
    n_d2 = 0

    do i = 1, size( d,1)
      do j = 1, size( d,2)
        if (d( i,j) == d1) then
          n_d1 = n_d1 + 1
          ij_d1( n_d1,:) = [i,j]
        elseif (d( i,j) == d2) then
          n_d2 = n_d2 + 1
          ij_d2( n_d2,:) = [i,j]
        end if
      end do
    end do

    ! Replace all (old) occurences of d1 with d2
    do ii = 1, n_d1
      i = ij_d1( ii,1)
      j = ij_d1( ii,2)
      d( i,j) = d2
    end do

    ! Replace all (old) occurences of d2 with d1
    do ii = 1, n_d2
      i = ij_d2( ii,1)
      j = ij_d2( ii,2)
      d( i,j) = d1
    end do

  end subroutine switch_entries_by_value_dp_2D

end module switch_array_elements
