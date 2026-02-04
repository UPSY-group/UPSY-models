module string_module

  use precisions, only: dp
  use mpi_basic, only: par

  implicit none

  private

  public :: type_string_utilities
  public :: separate_strings_by_double_vertical_bars, colour_string, &
    insert_val_into_string_int, insert_val_into_string_dp, capitalise_string, &
    remove_leading_spaces, str2int, int2str_with_leading_zeros, strrep

  type type_string_utilities
      private
    contains
      private
      procedure, public, nopass :: separate_strings_by_double_vertical_bars
      procedure, public, nopass :: colour_string
      procedure, public, nopass :: insert_val_into_string_int
      procedure, public, nopass :: insert_val_into_string_dp
      procedure, public, nopass :: capitalise_string
      procedure, public, nopass :: remove_leading_spaces
      procedure, public, nopass :: str2int
      procedure, public, nopass :: int2str_with_leading_zeros
      procedure, public, nopass :: strrep
  end type type_string_utilities

  logical :: do_colour_strings = .true.

contains

  pure subroutine separate_strings_by_double_vertical_bars( str_list, strs)
    !< Take a list of items separated by double vertical bars ("||"),
    !< and return them as separate strings

    character(len=*),                            intent(in   ) :: str_list
    character(len=*), dimension(:), allocatable, intent(  out) :: strs
    integer                       :: i, n_doublebars,ii
    character(len=:), allocatable :: str_list_redux

    ! Exception
    if (len_trim( str_list) == 0) then
      allocate( strs(0))
      return
    end if

    ! Count number of instances of "||"
    n_doublebars = 0
    do i = 1, len_trim( str_list)-1
      if (str_list( i:i+1) == '||') n_doublebars = n_doublebars + 1
    end do

    allocate( strs( n_doublebars+1))

    str_list_redux = str_list
    do i = 1, n_doublebars
      ii = index( str_list_redux,'||')
      strs(i) = str_list_redux( 1:ii-1)
      str_list_redux = str_list_redux( ii+2:len_trim( str_list_redux))
    end do
    strs( n_doublebars+1) = str_list_redux

  end subroutine separate_strings_by_double_vertical_bars

  pure function colour_string( str, col) result( str_col)
    !< Add colour to a string for writing to the terminal

    character(len=*),  intent(in) :: str, col
    character(len=:), allocatable :: str_col

    ! Safety: optionally don't do colours
    if (.not. do_colour_strings) then
      str_col = str
      return
    end if

    select case (col)
    case ('gray')
      str_col = achar(27) // '[90m' // str // achar(27) // '[0m'
    case ('red')
      str_col = achar(27) // '[91m' // str // achar(27) // '[0m'
    case ('green')
      str_col = achar(27) // '[92m' // str // achar(27) // '[0m'
    case ('yellow')
      str_col = achar(27) // '[93m' // str // achar(27) // '[0m'
    case ('blue')
      str_col = achar(27) // '[94m' // str // achar(27) // '[0m'
    case ( 'pink')
      str_col = achar(27) // '[95m' // str // achar(27) // '[0m'
    case ('light blue','light_blue')
      str_col = achar(27) // '[96m' // str // achar(27) // '[0m'
    case default
      str_col = str
    end select

  end function colour_string

  pure function insert_val_into_string_int( str, marker, val) result( str_new)
    !< Replace marker in str with val (where val is an integer)

    ! Example: str    = 'Johnny has {int_01} apples.'
    !          marker = '{int_01}'
    !          val    = 5
    !
    ! This returns: str = 'Johnny has 5 apples'

    character(len=*),  intent(in) :: str
    character(len=*),  intent(in) :: marker
    integer ,          intent(in) :: val
    character(len=:), allocatable :: str_new

    integer                       :: ci
    integer                       :: nc
    character(len=4)              :: fmt
    character(:), allocatable     :: val_str
    integer                       :: len_str, len_marker

    ! Find position ci in str where i_str occurs
    ci = index( str, marker)

    ! Safety
    if (ci == 0) error stop 'insert_val_into_string_int: couldnt find marker "' &
      // trim( marker) // '" in string "' // trim( str) // '"!'

    ! Write val to a string
    if     (abs( val) < 10) then
      nc = 1
    elseif (abs( val) < 100) then
      nc = 2
    elseif (abs( val) < 1000) then
      nc = 3
    elseif (abs( val) < 10000) then
      nc = 4
    elseif (abs( val) < 100000) then
      nc = 5
    elseif (abs( val) < 1000000) then
      nc = 6
    elseif (abs( val) < 10000000) then
      nc = 7
    elseif (abs( val) < 100000000) then
      nc = 8
    else
      error stop 'insert_val_into_string_int: only accepts integers up to 99999999'
    end if

    ! Add room for a minus sign if needed
    if (val < 0) nc = nc + 1

    write( fmt,'(A,I1,A)') '(I', nc, ')'
    allocate( character( nc) :: val_str)
    write( val_str, fmt) val

    ! Find total string length right now
    len_str    = len( str)
    len_marker = len( marker)

    ! Insert the integer string into the string
    str_new = str( 1:ci-1) // val_str // str( ci+len_marker:len_str)

  end function insert_val_into_string_int

  pure function insert_val_into_string_dp( str, marker, val) result( str_new)
    !< Replace marker in str with val (where val is a double-precision number)

    ! Example: str    = 'Johnny weighs {dp_01} kg.'
    !          marker = '{dp_01}'
    !          val    = 57.098
    !
    ! This returns: str = 'Johnny weighs 57.098 kg'

    character(len=*),  intent(in) :: str
    character(len=*),  intent(in) :: marker
    real(dp),          intent(in) :: val
    character(len=:), allocatable :: str_new

    integer                       :: ci
    character(len=20)             :: val_str
    integer                       :: len_str, len_marker

    ! Find position ci in str where i_str occurs
    ci = index( str, marker)

    ! Safety
    if (ci == 0) error stop 'insert_val_into_string_dp: couldnt find marker "' // &
      trim( marker) // '" in string "' // trim( str) // '"!'

    ! Write val to a string
    write( val_str,'(E20.14)') val

    ! Find total string length right now
    len_str    = len( str)
    len_marker = len( marker)

    ! Insert the integer string into the string
    str_new = str(1:ci-1) // val_str // str(ci+len_marker:len_str)

  end function insert_val_into_string_dp

  pure function capitalise_string( str) result( str_new)
    !< Capitalise all letters in a string

    character(len=*), intent(in)  :: str
    character(len=:), allocatable :: str_new
    integer                       :: i, index_cap
    character(len=*), parameter   :: upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    character(len=*), parameter   :: lower = 'abcdefghijklmnopqrstuvwxyz'

    str_new = str
    do i = 1, len_trim( str_new)
      index_cap = index( lower, str_new( i:i))
      if (index_cap > 0) str_new( i:i) = upper( index_cap:index_cap)
    end do

  end function capitalise_string

  pure function remove_leading_spaces( str) result( str_new)
    !< Remove leading spaces from a character string

    character(len=*),  intent(in) :: str
    character(len=:), allocatable :: str_new

    str_new = str
    do while (str_new( 1:1) == ' ')
      str_new = str_new( 2:len( str_new))
      if (len( str_new) == 0) exit
    end do

  end function remove_leading_spaces

  function str2int( str, stat) result( int)
    !< Convert a string containing an integer number to an actual integer
    character(len=*), intent(in   ) :: str
    integer,          intent(  out) :: stat
    integer                         :: int
    read( str, *, iostat = stat) int
  end function str2int

  function int2str_with_leading_zeros( int, n) result( str)
    !< Convert integer to a character string (with leading zeros)
    integer, intent(in) :: int
    integer, intent(in) :: n
    character(len=n)    :: str
    integer             :: i
    character(len=4)    :: fmt

    write( fmt,'(A,I1,A)') '(I', n, ')'
    write( str,fmt) int

    do i = 1, n
      if (str( i:i) == ' ') str( i:i) = '0'
    end do

  end function int2str_with_leading_zeros

  pure function strrep( str, old, new) result( str_new)
    !< Replace all occurences in [str] of [old] with [new]

    character(len=*),  intent(in) :: str
    character(len=*),  intent(in) :: old, new
    character(len=:), allocatable :: str_new
    integer                       :: i, j, nit

    str_new = str
    nit = 0
    do while (index( str_new, old) > 0 .and. nit < len( str))
      nit = nit + 1
      i = index( str_new, old)
      j = i + len( old)
      str_new = str_new( 1:i-1) // trim(new) // str( j: len_trim( str_new))
    end do

  end function strrep

end module string_module