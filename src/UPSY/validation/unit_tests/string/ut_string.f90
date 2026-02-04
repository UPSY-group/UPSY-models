module ut_string

  use UPSY_main, only: UPSY
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use precisions, only: dp
  use ut_basic, only: unit_test
  use mpi_basic, only: par
  use parameters, only: pi, NaN

  implicit none

  private

  public :: unit_tests_string_main

contains

  subroutine unit_tests_string_main( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'unit_tests_string_main'
    character(len=*), parameter   :: test_name_local = 'string'
    character(len=:), allocatable :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_separate_strings_by_double_vertical_bars( test_name)
    call test_colour_string                           ( test_name)
    call test_insert_val_into_string_int              ( test_name)
    call test_insert_val_into_string_dp               ( test_name)
    call test_capitalise_string                       ( test_name)
    call test_remove_leading_spaces                   ( test_name)
    call test_str2int                                 ( test_name)
    call test_int2str_with_leading_zeros              ( test_name)
    call test_strrep                                  ( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_string_main

  subroutine test_separate_strings_by_double_vertical_bars( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter                    :: routine_name = 'test_separate_strings_by_double_vertical_bars'
    character(len=*), parameter                    :: test_name_local = 'separate_strings_by_double_vertical_bars'
    character(len=:), allocatable                  :: test_name
    character(len=:), allocatable                  :: str
    character(len=1024), dimension(:), allocatable :: strs

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Some test cases
    str = 'john||william||dede||a||ladida,di+da&  di|da'
    call UPSY%stru%separate_strings_by_double_vertical_bars( str, strs)
    call unit_test( &
      size(strs,1) == 5 .and. &
      strs(1)=='john' .and. &
      strs(2)=='william' .and. &
      strs(3)=='dede' .and. &
      strs(4)=='a' .and. &
      strs(5)=='ladida,di+da&  di|da', test_name // '/multiple_strings')
    deallocate( str, strs)

    str = 'pete'
    call UPSY%stru%separate_strings_by_double_vertical_bars( str, strs)
    call unit_test( &
      size(strs,1) == 1 .and. &
      strs(1)=='pete', test_name // '/single_string')
    deallocate( str, strs)

    str = ''
    call UPSY%stru%separate_strings_by_double_vertical_bars( str, strs)
    call unit_test( size(strs,1) == 0, test_name // '/empty_string')
    deallocate( str, strs)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_separate_strings_by_double_vertical_bars

  subroutine test_colour_string( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'test_colour_string'
    character(len=*), parameter   :: test_name_local = 'colour_string'
    character(len=:), allocatable :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be gray', 'gray')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be red', 'red')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be green', 'green')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be yellow', 'yellow')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be blue', 'blue')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be pink', 'pink')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be light blue', 'light blue')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be white', 'anything else')

    call unit_test( .true., test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_colour_string

  subroutine test_insert_val_into_string_int( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'test_insert_val_into_string_int'
    character(len=*), parameter   :: test_name_local = 'insert_val_into_string_int'
    character(len=:), allocatable :: test_name
    character(len=:), allocatable :: str

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    str = 'Test: {int_01}' // &
      ', {int_02}' // &
      ', john' // &
      ', {int_04}' // &
      ', {int_05}' // &
      ', a_very_long_and_nonsensical_marker_should_work_too&&*@#(#(#))' // &
      ', {int_07}' // &
      ', {int_08}' // &
      ', {int_09}'

    str = UPSY%stru%insert_val_into_string_int( str, '{int_01}', 1)
    str = UPSY%stru%insert_val_into_string_int( str, '{int_02}', 2)
    str = UPSY%stru%insert_val_into_string_int( str, 'john', 3)
    str = UPSY%stru%insert_val_into_string_int( str, '{int_04}', 4)
    str = UPSY%stru%insert_val_into_string_int( str, '{int_05}', 5)
    str = UPSY%stru%insert_val_into_string_int( str, 'a_very_long_and_nonsensical_marker_should_work_too&&*@#(#(#))', 6)
    str = UPSY%stru%insert_val_into_string_int( str, '{int_07}', 7)
    str = UPSY%stru%insert_val_into_string_int( str, '{int_08}', 8)
    str = UPSY%stru%insert_val_into_string_int( str, '{int_09}', 9)

    call unit_test( str == 'Test: 1, 2, 3, 4, 5, 6, 7, 8, 9', test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_insert_val_into_string_int

  subroutine test_insert_val_into_string_dp( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'test_insert_val_into_string_dp'
    character(len=*), parameter   :: test_name_local = 'insert_val_into_string_dp'
    character(len=:), allocatable :: test_name
    character(len=:), allocatable :: str

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    str = 'Test: {int_01}' // &
      ', {dp_02}' // &
      ', john' // &
      ', a_very_long_and_nonsensical_marker_should_work_too&&*@#(#(#))'

    str = UPSY%stru%insert_val_into_string_dp( str, '{int_01}', 1._dp)
    str = UPSY%stru%insert_val_into_string_dp( str, '{dp_02}', pi)
    str = UPSY%stru%insert_val_into_string_dp( str, 'john', 3e12_dp)
    str = UPSY%stru%insert_val_into_string_dp( str, 'a_very_long_and_nonsensical_marker_should_work_too&&*@#(#(#))', NaN)

    call unit_test( str == 'Test: 0.10000000000000E+01, 0.31415926535898E+01, 0.30000000000000E+13,                  NaN', test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_insert_val_into_string_dp

  subroutine test_capitalise_string( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'test_capitalise_string'
    character(len=*), parameter   :: test_name_local = 'capitalise_string'
    character(len=:), allocatable :: test_name
    character(len=:), allocatable :: str

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    str = 'This is a test STRING with mostly lowercase characters and also some numbers 1234567892#$%^&'
    str = UPSY%stru%capitalise_string( str)
    call unit_test( str == 'THIS IS A TEST STRING WITH MOSTLY LOWERCASE CHARACTERS AND ALSO SOME NUMBERS 1234567892#$%^&', test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_capitalise_string

  subroutine test_remove_leading_spaces( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'test_remove_leading_spaces'
    character(len=*), parameter   :: test_name_local = 'remove_leading_spaces'
    character(len=:), allocatable :: test_name
    character(len=:), allocatable :: str

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    str = '      bla'
    str = UPSY%stru%remove_leading_spaces( str)
    call unit_test( str == 'bla', test_name // '/string')

    str = '      '
    str = UPSY%stru%remove_leading_spaces( str)
    call unit_test( str == '', test_name // '/empty_string')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remove_leading_spaces

  subroutine test_str2int( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'test_str2int'
    character(len=*), parameter   :: test_name_local = 'str2int'
    character(len=:), allocatable :: test_name
    character(len=:), allocatable :: str
    integer                       :: int, stat

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    str = '1337'
    int = UPSY%stru%str2int( str, stat)
    call unit_test( int == 1337 .and. stat == 0, test_name // '/int')

    str = '0001337'
    int = UPSY%stru%str2int( str, stat)
    call unit_test( int == 1337 .and. stat == 0, test_name // '/int_leading_zeros')

    str = '   1337'
    int = UPSY%stru%str2int( str, stat)
    call unit_test( int == 1337 .and. stat == 0, test_name // '/int_leading_spaces')

    str = 'this shouldnt work'
    int = UPSY%stru%str2int( str, stat)
    call unit_test( stat /= 0, test_name // '/invalid')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_str2int

  subroutine test_int2str_with_leading_zeros( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'test_int2str_with_leading_zeros'
    character(len=*), parameter   :: test_name_local = 'int2str_with_leading_zeros'
    character(len=:), allocatable :: test_name
    character(len=:), allocatable :: str

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    str = UPSY%stru%int2str_with_leading_zeros( 1337, 5)
    call unit_test( str == '01337', test_name // '/int')

    str = UPSY%stru%int2str_with_leading_zeros( 1337, 3)
    call unit_test( str == '***', test_name // '/invalid')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_int2str_with_leading_zeros

  subroutine test_strrep( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'test_strrep'
    character(len=*), parameter   :: test_name_local = 'strrep'
    character(len=:), allocatable :: test_name
    character(len=:), allocatable :: str

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    str = 'Please replace John with Pete'
    str = UPSY%stru%strrep( str, 'John', 'Pete')
    call unit_test( str == 'Please replace Pete with Pete', test_name // '/one')

    str = 'Please replace John, John, John-Bartholomew and Johnnieboy with Pete'
    str = UPSY%stru%strrep( str, 'John', 'Pete')
    call unit_test( str == 'Please replace Pete, Pete, Pete-Bartholomew and Petenieboy with Pete', test_name // '/two')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_strrep

end module ut_string