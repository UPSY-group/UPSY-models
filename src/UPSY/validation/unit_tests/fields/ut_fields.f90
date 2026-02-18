module ut_fields

  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning
  use ut_fields_create_field, only: test_create_field
  use ut_fields_io, only: test_io
  use ut_fields_reallocate, only: test_reallocate_field
  use ut_fields_remap, only: test_remap_field

  implicit none

  private

  public :: unit_tests_fields_main

contains

  subroutine unit_tests_fields_main( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_fields_main'
    character(len=1024), parameter :: test_name_local = 'fields'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Run all unit tests
    call test_create_field( test_name)
    call test_io( test_name)
    call test_reallocate_field( test_name)
    call test_remap_field( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_fields_main

end module ut_fields