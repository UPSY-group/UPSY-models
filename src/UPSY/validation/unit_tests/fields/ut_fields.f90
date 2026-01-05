module ut_fields

  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use ut_fields_create_field, only: test_create_field

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

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_fields_main

end module ut_fields