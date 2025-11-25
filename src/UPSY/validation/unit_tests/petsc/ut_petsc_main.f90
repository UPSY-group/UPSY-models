module ut_petsc_main

  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use ut_petsc_matrix_vector_conversion, only: test_multiply_PETSc_matrix_with_vector_1D, test_matrix_PETSc_CSR_conversion
  use ut_petsc_dmplex, only: unit_tests_petsc_dmplex_main

  implicit none

  private

  public :: unit_tests_petsc_main

contains

  subroutine unit_tests_petsc_main( test_name_parent)
    ! Run all unit tests for the PETSc subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_petsc_main'
    character(len=1024), parameter :: test_name_local = 'petsc'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Run all unit tests for the PETSc subroutines
    call test_multiply_PETSc_matrix_with_vector_1D( test_name)
    call test_matrix_PETSc_CSR_conversion( test_name)
    call unit_tests_petsc_dmplex_main( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_petsc_main

end module ut_petsc_main
