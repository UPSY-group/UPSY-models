module ut_petsc_dmplex

#include <petsc/finclude/petscksp.h>
  use petscksp
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use tests_main
  use assertions_basic
  use ut_basic

  implicit none

  private

  public :: unit_tests_petsc_dmplex_main

contains

  subroutine unit_tests_petsc_dmplex_main( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'unit_tests_petsc_dmplex_main'
    character(len=1024), parameter          :: test_name_local = 'dmplex'
    character(len=1024)                     :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_petsc_dmplex_main

end module ut_petsc_dmplex
