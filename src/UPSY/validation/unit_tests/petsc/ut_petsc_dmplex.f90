module ut_petsc_dmplex

#include <petsc/finclude/petscdmplex.h>
  use petscdmplex
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

    call test_setup_simple_petsc_dmplex_unstructured_grid( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_petsc_dmplex_main

  subroutine test_setup_simple_petsc_dmplex_unstructured_grid( test_name_parent)

    ! See also: https://petsc.org/release/manual/dmplex/

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'test_setup_simple_petsc_dmplex_unstructured_grid'
    character(len=1024), parameter          :: test_name_local = 'setup_simple_petsc_dmplex_unstructured_grid'
    character(len=1024)                     :: test_name
    type(tDM)                               :: dm
    integer                                 :: perr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Create a DMPLEX object
    call DMPlexCreate( MPI_COMM_WORLD, dm, perr)

    ! First, we declare the set of points present in a mesh
    call DMPlexSetChart( dm, 0, 11, perr)

    ! Note that a chart here corresponds to a semi-closed interval
    ! (e.g [0,11)={0,1,…,10}[0,11)={0,1,…,10}) specifying the range of indices we’d like
    ! to use to define points on the current rank. We then define the covering relation,
    ! which we call the cone, which are also the in-edges in the DAG. In order to
    ! preallocate correctly, we first provide sizes,

    ! DMPlexSetConeSize(dm, point, number of points that cover the point)
    call DMPlexSetConeSize( dm,  0, 3, perr)
    call DMPlexSetConeSize( dm,  1, 3, perr)
    call DMPlexSetConeSize( dm,  6, 2, perr)
    call DMPlexSetConeSize( dm,  7, 2, perr)
    call DMPlexSetConeSize( dm,  8, 2, perr)
    call DMPlexSetConeSize( dm,  9, 2, perr)
    call DMPlexSetConeSize( dm, 10, 2, perr)
    call DMSetUp( dm, perr)

    ! and then point values (recall each point is an integer that represents a single
    ! geometric entity, a cell, face, edge, or vertex),

    ! DMPlexSetCone(dm, point, [points that cover the point])
    call DMPlexSetCone( dm,  0, [6, 7, 8], perr)
    call DMPlexSetCone( dm,  1, [7, 9, 10], perr)
    call DMPlexSetCone( dm,  6, [2, 3], perr)
    call DMPlexSetCone( dm,  7, [3, 4], perr)
    call DMPlexSetCone( dm,  8, [4, 2], perr)
    call DMPlexSetCone( dm,  9, [4, 5], perr)
    call DMPlexSetCone( dm, 10, [5, 3], perr)

    ! calculate the dual relations
    call DMPlexSymmetrize( dm, perr)

    ! The “symmetrization” is in the sense of the DAG. Each point knows its covering (cone)
    ! and each point knows what it covers (support). Note that when using automatic
    ! symmetrization, cones will be ordered but supports will not. The user can enforce
    ! an ordering on supports by rewriting them after symmetrization using DMPlexSetSupport().

    ! In order to support efficient queries, we construct fast search structures and indices
    ! for the different types of points using

    call DMPlexStratify( dm, perr)

    ! Yay, we made it this far!
    call unit_test( .true., test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_setup_simple_petsc_dmplex_unstructured_grid

end module ut_petsc_dmplex
