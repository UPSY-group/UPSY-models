module ut_demo_model

  use precisions, only: dp
  use parameters, only: pi
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use models_demo, only: type_demo_model
  use ut_basic, only: unit_test, foldername_unit_tests_output

  implicit none

  private

  public :: test_demo_model

contains

  subroutine test_demo_model( test_name_parent, mesh1, mesh2)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh1, mesh2

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_demo_model'
    character(len=1024), parameter :: test_name_local = 'demo_model'
    character(len=1024)            :: test_name
    type(type_demo_model)          :: demo_model1, demo_model2, demo_model3
    integer, parameter             :: nz = 15
    character(:), allocatable      :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate and initialise an instance of the demo model
    call demo_model1%allocate  ( demo_model1%allocate_ct( mesh1))
    call demo_model1%initialise( demo_model1%initialise_ct( 1, 2))

    ! Run the demo model
    call demo_model1%run( demo_model1%run_ct( 3, 4))

    ! Write the current state of the demo model to a restart file
    call demo_model1%write_to_restart_file( foldername_unit_tests_output, filename)

    ! Create, but do not initialise, a second instance of the demo model
    call demo_model2%allocate( demo_model2%allocate_ct( mesh1))

    ! Initialise this instance's data from the restart file
    call demo_model2%read_from_restart_file( filename)

    ! Test if the two instances are identical
    call unit_test( demo_model1 == demo_model2, trim( test_name) // '_restart')

    ! Remap the first instance to another mesh, and re-initialise its data
    ! (so that the test won't fail due to remapping inaccuracies)
    call demo_model1%remap( demo_model1%remap_ct( mesh2))
    call demo_model1%initialise( demo_model1%initialise_ct( 1, 2))

    ! Create and initialise a third instance of the demo model directly on the new mesh
    call demo_model3%allocate  ( demo_model3%allocate_ct( mesh1))
    call demo_model3%initialise( demo_model3%initialise_ct( 1, 2))

    ! Test if the two instances are identical
    call unit_test( demo_model1 == demo_model3, trim( test_name) // '_remap')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_demo_model

end module ut_demo_model