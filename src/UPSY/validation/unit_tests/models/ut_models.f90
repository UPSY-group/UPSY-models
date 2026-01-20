module ut_models

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use parameters, only: pi
  use mesh_memory, only: allocate_mesh_primary, crop_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_Lloyds_algorithm, only: Lloyds_algorithm_single_iteration
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_disc_calc_matrix_operators_2D, only: calc_all_matrix_operators_mesh
  use ut_demo_model, only: test_demo_model

  implicit none

  private

  public :: unit_tests_models_main

contains

  subroutine unit_tests_models_main( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_models_main'
    character(len=1024), parameter :: test_name_local = 'models'
    character(len=1024)            :: test_name
    real(dp)                       :: alpha_min, res_max
    real(dp), parameter            :: xmin = -1._dp
    real(dp), parameter            :: xmax =  1._dp
    real(dp), parameter            :: ymin = -1._dp
    real(dp), parameter            :: ymax =  1._dp
    type(type_mesh)                :: mesh1, mesh2

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Source mesh
    ! ===========

    call allocate_mesh_primary( mesh1, 'dummy_mesh_1', 100, 200)
    call initialise_dummy_mesh_5( mesh1, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 23.2_dp
    call refine_mesh_uniform( mesh1, res_max, alpha_min)
    call Lloyds_algorithm_single_iteration( mesh1, alpha_min)
    call crop_mesh_primary( mesh1)
    call calc_all_secondary_mesh_data( mesh1, 0._dp, -90._dp, 71._dp)
    call calc_all_matrix_operators_mesh( mesh1)

    ! Destination mesh
    ! ================

    call allocate_mesh_primary( mesh2, 'dummy_mesh_2', 100, 200)
    call initialise_dummy_mesh_5( mesh2, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 18._dp
    call refine_mesh_uniform( mesh2, res_max, alpha_min)
    call Lloyds_algorithm_single_iteration( mesh2, alpha_min)
    call crop_mesh_primary( mesh2)
    call calc_all_secondary_mesh_data( mesh2, 0._dp, -90._dp, 71._dp)
    call calc_all_matrix_operators_mesh( mesh2)

    ! Run all unit tests
    ! ==================

    call test_demo_model( test_name, mesh1, mesh2)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_models_main

end module ut_models