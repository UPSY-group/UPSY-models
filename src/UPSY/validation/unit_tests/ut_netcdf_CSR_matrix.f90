module ut_netcdf_CSR_matrix

  ! Unit tests for the netcdf i/o - CSR matrix routines

  use tests_main
  use ut_basic
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use precisions, only: dp
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use ut_mpi_CSR_matrix_vector_multiplication, only: initialise_simple_matrix_equation_1, &
    initialise_simple_matrix_equation_2
  use model_configuration, only: C

  implicit none

  private

  public :: test_netcdf_CSR_matrix

contains

  subroutine test_netcdf_CSR_matrix( test_name_parent)
    ! Test the netcdf i/o subroutines for x/y-gridded data

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'test_netcdf_CSR_matrix'
    character(len=*), parameter   :: test_name_local = 'CSR_matrix_NetCDF'
    character(len=:), allocatable :: test_name
    type(type_CSR_matrix_dp)      :: A1, A2
    real(dp), dimension(1)        :: x1, x2, y1, y2

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call initialise_simple_matrix_equation_1( A1, x1, y1)
    call initialise_simple_matrix_equation_2( A2, x2, y2)

    call test_netcdf_CSR_single_matrix( test_name, A1, '1')
    call test_netcdf_CSR_single_matrix( test_name, A2, '2')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_netcdf_CSR_matrix

  subroutine test_netcdf_CSR_single_matrix( test_name_parent, A, test_number)
    ! Test the netcdf i/o subroutines for x/y-gridded data

    ! In/output variables:
    character(len=*),         intent(in) :: test_name_parent
    type(type_CSR_matrix_dp), intent(in) :: A
    character(len=*),         intent(in) :: test_number

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'test_netcdf_CSR_single_matrix'
    character(len=*), parameter   :: test_name_local = 'CSR_matrix_NetCDF'
    character(len=:), allocatable :: test_name
    character(len=:), allocatable :: filename_base
    type(type_CSR_matrix_dp)      :: A_from_file

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    filename_base = 'A_' // test_number // '.nc'
    call A%write_to_dist_NetCDFs( C%output_dir, filename_base)
    call A_from_file%read_from_dist_NetCDFs( C%output_dir, filename_base)

    call unit_test( test_eq( A%m          , A_from_file%m          ), trim( test_name) // '_' // test_number // '_m')
    call unit_test( test_eq( A%n          , A_from_file%n          ), trim( test_name) // '_' // test_number // '_n')
    call unit_test( test_eq( A%nnz        , A_from_file%nnz        ), trim( test_name) // '_' // test_number // '_nnz')

    call unit_test( test_eq( A%m_loc      , A_from_file%m_loc      ), trim( test_name) // '_' // test_number // '_m_loc')
    call unit_test( test_eq( A%i1         , A_from_file%i1         ), trim( test_name) // '_' // test_number // '_i1')
    call unit_test( test_eq( A%i2         , A_from_file%i2         ), trim( test_name) // '_' // test_number // '_i2')

    call unit_test( test_eq( A%m_node     , A_from_file%m_node     ), trim( test_name) // '_' // test_number // '_m_node')
    call unit_test( test_eq( A%i1_node    , A_from_file%i1_node    ), trim( test_name) // '_' // test_number // '_i1_node')
    call unit_test( test_eq( A%i2_node    , A_from_file%i2_node    ), trim( test_name) // '_' // test_number // '_i2_node')

    call unit_test( test_eq( A%n_loc      , A_from_file%n_loc      ), trim( test_name) // '_' // test_number // '_n_loc')
    call unit_test( test_eq( A%j1         , A_from_file%j1         ), trim( test_name) // '_' // test_number // '_j1')
    call unit_test( test_eq( A%j2         , A_from_file%j2         ), trim( test_name) // '_' // test_number // '_j2')

    call unit_test( test_eq( A%n_node     , A_from_file%n_node     ), trim( test_name) // '_' // test_number // '_n_node')
    call unit_test( test_eq( A%j1_node    , A_from_file%j1_node    ), trim( test_name) // '_' // test_number // '_j1_node')
    call unit_test( test_eq( A%j2_node    , A_from_file%j2_node    ), trim( test_name) // '_' // test_number // '_j2_node')

    call unit_test( test_eq( A%j_min_node , A_from_file%j_min_node ), trim( test_name) // '_' // test_number // '_j_min_node')
    call unit_test( test_eq( A%j_max_node , A_from_file%j_max_node ), trim( test_name) // '_' // test_number // '_j_max_node')
    call unit_test( test_eq( A%needs_x_tot, A_from_file%needs_x_tot), trim( test_name) // '_' // test_number // '_needs_x_tot')

    call unit_test( A%pai_x == A_from_file%pai_x, trim( test_name) // '_' // test_number // '_pai_x')
    call unit_test( A%pai_y == A_from_file%pai_y, trim( test_name) // '_' // test_number // '_pai_y')

    call unit_test( test_eq( A%ptr(A%i1:A%i2), A_from_file%ptr( A%i1:A%i2)), trim( test_name) // '_' // test_number // '_ptr')
    call unit_test( test_eq( A%ind           , A_from_file%ind            ), trim( test_name) // '_' // test_number // '_ind')
    call unit_test( test_eq( A%val           , A_from_file%val            ), trim( test_name) // '_' // test_number // '_val')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_netcdf_CSR_single_matrix

end module ut_netcdf_CSR_matrix