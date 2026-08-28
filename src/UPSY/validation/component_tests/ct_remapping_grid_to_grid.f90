module ct_remapping_grid_to_grid

  ! Test everything related to remapping

  use UPSY_main, only: UPSY
  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_CHAR
  use precisions, only: dp
  use mpi_basic, only: par
  use parameters
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning
  use grid_types, only: type_grid
  use ct_remapping_basic, only: calc_test_function_on_grid
  use remapping_grid_to_grid, only: map_from_xy_grid_to_xy_grid_2D
  use netcdf_io_main
  use netcdf, only: NF90_DOUBLE, NF90_PUT_VAR

  implicit none

  private

  public :: run_all_grid_to_grid_remapping_tests

contains

  !> Run all the grid-to-grid remapping tests
  subroutine run_all_grid_to_grid_remapping_tests( output_dir, test_grid_filenames)

    ! In/output variables:
    character(len=1024)           , intent(in) :: output_dir
    character(len=*), dimension(:), intent(in) :: test_grid_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_all_grid_to_grid_remapping_tests'
    integer                        :: i1, i2
    character(len=1024)            :: filename_grid1, filename_grid2

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '    Running grid-to-grid remapping component tests...'
    if (par%primary) write(0,*) ''

    do i1 = 1, size( test_grid_filenames)
      filename_grid1 = test_grid_filenames( i1)
      do i2 = 1, size( test_grid_filenames)
        filename_grid2 = test_grid_filenames( i2)
        call run_grid_to_grid_remapping_tests_on_grid_grid_combo( output_dir, filename_grid1, filename_grid2)
      end do
    end do

    if (par%primary) write(0,*) ''

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_grid_to_grid_remapping_tests

  !> Run all the grid-to-grid remapping tests on one grid/grid combination
  subroutine run_grid_to_grid_remapping_tests_on_grid_grid_combo( output_dir, filename_grid1, filename_grid2)

    ! In/output variables:
    character(len=*), intent(in) :: output_dir
    character(len=*), intent(in) :: filename_grid1
    character(len=*), intent(in) :: filename_grid2

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'run_grid_to_grid_remapping_tests_on_grid_grid_combo'
    character(len=1024)                   :: grid_name1, grid_name2
    integer                               :: ncid
    type(type_grid)                       :: grid1, grid2
    real(dp), dimension(:), allocatable   :: d_grid1_ex, d_grid2_ex
    real(dp), dimension(:), allocatable   :: d_grid2_cons
    real(dp), dimension(:), allocatable   :: d_grid1_ex_tot, d_grid2_cons_tot
    character(len=1024)                   :: filename
    integer                               :: nerr

    ! Add routine to call stack
    call init_routine( routine_name)

    grid_name1 = filename_grid1( index( filename_grid1, '/', back = .true.)+1 : len_trim( filename_grid1)-3)
    grid_name2 = filename_grid2( index( filename_grid2, '/', back = .true.)+1 : len_trim( filename_grid2)-3)
    filename = trim( output_dir) // '/res_' // &
      grid_name1( 1:len_trim(grid_name1)) // '_TO_' // grid_name2( 1:len_trim(grid_name2)) // '.nc'

    if (par%primary) write(0,*) '      Running grid-to-grid remapping tests on grid-grid combination:'
    if (par%primary) write(0,*) '        from grid: ', UPSY%stru%colour_string( trim( grid_name1),'light blue')
    if (par%primary) write(0,*) '          to grid: ', UPSY%stru%colour_string( trim( grid_name2),'light blue')

    ! Set up the grid and the grid from the provided files
    call open_existing_netcdf_file_for_reading( filename_grid1, ncid)
    call setup_xy_grid_from_file( filename_grid1, ncid, grid1)
    call close_netcdf_file( ncid)

    call open_existing_netcdf_file_for_reading( filename_grid2, ncid)
    call setup_xy_grid_from_file( filename_grid2, ncid, grid2)
    call close_netcdf_file( ncid)

    ! Calculate exact solution on both grides
    call calc_test_function_on_grid( grid1, d_grid1_ex)
    call calc_test_function_on_grid( grid2, d_grid2_ex)

    ! Map data to the new grid
    allocate( d_grid2_cons  ( grid2%n_loc))

    call map_from_xy_grid_to_xy_grid_2D( grid1, grid2, d_grid1_ex, d_grid2_cons)

    ! Write results to NetCDF
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid2)
    call add_field_grid_dp_2D_notime( filename, ncid, 'd_grid2_ex')
    call add_field_grid_dp_2D_notime( filename, ncid, 'd_grid2_cons')

    call write_to_field_multopt_grid_dp_2D_notime( grid2, filename, ncid, 'd_grid2_ex', d_grid2_ex)
    call write_to_field_multopt_grid_dp_2D_notime( grid2, filename, ncid, 'd_grid2_cons', d_grid2_cons)

    call close_netcdf_file( ncid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_grid_to_grid_remapping_tests_on_grid_grid_combo

end module ct_remapping_grid_to_grid
