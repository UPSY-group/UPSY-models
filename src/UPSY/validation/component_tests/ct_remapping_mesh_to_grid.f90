module ct_remapping_mesh_to_grid

  ! Test everything related to remapping

  use UPSY_main, only: UPSY
  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_CHAR
  use precisions, only: dp
  use mpi_basic, only: par
  use parameters
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use netcdf_io_main
  use apply_maps, only: clear_all_maps_involving_this_mesh
  use remapping_main, only: map_from_mesh_vertices_to_xy_grid_2D, map_from_mesh_triangles_to_xy_grid_2D
  use ct_remapping_basic, only: calc_test_function_on_grid, calc_test_function_on_mesh, &
    calc_test_function_on_mesh_triangles

  implicit none

  private

  public :: run_all_mesh_to_grid_remapping_tests

contains

  !> Run all the mesh-to-grid remapping tests
  subroutine run_all_mesh_to_grid_remapping_tests( output_dir, test_mesh_filenames, test_grid_filenames)

    ! In/output variables:
    character(len=1024)           , intent(in) :: output_dir
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames
    character(len=*), dimension(:), intent(in) :: test_grid_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_all_mesh_to_grid_remapping_tests'
    integer                        :: i_mesh, i_grid
    character(len=1024)            :: filename_mesh, filename_grid

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '    Running mesh-to-grid remapping component tests...'
    if (par%primary) write(0,*) ''

    do i_mesh = 1, size( test_mesh_filenames)
      filename_mesh = test_mesh_filenames( i_mesh)
      do i_grid = 1, size( test_grid_filenames)
        filename_grid = test_grid_filenames( i_grid)
        call run_mesh_vertices_to_grid_remapping_tests_on_mesh_grid_combo ( output_dir, filename_mesh, filename_grid)
        call run_mesh_triangles_to_grid_remapping_tests_on_mesh_grid_combo( output_dir, filename_mesh, filename_grid)
      end do
    end do

    if (par%primary) write(0,*) ''

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_mesh_to_grid_remapping_tests

  !> Run all the mesh-to-grid remapping tests on one mesh/grid combination
  subroutine run_mesh_vertices_to_grid_remapping_tests_on_mesh_grid_combo( output_dir, filename_mesh, filename_grid)

    ! In/output variables:
    character(len=*), intent(in) :: output_dir
    character(len=*), intent(in) :: filename_mesh
    character(len=*), intent(in) :: filename_grid

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'run_mesh_to_grid_remapping_tests_on_mesh_grid_combo'
    character(len=1024)                   :: mesh_name, grid_name
    integer                               :: ncid
    type(type_grid)                       :: grid
    type(type_mesh)                       :: mesh
    real(dp), dimension(:), allocatable   :: d_grid_ex, d_mesh_ex, d_grid
    character(len=1024)                   :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    mesh_name = filename_mesh( index( filename_mesh, '/', back = .true.)+1 : len_trim( filename_mesh)-3)
    grid_name = filename_grid( index( filename_grid, '/', back = .true.)+1 : len_trim( filename_grid)-3)
    filename = trim( output_dir) // '/res_' // &
      trim( mesh_name) // '_vertices_TO_' // trim( grid_name) // '.nc'

    if (par%primary) write(0,*) '      Running mesh-vertices-to-grid remapping tests on mesh-grid combination:'
    if (par%primary) write(0,*) '        mesh: ', UPSY%stru%colour_string( trim( mesh_name),'light blue')
    if (par%primary) write(0,*) '        grid: ', UPSY%stru%colour_string( trim( grid_name),'light blue')

    ! Set up the mesh and the grid from the provided files
    call open_existing_netcdf_file_for_reading( filename_mesh, ncid)
    call setup_mesh_from_file( filename_mesh, ncid, mesh)
    call close_netcdf_file( ncid)

    call open_existing_netcdf_file_for_reading( filename_grid, ncid)
    call setup_xy_grid_from_file( filename_grid, ncid, grid)
    call close_netcdf_file( ncid)

    ! Calculate exact solution on the grid and the mesh
    call calc_test_function_on_grid( grid, d_grid_ex)
    call calc_test_function_on_mesh( mesh, d_mesh_ex)

    ! Map gridded data to the mesh
    allocate( d_grid( grid%n_loc))
    call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, output_dir, d_mesh_ex, d_grid)

    ! Write results to NetCDF
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call add_field_mesh_dp_2D_notime( filename, ncid, 'd_mesh_ex')
    call add_field_grid_dp_2D_notime( filename, ncid, 'd_grid_ex')
    call add_field_grid_dp_2D_notime( filename, ncid, 'd_grid')

    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd_mesh_ex', d_mesh_ex)
    call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'd_grid_ex', d_grid_ex)
    call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'd_grid'   , d_grid   )

    call close_netcdf_file( ncid)

    ! Clean up after yourself
    call clear_all_maps_involving_this_mesh( mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_mesh_vertices_to_grid_remapping_tests_on_mesh_grid_combo

  subroutine run_mesh_triangles_to_grid_remapping_tests_on_mesh_grid_combo( output_dir, filename_mesh, filename_grid)

    ! In/output variables:
    character(len=*), intent(in) :: output_dir
    character(len=*), intent(in) :: filename_mesh
    character(len=*), intent(in) :: filename_grid

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'run_mesh_triangles_to_grid_remapping_tests_on_mesh_grid_combo'
    character(len=1024)                   :: mesh_name, grid_name
    integer                               :: ncid
    type(type_grid)                       :: grid
    type(type_mesh)                       :: mesh
    real(dp), dimension(:), allocatable   :: d_grid_ex, d_mesh_ex, d_grid
    character(len=1024)                   :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    mesh_name = filename_mesh( index( filename_mesh, '/', back = .true.)+1 : len_trim( filename_mesh)-3)
    grid_name = filename_grid( index( filename_grid, '/', back = .true.)+1 : len_trim( filename_grid)-3)
    filename = trim( output_dir) // '/res_' // &
      trim( mesh_name) // '_triangles_TO_' // trim( grid_name) // '.nc'

    if (par%primary) write(0,*) '      Running mesh-triangles-to-grid remapping tests on mesh-grid combination:'
    if (par%primary) write(0,*) '        mesh: ', UPSY%stru%colour_string( trim( mesh_name),'light blue')
    if (par%primary) write(0,*) '        grid: ', UPSY%stru%colour_string( trim( grid_name),'light blue')

    ! Set up the mesh and the grid from the provided files
    call open_existing_netcdf_file_for_reading( filename_mesh, ncid)
    call setup_mesh_from_file( filename_mesh, ncid, mesh)
    call close_netcdf_file( ncid)

    call open_existing_netcdf_file_for_reading( filename_grid, ncid)
    call setup_xy_grid_from_file( filename_grid, ncid, grid)
    call close_netcdf_file( ncid)

    ! Calculate exact solution on the grid and the mesh
    call calc_test_function_on_grid( grid, d_grid_ex)
    call calc_test_function_on_mesh_triangles( mesh, d_mesh_ex)

    ! Map gridded data to the mesh
    allocate( d_grid( grid%n_loc))
    call map_from_mesh_triangles_to_xy_grid_2D( mesh, grid, output_dir, d_mesh_ex, d_grid)

    ! Write results to NetCDF
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call add_field_mesh_dp_2D_b_notime( filename, ncid, 'd_mesh_ex')
    call add_field_grid_dp_2D_notime( filename, ncid, 'd_grid_ex')
    call add_field_grid_dp_2D_notime( filename, ncid, 'd_grid')

    call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, filename, ncid, 'd_mesh_ex', d_mesh_ex)
    call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'd_grid_ex', d_grid_ex)
    call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'd_grid'   , d_grid   )

    call close_netcdf_file( ncid)

    ! Clean up after yourself
    call clear_all_maps_involving_this_mesh( mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_mesh_triangles_to_grid_remapping_tests_on_mesh_grid_combo

end module ct_remapping_mesh_to_grid
