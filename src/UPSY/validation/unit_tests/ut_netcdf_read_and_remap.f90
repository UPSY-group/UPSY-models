module ut_netcdf_read_and_remap

  use tests_main
  use assertions_basic
  use ut_basic
  use netcdf_io_main
  use mpi_basic, only: par
  use precisions, only: dp
  use parameters, only: pi
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use string_module, only: strrep
  use grid_types, only: type_grid
  use grid_basic, only: setup_square_grid
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_disc_calc_matrix_operators_2D, only: calc_all_matrix_operators_mesh
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, &
    MPI_MIN, MPI_MAX, MPI_COMM_WORLD
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use deallocate_dist_shared_mod, only: deallocate_dist_shared

  implicit none

  private

  public :: unit_tests_netcdf_read_and_remap

contains

  subroutine unit_tests_netcdf_read_and_remap( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'unit_tests_netcdf_read_and_remap'
    character(len=1024), parameter                 :: test_name_local = 'read_and_remap'
    character(len=1024)                            :: test_name
    character(len=1024)                            :: name
    real(dp)                                       :: xmin, xmax, ymin, ymax, alpha_min, res_max
    type(type_mesh)                                :: mesh
    character(len=1024), dimension(:), allocatable :: list_of_filenames
    integer                                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Create a simple test mesh
    name = 'test_mesh'
    xmin = -500e3_dp
    xmax =  500e3_dp
    ymin = -500e3_dp
    ymax =  500e3_dp
    alpha_min = 25._dp * pi / 180._dp
    res_max = 50e3_dp

    call allocate_mesh_primary( mesh, name, 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)
    call calc_all_matrix_operators_mesh( mesh)

    call create_dummy_input_file_xy_grid( mesh, list_of_filenames)
    call create_dummy_input_file_mesh   ( mesh, list_of_filenames)

    do i = 1, size( list_of_filenames,1)
      call test_read_from_file( test_name, mesh, list_of_filenames(i))
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_netcdf_read_and_remap

  subroutine test_read_from_file( test_name_parent, mesh, filename)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    character(len=*), intent(in) :: filename

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_read_from_file'
    character(len=1024)                         :: test_name
    real(dp), dimension(:), allocatable         :: d_dist
    real(dp), dimension(:), contiguous, pointer :: d_hybrid     => null()
    real(dp), dimension(:), contiguous, pointer :: d_hybrid_loc => null()
    type(MPI_WIN)                               :: wd_hybrid
    real(dp)                                    :: d_max, d_min
    integer                                     :: ierr
    logical                                     :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // &
      trim( filename( index( filename,'/',back=.true.)+11 : len_trim( filename)-3))

    ! Distributed data
    ! ================

    allocate( d_dist( mesh%vi1:mesh%vi2))

    call read_field_from_file_2D( filename, 'd', mesh, trim( foldername_unit_tests_output), d_dist)

    d_min = minval( d_dist)
    d_max = maxval( d_dist)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    test_result = &
      test_ge_le( d_min, -1.00_dp, -0.95_dp) .and. &
      test_ge_le( d_max,  0.95_dp,  1.00_dp)
    call unit_test( test_result, trim( test_name) // '_dist')

    ! Hybrid distributed/shared data
    ! ==============================

    call allocate_dist_shared( d_hybrid, wd_hybrid, mesh%pai_V%n_nih)

    call read_field_from_file_2D( filename, 'd', mesh, trim( foldername_unit_tests_output), d_hybrid)

    d_hybrid_loc( mesh%pai_V%i1: mesh%pai_V%i2) => d_hybrid( mesh%pai_V%i1: mesh%pai_V%i2)
    d_min = minval( d_hybrid_loc)
    d_max = maxval( d_hybrid_loc)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    test_result = &
      test_ge_le( d_min, -1.00_dp, -0.95_dp) .and. &
      test_ge_le( d_max,  0.95_dp,  1.00_dp)
    call unit_test( test_result, trim( test_name) // '_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_hybrid, wd_hybrid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_read_from_file

  subroutine create_dummy_input_file_xy_grid( mesh_mod, list_of_filenames)

    ! In/output variables:
    type(type_mesh),                                intent(in   ) :: mesh_mod
    character(len=1024), dimension(:), allocatable, intent(inout) :: list_of_filenames

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'create_dummy_input_file_xy_grid'
    character(len=1024)                 :: name
    real(dp)                            :: xmin, xmax, ymin, ymax, dx, x, y
    type(type_grid)                     :: grid
    real(dp), dimension(:), allocatable :: d
    integer                             :: i,j,n
    character(len=1024)                 :: filename
    integer                             :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Create test grid
    name = 'test_grid'
    xmin = mesh_mod%xmin
    xmax = mesh_mod%xmax
    ymin = mesh_mod%ymin
    ymax = mesh_mod%ymax
    dx = (xmax - xmin) / 50._dp
    call setup_square_grid( name, xmin, xmax, ymin, ymax, dx, grid)

    ! Create test data
    allocate( d( grid%n1: grid%n2))

    do n = grid%n1, grid%n2
      i = grid%n2ij( n,1)
      j = grid%n2ij( n,2)
      x = grid%x(i)
      y = grid%y(j)
      d( n) = test_function( x,y,xmin,xmax,ymin,ymax)
    end do

    filename = trim( foldername_unit_tests_output) // '/test_file_xy_grid.nc'
    call add_filename_to_list_of_filenames( list_of_filenames, filename)

    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    call add_field_grid_dp_2D_notime( filename, ncid, 'd')
    call write_to_field_multopt_grid_dp_2D_notime( grid, filename, ncid, 'd', d)
    call close_netcdf_file( ncid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_dummy_input_file_xy_grid

  subroutine create_dummy_input_file_mesh( mesh_mod, list_of_filenames)

    ! In/output variables:
    type(type_mesh),                                intent(in   ) :: mesh_mod
    character(len=1024), dimension(:), allocatable, intent(inout) :: list_of_filenames

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'create_dummy_input_file_mesh'
    character(len=1024)                 :: name
    real(dp)                            :: xmin, xmax, ymin, ymax, alpha_min, res_max
    type(type_mesh)                     :: mesh
    real(dp), dimension(:), allocatable :: d
    integer                             :: vi
    real(dp)                            :: x,y
    character(len=1024)                 :: filename
    integer                             :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Create a simple test mesh
    name = 'test_mesh'
    xmin = mesh_mod%xmin
    xmax = mesh_mod%xmax
    ymin = mesh_mod%ymin
    ymax = mesh_mod%ymax
    alpha_min = 25._dp * pi / 180._dp
    res_max = 60e3_dp

    call allocate_mesh_primary( mesh, name, 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)

    ! Create test data
    allocate( d( mesh%vi1:mesh%vi2))
    do vi = mesh%vi1, mesh%vi2
      x = mesh%V( vi,1)
      y = mesh%V( vi,2)
      d( vi) = test_function( x,y,xmin,xmax,ymin,ymax)
    end do

    filename = trim( foldername_unit_tests_output) // '/test_file_mesh.nc'
    call add_filename_to_list_of_filenames( list_of_filenames, filename)

    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call add_field_mesh_dp_2D_notime( filename, ncid, 'd')
    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'd', d)
    call close_netcdf_file( ncid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_dummy_input_file_mesh

  function test_function( x,y,xmin,xmax,ymin,ymax) result( res)
    real(dp), intent(in) :: x,y,xmin,xmax,ymin,ymax
    real(dp)             :: res
    res = &
      cos( 2.35_dp * pi * (x-xmin) / (xmax - xmin)) * &
      cos( 3.65_dp * pi * (y-ymin) / (ymax - ymin))
  end function test_function

  !> Add a filename to an allocatable list of filenames
  subroutine add_filename_to_list_of_filenames( list_of_filenames, filename)

    ! In/output variables:
    character(len=1024), dimension(:), allocatable, intent(inout) :: list_of_filenames
    character(len=1024),                            intent(in)    :: filename

    ! Local variables:
    character(len=1024), dimension(:), allocatable :: list_of_filenames_temp

    if (.not. allocated( list_of_filenames)) then
      allocate( list_of_filenames( 1))
      list_of_filenames( 1) = trim( adjustl( filename))
    else
      list_of_filenames_temp = list_of_filenames
      deallocate( list_of_filenames)
      allocate( list_of_filenames( size( list_of_filenames_temp)+1))
      list_of_filenames( 1:size( list_of_filenames_temp)  ) = list_of_filenames_temp
      list_of_filenames(   size( list_of_filenames_temp)+1) = trim( adjustl( filename))
    end if

  end subroutine add_filename_to_list_of_filenames

end module ut_netcdf_read_and_remap