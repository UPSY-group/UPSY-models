module ut_fields_io

  use precisions, only: dp
  use mpi_basic, only: par
  use parameters, only: pi
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use fields_main, only: third_dimension, &
    atype_field, atype_field_2D, atype_field_3D, &
    type_field_logical_2D, type_field_int_2D, type_field_dp_2D, &
    type_field_logical_3D, type_field_int_3D, type_field_dp_3D, &
    type_fields_registry
  use ut_basic
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use grid_basic, only: setup_square_grid
  use Arakawa_grid_mod, only: Arakawa_grid
  use mpi_f08, only: MPI_WIN
  use mesh_memory, only: allocate_mesh_primary, crop_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use netcdf_io_main

  implicit none

  private

  public :: test_io

contains

  subroutine test_io( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_io'
    character(len=1024), parameter :: test_name_local = 'io'
    character(len=1024)            :: test_name
    type(type_grid)                :: grid
    real(dp)                       :: alpha_min, res_max
    real(dp), parameter            :: xmin = -1._dp
    real(dp), parameter            :: xmax =  1._dp
    real(dp), parameter            :: ymin = -1._dp
    real(dp), parameter            :: ymax =  1._dp
    type(type_mesh)                :: mesh

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call setup_square_grid( 'dummy_grid', 0._dp, 1._dp, 0._dp, 1._dp, 0.1_dp, grid)

    call allocate_mesh_primary( mesh, 'dummy_mesh', 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 23.2_dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call crop_mesh_primary( mesh)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)

    call test_io_grid_logical_2D      ( test_name, grid)
    call test_io_grid_logical_3D_zeta ( test_name, grid)
    call test_io_grid_logical_3D_month( test_name, grid)
    call test_io_grid_logical_3D_ocean( test_name, grid)

    call test_io_grid_int_2D      ( test_name, grid)
    call test_io_grid_int_3D_zeta ( test_name, grid)
    call test_io_grid_int_3D_month( test_name, grid)
    call test_io_grid_int_3D_ocean( test_name, grid)

    call test_io_grid_dp_2D      ( test_name, grid)
    call test_io_grid_dp_3D_zeta ( test_name, grid)
    call test_io_grid_dp_3D_month( test_name, grid)
    call test_io_grid_dp_3D_ocean( test_name, grid)



    call test_io_mesh_a_logical_2D      ( test_name, mesh)
    call test_io_mesh_a_logical_3D_zeta ( test_name, mesh)
    call test_io_mesh_a_logical_3D_month( test_name, mesh)
    call test_io_mesh_a_logical_3D_ocean( test_name, mesh)

    call test_io_mesh_a_int_2D      ( test_name, mesh)
    call test_io_mesh_a_int_3D_zeta ( test_name, mesh)
    call test_io_mesh_a_int_3D_month( test_name, mesh)
    call test_io_mesh_a_int_3D_ocean( test_name, mesh)

    call test_io_mesh_a_dp_2D      ( test_name, mesh)
    call test_io_mesh_a_dp_3D_zeta ( test_name, mesh)
    call test_io_mesh_a_dp_3D_month( test_name, mesh)
    call test_io_mesh_a_dp_3D_ocean( test_name, mesh)



    call test_io_mesh_b_logical_2D      ( test_name, mesh)
    call test_io_mesh_b_logical_3D_zeta ( test_name, mesh)
    call test_io_mesh_b_logical_3D_month( test_name, mesh)
    call test_io_mesh_b_logical_3D_ocean( test_name, mesh)

    call test_io_mesh_b_int_2D      ( test_name, mesh)
    call test_io_mesh_b_int_3D_zeta ( test_name, mesh)
    call test_io_mesh_b_int_3D_month( test_name, mesh)
    call test_io_mesh_b_int_3D_ocean( test_name, mesh)

    call test_io_mesh_b_dp_2D      ( test_name, mesh)
    call test_io_mesh_b_dp_3D_zeta ( test_name, mesh)
    call test_io_mesh_b_dp_3D_month( test_name, mesh)
    call test_io_mesh_b_dp_3D_ocean( test_name, mesh)



    call test_io_mesh_c_logical_2D      ( test_name, mesh)
    call test_io_mesh_c_logical_3D_zeta ( test_name, mesh)
    call test_io_mesh_c_logical_3D_month( test_name, mesh)
    call test_io_mesh_c_logical_3D_ocean( test_name, mesh)

    call test_io_mesh_c_int_2D      ( test_name, mesh)
    call test_io_mesh_c_int_3D_zeta ( test_name, mesh)
    call test_io_mesh_c_int_3D_month( test_name, mesh)
    call test_io_mesh_c_int_3D_ocean( test_name, mesh)

    call test_io_mesh_c_dp_2D      ( test_name, mesh)
    call test_io_mesh_c_dp_3D_zeta ( test_name, mesh)
    call test_io_mesh_c_dp_3D_month( test_name, mesh)
    call test_io_mesh_c_dp_3D_ocean( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io




  subroutine test_io_grid_logical_2D( test_name_parent, grid)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_grid), target,    intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_grid_logical_2D'
    character(len=1024), parameter               :: test_name_local = 'grid/logical_2D'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    logical, dimension(:  ), contiguous, pointer :: d1 => null()
    logical, dimension(:  ), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_grid_logical_2D'
    long_name = 'd_grid_logical_2D_long_name'
    units     = 'd_grid_logical_2D_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call flds_reg1%create_field( d1, w1, &
      grid, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( grid%n1) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      grid, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_grid_logical_2D

  subroutine test_io_grid_logical_3D_zeta( test_name_parent, grid)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_grid), target,    intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_grid_logical_3D_zeta'
    character(len=1024), parameter               :: test_name_local = 'grid/logical_3D_zeta'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 10
    logical, dimension(:,:), contiguous, pointer :: d1 => null()
    logical, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_grid_logical_3D_zeta'
    long_name = 'd_grid_logical_3D_zeta_long_name'
    units     = 'd_grid_logical_3D_zeta_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call flds_reg1%create_field( d1, w1, &
      grid, Arakawa_grid%a(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( grid%n1,3) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      grid, Arakawa_grid%a(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_grid_logical_3D_zeta

  subroutine test_io_grid_logical_3D_month( test_name_parent, grid)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_grid), target,    intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_grid_logical_3D_month'
    character(len=1024), parameter               :: test_name_local = 'grid/logical_3D_month'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    logical, dimension(:,:), contiguous, pointer :: d1 => null()
    logical, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_grid_logical_3D_month'
    long_name = 'd_grid_logical_3D_month_long_name'
    units     = 'd_grid_logical_3D_month_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call flds_reg1%create_field( d1, w1, &
      grid, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( grid%n1,3) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      grid, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_grid_logical_3D_month

  subroutine test_io_grid_logical_3D_ocean( test_name_parent, grid)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_grid), target,    intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_grid_logical_3D_ocean'
    character(len=1024), parameter               :: test_name_local = 'grid/logical_3D_ocean'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 20
    logical, dimension(:,:), contiguous, pointer :: d1 => null()
    logical, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_grid_logical_3D_ocean'
    long_name = 'd_grid_logical_3D_ocean_long_name'
    units     = 'd_grid_logical_3D_ocean_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call flds_reg1%create_field( d1, w1, &
      grid, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( grid%n1,3) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      grid, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_grid_logical_3D_ocean

  subroutine test_io_grid_int_2D( test_name_parent, grid)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_grid), target,    intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_grid_int_2D'
    character(len=1024), parameter               :: test_name_local = 'grid/int_2D'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer, dimension(:  ), contiguous, pointer :: d1 => null()
    integer, dimension(:  ), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_grid_int_2D'
    long_name = 'd_grid_int_2D_long_name'
    units     = 'd_grid_int_2D_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call flds_reg1%create_field( d1, w1, &
      grid, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( grid%n1) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      grid, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_grid_int_2D

  subroutine test_io_grid_int_3D_zeta( test_name_parent, grid)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_grid), target,    intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_grid_int_3D_zeta'
    character(len=1024), parameter               :: test_name_local = 'grid/int_3D_zeta'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 10
    integer, dimension(:,:), contiguous, pointer :: d1 => null()
    integer, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_grid_int_3D_zeta'
    long_name = 'd_grid_int_3D_zeta_long_name'
    units     = 'd_grid_int_3D_zeta_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call flds_reg1%create_field( d1, w1, &
      grid, Arakawa_grid%a(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( grid%n1,3) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      grid, Arakawa_grid%a(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_grid_int_3D_zeta

  subroutine test_io_grid_int_3D_month( test_name_parent, grid)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_grid), target,    intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_grid_int_3D_month'
    character(len=1024), parameter               :: test_name_local = 'grid/int_3D_month'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer, dimension(:,:), contiguous, pointer :: d1 => null()
    integer, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_grid_int_3D_month'
    long_name = 'd_grid_int_3D_month_long_name'
    units     = 'd_grid_int_3D_month_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call flds_reg1%create_field( d1, w1, &
      grid, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( grid%n1,3) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      grid, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_grid_int_3D_month

  subroutine test_io_grid_int_3D_ocean( test_name_parent, grid)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_grid), target,    intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_grid_int_3D_ocean'
    character(len=1024), parameter               :: test_name_local = 'grid/int_3D_ocean'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 20
    integer, dimension(:,:), contiguous, pointer :: d1 => null()
    integer, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_grid_int_3D_ocean'
    long_name = 'd_grid_int_3D_ocean_long_name'
    units     = 'd_grid_int_3D_ocean_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call flds_reg1%create_field( d1, w1, &
      grid, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( grid%n1,3) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      grid, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_grid_int_3D_ocean

  subroutine test_io_grid_dp_2D( test_name_parent, grid)

    ! In/output variables:
    character(len=*),        intent(in   ) :: test_name_parent
    type(type_grid), target, intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_io_grid_dp_2D'
    character(len=1024), parameter                :: test_name_local = 'grid/dp_2D'
    character(len=1024)                           :: test_name
    type(type_fields_registry)                    :: flds_reg1, flds_reg2
    character(:), allocatable                     :: name, long_name, units, filename
    real(dp), dimension(:  ), contiguous, pointer :: d1 => null()
    real(dp), dimension(:  ), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                 :: w1, w2
    integer                                       :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_grid_dp_2D'
    long_name = 'd_grid_dp_2D_long_name'
    units     = 'd_grid_dp_2D_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call flds_reg1%create_field( d1, w1, &
      grid, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( grid%n1) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      grid, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_grid_dp_2D

  subroutine test_io_grid_dp_3D_zeta( test_name_parent, grid)

    ! In/output variables:
    character(len=*),        intent(in   ) :: test_name_parent
    type(type_grid), target, intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_io_grid_dp_3D_zeta'
    character(len=1024), parameter                :: test_name_local = 'grid/dp_3D_zeta'
    character(len=1024)                           :: test_name
    type(type_fields_registry)                    :: flds_reg1, flds_reg2
    character(:), allocatable                     :: name, long_name, units, filename
    integer                                       :: nz = 10
    real(dp), dimension(:,:), contiguous, pointer :: d1 => null()
    real(dp), dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                 :: w1, w2
    integer                                       :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_grid_dp_3D_zeta'
    long_name = 'd_grid_dp_3D_zeta_long_name'
    units     = 'd_grid_dp_3D_zeta_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call flds_reg1%create_field( d1, w1, &
      grid, Arakawa_grid%a(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( grid%n1,3) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      grid, Arakawa_grid%a(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_grid_dp_3D_zeta

  subroutine test_io_grid_dp_3D_month( test_name_parent, grid)

    ! In/output variables:
    character(len=*),        intent(in   ) :: test_name_parent
    type(type_grid), target, intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_io_grid_dp_3D_month'
    character(len=1024), parameter                :: test_name_local = 'grid/dp_3D_month'
    character(len=1024)                           :: test_name
    type(type_fields_registry)                    :: flds_reg1, flds_reg2
    character(:), allocatable                     :: name, long_name, units, filename
    real(dp), dimension(:,:), contiguous, pointer :: d1 => null()
    real(dp), dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                 :: w1, w2
    integer                                       :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_grid_dp_3D_month'
    long_name = 'd_grid_dp_3D_month_long_name'
    units     = 'd_grid_dp_3D_month_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call flds_reg1%create_field( d1, w1, &
      grid, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( grid%n1,3) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      grid, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_grid_dp_3D_month

  subroutine test_io_grid_dp_3D_ocean( test_name_parent, grid)

    ! In/output variables:
    character(len=*),        intent(in   ) :: test_name_parent
    type(type_grid), target, intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_io_grid_dp_3D_ocean'
    character(len=1024), parameter                :: test_name_local = 'grid/dp_3D_ocean'
    character(len=1024)                           :: test_name
    type(type_fields_registry)                    :: flds_reg1, flds_reg2
    character(:), allocatable                     :: name, long_name, units, filename
    integer                                       :: nz = 20
    real(dp), dimension(:,:), contiguous, pointer :: d1 => null()
    real(dp), dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                 :: w1, w2
    integer                                       :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_grid_dp_3D_ocean'
    long_name = 'd_grid_dp_3D_ocean_long_name'
    units     = 'd_grid_dp_3D_ocean_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_xy_grid_in_netcdf_file( filename, ncid, grid)

    call flds_reg1%create_field( d1, w1, &
      grid, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( grid%n1,3) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      grid, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_grid_dp_3D_ocean




  subroutine test_io_mesh_a_logical_2D( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_a_logical_2D'
    character(len=1024), parameter               :: test_name_local = 'mesh/a/logical_2D'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    logical, dimension(:  ), contiguous, pointer :: d1 => null()
    logical, dimension(:  ), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_a_logical_2D'
    long_name = 'd_mesh_a_logical_2D_long_name'
    units     = 'd_mesh_a_logical_2D_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_a_logical_2D

  subroutine test_io_mesh_a_logical_3D_zeta( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_a_logical_3D_zeta'
    character(len=1024), parameter               :: test_name_local = 'mesh/a/logical_3D_zeta'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 10
    logical, dimension(:,:), contiguous, pointer :: d1 => null()
    logical, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_a_logical_3D_zeta'
    long_name = 'd_mesh_a_logical_3D_zeta_long_name'
    units     = 'd_mesh_a_logical_3D_zeta_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%a(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%a(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_a_logical_3D_zeta

  subroutine test_io_mesh_a_logical_3D_month( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_a_logical_3D_month'
    character(len=1024), parameter               :: test_name_local = 'mesh/a/logical_3D_month'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    logical, dimension(:,:), contiguous, pointer :: d1 => null()
    logical, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_a_logical_3D_month'
    long_name = 'd_mesh_a_logical_3D_month_long_name'
    units     = 'd_mesh_a_logical_3D_month_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_a_logical_3D_month

  subroutine test_io_mesh_a_logical_3D_ocean( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_a_logical_3D_ocean'
    character(len=1024), parameter               :: test_name_local = 'mesh/a/logical_3D_ocean'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 20
    logical, dimension(:,:), contiguous, pointer :: d1 => null()
    logical, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_a_logical_3D_ocean'
    long_name = 'd_mesh_a_logical_3D_ocean_long_name'
    units     = 'd_mesh_a_logical_3D_ocean_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_a_logical_3D_ocean

  subroutine test_io_mesh_a_int_2D( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_a_int_2D'
    character(len=1024), parameter               :: test_name_local = 'mesh/a/int_2D'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer, dimension(:  ), contiguous, pointer :: d1 => null()
    integer, dimension(:  ), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_a_int_2D'
    long_name = 'd_mesh_a_int_2D_long_name'
    units     = 'd_mesh_a_int_2D_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_a_int_2D

  subroutine test_io_mesh_a_int_3D_zeta( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_a_int_3D_zeta'
    character(len=1024), parameter               :: test_name_local = 'mesh/a/int_3D_zeta'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 10
    integer, dimension(:,:), contiguous, pointer :: d1 => null()
    integer, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_a_int_3D_zeta'
    long_name = 'd_mesh_a_int_3D_zeta_long_name'
    units     = 'd_mesh_a_int_3D_zeta_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%a(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%a(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_a_int_3D_zeta

  subroutine test_io_mesh_a_int_3D_month( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_a_int_3D_month'
    character(len=1024), parameter               :: test_name_local = 'mesh/a/int_3D_month'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer, dimension(:,:), contiguous, pointer :: d1 => null()
    integer, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_a_int_3D_month'
    long_name = 'd_mesh_a_int_3D_month_long_name'
    units     = 'd_mesh_a_int_3D_month_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_a_int_3D_month

  subroutine test_io_mesh_a_int_3D_ocean( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_a_int_3D_ocean'
    character(len=1024), parameter               :: test_name_local = 'mesh/a/int_3D_ocean'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 20
    integer, dimension(:,:), contiguous, pointer :: d1 => null()
    integer, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_a_int_3D_ocean'
    long_name = 'd_mesh_a_int_3D_ocean_long_name'
    units     = 'd_mesh_a_int_3D_ocean_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_a_int_3D_ocean

  subroutine test_io_mesh_a_dp_2D( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_a_dp_2D'
    character(len=1024), parameter               :: test_name_local = 'mesh/a/dp_2D'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    real(dp), dimension(:  ), contiguous, pointer :: d1 => null()
    real(dp), dimension(:  ), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_a_dp_2D'
    long_name = 'd_mesh_a_dp_2D_long_name'
    units     = 'd_mesh_a_dp_2D_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_a_dp_2D

  subroutine test_io_mesh_a_dp_3D_zeta( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_a_dp_3D_zeta'
    character(len=1024), parameter               :: test_name_local = 'mesh/a/dp_3D_zeta'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 10
    real(dp), dimension(:,:), contiguous, pointer :: d1 => null()
    real(dp), dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_a_dp_3D_zeta'
    long_name = 'd_mesh_a_dp_3D_zeta_long_name'
    units     = 'd_mesh_a_dp_3D_zeta_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%a(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%a(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_a_dp_3D_zeta

  subroutine test_io_mesh_a_dp_3D_month( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_a_dp_3D_month'
    character(len=1024), parameter               :: test_name_local = 'mesh/a/dp_3D_month'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    real(dp), dimension(:,:), contiguous, pointer :: d1 => null()
    real(dp), dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_a_dp_3D_month'
    long_name = 'd_mesh_a_dp_3D_month_long_name'
    units     = 'd_mesh_a_dp_3D_month_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_a_dp_3D_month

  subroutine test_io_mesh_a_dp_3D_ocean( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_a_dp_3D_ocean'
    character(len=1024), parameter               :: test_name_local = 'mesh/a/dp_3D_ocean'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 20
    real(dp), dimension(:,:), contiguous, pointer :: d1 => null()
    real(dp), dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_a_dp_3D_ocean'
    long_name = 'd_mesh_a_dp_3D_ocean_long_name'
    units     = 'd_mesh_a_dp_3D_ocean_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_a_dp_3D_ocean




  subroutine test_io_mesh_b_logical_2D( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_b_logical_2D'
    character(len=1024), parameter               :: test_name_local = 'mesh/b/logical_2D'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    logical, dimension(:  ), contiguous, pointer :: d1 => null()
    logical, dimension(:  ), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_b_logical_2D'
    long_name = 'd_mesh_b_logical_2D_long_name'
    units     = 'd_mesh_b_logical_2D_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%b(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%b(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_b_logical_2D

  subroutine test_io_mesh_b_logical_3D_zeta( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_b_logical_3D_zeta'
    character(len=1024), parameter               :: test_name_local = 'mesh/b/logical_3D_zeta'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 10
    logical, dimension(:,:), contiguous, pointer :: d1 => null()
    logical, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_b_logical_3D_zeta'
    long_name = 'd_mesh_b_logical_3D_zeta_long_name'
    units     = 'd_mesh_b_logical_3D_zeta_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_b_logical_3D_zeta

  subroutine test_io_mesh_b_logical_3D_month( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_b_logical_3D_month'
    character(len=1024), parameter               :: test_name_local = 'mesh/b/logical_3D_month'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    logical, dimension(:,:), contiguous, pointer :: d1 => null()
    logical, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_b_logical_3D_month'
    long_name = 'd_mesh_b_logical_3D_month_long_name'
    units     = 'd_mesh_b_logical_3D_month_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%b(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%b(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_b_logical_3D_month

  subroutine test_io_mesh_b_logical_3D_ocean( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_b_logical_3D_ocean'
    character(len=1024), parameter               :: test_name_local = 'mesh/b/logical_3D_ocean'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 20
    logical, dimension(:,:), contiguous, pointer :: d1 => null()
    logical, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_b_logical_3D_ocean'
    long_name = 'd_mesh_b_logical_3D_ocean_long_name'
    units     = 'd_mesh_b_logical_3D_ocean_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%b(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%b(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_b_logical_3D_ocean

  subroutine test_io_mesh_b_int_2D( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_b_int_2D'
    character(len=1024), parameter               :: test_name_local = 'mesh/b/int_2D'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer, dimension(:  ), contiguous, pointer :: d1 => null()
    integer, dimension(:  ), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_b_int_2D'
    long_name = 'd_mesh_b_int_2D_long_name'
    units     = 'd_mesh_b_int_2D_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%b(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%b(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_b_int_2D

  subroutine test_io_mesh_b_int_3D_zeta( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_b_int_3D_zeta'
    character(len=1024), parameter               :: test_name_local = 'mesh/b/int_3D_zeta'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 10
    integer, dimension(:,:), contiguous, pointer :: d1 => null()
    integer, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_b_int_3D_zeta'
    long_name = 'd_mesh_b_int_3D_zeta_long_name'
    units     = 'd_mesh_b_int_3D_zeta_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_b_int_3D_zeta

  subroutine test_io_mesh_b_int_3D_month( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_b_int_3D_month'
    character(len=1024), parameter               :: test_name_local = 'mesh/b/int_3D_month'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer, dimension(:,:), contiguous, pointer :: d1 => null()
    integer, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_b_int_3D_month'
    long_name = 'd_mesh_b_int_3D_month_long_name'
    units     = 'd_mesh_b_int_3D_month_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%b(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%b(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_b_int_3D_month

  subroutine test_io_mesh_b_int_3D_ocean( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_b_int_3D_ocean'
    character(len=1024), parameter               :: test_name_local = 'mesh/b/int_3D_ocean'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 20
    integer, dimension(:,:), contiguous, pointer :: d1 => null()
    integer, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_b_int_3D_ocean'
    long_name = 'd_mesh_b_int_3D_ocean_long_name'
    units     = 'd_mesh_b_int_3D_ocean_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%b(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%b(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_b_int_3D_ocean

  subroutine test_io_mesh_b_dp_2D( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_b_dp_2D'
    character(len=1024), parameter               :: test_name_local = 'mesh/b/dp_2D'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    real(dp), dimension(:  ), contiguous, pointer :: d1 => null()
    real(dp), dimension(:  ), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_b_dp_2D'
    long_name = 'd_mesh_b_dp_2D_long_name'
    units     = 'd_mesh_b_dp_2D_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%b(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%b(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_b_dp_2D

  subroutine test_io_mesh_b_dp_3D_zeta( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_b_dp_3D_zeta'
    character(len=1024), parameter               :: test_name_local = 'mesh/b/dp_3D_zeta'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 10
    real(dp), dimension(:,:), contiguous, pointer :: d1 => null()
    real(dp), dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_b_dp_3D_zeta'
    long_name = 'd_mesh_b_dp_3D_zeta_long_name'
    units     = 'd_mesh_b_dp_3D_zeta_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_b_dp_3D_zeta

  subroutine test_io_mesh_b_dp_3D_month( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_b_dp_3D_month'
    character(len=1024), parameter               :: test_name_local = 'mesh/b/dp_3D_month'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    real(dp), dimension(:,:), contiguous, pointer :: d1 => null()
    real(dp), dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_b_dp_3D_month'
    long_name = 'd_mesh_b_dp_3D_month_long_name'
    units     = 'd_mesh_b_dp_3D_month_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%b(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%b(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_b_dp_3D_month

  subroutine test_io_mesh_b_dp_3D_ocean( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_b_dp_3D_ocean'
    character(len=1024), parameter               :: test_name_local = 'mesh/b/dp_3D_ocean'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 20
    real(dp), dimension(:,:), contiguous, pointer :: d1 => null()
    real(dp), dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_b_dp_3D_ocean'
    long_name = 'd_mesh_b_dp_3D_ocean_long_name'
    units     = 'd_mesh_b_dp_3D_ocean_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%b(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%b(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_b_dp_3D_ocean




  subroutine test_io_mesh_c_logical_2D( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_c_logical_2D'
    character(len=1024), parameter               :: test_name_local = 'mesh/c/logical_2D'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    logical, dimension(:  ), contiguous, pointer :: d1 => null()
    logical, dimension(:  ), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_c_logical_2D'
    long_name = 'd_mesh_c_logical_2D_long_name'
    units     = 'd_mesh_c_logical_2D_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%c(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%c(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_c_logical_2D

  subroutine test_io_mesh_c_logical_3D_zeta( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_c_logical_3D_zeta'
    character(len=1024), parameter               :: test_name_local = 'mesh/c/logical_3D_zeta'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 10
    logical, dimension(:,:), contiguous, pointer :: d1 => null()
    logical, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_c_logical_3D_zeta'
    long_name = 'd_mesh_c_logical_3D_zeta_long_name'
    units     = 'd_mesh_c_logical_3D_zeta_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%c(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%c(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_c_logical_3D_zeta

  subroutine test_io_mesh_c_logical_3D_month( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_c_logical_3D_month'
    character(len=1024), parameter               :: test_name_local = 'mesh/c/logical_3D_month'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    logical, dimension(:,:), contiguous, pointer :: d1 => null()
    logical, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_c_logical_3D_month'
    long_name = 'd_mesh_c_logical_3D_month_long_name'
    units     = 'd_mesh_c_logical_3D_month_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%c(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%c(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_c_logical_3D_month

  subroutine test_io_mesh_c_logical_3D_ocean( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_c_logical_3D_ocean'
    character(len=1024), parameter               :: test_name_local = 'mesh/c/logical_3D_ocean'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 20
    logical, dimension(:,:), contiguous, pointer :: d1 => null()
    logical, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_c_logical_3D_ocean'
    long_name = 'd_mesh_c_logical_3D_ocean_long_name'
    units     = 'd_mesh_c_logical_3D_ocean_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%c(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = .true.

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%c(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_c_logical_3D_ocean

  subroutine test_io_mesh_c_int_2D( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_c_int_2D'
    character(len=1024), parameter               :: test_name_local = 'mesh/c/int_2D'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer, dimension(:  ), contiguous, pointer :: d1 => null()
    integer, dimension(:  ), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_c_int_2D'
    long_name = 'd_mesh_c_int_2D_long_name'
    units     = 'd_mesh_c_int_2D_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%c(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%c(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_c_int_2D

  subroutine test_io_mesh_c_int_3D_zeta( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_c_int_3D_zeta'
    character(len=1024), parameter               :: test_name_local = 'mesh/c/int_3D_zeta'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 10
    integer, dimension(:,:), contiguous, pointer :: d1 => null()
    integer, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_c_int_3D_zeta'
    long_name = 'd_mesh_c_int_3D_zeta_long_name'
    units     = 'd_mesh_c_int_3D_zeta_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%c(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%c(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_c_int_3D_zeta

  subroutine test_io_mesh_c_int_3D_month( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_c_int_3D_month'
    character(len=1024), parameter               :: test_name_local = 'mesh/c/int_3D_month'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer, dimension(:,:), contiguous, pointer :: d1 => null()
    integer, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_c_int_3D_month'
    long_name = 'd_mesh_c_int_3D_month_long_name'
    units     = 'd_mesh_c_int_3D_month_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%c(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%c(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_c_int_3D_month

  subroutine test_io_mesh_c_int_3D_ocean( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_c_int_3D_ocean'
    character(len=1024), parameter               :: test_name_local = 'mesh/c/int_3D_ocean'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 20
    integer, dimension(:,:), contiguous, pointer :: d1 => null()
    integer, dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_c_int_3D_ocean'
    long_name = 'd_mesh_c_int_3D_ocean_long_name'
    units     = 'd_mesh_c_int_3D_ocean_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%c(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 1337

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%c(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_c_int_3D_ocean

  subroutine test_io_mesh_c_dp_2D( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_c_dp_2D'
    character(len=1024), parameter               :: test_name_local = 'mesh/c/dp_2D'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    real(dp), dimension(:  ), contiguous, pointer :: d1 => null()
    real(dp), dimension(:  ), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_c_dp_2D'
    long_name = 'd_mesh_c_dp_2D_long_name'
    units     = 'd_mesh_c_dp_2D_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%c(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%c(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_c_dp_2D

  subroutine test_io_mesh_c_dp_3D_zeta( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_c_dp_3D_zeta'
    character(len=1024), parameter               :: test_name_local = 'mesh/c/dp_3D_zeta'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 10
    real(dp), dimension(:,:), contiguous, pointer :: d1 => null()
    real(dp), dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_c_dp_3D_zeta'
    long_name = 'd_mesh_c_dp_3D_zeta_long_name'
    units     = 'd_mesh_c_dp_3D_zeta_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%c(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%c(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_c_dp_3D_zeta

  subroutine test_io_mesh_c_dp_3D_month( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_c_dp_3D_month'
    character(len=1024), parameter               :: test_name_local = 'mesh/c/dp_3D_month'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    real(dp), dimension(:,:), contiguous, pointer :: d1 => null()
    real(dp), dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_c_dp_3D_month'
    long_name = 'd_mesh_c_dp_3D_month_long_name'
    units     = 'd_mesh_c_dp_3D_month_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%c(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%c(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_c_dp_3D_month

  subroutine test_io_mesh_c_dp_3D_ocean( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_io_mesh_c_dp_3D_ocean'
    character(len=1024), parameter               :: test_name_local = 'mesh/c/dp_3D_ocean'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg1, flds_reg2
    character(:), allocatable                    :: name, long_name, units, filename
    integer                                      :: nz = 20
    real(dp), dimension(:,:), contiguous, pointer :: d1 => null()
    real(dp), dimension(:,:), contiguous, pointer :: d2 => null()
    type(MPI_WIN)                                :: w1, w2
    integer                                      :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd_mesh_c_dp_3D_ocean'
    long_name = 'd_mesh_c_dp_3D_ocean_long_name'
    units     = 'd_mesh_c_dp_3D_ocean_units'

    filename  = trim( foldername_unit_tests_output) // '/' // trim( name) // '.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)

    call flds_reg1%create_field( d1, w1, &
      mesh, Arakawa_grid%c(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    d1( mesh%vi1,3) = 13.37_dp

    call flds_reg1%items(1)%p%write_to_netcdf ( filename, ncid)
    call close_netcdf_file( ncid)

    call flds_reg2%create_field( d2, w2, &
      mesh, Arakawa_grid%c(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = 'reallocate')

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call flds_reg2%items(1)%p%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    call unit_test( flds_reg1%items(1)%p == flds_reg2%items(1)%p, test_name)

    ! Clean up after yourself
    call flds_reg1%destroy
    call flds_reg2%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_io_mesh_c_dp_3D_ocean

end module ut_fields_io