module ct_ugrid

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use mpi_basic, only: par
  use mesh_basic, only: type_mesh
  use UPSY_main, only: UPSY
  use netcdf_io_main, only: open_existing_netcdf_file_for_reading,  setup_mesh_from_file, &
    close_netcdf_file
  use assertions_basic, only: assert
  use tests_main, only: test_mesh_is_self_consistent
  use ugrid_basic, only: type_ugrid
  use Arakawa_grid_mod, only: Arakawa_grid
  use netcdf_io_main, only: save_ugrid_as_netcdf

  implicit none

  private

  public :: run_all_ugrid_component_tests

contains

  subroutine run_all_ugrid_component_tests( output_dir, test_mesh_filenames)

    ! In/output variables:
    character(len=*),               intent(in) :: output_dir
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'run_all_ugrid_component_tests'
    integer                       :: i
    character(len=:), allocatable :: test_mesh_filename

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '  Running ugrid component tests...'
    if (par%primary) write(0,*) ''

    do i = 1, size( test_mesh_filenames)
      test_mesh_filename = trim( test_mesh_filenames( i))
      call run_ugrid_tests_on_mesh( output_dir, test_mesh_filename)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_ugrid_component_tests

  subroutine run_ugrid_tests_on_mesh( output_dir, test_mesh_filename)

    ! In/output variables:
    character(len=*), intent(in) :: output_dir
    character(len=*), intent(in) :: test_mesh_filename

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'run_ugrid_tests_on_mesh'
    type(type_mesh)               :: mesh
    character(len=:), allocatable :: mesh_raw_name
    integer                       :: ncid
    type(type_ugrid)              :: ugrid_a
    character(len=:), allocatable :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    mesh_raw_name = test_mesh_filename( index( test_mesh_filename,'/', back = .true.)+1: &
      len_trim( test_mesh_filename))
    if (par%primary) write(0,*) '      Running ugrid tests on mesh ', &
      UPSY%stru%colour_string( trim( mesh_raw_name),'light blue'), '...'

    ! Set up the mesh from the file (includes calculating secondary geometry data and matrix operators)
    call open_existing_netcdf_file_for_reading( trim( test_mesh_filename), ncid)
    call setup_mesh_from_file( test_mesh_filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Check mesh self-consistency
    call assert( test_mesh_is_self_consistent( mesh), 'mesh is not self-consistent')

    ! Convert mesh to ugrid and save to NetCDF
    call ugrid_a%build_from_grid( mesh, Arakawa_grid%a())
    filename = trim( output_dir) // '/ugrid_' // trim( mesh_raw_name) // '_a.nc'
    call save_ugrid_as_netcdf( filename, ugrid_a)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_ugrid_tests_on_mesh

end module ct_ugrid
