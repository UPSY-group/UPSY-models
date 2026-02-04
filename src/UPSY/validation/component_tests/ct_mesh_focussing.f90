module ct_mesh_focussing

  ! Test everything related to mesh focussing

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_CHAR
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use assertions_basic, only: assert
  use tests_main, only: test_mesh_is_self_consistent
  use netcdf_io_main
  use polyline_types, only: type_polyline
  use parameters, only: pi
  use mesh_focussing, only: focus_mesh_on_polyline

  implicit none

  private

  public :: run_all_mesh_focussing_component_tests

contains

  !> Run all mesh focussing component tests.
  subroutine run_all_mesh_focussing_component_tests( output_dir, test_mesh_filenames)

    ! In/output variables:
    character(len=*),               intent(in) :: output_dir
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_all_mesh_focussing_component_tests'
    character(len=1024)            :: foldername_mesh_focussing
    integer                        :: i
    character(len=1024)            :: test_mesh_filename

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '  Running mesh_focussing component tests...'
    if (par%primary) write(0,*) ''

    call create_mesh_focussing_component_tests_output_folder( output_dir, foldername_mesh_focussing)

    do i = 1, size(test_mesh_filenames)
      test_mesh_filename = trim(test_mesh_filenames( i))
      call run_mesh_focussing_tests_on_mesh( foldername_mesh_focussing, test_mesh_filename)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_mesh_focussing_component_tests

  !> Create the output folder for the mesh_focussing component tests
  subroutine create_mesh_focussing_component_tests_output_folder( output_dir, foldername_mesh_focussing)

    ! In/output variables:
    character(len=*), intent(in   ) :: output_dir
    character(len=*), intent(  out) :: foldername_mesh_focussing

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_mesh_focussing_component_tests_output_folder'
    logical                        :: ex
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    foldername_mesh_focussing = trim( output_dir) // '/mesh_focussing'

    if (par%primary) then

      ! Remove existing folder if necessary
      inquire( file = trim( foldername_mesh_focussing) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( foldername_mesh_focussing))
      end if

      ! Create the directory
      call system('mkdir ' // trim( foldername_mesh_focussing))

    end if
    call MPI_BCAST( foldername_mesh_focussing, len(foldername_mesh_focussing), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_mesh_focussing_component_tests_output_folder

  !> Run all the mapping/derivative tests on a particular mesh
  subroutine run_mesh_focussing_tests_on_mesh( foldername_mesh_focussing, test_mesh_filename)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_mesh_focussing
    character(len=*), intent(in) :: test_mesh_filename

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'run_mesh_focussing_tests_on_mesh'
    type(type_mesh)                    :: mesh, mesh_focused
    integer                            :: ncid
    real(dp), parameter                :: rr_min = 0.25_dp
    real(dp), parameter                :: rr_max = 0.75_dp
    integer,  parameter                :: nr     = 25
    integer,  parameter                :: ntheta = 100
    integer                            :: ri, thetai
    real(dp)                           :: rr, r, theta, xmid, ymid, x, y
    type(type_polyline)                :: ll
    integer, dimension(:), allocatable :: vi_new2vi_old, vi_old2vi_new
    integer, dimension(:), allocatable :: vi_ll2vi, vi2vi_ll
    character(len=3)                   :: i_str
    character(len=12)                  :: r_str
    character(len=1024)                :: mesh_raw_name, filename
    integer                            :: pprev, prog

    ! Add routine to call stack
    call init_routine( routine_name)

    mesh_raw_name = test_mesh_filename( index( test_mesh_filename,'/', back = .true.)+1: &
      len_trim( test_mesh_filename))
    if (par%primary) write(0,*) '      Running mesh focussing tests on mesh ', &
      UPSY%stru%colour_string( trim( mesh_raw_name),'light blue'), '...'

    ! Set up the mesh from the file (includes calculating secondary geometry data and matrix operators)
    call open_existing_netcdf_file_for_reading( trim(test_mesh_filename), ncid)
    call setup_mesh_from_file( test_mesh_filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Check mesh self-consistency
    call assert( test_mesh_is_self_consistent( mesh), 'mesh is not self-consistent')

    ! Focus this mesh on a bunch of different polylines
    pprev = 0
    do ri = 1, nr

      ! Display progress
      prog = floor( 100._dp * real( ri-1,dp) / real( nr-1,dp))
      if (prog >= pprev+10) then
        pprev = prog
        if (par%primary) write(0,*) '       Progress: ', prog, ' %...'
      end if

      rr = rr_min + (rr_max - rr_min) * real( ri-1,dp) / real( nr-1,dp)
      r = (mesh%xmax - mesh%xmin) * 0.5_dp * rr

      ! Create polyline
      ll%is_closed = .true.
      ll%n = ntheta
      allocate( ll%p( ll%n,2))

      xmid = (mesh%xmin + mesh%xmax) / 2._dp
      ymid = (mesh%ymin + mesh%ymax) / 2._dp

      do thetai = 1, ntheta
        theta = 2._dp * pi * real( thetai,dp) / real( ntheta,dp)
        x = xmid + r * cos( theta)
        y = ymid + r * sin( theta)
        ll%p( thetai,:) = [x,y]
      end do

      ! Focus mesh on the polyline
      call focus_mesh_on_polyline( mesh, ll, mesh_focused, &
        vi_new2vi_old, vi_old2vi_new, vi_ll2vi, vi2vi_ll)

      call assert( test_mesh_is_self_consistent( mesh_focused), 'focused mesh is not self-consistent')

      ! Clean up after yourself
      deallocate( ll%p, vi_new2vi_old, vi_old2vi_new, vi_ll2vi, vi2vi_ll)

      ! Save focused mesh to NetCDF
      write( i_str,'(I3.3)') ri
      write( r_str,'(E12.6)') r
      filename = trim( foldername_mesh_focussing) // '/' // &
        mesh_raw_name( 1: len_trim( mesh_raw_name) - 3) // '_' // i_str // '_r' // r_str // '.nc'
      call save_mesh_as_netcdf( trim( filename), mesh_focused)

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_mesh_focussing_tests_on_mesh

end module ct_mesh_focussing
