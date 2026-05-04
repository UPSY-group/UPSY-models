module ocean_ismip

  use precisions, only: dp
  use UPSY_main, only: UPSY
  use parameters, only: NaN
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ocean_model_types, only: type_ocean_model, type_ocean_model_ismip
  use netcdf_io_main
  use mpi_f08, only: MPI_BCAST, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
  use basic_model_utilities, only: list_files_in_folder

  implicit none

  private

  public :: initialise_ocean_model_ismip, run_ocean_model_ismip

contains

  subroutine run_ocean_model_ismip( mesh, ocean, time)

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    type(type_ocean_model), intent(inout) :: ocean
    real(dp),               intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_ocean_model_ismip'
    real(dp)                       :: w0, w1

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_ocean_model_ismip

  subroutine initialise_ocean_model_ismip( mesh, ismip)

    ! In/output variables:
    type(type_mesh),              intent(in   ) :: mesh
    type(type_ocean_model_ismip), intent(inout) :: ismip

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_ocean_model_ismip'
    character(len=1024)            :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Get filenames
    call gather_filenames( 'thetao', ismip%filenames_T, ismip%time_bnds_T)
    call gather_filenames( 'so',     ismip%filenames_S, ismip%time_bnds_S)

    ! Allocate memory

    allocate (ismip%T_t0( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate (ismip%S_t0( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)

    allocate (ismip%T_t1( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate (ismip%S_t1( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model_ismip

  subroutine gather_filenames( subfolder, filenames, time_bnds)

    ! In/output variables:
    character(len=*),                            intent(in   ) :: subfolder
    character(len=*), dimension(:), allocatable, intent(inout) :: filenames
    real(dp), dimension(:,:), allocatable,       intent(inout) :: time_bnds

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'gather_filenames'
    character(len=1024)            :: foldername, filename
    integer                        :: i, ncid, id_dim_time, nt, id_var_time_bnds
    integer                        :: ierr
    real(dp), dimension(:,:), allocatable :: tbnds

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Get full folder name
    foldername = trim(C%foldername_ocean_ismip) // '/' // trim(subfolder)

    ! Get filenames
    call list_files_in_folder( foldername, filenames)

    ! Get associated time ranges
    allocate (time_bnds( size(filenames), 2))

    do i = 1, size(filenames)

      ! Construct the full filename
      filename = trim(foldername) // '/' // trim(filenames( i))

      ! Open, get time length, and time_bnds dimension
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)
      call inquire_var_multopt( filename, ncid, 'time_bnds', id_var_time_bnds)

      ! Allocate and read timebounds
      allocate( tbnds( 2, nt))
      call read_var_primary( filename, ncid, id_var_time_bnds, tbnds)
      call close_netcdf_file( ncid)

      ! Copy tbnds into time_bounds
      call MPI_BCAST( tbnds(:,:), nt * 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      time_bnds( i, 1) = tbnds( 1, 1)
      time_bnds( i, 2) = tbnds( 2, nt)
      deallocate( tbnds)

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine gather_filenames

end module ocean_ismip
