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
  use calendar, only: convert_time_to_days

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
    real(dp)                       :: days
    character(len=1024)            :: filename_to_read_T, filename_to_read_S

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Convert model time (years) to days_since_1850
    call convert_time_to_days( time, days, calendar='noleap', allow_residual=.true.)

    ! Find filename which contains model time
    call find_filename_to_read( ocean%ismip%filenames_T, ocean%ismip%time_bnds_T, days, filename_to_read_T)
    call find_filename_to_read( ocean%ismip%filenames_S, ocean%ismip%time_bnds_S, days, filename_to_read_S)

    ! Interpolate timeframes
    call interpolate_timeframes( mesh, 'thetao', filename_to_read_T, days, ocean%T)
    call interpolate_timeframes( mesh, 'so'    , filename_to_read_S, days, ocean%S)

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

  subroutine find_filename_to_read( filenames, time_bnds, days, filename_to_read)

    ! In/output variables:
    character(len=*), dimension(:),      intent(in   ) :: filenames
    real(dp), dimension(:, :),           intent(in   ) :: time_bnds
    real(dp),                            intent(in   ) :: days
    character(len=1024),                 intent(  out) :: filename_to_read

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'find_filename_to_read'
    integer                        :: fidx

    ! Add routine to call stack
    call init_routine( routine_name)

    if (days < time_bnds( 1, 1)) then
      ! Model time before forcing time period, use first filename
      fidx = 1
    elseif (days > time_bnds( size(filenames), 2)) then
      ! Model time beyond forcing time period, use last filename
      fidx = size(filenames)
    else
      fidx = 1
      do while (days > time_bnds( fidx, 2))
        fidx = fidx + 1
      end do
    end if

    filename_to_read = filenames( fidx)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine find_filename_to_read

  subroutine interpolate_timeframes( mesh, varname, filename, days, ocean_field)

    ! In/output variables:
    type(type_mesh),              intent(in   ) :: mesh
    character(len=*),             intent(in   ) :: varname
    character(len=*),             intent(in   ) :: filename
    real(dp),                     intent(in   ) :: days
    real(dp), dimension(:,:),     intent(inout) :: ocean_field

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'interpolate_timeframes'
    real(dp), dimension(:,:), allocatable :: ocean_field_t0, ocean_field_t1
    integer                               :: i, ti0, ti1, w0, w1
    integer                               :: ncid, id_dim_time, nt, id_var_time, ierr
    real(dp), parameter                   :: tol = 1.e-8_dp
    real(dp), dimension(:), allocatable   :: days_from_file
    character(len=1024)                   :: foldername, filename_full

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Get full folder name and filename
    foldername = trim(C%foldername_ocean_ismip) // '/' // trim(varname)
    filename_full = trim(foldername) // '/' // trim(filename)
    
    ! Read time variable from the file
    call open_existing_netcdf_file_for_reading( filename_full, ncid)
    call check_time( filename_full, ncid)
    call inquire_dim_multopt( filename_full, ncid, field_name_options_time, id_dim_time, dim_length = nt)
    call inquire_var_multopt( filename_full, ncid, field_name_options_time, id_var_time)
    allocate( days_from_file( nt))
    call read_var_primary( filename_full, ncid, id_var_time, days_from_file)
    call MPI_BCAST( days_from_file(:), nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call close_netcdf_file( ncid)

    ! Check whether any single timeframe matches model time
    ti0 = -1
    ti1 = -1
    do i = 1, nt
      if (abs(days-days_from_file( i)) < tol) then
        ti0 = i
      end if
    end do

    ! If needed, find two bounding timeframes
    if (ti0 < 0) then
      ! No exact match found, so extracting two to interpolate
      ti0 = 1
      ti1 = 2
      do while (days > days_from_file(ti0))
        ti0 = ti1
        ti1 = ti1 + 1
      end do
    end if
 
    ! Read timeframe(s)
    if (ti1 < 0) then
      ! Found an exact timeframe, so just read that
      call read_field_from_file_3D_ocean( filename_full, trim(varname), mesh, C%output_dir, C%z_ocean, ocean_field, &
        time_to_read = days_from_file(ti0))
    else
      ! Interpolate between two fields

      ! Read bounding time frames
      allocate( ocean_field_t0( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
      allocate( ocean_field_t1( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)

      call read_field_from_file_3D_ocean( filename_full, trim(varname), mesh, C%output_dir, C%z_ocean, ocean_field_t0, &
        time_to_read = days_from_file(ti0))
      call read_field_from_file_3D_ocean( filename_full, trim(varname), mesh, C%output_dir, C%z_ocean, ocean_field_t1, &
        time_to_read = days_from_file(ti1))

      ! Get weights
      w0 = (days_from_file( ti1) - days) / (days_from_file( ti1) - days_from_file( ti0))
      w1 = 1._dp - w0

      ocean_field = w0 * ocean_field_t0 + w1 * ocean_field_t1 

    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine interpolate_timeframes


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
