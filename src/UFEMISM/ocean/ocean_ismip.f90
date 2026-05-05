module ocean_ismip

  use precisions, only: dp
  use UPSY_main, only: UPSY
  use parameters, only: NaN
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ocean_model_types, only: type_ocean_model, type_ocean_model_ismip, type_ocean_field_ismip
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

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Update timeframes if necessary
    call update_timeframes( mesh, ocean%ismip%T, time)
    call update_timeframes( mesh, ocean%ismip%S, time)

    ! TODO remove
    call convert_time_to_days( time, days, calendar='noleap', allow_residual=.true.)
    if (par%primary) print *, time, days, ocean%ismip%T%alltimes(ocean%ismip%T%ti0)

    ! Interpolate

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

    ! Define fields
    ismip%T%name = 'thetao'
    ismip%S%name = 'so'

    ! Get info from files
    call gather_fileinfo( ismip%T)
    call gather_fileinfo( ismip%S)

    ! Allocate memory
    allocate (ismip%T%val0( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate (ismip%S%val0( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate (ismip%T%val1( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate (ismip%S%val1( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)

    ! Update timeframes
    call update_timeframes( mesh, ismip%T, C%start_time_of_run)
    call update_timeframes( mesh, ismip%S, C%start_time_of_run)     

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model_ismip

  subroutine update_timeframes( mesh, field, time)

    ! In/output variables:
    type(type_mesh),              intent(in   ) :: mesh
    type(type_ocean_field_ismip), intent(inout) :: field
    real(dp),                     intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'update_timeframes'
    real(dp)                       :: days

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Convert model time (years) to days_since_1850
    call convert_time_to_days( time, days, calendar='noleap', allow_residual=.true.)

    ! Update the indices of time slices before and after current time
    call update_bracket_indices( field, days)

    ! Read timeframes


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine update_timeframes

  subroutine update_bracket_indices( field, days)

    ! In/output variables:
    type(type_ocean_field_ismip), intent(inout) :: field
    real(dp),                     intent(in   ) :: days

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'update_bracket_indices'
    integer                        :: i, n

    ! Add routine to call stack
    call init_routine( routine_name)

    n = size(field%alltimes)

    ! Model time before first available time value, return first two indices
    if (days <= field%alltimes(1)) then
      field%ti0 = 1
      field%ti1 = 2
      return
    end if

    ! Model time after last available time value, return last two indices
    if (days >= field%alltimes(n)) then
      field%ti0 = n-1
      field%ti1 = n
      return
    end if

    ! Model time within array, return bracketing indices
    do i = 1, n
      if (field%alltimes(i) <= days) then
        field%ti0 = i
        field%ti1 = i+1
      end if
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine update_bracket_indices


  !subroutine interpolate_timeframes( mesh, varname, filename, days, ocean_field)

  !  ! In/output variables:
  !  type(type_mesh),              intent(in   ) :: mesh
  !  character(len=*),             intent(in   ) :: varname
  !  character(len=*),             intent(in   ) :: filename
  !  real(dp),                     intent(in   ) :: days
  !  real(dp), dimension(:,:),     intent(inout) :: ocean_field

  !  ! Local variables:
  !  character(len=1024), parameter        :: routine_name = 'interpolate_timeframes'
  !  real(dp), dimension(:,:), allocatable :: ocean_field_t0, ocean_field_t1
  !  integer                               :: i, ti0, ti1, w0, w1
  !  integer                               :: ncid, id_dim_time, nt, id_var_time, ierr
  !  real(dp), parameter                   :: tol = 1.e-8_dp
  !  real(dp), dimension(:), allocatable   :: days_from_file
  !  character(len=1024)                   :: foldername, filename_full

  !  ! Add routine to call stack
  !  call init_routine( routine_name)

  !  ! Get full folder name and filename
  !  foldername = trim(C%foldername_ocean_ismip) // '/' // trim(varname)
  !  filename_full = trim(foldername) // '/' // trim(filename)
  !  
  !  ! Read time variable from the file
  !  call open_existing_netcdf_file_for_reading( filename_full, ncid)
  !  call check_time( filename_full, ncid)
  !  call inquire_dim_multopt( filename_full, ncid, field_name_options_time, id_dim_time, dim_length = nt)
  !  call inquire_var_multopt( filename_full, ncid, field_name_options_time, id_var_time)
  !  allocate( days_from_file( nt))
  !  call read_var_primary( filename_full, ncid, id_var_time, days_from_file)
  !  call MPI_BCAST( days_from_file(:), nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
  !  call close_netcdf_file( ncid)

  !  ! Check whether any single timeframe matches model time
  !  ti0 = -1
  !  ti1 = -1
  !  do i = 1, nt
  !    if (abs(days-days_from_file( i)) < tol) then
  !      ti0 = i
  !    end if
  !  end do

  !  ! If needed, find two bounding timeframes
  !  if (ti0 < 0) then
  !    ! No exact match found, so extracting two to interpolate
  !    ti0 = 1
  !    ti1 = 2
  !    do while (days > days_from_file(ti0))
  !      ti0 = ti1
  !      ti1 = ti1 + 1
  !    end do
  !  end if
 
  !  ! Read timeframe(s)
  !  if (ti1 < 0) then
  !    ! Found an exact timeframe, so just read that
  !    call read_field_from_file_3D_ocean( filename_full, trim(varname), mesh, C%output_dir, C%z_ocean, ocean_field, &
  !      time_to_read = days_from_file(ti0))
  !  else
  !    ! Interpolate between two fields

  !    ! Read bounding time frames
  !    allocate( ocean_field_t0( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
  !    allocate( ocean_field_t1( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)

  !    call read_field_from_file_3D_ocean( filename_full, trim(varname), mesh, C%output_dir, C%z_ocean, ocean_field_t0, &
  !      time_to_read = days_from_file(ti0))
  !    call read_field_from_file_3D_ocean( filename_full, trim(varname), mesh, C%output_dir, C%z_ocean, ocean_field_t1, &
  !      time_to_read = days_from_file(ti1))

  !    ! Get weights
  !    w0 = (days_from_file( ti1) - days) / (days_from_file( ti1) - days_from_file( ti0))
  !    w1 = 1._dp - w0

  !    ocean_field = w0 * ocean_field_t0 + w1 * ocean_field_t1 

  !  end if

  !  ! Remove routine from call stack
  !  call finalise_routine( routine_name)

  !end subroutine interpolate_timeframes


  subroutine gather_fileinfo( field)

    ! In/output variables:
    type(type_ocean_field_ismip), intent(inout) :: field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'gather_fileinfo'
    character(len=1024)            :: foldername, filename
    integer                        :: i, ncid, id_dim_time, nt, id_var_time
    integer                        :: ierr
    real(dp), dimension(:), allocatable :: time_tmp
    real(dp), dimension(1000)      :: time_buff
    integer, dimension(1000)       :: fi_buff
    integer                        :: cnt, t

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Get full folder name
    foldername = trim(C%foldername_ocean_ismip) // '/' // trim(field%name)

    ! Get filenames
    call list_files_in_folder( foldername, field%filenames)

    cnt = 1

    do i = 1, size(field%filenames)

      ! Construct the full filename
      filename = trim(foldername) // '/' // trim(field%filenames( i))

      ! Open, get time length, and time_bnds dimension
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)
      call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

      ! Allocate and read timebounds
      allocate( time_tmp( nt))
      call read_var_primary( filename, ncid, id_var_time, time_tmp)
      call close_netcdf_file( ncid)

      ! Copy time_tmptbnds into time_bounds
      call MPI_BCAST( time_tmp(:), nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      do t = 1, nt
        time_buff( cnt) = time_tmp( t)
        fi_buff( cnt) = i
        cnt = cnt + 1
      end do

      deallocate( time_tmp)

    end do

    allocate( field%alltimes( cnt-1))
    allocate( field%allfi   ( cnt-1))

    ! Copy all valid values from buffer
    field%alltimes = time_buff( 1:cnt-1)
    field%allfi    = fi_buff( 1:cnt-1)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine gather_fileinfo

end module ocean_ismip
