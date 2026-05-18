module ocean_ISMIP7

  use precisions, only: dp
  use UPSY_main, only: UPSY
  use parameters, only: NaN
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ocean_model_types, only: type_ocean_model, type_ocean_model_ISMIP7, type_ocean_field_ISMIP7
  use netcdf_io_main
  use mpi_f08, only: MPI_BCAST, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
  use basic_model_utilities, only: list_files_in_folder
  use calendar, only: convert_time_to_days

  implicit none

  private

  public :: initialise_ocean_model_ISMIP7, run_ocean_model_ISMIP7

contains

  subroutine run_ocean_model_ISMIP7( mesh, ocean, time)

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    type(type_ocean_model), intent(inout) :: ocean
    real(dp),               intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_ocean_model_ISMIP7'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Update timeframes if necessary
    call update_timeframes( mesh, ocean%ISMIP7%T, time)
    call update_timeframes( mesh, ocean%ISMIP7%S, time)

    ! Interpolate T and S between timeframes
    call interpolate_single_field( mesh, ocean%ISMIP7%T, ocean%T, time)
    call interpolate_single_field( mesh, ocean%ISMIP7%S, ocean%S, time)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_ocean_model_ISMIP7

  subroutine initialise_ocean_model_ISMIP7( mesh, ISMIP7)

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_ocean_model_ISMIP7), intent(inout) :: ISMIP7

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_ocean_model_ISMIP7'
    character(len=1024)            :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Define field names
    ISMIP7%T%name = 'thetao'
    ISMIP7%S%name = 'so'

    ! Get info from files
    call gather_fileinfo( ISMIP7%T)
    call gather_fileinfo( ISMIP7%S)

    ! Deallocate if necessary
    if (allocated( ISMIP7%T%val0 )) deallocate( ISMIP7%T%val0 )
    if (allocated( ISMIP7%S%val0 )) deallocate( ISMIP7%S%val0 )
    if (allocated( ISMIP7%T%val1 )) deallocate( ISMIP7%T%val1 )
    if (allocated( ISMIP7%S%val1 )) deallocate( ISMIP7%S%val1 )

    ! Allocate memory for timeframes
    allocate (ISMIP7%T%val0( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate (ISMIP7%S%val0( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate (ISMIP7%T%val1( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate (ISMIP7%S%val1( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)

    ! Update timeframes to the current model time
    call update_timeframes( mesh, ISMIP7%T, C%start_time_of_run)
    call update_timeframes( mesh, ISMIP7%S, C%start_time_of_run)     

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model_ISMIP7

  subroutine update_timeframes( mesh, field, time)

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_ocean_field_ISMIP7), intent(inout) :: field
    real(dp),                      intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'update_timeframes'
    real(dp)                       :: days
    integer                        :: ti0_old, ti1_old

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Get current bracket indices
    ti0_old = field%ti0
    ti1_old = field%ti1

    ! Convert model time (years) to days_since_1850
    call convert_time_to_days( time, days, calendar='noleap', allow_residual=.true.)

    ! Update the indices of time slices before and after current time
    call update_bracket_indices( field, days)

    ! Update timeframes if necessary
    if (field%ti0 /= ti0_old) then
      call update_single_timeframe( mesh, field, field%ti0, field%val0)
    end if

    if (field%ti1 /= ti1_old) then
      call update_single_timeframe( mesh, field, field%ti1, field%val1)
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine update_timeframes

  subroutine update_bracket_indices( field, days)

    ! In/output variables:
    type(type_ocean_field_ISMIP7), intent(inout) :: field
    real(dp),                      intent(in   ) :: days

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'update_bracket_indices'
    integer                        :: i, n

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Determine total numer of timeframes available for this field
    n = size(field%alltimes)

    if (days <= field%alltimes(1)) then
      ! Model time before first available time value, return first two indices
      field%ti0 = 1
      field%ti1 = 2

    elseif (days >= field%alltimes(n)) then
      ! Model time after last available time value, return last two indices
      field%ti0 = n-1
      field%ti1 = n

    else
      ! Model time within array, return bracketing indices
      do i = 1, n
        if (field%alltimes(i) <= days) then
          field%ti0 = i
          field%ti1 = i+1
        end if
      end do

    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine update_bracket_indices

  subroutine update_single_timeframe( mesh, field, ti, val)

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_ocean_field_ISMIP7), intent(in   ) :: field
    integer,                       intent(in   ) :: ti
    real(dp), dimension(:,:),      intent(inout) :: val

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'update_single_timeframe'
    character(len=1024)            :: filename
    integer                        :: fi
    integer                        :: ncid, id_dim_time, nt, id_var_time, ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Determine file index from time index
    fi = field%allfi(ti)

    ! Define the full filename of the file that contains the required timeframe
    filename = trim(field%foldername) // '/' // trim(field%filenames(fi))

    if (par%primary) then
      write(0,*) '   Reading ISMIP7 ocean forcing from file: ', &
        UPSY%stru%colour_string( trim( filename), 'light blue')
    end if

    ! Read ocean field from that timeframe
    call read_field_from_file_3D_ocean( filename, trim(field%name), mesh, C%output_dir, C%z_ocean, val, &
        time_to_read = field%alltimes( ti))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine update_single_timeframe

  subroutine interpolate_single_field( mesh, field, val_out, time)

    ! In/output variables:
    type(type_mesh),                                       intent(in   ) :: mesh
    type(type_ocean_field_ISMIP7),                         intent(in   ) :: field
    real(dp), dimension( mesh%vi1:mesh%vi2, 1:C%nz_ocean), intent(inout) :: val_out
    real(dp),                                              intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'interpolate_single_field'
    real(dp)                       :: w0, w1
    real(dp)                       :: days

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Convert model time (years) to days_since_1850
    call convert_time_to_days( time, days, calendar='noleap', allow_residual=.true.)

    ! Get weights
    if (days < field%alltimes( field%ti0)) then
      ! Model time before bracket times, put full weight on the first timeframe
      w0 = 1._dp
    elseif (days > field%alltimes( field%ti1)) then
      ! Model time after bracket times, put full weight on the last timeframe
      w0 = 0._dp
    else
      ! Model time between bracket times, determine interpolated weight
      w0 = (field%alltimes( field%ti1) - days) / (field%alltimes( field%ti1) - field%alltimes( field%ti0))
    end if

    w1 = 1._dp - w0

    ! Apply interpolation to values (T or S) of both timeframes, and write to main ocean%T or ocean%S
    val_out( mesh%vi1:mesh%vi2, :) = w0 * field%val0( mesh%vi1:mesh%vi2, :) + w1 * field%val1( mesh%vi1:mesh%vi2, :)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine interpolate_single_field

  subroutine gather_fileinfo( field)

    ! In/output variables:
    type(type_ocean_field_ISMIP7), intent(inout) :: field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'gather_fileinfo'
    character(len=1024)            :: filename
    integer                        :: i, ncid, id_dim_time, nt, id_var_time
    integer                        :: ierr
    real(dp), dimension(:), allocatable :: time_tmp
    real(dp), dimension(1000)      :: time_buff
    integer, dimension(1000)       :: fi_buff
    integer                        :: cnt, t

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Get full folder name
    field%foldername = trim(C%ocean_ISMIP7_forcing_foldername) // '/' // trim(field%name) // &
      '/' // trim(C%ocean_ISMIP7_forcing_version)

    ! Get all filenames in this folder, assuming this folder only contains files for this specific field (thetao or so)
    call list_files_in_folder( field%foldername, field%filenames, trim(field%name))

    ! Initialise counter for timeframes
    cnt = 1

    ! Loop over all files of this field
    do i = 1, size(field%filenames)

      ! Construct the full filename
      filename = trim(field%foldername) // '/' // trim(field%filenames( i))

      ! Open file and extract time variable
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)
      call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

      ! Copy time variable into a temporary array and close netcdf
      allocate( time_tmp( nt))
      call read_var_primary( filename, ncid, id_var_time, time_tmp)
      call close_netcdf_file( ncid)

      ! Copy temporary time array to all processes
      call MPI_BCAST( time_tmp(:), nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

      ! Loop over time array
      do t = 1, nt
        ! Add time value to buffer
        time_buff( cnt) = time_tmp( t)
        ! Add file index to buffer
        fi_buff( cnt) = i
        ! Increase count
        cnt = cnt + 1
      end do

      deallocate( time_tmp)

    end do

    ! Allocate arrays for times and file indices, combined for all available files in this folder
    if (allocated( field%alltimes )) deallocate( field%alltimes )
    if (allocated( field%allfi    )) deallocate( field%allfi    )

    allocate( field%alltimes( cnt-1))
    allocate( field%allfi   ( cnt-1))

    ! Copy all valid values from buffer, so the length of these arrays is equal to the actual amount of timeframes
    field%alltimes = time_buff( 1:cnt-1)
    field%allfi    = fi_buff( 1:cnt-1)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine gather_fileinfo

end module ocean_ISMIP7
