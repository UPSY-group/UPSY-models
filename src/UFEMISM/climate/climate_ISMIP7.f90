module climate_ISMIP7

  ! The ISMIP7 protocol provides separate NetCDF files for each year, with a single file contaning 12 monthly fields
  ! for that year. There are separate files for each variable, namely:
  !
  ! - tas        : the absolute 2m air temp (in K)
  ! - tas-anomaly: the 2m air temp anomaly relative to the 1960-1989 baseline
  ! - dtsdz      : the vertical surface temp gradient, to correct for differences in elevation between
  !                the observed and modelled ice-sheet geometry
  ! - pr         : the absolute precipitation (in kg m^-2 s^-1)
  ! - pr-anomaly : the precipitation anomaly relative to the 1960-1989 baseline
  !
  ! The directory structure for the forcing data is:
  !
  !  ../
  !    path/
  !      to/
  !        base_folder/
  !          tas/
  !            version/
  !              tas_somethingsomethingsomething_2015.nc
  !              tas_somethingsomethingsomething_2016.nc
  !              tas_somethingsomethingsomething_2017.nc
  !              ...
  !          tas-anomaly/
  !            version/
  !              tas-anomaly_somethingsomethingsomething_2015.nc
  !              tas-anomaly_somethingsomethingsomething_2016.nc
  !              tas-anomaly_somethingsomethingsomething_2017.nc
  !              ...
  !          dtsdz/
  !            version/
  !              dtsdz_somethingsomethingsomething_2015.nc
  !              dtsdz_somethingsomethingsomething_2016.nc
  !              dtsdz_somethingsomethingsomething_2017.nc
  !              ...
  !          pr/
  !            version/
  !              pr_somethingsomethingsomething_2015.nc
  !              pr_somethingsomethingsomething_2016.nc
  !              pr_somethingsomethingsomething_2017.nc
  !              ...
  !          pr-anomaly/
  !            version/
  !              pr-anomaly_somethingsomethingsomething_2015.nc
  !              pr-anomaly_somethingsomethingsomething_2016.nc
  !              pr-anomaly_somethingsomethingsomething_2017.nc
  !              ...
  !
  ! In the config, you only need to provide the path/to/base_folder and the version; UFEMISM will take
  ! care of the rest, assuming the directory structure is as expected.

  use mpi_basic, only: par
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use climate_model_types, only: type_climate_model, type_climate_model_ISMIP7, &
    type_climate_field_ISMIP7_monthly, type_climate_field_ISMIP7_yearly
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use mpi_f08, only: MPI_WIN
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use netcdf_io_main, only: read_field_from_file_2D_monthly, read_field_from_file_2D, read_time_from_file
  use basic_model_utilities, only: list_files_in_folder
  use parameters, only: NaN, sec_per_year, freshwater_density

  implicit none

  private

  public :: initialise_climate_model_ISMIP7, run_climate_model_ISMIP7

  interface initialise_climate_field
    module procedure initialise_climate_field_monthly
    module procedure initialise_climate_field_yearly
  end interface

  interface update_timeframes
    module procedure update_timeframes_monthly
    module procedure update_timeframes_yearly
  end interface

  interface interpolate_single_field
    module procedure interpolate_single_field_monthly
    module procedure interpolate_single_field_yearly
  end interface

contains

  subroutine run_climate_model_ISMIP7( mesh, ice, climate, time)

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    type(type_ice_model),         intent(in   ) :: ice
    type(type_climate_model), intent(inout) :: climate
    real(dp),                 intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'run_climate_model_ISMIP7'
    integer                                 :: vi, mi
    real(dp)                                :: delta_z
    real(dp), dimension( mesh%vi1:mesh%vi2) :: delta_ts

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Calculate elevation-based T2m correction
    call update_timeframes( mesh, climate%ISMIP7%dtsdz, time)
    call interpolate_single_field( mesh, climate%ISMIP7%dtsdz, time)

    do vi = mesh%vi1, mesh%vi2
      delta_z = ice%Hs( vi) - climate%ISMIP7%Hs_baseline ( vi)
      delta_ts( vi) = delta_z * climate%ISMIP7%dtsdz%val_interp( vi)
    end do

    ! Calculate monthly climate
    select case (C%climate_ISMIP7_choice_baseline)
    case default
      call crash('invalid climate_ISMIP7_choice_baseline "' // trim( C%climate_ISMIP7_choice_baseline) // '"')
    case ('yearly')
      ! Update timeframes
      call update_timeframes( mesh, climate%ISMIP7%tas, time)
      call update_timeframes( mesh, climate%ISMIP7%pr, time)

      ! Interpolate between timeframes
      call interpolate_single_field( mesh, climate%ISMIP7%tas, time)
      call interpolate_single_field( mesh, climate%ISMIP7%pr, time)

      ! Calculate monthly climate
      do vi = mesh%vi1, mesh%vi2
        do mi = 1, 12
          climate%T2m( vi, mi) = climate%ISMIP7%tas%val_interp( vi, mi) + delta_ts( vi)
          climate%Precip( vi, mi) = climate%ISMIP7%pr%val_interp( vi, mi) * sec_per_year / freshwater_density ! [m.w.e. yr^-1] 
        end do
      end do

    case ('fixed')
      ! Update timeframes
      call update_timeframes( mesh, climate%ISMIP7%tas_anomaly, time)
      call update_timeframes( mesh, climate%ISMIP7%pr_anomaly, time)

      ! Interpolate between timeframes
      call interpolate_single_field( mesh, climate%ISMIP7%tas_anomaly, time)
      call interpolate_single_field( mesh, climate%ISMIP7%pr_anomaly, time)

      ! Calculate monthly climate
      do vi = mesh%vi1, mesh%vi2
        do mi = 1, 12
          climate%T2m( vi, mi) = climate%ISMIP7%T2m_baseline( vi, mi) + climate%ISMIP7%tas_anomaly%val_interp( vi, mi) + delta_ts( vi)
          climate%Precip( vi, mi) = climate%ISMIP7%Precip_baseline( vi, mi) &
            + climate%ISMIP7%pr_anomaly%val_interp( vi, mi) * sec_per_year / freshwater_density ! [m.w.e. yr^-1] 
        end do
      end do

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_climate_model_ISMIP7

  subroutine initialise_climate_model_ISMIP7( mesh, refgeo_PD, refgeo_init, region_name, ISMIP7)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_reference_geometry),   intent(in   ) :: refgeo_PD
    type(type_reference_geometry),   intent(in   ) :: refgeo_init
    character(len=3),                intent(in   ) :: region_name
    type(type_climate_model_ISMIP7), intent(inout) :: ISMIP7

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_climate_model_ISMIP7'
    character(len=1024)            :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise fields and baseline climate
    select case (C%climate_ISMIP7_choice_baseline)
    case default
      call crash('invalid climate_ISMIP7_choice_baseline "' // trim( C%climate_ISMIP7_choice_baseline) // '"')
    case ('yearly')
      ! Initialise monthly fields
      call initialise_climate_field( mesh, ISMIP7, ISMIP7%tas, 'tas')
      call initialise_climate_field( mesh, ISMIP7, ISMIP7%pr, 'pr')
    case ('fixed')
      ! Initialise monthly fields
      call initialise_climate_field( mesh, ISMIP7, ISMIP7%tas_anomaly, 'tas-anomaly')
      call initialise_climate_field( mesh, ISMIP7, ISMIP7%pr_anomaly, 'pr-anomaly')

      ! Initialise baseline
      call initialise_climate_baseline_fixed( mesh, ISMIP7)
    end select

    ! Initialise vertical gradient
    call initialise_climate_field( mesh, ISMIP7, ISMIP7%dtsdz, 'dtsdz')

    ! Initialise the baseline surface elevation
    select case (C%climate_ISMIP7_choice_refgeo)
    case default
      call crash('invalid climate_ISMIP7_choice_refgeo "' // trim( C%climate_ISMIP7_choice_refgeo) // '"')
    case ('init')
      call initialise_Hs_baseline( mesh, ISMIP7, refgeo_init)
    case ('PD')
      call initialise_Hs_baseline( mesh, ISMIP7, refgeo_PD)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_climate_model_ISMIP7

  subroutine initialise_climate_field_monthly( mesh, ISMIP7, field, name)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_climate_model_ISMIP7),         intent(inout) :: ISMIP7
    type(type_climate_field_ISMIP7_monthly), intent(inout) :: field
    character(len=*),                        intent(in   ) :: name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_climate_field_monthly'
    character(len=1024)            :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Define name
    field%name = name

    ! Get info from files
    call gather_fileinfo( ISMIP7, field%filenames, field%timestamps, trim(field%name))

    ! Allocate memory for timeframes
    allocate (field%val0      ( mesh%vi1:mesh%vi2, 12), source = NaN)
    allocate (field%val1      ( mesh%vi1:mesh%vi2, 12), source = NaN)
    allocate (field%val_interp( mesh%vi1:mesh%vi2, 12), source = NaN)

    ! Update timeframes to the current model time
    call update_timeframes( mesh, field, C%start_time_of_run)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_climate_field_monthly

  subroutine initialise_climate_field_yearly( mesh, ISMIP7, field, name)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_climate_model_ISMIP7),         intent(inout) :: ISMIP7
    type(type_climate_field_ISMIP7_yearly),  intent(inout) :: field
    character(len=*),                        intent(in   ) :: name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_climate_field_yearly'
    character(len=1024)            :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Define name
    field%name = name

    ! Get info from files
    call gather_fileinfo( ISMIP7, field%filenames, field%timestamps, trim(field%name))

    ! Allocate memory for timeframes
    allocate (field%val0      ( mesh%vi1:mesh%vi2), source = NaN)
    allocate (field%val1      ( mesh%vi1:mesh%vi2), source = NaN)
    allocate (field%val_interp( mesh%vi1:mesh%vi2), source = NaN)

    ! Update timeframes to the current model time
    call update_timeframes( mesh, field, C%start_time_of_run)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_climate_field_yearly

  subroutine initialise_climate_baseline_fixed( mesh, ISMIP7)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_climate_model_ISMIP7), intent(inout) :: ISMIP7

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_climate_baseline_fixed'
    character(len=1024)            :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Read the fixed baseline climate
    call read_field_from_file_2D_monthly( C%climate_ISMIP7_filename_baseline, 'T2m'   , mesh, C%output_dir, ISMIP7%T2m_baseline)
    call read_field_from_file_2D_monthly( C%climate_ISMIP7_filename_baseline, 'Precip', mesh, C%output_dir, ISMIP7%Precip_baseline)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_climate_baseline_fixed

  subroutine initialise_Hs_baseline( mesh, ISMIP7, refgeo)

    ! In/output variables
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_climate_model_ISMIP7), intent(inout) :: ISMIP7
    type(type_reference_geometry),   intent(in   ) :: refgeo

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_Hs_baseline'

    ! Add routine to path
    call init_routine( routine_name)

    ISMIP7%Hs_baseline( mesh%vi1:mesh%vi2) = refgeo%Hs( mesh%vi1: mesh%vi2)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_Hs_baseline

  subroutine gather_fileinfo( ISMIP7, filenames, timestamps, var_name)

    ! In/output variables
    class(type_climate_model_ISMIP7),               intent(inout) :: ISMIP7
    character(len=1024), dimension(:), allocatable, intent(inout) :: filenames
    real(dp),            dimension(:), allocatable, intent(inout) :: timestamps
    character(len=*),                               intent(in   ) :: var_name

    ! Local variables:
    character(len=*), parameter                    :: routine_name = 'gather_fileinfo'
    character(len=:), allocatable                  :: foldername
    character(len=1024), dimension(:), allocatable :: list_of_filenames
    integer                                        :: i
    real(dp)                                       :: year

    ! Add routine to path
    call init_routine( routine_name)

    if (allocated( filenames )) deallocate( filenames)
    if (allocated( timestamps)) deallocate( timestamps)

    ! Construct foldernames
    foldername = trim( C%SMB_ISMIP7_forcing_foldername) // '/' // var_name // '/' // trim( C%SMB_ISMIP7_forcing_version)

    call list_files_in_folder( foldername, list_of_filenames, var_name)
    if (size( filenames,1) == 0) call crash('could not find any valid NetCDF files in directory "' // trim( foldername) // '"')

    ! Read timestamps
    allocate( timestamps( size( filenames,1)))
    do i = 1, size( filenames,1)
      call read_year_from_netcdf_filename( filenames( i), year)
      ! Add half a year, so that we have the timestamp at the middle of the year rather than the start
      timestamps( i) = year + 0.5_dp
    end do

    ! Append foldername to filenames
    do i = 1, size( filenames,1)
      filenames( i) = trim( foldername) // '/' // trim( filenames( i))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_fileinfo

  subroutine read_year_from_netcdf_filename( filename, year)
    ! Rather than trying to read from the file's time dimension, which would mean
    ! dealing with whatever calendar it uses, just read it from the filename

    ! In/output variables
    character(len=*), intent(in   ) :: filename
    real(dp),         intent(  out) :: year

    ! Local variables
    character(len=*), parameter :: routine_name = 'read_year_from_netcdf_filename'
    character(len=4) :: year_str
    integer          :: year_int, stat

    ! Add routine to path
    call init_routine( routine_name)

    ! Since the files are called e.g. 'acabf_AIS_CESM2-WACCM_ssp585_SDBN1-8000m_v2_2015.nc',
    ! just assumed that the last few characters spell out the year

    year_str = filename( len_trim( filename) - 6 : len_trim( filename) - 3)
    year_int = UPSY%stru%str2int( year_str, stat)
    if (stat /= 0) call crash('could not read year from filename "' // trim( filename) // '"')

    year = real( year_int,dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_year_from_netcdf_filename

  subroutine update_timeframes_monthly( mesh, field, time)

    ! In/output variables:
    type(type_mesh),                          intent(in   ) :: mesh
    type(type_climate_field_ISMIP7_monthly), intent(inout) :: field
    real(dp),                                 intent(in   ) :: time

    ! Local variables
    character(len=*), parameter :: routine_name = 'update_timeframes_monthly'
    integer                     :: ti0_old, ti1_old

    ! Add routine to path
    call init_routine( routine_name)

    ! Get current bracket indices
    ti0_old = field%ti0
    ti1_old = field%ti1

    ! Update the indices of time slices before and after current time
    call update_bracket_indices( field%timestamps, field%ti0, field%ti1, time)

    ! Update timeframes if necessary
    if (field%ti0 /= ti0_old) then
      call update_single_timeframe_monthly( mesh, field, field%ti0, field%val0)
    end if

    if (field%ti1 /= ti1_old) then
      call update_single_timeframe_monthly( mesh, field, field%ti1, field%val1)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_timeframes_monthly

  subroutine update_timeframes_yearly( mesh, field, time)

    ! In/output variables:
    type(type_mesh),                          intent(in   ) :: mesh
    type(type_climate_field_ISMIP7_yearly), intent(inout) :: field
    real(dp),                                 intent(in   ) :: time

    ! Local variables
    character(len=*), parameter :: routine_name = 'update_timeframes_yearly'
    integer                     :: ti0_old, ti1_old

    ! Add routine to path
    call init_routine( routine_name)

    ! Get current bracket indices
    ti0_old = field%ti0
    ti1_old = field%ti1

    ! Update the indices of time slices before and after current time
    call update_bracket_indices( field%timestamps, field%ti0, field%ti1, time)

    ! Update timeframes if necessary
    if (field%ti0 /= ti0_old) then
      call update_single_timeframe_yearly( mesh, field, field%ti0, field%val0)
    end if

    if (field%ti1 /= ti1_old) then
      call update_single_timeframe_yearly( mesh, field, field%ti1, field%val1)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_timeframes_yearly

  subroutine update_single_timeframe_monthly( mesh, field, ti, val)

    ! In/output variables:
    type(type_mesh),                          intent(in   ) :: mesh
    type(type_climate_field_ISMIP7_monthly),  intent(in   ) :: field
    integer,                                  intent(in   ) :: ti
    real(dp), dimension(:,:),                 intent(inout) :: val

    ! Local variables
    character(len=*), parameter             :: routine_name = 'update_single_timeframe_monthly'
    character(len=:), allocatable           :: filename
    real(dp), dimension(:), allocatable     :: time_from_file
    integer                                 :: mi
    real(dp), dimension( mesh%vi1:mesh%vi2) :: d_month

    ! Add routine to path
    call init_routine( routine_name)

    filename = trim(field%filenames( ti))

    ! Read time dimension from file
    call read_time_from_file( filename, time_from_file)
    if (size( time_from_file,1) /= 12) call crash('file "' // trim( filename) // '" doesnt have 12 months')

    ! Read all 12 months individually with extrapolation
    do mi = 1, 12
      call read_field_from_file_2D( filename, trim(field%name), mesh, C%output_dir, d_month, &
      time_to_read = time_from_file( mi), extrapolate_fillvalues = .true.)
      val( mesh%vi1:mesh%vi2, mi) = d_month( mesh%vi1:mesh%vi2)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_single_timeframe_monthly

  subroutine update_single_timeframe_yearly( mesh, field, ti, val)

    ! In/output variables:
    type(type_mesh),                          intent(in   ) :: mesh
    type(type_climate_field_ISMIP7_yearly),   intent(in   ) :: field
    integer,                                  intent(in   ) :: ti
    real(dp), dimension(:),                   intent(inout) :: val

    ! Local variables
    character(len=*), parameter             :: routine_name = 'update_single_timeframe_yearly'
    character(len=:), allocatable           :: filename
    real(dp), dimension(:), allocatable     :: time_from_file

    ! Add routine to path
    call init_routine( routine_name)

    filename = trim(field%filenames( ti))

    ! Read time dimension from file
    call read_time_from_file( filename, time_from_file)
    if (size( time_from_file,1) /= 1) call crash('file "' // trim( filename) // '" doesnt have 1 timeframe')

    ! Read with extrapolation
    call read_field_from_file_2D( filename, trim(field%name), mesh, C%output_dir, val, &
      time_to_read = time_from_file( 1), extrapolate_fillvalues = .true.)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_single_timeframe_yearly

  subroutine update_bracket_indices( timestamps, ti0, ti1, time)

    ! In/output variables:
    real(dp), dimension(:), intent(in   ) :: timestamps
    integer,                intent(inout) :: ti0
    integer,                intent(inout) :: ti1
    real(dp),               intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'update_bracket_indices'
    integer                        :: i, n

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Determine total numer of timeframes available for this field
    n = size(timestamps)

    if (time <= timestamps(1)) then
      ! Model time before first available time value, return first two indices
      ti0 = 1
      ti1 = 2

    elseif (time >= timestamps(n)) then
      ! Model time after last available time value, return last two indices
      ti0 = n-1
      ti1 = n

    else
      ! Model time within array, return bracketing indices
      do i = 1, n
        if (timestamps(i) <= time) then
          ti0 = i
          ti1 = i+1
        end if
      end do

    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine update_bracket_indices

  subroutine interpolate_single_field_monthly( mesh, field, time)

    ! In/output variables:
    type(type_mesh),                             intent(in   ) :: mesh
    type(type_climate_field_ISMIP7_monthly),     intent(inout) :: field
    real(dp),                                    intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'interpolate_single_field_monthly'
    real(dp)                       :: w0, w1
    real(dp)                       :: days

    ! Add routine to call stack
    call init_routine( routine_name)

    call get_interpolation_weights( field%timestamps( field%ti0), field%timestamps( field%ti1), time, w0, w1)

    ! Apply interpolation to values
    field%val_interp( mesh%vi1:mesh%vi2, :) = w0 * field%val0( mesh%vi1:mesh%vi2, :) + w1 * field%val1( mesh%vi1:mesh%vi2, :)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine interpolate_single_field_monthly

  subroutine interpolate_single_field_yearly( mesh, field, time)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_climate_field_ISMIP7_yearly),  intent(inout) :: field
    real(dp),                                intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'interpolate_single_field_yearly'
    real(dp)                       :: w0, w1
    real(dp)                       :: days

    ! Add routine to call stack
    call init_routine( routine_name)

    call get_interpolation_weights( field%timestamps( field%ti0), field%timestamps( field%ti1), time, w0, w1)

    ! Apply interpolation to values
    field%val_interp( mesh%vi1:mesh%vi2) = w0 * field%val0( mesh%vi1:mesh%vi2) + w1 * field%val1( mesh%vi1:mesh%vi2)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine interpolate_single_field_yearly

  subroutine get_interpolation_weights( timestamp0, timestamp1, time, w0, w1)

    ! In/output variables:
    real(dp), intent(in   ) :: timestamp0
    real(dp), intent(in   ) :: timestamp1
    real(dp), intent(in   ) :: time
    real(dp), intent(  out) :: w0
    real(dp), intent(  out) :: w1

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'get_interpolation_weights'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Get weights
    if (time < timestamp0) then
      ! Model time before bracket times, put full weight on the first timeframe
      w0 = 1._dp
    elseif (time > timestamp1) then
      ! Model time after bracket times, put full weight on the last timeframe
      w0 = 0._dp
    else
      ! Model time between bracket times, determine interpolated weight
      w0 = (timestamp1 - time) / (timestamp1 - timestamp0)
    end if

    w1 = 1._dp - w0

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine get_interpolation_weights

end module climate_ISMIP7
