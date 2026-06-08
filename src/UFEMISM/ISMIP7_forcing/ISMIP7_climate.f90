module ISMIP7_climate

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
  use climate_model_types, only: type_climate_model
  use netcdf_io_main, only: read_field_from_file_2D_monthly, read_field_from_file_2D
  use parameters, only: NaN, sec_per_year, freshwater_density
  use ISMIP7_climate_model_type, only: type_climate_model_ISMIP7
  use ISMIP7_forcing_field_types, only: update_timeframes, interpolate_single_field

  implicit none

  private

  public :: initialise_climate_model_ISMIP7, run_climate_model_ISMIP7

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
          climate%Precip( vi, mi) = max(0._dp, climate%ISMIP7%Precip_baseline( vi, mi) &
            + climate%ISMIP7%pr_anomaly%val_interp( vi, mi) * sec_per_year / freshwater_density) ! [m.w.e. yr^-1], must be positive
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
      call ISMIP7%tas%initialise( C%climate_ISMIP7_forcing_foldername, C%climate_ISMIP7_forcing_version, &
        mesh, 'tas')
      call ISMIP7%pr%initialise( C%climate_ISMIP7_forcing_foldername, C%climate_ISMIP7_forcing_version, &
        mesh, 'pr')
    case ('fixed')
      ! Initialise monthly fields
      call ISMIP7%tas_anomaly%initialise( C%climate_ISMIP7_forcing_foldername, C%climate_ISMIP7_forcing_version, &
        mesh, 'tas-anomaly')
      call ISMIP7%pr_anomaly%initialise( C%climate_ISMIP7_forcing_foldername, C%climate_ISMIP7_forcing_version, &
        mesh, 'pr-anomaly')

      ! Initialise baseline
      call initialise_climate_baseline_fixed( mesh, ISMIP7)
    end select

    ! Initialise vertical gradient
    call ISMIP7%dtsdz%initialise( C%climate_ISMIP7_forcing_foldername, C%climate_ISMIP7_forcing_version, &
      mesh, 'dtsdz')

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

  subroutine initialise_climate_baseline_fixed( mesh, ISMIP7)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_climate_model_ISMIP7), intent(inout) :: ISMIP7

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_climate_baseline_fixed'
    character(len=1024)            :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate baseline climate
    if (allocated( ISMIP7%T2m_baseline    )) deallocate( ISMIP7%T2m_baseline   )
    if (allocated( ISMIP7%Precip_baseline )) deallocate( ISMIP7%Precip_baseline)

    allocate (ISMIP7%T2m_baseline    ( mesh%vi1:mesh%vi2, 12), source = NaN)
    allocate (ISMIP7%Precip_baseline ( mesh%vi1:mesh%vi2, 12), source = NaN)

    if (par%primary) then
      write(0,*) '   Reading ISMIP7 climate baseline from file: ', &
        UPSY%stru%colour_string( trim( C%climate_ISMIP7_filename_baseline), 'light blue')
    end if

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

    if (allocated( ISMIP7%Hs_baseline    )) deallocate( ISMIP7%Hs_baseline   )
    allocate (ISMIP7%Hs_baseline ( mesh%vi1:mesh%vi2), source = NaN)

    ISMIP7%Hs_baseline( mesh%vi1:mesh%vi2) = refgeo%Hs( mesh%vi1: mesh%vi2)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_Hs_baseline

end module ISMIP7_climate
