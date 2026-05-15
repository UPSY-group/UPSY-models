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
  use reference_geometry_types, only: type_reference_geometry
  use climate_model_types, only: type_climate_model, type_climate_model_ISMIP7, &
    type_climate_field_ISMIP7_monthly, type_climate_field_ISMIP7_annual
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use mpi_f08, only: MPI_WIN
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use netcdf_io_main, only: read_field_from_file_2D, read_time_from_file
  use basic_model_utilities, only: list_files_in_folder
  use parameters, only: sec_per_year, ice_density

  implicit none

  private

  public :: initialise_climate_model_ISMIP7, run_climate_model_ISMIP7

contains

  subroutine run_climate_model_ISMIP7( mesh, climate, time)

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    type(type_climate_model), intent(inout) :: climate
    real(dp),                 intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_climate_model_ISMIP7'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Update timeframes if necessary

    ! Interpolate between timeframes

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

    ! Get info from files

    ! Allocate memory for timeframes

    ! Update timeframes to the current model time

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_climate_model_ISMIP7

end module climate_ISMIP7
