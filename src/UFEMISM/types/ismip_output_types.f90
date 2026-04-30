module ismip_output_types

  use precisions, only: dp
  use grid_types, only: type_grid

  implicit none

  private

  public :: type_ismip_grid_output, type_ismip_gridded_field

  type type_ismip_gridded_field
    ! A gridded output field

    real(dp), dimension(:), allocatable :: accum             ! [..] Accumulation of values per timestep
    logical                             :: is_initial        ! Whether this is the first time to write

    character(len=1024)                 :: name              ! Variable name
    character(len=1024)                 :: filename          ! Filename of output
    character(len=1024)                 :: long_name         ! Long name
    character(len=1024)                 :: standard_name     ! Standard name
    character(len=1024)                 :: units             ! Units
    character(len=2)                    :: fieldtype         ! ST or FL      

    integer                             :: nt                ! Counter for number of time values written

  end type type_ismip_gridded_field

  type type_ismip_grid_output
    ! Data fields for storing the regional ISMIP data

    character(len=3)                  :: IS_name      ! Ice sheet name (AIS or GIS)
    character(len=126)                :: crs          ! Associated crs (EPSG)
    type(type_grid)                   :: grid         ! Output grid
    character(len=1024)               :: folder       ! Subfolder exp_RES

    real(dp)                          :: t_prev       ! [yr] Time value at previous writing
    real(dp)                          :: t_curr       ! [yr] Current time value

    ! Basic topography
    type(type_ismip_gridded_field)    :: lithk        ! [m]        land_ice_thickness
    type(type_ismip_gridded_field)    :: orog         ! [m]        surface_altitude
    type(type_ismip_gridded_field)    :: topg         ! [m]        bedrock_altitude
    type(type_ismip_gridded_field)    :: base         ! [m]        base_altitude

    ! Geothermal heat flux
    type(type_ismip_gridded_field)    :: hfgeoubed    ! [W m-2]    upward_geothermal_heat_flux_in_land_ice

    ! Surface and basal mass balance
    type(type_ismip_gridded_field)    :: acabf        ! [kg m-2 s-1] land_ice_surface_specific_mass_balance_flux
    type(type_ismip_gridded_field)    :: libmassbfgr  ! [kg m-2 s-1] land_ice_basal_specific_mass_balance_flux (gr)
    type(type_ismip_gridded_field)    :: libmassbffl  ! [kg m-2 s-1] land_ice_basal_specific_mass_balance_flux (fl)

    ! Thickness tendency
    type(type_ismip_gridded_field)    :: dlithkdt     ! [m s-1]    tendency_of_land_ice_thickness

    ! Velocities
    type(type_ismip_gridded_field)    :: xvelsurf     ! [m s-1]    land_ice_surface_x_velocity
    type(type_ismip_gridded_field)    :: yvelsurf     ! [m s-1]    land_ice_surface_y_velocity
    type(type_ismip_gridded_field)    :: zvelsurf     ! [m s-1]    land_ice_surface_upward_velocity
    type(type_ismip_gridded_field)    :: xvelbase     ! [m s-1]    land_ice_basal_x_velocity
    type(type_ismip_gridded_field)    :: yvelbase     ! [m s-1]    land_ice_basal_y_velocity
    type(type_ismip_gridded_field)    :: zvelbase     ! [m s-1]    land_ice_base_upward_velocity
    type(type_ismip_gridded_field)    :: xvelmean     ! [m s-1]    land_ice_vertical_mean_x_velocity
    type(type_ismip_gridded_field)    :: yvelmean     ! [m s-1]    land_ice_vertical_mean_y_velocity

    ! Temperatures
    type(type_ismip_gridded_field)    :: litemptop    ! [K]        temperature_at_top_of_ice_sheet_model
    type(type_ismip_gridded_field)    :: litempavg    ! [K]        land_ice_temperature
    type(type_ismip_gridded_field)    :: litempgradgr ! [K m-1]    Vertical basal temperature gradient beneath grounded ice
    type(type_ismip_gridded_field)    :: litempgradfl ! [K m-1]    Vertical basal temperature gradient beneath floating ice
    type(type_ismip_gridded_field)    :: litempbotgr  ! [K]        temperature_at_base_of_ice_sheet_model (gr)
    type(type_ismip_gridded_field)    :: litempbotfl  ! [K]        temperature_at_base_of_ice_sheet_model (fl)

    ! Basal drag
    type(type_ismip_gridded_field)    :: strbasemag   ! [Pa]       land_ice_basal_drag

    ! Lateral mass balance
    type(type_ismip_gridded_field)    :: licalvf      ! [kg m-2 s-1] land_ice_specific_mass_flux_due_to_calving
    type(type_ismip_gridded_field)    :: lifmassbf    ! [kg m-2 s-1] Loss if ice mass resulting from ice front melting

    ! Area fractions
    type(type_ismip_gridded_field)    :: sftgif       ! []         land_ice_area_fraction
    type(type_ismip_gridded_field)    :: sftgrf       ! []         grounded_ice_sheet_area_fraction
    type(type_ismip_gridded_field)    :: sftflf       ! []         floating_ice_shelf_area_fraction

  end type type_ismip_grid_output
end module ismip_output_types
