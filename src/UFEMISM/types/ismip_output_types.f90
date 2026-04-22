module ismip_output_types

  use precisions, only: dp
  use grid_types, only: type_grid

  implicit none

  private

  public :: type_ismip_grid_output, type_ismip_gridded_field

  type type_ismip_gridded_field
    ! A gridded output field

    !real(dp), dimension(:), allocatable :: val              ! The value of this field on the mesh
    character(len=1024)                 :: name              ! Variable name
    character(len=1024)                 :: filename          ! Filename of output
    character(len=1024)                 :: long_name         ! Long name
    character(len=1024)                 :: units             ! Units
    character(len=2)                    :: fieldtype         ! ST or FL      

  end type type_ismip_gridded_field

  type type_ismip_grid_output
    ! Data fields for storing the regional ISMIP data

    character(len=3)                  :: IS_name      ! Ice sheet name (AIS or GIS)
    type(type_grid)                   :: grid         ! Output grid

    type(type_ismip_gridded_field)    :: lithk        ! [m]    land_ice_thickness
    type(type_ismip_gridded_field)    :: orog         ! [m]    surface_altitude
    type(type_ismip_gridded_field)    :: topg         ! [m]    bedrock_altitude
    type(type_ismip_gridded_field)    :: xvelsurf     ! [m]    land_ice_surface_x_velocity
    type(type_ismip_gridded_field)    :: yvelsurf     ! [m]    land_ice_surface_y_velocity
    type(type_ismip_gridded_field)    :: zvelsurf     ! [m]    land_ice_surface_z_velocity


  end type type_ismip_grid_output
end module ismip_output_types
