module ice_model_main

  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_basic, only: type_ice_geometry_model
  use ice_velocity_model_basic, only: atype_ice_velocity_model

  implicit none

  private

  public :: type_ice_model

  type, extends(atype_ice_model_data) :: type_ice_model

    ! Geometry
    type(type_ice_geometry_model), allocatable :: geom

    ! Velocity
    class(atype_ice_velocity_model), allocatable :: vel

  end type type_ice_model

end module ice_model_main
