module ice_model_main

  use ice_model_data, only: atype_ice_model_data

  implicit none

  private

  public :: type_ice_model

  type, extends(atype_ice_model_data) :: type_ice_model
  end type type_ice_model

end module ice_model_main
