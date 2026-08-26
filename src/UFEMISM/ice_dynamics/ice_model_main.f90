module ice_model_main

  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_basic, only: type_ice_geometry_model
  use ice_velocity_model_basic, only: atype_ice_velocity_model
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use ISMIP7_fracture, only: type_ISMIP7_fracture_model

  implicit none

  private

  public :: type_ice_model

  type, extends(atype_ice_model_data) :: type_ice_model

      class(type_ice_geometry_model),       allocatable :: geom                      !< Ice-sheet geometry
      class(atype_ice_velocity_model),      allocatable :: vel                       !< Ice-sheet velocity
      class(atype_momentum_balance_solver), allocatable :: momentum_balance_solver   !< The momentum balance solver
      class(type_ISMIP7_fracture_model),    allocatable :: ISMIP7_fracture           !< The ISMIP7 mask-based hydrofracture model

    contains

      procedure, public :: get_model_name

  end type type_ice_model

contains

  function get_model_name( self) result( model_name)
    class(type_ice_model), intent(in) :: self
    character(len=:), allocatable     :: model_name
    model_name = 'ice'
  end function get_model_name

end module ice_model_main
