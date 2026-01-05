module fields_dimensions

  use precisions, only: dp
  use control_resources_and_error_messaging, only: crash

  implicit none

  private

  public :: third_dimension, type_third_dimension

  type type_third_dimension_factory
  contains
    procedure :: ice_zeta    => third_dimension_factory_ice_zeta
    procedure :: month       => third_dimension_factory_month
    procedure :: ocean_depth => third_dimension_factory_ocean_depth
  end type type_third_dimension_factory

  type type_third_dimension
    character(len=1024) :: name
    integer             :: n
  contains
    procedure, private :: eq => test_third_dimension_equality
    generic :: operator(==) => eq
  end type type_third_dimension

  type(type_third_dimension_factory) :: third_dimension

contains

  function third_dimension_factory_ice_zeta( self, n) result( res)
    class(type_third_dimension_factory), intent(in) :: self
    integer,                             intent(in) :: n
    type(type_third_dimension)                      :: res
    res%name = 'ice_zeta'
    res%n    = n
  end function third_dimension_factory_ice_zeta

  function third_dimension_factory_month( self) result( res)
    class(type_third_dimension_factory), intent(in) :: self
    type(type_third_dimension)                      :: res
    res%name = 'month'
    res%n    = 12
  end function third_dimension_factory_month

  function third_dimension_factory_ocean_depth( self, n) result( res)
    class(type_third_dimension_factory), intent(in) :: self
    integer,                             intent(in) :: n
    type(type_third_dimension)                      :: res
    res%name = 'ocean_depth'
    res%n    = n
  end function third_dimension_factory_ocean_depth

  pure function test_third_dimension_equality( self, other) result( res)
    class(type_third_dimension), intent(in) :: self
    class(type_third_dimension), intent(in) :: other
    logical                                 :: res
    res = self%name == other%name .and. self%n == other%n
  end function test_third_dimension_equality

end module fields_dimensions