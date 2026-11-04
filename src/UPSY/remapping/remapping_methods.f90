module remapping_methods

  use call_stack_and_comp_time_tracking, only: crash

  implicit none

  private

  public :: remapping_method, type_remapping_method, type_remapping_method_factory

  type type_remapping_method_factory
  contains
    procedure :: nearest_neighbour         => remapping_method_factory_nearest_neighbour
    procedure :: linear                    => remapping_method_factory_linear
    procedure :: first_order_conservative  => remapping_method_factory_first_order_conservative
    procedure :: second_order_conservative => remapping_method_factory_second_order_conservative
  end type type_remapping_method_factory

  type type_remapping_method
    integer, private :: g
  contains
    procedure, private :: eq => test_remapping_method_equality
    generic :: operator(==) => eq
    procedure, public :: str => remapping_method_str
  end type type_remapping_method

  type(type_remapping_method_factory) :: remapping_method

contains

  function remapping_method_factory_nearest_neighbour( self) result( res)
    class(type_remapping_method_factory), intent(in) :: self
    type(type_remapping_method) :: res
    res%g = 1
  end function remapping_method_factory_nearest_neighbour

  function remapping_method_factory_linear( self) result( res)
    class(type_remapping_method_factory), intent(in) :: self
    type(type_remapping_method) :: res
    res%g = 2
  end function remapping_method_factory_linear

  function remapping_method_factory_first_order_conservative( self) result( res)
    class(type_remapping_method_factory), intent(in) :: self
    type(type_remapping_method) :: res
    res%g = 3
  end function remapping_method_factory_first_order_conservative

  function remapping_method_factory_second_order_conservative( self) result( res)
    class(type_remapping_method_factory), intent(in) :: self
    type(type_remapping_method) :: res
    res%g = 4
  end function remapping_method_factory_second_order_conservative

  pure function test_remapping_method_equality( self, other) result( res)
    class(type_remapping_method), intent(in) :: self
    class(type_remapping_method), intent(in) :: other
    logical                              :: res
    res = self%g == other%g
  end function test_remapping_method_equality

  function remapping_method_str( self) result( str)
    class(type_remapping_method), intent(in) :: self
    character(len=:), allocatable :: str
    if     (self == remapping_method%nearest_neighbour()) then
      str = 'nearest_neighbour'
    elseif (self == remapping_method%linear()) then
      str = 'linear'
    elseif (self == remapping_method%first_order_conservative()) then
      str = 'first_order_conservative'
    elseif (self == remapping_method%second_order_conservative()) then
      str = 'second_order_conservative'
    else
      call crash('invalid remapping method in remapping_method_str')
    end if
  end function remapping_method_str

end module remapping_methods