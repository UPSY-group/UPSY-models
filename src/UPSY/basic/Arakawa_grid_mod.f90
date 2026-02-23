module Arakawa_grid_mod

  use call_stack_and_comp_time_tracking, only: crash

  implicit none

  private

  public :: Arakawa_grid, type_Arakawa_grid

  type type_Arakawa_grid_factory
  contains
    procedure :: a  => Arakawa_grid_factory_a
    procedure :: b  => Arakawa_grid_factory_b
    procedure :: c  => Arakawa_grid_factory_c
    procedure :: cx => Arakawa_grid_factory_cx
    procedure :: cy => Arakawa_grid_factory_cy
  end type type_Arakawa_grid_factory

  type type_Arakawa_grid
    integer, private :: g
  contains
    procedure, private :: eq => test_Arakawa_grid_equality
    generic :: operator(==) => eq
    procedure, public :: str => Arakawa_grid_str
  end type type_Arakawa_grid

  type(type_Arakawa_grid_factory) :: Arakawa_grid

contains

  function Arakawa_grid_factory_a( self) result( res)
    class(type_Arakawa_grid_factory), intent(in) :: self
    type(type_Arakawa_grid) :: res
    res%g = 1
  end function Arakawa_grid_factory_a

  function Arakawa_grid_factory_b( self) result( res)
    class(type_Arakawa_grid_factory), intent(in) :: self
    type(type_Arakawa_grid) :: res
    res%g = 2
  end function Arakawa_grid_factory_b

  function Arakawa_grid_factory_c( self) result( res)
    class(type_Arakawa_grid_factory), intent(in) :: self
    type(type_Arakawa_grid) :: res
    res%g = 3
  end function Arakawa_grid_factory_c

  function Arakawa_grid_factory_cx( self) result( res)
    class(type_Arakawa_grid_factory), intent(in) :: self
    type(type_Arakawa_grid) :: res
    res%g = 4
  end function Arakawa_grid_factory_cx

  function Arakawa_grid_factory_cy( self) result( res)
    class(type_Arakawa_grid_factory), intent(in) :: self
    type(type_Arakawa_grid) :: res
    res%g = 5
  end function Arakawa_grid_factory_cy

  pure function test_Arakawa_grid_equality( self, other) result( res)
    class(type_Arakawa_grid), intent(in) :: self
    class(type_Arakawa_grid), intent(in) :: other
    logical                              :: res
    res = self%g == other%g
  end function test_Arakawa_grid_equality

  function Arakawa_grid_str( self) result( str)
    class(type_Arakawa_grid), intent(in) :: self
    character(len=2)                     :: str
    if (self == Arakawa_grid%a()) then
      str = 'a'
    elseif (self == Arakawa_grid%b()) then
      str = 'b'
    elseif (self == Arakawa_grid%c()) then
      str = 'c'
    elseif (self == Arakawa_grid%cx()) then
      str = 'cx'
    elseif (self == Arakawa_grid%cy()) then
      str = 'cy'
    else
      call crash('invalid Arakawa grid in Arakawa_grid_str')
    end if
  end function Arakawa_grid_str

end module Arakawa_grid_mod