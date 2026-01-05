module fields_basic

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid, type_Arakawa_grid
  use fields_dimensions, only: third_dimension, type_third_dimension
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: &
    atype_field, atype_field_2D, atype_field_3D, &
    type_field_logical_2D, type_field_int_2D, type_field_dp_2D, &
    type_field_logical_3D, type_field_int_3D, type_field_dp_3D

  ! Abstract basic field type
  ! =========================

  type, abstract :: atype_field

    ! Metadata
    character(len=1024), private :: name_val
    character(len=1024), private :: long_name_val
    character(len=1024), private :: units_val

    ! Grid
    class(*), pointer,          private :: grid_val
    type(type_Arakawa_grid),    private :: Arakawa_grid_val

    ! Pointer to array containing the actual field data
    ! class(*), dimension(:), pointer, public :: d
    type(MPI_WIN),                   public :: w

    contains

    procedure, public :: lbound => field_lbound
    procedure, public :: ubound => field_ubound
    procedure, public :: print_info

    ! ===== Set/get functions

    ! Metadata
    procedure, public :: set_name
    procedure, public :: set_long_name
    procedure, public :: set_units

    procedure, public :: name      => get_name
    procedure, public :: long_name => get_long_name
    procedure, public :: units     => get_units

    procedure, public :: is_name
    procedure, public :: is_long_name
    procedure, public :: is_units

    ! Grid
    procedure, public :: set_grid
    procedure, public :: set_Arakawa_grid

    procedure, public :: grid         => get_grid
    procedure, public :: Arakawa_grid => get_Arakawa_grid

    procedure, public :: is_grid
    procedure, public :: is_Arakawa_grid
    procedure, public :: is_third_dimension

  end type atype_field

  ! Concrete field types types for different data types/ranks
  ! =========================================================

  type, abstract, extends(atype_field) :: atype_field_2D
  end type atype_field_2D

  type, abstract, extends(atype_field) :: atype_field_3D
    type(type_third_dimension), private :: third_dimension_val
    contains
    procedure, public :: set_third_dimension
    procedure, public :: third_dimension => get_third_dimension
  end type atype_field_3D

  type, extends( atype_field_2D) :: type_field_logical_2D
    logical, dimension(:), contiguous, pointer :: d
  end type type_field_logical_2D

  type, extends( atype_field_2D) :: type_field_int_2D
    integer, dimension(:), contiguous, pointer :: d
  end type type_field_int_2D

  type, extends( atype_field_2D) :: type_field_dp_2D
    real(dp), dimension(:), contiguous, pointer :: d
  end type type_field_dp_2D

  type, extends( atype_field_3D) :: type_field_logical_3D
    logical, dimension(:,:), contiguous, pointer :: d
  end type type_field_logical_3D

  type, extends( atype_field_3D) :: type_field_int_3D
    integer, dimension(:,:), contiguous, pointer :: d
  end type type_field_int_3D

  type, extends( atype_field_3D) :: type_field_dp_3D
    real(dp), dimension(:,:), contiguous, pointer :: d
  end type type_field_dp_3D

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  interface

    module function field_lbound( field, dim) result( lb)
      class(atype_field), intent(in) :: field
      integer,            intent(in) :: dim
      integer                        :: lb
    end function field_lbound

    module function field_ubound( field, dim) result( ub)
      class(atype_field), intent(in) :: field
      integer,            intent(in) :: dim
      integer                        :: ub
    end function field_ubound

    module subroutine print_info( field)
      class(atype_field), intent(in) :: field
    end subroutine print_info

    ! ===== Set/get functions

    ! Metadata

    module subroutine set_name( field, name)
      class(atype_field), intent(inout) :: field
      character(len=*),   intent(in   ) :: name
    end subroutine set_name

    module subroutine set_long_name( field, long_name)
      class(atype_field), intent(inout) :: field
      character(len=*),   intent(in   ) :: long_name
    end subroutine set_long_name

    module subroutine set_units( field, units)
      class(atype_field), intent(inout) :: field
      character(len=*),   intent(in   ) :: units
    end subroutine set_units

    module function get_name( field) result( name)
      class(atype_field), intent(in) :: field
      character(:), allocatable      :: name
    end function get_name

    module function get_long_name( field) result( long_name)
      class(atype_field), intent(in) :: field
      character(:), allocatable      :: long_name
    end function get_long_name

    module function get_units( field) result( units)
      class(atype_field), intent(in) :: field
      character(:), allocatable      :: units
    end function get_units

    module function is_name( field, name) result( res)
      class(atype_field), intent(in) :: field
      character(len=*),   intent(in) :: name
      logical                        :: res
    end function is_name

    module function is_long_name( field, long_name) result( res)
      class(atype_field), intent(in) :: field
      character(len=*),   intent(in) :: long_name
      logical                        :: res
    end function is_long_name

    module function is_units( field, units) result( res)
      class(atype_field), intent(in) :: field
      character(len=*),   intent(in) :: units
      logical                        :: res
    end function is_units

    ! Grid

    module subroutine set_grid( field, grid)
      class(atype_field), intent(inout) :: field
      class(*), target,   intent(in   ) :: grid
    end subroutine set_grid

    module subroutine set_Arakawa_grid( field, field_Arakawa_grid)
      class(atype_field),      intent(inout) :: field
      type(type_Arakawa_grid), intent(in   ) :: field_Arakawa_grid
    end subroutine set_Arakawa_grid

    module function get_grid( field) result( grid)
      class(atype_field), intent(in) :: field
      class(*), pointer              :: grid
    end function get_grid

    module function get_Arakawa_grid( field) result( field_Arakawa_grid)
      class(atype_field), intent(in) :: field
      type(type_Arakawa_grid)        :: field_Arakawa_grid
    end function get_Arakawa_grid

    module function is_grid( field, grid) result( res)
      class(atype_field), intent(in) :: field
      class(*),           intent(in) :: grid
      logical                        :: res
    end function is_grid

    module function is_Arakawa_grid( field, field_Arakawa_grid) result( res)
      class(atype_field),      intent(in) :: field
      type(type_Arakawa_grid), intent(in) :: field_Arakawa_grid
      logical                             :: res
    end function is_Arakawa_grid

    ! Third dimension, only for concrete type type_field_3D

    module subroutine set_third_dimension( field, field_third_dimension)
      class(atype_field_3D),      intent(inout) :: field
      type(type_third_dimension), intent(in   ) :: field_third_dimension
    end subroutine set_third_dimension

    module function get_third_dimension( field) result( field_third_dimension)
      class(atype_field_3D), intent(in) :: field
      type(type_third_dimension)        :: field_third_dimension
    end function get_third_dimension

    module function is_third_dimension( field, field_third_dimension) result( res)
      ! ..except for this one, which also works on type_field_2D - it just always returns .false. there
      class(atype_field),         intent(in) :: field
      type(type_third_dimension), intent(in) :: field_third_dimension
      logical                                :: res
    end function is_third_dimension

  end interface

contains

end module fields_basic