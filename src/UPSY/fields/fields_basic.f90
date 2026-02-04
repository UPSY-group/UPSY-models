module fields_basic

  use precisions, only: dp
  use parameters, only: NaN
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, &
    warning, crash, colour_string
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid, type_Arakawa_grid
  use fields_dimensions, only: third_dimension, type_third_dimension
  use parallel_array_info_type, only: type_par_arr_info
  use mpi_f08, only: MPI_WIN
  use remapping_main, only: &
    map_from_mesh_to_mesh_with_reallocation_2D, &
    map_from_mesh_to_mesh_with_reallocation_3D, &
    map_from_mesh_tri_to_mesh_tri_with_reallocation_2D, &
    map_from_mesh_tri_to_mesh_tri_with_reallocation_3D
  use model_configuration, only: C
  use mpi_distributed_shared_memory, only: hybrid_to_dist, dist_to_hybrid
  use reallocate_dist_shared_mod, only: reallocate_dist_shared

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
    character(len=1024), private :: remap_method_val

    ! Grid
    class(*), pointer,       private :: grid_val
    type(type_Arakawa_grid), private :: Arakawa_grid_val
    type(type_par_arr_info), private :: pai_val

    ! Pointer to array containing the actual field data
    ! class(*), dimension(:), pointer, public :: d
    type(MPI_WIN), pointer, public :: w

    contains

    generic,   public  :: operator(==) => eq
    procedure, private :: eq => test_field_equality

    generic,   public  :: reallocate => &
      reallocate_logical_2D, &
      reallocate_logical_3D, &
      reallocate_int_2D, &
      reallocate_int_3D, &
      reallocate_dp_2D, &
      reallocate_dp_3D
    procedure, private :: reallocate_logical_2D
    procedure, private :: reallocate_logical_3D
    procedure, private :: reallocate_int_2D
    procedure, private :: reallocate_int_3D
    procedure, private :: reallocate_dp_2D
    procedure, private :: reallocate_dp_3D

    generic,   public  :: remap => &
      reallocate_logical_2D, &
      reallocate_logical_3D, &
      reallocate_int_2D, &
      reallocate_int_3D, &
      remap_dp_2D, &
      remap_dp_3D
    procedure, private :: remap_dp_2D
    procedure, private :: remap_dp_3D

    procedure, public :: lbound => field_lbound
    procedure, public :: ubound => field_ubound
    procedure, public :: print_info

    ! ===== Set/get functions

    ! Metadata
    procedure, public :: set_name
    procedure, public :: set_long_name
    procedure, public :: set_units
    procedure, public :: set_remap_method

    procedure, public :: name         => get_name
    procedure, public :: long_name    => get_long_name
    procedure, public :: units        => get_units
    procedure, public :: remap_method => get_remap_method

    procedure, public :: is_name
    procedure, public :: is_long_name
    procedure, public :: is_units
    procedure, public :: is_remap_method

    ! Grid
    procedure, public :: set_grid
    procedure, public :: set_Arakawa_grid
    procedure, public :: set_pai

    procedure, public :: grid         => get_grid
    procedure, public :: Arakawa_grid => get_Arakawa_grid
    procedure, public :: pai          => get_pai

    procedure, public :: is_grid
    procedure, public :: is_Arakawa_grid
    procedure, public :: is_third_dimension
    procedure, public :: is_pai

    ! ===== i/o

    procedure, public :: write_to_netcdf
    procedure, public :: read_from_netcdf

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
    logical, dimension(:), contiguous, pointer :: d_nih
  end type type_field_logical_2D

  type, extends( atype_field_2D) :: type_field_int_2D
    integer, dimension(:), contiguous, pointer :: d_nih
  end type type_field_int_2D

  type, extends( atype_field_2D) :: type_field_dp_2D
    real(dp), dimension(:), contiguous, pointer :: d_nih
  end type type_field_dp_2D

  type, extends( atype_field_3D) :: type_field_logical_3D
    logical, dimension(:,:), contiguous, pointer :: d_nih
  end type type_field_logical_3D

  type, extends( atype_field_3D) :: type_field_int_3D
    integer, dimension(:,:), contiguous, pointer :: d_nih
  end type type_field_int_3D

  type, extends( atype_field_3D) :: type_field_dp_3D
    real(dp), dimension(:,:), contiguous, pointer :: d_nih
  end type type_field_dp_3D

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  interface

    ! ===== Basics

    module function test_field_equality( field1, field2) result( res)
      class(atype_field), intent(in) :: field1, field2
      logical                        :: res
    end function test_field_equality

    module function field_lbound( self, dim) result( lb)
      class(atype_field), intent(in) :: self
      integer,            intent(in) :: dim
      integer                        :: lb
    end function field_lbound

    module function field_ubound( self, dim) result( ub)
      class(atype_field), intent(in) :: self
      integer,            intent(in) :: dim
      integer                        :: ub
    end function field_ubound

    module subroutine print_info( self)
      class(atype_field), intent(in) :: self
    end subroutine print_info

    ! ===== Set/get functions

    ! Metadata

    module subroutine set_name( self, name)
      class(atype_field), intent(inout) :: self
      character(len=*),   intent(in   ) :: name
    end subroutine set_name

    module subroutine set_long_name( self, long_name)
      class(atype_field), intent(inout) :: self
      character(len=*),   intent(in   ) :: long_name
    end subroutine set_long_name

    module subroutine set_units( self, units)
      class(atype_field), intent(inout) :: self
      character(len=*),   intent(in   ) :: units
    end subroutine set_units

    module subroutine set_remap_method( self, remap_method)
      class(atype_field), intent(inout) :: self
      character(len=*),   intent(in   ) :: remap_method
    end subroutine set_remap_method

    module function get_name( self) result( name)
      class(atype_field), intent(in) :: self
      character(:), allocatable      :: name
    end function get_name

    module function get_long_name( self) result( long_name)
      class(atype_field), intent(in) :: self
      character(:), allocatable      :: long_name
    end function get_long_name

    module function get_units( self) result( units)
      class(atype_field), intent(in) :: self
      character(:), allocatable      :: units
    end function get_units

    module function get_remap_method( self) result( remap_method)
      class(atype_field), intent(in) :: self
      character(:), allocatable      :: remap_method
    end function get_remap_method

    module function is_name( self, name) result( res)
      class(atype_field), intent(in) :: self
      character(len=*),   intent(in) :: name
      logical                        :: res
    end function is_name

    module function is_long_name( self, long_name) result( res)
      class(atype_field), intent(in) :: self
      character(len=*),   intent(in) :: long_name
      logical                        :: res
    end function is_long_name

    module function is_units( self, units) result( res)
      class(atype_field), intent(in) :: self
      character(len=*),   intent(in) :: units
      logical                        :: res
    end function is_units

    module function is_remap_method( self, remap_method) result( res)
      class(atype_field), intent(in) :: self
      character(len=*),   intent(in) :: remap_method
      logical                        :: res
    end function is_remap_method

    ! Grid

    module subroutine set_grid( self, grid)
      class(atype_field), intent(inout) :: self
      class(*), target,   intent(in   ) :: grid
    end subroutine set_grid

    module subroutine set_Arakawa_grid( self, field_Arakawa_grid)
      class(atype_field),      intent(inout) :: self
      type(type_Arakawa_grid), intent(in   ) :: field_Arakawa_grid
    end subroutine set_Arakawa_grid

    module subroutine set_pai( self, field_pai)
      class(atype_field),      intent(inout) :: self
      type(type_par_arr_info), intent(in   ) :: field_pai
    end subroutine set_pai

    module subroutine set_third_dimension( self, field_third_dimension)
      class(atype_field_3D),      intent(inout) :: self
      type(type_third_dimension), intent(in   ) :: field_third_dimension
    end subroutine set_third_dimension

    module function get_grid( self) result( grid)
      class(atype_field), intent(in) :: self
      class(*), pointer              :: grid
    end function get_grid

    module function get_Arakawa_grid( self) result( field_Arakawa_grid)
      class(atype_field), intent(in) :: self
      type(type_Arakawa_grid)        :: field_Arakawa_grid
    end function get_Arakawa_grid

    module function get_pai( self) result( field_pai)
      class(atype_field), intent(in) :: self
      type(type_par_arr_info)        :: field_pai
    end function get_pai

    module function get_third_dimension( self) result( field_third_dimension)
      class(atype_field_3D), intent(in) :: self
      type(type_third_dimension)        :: field_third_dimension
    end function get_third_dimension

    module function is_grid( self, grid) result( res)
      class(atype_field), intent(in) :: self
      class(*),           intent(in) :: grid
      logical                        :: res
    end function is_grid

    module function is_Arakawa_grid( self, field_Arakawa_grid) result( res)
      class(atype_field),      intent(in) :: self
      type(type_Arakawa_grid), intent(in) :: field_Arakawa_grid
      logical                             :: res
    end function is_Arakawa_grid

    module function is_pai( self, field_pai) result( res)
      class(atype_field),      intent(in) :: self
      type(type_par_arr_info), intent(in) :: field_pai
      logical                             :: res
    end function is_pai

    module function is_third_dimension( self, field_third_dimension) result( res)
      ! ..except for this one, which also works on type_field_2D - it just always returns .false. there
      class(atype_field),         intent(in) :: self
      type(type_third_dimension), intent(in) :: field_third_dimension
      logical                                :: res
    end function is_third_dimension

    ! ===== Reallocate

    module subroutine reallocate_logical_2D( self, mesh_new, d_nih)
      class(atype_field),                         intent(inout) :: self
      type(type_mesh),                            intent(in   ) :: mesh_new
      logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_logical_2D

    module subroutine reallocate_logical_3D( self, mesh_new, d_nih)
      class(atype_field),                           intent(inout) :: self
      type(type_mesh),                              intent(in   ) :: mesh_new
      logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_logical_3D

    module subroutine reallocate_int_2D( self, mesh_new, d_nih)
      class(atype_field),                         intent(inout) :: self
      type(type_mesh),                            intent(in   ) :: mesh_new
      integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_int_2D

    module subroutine reallocate_int_3D( self, mesh_new, d_nih)
      class(atype_field),                           intent(inout) :: self
      type(type_mesh),                              intent(in   ) :: mesh_new
      integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_int_3D

    module subroutine reallocate_dp_2D( self, mesh_new, d_nih)
      class(atype_field),                          intent(inout) :: self
      type(type_mesh),                             intent(in   ) :: mesh_new
      real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_dp_2D

    module subroutine reallocate_dp_3D( self, mesh_new, d_nih)
      class(atype_field),                            intent(inout) :: self
      type(type_mesh),                               intent(in   ) :: mesh_new
      real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_dp_3D

    ! ===== Remap

    module subroutine remap_dp_2D( self, mesh_new, d_nih)
      class(atype_field),                          intent(inout) :: self
      type(type_mesh),                             intent(in   ) :: mesh_new
      real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_dp_2D

    module subroutine remap_dp_3D( self, mesh_new, d_nih)
      class(atype_field),                            intent(inout) :: self
      type(type_mesh),                               intent(in   ) :: mesh_new
      real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_dp_3D

    ! ===== i/o

    module subroutine write_to_netcdf( self, filename, ncid)
      class(atype_field), intent(in) :: self
      character(len=*),   intent(in) :: filename
      integer,            intent(in) :: ncid
    end subroutine write_to_netcdf

    module subroutine read_from_netcdf( self, filename, ncid)
      class(atype_field), intent(inout) :: self
      character(len=*),   intent(in   ) :: filename
      integer,            intent(in   ) :: ncid
    end subroutine read_from_netcdf

  end interface

end module fields_basic