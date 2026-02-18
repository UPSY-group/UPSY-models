module fields_registry

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use fields_basic, only: atype_field, &
    type_field_logical_2D, type_field_int_2D, type_field_dp_2D, &
    type_field_logical_3D, type_field_int_3D, type_field_dp_3D
  use Arakawa_grid_mod, only: type_Arakawa_grid, Arakawa_grid
  use fields_dimensions, only: type_third_dimension
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use deallocate_dist_shared_mod, only: deallocate_dist_shared
  use fields_init_field, only: initialise_field
  use mpi_f08, only: MPI_WIN, MPI_GATHER, MPI_INTEGER, MPI_COMM_WORLD

  implicit none

  private

  public :: type_fields_registry

  ! Wrapper so we can have a mixed-type array
  ! =========================================

  type :: type_field_box
    class(atype_field), allocatable :: p
  end type type_field_box

  ! Field collection
  ! ================

  type :: type_fields_registry

    integer,                           public :: n     = 0
    integer,                           public :: n_max = 0
    type(type_field_box), allocatable, public :: items(:)

    contains

    generic,   public  :: create_field => &
      create_field_logical_2D, create_field_int_2D, create_field_dp_2D, &
      create_field_logical_3D, create_field_int_3D, create_field_dp_3D
    procedure, private :: create_field_logical_2D
    procedure, private :: create_field_int_2D
    procedure, private :: create_field_dp_2D
    procedure, private :: create_field_logical_3D
    procedure, private :: create_field_int_3D
    procedure, private :: create_field_dp_3D

    generic,   public  :: reallocate_field => &
      reallocate_field_logical_2D, &
      reallocate_field_logical_3D, &
      reallocate_field_int_2D, &
      reallocate_field_int_3D, &
      reallocate_field_dp_2D, &
      reallocate_field_dp_3D
    procedure, private :: reallocate_field_logical_2D
    procedure, private :: reallocate_field_logical_3D
    procedure, private :: reallocate_field_int_2D
    procedure, private :: reallocate_field_int_3D
    procedure, private :: reallocate_field_dp_2D
    procedure, private :: reallocate_field_dp_3D

    generic,   public  :: remap_field => &
      remap_field_logical_2D, &
      remap_field_logical_3D, &
      remap_field_int_2D, &
      remap_field_int_3D, &
      remap_field_dp_2D, &
      remap_field_dp_3D
    procedure, private :: remap_field_logical_2D
    procedure, private :: remap_field_logical_3D
    procedure, private :: remap_field_int_2D
    procedure, private :: remap_field_int_3D
    procedure, private :: remap_field_dp_2D
    procedure, private :: remap_field_dp_3D

    procedure, public  :: write_to_netcdf
    procedure, public  :: read_from_netcdf

    generic,   public  :: operator(==) => eq
    procedure, private :: eq => test_fields_registry_equality
    procedure, private :: add => add_field_to_registry
    procedure, private :: extend => extend_field_registry
    procedure, public  :: find => find_field_by_name
    procedure, public  :: print_info
    procedure, public  :: destroy

  end type type_fields_registry

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  ! create_field
  interface

    module subroutine create_field_logical_2D( self, d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units, remap_method)

      class(type_fields_registry),                intent(inout) :: self
      logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                              intent(inout) :: w
      class(*), target,                           intent(in   ) :: field_grid
      type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
      character(len=*),                 optional, intent(in   ) :: name
      character(len=*),                 optional, intent(in   ) :: long_name
      character(len=*),                 optional, intent(in   ) :: units
      character(len=*),                 optional, intent(in   ) :: remap_method

    end subroutine create_field_logical_2D

    module subroutine create_field_logical_3D( self, d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)

      class(type_fields_registry),                  intent(inout) :: self
      logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                                intent(inout) :: w
      class(*), target,                             intent(in   ) :: field_grid
      type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
      type(type_third_dimension),                   intent(in   ) :: field_third_dimension
      character(len=*),                   optional, intent(in   ) :: name
      character(len=*),                   optional, intent(in   ) :: long_name
      character(len=*),                   optional, intent(in   ) :: units
      character(len=*),                   optional, intent(in   ) :: remap_method

    end subroutine create_field_logical_3D

    module subroutine create_field_int_2D( self, d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units, remap_method)

      class(type_fields_registry),                intent(inout) :: self
      integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                              intent(inout) :: w
      class(*), target,                           intent(in   ) :: field_grid
      type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
      character(len=*),                 optional, intent(in   ) :: name
      character(len=*),                 optional, intent(in   ) :: long_name
      character(len=*),                 optional, intent(in   ) :: units
      character(len=*),                 optional, intent(in   ) :: remap_method

    end subroutine create_field_int_2D

    module subroutine create_field_int_3D( self, d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)

      class(type_fields_registry),                  intent(inout) :: self
      integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                                intent(inout) :: w
      class(*), target,                             intent(in   ) :: field_grid
      type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
      type(type_third_dimension),                   intent(in   ) :: field_third_dimension
      character(len=*),                   optional, intent(in   ) :: name
      character(len=*),                   optional, intent(in   ) :: long_name
      character(len=*),                   optional, intent(in   ) :: units
      character(len=*),                   optional, intent(in   ) :: remap_method

    end subroutine create_field_int_3D

    module subroutine create_field_dp_2D( self, d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units, remap_method)

      class(type_fields_registry),                 intent(inout) :: self
      real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                               intent(inout) :: w
      class(*), target,                            intent(in   ) :: field_grid
      type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
      character(len=*),                  optional, intent(in   ) :: name
      character(len=*),                  optional, intent(in   ) :: long_name
      character(len=*),                  optional, intent(in   ) :: units
      character(len=*),                  optional, intent(in   ) :: remap_method

    end subroutine create_field_dp_2D

    module subroutine create_field_dp_3D( self, d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)

      class(type_fields_registry),                   intent(inout) :: self
      real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                                 intent(inout) :: w
      class(*), target,                              intent(in   ) :: field_grid
      type(type_Arakawa_grid),                       intent(in   ) :: field_Arakawa_grid
      type(type_third_dimension),                    intent(in   ) :: field_third_dimension
      character(len=*),                    optional, intent(in   ) :: name
      character(len=*),                    optional, intent(in   ) :: long_name
      character(len=*),                    optional, intent(in   ) :: units
      character(len=*),                    optional, intent(in   ) :: remap_method

    end subroutine create_field_dp_3D

  end interface

  ! reallocate
  interface

    module subroutine reallocate_field_logical_2D( self, mesh_new, field_name, d_nih)
      class(type_fields_registry),                intent(inout) :: self
      character(len=*),                           intent(in   ) :: field_name
      type(type_mesh),                            intent(in   ) :: mesh_new
      logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_logical_2D

    module subroutine reallocate_field_logical_3D( self, mesh_new, field_name, d_nih)
      class(type_fields_registry),                  intent(inout) :: self
      character(len=*),                             intent(in   ) :: field_name
      type(type_mesh),                              intent(in   ) :: mesh_new
      logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_logical_3D

    module subroutine reallocate_field_int_2D( self, mesh_new, field_name, d_nih)
      class(type_fields_registry),                intent(inout) :: self
      character(len=*),                           intent(in   ) :: field_name
      type(type_mesh),                            intent(in   ) :: mesh_new
      integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_int_2D

    module subroutine reallocate_field_int_3D( self, mesh_new, field_name, d_nih)
      class(type_fields_registry),                  intent(inout) :: self
      character(len=*),                             intent(in   ) :: field_name
      type(type_mesh),                              intent(in   ) :: mesh_new
      integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_int_3D

    module subroutine reallocate_field_dp_2D( self, mesh_new, field_name, d_nih)
      class(type_fields_registry),                 intent(inout) :: self
      character(len=*),                            intent(in   ) :: field_name
      type(type_mesh),                             intent(in   ) :: mesh_new
      real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_dp_2D

    module subroutine reallocate_field_dp_3D( self, mesh_new, field_name, d_nih)
      class(type_fields_registry),                   intent(inout) :: self
      character(len=*),                              intent(in   ) :: field_name
      type(type_mesh),                               intent(in   ) :: mesh_new
      real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_dp_3D

  end interface

  ! remap
  interface

    module subroutine remap_field_logical_2D( self, mesh_new, field_name, d_nih)
      class(type_fields_registry),                intent(inout) :: self
      character(len=*),                           intent(in   ) :: field_name
      type(type_mesh),                            intent(in   ) :: mesh_new
      logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_field_logical_2D

    module subroutine remap_field_logical_3D( self, mesh_new, field_name, d_nih)
      class(type_fields_registry),                  intent(inout) :: self
      character(len=*),                             intent(in   ) :: field_name
      type(type_mesh),                              intent(in   ) :: mesh_new
      logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_field_logical_3D

    module subroutine remap_field_int_2D( self, mesh_new, field_name, d_nih)
      class(type_fields_registry),                intent(inout) :: self
      character(len=*),                           intent(in   ) :: field_name
      type(type_mesh),                            intent(in   ) :: mesh_new
      integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_field_int_2D

    module subroutine remap_field_int_3D( self, mesh_new, field_name, d_nih)
      class(type_fields_registry),                  intent(inout) :: self
      character(len=*),                             intent(in   ) :: field_name
      type(type_mesh),                              intent(in   ) :: mesh_new
      integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_field_int_3D

    module subroutine remap_field_dp_2D( self, mesh_new, field_name, d_nih)
      class(type_fields_registry),                 intent(inout) :: self
      character(len=*),                            intent(in   ) :: field_name
      type(type_mesh),                             intent(in   ) :: mesh_new
      real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_field_dp_2D

    module subroutine remap_field_dp_3D( self, mesh_new, field_name, d_nih)
      class(type_fields_registry),                   intent(inout) :: self
      character(len=*),                              intent(in   ) :: field_name
      type(type_mesh),                               intent(in   ) :: mesh_new
      real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_field_dp_3D

  end interface

  ! i/o
  interface

    module subroutine write_to_netcdf( self, filename, ncid)
      class(type_fields_registry), intent(in) :: self
      character(len=*),            intent(in) :: filename
      integer,                     intent(in) :: ncid
    end subroutine write_to_netcdf

    module subroutine read_from_netcdf( self, filename, ncid)
      class(type_fields_registry), intent(inout) :: self
      character(len=*),            intent(in   ) :: filename
      integer,                     intent(in   ) :: ncid
    end subroutine read_from_netcdf

  end interface

  ! basics
  interface

    module function test_fields_registry_equality( flds_reg1, flds_reg2) result( res)
      class(type_fields_registry), intent(in) :: flds_reg1, flds_reg2
      logical                                 :: res
    end function test_fields_registry_equality

    module subroutine add_field_to_registry( self, field)
      class(type_fields_registry), intent(inout) :: self
      class(atype_field),          intent(in   ) :: field
    end subroutine add_field_to_registry

    module subroutine extend_field_registry( self)
      class(type_fields_registry), intent(inout) :: self
    end subroutine extend_field_registry

    module function find_field_by_name( self, name) result(i)
      class(type_fields_registry), intent(in) :: self
      character(len=*),            intent(in) :: name
      integer                                 :: i
    end function find_field_by_name

    module subroutine print_info( self)
      class(type_fields_registry), intent(in) :: self
    end subroutine print_info

    module subroutine destroy( self)
      class(type_fields_registry), intent(inout) :: self
    end subroutine destroy

  end interface

end module fields_registry