module fields_registry

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
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

    generic,   public  :: operator(==) => eq
    procedure, private :: eq => test_fields_registry_equality
    procedure, private :: add => add_field_to_registry
    procedure, private :: extend => extend_field_registry
    procedure, public  :: find => find_field_by_name
    procedure, public  :: print_info

    procedure, public  :: destroy

    procedure, public  :: write_to_netcdf
    procedure, public  :: read_from_netcdf

  end type type_fields_registry

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  ! basics
  interface

    module function test_fields_registry_equality( flds_reg1, flds_reg2) result( res)
      class(type_fields_registry), intent(in) :: flds_reg1, flds_reg2
      logical                                 :: res
    end function test_fields_registry_equality

    module subroutine add_field_to_registry( flds_reg, field)
      class(type_fields_registry), intent(inout) :: flds_reg
      class(atype_field),          intent(in   ) :: field
    end subroutine add_field_to_registry

    module subroutine extend_field_registry( flds_reg)
      class(type_fields_registry), intent(inout) :: flds_reg
    end subroutine extend_field_registry

    module subroutine print_info( flds_reg)
      class(type_fields_registry), intent(in) :: flds_reg
    end subroutine print_info

    module subroutine destroy( flds_reg)
      class(type_fields_registry), intent(inout) :: flds_reg
    end subroutine destroy

    module function find_field_by_name( flds_reg, name) result(i)
      class(type_fields_registry), intent(in) :: flds_reg
      character(len=*),            intent(in) :: name
      integer                                 :: i
    end function find_field_by_name

  end interface

  ! create_field
  interface

    module subroutine create_field_logical_2D( flds_reg, d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units)

      class(type_fields_registry),                intent(inout) :: flds_reg
      logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                              intent(inout) :: w
      class(*), target,                           intent(in   ) :: field_grid
      type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
      character(len=*),                           intent(in   ) :: name
      character(len=*),                           intent(in   ) :: long_name
      character(len=*),                           intent(in   ) :: units

    end subroutine create_field_logical_2D

    module subroutine create_field_logical_3D( flds_reg, d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units)

      class(type_fields_registry),                  intent(inout) :: flds_reg
      logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                                intent(inout) :: w
      class(*), target,                             intent(in   ) :: field_grid
      type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
      type(type_third_dimension),                   intent(in   ) :: field_third_dimension
      character(len=*),                             intent(in   ) :: name
      character(len=*),                             intent(in   ) :: long_name
      character(len=*),                             intent(in   ) :: units

    end subroutine create_field_logical_3D

    module subroutine create_field_int_2D( flds_reg, d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units)

      class(type_fields_registry),                intent(inout) :: flds_reg
      integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                              intent(inout) :: w
      class(*), target,                           intent(in   ) :: field_grid
      type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
      character(len=*),                           intent(in   ) :: name
      character(len=*),                           intent(in   ) :: long_name
      character(len=*),                           intent(in   ) :: units

    end subroutine create_field_int_2D

    module subroutine create_field_int_3D( flds_reg, d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units)

      class(type_fields_registry),                  intent(inout) :: flds_reg
      integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                                intent(inout) :: w
      class(*), target,                             intent(in   ) :: field_grid
      type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
      type(type_third_dimension),                   intent(in   ) :: field_third_dimension
      character(len=*),                             intent(in   ) :: name
      character(len=*),                             intent(in   ) :: long_name
      character(len=*),                             intent(in   ) :: units

    end subroutine create_field_int_3D

    module subroutine create_field_dp_2D( flds_reg, d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units)

      class(type_fields_registry),                 intent(inout) :: flds_reg
      real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                               intent(inout) :: w
      class(*), target,                            intent(in   ) :: field_grid
      type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
      character(len=*),                            intent(in   ) :: name
      character(len=*),                            intent(in   ) :: long_name
      character(len=*),                            intent(in   ) :: units

    end subroutine create_field_dp_2D

    module subroutine create_field_dp_3D( flds_reg, d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units)

      class(type_fields_registry),                   intent(inout) :: flds_reg
      real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                                 intent(inout) :: w
      class(*), target,                              intent(in   ) :: field_grid
      type(type_Arakawa_grid),                       intent(in   ) :: field_Arakawa_grid
      type(type_third_dimension),                    intent(in   ) :: field_third_dimension
      character(len=*),                              intent(in   ) :: name
      character(len=*),                              intent(in   ) :: long_name
      character(len=*),                              intent(in   ) :: units

    end subroutine create_field_dp_3D

  end interface

  ! i/o
  interface

    module subroutine write_to_netcdf( flds_reg, filename, ncid)
      class(type_fields_registry), intent(in) :: flds_reg
      character(len=*),            intent(in) :: filename
      integer,                     intent(in) :: ncid
    end subroutine write_to_netcdf

    module subroutine read_from_netcdf( flds_reg, filename, ncid)
      class(type_fields_registry), intent(inout) :: flds_reg
      character(len=*),            intent(in   ) :: filename
      integer,                     intent(in   ) :: ncid
    end subroutine read_from_netcdf

  end interface

contains

end module fields_registry