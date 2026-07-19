module models_basic

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use fields_main, only: type_fields_registry
  use Arakawa_grid_mod, only: type_Arakawa_grid
  use fields_dimensions, only: type_third_dimension
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_model, &
    atype_model_context_remap

  ! Abstract basic model type
  ! =========================

  type, abstract :: atype_model

      ! Metadata
      character(len=:), allocatable, private :: name_val
      character(len=3),              private :: region_name_val

      ! Mesh
      type(type_mesh), pointer :: mesh

      ! Fields registry
      type(type_fields_registry), private :: flds_reg

    contains

      ! These routines all consist of two parts: a 'common' part that is executed for
      ! all models inheriting from atype_model, and a 'specific' part that is
      ! only executed for each specific model class. The specific parts are defined
      ! in the deferred procedures 'allocate_model', 'initialise_model', etc.

      procedure, public :: allocate_model
      procedure, public :: deallocate_model
      procedure, public :: initialise_model
      procedure, public :: run_model
      procedure, public :: remap

      procedure(remap_model_ifc),      deferred :: remap_model

      ! i/o

      procedure, public :: write_to_restart_file
      procedure, public :: read_from_restart_file

      ! Memory management for fields

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

      ! ===== Basics

      generic,   public  :: operator(==) => eq
      procedure, private :: eq => test_model_equality

      ! Metadata
      procedure, public :: set_name
      procedure, public :: name => get_name
      procedure, public :: is_name

      procedure, public :: set_region_name
      procedure, public :: region_name => get_region_name
      procedure, public :: is_region_name

  end type atype_model

  ! Context classes for initialise/run/remap
  ! =================================================

  type, abstract :: atype_model_context_remap
    type(type_mesh), pointer :: mesh_new
  end type atype_model_context_remap

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine remap_model_ifc( self, context)
      import atype_model, atype_model_context_remap
      class(atype_model),                       intent(inout) :: self
      class(atype_model_context_remap), target, intent(in   ) :: context
    end subroutine remap_model_ifc

  end interface

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  interface

    module subroutine remap( self, context)
      class(atype_model),                       intent(inout) :: self
      class(atype_model_context_remap), target, intent(in   ) :: context
    end subroutine remap

  end interface

  ! create_field
  interface

    module subroutine create_field_logical_2D( self, d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units, remap_method)
      class(atype_model),                         intent(inout) :: self
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
      class(atype_model),                           intent(inout) :: self
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
      class(atype_model),                         intent(inout) :: self
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
      class(atype_model),                           intent(inout) :: self
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
      class(atype_model),                          intent(inout) :: self
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
      class(atype_model),                            intent(inout) :: self
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
      class(atype_model),                         intent(inout) :: self
      character(len=*),                           intent(in   ) :: field_name
      type(type_mesh),                            intent(in   ) :: mesh_new
      logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_logical_2D

    module subroutine reallocate_field_logical_3D( self, mesh_new, field_name, d_nih)
      class(atype_model),                           intent(inout) :: self
      character(len=*),                             intent(in   ) :: field_name
      type(type_mesh),                              intent(in   ) :: mesh_new
      logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_logical_3D

    module subroutine reallocate_field_int_2D( self, mesh_new, field_name, d_nih)
      class(atype_model),                         intent(inout) :: self
      character(len=*),                           intent(in   ) :: field_name
      type(type_mesh),                            intent(in   ) :: mesh_new
      integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_int_2D

    module subroutine reallocate_field_int_3D( self, mesh_new, field_name, d_nih)
      class(atype_model),                           intent(inout) :: self
      character(len=*),                             intent(in   ) :: field_name
      type(type_mesh),                              intent(in   ) :: mesh_new
      integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_int_3D

    module subroutine reallocate_field_dp_2D( self, mesh_new, field_name, d_nih)
      class(atype_model),                          intent(inout) :: self
      character(len=*),                            intent(in   ) :: field_name
      type(type_mesh),                             intent(in   ) :: mesh_new
      real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_dp_2D

    module subroutine reallocate_field_dp_3D( self, mesh_new, field_name, d_nih)
      class(atype_model),                            intent(inout) :: self
      character(len=*),                              intent(in   ) :: field_name
      type(type_mesh),                               intent(in   ) :: mesh_new
      real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_dp_3D

  end interface

  ! remap
  interface

    module subroutine remap_field_logical_2D( self, mesh_new, field_name, d_nih)
      class(atype_model),                         intent(inout) :: self
      character(len=*),                           intent(in   ) :: field_name
      type(type_mesh),                            intent(in   ) :: mesh_new
      logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_field_logical_2D

    module subroutine remap_field_logical_3D( self, mesh_new, field_name, d_nih)
      class(atype_model),                           intent(inout) :: self
      character(len=*),                             intent(in   ) :: field_name
      type(type_mesh),                              intent(in   ) :: mesh_new
      logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_field_logical_3D

    module subroutine remap_field_int_2D( self, mesh_new, field_name, d_nih)
      class(atype_model),                         intent(inout) :: self
      character(len=*),                           intent(in   ) :: field_name
      type(type_mesh),                            intent(in   ) :: mesh_new
      integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_field_int_2D

    module subroutine remap_field_int_3D( self, mesh_new, field_name, d_nih)
      class(atype_model),                           intent(inout) :: self
      character(len=*),                             intent(in   ) :: field_name
      type(type_mesh),                              intent(in   ) :: mesh_new
      integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_field_int_3D

    module subroutine remap_field_dp_2D( self, mesh_new, field_name, d_nih)
      class(atype_model),                          intent(inout) :: self
      character(len=*),                            intent(in   ) :: field_name
      type(type_mesh),                             intent(in   ) :: mesh_new
      real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_field_dp_2D

    module subroutine remap_field_dp_3D( self, mesh_new, field_name, d_nih)
      class(atype_model),                            intent(inout) :: self
      character(len=*),                              intent(in   ) :: field_name
      type(type_mesh),                               intent(in   ) :: mesh_new
      real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine remap_field_dp_3D

  end interface

  ! NetCDF i/o
  interface

    module subroutine write_to_restart_file( self, output_dir, filename)
      class(atype_model),                  intent(in   ) :: self
      character(len=*),                    intent(in   ) :: output_dir
      character(:), allocatable, optional, intent(  out) :: filename
    end subroutine write_to_restart_file

    module subroutine read_from_restart_file( self, filename)
      class(atype_model), intent(inout) :: self
      character(len=*),   intent(in   ) :: filename
    end subroutine read_from_restart_file

  end interface

  ! Basics
  interface

    module function test_model_equality( model1, model2) result( res)
      class(atype_model), intent(in) :: model1, model2
      logical                        :: res
    end function test_model_equality

    ! ===== Set/get functions

    ! Metadata

    module subroutine set_name( self, name)
      class(atype_model), intent(inout) :: self
      character(len=*),   intent(in   ) :: name
    end subroutine set_name

    module function get_name( self) result( name)
      class(atype_model), intent(in) :: self
      character(:), allocatable      :: name
    end function get_name

    module function is_name( self, name) result( res)
      class(atype_model), intent(in) :: self
      character(len=*),   intent(in) :: name
      logical                        :: res
    end function is_name

    module subroutine set_region_name( self, region_name)
      class(atype_model), intent(inout) :: self
      character(len=*),   intent(in   ) :: region_name
    end subroutine set_region_name

    module function get_region_name( self) result( region_name)
      class(atype_model), intent(in) :: self
      character(:), allocatable      :: region_name
    end function get_region_name

    module function is_region_name( self, region_name) result( res)
      class(atype_model), intent(in) :: self
      character(len=*),   intent(in) :: region_name
      logical                        :: res
    end function is_region_name

  end interface

contains

  subroutine allocate_model( self, name, region_name, mesh)
    !< Allocate stuff that is common to all models
    !< (call this from your model-specific allocate routine)

    ! In/output variables:
    class(atype_model),      intent(inout) :: self
    character(len=*),        intent(in   ) :: name
    character(len=*),        intent(in   ) :: region_name
    type(type_mesh), target, intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Set model metadata and mesh
    call self%set_name       ( name)
    call self%set_region_name( region_name)
    self%mesh => mesh

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_model

  subroutine deallocate_model( self)
    !< Deallocate stuff that is common to all models
    !< (call this from your model-specific deallocate routine)

    ! In/output variables:
    class(atype_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Un-set model metadata and mesh
    call self%set_name( 'empty_model')
    call self%set_region_name( '!!!')
    nullify( self%mesh)

    ! Deallocate shared memory for all the fields
    call self%flds_reg%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_model

  subroutine initialise_model( self)
    !< Initialise stuff that is common to all models
    !< (call this from your model-specific initialise routine)

    ! In/output variables:
    class(atype_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Right now, there's nothing that can be put here yet,
    ! but it's useful to have the allocate/deallocate/initialise/run/remap
    ! routines all follow the same general layout.

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_model

  subroutine run_model( self)
    !< Run stuff that is common to all models
    !< (call this from your model-specific run routine)

    ! In/output variables:
    class(atype_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Right now, there's nothing that can be put here yet,
    ! but it's useful to have the allocate/deallocate/initialise/run/remap
    ! routines all follow the same general layout.

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_model

end module models_basic
