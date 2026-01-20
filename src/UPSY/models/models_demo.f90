module models_demo

  use precisions, only: dp
  use parameters, only: pi, NaN
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use models_basic, only: atype_model, atype_model_context
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_demo_model

  type, extends( atype_model) :: type_demo_model

    ! Some ice-model-esque data fields
    real(dp), dimension(:  ), contiguous, pointer :: H        => null()
    real(dp), dimension(:,:), contiguous, pointer :: u_3D     => null()
    real(dp), dimension(:,:), contiguous, pointer :: v_3D     => null()
    logical,  dimension(:  ), contiguous, pointer :: mask_ice => null()
    real(dp), dimension(:,:), contiguous, pointer :: T2m      => null()
    type(MPI_WIN) :: wH, wu_3D, wv_3D, wmask_ice, wT2m

  contains

    procedure, public :: create    => create_demo_model_abs
    procedure, public :: init      => init_demo_model_abs
    procedure, public :: run       => run_demo_model_abs
    procedure, public :: remap     => remap_demo_model_abs

    procedure, public :: create_ct => demo_model_context_create
    procedure, public :: init_ct   => demo_model_context_init
    procedure, public :: run_ct    => demo_model_context_run
    procedure, public :: remap_ct  => demo_model_context_remap

  end type type_demo_model

  type, extends( atype_model_context) :: type_demo_model_context_create
    ! Whatever data your model requires to be created
    type(type_mesh), pointer, public :: mesh
  end type type_demo_model_context_create

  type, extends( atype_model_context) :: type_demo_model_context_init
    ! Whatever data your model requires to initialise itself
    integer, public :: a
    integer, public :: b
  end type type_demo_model_context_init

  type, extends( atype_model_context) :: type_demo_model_context_run
    ! Whatever data your model requires to run
    integer, public :: c
    integer, public :: d
  end type type_demo_model_context_run

  type, extends( atype_model_context) :: type_demo_model_context_remap
    ! Whatever data your model requires to be remapped
    type(type_mesh), pointer, public :: mesh_new
  end type type_demo_model_context_remap

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  interface

    module subroutine create_demo_model_abs( self, context)
      class(type_demo_model),     intent(inout) :: self
      class(atype_model_context), intent(in   ) :: context
    end subroutine create_demo_model_abs

    module subroutine init_demo_model_abs( self, context)
      class(type_demo_model),     intent(inout) :: self
      class(atype_model_context), intent(in   ) :: context
    end subroutine init_demo_model_abs

    module subroutine run_demo_model_abs( self, context)
      class(type_demo_model),     intent(inout) :: self
      class(atype_model_context), intent(in   ) :: context
    end subroutine run_demo_model_abs

    module subroutine remap_demo_model_abs( self, context)
      class(type_demo_model),     intent(inout) :: self
      class(atype_model_context), intent(in   ) :: context
    end subroutine remap_demo_model_abs


    module function demo_model_context_create( self, mesh) result( context)
      class(type_demo_model),  intent(in   ) :: self
      type(type_mesh), target, intent(in   ) :: mesh
      type(type_demo_model_context_create)  :: context
    end function demo_model_context_create

    module function demo_model_context_init( self, a, b) result( context)
      class(type_demo_model),  intent(in   ) :: self
      integer,                 intent(in   ) :: a, b
      type(type_demo_model_context_init)     :: context
    end function demo_model_context_init

    module function demo_model_context_run( self, c, d) result( context)
      class(type_demo_model),  intent(in   ) :: self
      integer,                 intent(in   ) :: c, d
      type(type_demo_model_context_run)      :: context
    end function demo_model_context_run

    module function demo_model_context_remap( self, mesh_new) result( context)
      class(type_demo_model),  intent(in   ) :: self
      type(type_mesh), target, intent(in   ) :: mesh_new
      type(type_demo_model_context_remap)   :: context
    end function demo_model_context_remap

  end interface

contains

end module models_demo
