module SMB_basic
  !< The basic surface mass balance model

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mpi_f08, only: MPI_WIN
  use UPSY_main, only: atype_model, atype_model_context
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid

  implicit none

  private

  public :: atype_SMB_model, type_SMB_model_context_allocate, type_SMB_model_context_initialise, &
    type_SMB_model_context_run, type_SMB_model_context_remap

  type, abstract, extends(atype_model) :: atype_SMB_model
    !< The basic surface mass balance model

      ! Main data fields
      real(dp), dimension(:), contiguous, pointer :: SMB   !< [m] Net annual  SMB
      type(MPI_WIN) :: wSMB

      ! Timestepping
      real(dp) :: t_next

    contains

      private

      ! procedure( allocate_SMB_model_ifc), deferred :: allocate_SMB_model

      ! The common parts of allocate/initialise/run/remap, executed by all extensions of atype_SMB_model
      procedure, public :: allocate   => allocate_SMB_model_common_abs
      procedure, public :: initialise => initialise_SMB_model_common_abs
      procedure, public :: run        => run_SMB_model_common_abs
      procedure, public :: remap      => remap_SMB_model_common_abs

      procedure, public :: ct_allocate   => SMB_model_context_allocate
      procedure, public :: ct_initialise => SMB_model_context_initialise
      procedure, public :: ct_run        => SMB_model_context_run
      procedure, public :: ct_remap      => SMB_model_context_remap

      ! ! The specific parts of allocate/initialise/run/remap, executed only be each concrete type of SMB_model
      procedure(allocate_SMB_model_ifc),   deferred :: allocate_SMB_model
      procedure(initialise_SMB_model_ifc), deferred :: initialise_SMB_model
      procedure(run_SMB_model_ifc),        deferred :: run_SMB_model
      procedure(remap_SMB_model_ifc),      deferred :: remap_SMB_model

  end type atype_SMB_model

  type, extends( atype_model_context) :: type_SMB_model_context_allocate
    ! The set of variables required by any SMB model in order to be allocated
    type(type_mesh), pointer, public :: mesh
  end type type_SMB_model_context_allocate

  type, extends( atype_model_context) :: type_SMB_model_context_initialise
    ! The set of variables required by any SMB model in order to be initialised
    character(len=3)  :: region_name
  end type type_SMB_model_context_initialise

  type, extends( atype_model_context) :: type_SMB_model_context_run
    ! The set of variables required by any SMB model in order to be run
    class(*),        pointer :: ice
    class(*),        pointer :: climate
    type(type_grid), pointer :: grid_smooth
    real(dp)                 :: time
    character(len=3)         :: region_name
  end type type_SMB_model_context_run

  type, extends( atype_model_context) :: type_SMB_model_context_remap
    ! The set of variables required by any SMB model in order to be remapped
    type(type_mesh), pointer, public :: mesh_new
  end type type_SMB_model_context_remap

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine allocate_SMB_model_ifc( self, context)
      import atype_SMB_model, type_SMB_model_context_allocate
      class(atype_SMB_model),                intent(inout) :: self
      type(type_SMB_model_context_allocate), intent(in   ) :: context
    end subroutine allocate_SMB_model_ifc

    subroutine initialise_SMB_model_ifc( self, context)
      import atype_SMB_model, type_SMB_model_context_initialise
      class(atype_SMB_model),                intent(inout) :: self
      type(type_SMB_model_context_initialise), intent(in   ) :: context
    end subroutine initialise_SMB_model_ifc

    subroutine run_SMB_model_ifc( self, context)
      import atype_SMB_model, type_SMB_model_context_run
      class(atype_SMB_model),           intent(inout) :: self
      type(type_SMB_model_context_run), intent(in   ) :: context
    end subroutine run_SMB_model_ifc

    subroutine remap_SMB_model_ifc( self, context)
      import atype_SMB_model, type_SMB_model_context_remap
      class(atype_SMB_model),             intent(inout) :: self
      type(type_SMB_model_context_remap), intent(in   ) :: context
    end subroutine remap_SMB_model_ifc

  end interface

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  interface

    module subroutine allocate_SMB_model_common_abs( self, context)
      class(atype_SMB_model),     intent(inout) :: self
      class(atype_model_context), intent(in   ) :: context
    end subroutine allocate_SMB_model_common_abs

    module subroutine initialise_SMB_model_common_abs( self, context)
      class(atype_SMB_model),     intent(inout) :: self
      class(atype_model_context), intent(in   ) :: context
    end subroutine initialise_SMB_model_common_abs

    module subroutine run_SMB_model_common_abs( self, context)
      class(atype_SMB_model),     intent(inout) :: self
      class(atype_model_context), intent(in   ) :: context
    end subroutine run_SMB_model_common_abs

    module subroutine remap_SMB_model_common_abs( self, context)
      class(atype_SMB_model),     intent(inout) :: self
      class(atype_model_context), intent(in   ) :: context
    end subroutine remap_SMB_model_common_abs


    module function SMB_model_context_allocate( self, mesh) result( context)
      class(atype_SMB_model),  intent(in   ) :: self
      type(type_mesh), target, intent(in   ) :: mesh
      type(type_SMB_model_context_allocate)  :: context
    end function SMB_model_context_allocate

    module function SMB_model_context_initialise( self) result( context)
      class(atype_SMB_model),   intent(in   ) :: self
      type(type_SMB_model_context_initialise) :: context
    end function SMB_model_context_initialise

    module function SMB_model_context_run( self) result( context)
      class(atype_SMB_model), intent(in   ) :: self
      type(type_SMB_model_context_run)      :: context
    end function SMB_model_context_run

    module function SMB_model_context_remap( self, mesh_new) result( context)
      class(atype_SMB_model),  intent(in   ) :: self
      type(type_mesh), target, intent(in   ) :: mesh_new
      type(type_SMB_model_context_remap)     :: context
    end function SMB_model_context_remap

  end interface

end module SMB_basic