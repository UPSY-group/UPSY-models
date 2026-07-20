module calving_model_basic

  use precisions, only: dp
  use parameters, only: pi, NaN
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use calving_model_data, only: atype_calving_model_data
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_calving_model

  type, abstract, extends(atype_calving_model_data) :: atype_calving_model
    !< Stuff that is common to all calving models
    !<
    !< (except for the variables that we want other models to
    !< be able to access, which are already defined in atype_calving_model_data)

    real(dp) :: t_next   !< Time when the calving model should be run next

    contains

      ! Type-bound procedures that apply to all demo models
      procedure, public :: allocate_calving_model
      procedure, public :: deallocate_calving_model
      procedure, public :: initialise_calving_model
      procedure, public :: run_calving_model
      procedure, public :: remap_calving_model

      ! Deferred procedures that must be defined by each individual demo model
      procedure(calving_model_allocate_ifc),   deferred :: allocate
      procedure(calving_model_deallocate_ifc), deferred :: deallocate
      procedure(calving_model_initialise_ifc), deferred :: initialise
      procedure(calving_model_run_ifc),        deferred :: run
      procedure(calving_model_remap_ifc),      deferred :: remap

  end type atype_calving_model

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine calving_model_allocate_ifc( self, name, region_name, mesh)
      import atype_calving_model, type_mesh
      class(atype_calving_model), intent(inout) :: self
      character(len=*),           intent(in   ) :: name
      character(len=*),           intent(in   ) :: region_name
      type(type_mesh), target,    intent(in   ) :: mesh
    end subroutine calving_model_allocate_ifc

    subroutine calving_model_deallocate_ifc( self)
      import atype_calving_model
      class(atype_calving_model), intent(inout) :: self
    end subroutine calving_model_deallocate_ifc

    subroutine calving_model_initialise_ifc( self)
      import atype_calving_model
      class(atype_calving_model), intent(inout) :: self
    end subroutine calving_model_initialise_ifc

    subroutine calving_model_run_ifc( self)
      import atype_calving_model
      class(atype_calving_model), intent(inout) :: self
    end subroutine calving_model_run_ifc

    subroutine calving_model_remap_ifc( self, mesh_new)
      import atype_calving_model, type_mesh
      class(atype_calving_model), intent(inout) :: self
      type(type_mesh), target,    intent(in   ) :: mesh_new
    end subroutine calving_model_remap_ifc

  end interface

contains

  subroutine allocate_calving_model( self, name, region_name, mesh, nz)
    !< Allocate stuff that is common to all calving models
    !< (call this from your calving model-specific allocate routine)

    ! In/output variables:
    class(atype_calving_model), intent(inout) :: self
    character(len=*),           intent(in   ) :: name
    character(len=*),           intent(in   ) :: region_name
    type(type_mesh), target,    intent(in   ) :: mesh
    integer,                    intent(in   ) :: nz

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_calving_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate stuff that is common to all models
    call self%allocate_model( name, region_name, mesh)

    ! Allocate stuff that is specific to calving models

    call self%create_field( self%Hi_calved, self%wHi_calved, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Hi_calved', &
      long_name = 'Ice thickness after calving', &
      units     = 'm', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_calving_model

  subroutine deallocate_calving_model( self)
    !< Deallocate stuff that is common to all calving models
    !< (call this from your calving model-specific deallocate routine)

    ! In/output variables:
    class(atype_calving_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_calving_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate stuff that is common to all models
    call self%deallocate_model()

    ! Deallocate stuff that is specific to calving models

    nullify( self%Hi_calved)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_calving_model

  subroutine initialise_calving_model( self)
    !< Initialise stuff that is common to all calving models
    !< (call this from your calving model-specific initiaise routine)

    ! In/output variables:
    class(atype_calving_model), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_calving_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise stuff that is common to all models
    call self%initialise_model()

    ! Initialise stuff that is specific to calving models

    ! Set time of next calculation to start time
    self%t_next = C%start_time_of_run

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_calving_model

  subroutine run_calving_model( self)
    !< Run stuff that is common to all calving models
    !< (call this from your calving model-specific run routine)

    ! In/output variables:
    class(atype_calving_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_calving_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run stuff that is common to all models
    call self%run_model()

    ! Run stuff that is specific to calving models

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_calving_model

  subroutine remap_calving_model( self, mesh_new)
    !< Remap stuff that is common to all calving models
    !< (call this from your calving model-specific remap routine)

    ! In/output variables:
    class(atype_calving_model), intent(inout) :: self
    type(type_mesh), target, intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_calving_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap stuff that is common to all models
    call self%remap_model( mesh_new)

    ! Remap stuff that is specific to calving models

    ! DENK DROM
    call crash('remapping calving models is not yet supported')
    ! think about what should happen with Hi_calved

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_calving_model

end module calving_model_basic
