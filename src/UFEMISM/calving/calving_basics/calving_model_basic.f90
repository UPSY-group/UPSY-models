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
      procedure, public :: allocate   => calving_model_allocate
      procedure, public :: deallocate => calving_model_deallocate
      procedure, public :: initialise => calving_model_initialise
      procedure, public :: run        => calving_model_run
      procedure, public :: remap      => calving_model_remap

      ! Deferred procedures that must be defined by each individual demo model
      procedure(calving_model_allocate_ifc),   deferred :: allocate_calving_model
      procedure(calving_model_deallocate_ifc), deferred :: deallocate_calving_model
      procedure(calving_model_initialise_ifc), deferred :: initialise_calving_model
      procedure(calving_model_run_ifc),        deferred :: run_calving_model
      procedure(calving_model_remap_ifc),      deferred :: remap_calving_model

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

  subroutine calving_model_allocate( self, name, region_name, mesh)

    ! In/output variables:
    class(atype_calving_model), intent(inout) :: self
    character(len=*),           intent(in   ) :: name
    character(len=*),           intent(in   ) :: region_name
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calving_model_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate stuff that is common to all models
    call self%allocate_model( name, region_name, mesh)

    ! Allocate stuff that is common to all calving models
    call self%create_field( self%Hi_calved, self%wHi_calved, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Hi_calved', &
      long_name = 'Ice thickness after calving', &
      units     = 'm', &
      remap_method = 'reallocate')

    ! Allocate stuff that is specific to each individual calving model implementation
    call self%allocate_calving_model( name, region_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calving_model_allocate

  subroutine calving_model_deallocate( self)

    ! In/output variables:
    class(atype_calving_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calving_model_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate stuff that is common to all models
    call self%deallocate_model()

    ! Deallocate stuff that is common to all calving models
    nullify( self%Hi_calved)

    ! Deallocate stuff that is specific to each individual calving model implementation
    call self%deallocate_calving_model()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calving_model_deallocate

  subroutine calving_model_initialise( self)

    ! In/output variables:
    class(atype_calving_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calving_model_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise stuff that is common to all models
    call self%initialise_model()

    ! Initialise stuff that is common to all calving models

    ! Set time of next calculation to start time
    self%t_next = C%start_time_of_run

    ! Initialise stuff that is specific to each individual calving model implementation
    call self%initialise_calving_model()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calving_model_initialise

  subroutine calving_model_run( self)

    ! In/output variables:
    class(atype_calving_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calving_model_run'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run stuff that is common to all models
    call self%run_model()

    ! Run stuff that is common to all calving models

    ! Run stuff that is specific to each individual calving model implementation
    call self%run_calving_model()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calving_model_run

  subroutine calving_model_remap( self, mesh_new)

    ! In/output variables:
    class(atype_calving_model), intent(inout) :: self
    type(type_mesh), target,    intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calving_model_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap stuff that is common to all models
    call self%remap_model( mesh_new)

    ! Remap stuff that is common to all calving models

    ! DENK DROM
    call crash('remapping calving models is not yet supported')
    ! think about what should happen with Hi_calved

    ! Remap stuff that is specific to each individual calving model implementation
    call self%remap_calving_model( mesh_new)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calving_model_remap

end module calving_model_basic
