module SMB_model_basic

  use precisions, only: dp
  use parameters, only: pi, NaN
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use models_basic, only: atype_model, atype_model_context_remap
  use SMB_model_data, only: atype_SMB_model_data
  use mpi_f08, only: MPI_WIN
  use ice_model_types, only: type_ice_model
  use climate_model_types, only: type_climate_model
  use grid_types, only: type_grid
  use reference_geometry_types, only: type_reference_geometry

  implicit none

  private

  public :: atype_SMB_model, type_SMB_model_context_remap

  type, abstract, extends(atype_SMB_model_data) :: atype_SMB_model
    !< Stuff that is common to all SMB models
    !<
    !< (except for the variables that we want other models to
    !< be able to access, which are already defined in atype_SMB_model_data)

    real(dp) :: t_next   !< Time when the SMB model should be run next

    contains

      ! These routines all consist of two parts: a 'common' part that is executed for
      ! all models inheriting from atype_SMB_model, and a 'specific' part that is
      ! only executed for each specific model class. The specific parts are defined
      ! in the deferred procedures 'allocate_SMB_model', 'initialise_SMB_model', etc.

      procedure, public :: allocate_SMB_model
      procedure, public :: deallocate_SMB_model
      procedure, public :: initialise_SMB_model
      procedure, public :: run_SMB_model
      procedure, public :: remap_model      => remap_model_abs

      procedure(SMB_model_allocate_ifc),   deferred :: allocate
      procedure(SMB_model_deallocate_ifc), deferred :: deallocate
      procedure(SMB_model_initialise_ifc), deferred :: initialise
      procedure(SMB_model_run_ifc),        deferred :: run
      procedure(remap_SMB_model_ifc),      deferred :: remap_SMB_model

      ! Factory functions to create model context objects
      procedure, nopass, public :: ct_remap

  end type atype_SMB_model

  ! Context classes for allocate/initialise/run/remap
  ! =================================================

  type, extends(atype_model_context_remap) :: type_SMB_model_context_remap
    real(dp)                               :: time
    character(len=3)                       :: region_name
    type(type_reference_geometry), pointer :: refgeo_init
    type(type_reference_geometry), pointer :: refgeo_PD
    type(type_ice_model),          pointer :: ice
  end type type_SMB_model_context_remap

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine SMB_model_allocate_ifc( self, region_name, mesh)
      import atype_SMB_model, type_mesh
      class(atype_SMB_model),  intent(inout) :: self
      character(len=*),        intent(in   ) :: region_name
      type(type_mesh), target, intent(in   ) :: mesh
    end subroutine SMB_model_allocate_ifc

    subroutine SMB_model_deallocate_ifc( self)
      import atype_SMB_model
      class(atype_SMB_model), intent(inout) :: self
    end subroutine SMB_model_deallocate_ifc

    subroutine SMB_model_initialise_ifc( self, ice, refgeo_init, refgeo_PD)
      import atype_SMB_model, type_ice_model, type_reference_geometry
      class(atype_SMB_model),        intent(inout) :: self
      type(type_ice_model),          intent(in   ) :: ice
      type(type_reference_geometry), intent(in   ) :: refgeo_init
      type(type_reference_geometry), intent(in   ) :: refgeo_PD
    end subroutine SMB_model_initialise_ifc

    subroutine SMB_model_run_ifc( self, time, ice, climate, grid_smooth)
      import atype_SMB_model, dp, type_ice_model, type_climate_model, type_grid
      class(atype_SMB_model),   intent(inout) :: self
      real(dp),                 intent(in   ) :: time
      type(type_ice_model),     intent(in   ) :: ice
      type(type_climate_model), intent(inout) :: climate
      type(type_grid),          intent(in   ) :: grid_smooth
    end subroutine SMB_model_run_ifc

    subroutine remap_SMB_model_ifc( self, context)
      import atype_SMB_model, type_SMB_model_context_remap
      class(atype_SMB_model),                     intent(inout) :: self
      type(type_SMB_model_context_remap), target, intent(in   ) :: context
    end subroutine remap_SMB_model_ifc

  end interface

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  interface

    module subroutine remap_model_abs( self, context)
      class(atype_SMB_model),                   intent(inout) :: self
      class(atype_model_context_remap), target, intent(in   ) :: context
    end subroutine remap_model_abs

    module function ct_remap( mesh_new, time, region_name, refgeo_init, refgeo_PD, ice) result( context)
      type(type_mesh),               target, intent(in) :: mesh_new
      real(dp),                              intent(in) :: time
      character(len=3),                      intent(in) :: region_name
      type(type_reference_geometry), target, intent(in) :: refgeo_init, refgeo_PD
      type(type_ice_model),          target, intent(in) :: ice
      type(type_SMB_model_context_remap)                :: context
    end function ct_remap

  end interface

contains

  subroutine allocate_SMB_model( self, name, region_name, mesh)
    !< Allocate stuff that is common to all SMB models
    !< (call this from your demo model-specific allocate routine)

    ! In/output variables:
    class(atype_SMB_model),  intent(inout) :: self
    character(len=*),        intent(in   ) :: name
    character(len=*),        intent(in   ) :: region_name
    type(type_mesh), target, intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_SMB_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate stuff that is common to all models
    call self%allocate_model( name, region_name, mesh)

    ! Allocate stuff that is specific to demo models

    ! Allocate generic fields
    call self%create_field( self%SMB, self%wSMB, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'SMB', &
      long_name = 'surface mass balance', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model

  subroutine deallocate_SMB_model( self)
    !< Deallocate stuff that is common to all SMB models
    !< (call this from your SMB model-specific deallocate routine)

    ! In/output variables:
    class(atype_SMB_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_SMB_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate stuff that is common to all models
    call self%deallocate_model()

    ! Deallocate stuff that is specific to SMB models

    nullify( self%SMB)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_SMB_model

  subroutine initialise_SMB_model( self)
    !< Initialise stuff that is common to all SMB models
    !< (call this from your SMB model-specific initialise routine)

    ! In/output variables:
    class(atype_SMB_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_SMB_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise stuff that is common to all models
    call self%initialise_model()

    ! Initialise stuff that is specific to SMB models

    ! Set time of next calculation to start time
    self%t_next = C%start_time_of_run

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model

  subroutine run_SMB_model( self, time, do_run_SMB_model)
    !< Run stuff that is common to all SMB models
    !< (call this from your SMB model-specific run routine)

    ! In/output variables:
    class(atype_SMB_model), intent(inout) :: self
    real(dp),                   intent(in   ) :: time
    logical,                    intent(  out) :: do_run_SMB_model

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_SMB_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run stuff that is common to all models
    call self%run_model()

    ! Run stuff that is specific to SMB models

    ! Check if we need to calculate a new SMB
    do_run_SMB_model = .false.
    if (C%do_asynchronous_SMB) then
      ! Asynchronous coupling: do not calculate a new SMB in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next SMB time step
      if (time == self%t_next) then
        ! Go on to calculate a new SMB
        do_run_SMB_model = .true.
        self%t_next = time + C%dt_SMB
      elseif (time > self%t_next) then
        ! This should not be possible
        call crash('overshot the SMB time step')
      else
        ! It is not yet time to calculate a new SMB
        do_run_SMB_model = .false.
      end if

    else
      ! Synchronous coupling: calculate a new SMB in every model loop
      do_run_SMB_model = .true.
      self%t_next = time + C%dt_SMB
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model

end module SMB_model_basic
