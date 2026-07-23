module SMB_model_basic

  use precisions, only: dp
  use parameters, only: pi, NaN
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use SMB_model_data, only: atype_SMB_model_data
  use mpi_f08, only: MPI_WIN
  use ice_model_types, only: type_ice_model
  use climate_model_types, only: type_climate_model
  use grid_types, only: type_grid
  use reference_geometry_types, only: type_reference_geometry

  implicit none

  private

  public :: atype_SMB_model

  type, abstract, extends(atype_SMB_model_data) :: atype_SMB_model
    !< Stuff that is common to all SMB models
    !<
    !< (except for the variables that we want other models to
    !< be able to access, which are already defined in atype_SMB_model_data)

    real(dp) :: t_next   !< Time when the SMB model should be run next

    contains

      ! Procedures for model memory management and operation
      procedure, public :: allocate   => SMB_model_allocate
      procedure, public :: deallocate => SMB_model_deallocate
      procedure, public :: initialise => SMB_model_initialise
      procedure, public :: run        => SMB_model_run
      procedure, public :: remap      => SMB_model_remap

      ! Deferred procedures that must be overridden by each individual SMB model implementation
      procedure(SMB_model_allocate_ifc),   deferred :: allocate_SMB_model
      procedure(SMB_model_deallocate_ifc), deferred :: deallocate_SMB_model
      procedure(SMB_model_initialise_ifc), deferred :: initialise_SMB_model
      procedure(SMB_model_run_ifc),        deferred :: run_SMB_model
      procedure(SMB_model_remap_ifc),      deferred :: remap_SMB_model

  end type atype_SMB_model

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine SMB_model_allocate_ifc( self)
      import atype_SMB_model
      class(atype_SMB_model),  intent(inout) :: self
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

    subroutine SMB_model_remap_ifc( self, mesh_new, time, refgeo_init, refgeo_PD, ice)
      import atype_SMB_model, type_mesh, dp, type_reference_geometry, type_ice_model
      class(atype_SMB_model),                intent(inout) :: self
      type(type_mesh), target,               intent(in   ) :: mesh_new
      real(dp),                              intent(in   ) :: time
      type(type_reference_geometry), target, intent(in   ) :: refgeo_init, refgeo_PD
      type(type_ice_model),          target, intent(in   ) :: ice
    end subroutine SMB_model_remap_ifc

  end interface

contains

  subroutine SMB_model_allocate( self, name, region_name, mesh)

    ! In/output variables:
    class(atype_SMB_model),  intent(inout) :: self
    character(len=*),        intent(in   ) :: name
    character(len=*),        intent(in   ) :: region_name
    type(type_mesh), target, intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate stuff that is common to all models
    call self%allocate_model( name, region_name, mesh)

    ! Allocate stuff that is common to all SMB models

    ! Allocate generic fields
    call self%create_field( self%SMB, self%wSMB, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'SMB', &
      long_name = 'surface mass balance', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    ! Allocate stuff that is specific to each individual SMB model implementation
    call self%allocate_SMB_model()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine SMB_model_allocate

  subroutine SMB_model_deallocate( self)

    ! In/output variables:
    class(atype_SMB_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate stuff that is common to all models
    call self%deallocate_model()

    ! Deallocate stuff that is common to all SMB models

    nullify( self%SMB)

    ! Deallocate stuff that is specific to each individual SMB model implementation
    call self%deallocate_SMB_model()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine SMB_model_deallocate

  subroutine SMB_model_initialise( self, ice, refgeo_init, refgeo_PD)

    ! In/output variables:
    class(atype_SMB_model),        intent(inout) :: self
    type(type_ice_model),          intent(in   ) :: ice
    type(type_reference_geometry), intent(in   ) :: refgeo_init
    type(type_reference_geometry), intent(in   ) :: refgeo_PD

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise stuff that is common to all models
    call self%initialise_model()

    ! Initialise stuff that is common to all SMB models

    ! Set time of next calculation to start time
    self%t_next = C%start_time_of_run

    ! Initialise stuff that is specific to each individual SMB model implementation
    call self%initialise_SMB_model( ice, refgeo_init, refgeo_PD)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine SMB_model_initialise

  subroutine SMB_model_run( self, time, ice, climate, grid_smooth)

    ! In/output variables:
    class(atype_SMB_model),   intent(inout) :: self
    real(dp),                 intent(in   ) :: time
    type(type_ice_model),     intent(in   ) :: ice
    type(type_climate_model), intent(inout) :: climate
    type(type_grid),          intent(in   ) :: grid_smooth

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_run'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run stuff that is common to all models
    call self%run_model()

    ! Run stuff that is common to all SMB models

    ! Check if we need to calculate a new SMB
    if (C%do_asynchronous_SMB) then
      ! Asynchronous coupling: do not calculate a new SMB in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next SMB time step
      if (time == self%t_next) then
        ! Go on to calculate a new SMB
        self%t_next = time + C%dt_SMB

        ! Run stuff that is specific to each individual SMB model implementation
        call self%run_SMB_model( time, ice, climate, grid_smooth)

      elseif (time > self%t_next) then
        ! This should not be possible
        call crash('overshot the SMB time step')
      else
        ! It is not yet time to calculate a new SMB
      end if

    else
      ! Synchronous coupling: calculate a new SMB in every model loop
      self%t_next = time + C%dt_SMB

      ! Run stuff that is specific to each individual SMB model implementation
      call self%run_SMB_model( time, ice, climate, grid_smooth)

    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine SMB_model_run

  subroutine SMB_model_remap( self, mesh_new, time, refgeo_init, refgeo_PD, ice)

    ! In/output variables:
    class(atype_SMB_model),                intent(inout) :: self
    type(type_mesh), target,               intent(in   ) :: mesh_new
    real(dp),                              intent(in   ) :: time
    type(type_reference_geometry), target, intent(in   ) :: refgeo_init, refgeo_PD
    type(type_ice_model),          target, intent(in   ) :: ice

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap stuff that is common to all models
    call self%remap_model( mesh_new)

    ! Remap stuff that is common to all SMB models

    call self%remap_field( mesh_new, 'SMB', self%SMB)

    ! Remap stuff that is specific to each individual SMB model implementation
    call self%remap_SMB_model( mesh_new, time, refgeo_init, refgeo_PD, ice)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine SMB_model_remap

end module SMB_model_basic
