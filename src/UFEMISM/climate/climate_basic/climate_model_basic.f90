module climate_model_basic

  use precisions, only: dp
  use parameters, only: pi, NaN
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use climate_model_data, only: atype_climate_model_data
  use mpi_f08, only: MPI_WIN
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry

  implicit none

  private

  public :: atype_climate_model

  type, abstract, extends(atype_climate_model_data) :: atype_climate_model
    !< Stuff that is common to all climate models
    !<
    !< (except for the variables that we want other models to
    !< be able to access, which are already defined in atype_climate_model_data)

      real(dp) :: t_next   !< Time when the climate model should be run next

    contains

      ! Type-bound procedures that apply to all climate models
      procedure, public :: allocate_climate_model
      procedure, public :: deallocate_climate_model
      procedure, public :: initialise_climate_model
      procedure, public :: run_climate_model
      procedure, public :: remap_climate_model

      ! Deferred procedures that must be defined by each individual climate model
      procedure(climate_model_allocate_ifc),   deferred :: allocate
      procedure(climate_model_deallocate_ifc), deferred :: deallocate
      procedure(climate_model_initialise_ifc), deferred :: initialise
      procedure(climate_model_run_ifc),        deferred :: run
      procedure(climate_model_remap_ifc),      deferred :: remap

  end type atype_climate_model

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine climate_model_allocate_ifc( self, region_name, mesh)
      import atype_climate_model, type_mesh
      class(atype_climate_model), intent(inout) :: self
      character(len=*),           intent(in   ) :: region_name
      type(type_mesh), target,    intent(in   ) :: mesh
    end subroutine climate_model_allocate_ifc

    subroutine climate_model_deallocate_ifc( self)
      import atype_climate_model
      class(atype_climate_model), intent(inout) :: self
    end subroutine climate_model_deallocate_ifc

    subroutine climate_model_initialise_ifc( self, refgeo_PD, refgeo_init)
      import atype_climate_model, type_reference_geometry
      class(atype_climate_model),    intent(inout) :: self
      type(type_reference_geometry), intent(in   ) :: refgeo_PD
      type(type_reference_geometry), intent(in   ) :: refgeo_init
    end subroutine climate_model_initialise_ifc

    subroutine climate_model_run_ifc( self, ice, time)
      import atype_climate_model, type_ice_model, dp
      class(atype_climate_model), intent(inout) :: self
      type(type_ice_model),       intent(in   ) :: ice
      real(dp),                   intent(in   ) :: time
    end subroutine climate_model_run_ifc

    subroutine climate_model_remap_ifc( self, mesh_new)
      import atype_climate_model, type_mesh
      class(atype_climate_model), intent(inout) :: self
      type(type_mesh), target,    intent(in   ) :: mesh_new
    end subroutine climate_model_remap_ifc

  end interface

contains

  subroutine allocate_climate_model( self, name, region_name, mesh)
    !< Allocate stuff that is common to all climate models
    !< (call this from your climate model-specific allocate routine)

    ! In/output variables:
    class(atype_climate_model), intent(inout) :: self
    character(len=*),           intent(in   ) :: name
    character(len=*),           intent(in   ) :: region_name
    type(type_mesh), target,    intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_climate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all models
    call self%allocate_model( name, region_name, mesh)

    ! Allocate all the stuff that is specific to climate models

    call self%create_field( self%T2m, self%wT2m, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m', &
      long_name = 'Monthly mean 2-m air temperature', &
      units     = 'K', &
      remap_method = 'reallocate')

    call self%create_field( self%Precip, self%wPrecip, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Precip', &
      long_name = 'Monthly total precipitation', &
      units     = 'm.w.e.', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_climate_model

  subroutine deallocate_climate_model( self)
    !< Deallocate stuff that is common to all climate models
    !< (call this from your climate model-specific deallocate routine)

    ! In/output variables:
    class(atype_climate_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_climate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate stuff that is common to all models
    call self%deallocate_model()

    ! Deallocate stuff that is specific to climate models

    nullify( self%T2m)
    nullify( self%Precip)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_climate_model

  subroutine initialise_climate_model( self)
    !< Initialise stuff that is common to all climate models
    !< (call this from your climate model-specific initialise routine)

    ! In/output variables:
    class(atype_climate_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_climate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise stuff that is common to all models
    call self%initialise_model()

    ! Initialise stuff that is specific to climate models

    ! Set time of next calculation to start time
    self%t_next = C%start_time_of_run

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_climate_model

  subroutine run_climate_model( self, time, do_run_climate_model)
    !< Run stuff that is common to all climate models
    !< (call this from your climate model-specific run routine)

    ! In/output variables:
    class(atype_climate_model), intent(inout) :: self
    real(dp),                   intent(in   ) :: time
    logical,                    intent(  out) :: do_run_climate_model

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_climate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run stuff that is common to all models
    call self%run_model()

    ! Run stuff that is specific to climate models

    ! Check if we need to calculate a new climate
    do_run_climate_model = .false.
    if (C%do_asynchronous_climate) then
      ! Asynchronous coupling: do not calculate a new climate in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next climate time step
      if (time == self%t_next) then
        ! Go on to calculate a new climate
        do_run_climate_model = .true.
        self%t_next = time + C%dt_climate
      elseif (time > self%t_next) then
        ! This should not be possible
        call crash('overshot the climate time step')
      else
        ! It is not yet time to calculate a new climate
        do_run_climate_model = .false.
      end if

    else
      ! Synchronous coupling: calculate a new climate in every model loop
      do_run_climate_model = .true.
      self%t_next = time + C%dt_climate
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_climate_model

  subroutine remap_climate_model( self, mesh_new)
    !< Remap stuff that is common to all climate models
    !< (call this from your climate model-specific remap routine)

    ! In/output variables:
    class(atype_climate_model), intent(inout) :: self
    type(type_mesh), target,    intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_climate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap stuff that is common to all models
    call self%remap_model( mesh_new)

    ! Remap stuff that is specific to climate models

    call self%remap_field( mesh_new, 'T2m'   , self%T2m)
    call self%remap_field( mesh_new, 'Precip', self%Precip)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_climate_model

end module climate_model_basic
