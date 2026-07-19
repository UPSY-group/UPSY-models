module ocean_model_basic

  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use ocean_model_data, only: atype_ocean_model_data
  use mpi_f08, only: MPI_WIN
  use ice_model_types, only: type_ice_model
  use climate_model_types, only: type_climate_model
  use grid_types, only: type_grid

  implicit none

  private

  public :: atype_ocean_model

  type, abstract, extends(atype_ocean_model_data) :: atype_ocean_model
    !< Stuff that is common to all ocean models
    !<
    !< (except for the variables that we want other models to
    !< be able to access, which are already defined in atype_ocean_model_data)

    real(dp) :: t_next   !< Time when the ocean model should be run next

    contains

      ! Type-bound procedures that apply to all ocean models
      procedure, public :: allocate_ocean_model
      procedure, public :: deallocate_ocean_model
      procedure, public :: initialise_ocean_model
      procedure, public :: run_ocean_model
      procedure, public :: remap_ocean_model

      ! Deferred procedures that must be defined by each individual ocean model
      procedure(ocean_model_allocate_ifc),   deferred :: allocate
      procedure(ocean_model_deallocate_ifc), deferred :: deallocate
      procedure(ocean_model_initialise_ifc), deferred :: initialise
      procedure(ocean_model_run_ifc),        deferred :: run
      procedure(ocean_model_remap_ifc),      deferred :: remap

  end type atype_ocean_model

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine ocean_model_allocate_ifc( self, name, region_name, mesh)
      import atype_ocean_model, type_mesh
      class(atype_ocean_model), intent(inout) :: self
      character(len=*),         intent(in   ) :: name
      character(len=*),         intent(in   ) :: region_name
      type(type_mesh), target,  intent(in   ) :: mesh
    end subroutine ocean_model_allocate_ifc

    subroutine ocean_model_deallocate_ifc( self)
      import atype_ocean_model
      class(atype_ocean_model), intent(inout) :: self
    end subroutine ocean_model_deallocate_ifc

    subroutine ocean_model_initialise_ifc( self)
      import atype_ocean_model
      class(atype_ocean_model), intent(inout) :: self
    end subroutine ocean_model_initialise_ifc

    subroutine ocean_model_run_ifc( self)
      import atype_ocean_model
      class(atype_ocean_model), intent(inout) :: self
    end subroutine ocean_model_run_ifc

    subroutine ocean_model_remap_ifc( self, mesh_new)
      import atype_ocean_model, type_mesh
      class(atype_ocean_model), intent(inout) :: self
      type(type_mesh), target,  intent(in   ) :: mesh_new
    end subroutine ocean_model_remap_ifc

  end interface

contains

  subroutine allocate_ocean_model( self, name, region_name, mesh)
    !< Allocate stuff that is common to all ocean models
    !< (call this from your ocean model-specific allocate routine)

    ! In/output variables:
    class(atype_ocean_model), intent(inout) :: self
    character(len=*),         intent(in   ) :: name
    character(len=*),         intent(in   ) :: region_name
    type(type_mesh), target,  intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_ocean_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate stuff that is common to all models
    call self%allocate_model( name, region_name, mesh)

    ! Allocate stuff that is specific to ocean models

    ! Main data fields
    ! ================

    call self%create_field( self%T, self%wT, &
      self%mesh, Arakawa_grid%a(), third_dimension%ocean_depth( C%nz_ocean), &
      name      = 'T_ocean', &
      long_name = 'Ocean temperature', &
      units     = 'degrees Celsius')

    call self%create_field( self%S, self%wS, &
      self%mesh, Arakawa_grid%a(), third_dimension%ocean_depth( C%nz_ocean), &
      name      = 'S_ocean', &
      long_name = 'Ocean salinity', &
      units     = 'PSU')

    ! Secondary data fields
    ! =====================

    call self%create_field( self%T_draft, self%wT_draft, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'T_draft', &
      long_name = 'Ocean temperature at ice draft', &
      units     = 'degrees Celsius')

    call self%create_field( self%T_freezing_point, self%wT_freezing_point, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'T_freezing_point', &
      long_name = 'Ocean freezing temperature at ice draft', &
      units     = 'degrees Celsius')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_ocean_model

  subroutine deallocate_ocean_model( self)
    !< Deallocate stuff that is common to all ocean models
    !< (call this from your ocean model-specific deallocate routine)

    ! In/output variables:
    class(atype_ocean_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_ocean_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate stuff that is common to all models
    call self%deallocate_model()

    ! Deallocate stuff that is specific to ocean models

    nullify( self%T)
    nullify( self%S)
    nullify( self%T_draft)
    nullify( self%T_freezing_point)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_ocean_model

  subroutine initialise_ocean_model( self)
    !< Initialise stuff that is common to all ocean models
    !< (call this from your ocean model-specific initialise routine)

    ! In/output variables:
    class(atype_ocean_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_ocean_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise stuff that is common to all models
    call self%initialise_model()

    ! Initialise stuff that is specific to ocean models

    ! Set time of next calculation to start time
    self%t_next = C%start_time_of_run

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model

  subroutine run_ocean_model( self)
    !< Run stuff that is common to all ocean models
    !< (call this from your ocean model-specific run routine)

    ! In/output variables:
    class(atype_ocean_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_ocean_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run stuff that is common to all models
    call self%run_model()

    ! Run stuff that is specific to ocean models

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_ocean_model

  subroutine remap_ocean_model( self, mesh_new)
    !< Remap stuff that is common to all ocean models
    !< (call this from your ocean model-specific remap routine)

    ! In/output variables:
    class(atype_ocean_model), intent(inout) :: self
    type(type_mesh), target,  intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_ocean_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap stuff that is common to all models
    call self%remap_model( mesh_new)

    ! Remap stuff that is specific to ocean models

    call self%remap_field( mesh_new, 'T'               , self%T)
    call self%remap_field( mesh_new, 'S'               , self%S)
    call self%remap_field( mesh_new, 'T_draft'         , self%T_draft)
    call self%remap_field( mesh_new, 'T_freezing_point', self%T_freezing_point)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_ocean_model

end module ocean_model_basic
