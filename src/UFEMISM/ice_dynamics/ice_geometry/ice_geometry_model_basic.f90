module ice_geometry_model_basic

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use Arakawa_grid_mod, only: Arakawa_grid
  use mesh_types, only: type_mesh

  implicit none

  private

  public :: type_ice_geometry_model

  type, extends(atype_ice_geometry_model_data) :: type_ice_geometry_model

    contains

      procedure, public :: allocate   => allocate_ice_geometry_model
      procedure, public :: deallocate => deallocate_ice_geometry_model
      procedure, public :: remap      => remap_ice_geometry_model

      final :: finalise_ice_geometry_model

  end type type_ice_geometry_model

contains

  subroutine allocate_ice_geometry_model( self, name, region_name, mesh)

    ! In/output variables:
    class(type_ice_geometry_model), intent(inout) :: self
    character(len=*),               intent(in   ) :: name
    character(len=*),               intent(in   ) :: region_name
    type(type_mesh), target,        intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_ice_geometry_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all models
    call self%allocate_model( name, region_name, mesh)

    ! Allocate all the stuff that is specific to the ice_geometry model

    call self%create_field( self%Hi, self%wHi, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Hi', &
      long_name = 'Ice thickness', &
      units     = 'm', &
      remap_method = 'reallocate')

    call self%create_field( self%Hb, self%wHb, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Hb', &
      long_name = 'Bedrock elevation (w.r.t. PD sea level)', &
      units     = 'm', &
      remap_method = 'reallocate')

    call self%create_field( self%SL, self%wSL, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'SL', &
      long_name = 'Geoid elevation (w.r.t. PD sea level)', &
      units     = 'm', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_ice_geometry_model

  subroutine deallocate_ice_geometry_model( self)

    ! In/output variables:
    class(type_ice_geometry_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_ice_geometry_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate stuff that is common to all models
    call self%deallocate_model()

    ! Deallocate stuff that is specific to the ice_geometry model
    nullify( self%Hi)
    nullify( self%Hb)
    nullify( self%SL)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_ice_geometry_model

  subroutine remap_ice_geometry_model( self, mesh_new)

    ! In/output variables:
    class(type_ice_geometry_model), intent(inout) :: self
    type(type_mesh), target,        intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_ice_geometry_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap stuff that is common to all models
    call self%remap_model( mesh_new)

    ! Remap stuff that is specific to ice_geometry models
    call self%remap_field( mesh_new, 'Hi', self%Hi)
    call self%remap_field( mesh_new, 'Hb', self%Hb)
    call self%remap_field( mesh_new, 'SL', self%SL)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_ice_geometry_model

  subroutine finalise_ice_geometry_model( self)

    ! In/output variables:
    type(type_ice_geometry_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'finalise_ice_geometry_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%deallocate()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine finalise_ice_geometry_model

end module ice_geometry_model_basic