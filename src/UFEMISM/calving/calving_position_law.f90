module calving_position

  use parameters, only: pi
  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use ice_geometry_basics, only: is_floating 
  use ice_geometry_calculations, only: type_ice_geometry 
  use calving_model_basic, only: atype_calving_model, type_calving_model_context_allocate, &
    type_calving_model_context_initialise, type_calving_model_context_run, &
    type_calving_model_context_remap

  implicit none

  private

  public :: type_calving_model_position

  type, extends(atype_calving_model) :: type_calving_model_position

    contains

      procedure, public :: allocate_calving_model   => allocate_calving_model_position_abs
      procedure, public :: deallocate_calving_model => deallocate_calving_model_position_abs
      procedure, public :: initialise_calving_model => initialise_calving_model_position_abs
      procedure, public :: run_calving_model        => run_calving_model_position_abs
      procedure, public :: remap_calving_model      => remap_calving_model_position_abs

      procedure, private :: run_calving_model_position
      procedure, private :: run_threshold_thickness_front_calving
      procedure, private :: run_threshold_thickness_shelf_calving
      
  end type type_calving_model_position

contains

  subroutine allocate_calving_model_position_abs( self, context)

    ! In/output variables:
    class(type_calving_model_position),                intent(inout) :: self
    type(type_calving_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_calving_model_position_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_calving_model_position_abs

  subroutine deallocate_calving_model_position_abs( self)

    ! In/output variables:
    class(type_calving_model_position), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'deallocate_calving_model_position_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_calving_model_position_abs

  subroutine initialise_calving_model_position_abs( self, context)

    ! In/output variables:
    class(type_calving_model_position),                  intent(inout) :: self
    type(type_calving_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_calving_model_position_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_calving_model_position_abs

  subroutine run_calving_model_position_abs( self, context)

    ! In/output variables:
    class(type_calving_model_position),           intent(inout) :: self
    type(type_calving_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_calving_model_position_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call self%run_calving_model_position( self%mesh, context%time, context%icegeom)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_calving_model_position_abs

  subroutine remap_calving_model_position_abs( self, context)

    ! In/output variables:
    class(type_calving_model_position),             intent(inout) :: self
    type(type_calving_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_calving_model_position_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_calving_model_position_abs

  subroutine run_calving_model_position( self, mesh, time, icegeom)

    ! In/output variables:
    class(type_calving_model_position), intent(inout) :: self
    type(type_mesh),                    intent(in   ) :: mesh
    real(dp),                           intent(in   ) :: time
    type(type_ice_geometry),            intent(in   ) :: icegeom

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_calving_model_position'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run the calving model that follow the chosen position law
    select case (C%choice_calving_law )
    case default
      call crash('unknown choice_calving_model_position "' // trim( C%choice_calving_law ) // '"')
    case ('threshold_thickness_front')
      call self%run_threshold_thickness_front_calving(mesh, icegeom)
    case ('threshold_thickness_all')
      call self%run_threshold_thickness_shelf_calving(mesh, icegeom)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_calving_model_position

  subroutine run_threshold_thickness_front_calving(self, mesh, icegeom)
    ! Calculate where ice thickness should be set to 0 because of calving, 
    ! when the ice thickness at the floating front is under a certain (user defined) thickness.

    ! In/output variables
    class(type_calving_model_position), intent(inout) :: self
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_geometry),            intent(in   ) :: icegeom

    ! Local variables:
    character(len=256), parameter :: routine_name = ' run_threshold_thickness_front_calving'
    integer                       :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      if ( icegeom%mask_cf_fl(vi) .and. icegeom%Hi_eff( vi) < C%calving_threshold_thickness_shelf) then
        self%Hi_calved(vi) = 0.0_dp
      end if
    end do
    
    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_threshold_thickness_front_calving

  subroutine run_threshold_thickness_shelf_calving(self, mesh, icegeom)
    ! Calculate where ice thickness should be set to 0 because of calving, 
    ! when the (entire) floating ice is under a certain user defined thickness.

    ! In/output variables
    class(type_calving_model_position), intent(inout) :: self
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_ice_geometry),            intent(in   ) :: icegeom

    ! Local variables:
    character(len=256), parameter :: routine_name = ' run_threshold_thickness_shelf_calving'
    integer                       :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      if ( is_floating( icegeom%Hi_eff( vi), icegeom%Hb( vi), icegeom%SL( vi)) .and. icegeom%Hi_eff( vi) < C%calving_threshold_thickness_shelf)  then
        self%Hi_calved(vi) = 0.0_dp
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_threshold_thickness_shelf_calving

end module calving_position