module ice_velocity_model_clean

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use ice_velocity_model_basic, only: atype_ice_velocity_model

  implicit none

  private

  public :: type_ice_velocity_model_clean

  ! Temporary concrete extension, to be used until all velocity models have been converted
  ! to extensions of the abstract base type
  type, extends(atype_ice_velocity_model) :: type_ice_velocity_model_clean
  end type type_ice_velocity_model_clean

  ! Interfaces for procedures defined in submodules
  interface

    ! module subroutine remap_ice_velocity_model( self, mesh_old, mesh_new)
    !   class(type_ice_velocity_model),       intent(inout) :: self
    !   type(type_mesh),                      intent(in   ) :: mesh_old
    !   type(type_mesh),                      intent(in   ) :: mesh_new
    ! end subroutine remap_ice_velocity_model

  end interface

contains

end module ice_velocity_model_clean