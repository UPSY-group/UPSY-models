module ice_shelf_base_slopes

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use mesh_disc_apply_operators, only: ddx_a_b_2D, ddy_a_b_2D

  implicit none

  private

  public :: calc_ice_shelf_base_slopes

contains

  subroutine calc_ice_shelf_base_slopes( mesh, ice)

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_ice_shelf_base_slopes'

    ! Add routine to path
    call init_routine( routine_name)

    ! Straightforward calculation everywhere
    call ddx_a_b_2D( mesh, ice%Hib, ice%dHib_dx_b)
    call ddy_a_b_2D( mesh, ice%Hib, ice%dHib_dy_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ice_shelf_base_slopes

end module ice_shelf_base_slopes
