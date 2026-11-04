module remapping_calc_CSR_operators

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use remapping_methods, only: type_remapping_method, remapping_method
  use ugrid_basic, only: type_ugrid
  use CSR_matrix_mod, only: type_CSR_matrix_dp

  implicit none

  private

  public :: calc_CSR_remapping_operator

contains

  subroutine calc_CSR_remapping_operator( ugrid_src, ugrid_dst, choice_remapping_method, M)

    ! In/output variables:
    type(type_ugrid),            intent(in   ) :: ugrid_src
    type(type_ugrid),            intent(in   ) :: ugrid_dst
    type(type_remapping_method), intent(in   ) :: choice_remapping_method
    type(type_CSR_matrix_dp),    intent(  out) :: M

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_CSR_remapping_operator'

    ! Add routine to path
    call init_routine( routine_name)

    if     (choice_remapping_method == remapping_method%nearest_neighbour()) then
      call calc_CSR_remapping_operator_nearest_neighbour( ugrid_src, ugrid_dst, M)
    elseif (choice_remapping_method == remapping_method%linear()) then
      call calc_CSR_remapping_operator_linear( ugrid_src, ugrid_dst, M)
    elseif (choice_remapping_method == remapping_method%first_order_conservative()) then
      call calc_CSR_remapping_operator_first_order_conservative( ugrid_src, ugrid_dst, M)
    elseif (choice_remapping_method == remapping_method%second_order_conservative()) then
      call calc_CSR_remapping_operator_second_order_conservative( ugrid_src, ugrid_dst, M)
    else
      call crash("unsupported remapping method '" // choice_remapping_method%str() // "'")
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_CSR_remapping_operator

  subroutine calc_CSR_remapping_operator_nearest_neighbour( ugrid_src, ugrid_dst, M)

    ! In/output variables:
    type(type_ugrid),         intent(in   ) :: ugrid_src
    type(type_ugrid),         intent(in   ) :: ugrid_dst
    type(type_CSR_matrix_dp), intent(  out) :: M

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_CSR_remapping_operator_nearest_neighbour'

    ! Add routine to path
    call init_routine( routine_name)

    call crash('fixme')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_CSR_remapping_operator_nearest_neighbour

  subroutine calc_CSR_remapping_operator_linear( ugrid_src, ugrid_dst, M)

    ! In/output variables:
    type(type_ugrid),         intent(in   ) :: ugrid_src
    type(type_ugrid),         intent(in   ) :: ugrid_dst
    type(type_CSR_matrix_dp), intent(  out) :: M

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_CSR_remapping_operator_linear'

    ! Add routine to path
    call init_routine( routine_name)

    call crash('fixme')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_CSR_remapping_operator_linear

  subroutine calc_CSR_remapping_operator_first_order_conservative( ugrid_src, ugrid_dst, M)

    ! In/output variables:
    type(type_ugrid),         intent(in   ) :: ugrid_src
    type(type_ugrid),         intent(in   ) :: ugrid_dst
    type(type_CSR_matrix_dp), intent(  out) :: M

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_CSR_remapping_operator_first_order_conservative'

    ! Add routine to path
    call init_routine( routine_name)

    call crash('fixme')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_CSR_remapping_operator_first_order_conservative

  subroutine calc_CSR_remapping_operator_second_order_conservative( ugrid_src, ugrid_dst, M)

    ! In/output variables:
    type(type_ugrid),         intent(in   ) :: ugrid_src
    type(type_ugrid),         intent(in   ) :: ugrid_dst
    type(type_CSR_matrix_dp), intent(  out) :: M

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_CSR_remapping_operator_second_order_conservative'

    ! Add routine to path
    call init_routine( routine_name)

    call crash('fixme')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_CSR_remapping_operator_second_order_conservative

end module remapping_calc_CSR_operators