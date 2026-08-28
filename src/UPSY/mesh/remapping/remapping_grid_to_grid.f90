module remapping_grid_to_grid

  use mpi_basic, only: par
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use grid_types, only: type_grid
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use remapping_types, only: type_map

  implicit none

  private

  public :: create_map_from_xy_grid_to_xy_grid

contains

  subroutine create_map_from_xy_grid_to_xy_grid( grid_src, grid_dst, map)
    !< Create a new mapping object from one x/y-grid to another x/y-grid

    ! By default uses 2nd-order conservative remapping.

    ! In/output variables
    type(type_grid), intent(in   ) :: grid_src
    type(type_grid), intent(in   ) :: grid_dst
    type(type_map),  intent(inout) :: map

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'create_map_from_xy_grid_to_xy_grid'

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise map metadata
    if (map%is_in_use) call crash('this map is already in use!')
    map%is_in_use = .true.
    map%name_src  = grid_src%name
    map%name_dst  = grid_dst%name
    map%method    = '2nd_order_conservative'

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map_from_xy_grid_to_xy_grid

end module remapping_grid_to_grid
