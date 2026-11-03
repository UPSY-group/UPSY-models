submodule( ugrid_basic) build_from_xy_grid_mod

contains

  subroutine build_from_xy_grid( ugrid, xy_grid, xy_grid_Arakawa_grid)

    ! In/output variables:
    type(type_ugrid),        intent(inout) :: ugrid
    type(type_xy_grid),      intent(in   ) :: xy_grid
    type(type_Arakawa_grid), intent(in   ) :: xy_grid_Arakawa_grid

    ! Local variables:
    character(len=*), parameter :: routine_name = 'build_from_xy_grid'

    ! Add routine to path
    call init_routine( routine_name)

    call crash('fixme!')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine build_from_xy_grid

end submodule build_from_xy_grid_mod
