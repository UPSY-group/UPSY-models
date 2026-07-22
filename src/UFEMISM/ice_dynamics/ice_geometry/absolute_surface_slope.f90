submodule( ice_geometry_model_basic) absolute_surface_slope

contains

  subroutine calc_absolute_surface_slope( self)
    !< Determine the ice-filled fraction and effective ice thickness of floating margin pixels

    ! In- and output variables
    class(type_ice_geometry_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'calc_absolute_surface_slope'
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: dHs_dx, dHs_dy
    integer                                          :: vi

    ! Add routine to path
    call init_routine( routine_name)

    call ddx_a_a_2D( self%mesh, self%Hs, dHs_dx)
    call ddy_a_a_2D( self%mesh, self%Hs, dHs_dy)

    do vi = self%mesh%vi1, self%mesh%vi2
      self%Hs_slope( vi) = sqrt( dHs_dx( vi)**2 + dHs_dy( vi)**2)
    end do

    call checksum( self%mesh%pai_V, self%Hs_slope, 'Hs_slope')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_absolute_surface_slope

end submodule absolute_surface_slope