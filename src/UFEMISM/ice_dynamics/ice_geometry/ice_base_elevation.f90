submodule( ice_geometry_model_basic) ice_base_elevation

contains

  subroutine calc_ice_base_elevation( self)

    ! In- and output variables
    class(type_ice_geometry_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_ice_base_elevation'
    integer                     :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = self%mesh%vi1, self%mesh%vi2
      self%Hib( vi) = ice_surface_elevation( self%Hi( vi), self%Hb( vi), self%SL( vi)) - self%Hi( vi)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ice_base_elevation

end submodule ice_base_elevation