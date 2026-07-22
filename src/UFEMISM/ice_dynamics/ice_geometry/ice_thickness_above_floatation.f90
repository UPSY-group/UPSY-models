submodule( ice_geometry_model_basic) ice_thickness_above_floatation

contains

  subroutine calc_thickness_above_floatation( self)

    ! In- and output variables
    class(type_ice_geometry_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_thickness_above_floatation'
    integer                     :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = self%mesh%vi1, self%mesh%vi2
      self%TAF( vi) = thickness_above_floatation( self%Hi( vi), self%Hb( vi), self%SL( vi))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_thickness_above_floatation

end submodule ice_thickness_above_floatation