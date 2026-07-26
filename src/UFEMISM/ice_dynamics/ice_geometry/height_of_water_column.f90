submodule( ice_geometry_model_basic) height_of_water_column

contains

  subroutine calc_height_of_water_column( self)

    ! In- and output variables
    class(type_ice_geometry_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_height_of_water_column'
    integer                     :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = self%mesh%vi1, self%mesh%vi2
      self%Ho( vi) = height_of_water_column_at_ice_front( self%Hi( vi), self%Hb( vi), self%SL( vi))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_height_of_water_column

end submodule height_of_water_column