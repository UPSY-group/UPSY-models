submodule( ice_geometry_model_basic) ice_base_slopes

contains

  subroutine calc_ice_base_slopes( self)

    ! In- and output variables
    class(type_ice_geometry_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_ice_base_slopes'

    ! Add routine to path
    call init_routine( routine_name)

    call ddx_a_b_2D( self%mesh, self%Hib, self%dHib_dx_b)
    call ddy_a_b_2D( self%mesh, self%Hib, self%dHib_dy_b)

    call checksum( self%mesh%pai_Tri, self%dHib_dx_b, 'self%dHib_dx_b')
    call checksum( self%mesh%pai_Tri, self%dHib_dy_b, 'self%dHib_dy_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_ice_base_slopes

end submodule ice_base_slopes