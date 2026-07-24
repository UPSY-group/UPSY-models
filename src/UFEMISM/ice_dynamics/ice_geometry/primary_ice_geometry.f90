submodule( ice_geometry_model_basic) primary_ice_geometry

contains

  subroutine set_Hi( self, Hi)

    ! In/output variables:
    class(type_ice_geometry_model),                   intent(inout) :: self
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: Hi

    ! Local variables:
    character(len=*), parameter :: routine_name = 'set_Hi'
    integer                     :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = self%mesh%vi1, self%mesh%vi2
      self%Hi( vi) = Hi( vi)
    end do

    call self%mark_all_secondary_fields_as_no_longer_uptodate()

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_Hi

  subroutine set_Hb( self, Hb)

    ! In/output variables:
    class(type_ice_geometry_model),                   intent(inout) :: self
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: Hb

    ! Local variables:
    character(len=*), parameter :: routine_name = 'set_Hb'
    integer                     :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = self%mesh%vi1, self%mesh%vi2
      self%Hb( vi) = Hb( vi)
    end do

    call self%mark_all_secondary_fields_as_no_longer_uptodate()

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_Hb

  subroutine set_SL( self, SL)

    ! In/output variables:
    class(type_ice_geometry_model),                   intent(inout) :: self
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: SL

    ! Local variables:
    character(len=*), parameter :: routine_name = 'set_SL'
    integer                     :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = self%mesh%vi1, self%mesh%vi2
      self%SL( vi) = SL( vi)
    end do

    call self%mark_all_secondary_fields_as_no_longer_uptodate()

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_SL

  subroutine mark_all_secondary_fields_as_no_longer_uptodate( self)

    ! In/output variables:
    class(type_ice_geometry_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'mark_all_secondary_fields_as_no_longer_uptodate'

    ! Add routine to path
    call init_routine( routine_name)

    self%is_uptodate_Hs                 = .false.
    self%is_uptodate_Hib                = .false.
    self%is_uptodate_TAF                = .false.
    self%is_uptodate_Hi_eff             = .false.
    self%is_uptodate_Hs_slope           = .false.
    self%is_uptodate_Ho                 = .false.
    self%is_uptodate_ice_base_slopes    = .false.
    self%is_uptodate_grounded_fractions = .false.
    self%is_uptodate_masks              = .false.

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine mark_all_secondary_fields_as_no_longer_uptodate

end submodule primary_ice_geometry