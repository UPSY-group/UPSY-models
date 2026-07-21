submodule( ice_geometry_model_basic) determine_masks_mod

contains

  subroutine determine_masks( self)
    !< Determine the different masks

    ! In/output variables:
    class(type_ice_geometry_model),intent(inout) :: self

    ! Local variables:
    character(len=*), parameter      :: routine_name = 'determine_masks'
    integer                          :: vi, ci, vj
    logical, dimension(self%mesh%nV) :: mask_icefree_land_tot
    logical, dimension(self%mesh%nV) :: mask_icefree_ocean_tot
    logical, dimension(self%mesh%nV) :: mask_grounded_ice_tot
    logical, dimension(self%mesh%nV) :: mask_floating_ice_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! === Basic masks ===
    ! ===================

    ! Initialise basic masks
    self%mask_icefree_land  = .false.
    self%mask_icefree_ocean = .false.
    self%mask_grounded_ice  = .false.
    self%mask_floating_ice  = .false.
    self%mask               = 0

    ! Calculate
    do vi = self%mesh%vi1, self%mesh%vi2

      if (is_floating( self%Hi( vi), self%Hb( vi), self%SL( vi))) then
        ! Ice thickness is below the floatation thickness; either floating ice, or ice-free ocean

        if (self%Hi( vi) > 0._dp) then
          ! Floating ice

          self%mask_floating_ice( vi) = .true.
          self%mask( vi) = C%type_floating_ice

        else
          ! Ice-free ocean

          self%mask_icefree_ocean( vi) = .true.
          self%mask( vi) = C%type_icefree_ocean

        end if

      else
        ! Ice thickness is above the floatation thickness; either grounded ice, or ice-free land

        if (self%Hi( vi) > 0._dp) then
          ! Grounded ice

          self%mask_grounded_ice( vi) = .true.
          self%mask( vi) = C%type_grounded_ice

        else
          ! Ice-free land

          self%mask_icefree_land( vi) = .true.
          self%mask( vi) = C%type_icefree_land

        end if

      end if

    end do

    call checksum( self%mesh%pai_V, self%mask_icefree_land , 'mask_icefree_land')
    call checksum( self%mesh%pai_V, self%mask_icefree_ocean, 'mask_icefree_ocean')
    call checksum( self%mesh%pai_V, self%mask_grounded_ice , 'mask_grounded_ice')
    call checksum( self%mesh%pai_V, self%mask_floating_ice , 'mask_floating_ice')
    call checksum( self%mesh%pai_V, self%mask              , 'mask')

    ! === Transitional masks ===
    ! ==========================

    ! Gather basic masks to all processes
    call gather_to_all( self%mask_icefree_land , mask_icefree_land_tot )
    call gather_to_all( self%mask_icefree_ocean, mask_icefree_ocean_tot)
    call gather_to_all( self%mask_grounded_ice , mask_grounded_ice_tot )
    call gather_to_all( self%mask_floating_ice , mask_floating_ice_tot )

    ! Initialise transitional masks
    self%mask_margin    = .false.
    self%mask_gl_gr     = .false.
    self%mask_gl_fl     = .false.
    self%mask_cf_gr     = .false.
    self%mask_cf_fl     = .false.
    self%mask_coastline = .false.

    do vi = self%mesh%vi1, self%mesh%vi2

      ! Ice margin
      if (mask_grounded_ice_tot( vi) .OR. mask_floating_ice_tot( vi)) then
        do ci = 1, self%mesh%nC( vi)
          vj = self%mesh%C( vi,ci)
          if (.not. (mask_grounded_ice_tot( vj) .OR. mask_floating_ice_tot( vj))) then
            self%mask_margin( vi) = .true.
            self%mask( vi) = C%type_margin
          end if
        end do
      end if

      ! Grounding line (grounded side)
      if (mask_grounded_ice_tot( vi)) then
        do ci = 1, self%mesh%nC( vi)
          vj = self%mesh%C( vi,ci)
          if (mask_floating_ice_tot( vj)) then
            self%mask_gl_gr( vi) = .true.
            self%mask( vi) = C%type_groundingline_gr
          end if
        end do
      end if

      ! Grounding line (floating side)
      if (mask_floating_ice_tot( vi)) then
        do ci = 1, self%mesh%nC( vi)
          vj = self%mesh%C( vi,ci)
          if (mask_grounded_ice_tot( vj)) then
            self%mask_gl_fl( vi) = .true.
            self%mask( vi) = C%type_groundingline_fl
          end if
        end do
      end if

      ! Calving front (grounded)
      if (mask_grounded_ice_tot( vi)) then
        do ci = 1, self%mesh%nC(vi)
          vj = self%mesh%C( vi,ci)
          if (mask_icefree_ocean_tot( vj)) then
            self%mask_cf_gr( vi) = .true.
            self%mask( vi) = C%type_calvingfront_gr
          end if
        end do
      end if

      ! Calving front (floating)
      if (mask_floating_ice_tot( vi)) then
        do ci = 1, self%mesh%nC(vi)
          vj = self%mesh%C( vi,ci)
          if (mask_icefree_ocean_tot( vj)) then
            self%mask_cf_fl( vi) = .true.
            self%mask( vi) = C%type_calvingfront_fl
          end if
        end do
      end if

      ! Coastline
      if (mask_icefree_land_tot( vi)) then
        do ci = 1, self%mesh%nC(vi)
          vj = self%mesh%C( vi,ci)
          if (mask_icefree_ocean_tot( vj)) then
            self%mask_coastline( vi) = .true.
            self%mask( vi) = C%type_coastline
          end if
        end do
      end if

    end do

    call checksum( self%mesh%pai_V, self%mask_margin   , 'mask_margin')
    call checksum( self%mesh%pai_V, self%mask_gl_gr    , 'mask_gl_gr')
    call checksum( self%mesh%pai_V, self%mask_gl_fl    , 'mask_gl_fl')
    call checksum( self%mesh%pai_V, self%mask_cf_gr    , 'mask_cf_gr')
    call checksum( self%mesh%pai_V, self%mask_cf_fl    , 'mask_cf_fl')
    call checksum( self%mesh%pai_V, self%mask_coastline, 'mask_coastline')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_masks

end submodule determine_masks_mod