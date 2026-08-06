submodule( ice_geometry_model_basic) subgrid_grounded_fractions_main
  !< Routines for calculating sub-grid grounded fractions

contains

  subroutine calc_grounded_fractions( self, dHb)
    !< Calculate the sub-grid grounded-area fractions

    ! In- and output variables
    class(type_ice_geometry_model),                   intent(inout) :: self
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: dHb

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'calc_grounded_fractions'
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: fraction_gr_TAF_a
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: fraction_gr_CDF_a
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2) :: fraction_gr_TAF_b
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2) :: fraction_gr_CDF_b
    logical,  dimension(self%mesh%nV)                :: mask_floating_ice_tot
    integer                                          :: ti, via, vib, vic

    ! Add routine to path
    call init_routine( routine_name)

    ! Use the specified way of calculating sub-grid grounded fractions
    select case (C%choice_subgrid_grounded_fraction)
    case default
      call crash('unknown choice_subgrid_grounded_fraction "' // &
        trim( C%choice_subgrid_grounded_fraction) // '"')
    case('bilin_interp_TAF')
      ! Bilinearly interpolate the thickness above floatation to calculate the grounded fractions

      call self%calc_grounded_fractions_bilin_interp_TAF_a( fraction_gr_TAF_a)
      call self%calc_grounded_fractions_bilin_interp_TAF_b( fraction_gr_TAF_b)

      self%fraction_gr   = fraction_gr_TAF_a
      self%fraction_gr_b = fraction_gr_TAF_b

    case ('bedrock_CDF')
      ! Use the sub-grid bedrock cumulative density functions to calculate the grounded fractions

      call self%calc_grounded_fractions_bedrock_CDF_a( dHb, fraction_gr_CDF_a)
      call self%calc_grounded_fractions_bedrock_CDF_b( dHb, fraction_gr_CDF_b)

      self%fraction_gr   = fraction_gr_CDF_a
      self%fraction_gr_b = fraction_gr_CDF_b

    case ('bilin_interp_TAF+bedrock_CDF')
      ! Use the TAF method at the grounding line, and the CDF method inland

      call self%calc_grounded_fractions_bilin_interp_TAF_a( fraction_gr_TAF_a)
      call self%calc_grounded_fractions_bilin_interp_TAF_b( fraction_gr_TAF_b)

      call self%calc_grounded_fractions_bedrock_CDF_a( dHb, fraction_gr_CDF_a)
      call self%calc_grounded_fractions_bedrock_CDF_b( dHb, fraction_gr_CDF_b)

      ! Gather global floating ice mask
      call gather_to_all( self%mask_floating_ice, mask_floating_ice_tot)

      ! a-grid (vertices): take the smallest value (used for basal melt?)
      self%fraction_gr = min( fraction_gr_TAF_a, fraction_gr_CDF_a)

      ! b-grid (triangles): take CDF inland, TAF at grounding line (used for basal friction)
      do ti = self%mesh%ti1, self%mesh%ti2

        ! The three vertices spanning triangle ti
        via = self%mesh%Tri( ti,1)
        vib = self%mesh%Tri( ti,2)
        vic = self%mesh%Tri( ti,3)

        if (mask_floating_ice_tot( via) .OR. mask_floating_ice_tot( vib) .OR. mask_floating_ice_tot( vic)) then
          ! At least one corner of this triangle is afloat; grounding line
          self%fraction_gr_b( ti) = fraction_gr_TAF_b( ti)
        else
          ! All three corners of the triangle are grounded: inland
          self%fraction_gr_b( ti) = fraction_gr_CDF_b( ti)
        end if

      end do

    end select

    call checksum( self%mesh%pai_V  , self%fraction_gr  , 'fraction_gr'  )
    call checksum( self%mesh%pai_Tri, self%fraction_gr_b, 'fraction_gr_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions

end submodule subgrid_grounded_fractions_main