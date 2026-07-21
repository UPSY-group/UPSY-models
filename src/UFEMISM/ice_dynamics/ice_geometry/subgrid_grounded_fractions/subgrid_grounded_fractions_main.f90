submodule( ice_geometry_model_basic) subgrid_grounded_fractions_main
  !< Routines for calculating sub-grid grounded fractions

  use subgrid_grounded_fractions_bedrock_CDF, only: &
    calc_grounded_fractions_bedrock_CDF_a, &
    calc_grounded_fractions_bedrock_CDF_b
  use subgrid_grounded_fractions_bilin_TAF, only: &
    calc_grounded_fractions_bilin_interp_TAF_a, &
    calc_grounded_fractions_bilin_interp_TAF_b

contains

  subroutine calc_grounded_fractions( self, dHb, fraction_gr, fraction_gr_b)
    !< Calculate the sub-grid grounded-area fractions

    ! In- and output variables
    class(type_ice_geometry_model),                   intent(in   ) :: self
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: dHb
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: fraction_gr
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(  out) :: fraction_gr_b

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'calc_grounded_fractions'
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: TAF
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: fraction_gr_TAF_a
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: fraction_gr_CDF_a
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2) :: fraction_gr_TAF_b
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2) :: fraction_gr_CDF_b
    logical,  dimension(self%mesh%nV)                :: mask_floating_ice_tot
    integer                                          :: ti, via, vib, vic, vi


    ! Add routine to path
    call init_routine( routine_name)

    do vi = self%mesh%vi1, self%mesh%vi2
      TAF(vi) = thickness_above_floatation( self%Hi( vi), self%Hb( vi), self%SL( vi))
    end do

    ! Use the specified way of calculating sub-grid grounded fractions
    select case (C%choice_subgrid_grounded_fraction)
    case default
      call crash('unknown choice_subgrid_grounded_fraction "' // &
        trim( C%choice_subgrid_grounded_fraction) // '"')
    case('bilin_interp_TAF')
      ! Bilinearly interpolate the thickness above floatation to calculate the grounded fractions

      call calc_grounded_fractions_bilin_interp_TAF_a( self%mesh, TAF, fraction_gr_TAF_a)
      call calc_grounded_fractions_bilin_interp_TAF_b( self%mesh, TAF, fraction_gr_TAF_b)

      fraction_gr   = fraction_gr_TAF_a
      fraction_gr_b = fraction_gr_TAF_b

    case ('bedrock_CDF')
      ! Use the sub-grid bedrock cumulative density functions to calculate the grounded fractions

      call calc_grounded_fractions_bedrock_CDF_a( self%mesh, self%Hi, self%SL, dHb,      self%bedrock_cdf  , fraction_gr_CDF_a)
      call calc_grounded_fractions_bedrock_CDF_b( self%mesh, self%Hi, self%SL, dHb, TAF, self%bedrock_cdf_b, fraction_gr_CDF_b)

      fraction_gr   = fraction_gr_CDF_a
      fraction_gr_b = fraction_gr_CDF_b

    case ('bilin_interp_TAF+bedrock_CDF')
      ! Use the TAF method at the grounding line, and the CDF method inland

      call calc_grounded_fractions_bilin_interp_TAF_a( self%mesh, TAF, fraction_gr_TAF_a)
      call calc_grounded_fractions_bilin_interp_TAF_b( self%mesh, TAF, fraction_gr_TAF_b)

      call calc_grounded_fractions_bedrock_CDF_a( self%mesh, self%Hi, self%SL, dHb,      self%bedrock_cdf  , fraction_gr_CDF_a)
      call calc_grounded_fractions_bedrock_CDF_b( self%mesh, self%Hi, self%SL, dHb, TAF, self%bedrock_cdf_b, fraction_gr_CDF_b)

      ! Gather global floating ice mask
      call gather_to_all( self%mask_floating_ice, mask_floating_ice_tot)

      ! a-grid (vertices): take the smallest value (used for basal melt?)
      fraction_gr = min( fraction_gr_TAF_a, fraction_gr_CDF_a)

      ! b-grid (triangles): take CDF inland, TAF at grounding line (used for basal friction)
      do ti = self%mesh%ti1, self%mesh%ti2

        ! The three vertices spanning triangle ti
        via = self%mesh%Tri( ti,1)
        vib = self%mesh%Tri( ti,2)
        vic = self%mesh%Tri( ti,3)

        if (mask_floating_ice_tot( via) .OR. mask_floating_ice_tot( vib) .OR. mask_floating_ice_tot( vic)) then
          ! At least one corner of this triangle is afloat; grounding line
          fraction_gr_b( ti) = fraction_gr_TAF_b( ti)
        else
          ! All three corners of the triangle are grounded: inland
          fraction_gr_b( ti) = fraction_gr_CDF_b( ti)
        end if

      end do

    end select

    call checksum( self%mesh%pai_V  , fraction_gr  , 'fraction_gr'  )
    call checksum( self%mesh%pai_Tri, fraction_gr_b, 'fraction_gr_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions

end submodule subgrid_grounded_fractions_main