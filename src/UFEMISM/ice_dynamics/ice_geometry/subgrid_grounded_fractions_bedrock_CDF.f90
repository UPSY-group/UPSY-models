submodule(ice_geometry_model_basic) subgrid_grounded_fractions_bedrock_CDF

contains

  subroutine calc_grounded_fractions_bedrock_CDF_a( self, dHb, fraction_gr)
    !< Calculate the sub-grid grounded fractions of the vertices,
    !< using the sub-grid bedrock cumulative density functions (CDFs)

    ! In- and output variables
    class(type_ice_geometry_model),                   intent(in   ) :: self
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: dHb
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: fraction_gr

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_grounded_fractions_bedrock_CDF_a'
    integer                     :: vi, il, iu
    real(dp)                    :: Hb_float, wl, wu

    ! Add routine to path
    call init_routine( routine_name)

    do vi = self%mesh%vi1, self%mesh%vi2

      ! Compute the bedrock depth at which the current ice thickness and sea level
      ! will make this point afloat. Account for GIA here so we don't have to do it in
      ! the computation of the cumulative density function (CDF).

      Hb_float = self%SL( vi) - self%Hi( vi) * ice_density/seawater_density - dHb( vi)

      ! Get the fraction of bedrock within vertex coverage that is below
      ! Hb_float as a linear interpolation of the numbers in the CDF.

      if     (Hb_float <= self%bedrock_cdf( vi,1)) then
        ! All sub-grid points are above the floating bedrock elevation

        fraction_gr( vi) = 1._dp

      elseif (Hb_float >= self%bedrock_cdf( vi,C%subgrid_bedrock_cdf_nbins)) then
        ! All sub-grid points are below the floating bedrock elevation

        fraction_gr( vi) = 0._dp

      else

        ! Find the 2 elements in the CDF surrounding Hb_float
        iu = 1
        do while (self%bedrock_cdf( vi,iu) < Hb_float)
          iu = iu+1
        end do
        il = iu-1

        ! Interpolate the two enveloping bedrock bins to find the grounded fraction
        wl = (self%bedrock_cdf( vi,iu) - Hb_float) / (self%bedrock_cdf( vi,iu) - self%bedrock_cdf( vi,il))
        wu = 1._dp - wl
        fraction_gr( vi) = 1._dp - (real( il-1,dp) * wl + real( iu-1) * wu) / real( C%subgrid_bedrock_cdf_nbins-1,dp)

        ! Safety
        fraction_gr( vi) = min( 1._dp, max( 0._dp, fraction_gr( vi)))

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions_bedrock_CDF_a

  subroutine calc_grounded_fractions_bedrock_CDF_b( self, dHb, fraction_gr_b)
    !< Calculate the sub-grid grounded fractions of the triangles,
    !< using the sub-grid bedrock cumulative density functions (CDFs)

    ! In- and output variables
    class(type_ice_geometry_model),                   intent(in   ) :: self
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: dHb
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(  out) :: fraction_gr_b

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'calc_grounded_fractions_bedrock_CDF_b'
    real(dp), dimension(self%mesh%nV)                :: TAF_tot
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2) :: Hi_b, SL_b, dHb_b
    integer                                          :: ti, via, vib, vic, il, iu
    real(dp)                                         :: Hb_float, wl, wu

    ! Add routine to path
    call init_routine( routine_name)

    ! Map ice thickness, sea level, and bedrock deformation to the b-grid (triangles)
    call map_a_b_2D( self%mesh, self%Hi, Hi_b )
    call map_a_b_2D( self%mesh, self%SL, SL_b )
    call map_a_b_2D( self%mesh, dHb    , dHb_b)

    ! Gather global thickness above floatation
    call gather_dist_shared_to_all( self%mesh%pai_V, self%TAF, TAF_tot)

    do ti = self%mesh%ti1, self%mesh%ti2

      ! On the domain border, remapping issues make this answer unreliable
      ! (NOTE: only relevant when there's ice at the domain border, which in
      !        realistic experiments should never be the case; only happens
      !        in idealised geometries (e.g. MISMIP+))
      if (self%mesh%TriBI( ti) > 0) then
        ! if any of the three vertices spanning this triangle are grounded, treat it as grounded
        via = self%mesh%Tri( ti,1)
        vib = self%mesh%Tri( ti,2)
        vic = self%mesh%Tri( ti,3)
        if (TAF_tot( via) > 0._dp .or. TAF_tot( vib) > 0._dp .or. TAF_tot( vic) > 0._dp) then
          fraction_gr_b( ti) = 1._dp
        else
          fraction_gr_b( ti) = 0._dp
        end if
        cycle
      end if

      ! Compute the bedrock depth at which the current ice thickness and sea level
      ! will make this point afloat. Account for GIA here so we don't have to do it in
      ! the computation of the cumulative density function (CDF).

      Hb_float = SL_b( ti) - Hi_b( ti) * ice_density/seawater_density - dHb_b( ti)

      ! Get the fraction of bedrock within vertex coverage that is below
      ! Hb_float as a linear interpolation of the numbers in the CDF.

      if     (Hb_float <= self%bedrock_cdf_b( ti,1)) then
        ! All sub-grid points are above the floating bedrock elevation

        fraction_gr_b( ti) = 1._dp

      elseif (Hb_float >= self%bedrock_cdf_b( ti,C%subgrid_bedrock_cdf_nbins)) then
        ! All sub-grid points are below the floating bedrock elevation

        fraction_gr_b( ti) = 0._dp

      else

        ! Find the 2 elements in the CDF surrounding Hb_float
        iu = 1
        do while (self%bedrock_cdf_b( ti,iu) < Hb_float)
          iu = iu+1
        end do
        il = iu-1

        ! Interpolate the two enveloping bedrock bins to find the grounded fraction
        wl = (self%bedrock_cdf_b( ti,iu) - Hb_float) / (self%bedrock_cdf_b( ti,iu) - self%bedrock_cdf_b( ti,il))
        wu = 1._dp - wl
        fraction_gr_b( ti) = 1._dp - (real( il-1,dp) * wl + real( iu-1) * wu) / real( C%subgrid_bedrock_cdf_nbins-1,dp)

        ! Safety
        fraction_gr_b( ti) = min( 1._dp, max( 0._dp, fraction_gr_b( ti)))

      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions_bedrock_CDF_b

end submodule subgrid_grounded_fractions_bedrock_CDF