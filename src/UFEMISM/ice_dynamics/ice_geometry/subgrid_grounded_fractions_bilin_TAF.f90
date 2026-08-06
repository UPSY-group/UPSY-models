submodule(ice_geometry_model_basic) subgrid_grounded_fractions_bilin_TAF

contains

  subroutine calc_grounded_fractions_bilin_interp_TAF_a( self, fraction_gr)
    !< Calculate the sub-grid grounded fractions of the vertices

    ! Bilinearly interpolate the thickness above floatation (the CISM/PISM approach)

    ! In- and output variables
    class(type_ice_geometry_model),                   intent(in   ) :: self
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: fraction_gr

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'calc_grounded_fractions_bilin_interp_TAF_a'
    real(dp), dimension(self%mesh%nV)                :: TAF_tot
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2) :: TAF_b
    real(dp), dimension(self%mesh%nTri)              :: TAF_b_tot
    integer                                          :: vi, ci, vj, iti1, iti2, ti1, ti2, iti, ti
    real(dp)                                         :: TAF_max, TAF_min
    real(dp), dimension(2)                           :: va, vb, vc
    real(dp)                                         :: TAFa, TAFb, TAFc, A_vor, A_tri_tot, A_tri_grnd, A_grnd

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather global thickness above floatation
    call gather_to_all( self%TAF, TAF_tot)

    ! Map thickness-above-floatation to the b-grid
    call map_a_b_2D(  self%mesh, self%TAF, TAF_b)

    ! Gather global thickness above floatation on the b-grid
    call gather_to_all( TAF_b, TAF_b_tot)

    do vi = self%mesh%vi1, self%mesh%vi2

      ! Determine maximum and minimum TAF of the local neighbourhood
      TAF_max = -1E6_dp
      TAF_min =  1E6_dp

      TAF_max = max( TAF_max, TAF_tot( vi))
      TAF_min = min( TAF_min, TAF_tot( vi))

      do ci = 1, self%mesh%nC( vi)
        vj = self%mesh%C( vi,ci)
        TAF_max = max( TAF_max, TAF_tot( vj))
        TAF_min = min( TAF_min, TAF_tot( vj))
      end do

      ! if the entire local neighbourhood is grounded, the answer is trivial
      if (TAF_min >= 0._dp) then
        fraction_gr( vi) = 1._dp
        cycle
      end if

      ! if the entire local neighbourhood is floating, the answer is trivial
      if (TAF_max <= 0._dp) then
        fraction_gr( vi) = 0._dp
        cycle
      end if

      ! The local neighbourhood contains both grounded and floating vertices.
      A_vor  = 0._dp
      A_grnd = 0._dp

      va   = self%mesh%V( vi,:)
      TAFa = TAF_tot( vi)

      if (self%mesh%VBI( vi) == 0) then
        ! Free vertex

        do iti1 = 1, self%mesh%niTri( vi)

          iti2 = iti1 + 1
          if (iti2 == self%mesh%niTri( vi) + 1) iti2 = 1

          ti1 = self%mesh%iTri( vi,iti1)
          ti2 = self%mesh%iTri( vi,iti2)

          vb = self%mesh%Tricc( ti1,:)
          vc = self%mesh%Tricc( ti2,:)

          TAFb = TAF_b_tot( ti1)
          TAFc = TAF_b_tot( ti2)

          ! Calculate total area of, and grounded area within, this subtriangle
          call calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

          A_vor  = A_vor  + A_tri_tot
          A_grnd = A_grnd + A_tri_grnd

        end do ! do vori1 = 1, nVor

      else
        ! Border vertex

        ! First subtriangle
        vj   = self%mesh%C( vi,1)
        vb   = 0.5_dp * (self%mesh%V(  vi,:) + self%mesh%V(  vj,:))
        TAFb = 0.5_dp * (TAF_tot( vi  ) + TAF_tot( vj  ))

        iti  = 1
        ti   = self%mesh%iTri( vi,iti)
        vc   = self%mesh%Tricc( ti,:)
        TAFc = TAF_b_tot( ti)

        ! Calculate total area of, and grounded area within, this subtriangle
        call calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

        A_vor  = A_vor  + A_tri_tot
        A_grnd = A_grnd + A_tri_grnd

        ! Middle subtriangles
        do iti1 = 1, self%mesh%niTri( vi)-1

          iti2 = iti1 + 1
          if (iti2 == self%mesh%niTri( vi) + 1) iti2 = 1

          ti1 = self%mesh%iTri( vi,iti1)
          ti2 = self%mesh%iTri( vi,iti2)

          vb = self%mesh%Tricc( ti1,:)
          vc = self%mesh%Tricc( ti2,:)

          TAFb = TAF_b_tot( ti1)
          TAFc = TAF_b_tot( ti2)

          ! Calculate total area of, and grounded area within, this subtriangle
          call calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

          A_vor  = A_vor  + A_tri_tot
          A_grnd = A_grnd + A_tri_grnd

        end do ! do vori1 = 1, nVor

        ! Last subtriangle
        iti  = self%mesh%niTri( vi)
        ti   = self%mesh%iTri( vi,iti)
        vb   = self%mesh%Tricc( ti,:)
        TAFb = TAF_b_tot( ti)

        vj   = self%mesh%C( vi, self%mesh%nC( vi))
        vc   = 0.5_dp * (self%mesh%V(  vi,:) + self%mesh%V(  vj,:))
        TAFc = 0.5_dp * (TAF_tot( vi  ) + TAF_tot( vj  ))

        ! Calculate total area of, and grounded area within, this subtriangle
        call calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

        A_vor  = A_vor  + A_tri_tot
        A_grnd = A_grnd + A_tri_grnd

      end if ! if (self%mesh%VBI( vi) == 0) then

      ! Calculate the sub-grid grounded fraction of this Voronoi cell
      fraction_gr( vi) = min( 1._dp, max( 0._dp, A_grnd / A_vor ))

    end do ! do vi = self%mesh%vi1, self%mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions_bilin_interp_TAF_a

  subroutine calc_grounded_fractions_bilin_interp_TAF_b( self, fraction_gr_b)
    !< Calculate the sub-grid grounded fractions of the triangles

    ! Bilinearly interpolate the thickness above floatation (the CISM/PISM approach)

    ! In- and output variables
    class(type_ice_geometry_model),                   intent(in   ) :: self
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(  out) :: fraction_gr_b

    ! Local variables:
    character(len=*), parameter       :: routine_name = 'calc_grounded_fractions_bilin_interp_TAF_b'
    real(dp), dimension(self%mesh%nV) :: TAF_tot
    integer                           :: ti, via, vib, vic
    real(dp)                          :: TAF_max, TAF_min
    real(dp), dimension(2)            :: va, vb, vc
    real(dp)                          :: TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather global thickness above floatation
    call gather_to_all( self%TAF, TAF_tot)

    do ti = self%mesh%ti1, self%mesh%ti2

      ! The indices of the three vertices spanning triangle ti
      via = self%mesh%Tri( ti,1)
      vib = self%mesh%Tri( ti,2)
      vic = self%mesh%Tri( ti,3)

      ! The coordinates of the three vertices spanning triangle ti
      va   = self%mesh%V( via,:)
      vb   = self%mesh%V( vib,:)
      vc   = self%mesh%V( vic,:)

      ! Thickness above floatation at the three corners of the triangle
      TAFa = TAF_tot( via)
      TAFb = TAF_tot( vib)
      TAFc = TAF_tot( vic)

      ! Determine maximum and minimum TAF of the local neighbourhood
      TAF_max = maxval([ TAFa, TAFb, TAFc])
      TAF_min = minval([ TAFa, TAFb, TAFc])

      ! if the entire local neighbourhood is grounded, the answer is trivial
      if (TAF_min >= 0._dp) then
        fraction_gr_b( ti) = 1._dp
        cycle
      end if

      ! if the entire local neighbourhood is floating, the answer is trivial
      if (TAF_max <= 0._dp) then
        fraction_gr_b( ti) = 0._dp
        cycle
      end if

      ! Calculate total area of, and grounded area within, this subtriangle
      call calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)

      ! Calculate the sub-grid grounded fraction of this Voronoi cell
      fraction_gr_b( ti) = min( 1._dp, max( 0._dp, A_tri_grnd / A_tri_tot ))

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions_bilin_interp_TAF_b

  subroutine calc_grounded_area_triangle( va, vb, vc, TAFa, TAFb, TAFc, A_tri_tot, A_tri_grnd)
    !< Calculate the grounded area of the triangle [va,vb,vc],
    !< where the thickness above floatation is given at all three corners

    ! In- and output variables
    real(dp), dimension(2), intent(in   ) :: va, vb, vc
    real(dp),               intent(in   ) :: TAFa, TAFb, TAFc
    real(dp),               intent(  out) :: A_tri_tot, A_tri_grnd

    ! Local variables:
    real(dp)            :: A_flt
    real(dp), parameter :: tol = 1E-9_dp
    real(dp)            :: TAFa_corr, TAFb_corr, TAFc_corr

    ! Determine total area of this subtriangle
    A_tri_tot = triangle_area( va, vb, vc)

    ! TAF of zero can cause problems, correct for this
    if (TAFa >= 0._dp) then
      TAFa_corr = max( tol, TAFa)
    else
      TAFa_corr = min( -tol, TAFa)
    end if
    if (TAFb >= 0._dp) then
      TAFb_corr = max( tol, TAFb)
    else
      TAFb_corr = min( -tol, TAFb)
    end if
    if (TAFc >= 0._dp) then
      TAFc_corr = max( tol, TAFc)
    else
      TAFc_corr = min( -tol, TAFc)
    end if


    if     (TAFa_corr >= 0._dp .and. TAFb_corr >= 0._dp .and. TAFc_corr >= 0._dp) then
      ! if all three corners are grounded, the answer is trivial
      A_tri_grnd = A_tri_tot
    elseif (TAFa_corr <= 0._dp .and. TAFb_corr <= 0._dp .and. TAFc_corr <= 0._dp) then
      ! if all three corners are floating, the answer is trivial
      A_tri_grnd = 0._dp
    else
      ! At least one corner is grounded and at least one corner is floating

      if     (TAFa_corr >= 0._dp .and. TAFb_corr <= 0._dp .and. TAFc_corr <= 0._dp) then
        ! a is grounded, b and c are floating
        call calc_grounded_area_triangle_1grnd_2flt( va, vb, vc, TAFa_corr, TAFb_corr, TAFc_corr, A_tri_grnd)
      elseif (TAFa_corr <= 0._dp .and. TAFb_corr >= 0._dp .and. TAFc_corr <= 0._dp) then
        ! b is grounded, a and c are floating
        call calc_grounded_area_triangle_1grnd_2flt( vb, vc, va, TAFb_corr, TAFc_corr, TAFa_corr, A_tri_grnd)
      elseif (TAFa_corr <= 0._dp .and. TAFb_corr <= 0._dp .and. TAFc_corr >= 0._dp) then
        ! c is grounded, a and b are floating
        call calc_grounded_area_triangle_1grnd_2flt( vc, va, vb, TAFc_corr, TAFa_corr, TAFb_corr, A_tri_grnd)
      elseif (TAFa_corr <= 0._dp .and. TAFb_corr >= 0._dp .and. TAFc_corr >= 0._dp) then
        ! a is floating, b and c are grounded
        call calc_grounded_area_triangle_1flt_2grnd( va, vb, vc, TAFa_corr, TAFb_corr, TAFc_corr, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      elseif (TAFa_corr >= 0._dp .and. TAFb_corr <= 0._dp .and. TAFc_corr >= 0._dp) then
        ! b is floating, c and a are grounded
        call calc_grounded_area_triangle_1flt_2grnd( vb, vc, va, TAFb_corr, TAFc_corr, TAFa_corr, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      elseif (TAFa_corr >= 0._dp .and. TAFb_corr >= 0._dp .and. TAFc_corr <= 0._dp) then
        ! c is floating, a and b are grounded
        call calc_grounded_area_triangle_1flt_2grnd( vc, va, vb, TAFc_corr, TAFa_corr, TAFb_corr, A_flt)
        A_tri_grnd = A_tri_tot - A_flt
      else
        A_tri_grnd = 0._dp
        call crash('TAF = [{dp_01},{dp_02},{dp_03}]', dp_01 = TAFa_corr, dp_02 = TAFb_corr, dp_03 = TAFc_corr)
      end if

    end if

  end subroutine calc_grounded_area_triangle

  subroutine calc_grounded_area_triangle_1grnd_2flt( va, vb, vc, TAFa, TAFb, TAFc, A_tri_grnd)
    !< Calculate the grounded area of the triangle [va,vb,vc], where vertex a is grounded
    !< and b and c are floating

    ! In- and output variables
    real(dp), dimension(2), intent(in   ) :: va, vb, vc
    real(dp),               intent(in   ) :: TAFa, TAFb, TAFc
    real(dp),               intent(  out) :: A_tri_grnd

    ! Local variables:
    real(dp)               :: lambda_ab, lambda_ac
    real(dp), dimension(2) :: pab, pac

    lambda_ab = TAFa / (TAFa - TAFb)
    pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)

    lambda_ac = TAFa / (TAFa - TAFc)
    pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)

    A_tri_grnd = triangle_area( va, pab, pac)

  end subroutine calc_grounded_area_triangle_1grnd_2flt

  subroutine calc_grounded_area_triangle_1flt_2grnd( va, vb, vc, TAFa, TAFb, TAFc, A_tri_flt)
    !< Calculate the grounded area of the triangle [va,vb,vc], where vertex a is floating
    !< and b and c are grounded

    ! In- and output variables
    real(dp), dimension(2), intent(in   ) :: va, vb, vc
    real(dp),               intent(in   ) :: TAFa, TAFb, TAFc
    real(dp),               intent(  out) :: A_tri_flt

    ! Local variables:
    real(dp)               :: lambda_ab, lambda_ac
    real(dp), dimension(2) :: pab, pac

    lambda_ab = TAFa / (TAFa - TAFb)
    pab = (va * (1._dp - lambda_ab)) + (vb * lambda_ab)

    lambda_ac = TAFa / (TAFa - TAFc)
    pac = (va * (1._dp - lambda_ac)) + (vc * lambda_ac)

    A_tri_flt = triangle_area( va, pab, pac)

  end subroutine calc_grounded_area_triangle_1flt_2grnd

end submodule subgrid_grounded_fractions_bilin_TAF