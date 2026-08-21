submodule(ice_velocity_model_basic) vertical_velocities

contains

  subroutine calc_vertical_velocities( self, ice, geom, BMB)
    !< Calculate vertical velocities w from conservation of mass

    ! NOTE: since the vertical velocities for floating ice depend on
    !       the thinning rate dH/dt, this routine must be called
    !       after having calculated dHi_dt!
    !
    ! Derivation:
    !
    ! Conservation of mass, combined with the incompressibility
    ! condition (i.e. constant density) of ice, is described by:
    !
    !   du/dx + dv/dy + dw/dz = 0
    !
    ! Applying the zeta coordinate transformation yields:
    !
    !   du/dxp + dzeta/dx du/dzeta + dv/dxp + dzeta/dy dv/dzeta + dzeta/dz dw/dzeta = 0
    !
    ! The terms du/dxp + dv/dyp describe the two-dimensional divergence in scaled coordinates:
    !
    !   grad uv = du/dxp + dv/dyp
    !
    ! The average value over a single grid cell (Voronoi cell) of this divergence is:
    !
    !   grad uv = intint_Voronoi (grad uv) dA / intint dA = 1/A intint_Voronoi (grad uv) dA
    !
    ! By applying the divergence theorem, the surface integral over the Voronoi cell
    ! can be transformed into a loop integral over the boundary of that Voronoi cell:
    !
    !   grad uv = 1/A cint (uv * n_hat) dS
    !
    ! Here, n_hat is the outward unit normal to the Voronoi cell boundary. Substituting
    ! this into the equation for conservation of mass yields:
    !
    !   dw/dzeta = -1 / dzeta/dz [ 1/A cint (uv * n_hat) dS + dzeta/dx du/zeta + dzeta/dy dv/dzeta]
    !
    ! The vertical velocity w at the ice base is equal to the horizontal motion along
    ! the sloping ice base, plus the vertical motion of the ice base itself, plus the
    ! vertical motion of an ice particle with respect to the ice base (i.e. the basal melt rate):
    !
    !   w( z=b) = u( z=b) * dH_base/dx + v( z=b) * dH_base/dy + dH_base/dt + M_base
    !
    ! With this boundary condition, dw/dzeta can be integrated over zeta to yield w( z).

    ! In- and output variables:
    class(atype_ice_velocity_model),                  intent(inout) :: self
    class(atype_ice_model_data),                      intent(in   ) :: ice
    class(atype_ice_geometry_model_data),             intent(in   ) :: geom
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: BMB

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'calc_vertical_velocities'
    integer                               :: vi,ks,ci,vj,ei
    real(dp), dimension(:  ), allocatable :: dHib_dx
    real(dp), dimension(:  ), allocatable :: dHib_dy
    real(dp), dimension(:  ), allocatable :: dHib_dt
    real(dp)                              :: dzeta
    real(dp), dimension(:,:), allocatable :: u_3D_c, u_3D_c_tot
    real(dp), dimension(:,:), allocatable :: v_3D_c, v_3D_c_tot
    real(dp)                              :: cint_un_dS, dS, u_ks, v_ks, un_dS, grad_uv_ks
    real(dp), dimension(2)                :: n_hat
    real(dp)                              :: du_dzeta_ks, dv_dzeta_ks
    real(dp)                              :: dzeta_dx_ks, dzeta_dy_ks, dzeta_dz_ks
    real(dp)                              :: dw_dzeta_ks

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate shared memory
    allocate( dHib_dx   ( self%mesh%vi1:self%mesh%vi2              ))
    allocate( dHib_dy   ( self%mesh%vi1:self%mesh%vi2              ))
    allocate( dHib_dt   ( self%mesh%vi1:self%mesh%vi2              ))
    allocate( u_3D_c    ( self%mesh%ei1:self%mesh%ei2, self%mesh%nz))
    allocate( v_3D_c    ( self%mesh%ei1:self%mesh%ei2, self%mesh%nz))
    allocate( u_3D_c_tot( self%mesh%nE,                self%mesh%nz))
    allocate( v_3D_c_tot( self%mesh%nE,                self%mesh%nz))

    do vi = self%mesh%vi1, self%mesh%vi2

      ! Calculate rate of change of ice base elevation
      if     (geom%mask_grounded_ice( vi)) then
        ! For grounded ice, the ice base simply moves with the bedrock
        dHib_dt( vi) =  ice%dHb_dt( vi)
      elseif (geom%mask_floating_ice( vi)) then
        ! For floating ice, the ice base moves according to the thinning rate times the density fraction
        dHib_dt( vi) = -ice%dHi_dt( vi) * ice_density / seawater_density
      else
        ! No ice, so no vertical velocity
        dHib_dt( vi) = 0._dp
      end if

    end do

    ! Calculate slopes of the ice base
    call ddx_a_a_2D( self%mesh, geom%Hib, dHib_dx)
    call ddy_a_a_2D( self%mesh, geom%Hib, dHib_dy)

    ! Calculate u,v on the c-grid (edges)
    call map_velocities_from_b_to_c_3D( self%mesh, self%u_3D_b, self%v_3D_b, u_3D_c, v_3D_c)
    call gather_to_all( u_3D_c, u_3D_c_tot)
    call gather_to_all( v_3D_c, v_3D_c_tot)

    ! Calculate vertical velocities by solving conservation of mass in each 3-D cell
    do vi = self%mesh%vi1, self%mesh%vi2

      ! No ice means no velocity
      if (.not. (geom%mask_grounded_ice( vi) .or. geom%mask_floating_ice( vi))) then
        self%w_3D( vi,:) = 0._dp
        cycle
      end if

      ! Calculate the vertical velocity at the ice base
      !
      ! NOTE: BMB is defined so that a positive number means accumulation of ice;
      !       at the ice base, that means that a positive BMB means a positive
      !       value of w

      if (geom%mask_floating_ice( vi)) then

        self%w_3D( vi,C%nz) = (self%u_3D( vi,C%nz) * dHib_dx( vi)) + &
                              (self%v_3D( vi,C%nz) * dHib_dy( vi)) + &
                              dHib_dt( vi) + MIN( 0._dp, BMB( vi))

      else

        self%w_3D( vi,C%nz) = (self%u_3D( vi,C%nz) * dHib_dx( vi)) + &
                              (self%v_3D( vi,C%nz) * dHib_dy( vi)) + &
                              dHib_dt( vi) + MIN( 0._dp, BMB( vi))

      end if


      ! Exception for very thin ice / ice margin: assume horizontal stretching
      ! is negligible, so that w( z) = w( z = b)
      if (geom%Hi( vi) < 10._dp) then
        self%w_3D( vi,:) = self%w_3D( vi,C%nz)
        cycle
      end if

      ! Calculate vertical velocities by integrating dw/dz over the vertical column

      do ks = self%mesh%nz-1, 1, -1

        dzeta = self%mesh%zeta( ks+1) - self%mesh%zeta( ks)

        ! Integrate u*n_hat around the Voronoi cell boundary
        cint_un_dS = 0._dp
        do ci = 1, self%mesh%nC( vi)
          vj = self%mesh%C(  vi,ci)
          ei = self%mesh%VE( vi,ci)
          ! Velocities at this section of the boundary
          u_ks = 0.5_dp * (u_3D_c_tot( ei,ks) + u_3D_c_tot( ei,ks+1))
          v_ks = 0.5_dp * (v_3D_c_tot( ei,ks) + v_3D_c_tot( ei,ks+1))
          ! Length of this section of the boundary
          dS = self%mesh%Cw( vi,ci)
          ! Outward normal vector to this section of the boundary
          n_hat = self%mesh%V( vj,:) - self%mesh%V( vi,:)
          n_hat = n_hat / NORM2( n_hat)
          ! Line integral over this section of the boundary
          un_dS = (u_ks * n_hat( 1) + v_ks * n_hat( 2)) * dS
          ! Add to loop integral
          cint_un_dS = cint_un_dS + un_dS
        end do

        ! Calculate grad uv from the divergence theorem
        grad_uv_ks = cint_un_dS / self%mesh%A( vi)

        ! Calculate du/dzeta, dv/dzeta
        du_dzeta_ks = (self%u_3D( vi,ks+1) - self%u_3D( vi,ks)) / dzeta
        dv_dzeta_ks = (self%v_3D( vi,ks+1) - self%v_3D( vi,ks)) / dzeta

        ! Calculate dzeta/dx, dzeta/dy, dzeta/dz
        dzeta_dx_ks = 0.5_dp * (ice%dzeta_dx_ak( vi,ks) + ice%dzeta_dx_ak( vi,ks+1))
        dzeta_dy_ks = 0.5_dp * (ice%dzeta_dy_ak( vi,ks) + ice%dzeta_dy_ak( vi,ks+1))
        dzeta_dz_ks = 0.5_dp * (ice%dzeta_dz_ak( vi,ks) + ice%dzeta_dz_ak( vi,ks+1))

        ! Calculate dw/dzeta
        dw_dzeta_ks = -1._dp / dzeta_dz_ks * (grad_uv_ks + dzeta_dx_ks * du_dzeta_ks + dzeta_dy_ks * dv_dzeta_ks)

        ! Calculate w
        self%w_3D( vi,ks) = self%w_3D( vi,ks+1) - dzeta * dw_dzeta_ks

      end do

    end do

    ! Also calculate dw/dz (inexpensive, no need to allow turning this off)
    call self%calc_dw_dz( geom)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_vertical_velocities

  subroutine calc_dw_dz( self, geom)

    ! In- and output variables:
    class(atype_ice_velocity_model),      intent(inout) :: self
    class(atype_ice_geometry_model_data), intent(in   ) :: geom

    ! Local variables:
    character(len=*), parameter       :: routine_name = 'calc_dw_dz'
    integer                           :: vi
    real(dp), dimension(self%mesh%nz) :: w_prof, dw_dzeta_prof

    ! Add routine to path
    call init_routine( routine_name)

    do vi = self%mesh%vi1, self%mesh%vi2

      ! No ice means undefined velocity
      if (.not. (geom%mask_grounded_ice( vi) .or. geom%mask_floating_ice( vi))) then
        self%dw_dz_3D( vi,:) = NaN
        cycle
      end if

      ! Calculate dw/dzeta in the vertical column, on the regular zeta-grid (use two-sided differencing)
      w_prof = self%w_3D( vi,:)
      call multiply_CSR_matrix_with_vector_local( self%mesh%M_ddzeta_k_k_1D, w_prof, dw_dzeta_prof)

      ! Calculate dw/dz from dw/dzeta (see Eqs. C2-C3 in Berends et al., 2024:
      ! Berends, C. J., van de Wal, R. S. W., and Zegeling, P. A.: Improvements
      ! on the discretisation of boundary conditions to the momentum balance
      ! for glacial ice, Journal of Glaciology, 1-15, doi: 10.1017/jog.2024.45, 2024)
      self%dw_dz_3D( vi,:) = -1._dp / geom%Hi( vi) * dw_dzeta_prof

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_dw_dz

end submodule vertical_velocities
