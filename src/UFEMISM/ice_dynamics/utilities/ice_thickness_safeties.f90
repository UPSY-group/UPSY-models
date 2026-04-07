module ice_thickness_safeties
  !< Different kinds of "safeties" to keep the ice sheet stable during nudging-based initialisation runs

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use parameters, only: ice_density, seawater_density
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use subgrid_ice_margin, only: calc_effective_thickness
  use ice_geometry_basics, only: is_floating
  use mpi_distributed_memory, only: gather_to_all
  use mpi_basic, only: par, sync
  use masks_mod, only: determine_masks
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_MIN, MPI_SUM, MPI_COMM_WORLD

  implicit none

  private

  public :: alter_ice_thickness, calc_and_apply_spill_over_flux

contains

  subroutine alter_ice_thickness( mesh, ice, Hi_old,Hb,SL, Hi_new, refgeo, time)
    !< Modify the predicted ice thickness in some sneaky way

    ! In- and output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hi_old
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hb
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: SL
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: Hi_new
    type(type_reference_geometry),          intent(in   ) :: refgeo
    real(dp),                               intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'alter_ice_thickness'
    integer                                :: vi
    real(dp)                               :: decay_start, decay_end
    real(dp)                               :: fixiness, limitness, fix_H_applied, limit_H_applied
    real(dp), dimension(mesh%vi1:mesh%vi2) :: modiness_up, modiness_down
    real(dp), dimension(mesh%vi1:mesh%vi2) :: Hi_save, Hi_eff_new, fraction_margin_new
    real(dp)                               :: floating_area, calving_area, mass_lost

    ! Add routine to path
    call init_routine( routine_name)

    ! Save predicted ice thickness for future reference
    Hi_save = Hi_new

    ! Calculate would-be effective thickness
    call calc_effective_thickness( mesh, Hi_new, Hb, SL, Hi_eff_new, fraction_margin_new)
    
    ! == Mask conservation
    ! ====================

    ! if desired, don't let grounded ice cross the floatation threshold
    if (C%do_protect_grounded_mask .and. time <= C%protect_grounded_mask_t_end) then
      do vi = mesh%vi1, mesh%vi2
        if (ice%mask_grounded_ice( vi)) then
          Hi_new( vi) = max( Hi_new( vi), (ice%SL( vi) - ice%Hb( vi)) * seawater_density/ice_density + .1_dp)
        end if
      end do
    end if

    ! General cases
    ! =============

    ! if so specified, remove very thin ice
    do vi = mesh%vi1, mesh%vi2
      if (Hi_eff_new( vi) < C%Hi_min) then
        Hi_new( vi) = 0._dp
      end if
    end do

    ! if so specified, remove thin floating ice
    if (C%choice_calving_law == 'threshold_thickness') then
      do vi = mesh%vi1, mesh%vi2
        if (is_floating( Hi_eff_new( vi), ice%Hb( vi), ice%SL( vi)) .and. Hi_eff_new( vi) < C%calving_threshold_thickness_shelf) then
          Hi_new( vi) = 0._dp
        end if
      end do
    end if

    ! DENK DROM
    if (C%remove_ice_absent_at_PD) then
      do vi = mesh%vi1, mesh%vi2
        if (refgeo%Hi( vi) == 0._dp) then
          Hi_new( vi) = 0._dp
        end if
      end do
    end if

    ! if so specified, remove all floating ice
    if (C%do_remove_shelves) then
      do vi = mesh%vi1, mesh%vi2
        if (is_floating( Hi_eff_new( vi), ice%Hb( vi), ice%SL( vi))) then
          Hi_new( vi) = 0._dp
        end if
      end do
    end if

    ! if so specified, remove all floating ice beyond the present-day calving front
    if (C%remove_shelves_larger_than_PD) then
      do vi = mesh%vi1, mesh%vi2
        if (refgeo%Hi( vi) == 0._dp .and. refgeo%Hb( vi) < 0._dp) then
          Hi_new( vi) = 0._dp
        end if
      end do
    end if

    ! if so specified, remove all floating ice crossing the continental shelf edge
    if (C%continental_shelf_calving) then
      do vi = mesh%vi1, mesh%vi2
        if (refgeo%Hi( vi) == 0._dp .and. refgeo%Hb( vi) < C%continental_shelf_min_height) then
          Hi_new( vi) = 0._dp
        end if
      end do
    end if

    ! === Fixiness ===
    ! ================

    ! Intial value
    fixiness = 1._dp

    ! Make sure that the start and end times make sense
    decay_start = C%fixiness_t_start
    decay_end   = C%fixiness_t_end

    ! Compute decaying fixiness
    if (decay_start >= decay_end) then
      ! Fix interval makes no sense: ignore fixiness
      fixiness = 0._dp
    elseif (time <= decay_start) then
      ! Before fix interval: check chosen option
      if (C%do_fixiness_before_start) then
        fixiness = 1._dp
      else
        fixiness = 0._dp
      end if
    elseif (time >= decay_end) then
      ! After fix interval: remove any fix/delay
      fixiness = 0._dp
    else
      ! We're within the fix interval: fixiness decreases with time
      fixiness = 1._dp - (time - decay_start) / (decay_end - decay_start)
    end if

    ! Just in case
    fixiness = min( 1._dp, max( 0._dp, fixiness))

    ! === Limitness ===
    ! =================

    ! Intial value
    limitness = 1._dp

    ! Make sure that the start and end times make sense
    decay_start = C%limitness_t_start
    decay_end   = C%limitness_t_end

    ! Compute decaying limitness
    if (decay_start >= decay_end) then
      ! Limit interval makes no sense: ignore limitness
      limitness = 0._dp
    elseif (time <= decay_start) then
      ! Before limit interval: check chosen option
      if (C%do_limitness_before_start) then
        limitness = 1._dp
      else
        limitness = 0._dp
      end if
    elseif (time >= decay_end) then
      ! After limit interval: remove any limits
      limitness = 0._dp
    else
      ! Limitness decreases with time
      limitness = 1._dp - (time - decay_start) / (decay_end - decay_start)
    end if

    ! Just in case
    limitness = min( 1._dp, max( 0._dp, limitness))

    ! === Modifier ===
    ! ================

    ! Intial value
    modiness_up   = 0._dp
    modiness_down = 0._dp

    select case (C%modiness_H_style)
    case default
      call crash('unknown modiness_H_choice "' // trim( C%modiness_H_style) // '"')
    case ('none')
      modiness_up   = 0._dp
      modiness_down = 0._dp

    case ('Ti_hom')
      modiness_up   = 1._dp - exp(ice%Ti_hom/C%modiness_T_hom_ref)
      modiness_down = 1._dp - exp(ice%Ti_hom/C%modiness_T_hom_ref)

    case ('Ti_hom_up')
      modiness_up   = 1._dp - exp(ice%Ti_hom/C%modiness_T_hom_ref)
      modiness_down = 0._dp

    case ('Ti_hom_down')
      modiness_up   = 0._dp
      modiness_down = 1._dp - exp(ice%Ti_hom/C%modiness_T_hom_ref)

    case ('no_thick_inland')
      do vi = mesh%vi1, mesh%vi2
        if (ice%mask_grounded_ice( vi) .and. .not. ice%mask_gl_gr( vi)) then
          modiness_up( vi) = 1._dp
        end if
      end do
      modiness_down = 0._dp

    case ('no_thin_inland')
      do vi = mesh%vi1, mesh%vi2
        if (ice%mask_grounded_ice( vi) .and. .not. ice%mask_gl_gr( vi)) then
          modiness_down( vi) = 1._dp
        end if
      end do
      modiness_up = 0._dp

    end select

    ! Just in case
    modiness_up   = min( 1._dp, max( 0._dp, modiness_up  ))
    modiness_down = min( 1._dp, max( 0._dp, modiness_down))

    ! === Fix, delay, limit ===
    ! =========================

    do vi = mesh%vi1, mesh%vi2

      ! Initialise
      fix_H_applied   = 0._dp
      limit_H_applied = 0._dp

      if (    ice%mask_gl_gr( vi)) then
        fix_H_applied   = C%fixiness_H_gl_gr * fixiness
        limit_H_applied = C%limitness_H_gl_gr * limitness

      elseif (ice%mask_gl_fl( vi)) then
        fix_H_applied   = C%fixiness_H_gl_fl * fixiness
        limit_H_applied = C%limitness_H_gl_fl * limitness

      elseif (ice%mask_grounded_ice( vi)) then
        fix_H_applied   = C%fixiness_H_grounded * fixiness
        limit_H_applied = C%limitness_H_grounded * limitness

      elseif (ice%mask_floating_ice( vi)) then
        fix_H_applied   = C%fixiness_H_floating * fixiness
        limit_H_applied = C%limitness_H_floating * limitness

      elseif (ice%mask_icefree_land( vi)) then
        if (C%fixiness_H_freeland .and. fixiness > 0._dp) fix_H_applied = 1._dp
        limit_H_applied = C%limitness_H_grounded * limitness

      elseif (ice%mask_icefree_ocean( vi)) then
        if (C%fixiness_H_freeocean .and. fixiness > 0._dp) fix_H_applied = 1._dp
        limit_H_applied = C%limitness_H_floating * limitness
      else
        ! if we reached this point, vertex is neither grounded, floating, nor ice free. That's a problem
        call crash('vertex neither grounded, floating, nor ice-free?')
      end if

      ! Apply fixiness
      Hi_new( vi) = Hi_old( vi) * fix_H_applied + Hi_new( vi) * (1._dp - fix_H_applied)

      ! Apply limitness
      Hi_new( vi) = min( Hi_new( vi), refgeo%Hi( vi) + (1._dp - modiness_up(   vi)) * limit_H_applied &
                                                    + (1._dp - limitness         ) * (Hi_new( vi) - refgeo%Hi( vi)) )

      Hi_new( vi) = max( Hi_new( vi), refgeo%Hi( vi) - (1._dp - modiness_down( vi)) * limit_H_applied &
                                                    - (1._dp - limitness         ) * (refgeo%Hi( vi) - Hi_new( vi)) )

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine alter_ice_thickness

  subroutine calc_and_apply_spill_over_flux( mesh, ice, Hi_new, dt)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(inout) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: Hi_new
    real(dp),                               intent(in   ) :: dt

    ! Local variables:
    character(len=1024), parameter                       :: routine_name = 'calc_and_apply_spill_over_flux'
    integer                                              :: vi, ci, vj, cj, cm, ierr
    real(dp)                                             :: Q_max, Q_min, Q_dsttot, Q_srctot
    real(dp), dimension(mesh%nC_mem)                     :: weight        ! [m2/y] Perpendicular outflow to ocean
    real(dp), dimension(mesh%vi1: mesh%vi2, mesh%nC_mem) :: relweight     ! [0-1] Relative outflow weight
    real(dp), dimension(mesh%nV, mesh%nC_mem)            :: relweight_tot ! [0-1] Relative outflow weight
    real(dp), dimension(mesh%vi1: mesh%vi2)              :: Q_src, Q_dst  ! [m3/y] Source and sink spill fluxes
    real(dp), dimension(mesh%nV)                         :: Q_src_tot     ! [m3/y] Source spill flux
    real(dp)                                             :: w_eps = 1.e-2 ! [m2/y] Small value
    logical, dimension(mesh%nV)                          :: mask_icefree_ocean_tot
    real(dp)                                             :: Hi_ups, Hs_ups ! [m] Upstream ice values
    real(dp), dimension(mesh%nV)                         :: Hi_new_tot, Hb_tot

    ! Initialise
    ice%Qspill = 0._dp
    Q_src      = 0._dp
    Q_dst      = 0._dp
    relweight  = 0._dp

    call gather_to_all( ice%mask_icefree_ocean, mask_icefree_ocean_tot)
    call gather_to_all( Hi_new, Hi_new_tot)
    call gather_to_all( ice%Hb, Hb_tot)

    ! Compute spill flux source
    do vi = mesh%vi1, mesh%vi2
  
      ! Determine upstream ice thickness
      if (ice%mask_cf_fl( vi) .or. ice%mask_cf_gr( vi)) then

        ! Find connection with the strongest inflow into this cell
        cm = minloc(ice%u_perp( vi, :), dim=1)

        ! If there is no inflow at all, for example during initialisation,
        ! use effective thickness
        if (ice%u_perp( vi, cm) >= 0._dp) then
          Hi_ups = ice%Hi_eff( vi)
        end if

        ! Find vertex index of strongest inflow cell
        vj = mesh%C( vi, cm)

        ! Skip if upstream cell is empty
        if (Hi_new_tot( vj) == 0._dp) cycle

        ! Determine upstream thickness from strongest inflow cell
        Hi_ups = Hi_new_tot( vj)

      else

        ! Default: use effective ice thickness
        Hi_ups = ice%Hi_eff( vi)

      end if

      ! Only compute spill-over when ice thickness exceeds effective ice thickness
      if ((ice%mask_cf_gr( vi) .or. ice%mask_cf_fl( vi)) .and. Hi_new( vi) > Hi_ups) then

        ! Determine source flux of spill-over ice in m^3/yr
        Q_src( vi) = - (Hi_new( vi) - Hi_ups) * mesh%A( vi) / dt

        weight = 0._dp
  
        ! Determine weights of surrounding ocean cells
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi, ci)
    
          if (mask_icefree_ocean_tot( vj)) then
            ! Define weight by outflow perpendicular velocity into ocean cells.
            ! Add small value to avoid division by 0 if no outflow velocity enters
            ! any ocean cell. In that case, weights will be equally distributed
            ! over all neighbouring ocean cells.
            weight( ci) = max(0._dp, ice%u_perp( vi, ci)) + w_eps
          end if
        end do

        ! Set spill source to zero if no surrounding ocean cells available
        if (sum(weight) < w_eps) then
          Q_src( vi) = 0._dp
        else
          ! Define the relative weight of Q_src to Q_dst of the downstream ocean cells
          do ci = 1, mesh%nC( vi)
            relweight( vi, ci) = weight( ci) / sum(weight)
          end do
        end if

      end if

    end do

    ! Only apply distribution during corrector step
    call gather_to_all( relweight, relweight_tot)
    call gather_to_all( Q_src, Q_src_tot)

    ! Compute spill flux destination
    do vi = mesh%vi1, mesh%vi2

      ! Skip if not an ocean cell
      if (.not. mask_icefree_ocean_tot( vi)) cycle

        do ci = 1, mesh%nC( vi)
          !Neighbouring cell which potentially has nonzero Q_src
          vj = mesh%C( vi, ci)
  
          ! Connections of neighbouring cell
          do cj = 1, mesh%nC( vj)
  
            if (mesh%C( vj, cj) == vi) then
            ! Yes, this connection is a match, receive fraction of Q_src
              if (Q_src_tot( vj) < -1.e-2_dp .and. relweight_tot( vj, cj) > 1.e-6_dp) then
                Q_dst( vi) = Q_dst( vi) - Q_src_tot( vj) * relweight_tot( vj, cj)
              end if  

            end if
  
          end do

      end do

    end do

    ! Apply spill rates
    do vi = mesh%vi1, mesh%vi2

      ! Combine source and destination and convert to thickness rate in m/yr
      ice%Qspill( vi) = (Q_src( vi) + Q_dst( vi)) / mesh%A( vi)

      ! Update ice thickness
      Hi_new( vi) = Hi_new( vi) + ice%Qspill( vi) * dt

    end do

    ! Q_max = maxval( ice%Qspill)
    ! call MPI_ALLREDUCE( MPI_IN_PLACE, Q_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    ! Q_min = minval( ice%Qspill)
    ! call MPI_ALLREDUCE( MPI_IN_PLACE, Q_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    ! Q_dsttot = sum( Q_dst)*1e-9
    ! call MPI_ALLREDUCE( MPI_IN_PLACE, Q_dsttot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    ! Q_srctot = sum( Q_src)*1e-9
    ! call MPI_ALLREDUCE( MPI_IN_PLACE, Q_srctot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! if (par%primary) write (*,*) Q_max, Q_min, Q_dsttot, Q_srctot

    ! Update masks

    ! call determine_masks( mesh, Hi_new, ice%Hb, ice%SL, ice%mask, ice%mask_icefree_land, & 
    !                       ice%mask_icefree_ocean, ice%mask_grounded_ice, ice%mask_floating_ice, &
    !                       ice%mask_margin, ice%mask_gl_fl, ice%mask_gl_gr,ice%mask_cf_gr, &
    !                       ice%mask_cf_fl, ice%mask_coastline)

  end subroutine calc_and_apply_spill_over_flux

end module ice_thickness_safeties
