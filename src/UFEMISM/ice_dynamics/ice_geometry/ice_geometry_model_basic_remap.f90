submodule(ice_geometry_model_basic) ice_geometry_model_basic_remap

contains

  subroutine remap_ice_geometry_model( self, mesh_old, mesh_new, refgeo_PD, GIA, forcing, time)
    !< Remap the basic ice geometry Hi,Hb,Hs,SL.

    ! In/output variables:
    class(type_ice_geometry_model),       intent(inout) :: self
    type(type_mesh),                      intent(in   ) :: mesh_old
    type(type_mesh),                      intent(in   ) :: mesh_new
    type(type_reference_geometry),        intent(in   ) :: refgeo_PD
    type(type_GIA_model),                 intent(in   ) :: GIA
    type(type_global_forcing),            intent(in   ) :: forcing
    real(dp),                             intent(in   ) :: time

    ! Local variables:
    character(len=*), parameter                     :: routine_name = 'remap_ice_geometry_model'
    real(dp), dimension( mesh_old%nV)               :: Hi_old_tot
    logical,  dimension( mesh_old%nV)               :: mask_floating_ice_tot
    logical,  dimension( mesh_old%nV)               :: mask_icefree_ocean_tot
    real(dp), dimension( mesh_new%vi1:mesh_new%vi2) :: Hi_new
    real(dp), dimension( mesh_new%vi1:mesh_new%vi2) :: Hs_new
    logical,  dimension( mesh_new%vi1:mesh_new%vi2) :: mask_noice_new
    integer                                         :: vi
    integer                                         :: mi, mi_used
    logical                                         :: found_map
    type(type_CSR_matrix_dp)                        :: M_CSR
    integer                                         :: vi_new, k1, k2, k, vi_old
    integer                                         :: n_ice, n_nonice
    integer                                         :: n_shelf, n_open_ocean
    real(dp)                                        :: sum_Hi_shelf

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap stuff that is common to all models
    call self%remap_model( mesh_new)

    ! Remap stuff that is specific to ice_geometry models

    ! == Basic: remap surface elevation Hs from the old mesh, remap bedrock elevation Hb
    !    from its (presumably high-resolution) source grid, define remapped ice thickness
    !    as the difference between the two. As surface elevation is typically much smoother
    !    then ice thickness, remapping works much better.
    ! =====================================================================================

    ! Remap bedrock from the original high-resolution grid, and add the (very smooth) modelled deformation to it
    ! Remapping of Hb in the refgeo structure has already happened, only need to copy the data
    if (par%primary) call warning('GIA model isnt finished yet - need to include dHb in mesh update!')
    call reallocate_bounds( self%Hb, mesh_new%vi1, mesh_new%vi2)  ! [m] Bedrock elevation (w.r.t. PD sea level)
    self%Hb = refgeo_PD%Hb

    ! Remap sea level
    call reallocate_bounds( self%SL, mesh_new%vi1, mesh_new%vi2)  ! [m] Sea level (geoid) elevation (w.r.t. PD sea level)
    select case (C%choice_sealevel_model)
    case default
      call crash('unknown choice_sealevel_model "' // trim( C%choice_sealevel_model) // '"')
    case ('fixed')
      self%SL = C%fixed_sealevel
    case ('prescribed')
      call update_sealevel_in_model(forcing, mesh_new, self, time)
    end select

    ! Gather global ice thickness and masks
    call gather_to_all( self%Hi, Hi_old_tot)
    call gather_to_all( self%mask_floating_ice , mask_floating_ice_tot )
    call gather_to_all( self%mask_icefree_ocean, mask_icefree_ocean_tot)

    ! First, naively remap ice thickness and surface elevation without any restrictions
    call map_from_mesh_to_mesh_2D( mesh_old, mesh_new, C%output_dir, self%Hi, Hi_new, '2nd_order_conservative')
    call map_from_mesh_to_mesh_2D( mesh_old, mesh_new, C%output_dir, self%Hs, Hs_new, '2nd_order_conservative')

    ! Calculate remapped ice thickness as the difference between new bedrock and remapped surface elevation
    do vi = mesh_new%vi1, mesh_new%vi2
      if (Hi_new( vi) > 0._dp) then
        if (Hs_new( vi) <= self%Hb( vi)) then
          Hi_new( vi) = 0._dp
        else
          Hi_new( vi) = Hi_from_Hb_Hs_and_SL( self%Hb( vi), Hs_new( vi), self%SL( vi))
        end if
      else
        Hi_new( vi) = 0._dp
      end if
    end do ! do vi = mesh_new%vi1, mesh_new%vi2

    ! Apply boundary conditions at the domain border
    call calc_mask_noice( mesh_new, mask_noice_new)
    call apply_ice_thickness_BC_explicit( mesh_new, mask_noice_new, self%Hb, self%SL, Hi_new)

    ! == Corrections
    ! ==============

    ! Browse the Atlas to find the map between mesh_old and mesh_new
    found_map = .false.
    mi_used   = 0
    do mi = 1, size( Atlas, 1)
      if (Atlas( mi)%name_src == mesh_old%name .and. Atlas( mi)%name_dst == mesh_new%name &
          .and. Atlas( mi)%method == '2nd_order_conservative') then
        found_map = .true.
        mi_used  = mi
        exit
      end if
    end do
    ! Safety
    if (.not. found_map) call crash('couldnt find which map was used')

    ! Convert the mapping matrix to CSR format
    call mat_petsc2CSR( Atlas( mi_used)%M, M_CSR)

    ! == For those vertices of the new mesh that overlap with both old-mesh ice and old-mesh
    !    non-ice, remove very thin remapped ice
    ! ======================================================================================

    do vi_new = mesh_new%vi1, mesh_new%vi2

      k1 = M_CSR%ptr( vi_new)
      k2 = M_CSR%ptr( vi_new+1) - 1

      n_ice    = 0
      n_nonice = 0

      do k = k1, k2

        vi_old = M_CSR%ind( k)

        if     (Hi_old_tot( vi_old) > 1._dp) then
          n_ice = n_ice + 1
        elseif (Hi_old_tot( vi_old) < 1._dp) then
          n_nonice = n_nonice + 1
        end if

      end do ! do k = k1, k2

      if (n_ice > 0 .and. n_nonice > 0) then
        ! This new-mesh vertex overlaps with both old-mesh ice vertices,
        ! and old-mesh non-ice vertices
        if (Hi_new( vi_new) < 1._dp) then
          ! Remove very thin remapped ice

          Hi_new( vi_new) = 0._dp
        end if
      end if

    end do ! do vi_new = mesh_new%vi1, mesh_new%vi2

    ! == For those vertices of the new mesh that overlap with both old-mesh shelf and old-mesh
    !    open ocean, average only over the contributing old-mesh shelf vertices
    ! ======================================================================================

    do vi_new = mesh_new%vi1, mesh_new%vi2

      k1 = M_CSR%ptr( vi_new)
      k2 = M_CSR%ptr( vi_new+1) - 1

      n_shelf      = 0
      n_open_ocean = 0
      sum_Hi_shelf = 0._dp

      do k = k1, k2

        vi_old = M_CSR%ind( k)

        if     (mask_floating_ice_tot(  vi_old)) then
          n_shelf      = n_shelf      + 1
          sum_Hi_shelf = sum_Hi_shelf + Hi_old_tot( vi_old)
        elseif (mask_icefree_ocean_tot( vi_old)) then
          n_open_ocean = n_open_ocean + 1
        end if

      end do

      if (n_shelf > 0 .and. n_open_ocean > 0) then
        ! This new-mesh vertex overlaps with both old-mesh shelf vertices,
        ! and old-mesh open-ocean vertices
        Hi_new( vi_new) = sum_Hi_shelf / real( n_shelf,dp)
      end if

    end do

    ! Recalculate Hs
    call reallocate_bounds( self%Hs, mesh_new%vi1, mesh_new%vi2)
    do vi = mesh_new%vi1, mesh_new%vi2
      self%Hs( vi) = ice_surface_elevation( Hi_new( vi), self%Hb( vi), self%SL( vi))
    end do

    ! Move Hi_new to geom%Hi
    deallocate( self%Hi)
    allocate( self%Hi( mesh_new%vi1: mesh_new%vi2))
    self%Hi = Hi_new

    call apply_mask_noice_direct( mesh_new, mask_noice_new, self%Hi)

    call reallocate_secondary_geometry_variables( self, mesh_new)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_ice_geometry_model

  subroutine reallocate_secondary_geometry_variables( self, mesh_new)

    ! In/output variables:
    class(type_ice_geometry_model), intent(inout) :: self
    type(type_mesh),                intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'reallocate_secondary_geometry_variables'

    ! Add routine to path
    call init_routine( routine_name)

    ! Secondary ice geometry fields
   !call reallocate_bounds( self%Hs                     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%Hib                    , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%TAF                    , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%Hi_eff                 , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%Hs_slope               , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%Ho                     , mesh_new%vi1, mesh_new%vi2)

    ! Horizontal derivatives
    call reallocate_bounds( self%dHib_dx_b              , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%dHib_dy_b              , mesh_new%ti1, mesh_new%ti2)

    ! Sub-grid bedrock cumulative density functions (CDFs)
    call reallocate_bounds( self%bedrock_cdf            , mesh_new%vi1, mesh_new%vi2, C%subgrid_bedrock_cdf_nbins)
    call reallocate_bounds( self%bedrock_cdf_b          , mesh_new%ti1, mesh_new%ti2, C%subgrid_bedrock_cdf_nbins)

    ! Area fractions
    call reallocate_bounds( self%fraction_gr            , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%fraction_gr_b          , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%fraction_margin        , mesh_new%vi1, mesh_new%vi2)

    ! Ice masks
    call reallocate_bounds( self%mask_icefree_land      , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%mask_icefree_ocean     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%mask_grounded_ice      , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%mask_floating_ice      , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%mask_margin            , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%mask_gl_gr             , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%mask_gl_fl             , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%mask_cf_gr             , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%mask_cf_fl             , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%mask_coastline         , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%mask                   , mesh_new%vi1, mesh_new%vi2)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine reallocate_secondary_geometry_variables

end submodule ice_geometry_model_basic_remap