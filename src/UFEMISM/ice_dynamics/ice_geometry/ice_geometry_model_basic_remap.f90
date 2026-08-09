submodule(ice_geometry_model_basic) ice_geometry_model_basic_remap

contains

  subroutine remap_ice_geometry_model( self, mesh_old, mesh_new, refgeo_PD, GIA, forcing, time)
    ! Remap the surface elevation Hs conservatively from the old mesh, and remap
    ! the bedrock elevation Hb from its (presumably high-resolution) source grid
    ! in the reference geometry object. Then, define the remapped ice thickness
    ! as the difference between the two. This works well, because the surface
    ! elevation is typically much smoother than the ice thickness.

    ! In/output variables:
    class(type_ice_geometry_model), intent(inout) :: self
    type(type_mesh),                intent(in   ) :: mesh_old
    type(type_mesh),                intent(in   ) :: mesh_new
    type(type_reference_geometry),  intent(in   ) :: refgeo_PD
    type(type_GIA_model),           intent(in   ) :: GIA
    type(type_global_forcing),      intent(in   ) :: forcing
    real(dp),                       intent(in   ) :: time

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'remap_ice_geometry_model'
    real(dp), dimension(:), allocatable :: Hi_old, Hi_new
    real(dp), dimension(:), allocatable :: Hb_old, Hb_new
    real(dp), dimension(:), allocatable :: SL_old, SL_new
    real(dp), dimension(:), allocatable :: Hs_old
    logical,  dimension(:), allocatable :: mask_floating_ice_old, mask_icefree_ocean_old

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap stuff that is common to all models
    call self%remap_model( mesh_new)

    ! Remap stuff that is specific to ice_geometry models

    allocate( Hi_old                ( mesh_old%vi1:mesh_old%vi2), source = self%Hi                ( mesh_old%vi1:mesh_old%vi2))
    allocate( Hb_old                ( mesh_old%vi1:mesh_old%vi2), source = self%Hb                ( mesh_old%vi1:mesh_old%vi2))
    allocate( SL_old                ( mesh_old%vi1:mesh_old%vi2), source = self%SL                ( mesh_old%vi1:mesh_old%vi2))
    allocate( Hs_old                ( mesh_old%vi1:mesh_old%vi2), source = self%Hs                ( mesh_old%vi1:mesh_old%vi2))
    allocate( mask_floating_ice_old ( mesh_old%vi1:mesh_old%vi2), source = self%mask_floating_ice ( mesh_old%vi1:mesh_old%vi2))
    allocate( mask_icefree_ocean_old( mesh_old%vi1:mesh_old%vi2), source = self%mask_icefree_ocean( mesh_old%vi1:mesh_old%vi2))

    allocate( Hi_new                ( mesh_new%vi1:mesh_new%vi2), source = NaN)
    allocate( Hb_new                ( mesh_new%vi1:mesh_new%vi2), source = NaN)
    allocate( SL_new                ( mesh_new%vi1:mesh_new%vi2), source = NaN)

    call remap_ice_geometry_model_bedrock      ( mesh_old, mesh_new, refgeo_PD, GIA, Hb_old, Hb_new)
    call remap_ice_geometry_model_sealevel     ( mesh_old, mesh_new, forcing, time, SL_old, SL_new)
    call remap_ice_geometry_model_ice_thickness( mesh_old, mesh_new, &
      Hi_old, Hs_old, mask_floating_ice_old, mask_icefree_ocean_old, Hb_new, SL_new, Hi_new)

    call reallocate_bounds( self%Hi, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%Hb, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%SL, mesh_new%vi1, mesh_new%vi2)

    self%Hi( mesh_new%vi1:mesh_new%vi2) = Hi_new( mesh_new%vi1:mesh_new%vi2)
    self%Hb( mesh_new%vi1:mesh_new%vi2) = Hb_new( mesh_new%vi1:mesh_new%vi2)
    self%SL( mesh_new%vi1:mesh_new%vi2) = SL_new( mesh_new%vi1:mesh_new%vi2)

    call reallocate_and_recalculate_secondary_geometry_variables( self, mesh_new, refgeo_PD)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_ice_geometry_model

  subroutine remap_ice_geometry_model_bedrock( mesh_old, mesh_new, refgeo_PD, GIA, Hb_old, Hb_new)

    ! In/output variables:
    type(type_mesh),                                intent(in   ) :: mesh_old
    type(type_mesh),                                intent(in   ) :: mesh_new
    type(type_reference_geometry),                  intent(in   ) :: refgeo_PD
    type(type_GIA_model),                           intent(in   ) :: GIA
    real(dp), dimension(mesh_old%vi1:mesh_old%vi2), intent(in   ) :: Hb_old
    real(dp), dimension(mesh_new%vi1:mesh_new%vi2), intent(  out) :: Hb_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_ice_geometry_model_bedrock'

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap bedrock from the original high-resolution grid, and add the (very smooth) modelled deformation to it
    ! Remapping of Hb in the refgeo structure has already happened, only need to copy the data
    if (par%primary) call warning('GIA model isnt finished yet - need to include dHb in mesh update!')
    Hb_new( mesh_new%vi1 : mesh_new%vi2) = refgeo_PD%Hb

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_ice_geometry_model_bedrock

  subroutine remap_ice_geometry_model_sealevel( mesh_old, mesh_new, forcing, time, SL_old, SL_new)

    ! In/output variables:
    type(type_mesh),                                intent(in   ) :: mesh_old
    type(type_mesh),                                intent(in   ) :: mesh_new
    type(type_global_forcing),                      intent(in   ) :: forcing
    real(dp),                                       intent(in   ) :: time
    real(dp), dimension(mesh_old%vi1:mesh_old%vi2), intent(in   ) :: SL_old
    real(dp), dimension(mesh_new%vi1:mesh_new%vi2), intent(  out) :: SL_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_ice_geometry_model_sealevel'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_sealevel_model)
    case default
      call crash('unknown choice_sealevel_model "' // trim( C%choice_sealevel_model) // '"')
    case ('fixed')
      SL_new( mesh_new%vi1 : mesh_new%vi2) = C%fixed_sealevel
    case ('prescribed')
      call update_sealevel_in_model( forcing, mesh_new, time, SL_new)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_ice_geometry_model_sealevel

  subroutine remap_ice_geometry_model_ice_thickness( mesh_old, mesh_new, &
    Hi_old, Hs_old, mask_floating_ice_old, mask_icefree_ocean_old, Hb_new, SL_new, Hi_new)

    ! In/output variables:
    type(type_mesh),                                intent(in   ) :: mesh_old
    type(type_mesh),                                intent(in   ) :: mesh_new
    real(dp), dimension(mesh_old%vi1:mesh_old%vi2), intent(in   ) :: Hi_old
    real(dp), dimension(mesh_old%vi1:mesh_old%vi2), intent(in   ) :: Hs_old
    logical,  dimension(mesh_old%vi1:mesh_old%vi2), intent(in   ) :: mask_floating_ice_old
    logical,  dimension(mesh_old%vi1:mesh_old%vi2), intent(in   ) :: mask_icefree_ocean_old
    real(dp), dimension(mesh_new%vi1:mesh_new%vi2), intent(in   ) :: Hb_new
    real(dp), dimension(mesh_new%vi1:mesh_new%vi2), intent(in   ) :: SL_new
    real(dp), dimension(mesh_new%vi1:mesh_new%vi2), intent(  out) :: Hi_new

    ! Local variables:
    character(len=*), parameter                     :: routine_name = 'remap_ice_geometry_model_ice_thickness'
    real(dp), dimension( mesh_new%vi1:mesh_new%vi2) :: Hs_new
    logical,  dimension( mesh_new%vi1:mesh_new%vi2) :: mask_noice_new
    integer                                         :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! First, naively remap ice thickness and surface elevation without any restrictions
    call map_from_mesh_to_mesh_2D( mesh_old, mesh_new, C%output_dir, Hi_old, Hi_new, '2nd_order_conservative')
    call map_from_mesh_to_mesh_2D( mesh_old, mesh_new, C%output_dir, Hs_old, Hs_new, '2nd_order_conservative')

    ! In grid cells where the remapping predicts there should be ice,
    ! recalculate the remapped ice thickness as the difference between
    ! the new bedrock elevationand the remapped surface elevation

    do vi = mesh_new%vi1, mesh_new%vi2
      if (Hi_new( vi) > 0._dp) then
        if (Hs_new( vi) <= Hb_new( vi)) then
          Hi_new( vi) = 0._dp
        else
          Hi_new( vi) = Hi_from_Hb_Hs_and_SL( Hb_new( vi), Hs_new( vi), SL_new( vi))
        end if
      else
        Hi_new( vi) = 0._dp
      end if
    end do

    ! Apply boundary conditions at the domain border
    call calc_mask_noice( mesh_new, mask_noice_new)
    call apply_ice_thickness_BC_explicit( mesh_new, mask_noice_new, Hb_new, SL_new, Hi_new)

    ! Apply corrections at the ice margin
    call correct_remapped_ice_margin( mesh_old, mesh_new, Hi_old, Hs_old, &
      mask_floating_ice_old, mask_icefree_ocean_old, mask_noice_new, Hi_new)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_ice_geometry_model_ice_thickness

  subroutine correct_remapped_ice_margin( mesh_old, mesh_new, &
    Hi_old, Hs_old, mask_floating_ice_old, mask_icefree_ocean_old, mask_noice_new, Hi_new)

    ! In/output variables:
    type(type_mesh),                                intent(in   ) :: mesh_old
    type(type_mesh),                                intent(in   ) :: mesh_new
    real(dp), dimension(mesh_old%vi1:mesh_old%vi2), intent(in   ) :: Hi_old
    real(dp), dimension(mesh_old%vi1:mesh_old%vi2), intent(in   ) :: Hs_old
    logical,  dimension(mesh_old%vi1:mesh_old%vi2), intent(in   ) :: mask_floating_ice_old
    logical,  dimension(mesh_old%vi1:mesh_old%vi2), intent(in   ) :: mask_icefree_ocean_old
    logical,  dimension(mesh_new%vi1:mesh_new%vi2), intent(in   ) :: mask_noice_new
    real(dp), dimension(mesh_new%vi1:mesh_new%vi2), intent(inout) :: Hi_new

    ! Local variables:
    character(len=*), parameter                     :: routine_name = 'remap_ice_geometry_model_ice_thickness'
    real(dp), dimension( mesh_old%nV)               :: Hi_old_tot
    logical,  dimension( mesh_old%nV)               :: mask_floating_ice_old_tot
    logical,  dimension( mesh_old%nV)               :: mask_icefree_ocean_old_tot
    real(dp), dimension( mesh_new%vi1:mesh_new%vi2) :: Hs_new
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

    call gather_to_all( Hi_old                , Hi_old_tot                )
    call gather_to_all( mask_floating_ice_old , mask_floating_ice_old_tot )
    call gather_to_all( mask_icefree_ocean_old, mask_icefree_ocean_old_tot)

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

      end do

      if (n_ice > 0 .and. n_nonice > 0) then
        ! This new-mesh vertex overlaps with both old-mesh ice vertices,
        ! and old-mesh non-ice vertices
        if (Hi_new( vi_new) < 1._dp) then
          ! Remove very thin remapped ice
          Hi_new( vi_new) = 0._dp
        end if
      end if

    end do

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

        if     (mask_floating_ice_old_tot(  vi_old)) then
          n_shelf      = n_shelf      + 1
          sum_Hi_shelf = sum_Hi_shelf + Hi_old_tot( vi_old)
        elseif (mask_icefree_ocean_old_tot( vi_old)) then
          n_open_ocean = n_open_ocean + 1
        end if

      end do

      if (n_shelf > 0 .and. n_open_ocean > 0) then
        ! This new-mesh vertex overlaps with both old-mesh shelf vertices,
        ! and old-mesh open-ocean vertices
        Hi_new( vi_new) = sum_Hi_shelf / real( n_shelf,dp)
      end if

    end do

    call apply_mask_noice_direct( mesh_new, mask_noice_new, Hi_new)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine correct_remapped_ice_margin

  subroutine reallocate_and_recalculate_secondary_geometry_variables( self, mesh_new, refgeo_PD)

    ! In/output variables:
    class(type_ice_geometry_model), intent(inout) :: self
    type(type_mesh),                intent(in   ) :: mesh_new
    type(type_reference_geometry),  intent(in   ) :: refgeo_PD

    ! Local variables:
    character(len=*), parameter                    :: routine_name = 'reallocate_and_recalculate_secondary_geometry_variables'
    real(dp), dimension(mesh_new%vi1:mesh_new%vi2) :: dHb_new

    ! Add routine to path
    call init_routine( routine_name)

    ! Secondary ice geometry fields
    call reallocate_bounds( self%Hs                     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%Hib                    , mesh_new%vi1, mesh_new%vi2)
    call self%remap_field( mesh_new, 'TAF'                    , self%TAF                    )
    call self%remap_field( mesh_new, 'Hi_eff'                 , self%Hi_eff                 )
    call self%remap_field( mesh_new, 'Hs_slope'               , self%Hs_slope               )
    call self%remap_field( mesh_new, 'Ho'                     , self%Ho                     )

    ! Horizontal derivatives
    call self%remap_field( mesh_new, 'dHib_dx_b', self%dHib_dx_b)
    call self%remap_field( mesh_new, 'dHib_dy_b', self%dHib_dy_b)

    ! Sub-grid bedrock cumulative density functions (CDFs)
    call self%remap_field( mesh_new, 'bedrock_cdf'  , self%bedrock_cdf)
    call self%remap_field( mesh_new, 'bedrock_cdf_b', self%bedrock_cdf_b)

    ! Area fractions
    call self%remap_field( mesh_new, 'fraction_gr'    , self%fraction_gr    )
    call self%remap_field( mesh_new, 'fraction_gr_b'  , self%fraction_gr_b  )
    call self%remap_field( mesh_new, 'fraction_margin', self%fraction_margin)

    ! Ice masks
    call self%remap_field( mesh_new, 'mask_icefree_land' , self%mask_icefree_land )
    call self%remap_field( mesh_new, 'mask_icefree_ocean', self%mask_icefree_ocean)
    call self%remap_field( mesh_new, 'mask_grounded_ice' , self%mask_grounded_ice )
    call self%remap_field( mesh_new, 'mask_floating_ice' , self%mask_floating_ice )
    call self%remap_field( mesh_new, 'mask_margin'       , self%mask_margin       )
    call self%remap_field( mesh_new, 'mask_gl_gr'        , self%mask_gl_gr        )
    call self%remap_field( mesh_new, 'mask_gl_fl'        , self%mask_gl_fl        )
    call self%remap_field( mesh_new, 'mask_cf_gr'        , self%mask_cf_gr        )
    call self%remap_field( mesh_new, 'mask_cf_fl'        , self%mask_cf_fl        )
    call self%remap_field( mesh_new, 'mask_coastline'    , self%mask_coastline    )
    call self%remap_field( mesh_new, 'mask'              , self%mask              )

    ! Only recalculate bedrock CDFs if they are really needed
    ! (as this is a rather time-consuming step)
    if (C%choice_subgrid_grounded_fraction == 'bedrock_CDF' .or. &
        C%choice_subgrid_grounded_fraction == 'bilin_interp_TAF+bedrock_CDF') then
      call self%calc_bedrock_CDFs( mesh_new, refgeo_PD)
    end if

    ! Recalculate all other secondary geometry variables
    ! FIXME: this should receive the already-remapped dHb field!
    dHb_new = 0._dp
    call self%calc_all_secondary_geometry_variables( dHb_new)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine reallocate_and_recalculate_secondary_geometry_variables

end submodule ice_geometry_model_basic_remap