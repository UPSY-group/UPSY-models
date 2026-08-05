module ice_geometry_model_basic

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use Arakawa_grid_mod, only: Arakawa_grid
  use mesh_types, only: type_mesh
  use parameters, only: NaN
  use checksum_mod, only: checksum
  use model_configuration, only: C
  use mpi_distributed_memory, only: gather_to_all
  use ice_geometry_basics, only: is_floating, thickness_above_floatation, &
    ice_surface_elevation, height_of_water_column_at_ice_front, hi_from_hb_hs_and_sl
  use crash_mod, only: crash, warning
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D, &
    ddx_a_b_2D, ddy_a_b_2D
  use reference_geometry_types, only: type_reference_geometry
  use GIA_model_types, only: type_GIA_model
  use global_forcing_types, only: type_global_forcing
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use mpi_basic, only: par
  use remapping_main, only: Atlas, map_from_mesh_to_mesh_2D
  use reallocate_mod, only: reallocate_bounds
  use global_forcings_main, only: update_sealevel_in_model
  use masks_mod, only: calc_mask_noice
  use ice_thickness_boundary_conditions, only: apply_ice_thickness_BC_explicit
  use petsc_basic, only: mat_petsc2csr

  implicit none

  private

  public :: type_ice_geometry_model, remap_basic_ice_geometry

  type, extends(atype_ice_geometry_model_data) :: type_ice_geometry_model

    contains

      procedure, public :: allocate   => allocate_ice_geometry_model
      procedure, public :: deallocate => deallocate_ice_geometry_model
      procedure, public :: remap      => remap_ice_geometry_model

      final :: finalise_ice_geometry_model

      procedure, public :: calc_surface_elevation
      procedure, public :: calc_ice_base_elevation
      procedure, public :: calc_thickness_above_floatation
      procedure, public :: calc_height_of_water_column
      procedure, public :: determine_masks
      procedure, public :: calc_effective_thickness
      procedure, public :: calc_grounded_fractions
      procedure, public :: calc_absolute_surface_slope
      procedure, public :: calc_ice_base_slopes

      procedure, public :: calc_all_secondary_geometry_variables

      procedure, public :: get_model_name

  end type type_ice_geometry_model

  ! Interfaces for procedures defined in submodules
  interface

    module subroutine calc_surface_elevation( self)
      class(type_ice_geometry_model),intent(inout) :: self
    end subroutine calc_surface_elevation

    module subroutine calc_ice_base_elevation( self)
      class(type_ice_geometry_model),intent(inout) :: self
    end subroutine calc_ice_base_elevation

    module subroutine calc_thickness_above_floatation( self)
      class(type_ice_geometry_model),intent(inout) :: self
    end subroutine calc_thickness_above_floatation

    module subroutine calc_height_of_water_column( self)
      class(type_ice_geometry_model),intent(inout) :: self
    end subroutine calc_height_of_water_column

    module subroutine determine_masks( self)
      class(type_ice_geometry_model),intent(inout) :: self
    end subroutine determine_masks

    module subroutine calc_effective_thickness( self)
      class(type_ice_geometry_model), intent(inout) :: self
    end subroutine calc_effective_thickness

    module subroutine calc_grounded_fractions( self, dHb)
      class(type_ice_geometry_model),                   intent(inout) :: self
      real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: dHb
    end subroutine calc_grounded_fractions

    module subroutine calc_absolute_surface_slope( self)
      class(type_ice_geometry_model), intent(inout) :: self
    end subroutine calc_absolute_surface_slope

    module subroutine calc_ice_base_slopes( self)
      class(type_ice_geometry_model), intent(inout) :: self
    end subroutine calc_ice_base_slopes

  end interface

contains

  subroutine allocate_ice_geometry_model( self, region_name, mesh)

    ! In/output variables:
    class(type_ice_geometry_model), intent(inout) :: self
    character(len=*),               intent(in   ) :: region_name
    type(type_mesh), target,        intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_ice_geometry_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all models
    call self%allocate_model( region_name, mesh)

    ! Allocate all the stuff that is specific to the ice_geometry model

    ! Primary ice geometry variables
    allocate( self%Hi( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%Hb( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%SL( mesh%vi1:mesh%vi2), source = NaN)

    ! Derived ice geometry variables
    allocate( self%Hs      ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%Hib     ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%TAF     ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%Hi_eff  ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%Hs_slope( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%Ho      ( mesh%vi1:mesh%vi2), source = NaN)

    ! Horizontal derivatives
    allocate( self%dHib_dx_b( mesh%ti1:mesh%ti2), source = NaN)
    allocate( self%dHib_dy_b( mesh%ti1:mesh%ti2), source = NaN)

    ! Sub-grid bedrock cumulative density functions (CDFs)
    allocate( self%bedrock_cdf  ( mesh%vi1:mesh%vi2, C%subgrid_bedrock_cdf_nbins), source = NaN)
    allocate( self%bedrock_cdf_b( mesh%ti1:mesh%ti2, C%subgrid_bedrock_cdf_nbins), source = NaN)

    ! Area fractions
    allocate( self%fraction_gr    ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%fraction_gr_b  ( mesh%ti1:mesh%ti2), source = NaN)
    allocate( self%fraction_margin( mesh%vi1:mesh%vi2), source = NaN)

    ! Ice masks
    allocate( self%mask_icefree_land ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( self%mask_icefree_ocean( mesh%vi1:mesh%vi2), source = .false.)
    allocate( self%mask_grounded_ice ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( self%mask_floating_ice ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( self%mask_margin       ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( self%mask_gl_gr        ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( self%mask_gl_fl        ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( self%mask_cf_gr        ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( self%mask_cf_fl        ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( self%mask_coastline    ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( self%mask              ( mesh%vi1:mesh%vi2), source = -42)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_ice_geometry_model

  subroutine deallocate_ice_geometry_model( self)

    ! In/output variables:
    class(type_ice_geometry_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_ice_geometry_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! DENK DROM

    ! Deallocate stuff that is common to all models
    ! call self%deallocate_model()

    ! Deallocate stuff that is specific to the ice_geometry model

    ! nullify( self%Hi)
    ! nullify( self%Hb)
    ! nullify( self%SL)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_ice_geometry_model

  subroutine remap_ice_geometry_model( self, mesh_new)

    ! In/output variables:
    class(type_ice_geometry_model), intent(inout) :: self
    type(type_mesh), target,        intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_ice_geometry_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap stuff that is common to all models
    call self%remap_model( mesh_new)

    ! Remap stuff that is specific to ice_geometry models

    ! call self%remap_field( mesh_new, 'Hi', self%Hi)
    ! call self%remap_field( mesh_new, 'Hb', self%Hb)
    ! call self%remap_field( mesh_new, 'SL', self%SL)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_ice_geometry_model

  subroutine finalise_ice_geometry_model( self)

    ! In/output variables:
    type(type_ice_geometry_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'finalise_ice_geometry_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%deallocate()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine finalise_ice_geometry_model

  function get_model_name( self) result( model_name)
    class(type_ice_geometry_model), intent(in) :: self
    character(len=:), allocatable              :: model_name
    model_name = 'ice_geometry'
  end function get_model_name

  subroutine calc_all_secondary_geometry_variables( self, dHb)

    ! In/output variables:
    class(type_ice_geometry_model), intent(inout) :: self
    real(dp), dimension(:),         intent(in   ) :: dHb

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_all_secondary_geometry_variables'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%calc_surface_elevation()
    call self%calc_ice_base_elevation()
    call self%calc_thickness_above_floatation()
    call self%calc_height_of_water_column()
    call self%determine_masks()
    call self%calc_effective_thickness()
    call self%calc_grounded_fractions( dHb)
    call self%calc_absolute_surface_slope()
    call self%calc_ice_base_slopes()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calc_all_secondary_geometry_variables

  subroutine remap_basic_ice_geometry( mesh_old, mesh_new, refgeo_PD, GIA, geom, mask_noice, forcing, time)
    !< Remap the basic ice geometry Hi,Hb,Hs,SL.

    ! In/output variables:
    type(type_mesh),                      intent(in   ) :: mesh_old
    type(type_mesh),                      intent(in   ) :: mesh_new
    type(type_reference_geometry),        intent(in   ) :: refgeo_PD
    type(type_GIA_model),                 intent(in   ) :: GIA
    class(atype_ice_geometry_model_data), intent(inout) :: geom
    logical, dimension(:), allocatable,   intent(inout) :: mask_noice
    type(type_global_forcing),            intent(in   ) :: forcing
    real(dp),                             intent(in   ) :: time

    ! Local variables:
    character(len=*), parameter                     :: routine_name = 'remap_basic_ice_geometry'
    real(dp), dimension( mesh_old%nV)               :: Hi_old_tot
    logical,  dimension( mesh_old%nV)               :: mask_floating_ice_tot
    logical,  dimension( mesh_old%nV)               :: mask_icefree_ocean_tot
    real(dp), dimension( mesh_new%vi1:mesh_new%vi2) :: Hi_new
    real(dp), dimension( mesh_new%vi1:mesh_new%vi2) :: Hs_new
    integer                                         :: vi
    integer                                         :: mi, mi_used
    logical                                         :: found_map
    type(type_CSR_matrix_dp)                 :: M_CSR
    integer                                         :: vi_new, k1, k2, k, vi_old
    integer                                         :: n_ice, n_nonice
    integer                                         :: n_shelf, n_open_ocean
    real(dp)                                        :: sum_Hi_shelf

    ! Add routine to path
    call init_routine( routine_name)

    ! == Basic: remap surface elevation Hs from the old mesh, remap bedrock elevation Hb
    !    from its (presumably high-resolution) source grid, define remapped ice thickness
    !    as the difference between the two. As surface elevation is typically much smoother
    !    then ice thickness, remapping works much better.
    ! =====================================================================================

    ! Remap bedrock from the original high-resolution grid, and add the (very smooth) modelled deformation to it
    ! Remapping of Hb in the refgeo structure has already happened, only need to copy the data
    if (par%primary) call warning('GIA model isnt finished yet - need to include dHb in mesh update!')
    call reallocate_bounds( geom%Hb, mesh_new%vi1, mesh_new%vi2)  ! [m] Bedrock elevation (w.r.t. PD sea level)
    geom%Hb = refgeo_PD%Hb

    ! Remap sea level
    call reallocate_bounds( geom%SL, mesh_new%vi1, mesh_new%vi2)  ! [m] Sea level (geoid) elevation (w.r.t. PD sea level)
    select case (C%choice_sealevel_model)
    case default
      call crash('unknown choice_sealevel_model "' // trim( C%choice_sealevel_model) // '"')
    case ('fixed')
      geom%SL = C%fixed_sealevel
    case ('prescribed')
      call update_sealevel_in_model(forcing, mesh_new, geom, time)
    end select

    ! Gather global ice thickness and masks
    call gather_to_all( geom%Hi, Hi_old_tot)
    call gather_to_all( geom%mask_floating_ice , mask_floating_ice_tot )
    call gather_to_all( geom%mask_icefree_ocean, mask_icefree_ocean_tot)

    ! First, naively remap ice thickness and surface elevation without any restrictions
    call map_from_mesh_to_mesh_2D( mesh_old, mesh_new, C%output_dir, geom%Hi, Hi_new, '2nd_order_conservative')
    call map_from_mesh_to_mesh_2D( mesh_old, mesh_new, C%output_dir, geom%Hs, Hs_new, '2nd_order_conservative')

    ! Calculate remapped ice thickness as the difference between new bedrock and remapped surface elevation
    do vi = mesh_new%vi1, mesh_new%vi2
      if (Hi_new( vi) > 0._dp) then
        if (Hs_new( vi) <= geom%Hb( vi)) then
          Hi_new( vi) = 0._dp
        else
          Hi_new( vi) = Hi_from_Hb_Hs_and_SL( geom%Hb( vi), Hs_new( vi), geom%SL( vi))
        end if
      else
        Hi_new( vi) = 0._dp
      end if
    end do ! do vi = mesh_new%vi1, mesh_new%vi2

    ! reallocate no-ice mask
    ! T: no ice is allowed here, F: ice is allowed here
    call reallocate_bounds( mask_noice, mesh_new%vi1, mesh_new%vi2)

    ! Apply boundary conditions at the domain border
    call calc_mask_noice( mesh_new, mask_noice)
    call apply_ice_thickness_BC_explicit( mesh_new, mask_noice, geom%Hb, geom%SL, Hi_new)

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

      end do ! do k = k1, k2

      if (n_shelf > 0 .and. n_open_ocean > 0) then
        ! This new-mesh vertex overlaps with both old-mesh shelf vertices,
        ! and old-mesh open-ocean vertices
        Hi_new( vi_new) = sum_Hi_shelf / real( n_shelf,dp)
      end if

    end do ! do vi_new = mesh_new%vi1, mesh_new%vi2

    ! Recalculate Hs
    call reallocate_bounds( geom%Hs, mesh_new%vi1, mesh_new%vi2)
    do vi = mesh_new%vi1, mesh_new%vi2
      geom%Hs( vi) = ice_surface_elevation( Hi_new( vi), geom%Hb( vi), geom%SL( vi))
    end do

    ! Move Hi_new to geom%Hi
    deallocate( geom%Hi)
    allocate( geom%Hi( mesh_new%vi1: mesh_new%vi2))
    geom%Hi = Hi_new

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_basic_ice_geometry

end module ice_geometry_model_basic