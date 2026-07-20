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
  use ice_geometry_basics, only: is_floating

  implicit none

  private

  public :: type_ice_geometry_model

  type, extends(atype_ice_geometry_model_data) :: type_ice_geometry_model

    contains

      procedure, public :: allocate   => allocate_ice_geometry_model
      procedure, public :: deallocate => deallocate_ice_geometry_model
      procedure, public :: remap      => remap_ice_geometry_model

      final :: finalise_ice_geometry_model

      procedure, public :: determine_masks

  end type type_ice_geometry_model

contains

  subroutine allocate_ice_geometry_model( self, name, region_name, mesh)

    ! In/output variables:
    class(type_ice_geometry_model), intent(inout) :: self
    character(len=*),               intent(in   ) :: name
    character(len=*),               intent(in   ) :: region_name
    type(type_mesh), target,        intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_ice_geometry_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all models
    call self%allocate_model( name, region_name, mesh)

    ! Allocate all the stuff that is specific to the ice_geometry model

    ! DENK DROM
    allocate( self%Hi( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%Hb( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%SL( mesh%vi1:mesh%vi2), source = NaN)

    ! call self%create_field( self%Hi, self%wHi, &
    !   self%mesh, Arakawa_grid%a(), &
    !   name      = 'Hi', &
    !   long_name = 'Ice thickness', &
    !   units     = 'm', &
    !   remap_method = 'reallocate')

    ! call self%create_field( self%Hb, self%wHb, &
    !   self%mesh, Arakawa_grid%a(), &
    !   name      = 'Hb', &
    !   long_name = 'Bedrock elevation (w.r.t. PD sea level)', &
    !   units     = 'm', &
    !   remap_method = 'reallocate')

    ! call self%create_field( self%SL, self%wSL, &
    !   self%mesh, Arakawa_grid%a(), &
    !   name      = 'SL', &
    !   long_name = 'Geoid elevation (w.r.t. PD sea level)', &
    !   units     = 'm', &
    !   remap_method = 'reallocate')

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


  subroutine determine_masks( self, &
    mask, mask_icefree_land, mask_icefree_ocean, mask_grounded_ice, mask_floating_ice, &
    mask_margin, mask_gl_fl, mask_gl_gr, mask_cf_gr, mask_cf_fl,mask_coastline)
    !< Determine the different masks

    ! In/output variables:
    class(type_ice_geometry_model),                   intent(inout) :: self
    logical,  dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: mask_icefree_land       ! T: ice-free land , F: otherwise
    logical,  dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: mask_icefree_ocean      ! T: ice-free ocean, F: otherwise
    logical,  dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: mask_grounded_ice       ! T: grounded ice  , F: otherwise
    logical,  dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: mask_floating_ice       ! T: floating ice  , F: otherwise
    logical,  dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: mask_margin             ! T: ice next to ice-free, F: otherwise
    logical,  dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: mask_gl_gr              ! T: grounded ice next to floating ice, F: otherwise
    logical,  dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: mask_gl_fl              ! T: floating ice next to grounded ice, F: otherwise
    logical,  dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: mask_cf_gr              ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    logical,  dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: mask_cf_fl              ! T: floating ice next to ice-free water (sea or lake), F: otherwise
    logical,  dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: mask_coastline          ! T: ice-free land next to ice-free ocean, F: otherwise
    integer,  dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: mask

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
    mask_icefree_land  = .false.
    mask_icefree_ocean = .false.
    mask_grounded_ice  = .false.
    mask_floating_ice  = .false.
    mask               = 0

    ! Calculate
    do vi = self%mesh%vi1, self%mesh%vi2

      if (is_floating( self%Hi( vi), self%Hb( vi), self%SL( vi))) then
        ! Ice thickness is below the floatation thickness; either floating ice, or ice-free ocean

        if (self%Hi( vi) > 0._dp) then
          ! Floating ice

          mask_floating_ice( vi) = .true.
          mask( vi) = C%type_floating_ice

        else
          ! Ice-free ocean

          mask_icefree_ocean( vi) = .true.
          mask( vi) = C%type_icefree_ocean

        end if

      else
        ! Ice thickness is above the floatation thickness; either grounded ice, or ice-free land

        if (self%Hi( vi) > 0._dp) then
          ! Grounded ice

          mask_grounded_ice( vi) = .true.
          mask( vi) = C%type_grounded_ice

        else
          ! Ice-free land

          mask_icefree_land( vi) = .true.
          mask( vi) = C%type_icefree_land

        end if

      end if

    end do

    call checksum( self%mesh%pai_V, mask_icefree_land , 'mask_icefree_land')
    call checksum( self%mesh%pai_V, mask_icefree_ocean, 'mask_icefree_ocean')
    call checksum( self%mesh%pai_V, mask_grounded_ice , 'mask_grounded_ice')
    call checksum( self%mesh%pai_V, mask_floating_ice , 'mask_floating_ice')
    call checksum( self%mesh%pai_V, mask              , 'mask')

    ! === Transitional masks ===
    ! ==========================

    ! Gather basic masks to all processes
    call gather_to_all( mask_icefree_land , mask_icefree_land_tot )
    call gather_to_all( mask_icefree_ocean, mask_icefree_ocean_tot)
    call gather_to_all( mask_grounded_ice , mask_grounded_ice_tot )
    call gather_to_all( mask_floating_ice , mask_floating_ice_tot )

    ! Initialise transitional masks
    mask_margin    = .false.
    mask_gl_gr     = .false.
    mask_gl_fl     = .false.
    mask_cf_gr     = .false.
    mask_cf_fl     = .false.
    mask_coastline = .false.

    do vi = self%mesh%vi1, self%mesh%vi2

      ! Ice margin
      if (mask_grounded_ice_tot( vi) .OR. mask_floating_ice_tot( vi)) then
        do ci = 1, self%mesh%nC( vi)
          vj = self%mesh%C( vi,ci)
          if (.not. (mask_grounded_ice_tot( vj) .OR. mask_floating_ice_tot( vj))) then
            mask_margin( vi) = .true.
            mask( vi) = C%type_margin
          end if
        end do
      end if

      ! Grounding line (grounded side)
      if (mask_grounded_ice_tot( vi)) then
        do ci = 1, self%mesh%nC( vi)
          vj = self%mesh%C( vi,ci)
          if (mask_floating_ice_tot( vj)) then
            mask_gl_gr( vi) = .true.
            mask( vi) = C%type_groundingline_gr
          end if
        end do
      end if

      ! Grounding line (floating side)
      if (mask_floating_ice_tot( vi)) then
        do ci = 1, self%mesh%nC( vi)
          vj = self%mesh%C( vi,ci)
          if (mask_grounded_ice_tot( vj)) then
            mask_gl_fl( vi) = .true.
            mask( vi) = C%type_groundingline_fl
          end if
        end do
      end if

      ! Calving front (grounded)
      if (mask_grounded_ice_tot( vi)) then
        do ci = 1, self%mesh%nC(vi)
          vj = self%mesh%C( vi,ci)
          if (mask_icefree_ocean_tot( vj)) then
            mask_cf_gr( vi) = .true.
            mask( vi) = C%type_calvingfront_gr
          end if
        end do
      end if

      ! Calving front (floating)
      if (mask_floating_ice_tot( vi)) then
        do ci = 1, self%mesh%nC(vi)
          vj = self%mesh%C( vi,ci)
          if (mask_icefree_ocean_tot( vj)) then
            mask_cf_fl( vi) = .true.
            mask( vi) = C%type_calvingfront_fl
          end if
        end do
      end if

      ! Coastline
      if (mask_icefree_land_tot( vi)) then
        do ci = 1, self%mesh%nC(vi)
          vj = self%mesh%C( vi,ci)
          if (mask_icefree_ocean_tot( vj)) then
            mask_coastline( vi) = .true.
            mask( vi) = C%type_coastline
          end if
        end do
      end if

    end do ! do vi = self%mesh%vi1, self%mesh%vi2

    call checksum( self%mesh%pai_V, mask_margin   , 'mask_margin')
    call checksum( self%mesh%pai_V, mask_gl_gr    , 'mask_gl_gr')
    call checksum( self%mesh%pai_V, mask_gl_fl    , 'mask_gl_fl')
    call checksum( self%mesh%pai_V, mask_cf_gr    , 'mask_cf_gr')
    call checksum( self%mesh%pai_V, mask_cf_fl    , 'mask_cf_fl')
    call checksum( self%mesh%pai_V, mask_coastline, 'mask_coastline')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_masks

end module ice_geometry_model_basic