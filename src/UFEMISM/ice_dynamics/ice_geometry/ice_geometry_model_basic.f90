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
  use ice_geometry_basics, only: is_floating, thickness_above_floatation
  use crash_mod, only: crash

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
      procedure, public :: calc_effective_thickness
      procedure, public :: calc_grounded_fractions

  end type type_ice_geometry_model

  ! Interfaces for procedures defined in submodules
  interface

    module subroutine determine_masks( self)
      class(type_ice_geometry_model),intent(inout) :: self
    end subroutine determine_masks

    module subroutine calc_effective_thickness( self, Hi_eff, fraction_margin)
      class(type_ice_geometry_model),                   intent(inout) :: self
      real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: Hi_eff
      real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: fraction_margin
    end subroutine calc_effective_thickness

    module subroutine calc_grounded_fractions( self, dHb)
      class(type_ice_geometry_model),                   intent(inout) :: self
      real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: dHb
    end subroutine calc_grounded_fractions

  end interface

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

    allocate( self%Hi( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%Hb( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%SL( mesh%vi1:mesh%vi2), source = NaN)

    allocate( self%Hs ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%Hib( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%TAF( mesh%vi1:mesh%vi2), source = NaN)

    ! Sub-grid bedrock cumulative density functions (CDFs)
    allocate( self%bedrock_cdf  ( mesh%vi1:mesh%vi2, C%subgrid_bedrock_cdf_nbins), source = NaN)
    allocate( self%bedrock_cdf_b( mesh%ti1:mesh%ti2, C%subgrid_bedrock_cdf_nbins), source = NaN)

    ! Area fractions
    allocate( self%fraction_gr  ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%fraction_gr_b( mesh%ti1:mesh%ti2), source = NaN)

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

end module ice_geometry_model_basic