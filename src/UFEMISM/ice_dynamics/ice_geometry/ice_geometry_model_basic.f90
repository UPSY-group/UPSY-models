module ice_geometry_model_basic

  use UPSY_main, only: UPSY
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use Arakawa_grid_mod, only: Arakawa_grid
  use mesh_types, only: type_mesh
  use parameters, only: NaN, ice_density, seawater_density
  use checksum_mod, only: checksum
  use model_configuration, only: C
  use mpi_distributed_memory, only: gather_to_all, distribute_from_primary
  use mpi_distributed_shared_memory, only: distribute_dist_shared_from_primary
  use ice_geometry_basics, only: is_floating, thickness_above_floatation, &
    ice_surface_elevation, height_of_water_column_at_ice_front, hi_from_hb_hs_and_sl
  use crash_mod, only: crash, warning
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D, &
    ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D
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
  use petsc_basic, only: mat_petsc2csr, mat_petsc2CSR
  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_DOUBLE_PRECISION
  use remapping_grid_to_mesh_vertices, only: create_map_from_xy_grid_to_mesh_vertices
  use remapping_grid_to_mesh_triangles, only: create_map_from_xy_grid_to_mesh_triangles
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary
  use netcdf_io_main
  use conservation_of_mass_utilities, only: apply_mask_noice_direct
  use plane_geometry, only: triangle_area
  use fields_dimensions, only: third_dimension

  implicit none

  private

  public :: type_ice_geometry_model

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
      procedure, public :: calc_absolute_surface_slope
      procedure, public :: calc_ice_base_slopes

      procedure, public  :: initialise_bedrock_CDFs
      procedure, private :: initialise_bedrock_CDFs_from_file
      procedure, public  :: calc_bedrock_CDFs
      procedure, private :: calc_bedrock_CDFs_a
      procedure, private :: calc_bedrock_CDFs_b

      procedure, public  :: calc_grounded_fractions
      procedure, private :: calc_grounded_fractions_bilin_interp_TAF_a
      procedure, private :: calc_grounded_fractions_bilin_interp_TAF_b
      procedure, private :: calc_grounded_fractions_bedrock_CDF_a
      procedure, private :: calc_grounded_fractions_bedrock_CDF_b

      procedure, public :: calc_all_secondary_geometry_variables

      procedure, public :: get_model_name

  end type type_ice_geometry_model

  ! Interfaces for procedures defined in submodules
  interface

    module subroutine remap_ice_geometry_model( self, mesh_old, mesh_new, refgeo_PD, GIA, forcing, time)
      class(type_ice_geometry_model),       intent(inout) :: self
      type(type_mesh),                      intent(in   ) :: mesh_old
      type(type_mesh),                      intent(in   ) :: mesh_new
      type(type_reference_geometry),        intent(in   ) :: refgeo_PD
      type(type_GIA_model),                 intent(in   ) :: GIA
      type(type_global_forcing),            intent(in   ) :: forcing
      real(dp),                             intent(in   ) :: time
    end subroutine remap_ice_geometry_model

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

    module subroutine calc_absolute_surface_slope( self)
      class(type_ice_geometry_model), intent(inout) :: self
    end subroutine calc_absolute_surface_slope

    module subroutine calc_ice_base_slopes( self)
      class(type_ice_geometry_model), intent(inout) :: self
    end subroutine calc_ice_base_slopes

    module subroutine initialise_bedrock_CDFs( self, mesh, refgeo, region_name)
      class(type_ice_geometry_model), intent(inout) :: self
      type(type_mesh),                intent(in   ) :: mesh
      type(type_reference_geometry),  intent(in   ) :: refgeo
      character(len=3),               intent(in   ) :: region_name
    end subroutine initialise_bedrock_CDFs

    module subroutine initialise_bedrock_CDFs_from_file( self, mesh, region_name)
      class(type_ice_geometry_model), intent(inout) :: self
      type(type_mesh),                intent(in   ) :: mesh
      character(len=3),               intent(in   ) :: region_name
    end subroutine initialise_bedrock_CDFs_from_file

    module subroutine calc_bedrock_CDFs( self, mesh, refgeo)
      class(type_ice_geometry_model), intent(inout) :: self
      type(type_mesh),                intent(in   ) :: mesh
      type(type_reference_geometry),  intent(in   ) :: refgeo
    end subroutine calc_bedrock_CDFs

    module subroutine calc_bedrock_CDFs_a( self, mesh, refgeo)
      class(type_ice_geometry_model), intent(inout) :: self
      type(type_mesh),                intent(in   ) :: mesh
      type(type_reference_geometry),  intent(in   ) :: refgeo
    end subroutine calc_bedrock_CDFs_a

    module subroutine calc_bedrock_CDFs_b( self, mesh, refgeo)
      class(type_ice_geometry_model), intent(inout) :: self
      type(type_mesh),                intent(in   ) :: mesh
      type(type_reference_geometry),  intent(in   ) :: refgeo
    end subroutine calc_bedrock_CDFs_b

    module subroutine calc_grounded_fractions( self, dHb)
      class(type_ice_geometry_model),                   intent(inout) :: self
      real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: dHb
    end subroutine calc_grounded_fractions

    module subroutine calc_grounded_fractions_bilin_interp_TAF_a( self, fraction_gr)
      class(type_ice_geometry_model),                   intent(in   ) :: self
      real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: fraction_gr
    end subroutine calc_grounded_fractions_bilin_interp_TAF_a

    module subroutine calc_grounded_fractions_bilin_interp_TAF_b( self, fraction_gr_b)
      class(type_ice_geometry_model),                   intent(in   ) :: self
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(  out) :: fraction_gr_b
    end subroutine calc_grounded_fractions_bilin_interp_TAF_b

    module subroutine calc_grounded_fractions_bedrock_CDF_a( self, dHb, fraction_gr)
      class(type_ice_geometry_model),                   intent(in   ) :: self
      real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: dHb
      real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(  out) :: fraction_gr
    end subroutine calc_grounded_fractions_bedrock_CDF_a

    module subroutine calc_grounded_fractions_bedrock_CDF_b( self, dHb, fraction_gr_b)
      class(type_ice_geometry_model),                   intent(in   ) :: self
      real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: dHb
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2), intent(  out) :: fraction_gr_b
    end subroutine calc_grounded_fractions_bedrock_CDF_b

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
    call self%create_field( self%dHib_dx_b, self%wdHib_dx_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'dHib_dx_b', &
      long_name = 'Horizontal derivative in x-direction of ice draft on b-grid', &
      units     = '-', &
      remap_method = 'reallocate')
    call self%create_field( self%dHib_dy_b, self%wdHib_dy_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'dHib_dy_b', &
      long_name = 'Horizontal derivative in y-direction of ice draft on b-grid', &
      units     = '-', &
      remap_method = 'reallocate')

    ! Sub-grid bedrock cumulative density functions (CDFs)
    call self%create_field( self%bedrock_cdf, self%wbedrock_cdf, &
      self%mesh, Arakawa_grid%a(), third_dimension%bedrock_CDF( C%subgrid_bedrock_cdf_nbins), &
      name      = 'bedrock_cdf', &
      long_name = 'Bedrock CDF of vertices', &
      units     = 'm', &
      remap_method = 'reallocate')
    call self%create_field( self%bedrock_cdf_b, self%wbedrock_cdf_b, &
      self%mesh, Arakawa_grid%b(), third_dimension%bedrock_CDF( C%subgrid_bedrock_cdf_nbins), &
      name      = 'bedrock_cdf_b', &
      long_name = 'Bedrock CDF of triangles', &
      units     = 'm', &
      remap_method = 'reallocate')

    ! Area fractions
    call self%create_field( self%fraction_gr, self%wfraction_gr, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'fraction_gr', &
      long_name = 'Grounded area fractions of vertices', &
      units     = '0-1', &
      remap_method = 'reallocate')
    call self%create_field( self%fraction_gr_b, self%wfraction_gr_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'fraction_gr_b', &
      long_name = 'Grounded area fractions of triangles', &
      units     = '0-1', &
      remap_method = 'reallocate')
    call self%create_field( self%fraction_margin, self%wfraction_margin, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'fraction_margin', &
      long_name = 'Ice-covered area fractions of ice margins', &
      units     = '0-1', &
      remap_method = 'reallocate')

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
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 7)

  end subroutine allocate_ice_geometry_model

  subroutine deallocate_ice_geometry_model( self)

    ! In/output variables:
    class(type_ice_geometry_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_ice_geometry_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate stuff that is common to all models
    call self%deallocate_model()

    ! Deallocate stuff that is specific to the ice_geometry model

    ! Primary ice geometry variables
    deallocate( self%Hi)
    deallocate( self%Hb)
    deallocate( self%SL)

    ! Derived ice geometry variables
    deallocate( self%Hs      )
    deallocate( self%Hib     )
    deallocate( self%TAF     )
    deallocate( self%Hi_eff  )
    deallocate( self%Hs_slope)
    deallocate( self%Ho      )

    ! Horizontal derivatives
    nullify( self%dHib_dx_b)
    nullify( self%dHib_dy_b)

    ! Sub-grid bedrock cumulative density functions (CDFs)
    nullify( self%bedrock_cdf  )
    nullify( self%bedrock_cdf_b)

    ! Area fractions
    nullify( self%fraction_gr    )
    nullify( self%fraction_gr_b  )
    nullify( self%fraction_margin)

    ! Ice masks
    deallocate( self%mask_icefree_land )
    deallocate( self%mask_icefree_ocean)
    deallocate( self%mask_grounded_ice )
    deallocate( self%mask_floating_ice )
    deallocate( self%mask_margin       )
    deallocate( self%mask_gl_gr        )
    deallocate( self%mask_gl_fl        )
    deallocate( self%mask_cf_gr        )
    deallocate( self%mask_cf_fl        )
    deallocate( self%mask_coastline    )
    deallocate( self%mask              )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_ice_geometry_model

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

end module ice_geometry_model_basic