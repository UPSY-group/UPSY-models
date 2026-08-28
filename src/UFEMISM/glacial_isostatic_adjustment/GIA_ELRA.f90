module GIA_ELRA

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: crash, init_routine, finalise_routine
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use GIA_model_types, only: type_GIA_model, type_ELRA_model
  use region_types, only: type_model_region
  use grid_basic, only: setup_square_grid
  use reallocate_mod, only: reallocate_bounds
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary, distribute_gridded_data_from_primary
  use grid_types, only: type_grid
  use reference_geometry_types, only: type_reference_geometry
  use ice_geometry_basics, only: is_floating
  use remapping_main, only: map_from_mesh_vertices_to_xy_grid_2D, map_from_xy_grid_to_mesh_2D
  use kelvin_function, only: kelvin

  implicit none

  private

  public :: run_ELRA_model
  public :: calculate_ELRA_bedrock_deformation_rate
  public :: initialise_ELRA_model
  public :: remap_ELRA_model

contains

  subroutine run_ELRA_model( region)
    ! Use the ELRA model to update bedrock elevation. Once every (dt_bedrock_ELRA) years,
    ! update deformation rates. In all other time steps, just incrementally add deformation.

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_ELRA_model'
    integer                     :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the bedrock deformation rate
    call calculate_ELRA_bedrock_deformation_rate( region%mesh, region%GIA%grid, region%ice, region%ice%geom, region%GIA, region%ELRA)

    ! Update bedrock with last calculated deformation rate
    do vi = region%mesh%vi1, region%mesh%vi2
      region%ice%dHb_dt( vi) = (region%GIA%dHb_next( vi) - region%GIA%dHb_prev( vi)) / C%dt_GIA
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_ELRA_model

  subroutine calculate_ELRA_bedrock_deformation_rate( mesh, grid, ice, geom, GIA, ELRA)
    ! Use the ELRA model to update bedrock deformation rates.

    ! In/output variables:
    type(type_mesh),                      intent(in   ) :: mesh
    type(type_grid),                      intent(in   ) :: grid
    class(atype_ice_model_data),          intent(inout) :: ice
    class(atype_ice_geometry_model_data), intent(in   ) :: geom
    type(type_GIA_model),                 intent(inout) :: GIA
    type(type_ELRA_model),                intent(inout) :: ELRA

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calculate_ELRA_bedrock_deformation_rate'
    integer                     :: vi,i,j,n,k,l,ii,jj
    real(dp)                    :: Lr

    ! Add routine to path
    call init_routine( routine_name)

    ! Influence radius of the lithospheric rigidity
    Lr = (C%ELRA_lithosphere_flex_rigidity / (C%ELRA_mantle_density * grav))**0.25_dp

    ! Calculate the absolute and relative surface loads on the mesh

    do vi = mesh%vi1, mesh%vi2

      ! Absolute surface load
      if (is_floating( geom%Hi( vi), geom%Hb( vi), geom%SL( vi))) then
        ELRA%surface_load_mesh( vi) = (geom%SL( vi) - geom%Hb( vi)) * grid%dx**2 * seawater_density
      elseif (geom%Hi( vi) > 0._dp) then
        ELRA%surface_load_mesh( vi) =  geom%Hi( vi) * grid%dx**2 * ice_density
      else
        ELRA%surface_load_mesh( vi) = 0._dp
      end if

      ! Relative surface load
      GIA%relative_surface_load_mesh( vi) = ELRA%surface_load_mesh( vi) - ELRA%surface_load_GIAeq( vi)

    end do

    ! Map relative surface load to the GIA grid
    call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, GIA%relative_surface_load_mesh, GIA%relative_surface_load_grid)

    !! Gather data to primary
    call gather_gridded_data_to_primary( grid, GIA%relative_surface_load_grid, ELRA%relative_surface_load_grid_tot)

    n = ELRA%flex_prof_rad

    ! Let the primary do the actual work
    if (par%primary) then

      do i = 1, grid%nx
      do j = 1, grid%ny
        ELRA%dHb_eq_grid( i, j) = 0._dp
        do k = -n, n
        do l = -n, n
          ii = max( 1, min( grid%nx, i+k ))
          jj = max( 1, min( grid%ny, j+l ))
          ELRA%dHb_eq_grid( i, j) = ELRA%dHb_eq_grid( i, j) + &
            (0.5_dp * grav * Lr**2 /(pi * C%ELRA_lithosphere_flex_rigidity) * ELRA%relative_surface_load_grid_tot( ii, jj) * ELRA%flex_prof_grid( k+n+1,l+n+1))
        end do
        end do
      end do
      end do

    end if

    ! Map the actual bedrock deformation from the mesh to the grid
    call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, ice%dHb, ELRA%dHb_grid_partial)

    ! gather data from all processors to primary, from partial grid vec to total 2D grid
    call gather_gridded_data_to_primary( grid, ELRA%dHb_grid_partial, ELRA%dHb_grid_tot)

	  ! Let the primary do the actual work
    if (par%primary) then

      ! Calculate the bedrock deformation rate from the difference between the current and the equilibrium deformation
      do i = 1, grid%nx
      do j = 1, grid%ny
        ELRA%dHb_dt_grid( i,j) = (ELRA%dHb_eq_grid( i,j) - ELRA%dHb_grid_tot( i,j)) / C%ELRA_bedrock_relaxation_time
      end do
      end do

    end if

    ! distribute from 2D grid data on primary to vector grid data on all processors
    call distribute_gridded_data_from_primary( grid, ELRA%dHb_dt_grid_partial, ELRA%dHb_dt_grid)

    ! remap from partial grid vec data to mesh model
    call map_from_xy_grid_to_mesh_2D( grid, mesh, C%output_dir, ELRA%dHb_dt_grid_partial, ELRA%dHb_dt_mesh)

    ! multiply the GIA time-step to calculate the bedrock deformation
    GIA%dHb_next = GIA%dHb_prev + ELRA%dHb_dt_mesh * C%dt_GIA

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calculate_ELRA_bedrock_deformation_rate

  subroutine initialise_ELRA_model( mesh, grid, ELRA, refgeo_GIAeq)
    ! Allocate and initialise the ELRA GIA model

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_grid),               intent(in   ) :: grid
    type(type_ELRA_model),         intent(inout) :: ELRA
    type(type_reference_geometry), intent(in   ) :: refgeo_GIAeq

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_ELRA_model'
    integer                     :: i,j,n,k,l
    real(dp)                    :: Lr, r

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) WRITE (0,*) '    Initialising ELRA GIA model...'

    ! Allocate memory
    allocate( ELRA%surface_load_GIAeq( mesh%vi1:mesh%vi2))
    allocate( ELRA%relative_surface_load_grid_tot( grid%nx, grid%ny))
    allocate( ELRA%surface_load_mesh( mesh%vi1:mesh%vi2))
    allocate( ELRA%dHb_eq_grid( grid%nx, grid%ny))
    allocate( ELRA%dHb_grid_partial( grid%n1:grid%n2))
    allocate( ELRA%dHb_grid_tot( grid%nx, grid%ny))
    allocate( ELRA%dHb_dt_grid( grid%nx, grid%ny))
    allocate( ELRA%dHb_dt_grid_partial( grid%n1:grid%n2))
    allocate( ELRA%dHb_dt_mesh( mesh%vi1:mesh%vi2))

    ! Fill in the 2D flexural profile (= Kelvin function), with which
    ! a surface load is convoluted to find the surface deformation
    ! ============================================================

    ! Influence radius of the lithospheric rigidity
    Lr = (C%ELRA_lithosphere_flex_rigidity / (C%ELRA_mantle_density * grav))**0.25_dp

    ! Calculate radius (in number of grid cells) of the flexural profile

    if (par%primary) then
      ELRA%flex_prof_rad = min( ceiling( grid%dx/2._dp), max( 1, int( 6._dp * Lr / grid%dx) - 1))
      n = 2 * ELRA%flex_prof_rad + 1
      allocate( ELRA%flex_prof_grid( n, n))

      ! Calculate flexural profile
      do i = -ELRA%flex_prof_rad, ELRA%flex_prof_rad
      do j = -ELRA%flex_prof_rad, ELRA%flex_prof_rad
        l = i+ELRA%flex_prof_rad+1
        k = j+ELRA%flex_prof_rad+1
        r = grid%dx * sqrt( (real(i,dp))**2 + (real(j,dp))**2)
        ELRA%flex_prof_grid( l,k) = kelvin(r / Lr)
      end do
      end do
    end if

    ! Calculate the reference load
    ! ===============================

    call initialise_ELRA_reference_load( mesh, grid, ELRA, refgeo_GIAeq)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ELRA_model

  subroutine initialise_ELRA_reference_load( mesh, grid, ELRA, refgeo_GIAeq)

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    type(type_grid),               intent(in   ) :: grid
    type(type_ELRA_model),         intent(inout) :: ELRA
    type(type_reference_geometry), intent(in   ) :: refgeo_GIAeq

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_ELRA_reference_load'
    integer                     :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate PD reference load on the mesh
    do vi = mesh%vi1, mesh%vi2
      if (is_floating( refgeo_GIAeq%Hi( vi), refgeo_GIAeq%Hb( vi), 0._dp)) then
        ELRA%surface_load_GIAeq( vi) = -refgeo_GIAeq%Hb( vi) * grid%dx**2 * seawater_density
      elseif (refgeo_GIAeq%Hi( vi) > 0._dp) then
        ELRA%surface_load_GIAeq( vi) =  refgeo_GIAeq%Hi( vi) * grid%dx**2 * ice_density
      else
        ELRA%surface_load_GIAeq( vi) = 0._dp
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ELRA_reference_load

  subroutine remap_ELRA_model( mesh_old, mesh_new, ELRA, refgeo_GIAeq, grid)
    ! Remap or reallocate all the data fields

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh_old
    type(type_mesh),                     intent(in   ) :: mesh_new
    type(type_ELRA_model),               intent(inout) :: ELRA
    type(type_reference_geometry),       intent(in   ) :: refgeo_GIAeq
    type(type_grid),                     intent(in   ) :: grid

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_ELRA_model'
    integer                     :: int_dummy

    ! Add routine to path
    call init_routine( routine_name)
    call reallocate_bounds( ELRA%surface_load_GIAeq, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( ELRA%surface_load_mesh, mesh_new%vi1, mesh_new%vi2)

    ! Recalculate the reference load on the GIA grid
    call initialise_ELRA_reference_load( mesh_new, grid, ELRA, refgeo_GIAeq)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_ELRA_model

end module GIA_ELRA
