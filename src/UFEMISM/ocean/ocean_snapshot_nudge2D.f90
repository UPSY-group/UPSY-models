module ocean_snapshot_nudge2D

  use precisions, only: dp
  use parameters, only: NaN
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning, colour_string
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use grid_types, only: type_grid
  use ice_model_types, only: type_ice_model
  use ocean_model_types, only: type_ocean_model, type_ocean_model_snapshot_nudge2D
  use netcdf_io_main
  use remapping_main
  use reference_geometry_types, only: type_reference_geometry
  use reference_geometries_main, only: reallocate_reference_geometry_on_mesh
  use ice_geometry_basics, only: is_floating
  use mpi_distributed_memory, only: gather_to_all
  use mesh_utilities, only: extrapolate_Gaussian
  use mesh_data_smoothing, only: smooth_Gaussian

  implicit none

  private

  public :: initialise_ocean_model_snapshot_nudge2D, run_ocean_model_snapshot_nudge2D

contains

  subroutine run_ocean_model_snapshot_nudge2D( mesh, grid_smooth, ice, ocean, time)

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    type(type_grid),        intent(in   ) :: grid_smooth
    type(type_ice_model),   intent(in   ) :: ice
    type(type_ocean_model), intent(inout) :: ocean
    real(dp),               intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_ocean_model_snapshot_nudge2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Only nudge during the user-defined time window
    if (time < C%BMB_inversion_t_start .or.  time > C%BMB_inversion_t_end) return

    call nudge_deltaT( mesh, grid_smooth, ice, ocean%snapshot_nudge2D)
    call map_deltaT_to_reference_grid( mesh, ocean%snapshot_nudge2D)
    call add_deltaT_to_snapshot( mesh, ocean%snapshot_nudge2D)
    call set_applied_ocean_data( mesh, ocean)
    call write_to_output_file( mesh, ocean%snapshot_nudge2D)

    ! call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_ocean_model_snapshot_nudge2D

  subroutine nudge_deltaT( mesh, grid_smooth, ice, snapshot_nudge2D)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_grid),                         intent(in   ) :: grid_smooth
    type(type_ice_model),                    intent(in   ) :: ice
    type(type_ocean_model_snapshot_nudge2D), intent(inout) :: snapshot_nudge2D

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'nudge_deltaT'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: Hi_target_corr
    real(dp), dimension(mesh%vi1:mesh%vi2) :: dTdt
    integer                                :: vi
    real(dp), parameter                    :: deltaT_max =  2._dp
    real(dp), parameter                    :: deltaT_min = -2._dp

    ! Add routine to call stack
    call init_routine( routine_name)

    call calc_corrected_target_thickness( mesh, ice, snapshot_nudge2D, Hi_target_corr)
    call calc_dTdt( mesh, grid_smooth, ice, Hi_target_corr, snapshot_nudge2D%target_mask_shelf, dTdt)

    do vi = mesh%vi1, mesh%vi2

      snapshot_nudge2D%deltaT_nudge( vi) = snapshot_nudge2D%deltaT_nudge( vi) + C%dt_ocean * dTdt( vi)

      ! Apply limits
      snapshot_nudge2D%deltaT_nudge( vi) = max( deltaT_min, min( deltaT_max, &
        snapshot_nudge2D%deltaT_nudge( vi)))

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine nudge_deltaT

  subroutine calc_corrected_target_thickness( mesh, ice, snapshot_nudge2D, Hi_target_corr)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ice_model),                    intent(in   ) :: ice
    type(type_ocean_model_snapshot_nudge2D), intent(in   ) :: snapshot_nudge2D
    real(dp), dimension(mesh%vi1:mesh%vi2),  intent(  out) :: Hi_target_corr

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_corrected_target_thickness'
    integer                        :: vi, ci, vj
    logical,  dimension(mesh%nV)   :: mask_floating_ice_tot, mask_cf_fl_tot
    real(dp), dimension(mesh%nV)   :: Hi_target_tot
    real(dp)                       :: w_sum, wH_sum

    ! Add routine to call stack
    call init_routine( routine_name)

    Hi_target_corr = snapshot_nudge2D%target_geometry%Hi

    ! Exception: target ice thickness at the floating calving front
    ! is often wrong (because of the difficulty of remapping a discontinuous
    ! field), so instead use the mean of the neighbouring non-front shelf
    ! vertices.
    call gather_to_all( ice%mask_floating_ice, mask_floating_ice_tot)
    call gather_to_all( ice%mask_cf_fl       , mask_cf_fl_tot)
    call gather_to_all( Hi_target_corr       , Hi_target_tot)

    do vi = mesh%vi1, mesh% vi2
      if (mask_cf_fl_tot( vi)) then
        w_sum  = 0._dp
        wH_sum = 0._dp
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          if (mask_floating_ice_tot( vj) .and. .not. mask_cf_fl_tot( vj)) then
            w_sum = w_sum + 1._dp
            wH_sum = wH_sum + Hi_target_tot( vj)
          end if
        end do
        if (w_sum > 0._dp) then
          Hi_target_corr( vi) = wH_sum / w_sum
        end if
      end if
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calc_corrected_target_thickness

  subroutine calc_dTdt( mesh, grid_smooth, ice, Hi_target_corr, target_mask_shelf, dTdt)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_grid),                        intent(in   ) :: grid_smooth
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hi_target_corr
    logical , dimension(mesh%vi1:mesh%vi2), intent(in   ) :: target_mask_shelf
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: dTdt

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_dTdt'
    integer,  dimension(mesh%vi1:mesh%vi2) :: mask_extrapolation
    integer                                :: vi
    real(dp)                               :: deltaH, dHdt, dBMBdt
    real(dp), parameter                    :: c_H     = 0.00001_dp
    real(dp), parameter                    :: c_dHdt  = 0.0003_dp

    ! Add routine to call stack
    call init_routine( routine_name)

    mask_extrapolation = 1
    dTdt               = 0._dp

    do vi = mesh%vi1, mesh%vi2

      ! Only apply nudging to fully floating shelf vertices,
      ! skipping the grounding line and calving front.
      if (ice%fraction_gr( vi) < 0.01_dp .and. ice%Hi( vi) > 0.1_dp .and. .not. ice%mask_margin( vi)) then

        mask_extrapolation( vi) = 2

        deltaH = ice%Hi( vi) - Hi_target_corr( vi)
        dHdt   = ice%dHi_dt( vi)

        dTdt( vi) = c_H * deltaH + c_dHdt * dHdt

      end if
    end do

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    call extrapolate_Gaussian( mesh, mask_extrapolation, dTdt, 10e3_dp)

    call smooth_dTdt( mesh, grid_smooth, dTdt)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calc_dTdt

  subroutine smooth_dTdt( mesh, grid_smooth, dTdt)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_grid),                        intent(in   ) :: grid_smooth
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: dTdt

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'smooth_dTdt'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: dTdt_smoothed
    real(dp), parameter                    :: w_smooth = 0.5_dp

    ! Add routine to path
    call init_routine( routine_name)

    dTdt_smoothed = dTdt
    call smooth_Gaussian( mesh, grid_smooth, C%output_dir, dTdt_smoothed, 20e3_dp)

    dTdt = (1._dp - w_smooth) * dTdt + w_smooth * dTdt_smoothed

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine smooth_dTdt

  subroutine map_deltaT_to_reference_grid( mesh, snapshot_nudge2D)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ocean_model_snapshot_nudge2D), intent(inout) :: snapshot_nudge2D

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_deltaT_to_reference_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Map deltaT to reference grid
    call map_from_mesh_vertices_to_xy_grid_2D( mesh, snapshot_nudge2D%grid_ref, trim( C%output_dir), &
      snapshot_nudge2D%deltaT_nudge, snapshot_nudge2D%deltaT_nudge_grid, '2nd_order_conservative', &
      d_mesh_is_hybrid = .false., d_grid_is_hybrid = .false.)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine map_deltaT_to_reference_grid

  subroutine add_deltaT_to_snapshot( mesh, snapshot_nudge2D)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ocean_model_snapshot_nudge2D), intent(inout) :: snapshot_nudge2D

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_deltaT_to_snapshot'
    integer                        :: vi,k,n

    ! Add routine to call stack
    call init_routine( routine_name)

    ! On the mesh...
    do vi = mesh%vi1, mesh%vi2
      do k = 1, C%nz_ocean
        snapshot_nudge2D%T( vi,k) = snapshot_nudge2D%T_ref( vi,k) + snapshot_nudge2D%deltaT_nudge( vi)
        snapshot_nudge2D%S( vi,k) = snapshot_nudge2D%S_ref( vi,k)
      end do
    end do

    ! ...and on the reference grid
    do n = snapshot_nudge2D%grid_ref%n1, snapshot_nudge2D%grid_ref%n2
      do k = 1, snapshot_nudge2D%ndepth_ref
        snapshot_nudge2D%T_grid( n,k) = snapshot_nudge2D%T_ref_grid( n,k) + snapshot_nudge2D%deltaT_nudge_grid( n)
        snapshot_nudge2D%S_grid( n,k) = snapshot_nudge2D%S_ref_grid( n,k)
      end do
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine add_deltaT_to_snapshot

  subroutine set_applied_ocean_data( mesh, ocean)

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    type(type_ocean_model), intent(inout) :: ocean

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_applied_ocean_data'

    ! Add routine to call stack
    call init_routine( routine_name)

    ocean%T = ocean%snapshot_nudge2D%T
    ocean%S = ocean%snapshot_nudge2D%S

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_applied_ocean_data

  subroutine write_to_output_file( mesh, snapshot_nudge2D)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ocean_model_snapshot_nudge2D), intent(inout) :: snapshot_nudge2D

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_output_file'
    integer                        :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    call open_existing_netcdf_file_for_writing( snapshot_nudge2D%output_filename, ncid)

    call write_to_field_multopt_grid_dp_3D_ocean_notime( snapshot_nudge2D%grid_ref, &
      snapshot_nudge2D%output_filename, ncid, 't_an', snapshot_nudge2D%T_grid)
    call write_to_field_multopt_grid_dp_3D_ocean_notime( snapshot_nudge2D%grid_ref, &
      snapshot_nudge2D%output_filename, ncid, 's_an', snapshot_nudge2D%S_grid)

    call close_netcdf_file( ncid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_output_file

  subroutine initialise_ocean_model_snapshot_nudge2D( mesh, snapshot_nudge2D, region_name, refgeo_PD, refgeo_init)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ocean_model_snapshot_nudge2D), intent(inout) :: snapshot_nudge2D
    character(len=3),                        intent(in   ) :: region_name
    type(type_reference_geometry),           intent(in   ) :: refgeo_PD, refgeo_init

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'initialise_ocean_model_snapshot_nudge2D'
    character(len=1024)                   :: filename
    integer                               :: ncid
    real(dp), dimension(:,:), allocatable :: T_ref_partial_raw_layers
    real(dp), dimension(:,:), allocatable :: S_ref_partial_raw_layers

    ! Add routine to call stack
    call init_routine( routine_name)

    select case (region_name)
      case ('NAM')
        filename = C%filename_ocean_snapshot_NAM
      case ('EAS')
        filename = C%filename_ocean_snapshot_EAS
      case ('GRL')
        filename = C%filename_ocean_snapshot_GRL
      case ('ANT')
        filename = C%filename_ocean_snapshot_ANT
      case default
        call crash('unknown region_name "' // region_name // '"')
    end select

    ! Allocate memory for meshed ocean data
    allocate( snapshot_nudge2D%T_ref       ( mesh%vi1:mesh%vi2, C%nz_ocean), source = 0._dp)
    allocate( snapshot_nudge2D%S_ref       ( mesh%vi1:mesh%vi2, C%nz_ocean), source = 0._dp)
    allocate( snapshot_nudge2D%deltaT_nudge( mesh%vi1:mesh%vi2            ), source = 0._dp)
    allocate( snapshot_nudge2D%T           ( mesh%vi1:mesh%vi2, C%nz_ocean), source = 0._dp)
    allocate( snapshot_nudge2D%S           ( mesh%vi1:mesh%vi2, C%nz_ocean), source = 0._dp)

    ! Set up the grid from the file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_xy_grid_from_file( filename, ncid, snapshot_nudge2D%grid_ref)
    call setup_depth_from_file( filename, ncid, snapshot_nudge2D%ndepth_ref, snapshot_nudge2D%depth_ref)
    call close_netcdf_file( ncid)

    ! Allocate memory for gridded 3-D ocean snapshot data
    allocate( snapshot_nudge2D%T_ref_grid       ( snapshot_nudge2D%grid_ref%n1: snapshot_nudge2D%grid_ref%n2, snapshot_nudge2D%ndepth_ref), source = 0._dp)
    allocate( snapshot_nudge2D%S_ref_grid       ( snapshot_nudge2D%grid_ref%n1: snapshot_nudge2D%grid_ref%n2, snapshot_nudge2D%ndepth_ref), source = 0._dp)
    allocate( snapshot_nudge2D%deltaT_nudge_grid( snapshot_nudge2D%grid_ref%n1: snapshot_nudge2D%grid_ref%n2                             ), source = 0._dp)
    allocate( snapshot_nudge2D%T_grid           ( snapshot_nudge2D%grid_ref%n1: snapshot_nudge2D%grid_ref%n2, snapshot_nudge2D%ndepth_ref), source = 0._dp)
    allocate( snapshot_nudge2D%S_grid           ( snapshot_nudge2D%grid_ref%n1: snapshot_nudge2D%grid_ref%n2, snapshot_nudge2D%ndepth_ref), source = 0._dp)

    ! Read gridded 3-D ocean snapshot data
    call read_field_from_xy_file_dp_3D_ocean( filename, field_name_options_T_ocean, snapshot_nudge2D%T_ref_grid)
    call read_field_from_xy_file_dp_3D_ocean( filename, field_name_options_S_ocean, snapshot_nudge2D%S_ref_grid)

    ! Allocate memory for 3-D ocean snapshot data on the mesh, but with the original vertical layers
    allocate( T_ref_partial_raw_layers( mesh%vi1:mesh%vi2, snapshot_nudge2D%ndepth_ref))
    allocate( S_ref_partial_raw_layers( mesh%vi1:mesh%vi2, snapshot_nudge2D%ndepth_ref))

    ! Remap data horizontally
    call map_from_xy_grid_to_mesh_3D( snapshot_nudge2D%grid_ref, mesh, trim( C%output_dir), snapshot_nudge2D%T_ref_grid, T_ref_partial_raw_layers)
    call map_from_xy_grid_to_mesh_3D( snapshot_nudge2D%grid_ref, mesh, trim( C%output_dir), snapshot_nudge2D%S_ref_grid, S_ref_partial_raw_layers)

    ! Remap data vertically
    call map_from_vertical_to_vertical_2D_ocean( mesh, snapshot_nudge2D%depth_ref, C%z_ocean, T_ref_partial_raw_layers, snapshot_nudge2D%T_ref)
    call map_from_vertical_to_vertical_2D_ocean( mesh, snapshot_nudge2D%depth_ref, C%z_ocean, S_ref_partial_raw_layers, snapshot_nudge2D%S_ref)

    call set_target_geometry( mesh, snapshot_nudge2D, refgeo_PD, refgeo_init)

    call create_output_file( snapshot_nudge2D)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model_snapshot_nudge2D

  subroutine set_target_geometry( mesh, snapshot_nudge2D, refgeo_PD, refgeo_init)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ocean_model_snapshot_nudge2D), intent(inout) :: snapshot_nudge2D
    type(type_reference_geometry),           intent(in   ) :: refgeo_PD, refgeo_init

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_target_geometry'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    call reallocate_reference_geometry_on_mesh( mesh, snapshot_nudge2D%target_geometry)

    select case (C%choice_inversion_target_geometry)
    case default
      call crash('unknown choice_inversion_target_geometry "' // trim(C%choice_inversion_target_geometry) // '"')
    case ('init')
      snapshot_nudge2D%target_geometry%Hi = refgeo_init%Hi
      snapshot_nudge2D%target_geometry%Hb = refgeo_init%Hb
      snapshot_nudge2D%target_geometry%Hs = refgeo_init%Hs
      snapshot_nudge2D%target_geometry%SL = refgeo_init%SL
    case ('PD')
      snapshot_nudge2D%target_geometry%Hi = refgeo_PD%Hi
      snapshot_nudge2D%target_geometry%Hb = refgeo_PD%Hb
      snapshot_nudge2D%target_geometry%Hs = refgeo_PD%Hs
      snapshot_nudge2D%target_geometry%SL = refgeo_PD%SL
    end select

    ! Determine the shelf mask of the target geometry
    allocate( snapshot_nudge2D%target_mask_shelf( mesh%vi1:mesh%vi2), source = .false.)
    do vi = mesh%vi1, mesh%vi2
      if (snapshot_nudge2D%target_geometry%Hi( vi) > 0.1_dp) then
        snapshot_nudge2D%target_mask_shelf( vi) = is_floating( &
          snapshot_nudge2D%target_geometry%Hi( vi), &
          snapshot_nudge2D%target_geometry%Hb( vi), &
          snapshot_nudge2D%target_geometry%SL( vi))
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_target_geometry

  subroutine create_output_file( snapshot_nudge2D)

    ! In/output variables:
    type(type_ocean_model_snapshot_nudge2D), intent(inout) :: snapshot_nudge2D

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'create_ocean_model_snapshot_nudge2D_output_file'
    character(len=1024)                   :: filename
    integer                               :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    filename = trim( C%output_dir) // 'ocean_snapshot_nudged.nc'
    snapshot_nudge2D%output_filename = filename

    ! Print to terminal
    if (par%primary) write(0,'(a)') '     Creating ocean snapshot+nudge2D output file "' // &
      colour_string( trim( filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( filename, ncid)

    ! Set up the grid in the file
    call setup_xy_grid_in_netcdf_file( filename, ncid, snapshot_nudge2D%grid_ref)

    ! Add a time dimension to the file
    call add_time_dimension_to_file( filename, ncid)

    ! Add a depth dimension to the file
    call add_depth_dimension_to_file( filename, ncid, snapshot_nudge2D%depth_ref)

    ! Add the data fields to the file
    call add_field_grid_dp_3D_ocean_notime( filename, ncid, 't_an' , long_name = 'Ocean temperature', units = 'degrees C')
    call add_field_grid_dp_3D_ocean_notime( filename, ncid, 's_an' , long_name = 'Ocean salinity', units = 'PSU')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_output_file

end module ocean_snapshot_nudge2D
