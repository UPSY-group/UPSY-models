module ocean_snapshot_nudge2D_plus_anomalies

  use precisions, only: dp
  use parameters, only: NaN
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning, colour_string
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use grid_types, only: type_grid
  use ice_model_types, only: type_ice_model
  use ocean_model_types, only: type_ocean_model, type_ocean_model_snapshot_nudge2D_plus_anomalies
  use netcdf_io_main
  use remapping_main
  use reference_geometry_types, only: type_reference_geometry
  use reference_geometries_main, only: reallocate_reference_geometry_on_mesh
  use ice_geometry_basics, only: is_floating
  use mpi_distributed_memory, only: gather_to_all
  use mesh_utilities, only: extrapolate_Gaussian
  use mesh_data_smoothing, only: smooth_Gaussian
  use mpi_f08, only: MPI_BCAST, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD

  implicit none

  private

  public :: initialise_ocean_model_snapshot_nudge2D_plus_anomalies, run_ocean_model_snapshot_nudge2D_plus_anomalies

contains

  subroutine run_ocean_model_snapshot_nudge2D_plus_anomalies( mesh, grid_smooth, ice, ocean, time)

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    type(type_grid),        intent(in   ) :: grid_smooth
    type(type_ice_model),   intent(in   ) :: ice
    type(type_ocean_model), intent(inout) :: ocean
    real(dp),               intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_ocean_model_snapshot_nudge2D_plus_anomalies'
    real(dp)                       :: w0, w1

    ! Add routine to call stack
    call init_routine( routine_name)

    if (time >= C%BMB_inversion_t_start .and.  time <= C%BMB_inversion_t_end) then
      call nudge_deltaT( mesh, grid_smooth, ice, ocean%snapshot_nudge2D_plus_anomalies)
      call map_deltaT_to_reference_grid( mesh, ocean%snapshot_nudge2D_plus_anomalies)
      call add_deltaT_to_snapshot( mesh, ocean%snapshot_nudge2D_plus_anomalies)
    end if

    if (time >= C%ocean_snp_p_anml_tstart_anomalies) then

      ! If the current model time falls outside the enveloping window
      ! of the two timeframes that have been read, update them
      if (time < ocean%snapshot_nudge2D_plus_anomalies%anomaly_t0 .or. &
          time > ocean%snapshot_nudge2D_plus_anomalies%anomaly_t1) then
        call update_timeframes( mesh, ocean%snapshot_nudge2D_plus_anomalies, time)
      end if

      ! Interpolate between the two timeframes to find the applied anomaly
      w0 = (ocean%snapshot_nudge2D_plus_anomalies%anomaly_t1 - time) / &
          (ocean%snapshot_nudge2D_plus_anomalies%anomaly_t1 - ocean%snapshot_nudge2D_plus_anomalies%anomaly_t0)
      w1 = 1._dp - w0

      ocean%snapshot_nudge2D_plus_anomalies%T_anomaly = &
        w0 * ocean%snapshot_nudge2D_plus_anomalies%T_anomaly_0 + &
        w1 * ocean%snapshot_nudge2D_plus_anomalies%T_anomaly_1
      ocean%snapshot_nudge2D_plus_anomalies%S_anomaly = &
        w0 * ocean%snapshot_nudge2D_plus_anomalies%S_anomaly_0 + &
        w1 * ocean%snapshot_nudge2D_plus_anomalies%S_anomaly_1

      if (time >= C%ocean_snp_p_anml_tstart_anomalies) then
        ! Add anomaly to snapshot to find the applied ocean state
        ocean%snapshot_nudge2D_plus_anomalies%T = &
          ocean%snapshot_nudge2D_plus_anomalies%T_baseline + &
          ocean%snapshot_nudge2D_plus_anomalies%T_anomaly
        ocean%snapshot_nudge2D_plus_anomalies%S = &
          ocean%snapshot_nudge2D_plus_anomalies%S_baseline + &
          ocean%snapshot_nudge2D_plus_anomalies%S_anomaly
      else
        ! Only apply the baseline
        ocean%snapshot_nudge2D_plus_anomalies%T = ocean%snapshot_nudge2D_plus_anomalies%T_baseline
        ocean%snapshot_nudge2D_plus_anomalies%S = ocean%snapshot_nudge2D_plus_anomalies%S_baseline
      end if

    end if

    ! Set applied ocean data
    ocean%T = ocean%snapshot_nudge2D_plus_anomalies%T
    ocean%S = ocean%snapshot_nudge2D_plus_anomalies%S

    ! call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_ocean_model_snapshot_nudge2D_plus_anomalies

  subroutine nudge_deltaT( mesh, grid_smooth, ice, snapshot_nudge2D_plus_anomalies)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_grid),                         intent(in   ) :: grid_smooth
    type(type_ice_model),                    intent(in   ) :: ice
    type(type_ocean_model_snapshot_nudge2D_plus_anomalies), intent(inout) :: snapshot_nudge2D_plus_anomalies

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'nudge_deltaT'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: Hi_target_corr
    real(dp), dimension(mesh%vi1:mesh%vi2) :: dTdt
    integer                                :: vi
    real(dp), parameter                    :: deltaT_max =  2._dp
    real(dp), parameter                    :: deltaT_min = -2._dp

    ! Add routine to call stack
    call init_routine( routine_name)

    call calc_corrected_target_thickness( mesh, ice, snapshot_nudge2D_plus_anomalies, Hi_target_corr)
    call calc_dTdt( mesh, grid_smooth, ice, Hi_target_corr, snapshot_nudge2D_plus_anomalies%target_mask_shelf, dTdt)

    do vi = mesh%vi1, mesh%vi2

      snapshot_nudge2D_plus_anomalies%deltaT_nudge( vi) = snapshot_nudge2D_plus_anomalies%deltaT_nudge( vi) + C%dt_ocean * dTdt( vi)

      ! Apply limits
      snapshot_nudge2D_plus_anomalies%deltaT_nudge( vi) = max( deltaT_min, min( deltaT_max, &
        snapshot_nudge2D_plus_anomalies%deltaT_nudge( vi)))

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine nudge_deltaT

  subroutine calc_corrected_target_thickness( mesh, ice, snapshot_nudge2D_plus_anomalies, Hi_target_corr)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ice_model),                    intent(in   ) :: ice
    type(type_ocean_model_snapshot_nudge2D_plus_anomalies), intent(in   ) :: snapshot_nudge2D_plus_anomalies
    real(dp), dimension(mesh%vi1:mesh%vi2),  intent(  out) :: Hi_target_corr

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_corrected_target_thickness'
    integer                        :: vi, ci, vj
    logical,  dimension(mesh%nV)   :: mask_floating_ice_tot, mask_cf_fl_tot
    real(dp), dimension(mesh%nV)   :: Hi_target_tot
    real(dp)                       :: w_sum, wH_sum

    ! Add routine to call stack
    call init_routine( routine_name)

    Hi_target_corr = snapshot_nudge2D_plus_anomalies%target_geometry%Hi

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

  subroutine map_deltaT_to_reference_grid( mesh, snapshot_nudge2D_plus_anomalies)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ocean_model_snapshot_nudge2D_plus_anomalies), intent(inout) :: snapshot_nudge2D_plus_anomalies

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_deltaT_to_reference_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Map deltaT to reference grid
    call map_from_mesh_vertices_to_xy_grid_2D( mesh, snapshot_nudge2D_plus_anomalies%grid_ref, trim( C%output_dir), &
      snapshot_nudge2D_plus_anomalies%deltaT_nudge, snapshot_nudge2D_plus_anomalies%deltaT_nudge_grid, '2nd_order_conservative', &
      d_mesh_is_hybrid = .false., d_grid_is_hybrid = .false.)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine map_deltaT_to_reference_grid

  subroutine add_deltaT_to_snapshot( mesh, snapshot_nudge2D_plus_anomalies)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ocean_model_snapshot_nudge2D_plus_anomalies), intent(inout) :: snapshot_nudge2D_plus_anomalies

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_deltaT_to_snapshot'
    integer                        :: vi,k,n

    ! Add routine to call stack
    call init_routine( routine_name)

    ! On the mesh...
    do vi = mesh%vi1, mesh%vi2
      do k = 1, C%nz_ocean
        snapshot_nudge2D_plus_anomalies%T_baseline( vi,k) = snapshot_nudge2D_plus_anomalies%T_ref( vi,k) + snapshot_nudge2D_plus_anomalies%deltaT_nudge( vi)
        snapshot_nudge2D_plus_anomalies%S_baseline( vi,k) = snapshot_nudge2D_plus_anomalies%S_ref( vi,k)
      end do
    end do

    ! ...and on the reference grid
    do n = snapshot_nudge2D_plus_anomalies%grid_ref%n1, snapshot_nudge2D_plus_anomalies%grid_ref%n2
      do k = 1, snapshot_nudge2D_plus_anomalies%ndepth_ref
        snapshot_nudge2D_plus_anomalies%T_grid( n,k) = snapshot_nudge2D_plus_anomalies%T_ref_grid( n,k) + snapshot_nudge2D_plus_anomalies%deltaT_nudge_grid( n)
        snapshot_nudge2D_plus_anomalies%S_grid( n,k) = snapshot_nudge2D_plus_anomalies%S_ref_grid( n,k)
      end do
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine add_deltaT_to_snapshot

  subroutine initialise_ocean_model_snapshot_nudge2D_plus_anomalies( mesh, &
    snapshot_nudge2D_plus_anomalies, region_name, refgeo_PD, refgeo_init)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ocean_model_snapshot_nudge2D_plus_anomalies), intent(inout) :: snapshot_nudge2D_plus_anomalies
    character(len=3),                        intent(in   ) :: region_name
    type(type_reference_geometry),           intent(in   ) :: refgeo_PD, refgeo_init

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'initialise_ocean_model_snapshot_nudge2D_plus_anomalies'
    character(len=1024)                   :: filename
    integer                               :: ncid
    real(dp), dimension(:,:), allocatable :: T_ref_partial_raw_layers
    real(dp), dimension(:,:), allocatable :: S_ref_partial_raw_layers

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety
    if (C%BMB_inversion_t_end > C%ocean_snp_p_anml_tstart_anomalies) then
      call crash('can only nudge ocean before applying anomalies')
    end if

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
    allocate( snapshot_nudge2D_plus_anomalies%T_ref       ( mesh%vi1:mesh%vi2, C%nz_ocean), source = 0._dp)
    allocate( snapshot_nudge2D_plus_anomalies%S_ref       ( mesh%vi1:mesh%vi2, C%nz_ocean), source = 0._dp)
    allocate( snapshot_nudge2D_plus_anomalies%deltaT_nudge( mesh%vi1:mesh%vi2            ), source = 0._dp)
    allocate( snapshot_nudge2D_plus_anomalies%T_baseline  ( mesh%vi1:mesh%vi2, C%nz_ocean), source = 0._dp)
    allocate( snapshot_nudge2D_plus_anomalies%S_baseline  ( mesh%vi1:mesh%vi2, C%nz_ocean), source = 0._dp)

    ! Set up the grid from the file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_xy_grid_from_file( filename, ncid, snapshot_nudge2D_plus_anomalies%grid_ref)
    call setup_depth_from_file( filename, ncid, snapshot_nudge2D_plus_anomalies%ndepth_ref, snapshot_nudge2D_plus_anomalies%depth_ref)
    call close_netcdf_file( ncid)

    ! Allocate memory for gridded 3-D ocean snapshot data
    allocate( snapshot_nudge2D_plus_anomalies%T_ref_grid       ( snapshot_nudge2D_plus_anomalies%grid_ref%n1: snapshot_nudge2D_plus_anomalies%grid_ref%n2, snapshot_nudge2D_plus_anomalies%ndepth_ref), source = 0._dp)
    allocate( snapshot_nudge2D_plus_anomalies%S_ref_grid       ( snapshot_nudge2D_plus_anomalies%grid_ref%n1: snapshot_nudge2D_plus_anomalies%grid_ref%n2, snapshot_nudge2D_plus_anomalies%ndepth_ref), source = 0._dp)
    allocate( snapshot_nudge2D_plus_anomalies%deltaT_nudge_grid( snapshot_nudge2D_plus_anomalies%grid_ref%n1: snapshot_nudge2D_plus_anomalies%grid_ref%n2                             ), source = 0._dp)
    allocate( snapshot_nudge2D_plus_anomalies%T_grid           ( snapshot_nudge2D_plus_anomalies%grid_ref%n1: snapshot_nudge2D_plus_anomalies%grid_ref%n2, snapshot_nudge2D_plus_anomalies%ndepth_ref), source = 0._dp)
    allocate( snapshot_nudge2D_plus_anomalies%S_grid           ( snapshot_nudge2D_plus_anomalies%grid_ref%n1: snapshot_nudge2D_plus_anomalies%grid_ref%n2, snapshot_nudge2D_plus_anomalies%ndepth_ref), source = 0._dp)

    ! Read gridded 3-D ocean snapshot data
    call read_field_from_xy_file_dp_3D_ocean( filename, field_name_options_T_ocean, snapshot_nudge2D_plus_anomalies%T_ref_grid)
    call read_field_from_xy_file_dp_3D_ocean( filename, field_name_options_S_ocean, snapshot_nudge2D_plus_anomalies%S_ref_grid)

    ! Allocate memory for 3-D ocean snapshot data on the mesh, but with the original vertical layers
    allocate( T_ref_partial_raw_layers( mesh%vi1:mesh%vi2, snapshot_nudge2D_plus_anomalies%ndepth_ref))
    allocate( S_ref_partial_raw_layers( mesh%vi1:mesh%vi2, snapshot_nudge2D_plus_anomalies%ndepth_ref))

    ! Remap data horizontally
    call map_from_xy_grid_to_mesh_3D( snapshot_nudge2D_plus_anomalies%grid_ref, mesh, trim( C%output_dir), snapshot_nudge2D_plus_anomalies%T_ref_grid, T_ref_partial_raw_layers)
    call map_from_xy_grid_to_mesh_3D( snapshot_nudge2D_plus_anomalies%grid_ref, mesh, trim( C%output_dir), snapshot_nudge2D_plus_anomalies%S_ref_grid, S_ref_partial_raw_layers)

    ! Remap data vertically
    call map_from_vertical_to_vertical_2D_ocean( mesh, snapshot_nudge2D_plus_anomalies%depth_ref, C%z_ocean, T_ref_partial_raw_layers, snapshot_nudge2D_plus_anomalies%T_ref)
    call map_from_vertical_to_vertical_2D_ocean( mesh, snapshot_nudge2D_plus_anomalies%depth_ref, C%z_ocean, S_ref_partial_raw_layers, snapshot_nudge2D_plus_anomalies%S_ref)

    call set_target_geometry( mesh, snapshot_nudge2D_plus_anomalies, refgeo_PD, refgeo_init)



    ! Stuff for anomalies
    ! ===================

    ! Two anomaly snapshots enveloping the current model time
    allocate( snapshot_nudge2D_plus_anomalies%T_anomaly_0( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate( snapshot_nudge2D_plus_anomalies%S_anomaly_0( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)

    allocate( snapshot_nudge2D_plus_anomalies%T_anomaly_1( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate( snapshot_nudge2D_plus_anomalies%S_anomaly_1( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)

    ! Time-weighted anomaly
    allocate( snapshot_nudge2D_plus_anomalies%T_anomaly( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate( snapshot_nudge2D_plus_anomalies%S_anomaly( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)

    ! Applied climate
    allocate( snapshot_nudge2D_plus_anomalies%T( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate( snapshot_nudge2D_plus_anomalies%S( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)

    ! Initialise anomaly timeframes
    snapshot_nudge2D_plus_anomalies%anomaly_t0 = C%start_time_of_run - 200._dp
    snapshot_nudge2D_plus_anomalies%anomaly_t1 = C%start_time_of_run - 100._dp
    call update_timeframes( mesh, snapshot_nudge2D_plus_anomalies, C%start_time_of_run)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model_snapshot_nudge2D_plus_anomalies

  subroutine set_target_geometry( mesh, snapshot_nudge2D_plus_anomalies, refgeo_PD, refgeo_init)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ocean_model_snapshot_nudge2D_plus_anomalies), intent(inout) :: snapshot_nudge2D_plus_anomalies
    type(type_reference_geometry),           intent(in   ) :: refgeo_PD, refgeo_init

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_target_geometry'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    call reallocate_reference_geometry_on_mesh( mesh, snapshot_nudge2D_plus_anomalies%target_geometry)

    select case (C%choice_inversion_target_geometry)
    case default
      call crash('unknown choice_inversion_target_geometry "' // trim(C%choice_inversion_target_geometry) // '"')
    case ('init')
      snapshot_nudge2D_plus_anomalies%target_geometry%Hi = refgeo_init%Hi
      snapshot_nudge2D_plus_anomalies%target_geometry%Hb = refgeo_init%Hb
      snapshot_nudge2D_plus_anomalies%target_geometry%Hs = refgeo_init%Hs
      snapshot_nudge2D_plus_anomalies%target_geometry%SL = refgeo_init%SL
    case ('PD')
      snapshot_nudge2D_plus_anomalies%target_geometry%Hi = refgeo_PD%Hi
      snapshot_nudge2D_plus_anomalies%target_geometry%Hb = refgeo_PD%Hb
      snapshot_nudge2D_plus_anomalies%target_geometry%Hs = refgeo_PD%Hs
      snapshot_nudge2D_plus_anomalies%target_geometry%SL = refgeo_PD%SL
    end select

    ! Determine the shelf mask of the target geometry
    allocate( snapshot_nudge2D_plus_anomalies%target_mask_shelf( mesh%vi1:mesh%vi2), source = .false.)
    do vi = mesh%vi1, mesh%vi2
      if (snapshot_nudge2D_plus_anomalies%target_geometry%Hi( vi) > 0.1_dp) then
        snapshot_nudge2D_plus_anomalies%target_mask_shelf( vi) = is_floating( &
          snapshot_nudge2D_plus_anomalies%target_geometry%Hi( vi), &
          snapshot_nudge2D_plus_anomalies%target_geometry%Hb( vi), &
          snapshot_nudge2D_plus_anomalies%target_geometry%SL( vi))
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_target_geometry

  subroutine update_timeframes( mesh, snapshot_nudge2D_plus_anomalies, time)

    ! In/output variables:
    type(type_mesh),                                intent(in   ) :: mesh
    type(type_ocean_model_snapshot_nudge2D_plus_anomalies), intent(inout) :: snapshot_nudge2D_plus_anomalies
    real(dp),                                       intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'update_timeframes'
    character(len=1024)                 :: filename
    integer                             :: ncid, id_dim_time, nt, id_var_time, ierr
    real(dp), dimension(:), allocatable :: time_from_file
    integer                             :: ti0, ti1

    ! Add routine to path
    call init_routine( routine_name)

    filename = trim( C%ocean_snp_p_anml_filename_anomalies)

    ! Read time variable from the file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call check_time( filename, ncid)
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)
    call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)
    allocate( time_from_file( nt))
    call read_var_primary( filename, ncid, id_var_time, time_from_file)
    call MPI_BCAST( time_from_file(:), nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call close_netcdf_file( ncid)

    ! Find the two timeframes
    if (time < time_from_file( 1)) then
      if (par%primary) call warning('model time before start of anomaly file; using first timeframe')
      ti0 = 1
      ti1 = 2
    elseif (time > time_from_file( size( time_from_file,1))) then
      if (par%primary) call warning('model time beyond end of anomaly file; using last timeframe')
      ti0 = size( time_from_file,1) - 1
      ti1 = size( time_from_file,1)
    else
      ti0 = 1
      ti1 = 2
      do while (time_from_file( ti1) < time .and. ti1 < size( time_from_file,1))
        ti0 = ti1
        ti1 = ti1 + 1
      end do
    end if

    snapshot_nudge2D_plus_anomalies%anomaly_t0 = time_from_file( ti0)
    snapshot_nudge2D_plus_anomalies%anomaly_t1 = time_from_file( ti1)

    ! Read the two timeframes
    call read_field_from_file_3D_ocean( filename, 'temperature_anomaly', &
      mesh, C%output_dir, C%z_ocean, snapshot_nudge2D_plus_anomalies%T_anomaly_0, &
      time_to_read = snapshot_nudge2D_plus_anomalies%anomaly_t0)
    call read_field_from_file_3D_ocean( filename, 'temperature_anomaly', &
      mesh, C%output_dir, C%z_ocean, snapshot_nudge2D_plus_anomalies%T_anomaly_1, &
      time_to_read = snapshot_nudge2D_plus_anomalies%anomaly_t1)
    call read_field_from_file_3D_ocean( filename, 'salinity_anomaly', &
      mesh, C%output_dir, C%z_ocean, snapshot_nudge2D_plus_anomalies%S_anomaly_0, &
      time_to_read = snapshot_nudge2D_plus_anomalies%anomaly_t0)
    call read_field_from_file_3D_ocean( filename, 'salinity_anomaly', &
      mesh, C%output_dir, C%z_ocean, snapshot_nudge2D_plus_anomalies%S_anomaly_1, &
      time_to_read = snapshot_nudge2D_plus_anomalies%anomaly_t1)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_timeframes

end module ocean_snapshot_nudge2D_plus_anomalies
