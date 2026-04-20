MODULE climate_snapshot_plus_anomalies

  ! Snapshot + dT(time) climate model

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  use mpi_f08                                                , only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, MPI_BCAST
  USE call_stack_and_comp_time_tracking                  , ONLY: crash, init_routine, finalise_routine, warning
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model, type_climate_model_snapshot
  USE global_forcing_types                                   , ONLY: type_global_forcing
  use climate_realistic                                      , only: initialise_climate_model_realistic, initialise_insolation_forcing, remap_insolation
  USE global_forcings_main
  USE netcdf_io_main
  USE netcdf_basic
  use reallocate_mod                                         , only: reallocate_bounds
  use mpi_distributed_memory                                 , only: distribute_from_primary
  use climate_model_utilities                                , only: get_insolation_at_time, apply_geometry_downscaling_corrections
  use series_utilities

  IMPLICIT NONE

  private

  public :: run_climate_model_snp_p_anml
  public :: initialise_climate_model_snp_p_anml
  public :: remap_climate_snp_p_anml


CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_climate_model_snp_p_anml( mesh, ice, climate, time)
    ! Calculate the climate
    !
    ! Use a snapshot plus a prescribed uniform deltaT

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_climate_model_snp_p_anml'
    INTEGER                                               :: vi, m
    REAL(dp)                                              :: w0, w1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen snapshot+deltaT(time) climate model
    ! If the current model time falls outside the enveloping window
    ! of the two timeframes that have been read, update them
    if (time < climate%snapshot_p_anml%anomaly_t0 .or. time > climate%snapshot_p_anml%anomaly_t1) then
      call update_climate_timeframes( mesh, climate, time)
    end if

    ! Interpolate between the two timeframes to find the applied anomaly
    call interpolate_between_climate_anomaly_timeframes(mesh, climate, time)

    ! Update temperature and precipitation fields based on the mismatch between
    ! the ice sheet surface elevation in the forcing climate and the model's ice sheet surface elevation
    call apply_geometry_downscaling_corrections(mesh, ice, climate, climate%snapshot_p_anml%snapshot_baseline, 0.0_dp)

    ! if needed for IMAU-ITM or climate matrix, we need to update insolation
    IF (climate%snapshot%has_insolation) THEN
      CALL get_insolation_at_time( mesh, time, climate%snapshot_p_anml%snapshot_baseline)
      climate%Q_TOA          = climate%snapshot_p_anml%snapshot_baseline%Q_TOA
      climate%snapshot%Q_TOA = climate%snapshot_p_anml%snapshot_baseline%Q_TOA
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_snp_p_anml

  SUBROUTINE initialise_climate_model_snp_p_anml( mesh, ice, climate, region_name)
    ! Initialise the climate model
    !
    ! Use a realistic climate scheme

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_climate_model_snp_p_anml'
    INTEGER                                               :: vi, m
    CHARACTER(LEN=256)                                    :: filename_climate_snapshot
    REAL(dp)                                              :: timeframe_init_insolation, w0, w1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate variables that will be computed
    ALLOCATE( climate%snapshot_p_anml%snapshot_baseline%Hs(     mesh%vi1:mesh%vi2))
    ALLOCATE( climate%snapshot_p_anml%snapshot_baseline%T2m(    mesh%vi1:mesh%vi2,12))
    ALLOCATE( climate%snapshot_p_anml%snapshot_baseline%Precip( mesh%vi1:mesh%vi2,12))
    ALLOCATE( climate%snapshot_p_anml%T2m_anomaly_0(    mesh%vi1:mesh%vi2,12))
    ALLOCATE( climate%snapshot_p_anml%Precip_anomaly_0( mesh%vi1:mesh%vi2,12))
    ALLOCATE( climate%snapshot_p_anml%T2m_anomaly_1(    mesh%vi1:mesh%vi2,12))
    ALLOCATE( climate%snapshot_p_anml%Precip_anomaly_1( mesh%vi1:mesh%vi2,12))
    ALLOCATE( climate%snapshot_p_anml%T2m_anomaly(      mesh%vi1:mesh%vi2,12))
    ALLOCATE( climate%snapshot_p_anml%Precip_anomaly(   mesh%vi1:mesh%vi2,12))
    ALLOCATE( climate%snapshot_p_anml%T2m(              mesh%vi1:mesh%vi2,12))
    ALLOCATE( climate%snapshot_p_anml%Precip(           mesh%vi1:mesh%vi2,12))
    climate%snapshot_p_anml%snapshot_baseline%Hs     = 0._dp
    climate%snapshot_p_anml%snapshot_baseline%T2m    = 0._dp
    climate%snapshot_p_anml%snapshot_baseline%Precip = 0._dp
    climate%snapshot_p_anml%T2m_anomaly              = 0._dp
    climate%snapshot_p_anml%Precip_anomaly           = 0._dp
    climate%snapshot_p_anml%T2m_anomaly_0            = 0._dp
    climate%snapshot_p_anml%Precip_anomaly_0         = 0._dp
    climate%snapshot_p_anml%T2m_anomaly_1            = 0._dp
    climate%snapshot_p_anml%Precip_anomaly_1         = 0._dp
    climate%snapshot_p_anml%T2m                      = 0._dp
    climate%snapshot_p_anml%Precip                   = 0._dp

    ! Run the chosen realistic climate model
    climate%snapshot_p_anml%snapshot_baseline%has_insolation = .FALSE.

    ! Read single-time data from external file
    ! Determine which climate model to initialise for this region
    select case( region_name)
    case ('NAM')
        filename_climate_snapshot                                      = C%climate_snp_p_anml_filename_snapshot_NAM
        climate%snapshot_p_anml%filename_climate_anomalies             = C%climate_snp_p_anml_filename_anomalies_NAM
        climate%snapshot_p_anml%snapshot_baseline%precip_CC_correction = C%precip_CC_correction_NAM
        climate%snapshot_p_anml%snapshot_baseline%do_lapse_rates       = C%do_lapse_rate_corrections_NAM
        climate%snapshot_p_anml%snapshot_baseline%lapse_rate_temp      = C%lapse_rate_temp_NAM
        climate%snapshot_p_anml%snapshot_baseline%has_insolation       = C%choice_SMB_model_NAM == 'IMAU-ITM'
    case ('EAS')
        filename_climate_snapshot                                      = C%climate_snp_p_anml_filename_snapshot_EAS
        climate%snapshot_p_anml%filename_climate_anomalies             = C%climate_snp_p_anml_filename_anomalies_EAS
        climate%snapshot_p_anml%snapshot_baseline%precip_CC_correction = C%precip_CC_correction_EAS
        climate%snapshot_p_anml%snapshot_baseline%do_lapse_rates       = C%do_lapse_rate_corrections_EAS
        climate%snapshot_p_anml%snapshot_baseline%lapse_rate_temp      = C%lapse_rate_temp_EAS
        climate%snapshot_p_anml%snapshot_baseline%has_insolation       = C%choice_SMB_model_EAS == 'IMAU-ITM'
    case ('GRL')
        filename_climate_snapshot                                      = C%climate_snp_p_anml_filename_snapshot_GRL
        climate%snapshot_p_anml%filename_climate_anomalies             = C%climate_snp_p_anml_filename_anomalies_GRL
        climate%snapshot_p_anml%snapshot_baseline%precip_CC_correction = C%precip_CC_correction_GRL
        climate%snapshot_p_anml%snapshot_baseline%do_lapse_rates       = C%do_lapse_rate_corrections_GRL
        climate%snapshot_p_anml%snapshot_baseline%lapse_rate_temp      = C%lapse_rate_temp_GRL
        climate%snapshot_p_anml%snapshot_baseline%has_insolation       = C%choice_SMB_model_GRL == 'IMAU-ITM'
    case ('ANT')
        filename_climate_snapshot                                      = C%climate_snp_p_anml_filename_snapshot_ANT
        climate%snapshot_p_anml%filename_climate_anomalies             = C%climate_snp_p_anml_filename_anomalies_ANT
        climate%snapshot_p_anml%snapshot_baseline%precip_CC_correction = C%precip_CC_correction_ANT
        climate%snapshot_p_anml%snapshot_baseline%do_lapse_rates       = C%do_lapse_rate_corrections_ANT
        climate%snapshot_p_anml%snapshot_baseline%lapse_rate_temp      = C%lapse_rate_temp_ANT
        climate%snapshot_p_anml%snapshot_baseline%has_insolation       = C%choice_SMB_model_ANT == 'IMAU-ITM'
    case default
      CALL crash('unknown region_name "' // region_name // '"')
    end select

    ! Read in baseline variables
    CALL read_field_from_file_2D(         filename_climate_snapshot, 'Hs'    , mesh, C%output_dir, climate%snapshot_p_anml%snapshot_baseline%Hs)
    CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'T2m'   , mesh, C%output_dir, climate%snapshot_p_anml%snapshot_baseline%T2m)
    CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'Precip', mesh, C%output_dir, climate%snapshot_p_anml%snapshot_baseline%Precip)

    ! Read in anomalies

    ! Initialise anomaly timeframes
    climate%snapshot_p_anml%anomaly_t0 = C%start_time_of_run - 200._dp
    climate%snapshot_p_anml%anomaly_t1 = C%start_time_of_run - 100._dp
    call update_climate_timeframes( mesh, climate, C%start_time_of_run)

    ! Interpolate between the two timeframes to find the applied anomaly
    call interpolate_between_climate_anomaly_timeframes(mesh, climate, C%start_time_of_run)

    call apply_geometry_downscaling_corrections( mesh, ice, climate, climate%snapshot_p_anml%snapshot_baseline, 0.0_dp)

    ! Initialises the insolation (if needed)
    IF (climate%snapshot_p_anml%snapshot_baseline%has_insolation) THEN
    IF (C%choice_insolation_forcing == 'none') THEN
        CALL crash('Chosen climate or SMB model cannot be used with choice_insolation_forcing = "none"!')
    ELSE
        CALL initialise_insolation_forcing( climate%snapshot_p_anml%snapshot_baseline, mesh)
        IF (C%start_time_of_run < 0._dp) THEN
            timeframe_init_insolation = C%start_time_of_run
        ELSE
            timeframe_init_insolation = 0._dp
        END IF
        CALL get_insolation_at_time( mesh, timeframe_init_insolation, climate%snapshot_p_anml%snapshot_baseline)
        climate%Q_TOA          = climate%snapshot_p_anml%snapshot_baseline%Q_TOA
        climate%snapshot%Q_TOA = climate%snapshot_p_anml%snapshot_baseline%Q_TOA
    END IF
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_snp_p_anml

  SUBROUTINE remap_climate_snp_p_anml(mesh_old, mesh_new, ice, climate, region_name, time)
  ! In/out variables
    type(type_mesh),                        intent(in)    :: mesh_old
    type(type_mesh),                        intent(in)    :: mesh_new
    type(type_ice_model),                   intent(in)    :: ice
    type(type_climate_model),               intent(inout) :: climate
    character(LEN=3),                       intent(in)    :: region_name
    real(dp),                               intent(in)    :: time

    ! Local variables
    character(LEN=256), parameter                         :: routine_name = 'remap_climate_snp_p_anml'
    character(LEN=256)                                    :: choice_climate_model
    character(LEN=256)                                    :: filename_climate_snapshot,filename_atm_dT
    character(LEN=256)                                    :: choice_SMB_model
    INTEGER                                               :: vi, m
    real(dp)                                              :: w0, w1

    ! Add routine to path
    call init_routine( routine_name)

    select case( region_name)
    case ('NAM')
      filename_climate_snapshot                                      = C%climate_snp_p_anml_filename_snapshot_NAM
      climate%snapshot_p_anml%filename_climate_anomalies             = C%climate_snp_p_anml_filename_anomalies_NAM
      climate%snapshot_p_anml%snapshot_baseline%precip_CC_correction = C%precip_CC_correction_NAM
      climate%snapshot_p_anml%snapshot_baseline%do_lapse_rates       = C%do_lapse_rate_corrections_NAM
      climate%snapshot_p_anml%snapshot_baseline%lapse_rate_temp      = C%lapse_rate_temp_NAM
      climate%snapshot_p_anml%snapshot_baseline%has_insolation       = C%choice_SMB_model_NAM == 'IMAU-ITM'
    case ('EAS')
      filename_climate_snapshot                                      = C%climate_snp_p_anml_filename_snapshot_EAS
      climate%snapshot_p_anml%filename_climate_anomalies             = C%climate_snp_p_anml_filename_anomalies_EAS
      climate%snapshot_p_anml%snapshot_baseline%precip_CC_correction = C%precip_CC_correction_EAS
      climate%snapshot_p_anml%snapshot_baseline%do_lapse_rates       = C%do_lapse_rate_corrections_EAS
      climate%snapshot_p_anml%snapshot_baseline%lapse_rate_temp      = C%lapse_rate_temp_EAS
      climate%snapshot_p_anml%snapshot_baseline%has_insolation       = C%choice_SMB_model_EAS == 'IMAU-ITM'
    case ('GRL')
      filename_climate_snapshot                                      = C%climate_snp_p_anml_filename_snapshot_GRL
      climate%snapshot_p_anml%filename_climate_anomalies             = C%climate_snp_p_anml_filename_anomalies_GRL
      climate%snapshot_p_anml%snapshot_baseline%precip_CC_correction = C%precip_CC_correction_GRL
      climate%snapshot_p_anml%snapshot_baseline%do_lapse_rates       = C%do_lapse_rate_corrections_GRL
      climate%snapshot_p_anml%snapshot_baseline%lapse_rate_temp      = C%lapse_rate_temp_GRL
      climate%snapshot_p_anml%snapshot_baseline%has_insolation       = C%choice_SMB_model_GRL == 'IMAU-ITM'
    case ('ANT')
      filename_climate_snapshot                                      = C%climate_snp_p_anml_filename_snapshot_ANT
      climate%snapshot_p_anml%filename_climate_anomalies             = C%climate_snp_p_anml_filename_anomalies_ANT
      climate%snapshot_p_anml%snapshot_baseline%precip_CC_correction = C%precip_CC_correction_ANT
      climate%snapshot_p_anml%snapshot_baseline%do_lapse_rates       = C%do_lapse_rate_corrections_ANT
      climate%snapshot_p_anml%snapshot_baseline%lapse_rate_temp      = C%lapse_rate_temp_ANT
      climate%snapshot_p_anml%snapshot_baseline%has_insolation       = C%choice_SMB_model_ANT == 'IMAU-ITM'
    case default
      call crash('unknown region_name "' // region_name // '"')
    end select

    call reallocate_bounds( climate%snapshot_p_anml%snapshot_baseline%Hs, mesh_new%vi1, mesh_new%vi2)

    CALL read_field_from_file_2D(         filename_climate_snapshot, 'Hs'    , mesh_new, C%output_dir, climate%snapshot_p_anml%snapshot_baseline%Hs)
    CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'T2m'   , mesh_new, C%output_dir, climate%snapshot_p_anml%snapshot_baseline%T2m)
    CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'Precip', mesh_new, C%output_dir, climate%snapshot_p_anml%snapshot_baseline%Precip)

    
    ! Update anomaly timeframes
    call update_climate_timeframes( mesh_new, climate, time)

    ! Interpolate between the two timeframes to find the applied anomaly
    call interpolate_between_climate_anomaly_timeframes(mesh_new, climate, time)
    
    call apply_geometry_downscaling_corrections( mesh_new, ice, climate, climate%snapshot_p_anml%snapshot_baseline, 0.0_dp)
    

    IF (climate%snapshot_p_anml%snapshot_baseline%has_insolation .eqv. .TRUE.) THEN
      call remap_insolation( climate%snapshot_p_anml%snapshot_baseline, mesh_new)
    END IF

    ! Finalise routine path
    call finalise_routine( routine_name)

  END SUBROUTINE remap_climate_snp_p_anml

  subroutine update_climate_timeframes(mesh, climate, time)

    ! In/output variables:
    type(type_mesh),                  intent(in   ) :: mesh
    type(type_climate_model),         intent(inout) :: climate
    real(dp),                         intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'update_climate_timeframes'
    integer                             :: ncid, id_dim_time, nt, id_var_time, ierr
    real(dp), dimension(:), allocatable :: time_from_file
    integer                             :: ti0, ti1

    ! Add routine to path
    call init_routine( routine_name)

    ! Read time variable from the file
    call open_existing_netcdf_file_for_reading( climate%snapshot_p_anml%filename_climate_anomalies, ncid)
    call check_time( climate%snapshot_p_anml%filename_climate_anomalies, ncid)
    call inquire_dim_multopt( climate%snapshot_p_anml%filename_climate_anomalies, ncid, field_name_options_time, id_dim_time, dim_length = nt)
    call inquire_var_multopt( climate%snapshot_p_anml%filename_climate_anomalies, ncid, field_name_options_time, id_var_time)
    allocate( time_from_file( nt))
    call read_var_primary( climate%snapshot_p_anml%filename_climate_anomalies, ncid, id_var_time, time_from_file)
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

    climate%snapshot_p_anml%anomaly_t0 = time_from_file( ti0)
    climate%snapshot_p_anml%anomaly_t1 = time_from_file( ti1)

    ! Read the two timeframes
    call read_field_from_file_2D_monthly( climate%snapshot_p_anml%filename_climate_anomalies, 'T2m_anomaly', &
      mesh, C%output_dir, climate%snapshot_p_anml%T2m_anomaly_0, &
      time_to_read = climate%snapshot_p_anml%anomaly_t0)
    call read_field_from_file_2D_monthly( climate%snapshot_p_anml%filename_climate_anomalies, 'T2m_anomaly', &
      mesh, C%output_dir, climate%snapshot_p_anml%T2m_anomaly_1, &
      time_to_read = climate%snapshot_p_anml%anomaly_t1)
    call read_field_from_file_2D_monthly( climate%snapshot_p_anml%filename_climate_anomalies, 'Precip_anomaly', &
      mesh, C%output_dir, climate%snapshot_p_anml%Precip_anomaly_0, &
      time_to_read = climate%snapshot_p_anml%anomaly_t0)
    call read_field_from_file_2D_monthly( climate%snapshot_p_anml%filename_climate_anomalies, 'Precip_anomaly', &
      mesh, C%output_dir, climate%snapshot_p_anml%Precip_anomaly_1, &
      time_to_read = climate%snapshot_p_anml%anomaly_t1)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_climate_timeframes

  subroutine interpolate_between_climate_anomaly_timeframes(mesh, climate, time)
  ! Calculates anomalies by interpolating between the two 
  ! timeframes to find the applied anomaly, then applies it to the climate variables


    ! In/output variables:
    type(type_mesh),                  intent(in   ) :: mesh
    type(type_climate_model),         intent(inout) :: climate
    real(dp),                         intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'interpolate_between_climate_anomaly_timeframes'
    real(dp)                            :: w0, w1
    integer                             :: vi, m

    ! Add routine to path
    call init_routine( routine_name)

    w0 = (climate%snapshot_p_anml%anomaly_t1 - time) / (climate%snapshot_p_anml%anomaly_t1 - climate%snapshot_p_anml%anomaly_t0)
    w1 = 1._dp - w0

    do m = 1,12
    do vi = mesh%vi1, mesh%vi2

      ! Note that the baseline and the applied temperature and precip are monthly, but the anomaly is annual
      climate%snapshot_p_anml%T2m_anomaly( vi,m)    = w0 * climate%snapshot_p_anml%T2m_anomaly_0( vi,m)    + w1 * climate%snapshot_p_anml%T2m_anomaly_1( vi,m)
      climate%snapshot_p_anml%Precip_anomaly( vi,m) = w0 * climate%snapshot_p_anml%Precip_anomaly_0( vi,m) + w1 * climate%snapshot_p_anml%Precip_anomaly_1( vi,m)

      ! Add anomaly to snapshots to find the applied temperature and precipitation
      climate%snapshot_p_anml%T2m( vi,m)    = climate%snapshot_p_anml%snapshot_baseline%T2m( vi,m)    + climate%snapshot_p_anml%T2m_anomaly( vi,m)
      climate%snapshot_p_anml%Precip( vi,m) = climate%snapshot_p_anml%snapshot_baseline%Precip( vi,m) + climate%snapshot_p_anml%Precip_anomaly( vi,m)

      ! Copy to climate model
      climate%T2m( vi,m)    = climate%snapshot_p_anml%T2m( vi,m)
      climate%Precip( vi,m) = climate%snapshot_p_anml%Precip( vi,m)

    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine interpolate_between_climate_anomaly_timeframes


END MODULE climate_snapshot_plus_anomalies
