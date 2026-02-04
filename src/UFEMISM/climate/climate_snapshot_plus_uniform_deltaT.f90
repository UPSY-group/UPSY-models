MODULE climate_snapshot_plus_uniform_deltaT

  ! Snapshot + dT climate model

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
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
  use climate_model_utilities                                , only: allocate_climate_snapshot, read_climate_snapshot, get_insolation_at_time, apply_precipitation_CC_correction, apply_geometry_downscaling_corrections

  IMPLICIT NONE

  private

  public :: run_climate_model_snapshot_plus_uniform_deltaT
  public :: initialise_climate_model_snapshot_plus_uniform_deltaT
  public :: remap_climate_snapshot_plus_uniform_deltaT


CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_climate_model_snapshot_plus_uniform_deltaT( mesh, ice, climate, time)
    ! Calculate the climate
    !
    ! Use a snapshot plus a prescribed uniform deltaT

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_climate_model_snapshot_plus_uniform_deltaT'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen snapshot+deltaT climate model

    ! Update temperature and precipitation fields based on the mismatch between
    ! the ice sheet surface elevation in the forcing climate and the model's ice sheet surface elevation
    call apply_geometry_downscaling_corrections( mesh, ice, climate, climate%snapshot_unif_dT%snapshot, climate%snapshot_unif_dT%deltaT)

    ! if needed for IMAU-ITM or climate matrix, we need to update insolation
    IF (climate%snapshot%has_insolation) THEN
      CALL get_insolation_at_time( mesh, time, climate%snapshot_unif_dT%snapshot)
      climate%Q_TOA = climate%snapshot%Q_TOA
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_snapshot_plus_uniform_deltaT

  SUBROUTINE initialise_climate_model_snapshot_plus_uniform_deltaT( mesh, ice, climate, region_name)
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
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_climate_model_snapshot_plus_uniform_deltaT'
    INTEGER                                               :: vi, m
    CHARACTER(LEN=256)                                    :: filename_climate_snapshot
    LOGICAL                                               :: do_lapse_rates
    REAL(dp)                                              :: timeframe_init_insolation

    ! Add routine to path
    CALL init_routine( routine_name)

    ALLOCATE( climate%snapshot_unif_dT%snapshot%Hs(     mesh%vi1:mesh%vi2))
    ALLOCATE( climate%snapshot_unif_dT%snapshot%T2m(    mesh%vi1:mesh%vi2,12))
    ALLOCATE( climate%snapshot_unif_dT%snapshot%Precip( mesh%vi1:mesh%vi2,12))

    ! Run the chosen realistic climate model
    climate%snapshot_unif_dT%snapshot%has_insolation = .FALSE.

    ! Read single-time data from external file
    ! Determine which climate model to initialise for this region
    IF     (region_name == 'NAM') THEN
        filename_climate_snapshot                              = C%filename_climate_snapshot_unif_dT_NAM
        climate%snapshot_unif_dT%snapshot%precip_CC_correction = C%precip_CC_correction_NAM
        climate%snapshot_unif_dT%snapshot%do_lapse_rates       = C%do_lapse_rate_corrections_NAM
        climate%snapshot_unif_dT%snapshot%lapse_rate_temp      = C%lapse_rate_temp_NAM
        climate%snapshot_unif_dT%deltaT                        = C%uniform_deltaT_NAM
        climate%snapshot_unif_dT%snapshot%has_insolation       = C%choice_SMB_model_NAM == 'IMAU-ITM'
    ELSEIF (region_name == 'EAS') THEN
        filename_climate_snapshot                              = C%filename_climate_snapshot_unif_dT_EAS
        climate%snapshot_unif_dT%snapshot%precip_CC_correction = C%precip_CC_correction_EAS
        climate%snapshot_unif_dT%snapshot%do_lapse_rates       = C%do_lapse_rate_corrections_EAS
        climate%snapshot_unif_dT%snapshot%lapse_rate_temp      = C%lapse_rate_temp_EAS
        climate%snapshot_unif_dT%deltaT                        = C%uniform_deltaT_EAS
        climate%snapshot_unif_dT%snapshot%has_insolation       = C%choice_SMB_model_EAS == 'IMAU-ITM'
    ELSEIF (region_name == 'GRL') THEN
        filename_climate_snapshot                              = C%filename_climate_snapshot_unif_dT_GRL
        climate%snapshot_unif_dT%snapshot%precip_CC_correction = C%precip_CC_correction_GRL
        climate%snapshot_unif_dT%snapshot%do_lapse_rates       = C%do_lapse_rate_corrections_GRL
        climate%snapshot_unif_dT%snapshot%lapse_rate_temp      = C%lapse_rate_temp_GRL
        climate%snapshot_unif_dT%deltaT                        = C%uniform_deltaT_GRL
        climate%snapshot_unif_dT%snapshot%has_insolation       = C%choice_SMB_model_GRL == 'IMAU-ITM'
    ELSEIF (region_name == 'ANT') THEN
        filename_climate_snapshot                              = C%filename_climate_snapshot_unif_dT_ANT
        climate%snapshot_unif_dT%snapshot%precip_CC_correction = C%precip_CC_correction_ANT
        climate%snapshot_unif_dT%snapshot%do_lapse_rates       = C%do_lapse_rate_corrections_ANT
        climate%snapshot_unif_dT%snapshot%lapse_rate_temp      = C%lapse_rate_temp_ANT
        climate%snapshot_unif_dT%deltaT                        = C%uniform_deltaT_ANT
        climate%snapshot_unif_dT%snapshot%has_insolation       = C%choice_SMB_model_ANT == 'IMAU-ITM'
    ELSE
      CALL crash('unknown region_name "' // region_name // '"')
    END IF

    CALL read_field_from_file_2D(         filename_climate_snapshot, 'Hs'    , mesh, C%output_dir, climate%snapshot_unif_dT%snapshot%Hs)
    CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'T2m'   , mesh, C%output_dir, climate%T2m)
    CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'Precip', mesh, C%output_dir, climate%Precip)
    climate%snapshot_unif_dT%snapshot%T2m    = climate%T2m
    climate%snapshot_unif_dT%snapshot%Precip = climate%Precip

    ! Adding deltaT to the temperature field (uniform in space and time)
    do vi = mesh%vi1, mesh%vi2
    do m = 1, 12
        climate%T2m( vi, m)   = climate%T2m( vi, m) + climate%snapshot_unif_dT%deltaT
    end do
    end do

    ! apply corrections (increase in Precip due to deltaT, plus downscaling correction)
    call apply_precipitation_CC_correction(mesh, climate, climate%snapshot_unif_dT%snapshot%precip_CC_correction, climate%snapshot_unif_dT%deltaT)
    call apply_geometry_downscaling_corrections( mesh, ice, climate, climate%snapshot_unif_dT%snapshot, climate%snapshot_unif_dT%deltaT)

    ! Initialises the insolation (if needed)
    IF (climate%snapshot_unif_dT%snapshot%has_insolation) THEN
    IF (C%choice_insolation_forcing == 'none') THEN
        CALL crash('Chosen climate or SMB model cannot be used with choice_insolation_forcing = "none"!')
    ELSE
        CALL initialise_insolation_forcing( climate%snapshot_unif_dT%snapshot, mesh)
        IF (C%start_time_of_run < 0._dp) THEN
            timeframe_init_insolation = C%start_time_of_run
        ELSE
            timeframe_init_insolation = 0._dp
        END IF
        CALL get_insolation_at_time( mesh, timeframe_init_insolation, climate%snapshot_unif_dT%snapshot)
        climate%Q_TOA          = climate%snapshot_unif_dT%snapshot%Q_TOA
        climate%snapshot%Q_TOA = climate%snapshot_unif_dT%snapshot%Q_TOA
    END IF
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_snapshot_plus_uniform_deltaT

  SUBROUTINE remap_climate_snapshot_plus_uniform_deltaT(mesh_old, mesh_new, ice, climate, region_name)
  ! In/out variables
    type(type_mesh),                        intent(in)    :: mesh_old
    type(type_mesh),                        intent(in)    :: mesh_new
    type(type_ice_model),                   intent(in)    :: ice
    type(type_climate_model),               intent(inout) :: climate
    character(LEN=3),                       intent(in)    :: region_name

    ! Local variables
    character(LEN=256), parameter                         :: routine_name = 'remap_climate_snapshot_plus_uniform_deltaT'
    character(LEN=256)                                    :: choice_climate_model
    character(LEN=256)                                    :: filename_climate_snapshot
    character(LEN=256)                                    :: choice_SMB_model
    INTEGER                                               :: vi, m

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine which climate model to initialise for this region
    select case( region_name)
    case ('NAM')
      filename_climate_snapshot                              = C%filename_climate_snapshot_unif_dT_NAM
      climate%snapshot_unif_dT%snapshot%precip_CC_correction = C%precip_CC_correction_NAM
      climate%snapshot_unif_dT%snapshot%lapse_rate_temp      = C%lapse_rate_temp_NAM
      climate%snapshot_unif_dT%deltaT                        = C%uniform_deltaT_NAM
      climate%snapshot_unif_dT%snapshot%has_insolation       = C%choice_SMB_model_NAM == 'IMAU-ITM'
    case ('EAS')
      filename_climate_snapshot                              = C%filename_climate_snapshot_unif_dT_EAS
      climate%snapshot_unif_dT%snapshot%precip_CC_correction = C%precip_CC_correction_EAS
      climate%snapshot_unif_dT%snapshot%lapse_rate_temp      = C%lapse_rate_temp_EAS
      climate%snapshot_unif_dT%deltaT                        = C%uniform_deltaT_EAS
      climate%snapshot_unif_dT%snapshot%has_insolation       = C%choice_SMB_model_EAS == 'IMAU-ITM'
    case ('GRL')
      filename_climate_snapshot                              = C%filename_climate_snapshot_unif_dT_GRL
      climate%snapshot_unif_dT%snapshot%precip_CC_correction = C%precip_CC_correction_GRL
      climate%snapshot_unif_dT%snapshot%lapse_rate_temp      = C%lapse_rate_temp_GRL
      climate%snapshot_unif_dT%deltaT                        = C%uniform_deltaT_GRL
      climate%snapshot_unif_dT%snapshot%has_insolation       = C%choice_SMB_model_GRL == 'IMAU-ITM'
    case ('ANT')
      filename_climate_snapshot                              = C%filename_climate_snapshot_unif_dT_ANT
      climate%snapshot_unif_dT%snapshot%precip_CC_correction = C%precip_CC_correction_ANT
      climate%snapshot_unif_dT%snapshot%lapse_rate_temp      = C%lapse_rate_temp_ANT
      climate%snapshot_unif_dT%deltaT                        = C%uniform_deltaT_ANT
      climate%snapshot_unif_dT%snapshot%has_insolation       = C%choice_SMB_model_ANT == 'IMAU-ITM'
    case default
      call crash('unknown region_name "' // region_name // '"')
    end select

    call reallocate_bounds( climate%snapshot%Hs, mesh_new%vi1, mesh_new%vi2)

    IF (climate%snapshot_unif_dT%snapshot%has_insolation .eqv. .TRUE.) THEN
      call remap_insolation( climate%snapshot, mesh_new)
    END IF

    ! Read single-time data from external file
    call read_field_from_file_2D( filename_climate_snapshot, 'Hs', mesh_new, C%output_dir, climate%snapshot%Hs)
    call read_field_from_file_2D_monthly( filename_climate_snapshot, 'T2m', mesh_new, C%output_dir, climate%T2m)
    call read_field_from_file_2D_monthly( filename_climate_snapshot, 'Precip', mesh_new, C%output_dir, climate%Precip)

    ! Adding deltaT to the temperature field (uniform in space and time)
    do vi = mesh_new%vi1, mesh_new%vi2
    do m = 1, 12
        climate%T2m( vi, m) = climate%T2m( vi, m) + climate%snapshot_unif_dT%deltaT
    end do
    end do

    ! apply corrections (increase in Precip due to deltaT, plus downscaling correction)
    call apply_precipitation_CC_correction(mesh_new, climate, climate%snapshot_unif_dT%snapshot%precip_CC_correction, climate%snapshot_unif_dT%deltaT)
    call apply_geometry_downscaling_corrections( mesh_new, ice, climate, climate%snapshot_unif_dT%snapshot, climate%snapshot_unif_dT%deltaT)

    ! Finalise routine path
    call finalise_routine( routine_name)

  END SUBROUTINE remap_climate_snapshot_plus_uniform_deltaT


END MODULE climate_snapshot_plus_uniform_deltaT
