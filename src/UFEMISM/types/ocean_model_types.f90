MODULE ocean_model_types

  ! The different data types used in the ocean modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  use grid_types, only: type_grid
  use reference_geometry_types, only: type_reference_geometry

  IMPLICIT NONE

! ===== Types =====
! =================
  TYPE type_ocean_model_transient_deltaT
    ! Main data fields to compute the ocean model T and S transiently
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T0                           ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: S0                           ! [PSU]             Salinity

    ! deltaT record
    REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: dT_series_time
    REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: dT_series
    REAL(dp),                   ALLOCATABLE     :: dT_t0, dT_t1, dT_at_t0, dT_at_t1

  END TYPE type_ocean_model_transient_deltaT

  TYPE type_ocean_model_deltaT
    ! Main data fields to compute the ocean model T and S transiently
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T0                           ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: S0                           ! [PSU]             Salinity

    ! deltaT value
    REAL(dp),                   ALLOCATABLE     :: dT

  END TYPE type_ocean_model_deltaT

  TYPE type_ocean_model_GlacialIndex
    ! Main data fields to compute the ocean model T and S transiently
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T0_warm                           ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: S0_warm                           ! [PSU]             Salinity
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T0_cold                           ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: S0_cold                           ! [PSU]             Salinity

    ! Glacial Index record
    REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: GI_series_time
    REAL(dp), DIMENSION(:    ), ALLOCATABLE     :: GI_series
    REAL(dp),                   ALLOCATABLE     :: GI_t0, GI_t1, GI_at_t0, GI_at_t1


  END TYPE type_ocean_model_GlacialIndex

  type type_ocean_model_snapshot_nudge2D

    type(type_reference_geometry)           :: target_geometry            !< The geometry that the BMB inversion should aim to reproduce
    logical,  dimension(:    ), allocatable :: target_mask_shelf          !< Shelf mask of the target geometry

    type(type_grid)                         :: grid_ref
    integer                                 :: ndepth_ref
    real(dp), dimension(:    ), allocatable :: depth_ref
    real(dp), dimension(:,:  ), allocatable :: T_ref_grid, S_ref_grid     !< Reference 3-D ocean snapshot on its original grid (in vectorised form)
    real(dp), dimension(:,:  ), allocatable :: T_ref, S_ref               !< Reference 3-D ocean snapshot on the model mesh
    real(dp), dimension(:    ), allocatable :: deltaT_nudge               !< 2-D temperature nudging term on the model mesh
    real(dp), dimension(:,:  ), allocatable :: T, S                       !< Applied 3-D ocean snapshot on the model mesh
    real(dp), dimension(:    ), allocatable :: deltaT_nudge_grid          !< 2-D temperature nudging term on the original grid (in vectorised form)
    real(dp), dimension(:,:  ), allocatable :: T_grid, S_grid             !< Applied 3-D ocean snapshot on the original grid (in vectorised form)
    character(len=1024)                     :: output_filename            !< Filename for output file with nudged snapshot on the original grid

  end type type_ocean_model_snapshot_nudge2D

  type type_ocean_model_snapshot_plus_anomalies

    ! Baseline ocean
    real(dp), dimension(:,:), allocatable :: T_baseline
    real(dp), dimension(:,:), allocatable :: S_baseline

    ! Two anomaly timeframes enveloping the current model time
    real(dp)                              :: anomaly_t0
    real(dp), dimension(:,:), allocatable :: T_anomaly_0
    real(dp), dimension(:,:), allocatable :: S_anomaly_0

    real(dp)                              :: anomaly_t1
    real(dp), dimension(:,:), allocatable :: T_anomaly_1
    real(dp), dimension(:,:), allocatable :: S_anomaly_1

    ! Time-weighted anomaly
    real(dp), dimension(:,:), allocatable :: T_anomaly
    real(dp), dimension(:,:), allocatable :: S_anomaly

    ! Applied ocean
    real(dp), dimension(:,:), allocatable :: T    ! = baseline + anomaly
    real(dp), dimension(:,:), allocatable :: S

  end type type_ocean_model_snapshot_plus_anomalies

  type type_ocean_model_snapshot_nudge2D_plus_anomalies

    type(type_reference_geometry)           :: target_geometry            !< The geometry that the BMB inversion should aim to reproduce
    logical,  dimension(:    ), allocatable :: target_mask_shelf          !< Shelf mask of the target geometry

    type(type_grid)                         :: grid_ref
    integer                                 :: ndepth_ref
    real(dp), dimension(:    ), allocatable :: depth_ref
    real(dp), dimension(:,:  ), allocatable :: T_ref_grid, S_ref_grid     !< Reference 3-D ocean snapshot on its original grid (in vectorised form)
    real(dp), dimension(:,:  ), allocatable :: T_ref, S_ref               !< Reference 3-D ocean snapshot on the model mesh
    real(dp), dimension(:    ), allocatable :: deltaT_nudge               !< 2-D temperature nudging term on the model mesh
    real(dp), dimension(:,:  ), allocatable :: T_baseline, S_baseline     !< Applied 3-D ocean snapshot on the model mesh
    real(dp), dimension(:    ), allocatable :: deltaT_nudge_grid          !< 2-D temperature nudging term on the original grid (in vectorised form)
    real(dp), dimension(:,:  ), allocatable :: T_grid, S_grid             !< Applied 3-D ocean snapshot on the original grid (in vectorised form)

    ! Two anomaly timeframes enveloping the current model time
    real(dp)                              :: anomaly_t0
    real(dp), dimension(:,:), allocatable :: T_anomaly_0
    real(dp), dimension(:,:), allocatable :: S_anomaly_0

    real(dp)                              :: anomaly_t1
    real(dp), dimension(:,:), allocatable :: T_anomaly_1
    real(dp), dimension(:,:), allocatable :: S_anomaly_1

    ! Time-weighted anomaly
    real(dp), dimension(:,:), allocatable :: T_anomaly
    real(dp), dimension(:,:), allocatable :: S_anomaly

    ! Applied ocean
    real(dp), dimension(:,:), allocatable :: T    ! = baseline + anomaly
    real(dp), dimension(:,:), allocatable :: S

  end type type_ocean_model_snapshot_nudge2D_plus_anomalies

  TYPE type_ocean_model
    ! The ocean model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: T                           ! [degrees Celsius] Temperature
    REAL(dp), DIMENSION(:,:  ), ALLOCATABLE :: S                           ! [PSU]             Salinity

    ! Secondary data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T_draft                     ! [degrees Celsius] Temperature at ice base
    REAL(dp), DIMENSION(:    ), ALLOCATABLE :: T_freezing_point            ! [degrees Celsius] Pressure freezing point of water

    ! Sub-models
    TYPE(type_ocean_model_transient_deltaT)        :: deltaT_transient
    TYPE(type_ocean_model_GlacialIndex)            :: GI
    TYPE(type_ocean_model_deltaT)                  :: deltaT
    type(type_ocean_model_snapshot_nudge2D)        :: snapshot_nudge2D
    type(type_ocean_model_snapshot_plus_anomalies) :: snapshot_plus_anomalies
    type(type_ocean_model_snapshot_nudge2D_plus_anomalies) :: snapshot_nudge2D_plus_anomalies

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next

  END TYPE type_ocean_model

CONTAINS

END MODULE ocean_model_types