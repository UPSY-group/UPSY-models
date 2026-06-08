module ISMIP7_climate_model_type

  use precisions, only: dp
  use ISMIP7_forcing_field_types, only: type_climate_field_ISMIP7_monthly, type_climate_field_ISMIP7_yearly

  implicit none

  private

  public :: type_climate_model_ISMIP7

  type type_climate_model_ISMIP7

    ! Baseline TS and surface elevation
    real(dp), dimension(:,:), allocatable :: T2m_baseline    ! [K]                     Baseline monthly T2m
    real(dp), dimension(:,:), allocatable :: Precip_baseline ! [m.w.e.]                Baseline monthly Precip
    real(dp), dimension(:)  , allocatable :: Hs_baseline     ! [m w.r.t. PD sea level] Baseline surface elevation

    ! Fields
    type(type_climate_field_ISMIP7_monthly) :: tas
    type(type_climate_field_ISMIP7_monthly) :: tas_anomaly
    type(type_climate_field_ISMIP7_monthly) :: pr
    type(type_climate_field_ISMIP7_monthly) :: pr_anomaly
    type(type_climate_field_ISMIP7_yearly)  :: dtsdz

  end type type_climate_model_ISMIP7

contains

end module ISMIP7_climate_model_type
