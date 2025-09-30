module parameters

  ! Some mathematical and physical constants

  use precisions, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_signaling_nan

  implicit none

! ===== Global variables =====
! ============================
  real(dp)            :: pi
  real(dp)            :: NaN

  real(dp), parameter :: sec_per_year                     = 31556943.36_dp            ! = 365.2424 * 24 * 3600
  real(dp), parameter :: sec_per_day                      = 86400.0_dp                ! = 24 * 3600
  real(dp), parameter :: T0                               = 273.16_dp                 ! [K]                 Triple point of water
  real(dp), parameter :: Clausius_Clapeyron_gradient      = 8.7E-04_dp                ! [K m^-1]            Clausius Clapeyron gradient
  real(dp), parameter :: grav                             = 9.81_dp                   ! [m s^-2]            Acceleration of gravity
  real(dp), parameter :: earth_radius                     = 6.371221E6_dp             ! [m] Earth           Radius
  real(dp), parameter :: L_fusion                         = 3.335E+5_dp               ! [J kg-1]            Latent heat of fusion
  real(dp), parameter :: ice_density                      =  917.0_dp                 ! [kg m^-3]           Ice density
  real(dp), parameter :: freshwater_density               = 1000.0_dp                 ! [kg m^-3]           Freshwater density
  real(dp), parameter :: seawater_density                 = 1027.0_dp                 ! [kg m^-3]           Seawater density
  real(dp), parameter :: earth_density                    = 5511.57_dp                ! [kg m^-3]           Total mean Earth density
  real(dp), parameter :: R_gas                            = 8.314_dp                  ! [J mol^-1 K^-1]     Gas constant
  real(dp), parameter :: cp_ocean                         = 3.974E3_dp                ! [J kg^-1 K^-1]      Specific heat capacity of ocean water
  real(dp), parameter :: ocean_area                       = 3.611E14_dp               ! [m^2]               World ocean area
  real(dp), parameter :: earth_rotation_rate              = 7.2921E-5_dp              ! [s^-1]              Earth's rotation rate

! ===== LADDIE parameters ====
! ============================

  real(dp), parameter :: freezing_lambda_1                = -5.73E-2_dp               ! [K PSU^-1]          Freezing point salinity coefficient
  real(dp), parameter :: freezing_lambda_2                =  8.32E-2_dp               ! [K]                 Freezing point offset
  real(dp), parameter :: freezing_lambda_3                =  7.61E-4_dp               ! [K m^-1]            Freezing point depth coefficient
  real(dp), parameter :: cp_ice                           = 2009.0_dp                 ! [J kg^-1 K^-1]      Specific heat capacity of ice
  real(dp), parameter :: Stanton_number                   = 5.9E-4_dp                 ! []                  Effective thermal Stanton number
  real(dp), parameter :: Prandtl_number                   = 13.8_dp                   ! []                  Prandtl number
  real(dp), parameter :: Schmidt_number                   = 2432.0_dp                 ! []                  Schmidt number
  real(dp), parameter :: molecular_viscosity              = 1.95E-6_dp                ! [m^2 s^-1]          Molecular viscosity

contains

  subroutine initialise_constants
    call initialise_pi
    call initialise_NaN
  end subroutine initialise_constants

  subroutine initialise_pi
    ! NOTE: yes, it seems idiotic, but Fortran does not have pi as an intrinsic constant/function.
    ! The definition below yields the best numerical accuracy (as both 4 and 1 can be represented
    ! exactly in binary)
    ! See also: https://stackoverflow.com/questions/2157920/why-define-pi-4atan1-d0
    pi = 4._dp * atan( 1._dp)
  end subroutine initialise_pi

  subroutine initialise_NaN
    ! Strangely, Fortran does have an isnan intrinsic, but not nan itself...
    NaN = ieee_value( NaN, ieee_signaling_nan)
    ! Check if it worked
    if (.not. isnan( NaN)) error stop 'Could not initialise NaN'
  end subroutine initialise_NaN

end module parameters
