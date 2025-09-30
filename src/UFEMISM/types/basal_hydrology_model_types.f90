MODULE basal_hydrology_model_types

  ! The different data types used in the basal hydrology modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE laddie_model_types                                     , ONLY: type_laddie_model
  use laddie_forcing_types, only: type_laddie_forcing
  use reference_geometry_types, only: type_reference_geometry

  IMPLICIT NONE

  private

  public :: type_basal_hydrology_model, type_basal_hydrology_model_inverted

! ===== Types =====
! =================

  TYPE type_basal_hydrology_model
    ! The basal hydrology model data structure.

    ! Main data fields
    ! Fill here what we need for the basal hydrology model
    real(dp), dimension(:), allocatable :: W                   ! Basal water depth
    real(dp), dimension(:), allocatable :: W_max               ! Maximum basal water depth
    real(dp), dimension(:), allocatable :: W_min               ! Minimum basal water depth

    real(dp), dimension(:), allocatable :: W_til               ! Basal water depth in till
    real(dp), dimension(:), allocatable :: W_max_til           ! Maximum basal water depth in till
    real(dp), dimension(:), allocatable :: W_min_til           ! Minimum basal water depth in till

    real(dp), dimension(:), allocatable :: P                   ! Pressure of the ice on the basal water
    real(dp), dimension(:), allocatable :: P_o                 ! Overburden pressure of the ice on the basal water
    real(dp), dimension(:), allocatable :: P_min               ! Minimum pressure of the ice on the basal water
    real(dp), dimension(:), allocatable :: N                   ! Effective pressure (P_o - P)

    real(dp), dimension(:), allocatable :: m                   ! Total input of water to the basal system

    real(dp), dimension(:), allocatable :: rho_w               ! Density of water
    real(dp), dimension(:), allocatable :: rho_i               ! Density of ice

    real(dp), dimension(:), allocatable :: D                   ! Diffusivity

    real(dp), dimension(:), allocatable :: K                   ! Effective conductivity
    real(dp), dimension(:), allocatable :: k_coef              ! Coefficient of effective conductivity
    real(dp), dimension(:), allocatable :: alpha               ! Exponent used in effective conductivity
    real(dp), dimension(:), allocatable :: beta                ! Exponent used in effective conductivity

    real(dp), dimension(:), allocatable :: u                   ! Velocity in x direction
    real(dp), dimension(:), allocatable :: v                   ! Velocity in y direction

    real(dp), dimension(:), allocatable :: Cd                  ! Gradual drain of water in till

    real(dp), dimension(:), allocatable :: phi_0               ! Englacial porosity

    real(dp), dimension(:), allocatable :: g                   ! Gravitational acceleration

    real(dp), dimension(:), allocatable :: u_b                 ! Ice sliding velocity in x direction
    real(dp), dimension(:), allocatable :: v_b                 ! Ice sliding velocity in y direction

    real(dp), dimension(:), allocatable :: c_1                 ! Scaling coefficient 1 (non-negative opening term)
    real(dp), dimension(:), allocatable :: c_2                 ! Scaling coefficient 2 (closing term)

    real(dp), dimension(:), allocatable :: A                   ! Softness of the ice

    real(dp), dimension(:), allocatable :: W_r                 ! Maximum roughness scale of basal topography

    real(dp), dimension(:), allocatable :: psi                 ! Hydraulic potential

    real(dp), dimension(:), allocatable :: b                   ! Bedrock elevation

    real(dp), dimension(:), allocatable :: H_i                 ! Ice thickness

    logical, dimension(:), allocatable  :: floating            ! Is the ice floating?
    logical, dimension(:), allocatable  :: ice_free            ! Is there ice at this location?

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp)                                :: t_next
    !real(dp)                                :: dt                  ! Timestep
    

    ! LADDIE
    type(type_laddie_model)                       :: laddie
    type(type_laddie_forcing)                     :: forcing

    type(type_BMB_model_inverted) :: inv

  END TYPE type_basal_hydrology_model

CONTAINS

END MODULE basal_hydrology_model_types
