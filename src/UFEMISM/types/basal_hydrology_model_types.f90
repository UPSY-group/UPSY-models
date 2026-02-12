MODULE basal_hydrology_model_types

  ! The different data types used in the basal hydrology modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE laddie_model_types                                     , ONLY: type_laddie_model
  use laddie_forcing_types, only: type_laddie_forcing
  use reference_geometry_types, only: type_reference_geometry
  USE mesh_types                                             , ONLY: type_mesh
  USE CSR_sparse_matrix_type                                 , ONLY: type_sparse_matrix_CSR_dp

  IMPLICIT NONE

  private

  public :: type_basal_hydrology_model

! ===== Types =====
! =================

  TYPE type_basal_hydrology_model
    ! The basal hydrology model data structure.

    ! Main data fields
    ! Fill here what we need for the basal hydrology model
    ! Part of this might need to go in ice_model_types?
    real(dp), dimension(:), allocatable :: W                   ! Basal water depth
    real(dp), dimension(:), allocatable :: dW_dx_b             ! Derivative of W to x on B grid.
    real(dp), dimension(:), allocatable :: W_b                 ! Basal water depth on B grid

    real(dp), dimension(:), allocatable :: W_til               ! Basal water depth in till
    real(dp), dimension(:), allocatable :: W_til_next          ! Basal water depth in till at next timestep

    real(dp), dimension(:), allocatable :: P                   ! Pressure of the ice on the basal water
    real(dp), dimension(:), allocatable :: P_o                 ! Overburden pressure of the ice on the basal water
    real(dp), dimension(:), allocatable :: N                   ! Effective pressure (P_o - P)

    real(dp), dimension(:), allocatable :: m                   ! Total input of water to the basal system
    real(dp), dimension(:), allocatable :: m_extra             ! Additional input of water to the basal system because of time that has passed since the last basal hydrology model update

    real(dp), dimension(:), allocatable :: D                   ! Diffusivity
    real(dp), dimension(:), allocatable :: dD_dx_b             ! Derivative of D to x on B grid.
    real(dp), dimension(:), allocatable :: D_b                 ! Diffusivity on B grid

    real(dp), dimension(:), allocatable :: K                   ! Effective conductivity
    real(dp), dimension(:), allocatable :: dK_dx_b             ! Derivative of K to x on B grid.
    real(dp), dimension(:), allocatable :: K_b                 ! Effective conductivity on B grid.

    real(dp), dimension(:), allocatable :: u                   ! Velocity in x direction
    real(dp), dimension(:), allocatable :: du_dx_b             ! Derivative of u to x on B grid.
    real(dp), dimension(:), allocatable :: u_b                 ! Velocity in x direction on B grid.
    real(dp), dimension(:), allocatable :: u_c                 ! Velocity in x direction on C grid.
    real(dp), dimension(:), allocatable :: v                   ! Velocity in y direction
    real(dp), dimension(:), allocatable :: dv_dx_b             ! Derivative of v to x on B grid.
    real(dp), dimension(:), allocatable :: v_b                 ! Velocity in y direction on B grid.
    real(dp), dimension(:), allocatable :: v_c                 ! Velocity in y direction on C grid.

    ! Is this how it goes in UFEMISM as well? Compass indices do not really make sense here?
    real(dp), dimension(:), allocatable :: Q_e                 ! Normal component (east) of the advective flux VW
    real(dp), dimension(:), allocatable :: Q_w                 ! Normal component (west) of the advective flux VW
    real(dp), dimension(:), allocatable :: Q_n                 ! Normal component (north) of the advective flux VW
    real(dp), dimension(:), allocatable :: Q_s                 ! Normal component (south) of the advective flux VW
    real(dp), dimension(:), allocatable :: divQ                ! Divergence of the advective flux

    real(dp), dimension(:), allocatable :: q_til               ! Water flux towards till
    real(dp), dimension(:), allocatable :: q_water_layer       ! Water flux towards water layer

    real(dp), dimension(:), allocatable :: phi_0               ! Englacial porosity

    ! How do we handle constants? Are those also types somewhere?
    ! Seems to be in model_configuration.f90 and parameters.f90?

    real(dp), dimension(:), allocatable :: W_r                 ! Maximum roughness scale of basal topography
    real(dp), dimension(:), allocatable :: Y                   ! bed seperation

    real(dp), dimension(:), allocatable :: A                   ! Softness of the ice

    real(dp), dimension(:), allocatable :: psi                 ! Hydraulic potential

    real(dp), dimension(:), allocatable :: Z                   ! Sum of zeroth-order terms
    real(dp), dimension(:), allocatable :: C                   ! Closing rates
    real(dp), dimension(:), allocatable :: O                   ! Opening rates

    real(dp), dimension(:), allocatable :: R                   ! Subglacial water layer pressure: P + rho_w*g*b
    real(dp), dimension(:), allocatable :: dR_dx               ! Derivative of R to x on A grid.
    real(dp), dimension(:), allocatable :: dr_dx_b             ! Derivative of R to x on B grid.
    real(dp), dimension(:), allocatable :: dR_dy               ! Derivative of R to y on A grid.
    real(dp), dimension(:), allocatable :: dR_dy_b             ! Derivative of R to y on B grid.

    real(dp), allocatable               :: old_time            ! Time at previous timestep
    real(dp), allocatable               :: diff_time           ! Time since previous leg was run

    type(type_sparse_matrix_CSR_dp)     :: M_a_b               ! Matrix for going from grid a to grid b
    TYPE(type_sparse_matrix_CSR_dp)     :: M_b_c               ! Matrix for going from grid b to grid c
    TYPE(type_sparse_matrix_CSR_dp)     :: M_a_c               ! Matrix for going from grid b to grid c

    logical,  dimension(:), allocatable :: mask_a              ! Mask on b-grid on which to apply computation
    logical,  dimension(:), allocatable :: mask_b              ! Mask on b-grid on which to apply computation
    logical,  dimension(:), allocatable :: mask_W              ! Mask showing where the water layer is present

    ! Metadata
    CHARACTER(LEN=256)                      :: restart_filename            ! Name for generated restart file

    ! Timestepping
    REAL(dp), allocatable                   :: t_next
    real(dp), allocatable                   :: dt

    ! Ice model variables converted to m/s
    real(dp), dimension(:), allocatable     :: ice_u_base
    real(dp), dimension(:), allocatable     :: ice_v_base
    real(dp), dimension(:), allocatable     :: ice_w_base

    ! Water leaking back from till to water layer parameter
    real(dp), allocatable                   :: Cd

    ! LADDIE
    type(type_laddie_model)                       :: laddie
    type(type_laddie_forcing)                     :: forcing


  END TYPE type_basal_hydrology_model

CONTAINS

END MODULE basal_hydrology_model_types
