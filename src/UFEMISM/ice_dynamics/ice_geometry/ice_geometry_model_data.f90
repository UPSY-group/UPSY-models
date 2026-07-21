module ice_geometry_model_data

  use models_basic, only: atype_model
  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_ice_geometry_model_data

  type, abstract, extends(atype_model) :: atype_ice_geometry_model_data

    ! Primary ice geometry fields
    real(dp), dimension(:  ), allocatable :: Hi                      ! [m]       Ice thickness
    real(dp), dimension(:  ), allocatable :: Hb                      ! [m]       Bedrock elevation (w.r.t. PD sea level)
    real(dp), dimension(:  ), allocatable :: SL                      ! [m]       Geoid elevation   (w.r.t. PD sea level)
    type(MPI_WIN) :: wHi, wHb, wSL

    ! Secondary ice geometry fields
    real(dp), dimension(:  ), allocatable :: Hs                      ! [m]       Ice surface elevation (w.r.t. PD sea level)
    real(dp), dimension(:  ), allocatable :: Hib                     ! [m]       Ice base elevation    (w.r.t. PD sea level)
    real(dp), dimension(:  ), allocatable :: TAF                     ! [m]       Thickness above floatation
    ! real(dp), dimension(:  ), allocatable :: Ho                      ! [m]       Height of water column at ice front
    ! real(dp), dimension(:  ), allocatable :: Hi_eff                  ! [m]       Effective ice thickness
    ! real(dp), dimension(:  ), allocatable :: dHb                     ! [m]       Change in bedrock elevation w.r.t. reference geometry
    ! real(dp), dimension(:  ), allocatable :: fraction_margin         ! [0-1]     Sub-grid ice-filled fraction
    ! real(dp), dimension(:  ), allocatable :: fraction_gr             ! [0-1]     Sub-grid grounded fraction
    ! real(dp), dimension(:  ), allocatable :: fraction_gr_b           ! [0-1]     Sub-grid grounded fraction on b-grid (triangles)
    ! real(dp), dimension(:,:), allocatable :: bedrock_cdf             ! [m]       Bedrock cumulative density function
    ! real(dp), dimension(:,:), allocatable :: bedrock_cdf_b           ! [m]       Bedrock cumulative density function on b-grid (triangles)
    ! type(MPI_WIN) :: wHs, wHib, wTAF, wHo, wHi_eff, wdHb
    ! type(MPI_WIN) :: wfraction_margin, wfraction_gr, wfraction_gr_b, wbedrock_cdf, wbedrock_cdf_b

    ! Ice masks
    logical,  dimension(:  ), allocatable :: mask_icefree_land       ! T: ice-free land , F: otherwise
    logical,  dimension(:  ), allocatable :: mask_icefree_ocean      ! T: ice-free ocean, F: otherwise
    logical,  dimension(:  ), allocatable :: mask_grounded_ice       ! T: grounded ice  , F: otherwise
    logical,  dimension(:  ), allocatable :: mask_floating_ice       ! T: floating ice  , F: otherwise
    logical,  dimension(:  ), allocatable :: mask_margin             ! T: ice next to ice-free, F: otherwise
    logical,  dimension(:  ), allocatable :: mask_gl_gr              ! T: grounded ice next to floating ice, F: otherwise
    logical,  dimension(:  ), allocatable :: mask_gl_fl              ! T: floating ice next to grounded ice, F: otherwise
    logical,  dimension(:  ), allocatable :: mask_cf_gr              ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    logical,  dimension(:  ), allocatable :: mask_cf_fl              ! T: floating ice next to ice-free water (sea or lake), F: otherwise
    logical,  dimension(:  ), allocatable :: mask_coastline          ! T: ice-free land next to ice-free ocean, F: otherwise
    integer,  dimension(:  ), allocatable :: mask
    type(MPI_WIN) :: wmask_icefree_land, wmask_icefree_ocean, wmask_grounded_ice, wmask_floating_ice
    type(MPI_WIN) :: wmask_margin, wmask_gl_gr, wmask_gl_fl, wmask_cf_gr, wmask_cf_fl, wmask_coastline, wmask

  end type atype_ice_geometry_model_data

end module ice_geometry_model_data