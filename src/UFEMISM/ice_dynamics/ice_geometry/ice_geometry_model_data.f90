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

    ! Secondary ice geometry fields
    real(dp), dimension(:  ), allocatable :: Hs                      ! [m]       Ice surface elevation (w.r.t. PD sea level)
    real(dp), dimension(:  ), allocatable :: Hib                     ! [m]       Ice base elevation    (w.r.t. PD sea level)
    real(dp), dimension(:  ), allocatable :: TAF                     ! [m]       Thickness above floatation
    real(dp), dimension(:  ), allocatable :: Hi_eff                  ! [m]       Effective ice thickness
    real(dp), dimension(:  ), allocatable :: Hs_slope                ! [-]       Absolute surface gradient
    real(dp), dimension(:  ), allocatable :: Ho                      ! [m]       Depth of ocean column adjacent to the ice front

    ! Horizontal derivatives
    real(dp), dimension(:  ), contiguous, pointer :: dHib_dx_b             => null()  ! [-]       Horizontal derivative of ice draft on b-grid
    real(dp), dimension(:  ), contiguous, pointer :: dHib_dy_b             => null()  ! [-]       Horizontal derivative of ice draft on b-grid
    type(MPI_WIN) :: wdHib_dx_b, wdHib_dy_b

    ! Sub-grid bedrock cumulative density functions (CDFs)
    real(dp), dimension(:,:), contiguous, pointer :: bedrock_cdf           => null()  ! [-]       Sub-grid bedrock cumulative density functions on the a-grid (vertices)
    real(dp), dimension(:,:), contiguous, pointer :: bedrock_cdf_b         => null()  ! [-]       Sub-grid bedrock cumulative density functions on the b-grid (triangles)
    type(MPI_WIN) :: wbedrock_cdf, wbedrock_cdf_b

    ! Area fractions
    real(dp), dimension(:  ), contiguous, pointer :: fraction_gr           => null()  ! [0-1]     Grounded area fractions of vertices
    real(dp), dimension(:  ), contiguous, pointer :: fraction_gr_b         => null()  ! [0-1]     Grounded area fractions of triangles
    real(dp), dimension(:  ), contiguous, pointer :: fraction_margin       => null()  ! [0-1]     Ice-covered area fractions of ice margins
    type(MPI_WIN) :: wfraction_gr, wfraction_gr_b, wfraction_margin

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

  end type atype_ice_geometry_model_data

end module ice_geometry_model_data