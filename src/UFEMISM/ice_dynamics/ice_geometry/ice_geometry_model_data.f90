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
    real(dp), dimension(:  ), contiguous, pointer :: Ho                    => null()  ! [m]       Depth of ocean column adjacent to the ice front
    type(MPI_WIN) :: wHo

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
    logical,  dimension(:  ), contiguous, pointer :: mask_icefree_land     => null()  ! T: ice-free land , F: otherwise
    logical,  dimension(:  ), contiguous, pointer :: mask_icefree_ocean    => null()  ! T: ice-free ocean, F: otherwise
    logical,  dimension(:  ), contiguous, pointer :: mask_grounded_ice     => null()  ! T: grounded ice  , F: otherwise
    logical,  dimension(:  ), contiguous, pointer :: mask_floating_ice     => null()  ! T: floating ice  , F: otherwise
    logical,  dimension(:  ), contiguous, pointer :: mask_margin           => null()  ! T: ice next to ice-free, F: otherwise
    logical,  dimension(:  ), contiguous, pointer :: mask_gl_gr            => null()  ! T: grounded ice next to floating ice, F: otherwise
    logical,  dimension(:  ), contiguous, pointer :: mask_gl_fl            => null()  ! T: floating ice next to grounded ice, F: otherwise
    logical,  dimension(:  ), contiguous, pointer :: mask_cf_gr            => null()  ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    logical,  dimension(:  ), contiguous, pointer :: mask_cf_fl            => null()  ! T: floating ice next to ice-free water (sea or lake), F: otherwise
    logical,  dimension(:  ), contiguous, pointer :: mask_coastline        => null()  ! T: ice-free land next to ice-free ocean, F: otherwise
    integer,  dimension(:  ), contiguous, pointer :: mask                  => null()
    type(MPI_WIN) :: wmask_icefree_land, wmask_icefree_ocean, wmask_grounded_ice, wmask_floating_ice, wmask_margin
    type(MPI_WIN) :: wmask_gl_gr, wmask_gl_fl, wmask_cf_gr, wmask_cf_fl, wmask_coastline, wmask

  end type atype_ice_geometry_model_data

end module ice_geometry_model_data