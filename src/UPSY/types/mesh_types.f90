module mesh_types

  ! The different data types used in the mesh modules

  use precisions, only: dp
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use parallel_array_info_type, only: type_par_arr_info
  use mpi_f08, only: MPI_WIN
  use deallocate_dist_shared_mod, only: deallocate_dist_shared

  implicit none

  private

  public :: type_mesh

  type type_mesh
    ! The unstructured triangular mesh.

  ! Mesh creation config parameters
  ! ===============================

    logical                                 :: do_singlecore_mesh_creation = .true.   ! Whether or not to let each process create the entire mesh themselves (as opposed to letting each process create a partial mesh and then merging them together - which is not yet supported...)
    real(dp)                                :: resolution_tolerance        = 1._dp

  ! Basic meta properties
  ! =====================

    character(len=256)                      :: name                          !           A nice name tag, e.g. mesh_ANT_00001
    real(dp)                                :: lambda_M                      ! [degrees] Oblique stereographic projection parameters
    real(dp)                                :: phi_M                         ! [degrees]
    real(dp)                                :: beta_stereo                   ! [degrees]
    real(dp)                                :: xmin                          ! [m]       x and y range of the square covered by the mesh
    real(dp)                                :: xmax                          ! [m]
    real(dp)                                :: ymin                          ! [m]
    real(dp)                                :: ymax                          ! [m]
    real(dp)                                :: tol_dist                      ! [m]       Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    integer                                 :: nV_mem                        !           Size of allocated memory for vertices
    integer                                 :: nTri_mem                      !           Size of allocated memory for triangles
    integer                                 :: nC_mem = 32                   !           Maximum allowed number of connections per vertex
    integer                                 :: nV                            !           Number of vertices
    integer                                 :: nTri                          !           Number of triangles

  ! Vertical grid (scaled coordinate)
  ! =================================

    character(len=1024)                     :: choice_zeta_grid      = 'irregular_log'  !           Choice of zeta grid ("regular", "irregular_log", "old_15_layer_zeta")
    integer                                 :: nz                    = 12               !           Number of vertical layers
    real(dp)                                :: zeta_irregular_log_R  = 10               !           Ratio between surface and base layer spacings
    real(dp), dimension(:    ), allocatable :: zeta                                     ! [0-1]     Scaled vertical coordinate
    real(dp), dimension(:    ), allocatable :: zeta_stag                                !           Staggered zeta grid

  ! Primary mesh data
  ! =================

    ! Vertex data
    integer                                 :: vi_SW, vi_SE, vi_NW, vi_NE    !           Indices of the vertices at the four corners of the domain
    real(dp), dimension(:,:  ), allocatable :: V                             ! [m]       The x,y-coordinates of all the vertices
    integer,  dimension(:    ), allocatable :: nC                            !           The number of other vertices this vertex is connected to
    integer,  dimension(:,:  ), allocatable :: C                             !           The list   of other vertices this vertex is connected to (ordered counter-clockwise, from edge to edge for edge vertices)
    integer,  dimension(:    ), allocatable :: niTri                         !           The number of triangles this vertex is a part of
    integer,  dimension(:,:  ), allocatable :: iTri                          !           The list   of triangles this vertex is a part of (ordered counter-clockwise)
    integer,  dimension(:    ), allocatable :: VBI                           ! [0-8]     Each vertex's border index; 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest

    ! Triangle data
    integer,  dimension(:,:  ), allocatable :: Tri                           !           The triangle array: Tri(ti) = [vi1, vi2, vi3] (vertices ordered counter-clockwise)
    real(dp), dimension(:,:  ), allocatable :: Tricc                         ! [m]       The X,Y-coordinates of each triangle's circumcenter
    integer,  dimension(:,:  ), allocatable :: TriC                          !           The (up to) three neighbour triangles (order across from 1st, 2nd and 3d vertex, respectively)

  ! Refinement data
  ! ===============

    logical,  dimension(:,:  ), allocatable :: check_Delaunay_map            !           Map    of triangle pairs that should be checked for the local Delaunay criterion
    integer,  dimension(:,:  ), allocatable :: check_Delaunay_stack          !           Stack  of triangle pairs that ...
    integer                                 :: check_Delaunay_stackN         !           Stack  of triangle pairs that ...
    integer,  dimension(:    ), allocatable :: refinement_map                !           Map    of triangles that should be checked for needing refinement
    integer,  dimension(:    ), allocatable :: refinement_stack              !           Stack  of triangles that ...
    integer                                 :: refinement_stackN             !           Number of triangles that ...
    integer,  dimension(:,:  ), allocatable :: Tri_li                        !           List of overlap ranges between triangles and line segments (for line-based mesh refinement)

  ! Secondary mesh data (everything that can be calculated after mesh creation is finished)
  ! =======================================================================================

    ! Derived geometry data
    real(dp), dimension(:    ), allocatable :: A                             ! [m^2]     The area             of each vertex's Voronoi cell
    real(dp), dimension(:,:  ), allocatable :: VorGC                         ! [m]       The geometric centre of each vertex's Voronoi cell
    real(dp), dimension(:    ), allocatable :: R                             ! [m]       The resolution of each vertex (defined as distance to nearest neighbour)
    real(dp), dimension(:,:  ), allocatable :: Cw                            ! [m]       The width of all vertex connections (= length of the shared Voronoi cell edge)
    real(dp), dimension(:,:  ), allocatable :: D_x                           ! [m]       x-component of vertex-vertex connections
    real(dp), dimension(:,:  ), allocatable :: D_y                           ! [m]       y-component of vertex-vertex connections
    real(dp), dimension(:,:  ), allocatable :: D                             ! [m]       absolute distance of vertex-vertex connections
    real(dp), dimension(:,:  ), allocatable :: TriCw                         ! [m]       The width of all triangle connections (= shared triangle edge)
    integer,  dimension(:    ), allocatable :: TriBI                         ! [0-8]     Each triangle's border index; 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest
    real(dp), dimension(:,:  ), allocatable :: TriGC                         ! [m]       The X,Y-coordinates of each triangle's geometric centre
    real(dp), dimension(:    ), allocatable :: TriA                          ! [m^2]     The area of each triangle
    real(dp), dimension(:,:  ), allocatable :: TriD_x                        ! [m]       x-component of triangle-triangle connections
    real(dp), dimension(:,:  ), allocatable :: TriD_y                        ! [m]       y-component of triangle-triangle connections
    real(dp), dimension(:,:  ), allocatable :: TriD                          ! [m]       absolute distance of triangle-triangle connections

    ! lon/lat coordinates
    real(dp), dimension(:    ), allocatable :: lat                           ! [degrees north] Latitude  of each vertex
    real(dp), dimension(:    ), allocatable :: lon                           ! [degrees east]  Longitude of each vertex

    ! Edges (c-grid)
    integer                                 :: nE                            !           Number of edges
    real(dp), dimension(:,:  ), allocatable :: E                             ! [m]       The x,y-coordinates of all the edges
    integer,  dimension(:,:  ), allocatable :: VE                            !           Vertex-to-edge connectivity list
    integer,  dimension(:,:  ), allocatable :: EV                            !           Edge-to-vertex connectivity list
    integer,  dimension(:,:  ), allocatable :: ETri                          !           Edge-to-triangle connectivity list
    integer,  dimension(:,:  ), allocatable :: TriE                          !           Triangle-to-edge connectivity list (order of TriC, across from 1st, 2nd, 3rd vertex)
    integer,  dimension(:    ), allocatable :: EBI                           ! [0-8]     Each edge's border index; 0 = free, 1 = north, 2 = northeast, ..., 8 = northwest
    real(dp), dimension(:    ), allocatable :: EA                            ! [m^2]     Area of each edge's "cell"

    ! Voronoi mesh
    integer                                 :: nVor                          !           Total number of Voronoi vertices
    integer,  dimension(:    ), allocatable :: vi2vori                       !           To which Voronoi vertex (if any) each vertex   corresponds (>0 only for the four corner vertices)
    integer,  dimension(:    ), allocatable :: ti2vori                       !           To which Voronoi vertex (if any) each triangle corresponds (>0 always)
    integer,  dimension(:    ), allocatable :: ei2vori                       !           To which Voronoi vertex (if any) each edge     corresponds (>0 only for border edges)
    integer,  dimension(:    ), allocatable :: vori2vi                       !           To which regular vertex (if any) each Voronoi vertex corresponds (vori2vi( vi2vori( vi)) = vi)
    integer,  dimension(:    ), allocatable :: vori2ti                       !           To which triangle       (if any) each Voronoi vertex corresponds (vori2ti( ti2vori( ti)) = ti)
    integer,  dimension(:    ), allocatable :: vori2ei                       !           To which edge           (if any) each Voronoi vertex corresponds (vori2ei( ei2vori( ei)) = ei)
    real(dp), dimension(:,:  ), allocatable :: Vor                           ! [m]       Coordinates of all the Voronoi vertices
    integer,  dimension(:    ), allocatable :: VornC                         !           Number  of other Voronoi vertices that each Voronoi vertex is connected to (2 for corners, 3 otherwise)
    integer,  dimension(:,:  ), allocatable :: VorC                          !           Indices of other Voronoi vertices that each Voronoi vertex is connected to
    integer,  dimension(:    ), allocatable :: nVVor                         !           For each regular vertex, the number of Voronoi vertices spanning its Voronoi cell
    integer,  dimension(:,:  ), allocatable :: VVor                          !           For each regular vertex, the indices of the Voronoi vertices spanning its Voronoi cell

    ! Parallelisation ranges
    integer                                         :: vi1, vi2, nV_loc              ! Each process "owns" nV_loc    vertices  vi1      - vi2     , so that nV_loc    = vi2      + 1 - vi1
    integer                                         :: ti1, ti2, nTri_loc            ! Each process "owns" nTri_loc  triangles ti1      - ti2     , so that nTri_loc  = ti2      + 1 - ti1
    integer                                         :: ei1, ei2, nE_loc              ! Each process "owns" nE_loc    edges     ei1      - ei2     , so that nE_loc    = ei2      + 1 - ei1

    integer,  dimension(:    ), contiguous, pointer :: V_owning_process   => null()  ! Which process owns each vertex
    integer,  dimension(:    ), contiguous, pointer :: Tri_owning_process => null()  ! Which process owns each triangle
    integer,  dimension(:    ), contiguous, pointer :: E_owning_process   => null()  ! Which process owns each edge
    integer,  dimension(:    ), contiguous, pointer :: V_owning_node      => null()  ! Which node owns each vertex
    integer,  dimension(:    ), contiguous, pointer :: Tri_owning_node    => null()  ! Which node owns each triangle
    integer,  dimension(:    ), contiguous, pointer :: E_owning_node      => null()  ! Which node owns each edge
    type(MPI_WIN) :: wV_owning_process, wTri_owning_process, wE_owning_process
    type(MPI_WIN) :: wV_owning_node, wTri_owning_node, wE_owning_node

    type(type_par_arr_info)                         :: pai_V                          ! Parallelisation info for vertex-based fields
    type(type_par_arr_info)                         :: pai_Tri                        ! Parallelisation info for triangle-based fields
    type(type_par_arr_info)                         :: pai_E                          ! Parallelisation info for edge-based fields

    ! Buffer shared memory
    real(dp), dimension(:    ), contiguous, pointer :: buffer1_d_a_nih  => null()     ! Pre-allocated buffer memory on the a-grid (vertices)
    real(dp), dimension(:    ), contiguous, pointer :: buffer2_d_a_nih  => null()
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer1_d_ak_nih => null()
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer2_d_ak_nih => null()
    real(dp), dimension(:    ), contiguous, pointer :: buffer1_d_b_nih  => null()     ! Pre-allocated buffer memory on the b-grid (triangles)
    real(dp), dimension(:    ), contiguous, pointer :: buffer2_d_b_nih  => null()
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer1_d_bk_nih => null()
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer2_d_bk_nih => null()
    real(dp), dimension(:    ), contiguous, pointer :: buffer1_d_c_nih  => null()     ! Pre-allocated buffer memory on the c-grid (edges)
    real(dp), dimension(:    ), contiguous, pointer :: buffer2_d_c_nih  => null()
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer1_d_ck_nih => null()
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer2_d_ck_nih => null()
    type(MPI_WIN) :: wbuffer1_d_a_nih, wbuffer2_d_a_nih, wbuffer1_d_ak_nih, wbuffer2_d_ak_nih
    type(MPI_WIN) :: wbuffer1_d_b_nih, wbuffer2_d_b_nih, wbuffer1_d_bk_nih, wbuffer2_d_bk_nih
    type(MPI_WIN) :: wbuffer1_d_c_nih, wbuffer2_d_c_nih, wbuffer1_d_ck_nih, wbuffer2_d_ck_nih

  ! Matrix operators
  ! ================

    ! Grid-cell-to-matrix-row translation tables

    ! a-grid (vertices)
    integer                                         :: nna, nnauv, nnak, nnaks, nnakuv, nnaksuv
    integer,  dimension(:    ), contiguous, pointer :: n2vi     => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2viuv   => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2vik    => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2vikuv  => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2viks   => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2viksuv => null()
    integer,  dimension(:    ), contiguous, pointer :: vi2n     => null()
    integer,  dimension(:,:  ), contiguous, pointer :: viuv2n   => null()
    integer,  dimension(:,:  ), contiguous, pointer :: vik2n    => null()
    integer,  dimension(:,:,:), contiguous, pointer :: vikuv2n  => null()
    integer,  dimension(:,:  ), contiguous, pointer :: viks2n   => null()
    integer,  dimension(:,:,:), contiguous, pointer :: viksuv2n => null()
    type(MPI_WIN) :: wn2vi, wn2viuv, wn2vik, wn2vikuv, wn2viks, wn2viksuv
    type(MPI_WIN) :: wvi2n, wviuv2n, wvik2n, wvikuv2n, wviks2n, wviksuv2n

    ! b-grid (triangles)
    integer                                         :: nnb, nnbuv, nnbk, nnbks, nnbkuv, nnbksuv
    integer,  dimension(:    ), contiguous, pointer :: n2ti     => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2tiuv   => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2tik    => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2tikuv  => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2tiks   => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2tiksuv => null()
    integer,  dimension(:    ), contiguous, pointer :: ti2n     => null()
    integer,  dimension(:,:  ), contiguous, pointer :: tiuv2n   => null()
    integer,  dimension(:,:  ), contiguous, pointer :: tik2n    => null()
    integer,  dimension(:,:,:), contiguous, pointer :: tikuv2n  => null()
    integer,  dimension(:,:  ), contiguous, pointer :: tiks2n   => null()
    integer,  dimension(:,:,:), contiguous, pointer :: tiksuv2n => null()
    type(MPI_WIN) :: wn2ti, wn2tiuv, wn2tik, wn2tikuv, wn2tiks, wn2tiksuv
    type(MPI_WIN) :: wti2n, wtiuv2n, wtik2n, wtikuv2n, wtiks2n, wtiksuv2n

    ! c-grid (edges)
    integer                                         :: nnc, nncuv, nnck, nncks, nnckuv, nncksuv
    integer,  dimension(:    ), contiguous, pointer :: n2ei     => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2eiuv   => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2eik    => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2eikuv  => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2eiks   => null()
    integer,  dimension(:,:  ), contiguous, pointer :: n2eiksuv => null()
    integer,  dimension(:    ), contiguous, pointer :: ei2n     => null()
    integer,  dimension(:,:  ), contiguous, pointer :: eiuv2n   => null()
    integer,  dimension(:,:  ), contiguous, pointer :: eik2n    => null()
    integer,  dimension(:,:,:), contiguous, pointer :: eikuv2n  => null()
    integer,  dimension(:,:  ), contiguous, pointer :: eiks2n   => null()
    integer,  dimension(:,:,:), contiguous, pointer :: eiksuv2n => null()
    type(MPI_WIN) :: wn2ei, wn2eiuv, wn2eik, wn2eikuv, wn2eiks, wn2eiksuv
    type(MPI_WIN) :: wei2n, weiuv2n, weik2n, weikuv2n, weiks2n, weiksuv2n

    ! Basic 2-D mapping and gradient operators

    ! a-grid (vertices) to a-grid (vertices)
    type(type_CSR_matrix_dp)         :: M_ddx_a_a
    type(type_CSR_matrix_dp)         :: M_ddy_a_a
    ! a-grid (vertices) to b-grid (triangles)
    type(type_CSR_matrix_dp)         :: M_map_a_b
    type(type_CSR_matrix_dp)         :: M_ddx_a_b
    type(type_CSR_matrix_dp)         :: M_ddy_a_b
    ! b-grid (triangles) to a-grid (vertices)
    type(type_CSR_matrix_dp)         :: M_map_b_a
    type(type_CSR_matrix_dp)         :: M_ddx_b_a
    type(type_CSR_matrix_dp)         :: M_ddy_b_a
    ! b-grid (triangles) to b-grid (triangles)
    type(type_CSR_matrix_dp)         :: M_ddx_b_b
    type(type_CSR_matrix_dp)         :: M_ddy_b_b

    ! b-grid (triangles) to b-grid (triangles), 2nd-order accurate
    type(type_CSR_matrix_dp)         :: M2_ddx_b_b
    type(type_CSR_matrix_dp)         :: M2_ddy_b_b
    type(type_CSR_matrix_dp)         :: M2_d2dx2_b_b
    type(type_CSR_matrix_dp)         :: M2_d2dxdy_b_b
    type(type_CSR_matrix_dp)         :: M2_d2dy2_b_b

    ! Operators on the zeta grids
    type(type_CSR_matrix_dp)         :: M_ddzeta_k_k_1D
    type(type_CSR_matrix_dp)         :: M_d2dzeta2_k_k_1D
    type(type_CSR_matrix_dp)         :: M_map_k_ks_1D
    type(type_CSR_matrix_dp)         :: M_ddzeta_k_ks_1D
    type(type_CSR_matrix_dp)         :: M_map_ks_k_1D
    type(type_CSR_matrix_dp)         :: M_ddzeta_ks_k_1D

    ! Zeta operators in tridiagonal form for efficient use in thermodynamics
    real(dp), dimension(:    ), allocatable :: M_ddzeta_k_k_ldiag
    real(dp), dimension(:    ), allocatable :: M_ddzeta_k_k_diag
    real(dp), dimension(:    ), allocatable :: M_ddzeta_k_k_udiag
    real(dp), dimension(:    ), allocatable :: M_d2dzeta2_k_k_ldiag
    real(dp), dimension(:    ), allocatable :: M_d2dzeta2_k_k_diag
    real(dp), dimension(:    ), allocatable :: M_d2dzeta2_k_k_udiag

    ! 3-D gradient operators

    ! bk to ak (for calculating the horizontal stretch/shear strain rates in the BPA)
    type(type_CSR_matrix_dp)         :: M_ddx_bk_ak
    type(type_CSR_matrix_dp)         :: M_ddy_bk_ak

    ! ak to bk (for calculating the horizontal gradients of the effective viscosity in the BPA)
    type(type_CSR_matrix_dp)         :: M_ddx_ak_bk
    type(type_CSR_matrix_dp)         :: M_ddy_ak_bk

    ! bk to bks (for calculating the vertical shear strain rates in the BPA)
    type(type_CSR_matrix_dp)         :: M_ddz_bk_bks

    ! bks to bk (for calculating (the vertical gradient of) the effective viscosity in the BPA)
    type(type_CSR_matrix_dp)         :: M_map_bks_bk
    type(type_CSR_matrix_dp)         :: M_ddz_bks_bk

    ! Map between the bks-grid and the ak-grid (for calculating strain rates in the BPA)
    type(type_CSR_matrix_dp)         :: M_map_bks_ak
    type(type_CSR_matrix_dp)         :: M_map_ak_bks

    ! bk to bk (for constructing the BPA stiffness matrix)
    type(type_CSR_matrix_dp)         :: M2_ddx_bk_bk
    type(type_CSR_matrix_dp)         :: M2_ddy_bk_bk
    type(type_CSR_matrix_dp)         :: M2_d2dx2_bk_bk
    type(type_CSR_matrix_dp)         :: M2_d2dxdy_bk_bk
    type(type_CSR_matrix_dp)         :: M2_d2dy2_bk_bk
    type(type_CSR_matrix_dp)         :: M2_ddz_bk_bk
    type(type_CSR_matrix_dp)         :: M2_d2dz2_bk_bk

  contains

    final :: finalise

  end type type_mesh

contains

  subroutine finalise( mesh)
    !< Deallocate hybrid distributed/shared components (as those are
    !< not captured by automatic finalisation)

    type(type_mesh) :: mesh

    ! Parallelisation ranges
    if (associated( mesh%V_owning_process  )) call deallocate_dist_shared( mesh%V_owning_process  , mesh%wV_owning_process  )
    if (associated( mesh%Tri_owning_process)) call deallocate_dist_shared( mesh%Tri_owning_process, mesh%wTri_owning_process)
    if (associated( mesh%E_owning_process  )) call deallocate_dist_shared( mesh%E_owning_process  , mesh%wE_owning_process  )
    if (associated( mesh%V_owning_node     )) call deallocate_dist_shared( mesh%V_owning_node     , mesh%wV_owning_node     )
    if (associated( mesh%Tri_owning_node   )) call deallocate_dist_shared( mesh%Tri_owning_node   , mesh%wTri_owning_node   )
    if (associated( mesh%E_owning_node     )) call deallocate_dist_shared( mesh%E_owning_node     , mesh%wE_owning_node     )

    ! Buffer shared memory
    if (associated( mesh%buffer1_d_a_nih )) call deallocate_dist_shared( mesh%buffer1_d_a_nih , mesh%wbuffer1_d_a_nih )
    if (associated( mesh%buffer2_d_a_nih )) call deallocate_dist_shared( mesh%buffer2_d_a_nih , mesh%wbuffer2_d_a_nih )
    if (associated( mesh%buffer1_d_ak_nih)) call deallocate_dist_shared( mesh%buffer1_d_ak_nih, mesh%wbuffer1_d_ak_nih)
    if (associated( mesh%buffer2_d_ak_nih)) call deallocate_dist_shared( mesh%buffer2_d_ak_nih, mesh%wbuffer2_d_ak_nih)
    if (associated( mesh%buffer1_d_b_nih )) call deallocate_dist_shared( mesh%buffer1_d_b_nih , mesh%wbuffer1_d_b_nih )
    if (associated( mesh%buffer2_d_b_nih )) call deallocate_dist_shared( mesh%buffer2_d_b_nih , mesh%wbuffer2_d_b_nih )
    if (associated( mesh%buffer1_d_bk_nih)) call deallocate_dist_shared( mesh%buffer1_d_bk_nih, mesh%wbuffer1_d_bk_nih)
    if (associated( mesh%buffer2_d_bk_nih)) call deallocate_dist_shared( mesh%buffer2_d_bk_nih, mesh%wbuffer2_d_bk_nih)
    if (associated( mesh%buffer1_d_c_nih )) call deallocate_dist_shared( mesh%buffer1_d_c_nih , mesh%wbuffer1_d_c_nih )
    if (associated( mesh%buffer2_d_c_nih )) call deallocate_dist_shared( mesh%buffer2_d_c_nih , mesh%wbuffer2_d_c_nih )
    if (associated( mesh%buffer1_d_ck_nih)) call deallocate_dist_shared( mesh%buffer1_d_ck_nih, mesh%wbuffer1_d_ck_nih)
    if (associated( mesh%buffer2_d_ck_nih)) call deallocate_dist_shared( mesh%buffer2_d_ck_nih, mesh%wbuffer2_d_ck_nih)

    ! Grid-cell-to-matrix-row translation tables

    ! a-grid (vertices)
    if (associated( mesh%n2vi    )) call deallocate_dist_shared( mesh%n2vi    , mesh%wn2vi    )
    if (associated( mesh%n2viuv  )) call deallocate_dist_shared( mesh%n2viuv  , mesh%wn2viuv  )
    if (associated( mesh%n2vik   )) call deallocate_dist_shared( mesh%n2vik   , mesh%wn2vik   )
    if (associated( mesh%n2vikuv )) call deallocate_dist_shared( mesh%n2vikuv , mesh%wn2vikuv )
    if (associated( mesh%n2viks  )) call deallocate_dist_shared( mesh%n2viks  , mesh%wn2viks  )
    if (associated( mesh%n2viksuv)) call deallocate_dist_shared( mesh%n2viksuv, mesh%wn2viksuv)
    if (associated( mesh%vi2n    )) call deallocate_dist_shared( mesh%vi2n    , mesh%wvi2n    )
    if (associated( mesh%viuv2n  )) call deallocate_dist_shared( mesh%viuv2n  , mesh%wviuv2n  )
    if (associated( mesh%vik2n   )) call deallocate_dist_shared( mesh%vik2n   , mesh%wvik2n   )
    if (associated( mesh%vikuv2n )) call deallocate_dist_shared( mesh%vikuv2n , mesh%wvikuv2n )
    if (associated( mesh%viks2n  )) call deallocate_dist_shared( mesh%viks2n  , mesh%wviks2n  )
    if (associated( mesh%viksuv2n)) call deallocate_dist_shared( mesh%viksuv2n, mesh%wviksuv2n)

    ! b-grid (triangles)
    if (associated( mesh%n2ti    )) call deallocate_dist_shared( mesh%n2ti    , mesh%wn2ti    )
    if (associated( mesh%n2tiuv  )) call deallocate_dist_shared( mesh%n2tiuv  , mesh%wn2tiuv  )
    if (associated( mesh%n2tik   )) call deallocate_dist_shared( mesh%n2tik   , mesh%wn2tik   )
    if (associated( mesh%n2tikuv )) call deallocate_dist_shared( mesh%n2tikuv , mesh%wn2tikuv )
    if (associated( mesh%n2tiks  )) call deallocate_dist_shared( mesh%n2tiks  , mesh%wn2tiks  )
    if (associated( mesh%n2tiksuv)) call deallocate_dist_shared( mesh%n2tiksuv, mesh%wn2tiksuv)
    if (associated( mesh%ti2n    )) call deallocate_dist_shared( mesh%ti2n    , mesh%wti2n    )
    if (associated( mesh%tiuv2n  )) call deallocate_dist_shared( mesh%tiuv2n  , mesh%wtiuv2n  )
    if (associated( mesh%tik2n   )) call deallocate_dist_shared( mesh%tik2n   , mesh%wtik2n   )
    if (associated( mesh%tikuv2n )) call deallocate_dist_shared( mesh%tikuv2n , mesh%wtikuv2n )
    if (associated( mesh%tiks2n  )) call deallocate_dist_shared( mesh%tiks2n  , mesh%wtiks2n  )
    if (associated( mesh%tiksuv2n)) call deallocate_dist_shared( mesh%tiksuv2n, mesh%wtiksuv2n)

    ! c-grid (edges)
    if (associated( mesh%n2ei    )) call deallocate_dist_shared( mesh%n2ei    , mesh%wn2ei    )
    if (associated( mesh%n2eiuv  )) call deallocate_dist_shared( mesh%n2eiuv  , mesh%wn2eiuv  )
    if (associated( mesh%n2eik   )) call deallocate_dist_shared( mesh%n2eik   , mesh%wn2eik   )
    if (associated( mesh%n2eikuv )) call deallocate_dist_shared( mesh%n2eikuv , mesh%wn2eikuv )
    if (associated( mesh%n2eiks  )) call deallocate_dist_shared( mesh%n2eiks  , mesh%wn2eiks  )
    if (associated( mesh%n2eiksuv)) call deallocate_dist_shared( mesh%n2eiksuv, mesh%wn2eiksuv)
    if (associated( mesh%ei2n    )) call deallocate_dist_shared( mesh%ei2n    , mesh%wei2n    )
    if (associated( mesh%eiuv2n  )) call deallocate_dist_shared( mesh%eiuv2n  , mesh%weiuv2n  )
    if (associated( mesh%eik2n   )) call deallocate_dist_shared( mesh%eik2n   , mesh%weik2n   )
    if (associated( mesh%eikuv2n )) call deallocate_dist_shared( mesh%eikuv2n , mesh%weikuv2n )
    if (associated( mesh%eiks2n  )) call deallocate_dist_shared( mesh%eiks2n  , mesh%weiks2n  )
    if (associated( mesh%eiksuv2n)) call deallocate_dist_shared( mesh%eiksuv2n, mesh%weiksuv2n)

  end subroutine finalise

end module mesh_types
