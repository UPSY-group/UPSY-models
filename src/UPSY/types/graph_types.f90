module graph_types

  use precisions, only: dp
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use parallel_array_info_type, only: type_par_arr_info
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_graph, type_graph_pair

  type type_graph

    ! Metadata
    character(len=1024)                     :: parent_mesh_name              !<     Name of this graph's parent mesh
    real(dp)                                :: xmin, xmax, ymin, ymax        !< [m] Domain boundaries
    integer                                 :: n                             !<     Number of nodes in the graph
    integer                                 :: nC_mem                        !<     Maximum allowed number of connections per node

    ! Mapping between parent mesh and graph
    integer,  dimension(:    ), allocatable :: ni2vi                         !<     Graph node (ni) to parent mesh vertex   (vi) translation table
    integer,  dimension(:    ), allocatable :: ni2ti                         !<     Graph node (ni) to parent mesh triangle (ti) translation table
    integer,  dimension(:    ), allocatable :: ni2ei                         !<     Graph node (ni) to parent mesh edge     (ei) translation table
    integer,  dimension(:    ), allocatable :: vi2ni                         !<     Parent mesh vertex   (vi) to graph node (ni) translation table
    integer,  dimension(:    ), allocatable :: ti2ni                         !<     Parent mesh triangle (ti) to graph node (ti) translation table
    integer,  dimension(:    ), allocatable :: ei2ni                         !<     Parent mesh triangle (ei) to graph node (ti) translation table

    ! Node coordinates and connectivity
    real(dp), dimension(:,:  ), allocatable :: V                             !< [m] Node coordinates
    integer,  dimension(:    ), allocatable :: nC                            !<     Number  of other nodes each node is connected to
    integer,  dimension(:,:  ), allocatable :: C                             !<     Indices of other nodes each node is connected to

    ! Border nodes
    logical,  dimension(:    ), allocatable :: is_border                     !<     Whether or not a node is a border node
    real(dp), dimension(:,:  ), allocatable :: border_nhat                   !<     Outward unit normal vector at each border node

    ! Parallelisation ranges
    integer                                 :: ni1, ni2, n_loc               !< Each process "owns" n_loc nodes ni1 - ni2, so that n_loc = ni2 + 1 - ni1

    integer,  dimension(:    ), allocatable :: owning_process                !< Which process owns each node
    integer,  dimension(:    ), allocatable :: owning_node                   !< Which HPC node owns each graph node (yes, this is confusing...)

    type(type_par_arr_info)                 :: pai                           !< Parallelisation info for graph-based fields

    real(dp), dimension(:    ), contiguous, pointer :: buffer1_g_nih  => null()     !< Pre-allocated buffer memory on the graph
    real(dp), dimension(:    ), contiguous, pointer :: buffer2_g_nih  => null()     !< Pre-allocated buffer memory on the graph
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer1_gk_nih => null()     !< Pre-allocated buffer memory on the graph
    real(dp), dimension(:,:  ), contiguous, pointer :: buffer2_gk_nih => null()     !< Pre-allocated buffer memory on the graph
    type(MPI_WIN) :: wbuffer1_g_nih, wbuffer2_g_nih                                 !< MPI window to pre-allocated buffer memory on the graph
    type(MPI_WIN) :: wbuffer1_gk_nih, wbuffer2_gk_nih                               !< MPI window to pre-allocated buffer memory on the graph

  end type type_graph

  type type_graph_pair

    type(type_graph)                :: graph_a, graph_b

    type(type_sparse_matrix_CSR_dp) :: M_ddx_b_b
    type(type_sparse_matrix_CSR_dp) :: M_ddy_b_b

    type(type_sparse_matrix_CSR_dp) :: M2_ddx_b_b
    type(type_sparse_matrix_CSR_dp) :: M2_ddy_b_b
    type(type_sparse_matrix_CSR_dp) :: M2_d2dx2_b_b
    type(type_sparse_matrix_CSR_dp) :: M2_d2dxdy_b_b
    type(type_sparse_matrix_CSR_dp) :: M2_d2dy2_b_b

    type(type_sparse_matrix_CSR_dp) :: M_map_a_b
    type(type_sparse_matrix_CSR_dp) :: M_ddx_a_b
    type(type_sparse_matrix_CSR_dp) :: M_ddy_a_b

    type(type_sparse_matrix_CSR_dp) :: M_map_b_a
    type(type_sparse_matrix_CSR_dp) :: M_ddx_b_a
    type(type_sparse_matrix_CSR_dp) :: M_ddy_b_a

    integer, dimension(:,:), allocatable :: biuv2n
    integer, dimension(:,:), allocatable :: n2biuv

  end type type_graph_pair

end module graph_types