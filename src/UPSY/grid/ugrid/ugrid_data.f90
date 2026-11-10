module ugrid_data

  use precisions, only: dp
  use basic_grid_mod, only: atype_grid
  use mpi_f08, only: MPI_WIN
  use parallel_array_info_type, only: type_par_arr_info

  implicit none

  private

  public :: atype_ugrid_data

  type, abstract, extends(atype_grid) :: atype_ugrid_data

    integer                                       :: nNodes                                        !           Number of nodes in the ugrid
    integer                                       :: nFaces                                        !           Number of faces in the ugrid

    real(dp), dimension(:,:), contiguous, pointer :: node_coordinates                   => null()  ! [m]       x/y-coordinates of each node
    real(dp), dimension(:,:), contiguous, pointer :: face_coordinates                   => null()  ! [m]       x/y-coordinates of each face center
    type(MPI_WIN) :: wnode_coordinates, wface_coordinates

    integer                                       :: face_node_connectivity_fill_value             !           Fill value for unused slots
    integer,  dimension(:),   contiguous, pointer :: face_node_connectivity_n           => null()  ! [nFaces]  Number of nodes spanning each face
    integer                                       :: face_node_connectivity_n_max                  !           Maximum number of nodes spanning a face
    integer,  dimension(:,:), contiguous, pointer :: face_node_connectivity             => null()  ! [nFaces-by-face_node_connectivity_n_max] Indices of nodes spanning each face
    type(MPI_WIN) :: wface_node_connectivity_n, wface_node_connectivity

    integer                                       :: face_face_connectivity_fill_value             !           Fill value for unused slots
    integer,  dimension(:),   contiguous, pointer :: face_face_connectivity_n           => null()  ! [nFaces]  Number of faces adjacent to each face
    integer                                       :: face_face_connectivity_n_max                  !           Maximum number of faces adjacent to a face
    integer,  dimension(:,:), contiguous, pointer :: face_face_connectivity             => null()  ! [nFaces-by-max(face_face_connectivity)] Indices of faces adjacent to each face
    type(MPI_WIN) :: wface_face_connectivity_n, wface_face_connectivity

    type(type_par_arr_info)                       :: pai                                           ! Parallelisation info

  end type atype_ugrid_data

end module ugrid_data
