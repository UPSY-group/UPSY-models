module ugrid_data

  use precisions, only: dp
  use basic_grid_mod, only: atype_grid
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_ugrid_data

  type, abstract, extends(atype_grid) :: atype_ugrid_data

    integer                                       :: nNodes
    integer                                       :: nFaces

    real(dp), dimension(:,:), contiguous, pointer :: node_coordinates
    real(dp), dimension(:,:), contiguous, pointer :: face_coordinates
    type(MPI_WIN) :: wnode_coordinates, wface_coordinates

    integer                                       :: face_node_connectivity_fill_value
    integer                                       :: face_node_connectivity_n
    integer,  dimension(:,:), contiguous, pointer :: face_node_connectivity
    type(MPI_WIN) :: wface_node_connectivity

    integer                                       :: face_face_connectivity_fill_value
    integer                                       :: face_face_connectivity_n
    integer,  dimension(:,:), contiguous, pointer :: face_face_connectivity
    type(MPI_WIN) :: wface_face_connectivity

  end type atype_ugrid_data

end module ugrid_data
