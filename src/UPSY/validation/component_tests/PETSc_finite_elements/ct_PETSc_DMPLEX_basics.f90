module ct_PETSc_DMPLEX_basics

#include <petsc/finclude/petscdm.h>
#include <petsc/finclude/petscdmplex.h>

  use precisions, only: dp
  use mpi, only: MPI_COMM_WORLD         ! The newer mpi_f08 uses a wrapper type for communicators, which PETSc doesn't support...
  use petscdmplex, only: tDM, DMPlexCreate, DMPlexSetChart, DMPlexSetConeSize, &
    DMSetUp, DMPlexSetCone, DMPlexSymmetrize, DMPlexStratify, DMDestroy, &
    tVec, DMSetCoordinates
  use petsc_basic, only: vec_double2petsc

  implicit none

  private

  public :: create_simple_DMPLEX

contains

  subroutine create_simple_DMPLEX
    ! Implements the example from https://petsc.org/release/manual/dmplex/#ch-unstructured
    !
    ! The dummy mesh looks like this:
    !
    !                   [v2;p4]
    !                  /   |   \
    !                /     |     \
    !              /       |       \
    !         [e2;p8]      |      [e3;p9]
    !          /           |           \
    !        /             |             \
    !      /               |               \
    ! [v0;p2]  [f0;p0]  [e1;p7]  [f1;p1]  [v3;p5]
    !      \               |               /
    !        \             |             /
    !          \           |           /
    !         [e0;p6]      |      [e4;p10]
    !              \       |       /
    !                \     |     /
    !                  \   |   /
    !                   [v1;p3]

    type(tDM)                           :: dm
    integer                             :: perr
    real(dp), dimension(0:5 ,2)         :: vertex_coords
    real(dp), dimension(0:10,2)         :: coords_nby2
    real(dp), dimension(:), allocatable :: coords_2n
    integer                             :: n,i
    real(dp)                            :: x,y
    type(tVec)                          :: coords

    ! Create a DMPlex object
    call DMPlexCreate( MPI_COMM_WORLD, dm, perr)

    ! The example mesh has 11 'points' (vertices, faces, and edges)
    call DMPlexSetChart( dm, 0, 11, perr)

    ! In 2 dimensions the convention is to first number faces, then vertices, and then edges.
    ! PETSc indexes from zero
    !    DMPlexSetConeSize( dm, point, number of points that cover the point)
    call DMPlexSetConeSize( dm, 0, 3, perr)  ! Points 0 and 1 are faces, each connected to 3 edges.
    call DMPlexSetConeSize( dm, 1, 3, perr)
    call DMPlexSetConeSize( dm, 6, 2, perr)  ! Points 6-10 are edges, each connected to 2 vertices.
    call DMPlexSetConeSize( dm, 7, 2, perr)
    call DMPlexSetConeSize( dm, 8, 2, perr)
    call DMPlexSetConeSize( dm, 9, 2, perr)
    call DMPlexSetConeSize( dm, 10, 2, perr)
    call DMSetUp( dm, perr)

    !    DMPlexSetCone( dm, point, [points that cover the point])
    call DMPlexSetCone( dm,  0, [6, 7,  8], perr) ! Point  0 (= face 0) is connected to points [6,7, 8] (= edges [0,1,2])
    call DMPlexSetCone( dm,  1, [7, 9, 10], perr) ! Point  1 (= face 1) is connected to points [7,9,10] (= edges [1,3,4])  ...so not necessarily counter-clockwise
    call DMPlexSetCone( dm,  6, [2, 3    ], perr) ! Point  6 (= edge 0) is connected to points [2,3] (= vertices [0,1])
    call DMPlexSetCone( dm,  7, [3, 4    ], perr) ! Point  7 (= edge 1) is connected to points [3,4] (= vertices [1,2])
    call DMPlexSetCone( dm,  8, [4, 2    ], perr) ! Point  8 (= edge 2) is connected to points [4,2] (= vertices [2,0])
    call DMPlexSetCone( dm,  9, [4, 5    ], perr) ! Point  9 (= edge 3) is connected to points [4,5] (= vertices [2,3])
    call DMPlexSetCone( dm, 10, [5, 3    ], perr) ! Point 10 (= edge 4) is connected to points [5,3] (= vertices [3,1])

    ! Let PETSc automatically figure out the 'supports', i.e. the backward connections (so each
    ! vertex knows which edges and faces it spans)
    call DMPlexSymmetrize( dm, perr)

    ! In order to support efficient queries, we construct fast search structures
    ! and indices for the different types of points
    call DMPlexStratify( dm, perr)

    ! ADDITION: set coordinates

    ! Basic vertex coordinates
    vertex_coords( 0,:) = [-1._dp,  0._dp]
    vertex_coords( 1,:) = [ 0._dp, -1._dp]
    vertex_coords( 2,:) = [ 0._dp,  1._dp]
    vertex_coords( 3,:) = [ 1._dp,  0._dp]

    ! DMPLEX point coordinates

    ! Point 0 = face 0, spanned by vertices [0,1,2]
    coords_nby2(  0,:) = (vertex_coords( 0,:) + vertex_coords( 1,:) + vertex_coords( 2,:)) / 3._dp
    ! Point 1 = face 1, spanned by vertices [0,1,2]
    coords_nby2(  1,:) = (vertex_coords( 1,:) + vertex_coords( 2,:) + vertex_coords( 3,:)) / 3._dp
    ! Point 2 = vertex 0
    coords_nby2(  2,:) = vertex_coords( 0,:)
    ! Point 3 = vertex 1
    coords_nby2(  3,:) = vertex_coords( 1,:)
    ! Point 4 = vertex 2
    coords_nby2(  4,:) = vertex_coords( 2,:)
    ! Point 5 = vertex 3
    coords_nby2(  5,:) = vertex_coords( 3,:)
    ! Point 6 = edge 0, spanned by vertices [0,1]
    coords_nby2(  6,:) = (vertex_coords( 0,:) + vertex_coords( 1,:)) / 2._dp
    ! Point 7 = edge 1, spanned by vertices [1,2]
    coords_nby2(  7,:) = (vertex_coords( 1,:) + vertex_coords( 2,:)) / 2._dp
    ! Point 8 = edge 2, spanned by vertices [0,2]
    coords_nby2(  8,:) = (vertex_coords( 0,:) + vertex_coords( 2,:)) / 2._dp
    ! Point 9 = edge 3, spanned by vertices [2,3]
    coords_nby2(  9,:) = (vertex_coords( 2,:) + vertex_coords( 3,:)) / 2._dp
    ! Point 10 = edge 4, spanned by vertices [1,3]
    coords_nby2( 10,:) = (vertex_coords( 1,:) + vertex_coords( 3,:)) / 2._dp

    ! Reshape from n-by-2 to 2n
    n = 22
    allocate( coords_2n( 0:n-1))
    do i = 0, 10
      x = coords_nby2( i,1)
      y = coords_nby2( i,2)
      coords_2n( 2*i  ) = x
      coords_2n( 2*i+1) = y
    end do

    ! Convert to PETSc vector
    call vec_double2petsc( coords_2n, coords)

    ! Set DMPLEX point coordinates from the PETSc vector
    call DMSetCoordinates( dm, coords, perr)

    ! Clean up after yourself
    call DMDestroy( dm, perr)

  end subroutine create_simple_DMPLEX

end module ct_PETSc_DMPLEX_basics