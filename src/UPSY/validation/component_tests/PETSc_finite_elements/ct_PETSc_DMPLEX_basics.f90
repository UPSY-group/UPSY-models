module ct_PETSc_DMPLEX_basics

#include <petsc/finclude/petscsys.h>
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: warning
  use petsc, only: tDM, tVec, tPetscViewer, tPetscSection, DMPlexCreate, DMSetDimension, &
    DMSetCoordinateDim, DMPlexSetChart, DMPlexSetConeSize, DMSetUp, DMPlexSetCone, &
    DMPlexSymmetrize, DMPlexStratify, DMGetCoordinateSection, DMGetCoordinateDM, &
    DMCreateLocalVector, DMSetCoordinatesLocal, PetscSectionSetChart, PetscSectionSetDof, &
    PetscSectionSetUp, VecSetValues, VecAssemblyBegin, VecAssemblyEnd, INSERT_VALUES, &
    VecDestroy, PetscViewerCreate, PetscViewerSetType, &
    PETSCVIEWERHDF5, PetscViewerFileSetMode, FILE_MODE_WRITE, &
    PetscViewerFileSetName, PetscViewerPushFormat, PetscViewerPopFormat, &
    DMView, PetscViewerDestroy, PETSC_COMM_WORLD, &
    DMDestroy, PetscObjectSetName, PETSC_VIEWER_HDF5_PETSC, PetscErrorF

  implicit none

  private

  public :: create_simple_DMPLEX

contains

  subroutine create_simple_DMPLEX( foldername_output)
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

    ! In/output variables:
    character(len=*), intent(in) :: foldername_output

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'create_simple_DMPLEX'
    type(tDM)                           :: dm
    integer                             :: ierr
    real(dp), dimension(0:3,2)          :: vertex_coords
    real(dp), dimension(0:3,2)          :: coords_nby2
    real(dp), dimension(:), allocatable :: coords_2n
    integer                             :: n,i
    real(dp)                            :: x,y
    type(tVec)                          :: coords
    type(tDM)                           :: coordinate_dm
    type(tPetscSection)                 :: coordinate_section
    integer, dimension(0:7)             :: coords_indices
    character(len=:), allocatable       :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Create a DMPlex object
    PetscCall( DMPlexCreate( PETSC_COMM_WORLD, dm, ierr))
    PetscCall( PetscObjectSetName( dm, 'demo_dmplex', ierr))
    PetscCall( DMSetDimension( dm, 2, ierr))
    PetscCall( DMSetCoordinateDim( dm, 2, ierr))

    ! The example mesh has 11 'points' (vertices, faces, and edges)
    PetscCall( DMPlexSetChart( dm, 0, 11, ierr))

    ! In 2 dimensions the convention is to first number faces, then vertices, and then edges.
    ! PETSc indexes from zero
    !    DMPlexSetConeSize( dm, point, number of points that cover the point)
    PetscCall( DMPlexSetConeSize( dm,  0, 3, ierr))  ! Points 0 and 1 are faces, each connected to 3 edges.
    PetscCall( DMPlexSetConeSize( dm,  1, 3, ierr))
    PetscCall( DMPlexSetConeSize( dm,  6, 2, ierr))  ! Points 6-10 are edges, each connected to 2 vertices.
    PetscCall( DMPlexSetConeSize( dm,  7, 2, ierr))
    PetscCall( DMPlexSetConeSize( dm,  8, 2, ierr))
    PetscCall( DMPlexSetConeSize( dm,  9, 2, ierr))
    PetscCall( DMPlexSetConeSize( dm, 10, 2, ierr))
    PetscCall( DMSetUp( dm, ierr))

    !    DMPlexSetCone( dm, point, [points that cover the point])
    PetscCall( DMPlexSetCone( dm,  0, [6, 7,  8], ierr))   ! Point  0 (= face 0) is connected to points [6,7, 8] (= edges [0,1,2])
    PetscCall( DMPlexSetCone( dm,  1, [7, 9, 10], ierr))   ! Point  1 (= face 1) is connected to points [7,9,10] (= edges [1,3,4])  ...so not necessarily counter-clockwise
    PetscCall( DMPlexSetCone( dm,  6, [2, 3    ], ierr))   ! Point  6 (= edge 0) is connected to points [2,3] (= vertices [0,1])
    PetscCall( DMPlexSetCone( dm,  7, [3, 4    ], ierr))   ! Point  7 (= edge 1) is connected to points [3,4] (= vertices [1,2])
    PetscCall( DMPlexSetCone( dm,  8, [4, 2    ], ierr))   ! Point  8 (= edge 2) is connected to points [4,2] (= vertices [2,0])
    PetscCall( DMPlexSetCone( dm,  9, [4, 5    ], ierr))   ! Point  9 (= edge 3) is connected to points [4,5] (= vertices [2,3])
    PetscCall( DMPlexSetCone( dm, 10, [5, 3    ], ierr))   ! Point 10 (= edge 4) is connected to points [5,3] (= vertices [3,1])

    ! Let PETSc automatically figure out the 'supports', i.e. the backward connections (so each
    ! vertex knows which edges and faces it spans)
    PetscCall( DMPlexSymmetrize( dm, ierr))

    ! In order to support efficient queries, we construct fast search structures
    ! and indices for the different types of points
    PetscCall( DMPlexStratify( dm, ierr))

    ! ADDITION: set coordinates

    ! Basic vertex coordinates
    vertex_coords( 0,:) = [-1._dp,  0._dp]
    vertex_coords( 1,:) = [ 0._dp, -1._dp]
    vertex_coords( 2,:) = [ 0._dp,  1._dp]
    vertex_coords( 3,:) = [ 1._dp,  0._dp]

    ! DMPLEX coordinates are stored only for the vertices of the mesh.
    ! The topology points for faces and edges do not have coordinate entries.
    coords_nby2( 0,:) = vertex_coords( 0,:)
    coords_nby2( 1,:) = vertex_coords( 1,:)
    coords_nby2( 2,:) = vertex_coords( 2,:)
    coords_nby2( 3,:) = vertex_coords( 3,:)

    ! Reshape from n-by-2 to 2n
    n = 8
    allocate( coords_2n( 0:n-1))
    do i = 0, 3
      x = coords_nby2( i,1)
      y = coords_nby2( i,2)
      coords_2n( 2*i  ) = x
      coords_2n( 2*i+1) = y
    end do

    ! Define two coordinate degrees of freedom for each vertex. The coordinate
    ! vector must use this DMPlex-owned layout rather than a general parallel Vec.
    PetscCall( DMGetCoordinateSection( dm, coordinate_section, ierr))
    PetscCall( PetscSectionSetChart( coordinate_section, 0, 11, ierr))
    do i = 2, 5
      PetscCall( PetscSectionSetDof( coordinate_section, i, 2, ierr))
    end do
    PetscCall( PetscSectionSetUp( coordinate_section, ierr))

    ! Create and fill the coordinate DM's local vector.
    PetscCall( DMGetCoordinateDM( dm, coordinate_dm, ierr))
    PetscCall( DMCreateLocalVector( coordinate_dm, coords, ierr))
    do i = 0, 7
      coords_indices( i) = i
    end do
    PetscCall( VecSetValues( coords, n, coords_indices, coords_2n, INSERT_VALUES, ierr))
    PetscCall( VecAssemblyBegin( coords, ierr))
    PetscCall( VecAssemblyEnd( coords, ierr))
    PetscCall( DMSetCoordinatesLocal( dm, coords, ierr))
    PetscCall( VecDestroy( coords, ierr))

    filename = trim( foldername_output) // '/PETSc_DMPLEX_output.h5'
    call write_dmplex_to_hdf5( dm, filename)

    ! Clean up after yourself
    PetscCall( DMDestroy( dm, ierr))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_simple_DMPLEX

  subroutine write_dmplex_to_hdf5( dm, filename)

    ! In/output variables:
    type(tDM),        intent(in) :: dm
    character(len=*), intent(in) :: filename

    ! Local variables:
    character(len=*), parameter :: routine_name = 'write_dmplex_to_hdf5'
    type(tPetscViewer)          :: viewer
    integer                     :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Export DMPLEX to an HDF5 file with the actual topology and coordinate arrays,
    ! not just the high-level metadata summary.

    ! Open the HDF5 file
    PetscCall( PetscViewerCreate( PETSC_COMM_WORLD, viewer, ierr))
    PetscCall( PetscViewerSetType( viewer, PETSCVIEWERHDF5, ierr))
    PetscCall( PetscViewerFileSetMode( viewer, FILE_MODE_WRITE, ierr))
    PetscCall( PetscViewerFileSetName( viewer, trim( filename), ierr))

    ! This format is important for saving a DMPlex
    PetscCall( PetscViewerPushFormat( viewer, PETSC_VIEWER_HDF5_PETSC, ierr))

    ! Write topology + coordinates + labels
    PetscCall( DMView( dm, viewer, ierr))

    PetscCall( PetscViewerPopFormat( viewer, ierr))
    PetscCall( PetscViewerDestroy( viewer, ierr))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_dmplex_to_hdf5

end module ct_PETSc_DMPLEX_basics