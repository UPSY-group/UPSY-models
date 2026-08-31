module ct_PETSc_DMPLEX_basics

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: warning
  use petsc, only: tDM, tVec, tPetscViewer, tPetscSection, DMPlexCreate, DMSetDimension, &
    DMSetCoordinateDim, DMPlexSetChart, DMPlexSetConeSize, DMSetUp, DMPlexSetCone, &
    DMPlexSymmetrize, DMPlexStratify, DMGetCoordinateSection, DMGetCoordinateDM, &
    DMCreateLocalVector, DMSetCoordinatesLocal, PetscSectionSetChart, PetscSectionSetDof, &
    PetscSectionSetUp, VecSetValues, VecAssemblyBegin, VecAssemblyEnd, INSERT_VALUES, &
    VecDestroy, DMPlexTopologyView, &
    DMPlexCoordinatesView, PetscViewerCreate, PetscViewerSetType, &
    PETSCVIEWERASCII, PETSCVIEWERHDF5, PetscViewerFileSetMode, FILE_MODE_WRITE, &
    PetscViewerFileSetName, PetscViewerPushFormat, PetscViewerPopFormat, &
    PETSC_VIEWER_ASCII_INFO_DETAIL, DMView, PetscViewerDestroy, PETSC_COMM_WORLD, &
    DMDestroy, PetscObjectSetName, PETSC_VIEWER_HDF5_PETSC

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
    integer                             :: perr
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
    call DMPlexCreate( PETSC_COMM_WORLD, dm, perr)
    call PetscObjectSetName( dm, 'demo_dmplex', perr)
    call DMSetDimension( dm, 2, perr)
    call DMSetCoordinateDim( dm, 2, perr)

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
    call DMGetCoordinateSection( dm, coordinate_section, perr)
    call PetscSectionSetChart( coordinate_section, 0, 11, perr)
    do i = 2, 5
      call PetscSectionSetDof( coordinate_section, i, 2, perr)
    end do
    call PetscSectionSetUp( coordinate_section, perr)

    ! Create and fill the coordinate DM's local vector.
    call DMGetCoordinateDM( dm, coordinate_dm, perr)
    call DMCreateLocalVector( coordinate_dm, coords, perr)
    do i = 0, 7
      coords_indices( i) = i
    end do
    call VecSetValues( coords, n, coords_indices, coords_2n, INSERT_VALUES, perr)
    call VecAssemblyBegin( coords, perr)
    call VecAssemblyEnd( coords, perr)
    call DMSetCoordinatesLocal( dm, coords, perr)
    call VecDestroy( coords, perr)

    filename = trim( foldername_output) // '/PETSc_DMPLEX_output.h5'
    call write_dmplex_to_hdf5( dm, filename)

    ! Clean up after yourself
    call DMDestroy( dm, perr)

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
    integer                     :: perr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Export DMPLEX to an HDF5 file with the actual topology and coordinate arrays,
    ! not just the high-level metadata summary.

    ! Open the HDF5 file
    call PetscViewerCreate( PETSC_COMM_WORLD, viewer, perr)
    call PetscViewerSetType( viewer, PETSCVIEWERHDF5, perr)
    call PetscViewerFileSetMode( viewer, FILE_MODE_WRITE, perr)
    call PetscViewerFileSetName( viewer, trim( filename), perr)

    ! This format is important for saving a DMPlex
    call PetscViewerPushFormat( viewer, PETSC_VIEWER_HDF5_PETSC, perr)

    ! Write topology + coordinates + labels
    call DMView( dm, viewer, perr)

    call PetscViewerPopFormat( viewer, perr)
    call PetscViewerDestroy( viewer, perr)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_dmplex_to_hdf5

end module ct_PETSc_DMPLEX_basics