module ct_PETSc_DMPLEX_basics

#include <petsc/finclude/petscsys.h>
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: warning
  use mpi_basic, only: par
  use mesh_types, only: type_mesh
  use netcdf_io_main, only: open_existing_netcdf_file_for_reading, setup_mesh_from_file, &
    close_netcdf_file
  use petsc, only: PetscErrorF, tDM, tVec, tPetscViewer, tPetscSection, DMPlexCreate, DMSetDimension, &
    DMSetCoordinateDim, DMPlexSetChart, DMPlexSetConeSize, DMSetUp, DMPlexSetCone, &
    DMPlexSymmetrize, DMPlexStratify, DMGetCoordinateSection, DMGetCoordinateDM, &
    DMCreateLocalVector, DMSetCoordinatesLocal, PetscSectionSetChart, PetscSectionSetDof, &
    PetscSectionSetUp, VecSetValues, VecAssemblyBegin, VecAssemblyEnd, INSERT_VALUES, &
    VecDestroy, PETSC_COMM_WORLD, DMDestroy, PetscObjectSetName
  use petsc_dmplex, only: write_dmplex_to_hdf5, mesh_to_dmplex
  use mesh_memory, only: allocate_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5, initialise_dummy_mesh_9, &
    initialise_dummy_mesh_16
  use mesh_edges, only: construct_mesh_edges
  use string_module, only: strrep

  implicit none

  private

  public :: create_simple_DMPLEX, ct_convert_meshes_to_dmplex

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

    if (par%primary) write(0,*) '   Running PETSc DMPLEX test: create a super simple DMPLEX object and write it to an HDF5 file'
    if (par%primary) write(0,*) ''

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

  subroutine ct_convert_meshes_to_dmplex( foldername_output, test_mesh_filenames, test_grid_filenames)

    ! In/output variables:
    character(len=*),               intent(in) :: foldername_output
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames
    character(len=*), dimension(:), intent(in) :: test_grid_filenames

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'ct_convert_meshes_to_dmplex'
    integer                       :: i_mesh
    type(type_mesh), allocatable  :: mesh
    character(len=:), allocatable :: filename_mesh
    integer                       :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '   Running PETSc DMPLEX test: converting UFEMISM meshes to DMPLEX objects'
    if (par%primary) write(0,*) ''

    ! Dummy meshes
    allocate( mesh)
    call allocate_mesh_primary( mesh, 'test_mesh_dummy_5', 100, 200)
    call initialise_dummy_mesh_5( mesh, 0._dp, 1._dp, 0._dp, 1._dp)
    call construct_mesh_edges( mesh)
    call ct_convert_mesh_to_dmplex_and_write_to_hdf5( foldername_output, mesh)
    deallocate( mesh)

    allocate( mesh)
    call allocate_mesh_primary( mesh, 'test_mesh_dummy_9', 100, 200)
    call initialise_dummy_mesh_9( mesh, 0._dp, 1._dp, 0._dp, 1._dp)
    call construct_mesh_edges( mesh)
    call ct_convert_mesh_to_dmplex_and_write_to_hdf5( foldername_output, mesh)
    deallocate( mesh)

    allocate( mesh)
    call allocate_mesh_primary( mesh, 'test_mesh_dummy_16', 100, 200)
    call initialise_dummy_mesh_16( mesh, 0._dp, 1._dp, 0._dp, 1._dp)
    call construct_mesh_edges( mesh)
    call ct_convert_mesh_to_dmplex_and_write_to_hdf5( foldername_output, mesh)
    deallocate( mesh)

    ! Meshes read from files
    do i_mesh = 1, size( test_mesh_filenames)
      filename_mesh = test_mesh_filenames( i_mesh)
      allocate( mesh)
      call open_existing_netcdf_file_for_reading( filename_mesh, ncid)
      call setup_mesh_from_file( filename_mesh, ncid, mesh)
      call close_netcdf_file( ncid)
      call ct_convert_mesh_to_dmplex_and_write_to_hdf5 ( foldername_output, mesh)
      deallocate( mesh)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ct_convert_meshes_to_dmplex

  subroutine ct_convert_mesh_to_dmplex_and_write_to_hdf5( foldername_output, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_output
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'ct_convert_mesh_to_dmplex_and_write_to_hdf5'
    type(tDM)                     :: dm
    integer                       :: ierr
    character(len=:), allocatable :: mesh_name_cleaned
    character(len=:), allocatable :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '    Converting UFEMISM mesh ', trim( mesh%name), ' to DMPLEX object...'

    call mesh_to_dmplex( mesh, dm)

    mesh_name_cleaned = trim( mesh%name)
    mesh_name_cleaned = strrep( mesh_name_cleaned, '"', '')
    mesh_name_cleaned = strrep( mesh_name_cleaned, '.', '_')
    mesh_name_cleaned = strrep( mesh_name_cleaned, '/', '_')
    filename = trim( foldername_output) // '/' // trim( mesh_name_cleaned) // '_dmplex.h5'
    call write_dmplex_to_hdf5( dm, filename)

    ! Clean up after yourself
    PetscCall( DMDestroy( dm, ierr))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ct_convert_mesh_to_dmplex_and_write_to_hdf5

end module ct_PETSc_DMPLEX_basics