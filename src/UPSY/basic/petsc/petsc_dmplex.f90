module petsc_dmplex

#include <petsc/finclude/petscsys.h>
  use precisions, only: dp
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use petsc, only: PetscErrorF, PETSC_COMM_WORLD, tDM, tPetscViewer, PetscViewerCreate, &
    PetscViewerSetType, PETSCVIEWERHDF5, PetscViewerFileSetMode, FILE_MODE_WRITE, &
    PetscViewerFileSetName, PetscViewerPushFormat, PETSC_VIEWER_HDF5_PETSC, DMView, &
    PetscViewerPopFormat, PetscViewerDestroy, tVec, tPetscSection, DMPlexCreate, &
    PetscObjectSetName, DMSetDimension, DMSetCoordinateDim, DMPlexSetChart, &
    DMPlexSetConeSize, DMSetUp, DMPlexSetCone, DMPlexSymmetrize, DMPlexStratify, &
    DMGetCoordinateSection, PetscSectionSetChart, PetscSectionSetDof, PetscSectionSetUp, &
    DMGetCoordinateDM, DMCreateLocalVector, VecSetValues, INSERT_VALUES, VecAssemblyBegin, &
    VecAssemblyEnd, DMSetCoordinatesLocal, VecDestroy
  use assertions_basic, only: assert
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use mpi_f08, only: MPI_ALLGATHER, MPI_INTEGER, MPI_COMM_WORLD
  use string_module, only: colour_string
  use mesh_types, only: type_mesh

  implicit none

  private

  public :: mesh_to_dmplex, write_dmplex_to_hdf5

contains

  subroutine mesh_to_dmplex( mesh, dm)

    ! In/output variables:
    type(type_mesh),  intent(in   ) :: mesh
    type(tDM),        intent(  out) :: dm

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'mesh_to_dmplex'
    integer                             :: ierr
    integer,  dimension(:), allocatable :: vi2p, p2vi
    integer,  dimension(:), allocatable :: ti2p, p2ti
    integer,  dimension(:), allocatable :: ei2p, p2ei
    integer                             :: np, n, vi, p, i
    real(dp), dimension(:), allocatable :: coords_2n
    type(tVec)                          :: coords
    type(tDM)                           :: coordinate_dm
    type(tPetscSection)                 :: coordinate_section
    integer,  dimension(:), allocatable :: coords_indices

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Create a DMPlex object
    PetscCall( DMPlexCreate( PETSC_COMM_WORLD, dm, ierr))
    PetscCall( PetscObjectSetName( dm, 'dmplex_' // trim( mesh%name), ierr))
    PetscCall( DMSetDimension( dm, 2, ierr))
    PetscCall( DMSetCoordinateDim( dm, 2, ierr))

    call calc_vertex_triangle_edge_point_translation_tables( mesh, np, vi2p, p2vi, ti2p, p2ti, ei2p, p2ei)
    call set_dmplex_topology( mesh, dm, np, vi2p, p2vi, ti2p, p2ti, ei2p, p2ei)

    ! Let PETSc automatically figure out the 'supports', i.e. the backward connections (so each
    ! vertex knows which edges and faces it spans)
    PetscCall( DMPlexSymmetrize( dm, ierr))

    ! In order to support efficient queries, we construct fast search structures
    ! and indices for the different types of points
    PetscCall( DMPlexStratify( dm, ierr))

    ! Set vertex coordinates

    ! Reshape from nV-by-2 to 2*nV
    n = mesh%nV * 2
    allocate( coords_2n( 0:n-1))
    do vi = 1, mesh%nV
      coords_2n( 2*vi-2) = mesh%V( vi,1)
      coords_2n( 2*vi-1) = mesh%V( vi,2)
    end do

    ! Define two coordinate degrees of freedom for each vertex. The coordinate
    ! vector must use this DMPlex-owned layout rather than a general parallel Vec.
    PetscCall( DMGetCoordinateSection( dm, coordinate_section, ierr))
    PetscCall( PetscSectionSetChart( coordinate_section, 0, np, ierr))
    do vi = 1, mesh%nV
      p = vi2p( vi)
      PetscCall( PetscSectionSetDof( coordinate_section, p, 2, ierr))
    end do
    PetscCall( PetscSectionSetUp( coordinate_section, ierr))

    ! Create and fill the coordinate DM's local vector.
    PetscCall( DMGetCoordinateDM( dm, coordinate_dm, ierr))
    PetscCall( DMCreateLocalVector( coordinate_dm, coords, ierr))
    allocate( coords_indices( 0:n-1))
    do i = 0, n-1
      coords_indices( i) = i
    end do
    PetscCall( VecSetValues( coords, n, coords_indices, coords_2n, INSERT_VALUES, ierr))
    PetscCall( VecAssemblyBegin( coords, ierr))
    PetscCall( VecAssemblyEnd( coords, ierr))
    PetscCall( DMSetCoordinatesLocal( dm, coords, ierr))
    PetscCall( VecDestroy( coords, ierr))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine mesh_to_dmplex

  subroutine calc_vertex_triangle_edge_point_translation_tables( mesh, np, vi2p, p2vi, ti2p, p2ti, ei2p, p2ei)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    integer,                            intent(  out) :: np
    integer, dimension(:), allocatable, intent(inout) :: vi2p, p2vi, ti2p, p2ti, ei2p, p2ei

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_vertex_triangle_edge_point_translation_tables'
    integer                     :: vi, ti, ei, p

    ! Add routine to call stack
    call init_routine( routine_name)

    ! The number of 'points' is the combined number of vertices, faces, and edges
    np = mesh%nV + mesh%nTri + mesh%nE

    allocate( vi2p( 1:mesh%nV  ), source = -1)
    allocate( ti2p( 1:mesh%nTri), source = -1)
    allocate( ei2p( 1:mesh%nE  ), source = -1)

    allocate( p2vi( 0:np-1), source = -1)
    allocate( p2ti( 0:np-1), source = -1)
    allocate( p2ei( 0:np-1), source = -1)

    ! In 2 dimensions the convention is to first number faces, then vertices, and then edges.
    p = -1

    do vi = 1, mesh%nV
      p = p+1
      vi2p( vi) = p
      p2vi( p ) = vi
    end do

    do ti = 1, mesh%nTri
      p = p+1
      ti2p( ti) = p
      p2ti( p ) = ti
    end do

    do ei = 1, mesh%nE
      p = p+1
      ei2p( ei) = p
      p2ei( p ) = ei
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calc_vertex_triangle_edge_point_translation_tables

  subroutine set_dmplex_topology( mesh, dm, np, vi2p, p2vi, ti2p, p2ti, ei2p, p2ei)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    type(tDM),                       intent(inout) :: dm
    integer,                         intent(in   ) :: np
    integer, dimension(1:mesh%nV  ), intent(in   ) :: vi2p
    integer, dimension(0:np-1     ), intent(in   ) :: p2vi
    integer, dimension(1:mesh%nTri), intent(in   ) :: ti2p
    integer, dimension(0:np-1     ), intent(in   ) :: p2ti
    integer, dimension(1:mesh%nE  ), intent(in   ) :: ei2p
    integer, dimension(0:np-1     ), intent(in   ) :: p2ei

    ! Local variables:
    character(len=*), parameter :: routine_name = 'set_dmplex_topology'
    integer                     :: ierr
    integer                     :: vi, ti, ei, p
    integer, dimension(3)       :: cone_triangle
    integer, dimension(2)       :: cone_edge

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Set the total number of points
    PetscCall( DMPlexSetChart( dm, 0, np, ierr))

    ! First set the 'cone size', i.e. to how many lower-level points
    ! each higher-level point is connected

    ! In 2 dimensions the convention is to first number faces, then vertices, and then edges.

    ! Faces = mesh triangles, which are each connected to 3 edges
    do ti = 1, mesh%nTri
      p = ti2p( ti)
      ! DMPlexSetConeSize( dm, point, number of points that cover the point)
      PetscCall( DMPlexSetConeSize( dm,  p, 3, ierr))
    end do

    ! Edges = mesh edges, which are each connected to 2 vertices
    do ei = 1, mesh%nE
      p = ei2p( ei)
      ! DMPlexSetConeSize( dm, point, number of points that cover the point)
      PetscCall( DMPlexSetConeSize( dm,  p, 2, ierr))
    end do

    ! Finish setting up the cone sizes
    PetscCall( DMSetUp( dm, ierr))

    ! Then, set the actual cones

    ! Triangle-edge connectivity
    do ti = 1, mesh%nTri
      p = ti2p( ti)
      cone_triangle = [ei2p( mesh%TriE( ti,1)), ei2p( mesh%TriE( ti,2)), ei2p( mesh%TriE( ti,3))]
      ! DMPlexSetCone( dm, point, [points that cover the point])
      PetscCall( DMPlexSetCone( dm, p, cone_triangle, ierr))
    end do

    ! Edge-vertex connectivity
    do ei = 1, mesh%nE
      p = ei2p( ei)
      cone_edge = [vi2p( mesh%EV( ei,1)), vi2p( mesh%EV( ei,2))]
      ! DMPlexSetCone( dm, point, [points that cover the point])
      PetscCall( DMPlexSetCone( dm, p, cone_edge, ierr))
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_dmplex_topology

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

end module petsc_dmplex
