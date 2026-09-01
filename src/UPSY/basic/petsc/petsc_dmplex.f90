module petsc_dmplex

#include <petsc/finclude/petscsys.h>
  use precisions, only: dp
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use petsc, only: PetscErrorF, PETSC_COMM_WORLD, tDM, tPetscViewer, PetscViewerCreate, &
    PetscViewerSetType, PETSCVIEWERHDF5, PetscViewerFileSetMode, FILE_MODE_WRITE, &
    PetscViewerFileSetName, PetscViewerPushFormat, PETSC_VIEWER_HDF5_PETSC, DMView, &
    PetscViewerPopFormat, PetscViewerDestroy
  use assertions_basic, only: assert
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use mpi_f08, only: MPI_ALLGATHER, MPI_INTEGER, MPI_COMM_WORLD
  use string_module, only: colour_string

  implicit none

  private

  public :: write_dmplex_to_hdf5

contains

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
