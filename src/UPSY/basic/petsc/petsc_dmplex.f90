submodule(petsc_basic) petsc_dmplex

#include <petsc/finclude/petscsys.h>

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

end submodule petsc_dmplex
