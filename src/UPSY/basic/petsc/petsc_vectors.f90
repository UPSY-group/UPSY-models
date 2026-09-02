submodule(petsc_basic) petsc_vectors

#include <petsc/finclude/petscsys.h>

contains

  subroutine vec_double2petsc( xx, x)
    ! Convert a regular 1-D Fortran double-precision array to a PETSc parallel vector

    ! In- and output variables:
    real(dp), dimension(:), intent(in)  :: xx
    type(tVec),             intent(out) :: x

    ! Local variables:
    character(len=*), parameter  :: routine_name = 'vec_double2petsc'
    integer                      :: ierr
    integer                      :: nF_loc
    integer, dimension(1:par%n)  :: nF_loc_all
    integer                      :: nF_glob, i1F_glob, i2F_glob, iF_loc
    integer, dimension(size(xx)) :: iP_glob

    ! Add routine to path
    call init_routine( routine_name)

    ! == Determine local and global size and local ownership range of Fortran vector

    ! Local size
    nF_loc = size( xx,1)

    ! Global size
    call MPI_ALLGATHER( nF_loc, 1, MPI_INTEGER, nF_loc_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    nF_glob = sum( nF_loc_all)

    ! Local ownership ranges
    i1F_glob = sum( nF_loc_all( 1:par%i)) + 1
    i2F_glob = i1F_glob + nF_loc - 1

    ! Create parallel vector
    PetscCall( VecCreate( PETSC_COMM_WORLD, x, ierr))
    PetscCall( VecSetSizes( x, nF_loc, nF_glob, ierr))
    PetscCall( VecSetFromOptions( x, ierr))

    ! Determine global PETSc indices
    do iF_loc = 1, nF_loc
      iP_glob( iF_loc) = i1F_glob - 2 + iF_loc
    end do

    ! Set PETSc vector values
    PetscCall( VecSetValues( x, nF_loc, iP_glob, xx, INSERT_VALUES, ierr))

    ! Assemble vectors, using the 2-step process:
    !   VecAssemblyBegin(), VecAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.

    PetscCall( VecAssemblyBegin( x, ierr))
    PetscCall( VecAssemblyEnd(   x, ierr))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine vec_double2petsc

  subroutine vec_petsc2double( x, xx)
    ! Convert a PETSc parallel vector to a regular 1-D Fortran double-precision array

    ! In- and output variables:
    type(tVec),             intent(in)  :: x
    real(dp), dimension(:), intent(out) :: xx

    ! Local variables:
    character(len=*), parameter  :: routine_name = 'vec_petsc2double'
    integer                      :: ierr
    integer                      :: nF_loc
    integer, dimension(1:par%n)  :: nF_loc_all
    integer                      :: nF_glob, i1F_glob, i2F_glob, iF_loc
    integer, dimension(size(xx)) :: iP_glob
    integer                      :: nP_loc, nP_glob, i1P_glob, i2P_glob

    ! Add routine to path
    call init_routine( routine_name)

    ! == Determine local and global size and local ownership range of Fortran vector

    ! Local size
    nF_loc = size( xx,1)

    ! Global size
    call MPI_ALLGATHER( nF_loc, 1, MPI_INTEGER, nF_loc_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    nF_glob = sum( nF_loc_all)

    ! Local ownership ranges
    i1F_glob = sum( nF_loc_all( 1:par%i)) + 1
    i2F_glob = i1F_glob + nF_loc - 1

    ! == Determine local and global size and local ownership range of PETSc vector

    ! Local size
    PetscCall( VecGetLocalSize( x, nP_loc, ierr))

    ! Global size
    PetscCall( VecGetSize( x, nP_glob, ierr))

    ! Safety
    PetscCall( VecGetOwnershipRange( x, i1P_glob, i2P_glob, ierr))

    ! Safety
    call assert( nF_loc == nP_loc, 'nF_loc /= nP_loc')
    call assert( nF_glob == nP_glob, 'nF_glob /= nP_glob')
    call assert( i1F_glob == i1P_glob+1, 'i1F_glob /= i1P_glob+1')
    call assert( i2F_glob == i2P_glob, 'i2F_glob /= i2P_glob')

    ! Determine global PETSc indices
    do iF_loc = 1, nF_loc
      iP_glob( iF_loc) = i1F_glob - 2 + iF_loc
    end do

    ! Get PETSc vector values
    PetscCall( VecGetValues( x, nF_loc, iP_glob, xx, ierr))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine vec_petsc2double

end submodule petsc_vectors
