MODULE petsc_basic

  ! Contains routines that use the PETSc matrix solvers
  !
  ! Convention: xx = Fortran, x = PETSc

#include <petsc/finclude/petscksp.h>
  USE petscksp
  use assertions_basic
  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par
  USE call_stack_and_comp_time_tracking                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE parameters
  USE reallocate_mod                                         , ONLY: reallocate
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use mpi_distributed_memory, only: partition_list, gather_to_all
  use petsc_vectors, only: vec_double2petsc, vec_petsc2double

  IMPLICIT NONE

  private

  public :: multiply_PETSc_matrix_with_vector_1D, multiply_PETSc_matrix_with_vector_2D

  INTEGER :: perr    ! Error flag for PETSc routines
  character(len=:), allocatable :: PETSc_KSPtype_printed, PETSc_PCtype_printed

CONTAINS

! == Matrix-vector multiplication

  SUBROUTINE multiply_PETSc_matrix_with_vector_1D( A, xx, yy)
    ! Multiply a PETSc matrix with a FORTRAN vector: y = A*x

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(tMat),                          INTENT(IN)    :: A
    REAL(dp), DIMENSION(:    ), TARGET,  INTENT(IN)    :: xx
    REAL(dp), DIMENSION(:    ), TARGET,  INTENT(OUT)   :: yy

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'multiply_PETSc_matrix_with_vector_1D'
    integer                                            :: ierr
    TYPE(PetscInt)                                     :: nFx_loc
    INTEGER,  DIMENSION(1:par%n)                       :: nFx_loc_all
    TYPE(PetscInt)                                     :: nFx_glob, i1Fx_glob, i2Fx_glob
    TYPE(PetscInt)                                     :: nFy_loc
    INTEGER,  DIMENSION(1:par%n)                       :: nFy_loc_all
    TYPE(PetscInt)                                     :: nFy_glob, i1Fy_glob, i2Fy_glob
    TYPE(PetscInt)                                     :: mA_glob, nA_glob, mA_loc, nA_loc
    TYPE(tVec)                                         :: x, y

    ! Add routine to path
    CALL init_routine( routine_name)

    ! == Determine local and global size and local ownership range of Fortran vectors

    ! == x

    ! Local size
    nFx_loc = SIZE( xx,1)

    ! Global size
    CALL MPI_ALLGATHER( nFx_loc, 1, MPI_INTEGER, nFx_loc_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    nFx_glob = SUM( nFx_loc_all)

    ! Local ownership ranges
    i1Fx_glob = SUM( nFx_loc_all( 1:par%i)) + 1
    i2Fx_glob = i1Fx_glob + nFx_loc - 1

    ! == y

    ! Local size
    nFy_loc = SIZE( yy,1)

    ! Global size
    CALL MPI_ALLGATHER( nFy_loc, 1, MPI_INTEGER, nFy_loc_all, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    nFy_glob = SUM( nFy_loc_all)

    ! Local ownership ranges
    i1Fy_glob = SUM( nFy_loc_all( 1:par%i)) + 1
    i2Fy_glob = i1Fy_glob + nFy_loc - 1

    ! == Determine local and global size and local ownership range of PETSc matrix

    CALL MatGetSize(      A, mA_glob, nA_glob, perr)
    CALL MatGetLocalSize( A, mA_loc , nA_loc , perr)

    ! Safety
    IF (nA_loc /= nFx_loc .OR. mA_loc /= nFy_loc) THEN
      CALL crash('matrix and vector sub-sizes dont match!')
    END IF

    ! Convert Fortran array xx to PETSc vector x
    CALL vec_double2petsc( xx, x)

    ! Create parallel vector
    CALL VecCreate( PETSC_COMM_WORLD, y, perr)
    CALL VecSetSizes( y, nFy_loc, nFy_glob, perr)
    CALL VecSetFromOptions( y, perr)

    ! Compute the matrix-vector product
    CALL MatMult( A, x, y, perr)

    ! Convert PETSc vector y to Fortran array yy
    CALL vec_petsc2double( y, yy)

    ! Clean up after yourself
    CALL VecDestroy( x, perr)
    CALL VecDestroy( y, perr)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE multiply_PETSc_matrix_with_vector_1D

  SUBROUTINE multiply_PETSc_matrix_with_vector_2D( A, xx, yy)
    ! Multiply a PETSc matrix with a FORTRAN vector: y = A*x

    IMPLICIT NONE

    ! In- and output variables:
    TYPE(tMat),                          INTENT(IN)    :: A
    REAL(dp), DIMENSION(:,:  ),          INTENT(IN)    :: xx
    REAL(dp), DIMENSION(:,:  ),          INTENT(OUT)   :: yy

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'multiply_PETSc_matrix_with_vector_2D'
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: xx_1D, yy_1D
    INTEGER                                            :: k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF (SIZE( xx,2) /= SIZE( yy,2)) THEN
      CALL crash('vector sizes dont match!')
    END IF

    ! Allocate shared memory
    ALLOCATE( xx_1D( SIZE( xx,1)), source = 0._dp)
    ALLOCATE( yy_1D( SIZE( yy,1)), source = 0._dp)

    ! Compute each column separately
    DO k = 1, SIZE( xx,2)

      ! Copy this column of x
      xx_1D = xx( :,k)

      ! Compute the matrix-vector product
      CALL multiply_PETSc_matrix_with_vector_1D( A, xx_1D, yy_1D)

      ! Copy the result back
      yy( :,k) = yy_1D

    END DO

    ! Clean up after yourself
    DEALLOCATE( xx_1D)
    DEALLOCATE( yy_1D)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE multiply_PETSc_matrix_with_vector_2D

! ===== Unit tests =====
! ======================

END MODULE petsc_basic
