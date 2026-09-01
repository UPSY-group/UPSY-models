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

  public :: mat_CSR2petsc, mat_petsc2CSR, multiply_PETSc_matrix_with_vector_1D, multiply_PETSc_matrix_with_vector_2D

  INTEGER :: perr    ! Error flag for PETSc routines
  character(len=:), allocatable :: PETSc_KSPtype_printed, PETSc_PCtype_printed

CONTAINS

  SUBROUTINE mat_petsc2CSR( A, AA)
    ! Convert a PETSC parallel matrix to a CSR-format matrix in regular Fortran arrays

    IMPLICIT NONE

    ! In/output variables:
    TYPE(tMat),                          INTENT(IN)    :: A
    TYPE(type_CSR_matrix_dp),     INTENT(OUT)   :: AA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'mat_petsc2CSR'
    INTEGER                                            :: m_glob, n_glob, m_loc, n_loc, istart, iend, row_glob, row_loc
    INTEGER                                            :: ncols
    INTEGER,  DIMENSION(:    ), pointer                :: cols_
    REAL(dp), DIMENSION(:    ), pointer                :: vals_
    INTEGER,  DIMENSION(:    ), ALLOCATABLE, target    :: cols
    REAL(dp), DIMENSION(:    ), ALLOCATABLE, target    :: vals
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: nnz_row_loc
    INTEGER                                            :: nnz_loc
    INTEGER                                            :: k

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Retrieve global and local matrix size and ownership range
    CALL MatGetSize(           A, m_glob, n_glob, perr)
    CALL MatGetLocalSize(      A, m_loc , n_loc , perr)
    CALL MatGetOwnershipRange( A, istart, iend  , perr)

    ! Find number of non-zeros in each row
    ALLOCATE( nnz_row_loc( m_loc ))
    ALLOCATE( cols(        n_glob))
    ALLOCATE( vals(        n_glob))

    cols_ => cols
    vals_ => vals

    DO row_glob = istart+1, iend ! +1 because PETSc indexes from 0
      row_loc = row_glob - istart
      CALL MatGetRow( A, row_glob-1, ncols, cols_, vals_, perr)
      nnz_row_loc( row_loc) = ncols
      CALL MatRestoreRow( A, row_glob-1, ncols, cols_, vals_, perr)
    END DO

    nnz_loc = SUM( nnz_row_loc)

    ! Allocate memory for CSR matrix
    call AA%allocate( m_glob, n_glob, m_loc, n_loc, nnz_loc)

    ! Copy data from the PETSc matrix to the CSR arrays
    DO row_glob = istart+1, iend ! +1 because PETSc indexes from 0
      CALL MatGetRow( A, row_glob-1, ncols, cols_, vals_, perr)
      DO k = 1, ncols
        call AA%add_entry( row_glob, cols_( k)+1, vals_( k))
      END DO
      CALL MatRestoreRow( A, row_glob-1, ncols, cols_, vals_, perr)
    END DO

    ! Crop memory
    call AA%finalise

    ! Clean up after yourself
    DEALLOCATE( nnz_row_loc)
    DEALLOCATE( cols       )
    DEALLOCATE( vals       )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE mat_petsc2CSR

  SUBROUTINE mat_CSR2petsc( AA, A)
    ! Convert a CSR-format matrix in regular Fortran arrays to a PETSC parallel matrix
    !
    ! NOTE: the PETSc documentation seems to advise against using the MatCreateMPIAIJWithArrays
    !       routine used here. However, for the advised way of using MatSetValues with preallocation
    !       I've not been able to find a way that is fast enough to be useful without having to
    !       preallocate -WAY- too much memory. Especially for the remapping matrices, which
    !       can have hundreds or even thousands of non-zero elements per row, this can make the
    !       model run hella slow, whereas the current solution seems to work perfectly. So there you go.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_CSR_matrix_dp),     INTENT(IN)    :: AA
    TYPE(tMat),                          INTENT(OUT)   :: A

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'mat_CSR2petsc'
    INTEGER                                            :: i, k1, k2, nnz_proc, ii, k, kk
    INTEGER,  DIMENSION(:    ), ALLOCATABLE            :: ptr_proc, ind_proc
    REAL(dp), DIMENSION(:    ), ALLOCATABLE            :: val_proc

    ! Add routine to path
    CALL init_routine( routine_name)

    if (.not. AA%is_finalised) call crash('A is not finalised')

    ! Determine number of non-zeros for this process
    nnz_proc = AA%nnz

    ! Allocate memory for local CSR-submatrix
    ALLOCATE( ptr_proc( 0:AA%m_loc      ))
    ALLOCATE( ind_proc( 0:nnz_proc   - 1))
    ALLOCATE( val_proc( 0:nnz_proc   - 1))

    ! Copy matrix data
    DO i = AA%i1, AA%i2

      ! ptr
      ii = i - AA%i1
      ptr_proc( ii) = AA%ptr( i) - AA%ptr( AA%i1)

      ! index and val
      k1 = AA%ptr( i)
      k2 = AA%ptr( i+1) - 1
      DO k = k1, k2
        kk = k - AA%ptr( AA%i1)
        ind_proc( kk) = AA%ind( k) - 1
        val_proc( kk) = AA%val( k)
      END DO

    END DO
    ! Last row
    ptr_proc( AA%m_loc) = AA%ptr( AA%i2+1) - AA%ptr( AA%i1)

    ! Create PETSc matrix
    ! IF (AA%balanced) then
      CALL MatCreateMPIAIJWithArrays( PETSC_COMM_WORLD, AA%m_loc, AA%n_loc, AA%m, AA%n, ptr_proc, ind_proc, val_proc, A, perr)
    ! ELSE
    !   ! Special treatment if the rows are not partitioned according to PETSC
    !   CALL MatCreateMPIAIJWithArrays( PETSC_COMM_WORLD, nrows_proc, nrows_proc, AA%m, AA%n, ptr_proc, ind_proc, val_proc, A, perr)
    ! END IF

    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.

    CALL MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY, perr)
    CALL MatAssemblyEnd(   A, MAT_FINAL_ASSEMBLY, perr)

    ! Clean up after yourself
    DEALLOCATE( ptr_proc)
    DEALLOCATE( ind_proc)
    DEALLOCATE( val_proc)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE mat_CSR2petsc

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
