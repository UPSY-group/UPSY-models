submodule(petsc_basic) petsc_matrix_vector_multiplication

contains

  subroutine multiply_PETSc_matrix_with_vector_1D( A, xx, yy)
    ! Multiply a PETSc matrix with a FORTRAN vector: y = A*x

    ! In- and output variables:
    type(tMat),                     intent(in   ) :: A
    real(dp), dimension(:), target, intent(in   ) :: xx
    real(dp), dimension(:), target, intent(  out) :: yy

    ! Local variables:
    character(len=*), parameter :: routine_name = 'multiply_PETSc_matrix_with_vector_1D'
    integer                     :: ierr
    integer                     :: nFx_loc
    integer, dimension(1:par%n) :: nFx_loc_all
    integer                     :: nFx_glob, i1Fx_glob, i2Fx_glob
    integer                     :: nFy_loc
    integer, dimension(1:par%n) :: nFy_loc_all
    integer                     :: nFy_glob, i1Fy_glob, i2Fy_glob
    integer                     :: mA_glob, nA_glob, mA_loc, nA_loc
    type(tVec)                  :: x, y

    ! Add routine to path
    call init_routine( routine_name)

    ! == Determine local and global size and local ownership range of Fortran vectors

    ! == x

    ! Local size
    nFx_loc = size( xx,1)

    ! Global size
    call MPI_ALLGATHER( nFx_loc, 1, MPI_integer, nFx_loc_all, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    nFx_glob = sum( nFx_loc_all)

    ! Local ownership ranges
    i1Fx_glob = sum( nFx_loc_all( 1:par%i)) + 1
    i2Fx_glob = i1Fx_glob + nFx_loc - 1

    ! == y

    ! Local size
    nFy_loc = size( yy,1)

    ! Global size
    call MPI_ALLGATHER( nFy_loc, 1, MPI_integer, nFy_loc_all, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    nFy_glob = sum( nFy_loc_all)

    ! Local ownership ranges
    i1Fy_glob = sum( nFy_loc_all( 1:par%i)) + 1
    i2Fy_glob = i1Fy_glob + nFy_loc - 1

    ! == Determine local and global size and local ownership range of PETSc matrix

    call MatGetSize(      A, mA_glob, nA_glob, ierr)
    call MatGetLocalSize( A, mA_loc , nA_loc , ierr)

    ! Safety
    if (nA_loc /= nFx_loc .or. mA_loc /= nFy_loc) then
      call crash('matrix and vector sub-sizes dont match!')
    end if

    ! Convert Fortran array xx to PETSc vector x
    call vec_double2petsc( xx, x)

    ! Create parallel vector
    call VecCreate( PETSC_COMM_WORLD, y, ierr)
    call VecSetSizes( y, nFy_loc, nFy_glob, ierr)
    call VecSetFromOptions( y, ierr)

    ! Compute the matrix-vector product
    call MatMult( A, x, y, ierr)

    ! Convert PETSc vector y to Fortran array yy
    call vec_petsc2double( y, yy)

    ! Clean up after yourself
    call VecDestroy( x, ierr)
    call VecDestroy( y, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_PETSc_matrix_with_vector_1D

  subroutine multiply_PETSc_matrix_with_vector_2D( A, xx, yy)
    ! Multiply a PETSc matrix with a FORTRAN vector: y = A*x

    ! In- and output variables:
    type(tMat),               intent(in   ) :: A
    real(dp), dimension(:,:), intent(in   ) :: xx
    real(dp), dimension(:,:), intent(  out) :: yy

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'multiply_PETSc_matrix_with_vector_2D'
    real(dp), dimension(:), allocatable :: xx_1D, yy_1D
    integer                             :: k

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( xx,2) /= size( yy,2)) then
      call crash('vector sizes dont match!')
    end if

    ! Allocate shared memory
    allocate( xx_1D( size( xx,1)), source = 0._dp)
    allocate( yy_1D( size( yy,1)), source = 0._dp)

    ! Compute each column separately
    do k = 1, size( xx,2)

      ! Copy this column of x
      xx_1D = xx( :,k)

      ! Compute the matrix-vector product
      call multiply_PETSc_matrix_with_vector_1D( A, xx_1D, yy_1D)

      ! Copy the result back
      yy( :,k) = yy_1D

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine multiply_PETSc_matrix_with_vector_2D

end submodule petsc_matrix_vector_multiplication
