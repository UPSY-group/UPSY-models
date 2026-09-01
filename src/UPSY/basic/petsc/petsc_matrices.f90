submodule(petsc_basic) petsc_matrices

contains

  subroutine mat_petsc2CSR( A, AA)
    ! Convert a PETSC parallel matrix to a CSR-format matrix

    ! In/output variables:
    type(tMat),               intent(in   ) :: A
    type(type_CSR_matrix_dp), intent(  out) :: AA

    ! Local variables:
    character(len=*), parameter                 :: routine_name = 'mat_petsc2CSR'
    integer                                     :: ierr
    integer                                     :: m_glob, n_glob, m_loc, n_loc, istart, iend, row_glob, row_loc
    integer                                     :: ncols
    integer,  dimension(:), pointer             :: cols_
    real(dp), dimension(:), pointer             :: vals_
    integer,  dimension(:), allocatable, target :: cols
    real(dp), dimension(:), allocatable, target :: vals
    integer,  dimension(:), allocatable         :: nnz_row_loc
    integer                                     :: nnz_loc
    integer                                     :: k

    ! Add routine to path
    call init_routine( routine_name)

    ! Retrieve global and local matrix size and ownership range
    call MatGetSize(           A, m_glob, n_glob, ierr)
    call MatGetLocalSize(      A, m_loc , n_loc , ierr)
    call MatGetOwnershipRange( A, istart, iend  , ierr)

    ! Find number of non-zeros in each row
    allocate( nnz_row_loc( m_loc ))
    allocate( cols(        n_glob))
    allocate( vals(        n_glob))

    cols_ => cols
    vals_ => vals

    do row_glob = istart+1, iend ! +1 because PETSc indexes from 0
      row_loc = row_glob - istart
      call MatGetRow( A, row_glob-1, ncols, cols_, vals_, ierr)
      nnz_row_loc( row_loc) = ncols
      call MatRestoreRow( A, row_glob-1, ncols, cols_, vals_, ierr)
    end do

    nnz_loc = sum( nnz_row_loc)

    call AA%allocate( m_glob, n_glob, m_loc, n_loc, nnz_loc)

    ! Copy data from the PETSc matrix to the CSR arrays
    do row_glob = istart+1, iend ! +1 because PETSc indexes from 0
      call MatGetRow( A, row_glob-1, ncols, cols_, vals_, ierr)
      do k = 1, ncols
        call AA%add_entry( row_glob, cols_( k)+1, vals_( k))
      end do
      call MatRestoreRow( A, row_glob-1, ncols, cols_, vals_, ierr)
    end do

    call AA%finalise()

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine mat_petsc2CSR

  subroutine mat_CSR2petsc( AA, A)
    ! Convert a CSR-format matrix in regular Fortran arrays to a PETSC parallel matrix
    !
    ! NOTE: the PETSc documentation seems to advise against using the MatCreateMPIAIJWithArrays
    !       routine used here. However, for the advised way of using MatSetValues with preallocation
    !       I've not been able to find a way that is fast enough to be useful without having to
    !       preallocate -WAY- too much memory. Especially for the remapping matrices, which
    !       can have hundreds or even thousands of non-zero elements per row, this can make the
    !       model run hella slow, whereas the current solution seems to work perfectly. So there you go.

    ! In/output variables:
    type(type_CSR_matrix_dp), intent(in   ) :: AA
    type(tMat),               intent(  out) :: A

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'mat_CSR2petsc'
    integer                             :: ierr
    integer                             :: i, k1, k2, nnz_proc, ii, k, kk
    integer,  dimension(:), allocatable :: ptr_proc, ind_proc
    real(dp), dimension(:), allocatable :: val_proc

    ! Add routine to path
    call init_routine( routine_name)

    if (.not. AA%is_finalised) call crash('A is not finalised')

    ! Determine number of non-zeros for this process
    nnz_proc = AA%nnz

    ! allocate memory for local CSR-submatrix
    allocate( ptr_proc( 0:AA%m_loc      ))
    allocate( ind_proc( 0:nnz_proc   - 1))
    allocate( val_proc( 0:nnz_proc   - 1))

    ! Copy matrix data
    do i = AA%i1, AA%i2

      ! ptr
      ii = i - AA%i1
      ptr_proc( ii) = AA%ptr( i) - AA%ptr( AA%i1)

      ! index and val
      k1 = AA%ptr( i)
      k2 = AA%ptr( i+1) - 1
      do k = k1, k2
        kk = k - AA%ptr( AA%i1)
        ind_proc( kk) = AA%ind( k) - 1
        val_proc( kk) = AA%val( k)
      end do

    end do
    ! Last row
    ptr_proc( AA%m_loc) = AA%ptr( AA%i2+1) - AA%ptr( AA%i1)

    ! Create PETSc matrix
    call MatCreateMPIAIJWithArrays( PETSC_COMM_WORLD, AA%m_loc, AA%n_loc, AA%m, AA%n, ptr_proc, ind_proc, val_proc, A, ierr)

    ! Assemble matrix and vectors, using the 2-step process:
    !   MatAssemblyBegin(), MatAssemblyEnd()
    ! Computations can be done while messages are in transition
    ! by placing code between these two statements.

    call MatAssemblyBegin( A, MAT_FINAL_ASSEMBLY, ierr)
    call MatAssemblyEnd(   A, MAT_FINAL_ASSEMBLY, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine mat_CSR2petsc

end submodule petsc_matrices
