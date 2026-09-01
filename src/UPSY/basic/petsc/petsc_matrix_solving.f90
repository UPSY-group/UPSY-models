module petsc_matrix_solving

  use precisions, only: dp
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use petsc, only: tMat, tVec, tKSP, tPC, MatGetSize, MatGetLocalSize, &
    KSPcreate, PETSC_COMM_WORLD, KSPSetOperators, KSPSetInitialGuessNonzero, &
    PETSC_TRUE, KSPSetTolerances, PETSC_DEFAULT_REAL, KSPSetFromOptions, &
    KSPGetIterationNumber, KSPDestroy, VecDestroy, KSPGetType, PCGetType, &
    KSPSolve, KSPGMRES, KSPLGMRES, KSPFGMRES, KSPPGMRES, KSPBICG, KSPBCGS, KSPIBCGS, &
    PCBJACOBI, PCASM, PCGASM, PCGAMG, PCNONE, KSPSetType, PCSetType, KSPGetPC, MatDestroy
  use petsc_matrices, only: mat_CSR2petsc
  use petsc_vectors, only: vec_double2petsc, vec_petsc2double
  use crash_mod, only: crash
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, colour_string

  implicit none

  private

  public :: solve_matrix_equation_CSR_PETSc

  character(len=:), allocatable :: PETSc_KSPtype_printed, PETSc_PCtype_printed

contains

  subroutine solve_matrix_equation_CSR_PETSc( AA, bb, xx, rtol, abstol, n_Axb_its, &
    PETSc_KSPtype, PETSc_PCtype)

    ! In/output variables:
    type(type_CSR_matrix_dp),            intent(in   ) :: AA
    real(dp), dimension(:),              intent(in   ) :: bb
    real(dp), dimension(:),              intent(inout) :: xx
    real(dp),                            intent(in   ) :: rtol, abstol
    integer,                             intent(out)   :: n_Axb_its
    character(len=*), optional,          intent(in   ) :: PETSc_KSPtype
    character(len=*), optional,          intent(in   ) :: PETSc_PCtype

    ! Local variables
    character(len=*), parameter :: routine_name = 'solve_matrix_equation_CSR_PETSc'
    type(tMat)                  :: A
    integer                     :: perr

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. AA%is_finalised) call crash('A is not finalised')

    call mat_CSR2petsc( AA, A)
    call solve_matrix_equation_PETSc( A, bb, xx, rtol, abstol, &
      n_Axb_its, PETSc_KSPtype, PETSc_PCtype)
    call MatDestroy( A, perr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_matrix_equation_CSR_PETSc

  subroutine solve_matrix_equation_PETSc( A, bb, xx, rtol, abstol, &
    n_Axb_its, PETSc_KSPtype, PETSc_PCtype)
    !< Solve the matrix equation using a Krylov solver from PETSc

    ! In/output variables:
    type(tMat),                 intent(in   ) :: A
    real(dp), dimension(:),     intent(in   ) :: bb
    real(dp), dimension(:),     intent(inout) :: xx
    real(dp),                   intent(in   ) :: rtol, abstol
    integer,                    intent(  out) :: n_Axb_its
    character(len=*), optional, intent(in   ) :: PETSc_KSPtype
    character(len=*), optional, intent(in   ) :: PETSc_PCtype

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'solve_matrix_equation_PETSc'
    integer                       :: perr
    integer                       :: m, n, m_local, n_local
    type(tVec)                    :: b
    type(tVec)                    :: x
    type(tKSP)                    :: KSP_solver
    type(tPC)                     :: precond
    character(len=:), allocatable :: PETSc_KSPtype_, PETSc_PCtype_

    ! Add routine to path
    call init_routine( routine_name)

    ! Optional arguments
    if (present( PETSc_KSPtype)) then
      PETSc_KSPtype_ = PETSc_KSPtype
    else
      PETSc_KSPtype_ = 'gmres'
    end if

    if (present( PETSc_PCtype)) then
      PETSc_PCtype_ = PETSc_PCtype
    else
      PETSc_PCtype_ = 'bjacobi'
    end if

    ! Safety
    call MatGetSize( A, m, n, perr)
    call MatGetLocalSize( A, m_local, n_local, perr)

    if (n_local /= size( xx,1) .or. m_local /= size( bb,1)) then
      call crash('matrix and vector sub-sizes dont match')
    end if

    ! == Set up right-hand side and solution vectors as PETSc data structures
    ! =======================================================================

    call vec_double2petsc( xx, x)
    call vec_double2petsc( bb, b)

    ! Set up the solver
    ! =================

    ! Set up the KSP solver
    call KSPcreate( PETSC_COMM_WORLD, KSP_solver, perr)

    ! Set operators. Here the matrix that defines the linear system
    ! also serves as the preconditioning matrix.
    call KSPSetOperators( KSP_solver, A, A, perr)

    ! Set the type of KSP solver
    ! NOTE: copied the list of options from ISSM. From some brief
    !       manual testing, it seems that 'bicg' WILDLY outperforms
    !       all other options, so let's leave it at that for now.
    select case (PETSc_KSPtype_)
    case default
      call crash('unknown PETSc_KSPtype "' // trim( PETSc_KSPtype_) // '"')
    case ('gmres')
      call KSPSetType( KSP_solver, KSPGMRES, perr)
    case ('lgmres')
      call KSPSetType( KSP_solver, KSPLGMRES, perr)
    case ('fgmres')
      call KSPSetType( KSP_solver, KSPFGMRES, perr)
    case ('pipegmres')
      call KSPSetType( KSP_solver, KSPPGMRES, perr)
    case ('bicg')
      call KSPSetType( KSP_solver, KSPBICG, perr)
    case ('bicgstab')
      call KSPSetType( KSP_solver, KSPBCGS, perr)
    case ('ibicgstab')
      call KSPSetType( KSP_solver, KSPIBCGS, perr)
    end select

    ! Set preconditioner
    call KSPGetPC( KSP_solver, precond, perr)
    select case (PETSc_PCtype_)
    case default
      call crash('unknown PETSc_PCtype "' // trim( PETSc_PCtype_) // '"')
    case ('bjacobi')
      call PCSetType( precond, PCBJACOBI, perr)
    case ('asm')
      call PCSetType( precond, PCASM, perr)
    case ('gasm')
      call PCSetType( precond, PCGASM, perr)
    case ('gamg')
      call PCSetType( precond, PCGAMG, perr)
    case ('none')
      call PCSetType( precond, PCNONE, perr)
    end select

    ! Make sure PETSc knows we're starting from an initial guess
    call KSPSetInitialGuessNonzero( KSP_solver, PETSC_TRUE, perr)

    ! Iterative solver tolerances
    call KSPSetTolerances( KSP_solver, rtol, abstol, PETSC_DEFAULT_REAL, 1000, perr)

    ! Set runtime options, e.g.,
    !     -ksp_type <type> -pc_type <type> -ksp_monitor -ksp_rtol <rtol>
    ! These options will override those specified above as long as
    ! KSPSetFromOptions() is called _after_ any other customization routines.
    call KSPSetFromOptions( KSP_solver, perr)
    call print_PETSc_KSP_PC_type( KSP_solver, precond)

    ! == Solve Ax=b
    ! =============

    ! Solve the linear system
    call solve_matrix_equation_PETSc_KSPSolve( KSP_solver, b, x)

    ! Find out how many iterations it took
    call KSPGetIterationNumber( KSP_solver, n_Axb_its, perr)

    ! Get the solution back to the native UFEMISM storage structure
    call vec_petsc2double( x, xx)

    ! Clean up after yourself
    call KSPDestroy( KSP_solver, perr)
    call VecDestroy( x, perr)
    call VecDestroy( b, perr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_matrix_equation_PETSc

  subroutine print_PETSc_KSP_PC_type( KSP_solver, precond)
    ! Query which KSP and PC we're actually using, and tell the user

    ! In/output variables:
    type(tKSP), intent(in) :: KSP_solver
    type(tPC),  intent(in) :: precond

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'solve_matrix_equation_PETSc'
    character(len=:), allocatable :: PETSc_KSPtype_applied, PETSc_PCtype_applied
    logical                       :: do_print_KSP, do_print_PC
    integer                       :: perr

    ! Add routine to path
    call init_routine( routine_name)

    ! KSP solver
    call KSPGetType( KSP_solver, PETSc_KSPtype_applied, perr)

    do_print_KSP = .false.
    if (.not. allocated( PETSc_KSPtype_printed)) then
      do_print_KSP = .true.
      PETSc_KSPtype_printed = trim( PETSc_KSPtype_applied)
    else
      ! Maybe the KSP type has changed; check for this
      if (trim( PETSc_KSPtype_applied) /= PETSc_KSPtype_printed) then
        deallocate( PETSc_KSPtype_printed)
        do_print_KSP = .true.
        PETSc_KSPtype_printed = trim( PETSc_KSPtype_applied)
      end if
    end if

    if (do_print_KSP .and. par%primary) then
      write(0,*) '  Using PETSc KSP solver ', colour_string( PETSc_KSPtype_printed, 'light blue')
    end if

    ! Preconditioner
    call PCGetType( precond, PETSc_PCtype_applied, perr)

    do_print_PC = .false.
    if (.not. allocated( PETSc_PCtype_printed)) then
      do_print_PC = .true.
      PETSc_PCtype_printed = trim( PETSc_PCtype_applied)
    else
      ! Maybe the PC type has changed; check for this
      if (trim( PETSc_PCtype_applied) /= PETSc_PCtype_printed) then
        deallocate( PETSc_PCtype_printed)
        do_print_PC = .true.
        PETSc_PCtype_printed = trim( PETSc_PCtype_applied)
      end if
    end if

    if (do_print_PC .and. par%primary) then
      write(0,*) '  Using PETSc preconditioner ', colour_string( PETSc_PCtype_printed, 'light blue')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine print_PETSc_KSP_PC_type

  subroutine solve_matrix_equation_PETSc_KSPSolve( KSP_solver, b, x)
    ! Solve the matrix equation using a Krylov solver from PETSc
    !
    ! Wrapper for KSPSolve, so we can determine how much computation is spent
    ! on that relative to the initialisation and format conversion stuff.

    ! In/output variables:
    type(tKSP), intent(inout) :: KSP_solver
    type(tVec), intent(inout) :: b
    type(tVec), intent(inout) :: x

    ! Local variables:
    character(len=*), parameter :: routine_name = 'solve_matrix_equation_PETSc_KSPSolve'
    integer                     :: perr

    ! Add routine to path
    call init_routine( routine_name)

    ! Solve the linear system
    call KSPSolve( KSP_solver, b, x, perr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_matrix_equation_PETSc_KSPSolve

end module petsc_matrix_solving
