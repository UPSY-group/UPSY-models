module ct_PETSc_matrix_solving

  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: warning, crash
  use basic_model_utilities, only: list_files_in_folder
  use mpi_basic, only: par
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use parallel_array_info_type, only: type_par_arr_info
  use netcdf_io_main, only: open_existing_netcdf_file_for_reading, inquire_var, read_var_primary, &
    close_netcdf_file
  use mpi_distributed_memory, only: distribute_from_primary
  use petsc_basic, only: solve_matrix_equation_CSR_PETSc
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPI_WTIME

  implicit none

  private

  public :: run_all_PETSc_matrix_solving_tests

contains

  subroutine run_all_PETSc_matrix_solving_tests( foldername_output, foldername_test_matrix_equations)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_output, foldername_test_matrix_equations

    ! Local variables:
    character(len=*), parameter                    :: routine_name = 'run_all_PETSc_matrix_solving_tests'
    character(len=1024), dimension(7)              :: PETSc_KSPtypes
    character(len=1024), dimension(4)              :: PETSc_PCtypes
    character(len=1024), dimension(:), allocatable :: list_of_test_matrix_names
    character(len=:), allocatable                  :: filename_output
    integer                                        :: i
    character(len=:), allocatable                  :: test_matrix_equation_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Also tried:
    !   - KSPCG, KSPPIPECG, KSPCR, KSPPIPECR, KSPMINRES (don't work, require symmetric of even symmetric positive-definite matrices)
    !   - KSPGCR (doesn't converge)
    PETSc_KSPtypes( 1) = 'bicg'
    PETSc_KSPtypes( 2) = 'bicgstab'
    PETSc_KSPtypes( 3) = 'ibicgstab'
    PETSc_KSPtypes( 4) = 'gmres'
    PETSc_KSPtypes( 5) = 'pipegmres'
    PETSc_KSPtypes( 6) = 'fgmres'
    PETSc_KSPtypes( 7) = 'lgmres'

    ! Also tried:
    !   - PCILU (returns errors)
    !   - PCHYPRE+PCHYPRESetType( precond, 'boomeramg', perr) (very inaccurate+slow)
    !   - PCJACOBI (very inaccurate)
    !   - PCNONE   (very slow)
    PETSc_PCtypes( 1) = 'bjacobi'
    PETSc_PCtypes( 3) = 'gamg'
    PETSc_PCtypes( 2) = 'asm'
    PETSc_PCtypes( 4) = 'gasm'

    call list_test_matrix_equations( foldername_test_matrix_equations, list_of_test_matrix_names)

    call create_output_text_file( foldername_output, PETSc_KSPtypes, PETSc_PCtypes, list_of_test_matrix_names, filename_output)

    do i = 1, size( list_of_test_matrix_names,1)
      test_matrix_equation_name = trim( list_of_test_matrix_names( i))
      call run_all_PETSc_matrix_solving_tests_on_matrix_equation( filename_output, &
        PETSc_KSPtypes, PETSc_PCtypes, foldername_test_matrix_equations, test_matrix_equation_name, 'init')
      call run_all_PETSc_matrix_solving_tests_on_matrix_equation( filename_output, &
        PETSc_KSPtypes, PETSc_PCtypes, foldername_test_matrix_equations, test_matrix_equation_name, 'zero')
    end do

    call finalise_output_text_file( filename_output)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_PETSc_matrix_solving_tests

  subroutine list_test_matrix_equations( foldername_test_matrix_equations, list_of_test_matrix_names)

    ! In/output variables
    character(len=*),                            intent(in   ) :: foldername_test_matrix_equations
    character(len=*), dimension(:), allocatable, intent(inout) :: list_of_test_matrix_names

    ! Local variables:
    character(len=*), parameter                    :: routine_name = 'list_test_matrix_equations'
    character(len=1024), dimension(:), allocatable :: list_of_filenames
    character(len=:), allocatable                  :: filename
    integer                                        :: n_valid_matrix_equations, i, j

    ! Add routine to call stack
    call init_routine( routine_name)

    ! List all files in the folder
    call list_files_in_folder( foldername_test_matrix_equations, list_of_filenames)

    n_valid_matrix_equations = 0
    do i = 1, size( list_of_filenames,1)
      filename = trim( list_of_filenames( i))
      if (UPSY%stru%endswith( filename, '_b.nc')) then
        n_valid_matrix_equations = n_valid_matrix_equations + 1
      end if
    end do

    allocate( list_of_test_matrix_names( n_valid_matrix_equations))

    j = 0
    do i = 1, size( list_of_filenames,1)
      filename = trim( list_of_filenames( i))
      if (UPSY%stru%endswith( filename, '_b.nc')) then
        j = j+1
        list_of_test_matrix_names( j) = filename( 1: len_trim( filename)-5)
      end if
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine list_test_matrix_equations

  subroutine create_output_text_file( foldername_output, PETSc_KSPtypes, PETSc_PCtypes, list_of_test_matrix_names, filename_output)

    ! In/output variables
    character(len=*),               intent(in   ) :: foldername_output
    character(len=*), dimension(:), intent(in   ) :: PETSc_KSPtypes
    character(len=*), dimension(:), intent(in   ) :: PETSc_PCtypes
    character(len=*), dimension(:), intent(in   ) :: list_of_test_matrix_names
    character(len=:), allocatable,  intent(inout) :: filename_output

    ! Local variables:
    character(len=*), parameter :: routine_name = 'create_output_text_file'
    integer                     :: u, kspi, pci, mi

    ! Add routine to call stack
    call init_routine( routine_name)

    filename_output = trim( foldername_output) // '/PETSc_solver_results.m'

    if (par%primary) then

      open( newunit = u, status = 'new', action = 'write', file = filename_output)

      write( u, '(a)') 'function [KSPs, PCs, matrix_names, results] = PETSc_solver_results'
      write( u, '(a)') '% Results of solving test matrix equations with different PETSc KSP and PC options'

      write( u, '(a)') ''
      write( u, '(a)') 'KSPs = {'
      do kspi = 1, size( PETSc_KSPtypes,1)
      write( u, '(a)') "  '" // trim( PETSc_KSPtypes( kspi)) // "'"
      end do
      write( u, '(a)') '  };'

      write( u, '(a)') ''
      write( u, '(a)') 'PCs = {'
      do pci = 1, size( PETSc_PCtypes,1)
      write( u, '(a)') "  '" // trim( PETSc_PCtypes( pci)) // "'"
      end do
      write( u, '(a)') '  };'

      write( u, '(a)') ''
      write( u, '(a)') 'matrix_names = {'
      do mi = 1, size( list_of_test_matrix_names,1)
      write( u, '(a)') "  '" // trim( list_of_test_matrix_names( mi)) // "_init'"
      write( u, '(a)') "  '" // trim( list_of_test_matrix_names( mi)) // "_zero'"
      end do
      write( u, '(a)') '  };'

      write( u, '(a)') ''
      write( u, '(a)') 'results = [];'
      write( u, '(a)') ''

      close( u)

    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_output_text_file

  subroutine finalise_output_text_file( filename_output)

    ! In/output variables
    character(len=*), intent(in) :: filename_output

    ! Local variables:
    character(len=*), parameter :: routine_name = 'finalise_output_text_file'
    integer                     :: u

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) then

      open( file = filename_output, newunit = u, status = 'old', action = 'write', access = 'append')

      write( u,'(a)') ''
      write( u,'(a)') 'end'

      close( u)

    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine finalise_output_text_file

  subroutine run_all_PETSc_matrix_solving_tests_on_matrix_equation( filename_output, &
    PETSc_KSPtypes, PETSc_PCtypes, foldername_test_matrix_equations, test_matrix_equation_name, choice_x_init)

    ! In/output variables:
    character(len=*),               intent(in) :: filename_output
    character(len=*), dimension(:), intent(in) :: PETSc_KSPtypes
    character(len=*), dimension(:), intent(in) :: PETSc_PCtypes
    character(len=*),               intent(in) :: foldername_test_matrix_equations
    character(len=*),               intent(in) :: test_matrix_equation_name
    character(len=*),               intent(in) :: choice_x_init

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'run_all_PETSc_matrix_solving_tests_on_matrix_equation'
    type(type_CSR_matrix_dp)            :: A_CSR
    real(dp), dimension(:), allocatable :: b, x_init, x, x_PETSc
    integer                             :: iksp, ipc
    character(len=:), allocatable       :: PETSC_KSPtype, PETSc_PCtype
    character(len=:), allocatable       :: test_matrix_equation_name_init

    ! Add routine to call stack
    call init_routine( routine_name)

    call read_test_matrix_equation( foldername_test_matrix_equations, test_matrix_equation_name, &
      A_CSR, b, x_init, x)

    select case (choice_x_init)
    case default
      call crash('invalid choice_x_init "' // trim( choice_x_init) // '"')
    case ('init')
      ! Use the initial guess read from the NetCDF file
    case ('zero')
      x_init = 0._dp
    end select
    test_matrix_equation_name_init = trim( test_matrix_equation_name) // '_' // trim( choice_x_init)

    do iksp = 1, size( PETSc_KSPtypes,1)
      do ipc = 1, size( PETSc_PCtypes,1)

        PETSC_KSPtype = trim( PETSc_KSPtypes( iksp))
        PETSc_PCtype  = trim( PETSc_PCtypes ( ipc))

        call run_PETSc_matrix_solving_test_on_matrix_equation_KSP_PC( filename_output, &
          test_matrix_equation_name_init, A_CSR, b, x_init, x, PETSc_KSPtype, PETSc_PCtype)

      end do
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_PETSc_matrix_solving_tests_on_matrix_equation

  subroutine read_test_matrix_equation( foldername, test_matrix_equation_name, A_CSR, b, x_init, x)

    ! In/output variables:
    character(len=*),                    intent(in   ) :: foldername, test_matrix_equation_name
    type(type_CSR_matrix_dp),            intent(inout) :: A_CSR
    real(dp), dimension(:), allocatable, intent(inout) :: b, x_init, x

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_all_PETSc_matrix_solving_tests'

    ! Add routine to call stack
    call init_routine( routine_name)

    call A_CSR%read_from_dist_NetCDFs( foldername, test_matrix_equation_name)
    call read_test_matrix_equation_vector( foldername, test_matrix_equation_name, A_CSR%pai_x, 'b'     , b     )
    call read_test_matrix_equation_vector( foldername, test_matrix_equation_name, A_CSR%pai_x, 'x_init', x_init)
    call read_test_matrix_equation_vector( foldername, test_matrix_equation_name, A_CSR%pai_x, 'x'     , x     )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_test_matrix_equation

  subroutine read_test_matrix_equation_vector( foldername, test_matrix_equation_name, pai, vec_name, d)

    ! In/output variables:
    character(len=*),                    intent(in   ) :: foldername, test_matrix_equation_name
    type(type_par_arr_info),             intent(in   ) :: pai
    character(len=*),                    intent(in   ) :: vec_name
    real(dp), dimension(:), allocatable, intent(inout) :: d

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'run_all_PETSc_matrix_solving_tests'
    character(len=:), allocatable       :: filename, var_name
    integer                             :: ncid, id_var
    real(dp), dimension(:), allocatable :: d_tot

    ! Add routine to call stack
    call init_routine( routine_name)

    filename = trim( foldername) // '/' // trim( test_matrix_equation_name) // '_' // trim( vec_name) // '.nc'

    if (par%primary) then
      allocate( d_tot( 1: pai%n))
    else
      allocate( d_tot( 0))
    end if

    call open_existing_netcdf_file_for_reading( filename, ncid)
    var_name = trim( test_matrix_equation_name) // '_' // trim( vec_name)
    call inquire_var( filename, ncid, var_name, id_var)
    call read_var_primary( filename, ncid, id_var, d_tot)
    call close_netcdf_file( ncid)

    allocate( d( pai%i1 : pai%i2))
    call distribute_from_primary( d, d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_test_matrix_equation_vector

  subroutine run_PETSc_matrix_solving_test_on_matrix_equation_KSP_PC( filename_output, &
    test_matrix_equation_name, A_CSR, b, x_init, x, PETSc_KSPtype, PETSc_PCtype)

    ! In/output variables:
    character(len=*),         intent(in) :: filename_output
    character(len=*),         intent(in) :: test_matrix_equation_name
    type(type_CSR_matrix_dp), intent(in) :: A_CSR
    real(dp), dimension(:),   intent(in) :: b, x_init, x
    character(len=*),         intent(in) :: PETSc_KSPtype, PETSc_PCtype

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'run_PETSc_matrix_solving_test_on_matrix_equation_KSP_PC'
    real(dp), dimension(:), allocatable :: x_PETSc
    real(dp)                            :: rtol, abstol
    real(dp)                            :: tstart, tstop, tcomp
    integer                             :: n_Axb_its
    real(dp)                            :: sum_sq, n_sq, rmse
    integer                             :: ierr
    character(len=60)                   :: matrix_name_str
    character(len=20)                   :: KSP_str
    character(len=20)                   :: PC_str
    character(len=:), allocatable       :: str

    ! Add routine to call stack
    call init_routine( routine_name)

    allocate( x_PETSc( A_CSR%pai_x%i1 : A_CSR%pai_x%i2), source = x_init)

    rtol   = 1E-6  ! Values taken from ISMIP-HOM
    abstol = 1E-4
    tstart = MPI_WTIME()
    call solve_matrix_equation_CSR_PETSc( A_CSR, b, x_PETSc, rtol, abstol, n_Axb_its, &
      PETSc_KSPtype, PETSc_PCtype)
    tstop = MPI_WTIME()
    tcomp = tstop - tstart

    ! Calculate RMSE with respect to exact solution
    sum_sq = sum( (x_PETSc - x)**2)
    n_sq   = real( size( x_PETSc,1),dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, sum_sq, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, n_sq  , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    rmse = sqrt( sum_sq / n_sq)

    ! If n_Axb_its = 0, the solver didn't work
    if (n_Axb_its == 0) n_Axb_its = 1000 ! maximum allowed value used inside solve_matrix_equation_CSR_PETSc
    if (isnan( rmse)) rmse = 1e36_dp

    matrix_name_str = UPSY%stru%colour_string( trim( test_matrix_equation_name), 'light blue')
    KSP_str         = UPSY%stru%colour_string( trim( PETSc_KSPtype            ), 'light blue')
    PC_str          = UPSY%stru%colour_string( trim( PETSc_PCtype             ), 'light blue')
    str = '  Matrix ' // matrix_name_str // ' - KSP: ' // KSP_str // ', PC: ' // PC_str
    if (par%primary) write(0,'(a,i8,a,e20.14,a,e20.14)') str // &
      ' - solved in ', n_Axb_its, ' its, rmse = ', rmse, ', tcomp = ', tcomp

    ! Write test results to output
    call write_test_results_to_text_output( filename_output, test_matrix_equation_name, &
      PETSc_KSPtype, PETSc_PCtype, n_Axb_its, rmse, tcomp)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_PETSc_matrix_solving_test_on_matrix_equation_KSP_PC

  subroutine write_test_results_to_text_output( filename_output, test_matrix_equation_name, &
      PETSc_KSPtype, PETSc_PCtype, n_Axb_its, rmse, tcomp)

    ! In/output variables:
    character(len=*), intent(in) :: filename_output, test_matrix_equation_name
    character(len=*), intent(in) :: PETSc_KSPtype, PETSc_PCtype
    integer,          intent(in) :: n_Axb_its
    real(dp),         intent(in) :: rmse, tcomp

    ! Local variables:
    character(len=*), parameter :: routine_name = 'write_test_results_to_text_output'
    integer                     :: u

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) then

      open( file = filename_output, newunit = u, status = 'old', action = 'write', access = 'append')

      write( u,'(a)') ''
      write( u,'(a,a     ,a)') "result.matrix_name = '", trim( test_matrix_equation_name), "';"
      write( u,'(a,a     ,a)') "result.KSP         = '", trim( PETSc_KSPtype), "';"
      write( u,'(a,a     ,a)') "result.PC          = '", trim( PETSc_PCtype), "';"
      write( u,'(a,i8    ,a)') "result.n_Axb_its   = ", n_Axb_its, ";"
      write( u,'(a,e20.14,a)') "result.rmse        = ", rmse, ";"
      write( u,'(a,e20.14,a)') "result.tcomp       = ", tcomp, ";"

      write( u,'(a)') 'if isempty( results)'
      write( u,'(a)') '  results = result;'
      write( u,'(a)') 'else'
      write( u,'(a)') '  results( end+1) = result;'
      write( u,'(a)') 'end'

      close( u)

    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_test_results_to_text_output

end module ct_PETSc_matrix_solving
