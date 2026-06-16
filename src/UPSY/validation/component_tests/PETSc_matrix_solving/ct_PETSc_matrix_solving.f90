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

  implicit none

  private

  public :: run_all_PETSc_matrix_solving_tests

contains

  subroutine run_all_PETSc_matrix_solving_tests( foldername_output, foldername_test_matrix_equations)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_output, foldername_test_matrix_equations

    ! Local variables:
    character(len=*), parameter                    :: routine_name = 'run_all_PETSc_matrix_solving_tests'
    character(len=1024), dimension(:), allocatable :: list_of_test_matrix_names
    integer                                        :: i
    character(len=:), allocatable                  :: test_matrix_equation_name

    ! Add routine to call stack
    call init_routine( routine_name)

    call list_test_matrix_equations( foldername_test_matrix_equations, list_of_test_matrix_names)

    do i = 1, size( list_of_test_matrix_names,1)
      test_matrix_equation_name = trim( list_of_test_matrix_names( i))
      call run_all_PETSc_matrix_solving_tests_on_matrix_equation( foldername_output, &
        foldername_test_matrix_equations, test_matrix_equation_name)
    end do

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

  subroutine run_all_PETSc_matrix_solving_tests_on_matrix_equation( foldername_output, &
    foldername_test_matrix_equations, test_matrix_equation_name)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_output, foldername_test_matrix_equations, test_matrix_equation_name

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'run_all_PETSc_matrix_solving_tests_on_matrix_equation'
    type(type_CSR_matrix_dp)            :: A_CSR
    real(dp), dimension(:), allocatable :: b, x_init, x, x_PETSc
    character(len=1024), dimension(10)  :: PETSc_KSPtypes
    character(len=1024), dimension(6)   :: PETSc_PCtypes
    integer                             :: iksp, ipc
    character(len=:), allocatable       :: PETSC_KSPtype, PETSc_PCtype

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '  Running PETSc matrix solving tests on matrix equation ', &
      UPSY%stru%colour_string( test_matrix_equation_name, 'light blue'), '...'

    call read_test_matrix_equation( foldername_test_matrix_equations, test_matrix_equation_name, &
      A_CSR, b, x_init, x)

    PETSc_KSPtypes(  1) = 'gmres'
    PETSc_KSPtypes(  2) = 'pipegmres'
    PETSc_KSPtypes(  3) = 'cg'
    PETSc_KSPtypes(  4) = 'pipecg'
    PETSc_KSPtypes(  5) = 'bicg'
    PETSc_KSPtypes(  6) = 'bicgstab'
    PETSc_KSPtypes(  7) = 'ibicgstab'
    PETSc_KSPtypes(  8) = 'minres'
    PETSc_KSPtypes(  9) = 'cr'
    PETSc_KSPtypes( 10) = 'pipecr'

    PETSc_PCtypes( 1) = 'bjacobi'
    PETSc_PCtypes( 2) = 'asm'
    PETSc_PCtypes( 3) = 'gamg'
    PETSc_PCtypes( 4) = 'gasm'
    PETSc_PCtypes( 5) = 'jacobi'
    PETSc_PCtypes( 6) = 'none'

    do iksp = 1, size( PETSc_KSPtypes,1)
      do ipc = 1, size( PETSc_PCtypes,1)

        PETSC_KSPtype = trim( PETSc_KSPtypes( iksp))
        PETSc_PCtype  = trim( PETSc_PCtypes ( ipc))

        call run_all_PETSc_matrix_solving_tests_on_matrix_equation_KSP_PC( foldername_output, &
          A_CSR, b, x_init, x, PETSc_KSPtype, PETSc_PCtype)

      end do
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_PETSc_matrix_solving_tests_on_matrix_equation

  subroutine run_all_PETSc_matrix_solving_tests_on_matrix_equation_KSP_PC( foldername_output, &
    A_CSR, b, x_init, x, PETSc_KSPtype, PETSc_PCtype)

    ! In/output variables:
    character(len=*),         intent(in) :: foldername_output
    type(type_CSR_matrix_dp), intent(in) :: A_CSR
    real(dp), dimension(:),   intent(in) :: b, x_init, x
    character(len=*),         intent(in) :: PETSc_KSPtype, PETSc_PCtype

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_all_PETSc_matrix_solving_tests_on_matrix_equation_KSP_PC'

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '   KSP: ', &
      UPSY%stru%colour_string( PETSc_KSPtype, 'light blue'), ', PC: ', &
      UPSY%stru%colour_string( PETSc_PCtype, 'light blue')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_PETSc_matrix_solving_tests_on_matrix_equation_KSP_PC

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

end module ct_PETSc_matrix_solving
