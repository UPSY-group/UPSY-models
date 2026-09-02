module ct_PETSc_SNES_Poisson

#include <petsc/finclude/petscsys.h>
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: warning
  use mpi_basic, only: par
  use mesh_types, only: type_mesh
  use netcdf_io_main, only: open_existing_netcdf_file_for_reading, setup_mesh_from_file, &
    close_netcdf_file
  use petsc, only: PetscErrorF, tDM, DMDestroy
  use petsc_dmplex, only: mesh_to_dmplex
  use string_module, only: strrep

  implicit none

  private

  public :: ct_solve_Poisson_eq_with_PETSc_SNES

contains

  subroutine ct_solve_Poisson_eq_with_PETSc_SNES( foldername_output, test_mesh_filenames, test_grid_filenames)

    ! In/output variables:
    character(len=*),               intent(in) :: foldername_output
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames
    character(len=*), dimension(:), intent(in) :: test_grid_filenames

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'ct_solve_Poisson_eq_with_PETSc_SNES'
    integer                       :: i_mesh
    type(type_mesh), allocatable  :: mesh
    character(len=:), allocatable :: filename_mesh
    integer                       :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '   Running PETSc SNES test: solve a simple Poisson equation'
    if (par%primary) write(0,*) ''

    ! Meshes read from files
    do i_mesh = 1, size( test_mesh_filenames)
      filename_mesh = test_mesh_filenames( i_mesh)
      allocate( mesh)
      call open_existing_netcdf_file_for_reading( filename_mesh, ncid)
      call setup_mesh_from_file( filename_mesh, ncid, mesh)
      call close_netcdf_file( ncid)
      call ct_solve_Poisson_eq_with_PETSc_SNES_on_mesh( foldername_output, mesh)
      deallocate( mesh)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ct_solve_Poisson_eq_with_PETSc_SNES

  subroutine ct_solve_Poisson_eq_with_PETSc_SNES_on_mesh( foldername_output, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_output
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'ct_solve_Poisson_eq_with_PETSc_SNES_on_mesh'
    type(tDM)                     :: dm
    integer                       :: ierr
    character(len=:), allocatable :: mesh_name_cleaned
    character(len=:), allocatable :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '    Using PETSc SNES to solve the Poisson equation on mesh ', trim( mesh%name), '...'

    call mesh_to_dmplex( mesh, dm)


    ! Clean up after yourself
    PetscCall( DMDestroy( dm, ierr))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ct_solve_Poisson_eq_with_PETSc_SNES_on_mesh

end module ct_PETSc_SNES_Poisson