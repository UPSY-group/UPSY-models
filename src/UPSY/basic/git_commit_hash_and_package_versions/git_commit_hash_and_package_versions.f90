module git_commit_hash_and_package_versions

  use mpi_basic, only: par
  use string_module, only: colour_string
  use crash_mod, only: crash

  implicit none

  private

  public :: print_git_commit_hash_and_package_versions
  public :: git_commit_hash
  public :: has_uncommitted_changes
  public :: petsc_version
  public :: netcdf_version
  public :: openmpi_version
  public :: compiler_version
  public :: compiler_flags

  ! These parameters will be set automatically when compiling the code!
  character(len=*), parameter :: git_commit_hash         = '1a57eeea77fa306d5ac060c108cdb5ac036de765'
  logical,          parameter :: has_uncommitted_changes = .true.
  character(len=*), parameter :: petsc_version           = '3.25.1'
  character(len=*), parameter :: netcdf_version          = '4.10.0'
  character(len=*), parameter :: openmpi_version         = '5.0.10'
  character(len=*), parameter :: compiler_version        = 'GNU Fortran (conda-forge gcc 14.3.0-19) 14.3.0'
  character(len=*), parameter :: compiler_flags          = '-fdiagnostics-color=always -O3 -Wall -ffree-line-length-none -cpp -fimplicit-none -g -march=native'

contains

  subroutine print_git_commit_hash_and_package_versions

    ! Safety
    if (git_commit_hash  == 'INVALID') call crash('Invalid git commit hash - check the compile script!')
    if (petsc_version    == 'INVALID') call crash('Invalid PETSc version number - check the compile script!')
    if (netcdf_version   == 'INVALID') call crash('Invalid NetCDF version number - check the compile script!')
    if (openmpi_version  == 'INVALID') call crash('Invalid OpenMPI version number - check the compile script!')
    if (compiler_version == 'INVALID') call crash('Invalid compiler version - check the compile script!')
    if (compiler_flags   == 'INVALID') call crash('Invalid compiler flags - check the compile script!')

    if (par%primary) then

      write(0,'(A)') ''
      write(0,'(A)') ' This program was compiled from git commit ' // colour_string( trim( git_commit_hash), 'pink')

      if (has_uncommitted_changes) then
        write(0,'(A)') colour_string( ' WARNING: You have uncommitted changes; the current simulation might not be reproducible!', &
          'yellow')
      end if

      write(0,'(A)') ' The following package versions were used:'
      write(0,'(A)') colour_string( '  PETSc  ', 'pink') // ' - version: ' // colour_string( trim( petsc_version), 'pink')
      write(0,'(A)') colour_string( '  NetCDF ', 'pink') // ' - version: ' // colour_string( trim( netcdf_version), 'pink')
      write(0,'(A)') colour_string( '  OPenMPI', 'pink') // ' - version: ' // colour_string( trim( openmpi_version), 'pink')
      write(0,'(A)') '  Fortran compiler:' // ' ' // colour_string( trim( compiler_version), 'pink')
      write(0,'(A)') '  Compiler flags  :' // ' ' // colour_string( trim( compiler_flags), 'pink')

    end if

  end subroutine print_git_commit_hash_and_package_versions

end module git_commit_hash_and_package_versions