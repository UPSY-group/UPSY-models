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

  ! These parameters will be set automatically when compiling the code!
  character(len=*), parameter :: git_commit_hash         = 'INVALID'
  logical,          parameter :: has_uncommitted_changes = .false.
  character(len=*), parameter :: petsc_version           = 'INVALID'
  character(len=*), parameter :: netcdf_version          = 'INVALID'
  character(len=*), parameter :: openmpi_version         = 'INVALID'

contains

  subroutine print_git_commit_hash_and_package_versions

    ! Safety
    if (git_commit_hash == 'INVALID') call crash('Invalid git commit hash - check the compile script!')

    if (par%primary) then

      write(0,'(A)') ''
      write(0,'(A)') ' Running UFEMISM from git commit ' // colour_string( trim( git_commit_hash), 'pink')
      write(0,'(A)') ' PETSc version:   ' // trim( petsc_version)
      write(0,'(A)') ' NetCDF version:  ' // trim( netcdf_version)
      write(0,'(A)') ' OpenMPI version: ' // trim( openmpi_version)
      if (has_uncommitted_changes) then
        write(0,'(A)') colour_string( ' WARNING: You have uncommitted changes; the current simulation might not be reproducible!', &
          'yellow')
      end if

    end if

  end subroutine print_git_commit_hash_and_package_versions

end module git_commit_hash_and_package_versions