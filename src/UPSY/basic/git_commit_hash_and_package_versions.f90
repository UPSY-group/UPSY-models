module git_commit_hash_and_package_versions

  use mpi_basic, only: par
  use string_module, only: colour_string

  implicit none

  private

  public :: print_git_commit_hash_and_package_versions
  public :: git_commit_hash
  public :: has_uncommitted_changes

  character(len=*), parameter :: git_commit_hash         = '0d64f4389e2eb643fb84af3e08abe0142cd9d0d8'
  logical,          parameter :: has_uncommitted_changes = .false.

contains

  subroutine print_git_commit_hash_and_package_versions

    if (par%primary) then

      write(0,'(A)') ''
      write(0,'(A)') ' Running UFEMISM from git commit ' // colour_string( trim( git_commit_hash), 'pink')
      if (has_uncommitted_changes) then
        write(0,'(A)') colour_string( ' WARNING: You have uncommitted changes; the current simulation might not be reproducible!', &
          'yellow')
      end if

    end if

  end subroutine print_git_commit_hash_and_package_versions

end module git_commit_hash_and_package_versions