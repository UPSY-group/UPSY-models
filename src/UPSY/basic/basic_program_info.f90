module basic_program_info

  private

  public :: program_name, routine_path, n_MPI_windows_used

  character(len=:), allocatable :: program_name
  character(len=:), allocatable :: routine_path
  integer                       :: n_MPI_windows_used

end module basic_program_info