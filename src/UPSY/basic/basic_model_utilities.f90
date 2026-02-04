module basic_model_utilities

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use basic_program_info, only: program_name
  use string_module, only: colour_string, insert_val_into_string_dp, insert_val_into_string_int
  use crash_mod, only: crash

  implicit none

  private

  public :: get_git_commit_hash, git_commit_hash
  public :: check_for_uncommitted_changes, has_uncommitted_changes
  public :: generate_procedural_output_dir_name
  public :: print_model_start
  public :: print_model_end

  ! The hash of the current git commit
  character(len=1024) :: git_commit_hash
  logical             :: has_uncommitted_changes = .false.

contains

  subroutine get_git_commit_hash( git_commit_hash)

    ! In/output variables:
    character(len=*), intent(out) :: git_commit_hash

    ! Local variables:
    character(len=256), parameter :: routine_name = 'get_git_commit_hash'
    character(len=256), parameter :: filename_git_commit_hash = 'git_commit_hash.txt'
    integer                       :: ierr, ios
    integer, parameter            :: git_commit_hash_file_unit = 1847

    ! Add routine to path
    call init_routine( routine_name)

    ! Create a text file containing the hash of the current git commit
    call system( 'git rev-parse HEAD > ' // trim(filename_git_commit_hash), ierr)
    if (ierr /= 0) call crash('failed to obtain hash of current git commit')

    ! Read the hash from the temporary commit hash file
    open( unit = git_commit_hash_file_unit, file = filename_git_commit_hash, iostat = ios)
    if (ios /= 0) call crash('couldnt open temporary commit hash file "' // trim( filename_git_commit_hash) // '"!')
    read( unit = git_commit_hash_file_unit, fmt = '(A)', iostat = ios) git_commit_hash
    if (ios < 0) call crash('couldnt read commit hash from the temporary commit hash file')
    close( unit = git_commit_hash_file_unit)

    ! Delete the temporary commit hash file
    call system( 'rm -f ' // trim( filename_git_commit_hash), ierr)
    if (ierr /= 0) call crash('failed to delete temporary commit hash file')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine get_git_commit_hash

  subroutine check_for_uncommitted_changes

    ! Local variables:
    character(len=256), parameter :: routine_name = 'check_for_uncommitted_changes'
    character(len=256), parameter :: filename_git_status = 'git_status.txt'
    integer                       :: ierr, ios
    integer, parameter            :: git_status_file_unit = 1847
    character(len=1024)           :: single_line

    ! Add routine to path
    call init_routine( routine_name)

    ! Create a text file containing the output of git status
    call system( 'git status > ' // trim( filename_git_status), ierr)
    if (ierr /= 0) call crash('failed to write git status to text file')

    ! Check the temporary git status file for uncommitted changes
    open( unit = git_status_file_unit, file = filename_git_status, iostat = ios)
    if (ios /= 0) call crash('couldnt open temporary git status file "' // trim( filename_git_status) // '"!')

    do while (.true.)
        ! Read a single line from the temporary git status file
        read( unit = git_status_file_unit, fmt = '(A)', iostat = ios) single_line
        ! If we've reached the end of the file, stop reading.
        if (ios < 0) exit
        ! Check if the temporary git status file mentions any uncommitted changes
        if (single_line == 'Changes not staged for commit:') has_uncommitted_changes = .true.
    end do

    close( unit = git_status_file_unit)

    ! Mention uncommitted changes in the commit hash (done after writing the commit hash to the terminal,
    ! but still useful for the version that ends up in the NetCDF output files)
    if (has_uncommitted_changes) git_commit_hash = trim( git_commit_hash) // ' (with uncommitted changes!)'

    ! Delete the temporary git status file
    call system( 'rm -f ' // trim( filename_git_status), ierr)
    if (ierr /= 0) call crash('failed to delete temporary git status file')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_for_uncommitted_changes

  subroutine generate_procedural_output_dir_name( output_dir)
    ! Generate a procedural output directory for the current date (e.g. results_20210721_001)
    ! Keep increasing the counter at the end until a directory is available.

    ! In/output variables:
    character(len=*), intent(out) :: output_dir

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'generate_procedural_output_dir_name'
    integer, dimension(8)          :: values
    logical                        :: ex

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    output_dir = ' '

    ! Get current date and time
    call date_and_time( values = values)

    ! Get proper year (assume we're still in the 21st century...)
    output_dir(1:10) = 'results_20'
    select case( floor( real( values(1)) / 10._dp) - 200)
    case(0)
      output_dir(11:11) = '0'
    case(1)
      output_dir(11:11) = '1'
    case(2)
      output_dir(11:11) = '2'
    case(3)
      output_dir(11:11) = '3'
    case(4)
      output_dir(11:11) = '4'
    case(5)
      output_dir(11:11) = '5'
    case(6)
      output_dir(11:11) = '6'
    case(7)
      output_dir(11:11) = '7'
    case(8)
      output_dir(11:11) = '8'
    case(9)
      output_dir(11:11) = '9'
    case default
      call crash('error retrieving date and time!')
    end select

    select case( mod(values(1), 10))
    case(0)
      output_dir(12:12) = '0'
    case(1)
      output_dir(12:12) = '1'
    case(2)
      output_dir(12:12) = '2'
    case(3)
      output_dir(12:12) = '3'
    case(4)
      output_dir(12:12) = '4'
    case(5)
      output_dir(12:12) = '5'
    case(6)
      output_dir(12:12) = '6'
    case(7)
      output_dir(12:12) = '7'
    case(8)
      output_dir(12:12) = '8'
    case(9)
      output_dir(12:12) = '9'
    case default
      call crash('error retrieving date and time!')
    end select

    select case( values(2))
    case(1)
      output_dir(13:14) = '01'
    case(2)
      output_dir(13:14) = '02'
    case(3)
      output_dir(13:14) = '03'
    case(4)
      output_dir(13:14) = '04'
    case(5)
      output_dir(13:14) = '05'
    case(6)
      output_dir(13:14) = '06'
    case(7)
      output_dir(13:14) = '07'
    case(8)
      output_dir(13:14) = '08'
    case(9)
      output_dir(13:14) = '09'
    case(10)
      output_dir(13:14) = '10'
    case(11)
      output_dir(13:14) = '11'
    case(12)
      output_dir(13:14) = '12'
    case default
      call crash('error retrieving date and time!')
    end select

    select case( floor( real( values(3)) / 10._dp))
    case(0)
      output_dir(15:15) = '0'
    case(1)
      output_dir(15:15) = '1'
    case(2)
      output_dir(15:15) = '2'
    case(3)
      output_dir(15:15) = '3'
    case default
      call crash('error retrieving date and time!')
    end select

    select case( MOD(values(3),10))
    case(0)
      output_dir(16:16) = '0'
    case(1)
      output_dir(16:16) = '1'
    case(2)
      output_dir(16:16) = '2'
    case(3)
      output_dir(16:16) = '3'
    case(4)
      output_dir(16:16) = '4'
    case(5)
      output_dir(16:16) = '5'
    case(6)
      output_dir(16:16) = '6'
    case(7)
      output_dir(16:16) = '7'
    case(8)
      output_dir(16:16) = '8'
    case(9)
      output_dir(16:16) = '9'
    case default
      call crash('error retrieving date and time!')
    end select

    output_dir(17:20) = '_001'

    inquire( file = trim( output_dir) // '/.', exist = ex)

    do while (ex)

     if     (output_dir(20:20) == '0') then
       output_dir(20:20) = '1'
     elseif (output_dir(20:20) == '1') then
       output_dir(20:20) = '2'
     elseif (output_dir(20:20) == '2') then
       output_dir(20:20) = '3'
     elseif (output_dir(20:20) == '3') then
       output_dir(20:20) = '4'
     elseif (output_dir(20:20) == '4') then
       output_dir(20:20) = '5'
     elseif (output_dir(20:20) == '5') then
       output_dir(20:20) = '6'
     elseif (output_dir(20:20) == '6') then
       output_dir(20:20) = '7'
     elseif (output_dir(20:20) == '7') then
       output_dir(20:20) = '8'
     elseif (output_dir(20:20) == '8') then
       output_dir(20:20) = '9'
     elseif (output_dir(20:20) == '9') then
       output_dir(20:20) = '0'

       if     (output_dir(19:19) == '0') then
         output_dir(19:19) = '1'
       elseif (output_dir(19:19) == '1') then
         output_dir(19:19) = '2'
       elseif (output_dir(19:19) == '2') then
         output_dir(19:19) = '3'
       elseif (output_dir(19:19) == '3') then
         output_dir(19:19) = '4'
       elseif (output_dir(19:19) == '4') then
         output_dir(19:19) = '5'
       elseif (output_dir(19:19) == '5') then
         output_dir(19:19) = '6'
       elseif (output_dir(19:19) == '6') then
         output_dir(19:19) = '7'
       elseif (output_dir(19:19) == '7') then
         output_dir(19:19) = '8'
       elseif (output_dir(19:19) == '8') then
         output_dir(19:19) = '9'
       elseif (output_dir(19:19) == '9') then
         output_dir(19:19) = '0'

         if     (output_dir(18:18) == '0') then
           output_dir(18:18) = '1'
         elseif (output_dir(18:18) == '1') then
           output_dir(18:18) = '2'
         elseif (output_dir(18:18) == '2') then
           output_dir(18:18) = '3'
         elseif (output_dir(18:18) == '3') then
           output_dir(18:18) = '4'
         elseif (output_dir(18:18) == '4') then
           output_dir(18:18) = '5'
         elseif (output_dir(18:18) == '5') then
           output_dir(18:18) = '6'
         elseif (output_dir(18:18) == '6') then
           output_dir(18:18) = '7'
         elseif (output_dir(18:18) == '7') then
           output_dir(18:18) = '8'
         elseif (output_dir(18:18) == '8') then
           output_dir(18:18) = '9'
         elseif (output_dir(18:18) == '9') then
           output_dir(18:18) = '0'
         end if

       end if

     end if

     inquire( file = trim( output_dir) // '/.', exist = ex)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine generate_procedural_output_dir_name

  subroutine print_model_start

    ! Local variables:
    character(len=1024) :: str1, str2
    integer             :: i

    str1 = ' '
    if (par%n_nodes == 1) then
      str1 = '===== Running ' // trim( program_name) // ' on {int_01} cores ====='
      str1 = insert_val_into_string_int( str1, '{int_01}', par%n)
    else
      str1 = '===== Running ' // trim( program_name) // ' on {int_01} cores ({int_02} nodes) ====='
      str1 = insert_val_into_string_int( str1, '{int_01}', par%n)
      str1 = insert_val_into_string_int( str1, '{int_02}', par%n_nodes)
    end if

    str2 = ' '
    do i = 1, len_trim( str1)
      str2( i:i) = '='
    end do

    if (par%primary) then
      write(0,'(A)') ''
      write(0,'(A)') trim( colour_string( trim( str2),'green'))
      write(0,'(A)') trim( colour_string( trim( str1),'green'))
      write(0,'(A)') trim( colour_string( trim( str2),'green'))
    end if
    call sync

  end subroutine print_model_start

  subroutine print_model_end( tcomp)

    ! In/output variables:
    real(dp), intent(in) :: tcomp

    ! Local variables:
    character(len=1024) :: str1, str2
    integer             :: n,i
    integer             :: nr, ns, nm, nh, nd

    ! Calculate number of elapsed days, hours, minutes, and seconds since this run started
    ns = ceiling( tcomp)

    nr = mod( ns, 60*60*24)
    nd = (ns - nr) / (60*60*24)
    ns = ns - (nd*60*60*24)

    nr = mod( ns, 60*60)
    nh = (ns - nr) / (60*60)
    ns = ns - (nh*60*60)

    nr = mod( ns, 60)
    nm = (ns - nr) / (60)
    ns = ns - (nm*60)

    ! Print to screen
    str1 = '===== Finished running ' // trim( program_name) // &
      ' in {int_01} days, {int_02} hours, {int_03} minutes, and {int_04} seconds ====='
    str1 = insert_val_into_string_int( str1, '{int_01}', nd)
    str1 = insert_val_into_string_int( str1, '{int_02}', nh)
    str1 = insert_val_into_string_int( str1, '{int_03}', nm)
    str1 = insert_val_into_string_int( str1, '{int_04}', ns)

    n = len_trim( str1)
    str2 = ' '
    do i = 1, n
      str2( i:i) = '='
    end do

    if (par%primary) write(0,'(A)') ''
    if (par%primary) write(0,'(A)') trim( colour_string( trim( str2),'green'))
    if (par%primary) write(0,'(A)') trim( colour_string( trim( str1),'green'))
    if (par%primary) write(0,'(A)') trim( colour_string( trim( str2),'green'))
    if (par%primary) write(0,'(A)') ''
    call sync

  end subroutine print_model_end

end module basic_model_utilities