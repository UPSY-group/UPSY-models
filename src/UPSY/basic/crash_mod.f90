module crash_mod

  use precisions, only: dp
  use mpi_basic, only: par
  use string_module, only: insert_val_into_string_int, insert_val_into_string_dp, colour_string
  use basic_program_info, only: routine_path
  use mpi_f08, only: MPI_ABORT, MPI_COMM_WORLD

  implicit none

  private

  public :: crash, warning, happy

contains

  subroutine crash( err_msg, &
    int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
     dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)
    !< Crash the model, write the error message to the screen

    character(len=*),   intent(in) :: err_msg
    integer,  optional, intent(in) :: int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10
    real(dp), optional, intent(in) ::  dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10

    character(len=*), parameter :: print_colour     = 'red'
    character(len=*), parameter :: prefix           = 'ERROR'
    logical,          parameter :: do_crash_program = .true.

    call crash_warning_happy( err_msg, prefix, print_colour, do_crash_program, &
      int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
      dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)

  end subroutine crash

  subroutine warning( err_msg, &
    int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
     dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)
    !< Write the warning message to the screen, but don't crash the model

    character(len=*),   intent(in) :: err_msg
    integer,  optional, intent(in) :: int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10
    real(dp), optional, intent(in) ::  dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10

    character(len=*), parameter :: print_colour     = 'yellow'
    character(len=*), parameter :: prefix           = 'WARNING'
    logical,          parameter :: do_crash_program = .false.

    call crash_warning_happy( err_msg, prefix, print_colour, do_crash_program, &
      int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
      dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)

  end subroutine warning

  subroutine happy( err_msg, &
    int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
    dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)
    !< Write a happy message to the screen

    character(len=*),   intent(in) :: err_msg
    integer,  optional, intent(in) :: int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10
    real(dp), optional, intent(in) ::  dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10

    character(len=*), parameter :: print_colour     = 'yellow'
    character(len=*), parameter :: prefix           = 'SUCCESS'
    logical,          parameter :: do_crash_program = .false.

    call crash_warning_happy( err_msg, prefix, print_colour, do_crash_program, &
      int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
      dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)

  end subroutine happy

  subroutine crash_warning_happy( err_msg, prefix, print_colour, do_crash_program, &
    int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10, &
     dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10)
    !< Print an error message (including inserted numbers) to the terminal
    !< in an appropriate colour, and optionally crash the program

    character(len=*),   intent(in) :: err_msg
    character(len=*),   intent(in) :: prefix
    character(len=*),   intent(in) :: print_colour
    logical,            intent(in) :: do_crash_program
    integer,  optional, intent(in) :: int_01, int_02, int_03, int_04, int_05, int_06, int_07, int_08, int_09, int_10
    real(dp), optional, intent(in) ::  dp_01,  dp_02,  dp_03,  dp_04,  dp_05,  dp_06,  dp_07,  dp_08,  dp_09,  dp_10

    ! Local variables:
    character(len=:), allocatable :: err_msg_loc
    integer                       :: nc
    character(len=9)              :: fmt
    character(len=:), allocatable :: process_str
    integer                       :: ierr

    ! Get local, edit-able copy of error message string
    err_msg_loc = err_msg

    ! Set the process string (e.g. "05/16")
    if     (par%n < 10) then
      nc = 1
    elseif (par%n < 100) then
      nc = 2
    elseif (par%n < 1000) then
      nc = 3
    elseif (par%n < 10000) then
      nc = 4
    else
      nc = 5
    end if

    write( fmt,'(A,I1,A,I1,A)') '(I', nc, ',A,I', nc, ')'
    allocate( character( 2*nc+1) :: process_str)
    write( process_str,fmt) par%i, '/', par%n

    ! Insert numbers into string if needed
    if (present( int_01)) err_msg_loc = insert_val_into_string_int( err_msg_loc, '{int_01}', int_01)
    if (present( int_02)) err_msg_loc = insert_val_into_string_int( err_msg_loc, '{int_02}', int_02)
    if (present( int_03)) err_msg_loc = insert_val_into_string_int( err_msg_loc, '{int_03}', int_03)
    if (present( int_04)) err_msg_loc = insert_val_into_string_int( err_msg_loc, '{int_04}', int_04)
    if (present( int_05)) err_msg_loc = insert_val_into_string_int( err_msg_loc, '{int_05}', int_05)
    if (present( int_06)) err_msg_loc = insert_val_into_string_int( err_msg_loc, '{int_06}', int_06)
    if (present( int_07)) err_msg_loc = insert_val_into_string_int( err_msg_loc, '{int_07}', int_07)
    if (present( int_08)) err_msg_loc = insert_val_into_string_int( err_msg_loc, '{int_08}', int_08)
    if (present( int_09)) err_msg_loc = insert_val_into_string_int( err_msg_loc, '{int_09}', int_09)
    if (present( int_10)) err_msg_loc = insert_val_into_string_int( err_msg_loc, '{int_10}', int_10)

    if (present( dp_01 )) err_msg_loc = insert_val_into_string_dp(  err_msg_loc, '{dp_01}' , dp_01 )
    if (present( dp_02 )) err_msg_loc = insert_val_into_string_dp(  err_msg_loc, '{dp_02}' , dp_02 )
    if (present( dp_03 )) err_msg_loc = insert_val_into_string_dp(  err_msg_loc, '{dp_03}' , dp_03 )
    if (present( dp_04 )) err_msg_loc = insert_val_into_string_dp(  err_msg_loc, '{dp_04}' , dp_04 )
    if (present( dp_05 )) err_msg_loc = insert_val_into_string_dp(  err_msg_loc, '{dp_05}' , dp_05 )
    if (present( dp_06 )) err_msg_loc = insert_val_into_string_dp(  err_msg_loc, '{dp_06}' , dp_06 )
    if (present( dp_07 )) err_msg_loc = insert_val_into_string_dp(  err_msg_loc, '{dp_07}' , dp_07 )
    if (present( dp_08 )) err_msg_loc = insert_val_into_string_dp(  err_msg_loc, '{dp_08}' , dp_08 )
    if (present( dp_09 )) err_msg_loc = insert_val_into_string_dp(  err_msg_loc, '{dp_09}' , dp_09 )
    if (present( dp_10 )) err_msg_loc = insert_val_into_string_dp(  err_msg_loc, '{dp_10}' , dp_10 )

    ! Write the error to the screen
    write(0,'(A,A,A,A,A,A)') colour_string(' ' // prefix // ': ' // trim( err_msg_loc), print_colour) // &
      ' in ' // colour_string( trim( routine_path),'light blue') // &
      ' on process ', colour_string( process_str,'light blue'), ' (0 = primary)'

    ! Stop the program
    if (do_crash_program) error stop

  end subroutine crash_warning_happy

end module crash_mod