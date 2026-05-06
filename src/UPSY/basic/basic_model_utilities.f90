module basic_model_utilities

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use basic_program_info, only: program_name
  use string_module, only: colour_string, insert_val_into_string_dp, insert_val_into_string_int
  use crash_mod, only: crash
  use mpi_f08, only: MPI_BCAST, MPI_INTEGER, MPI_CHAR, MPI_COMM_WORLD

  implicit none

  private

  public :: generate_procedural_output_dir_name
  public :: print_model_start
  public :: print_model_end
  public :: list_files_in_folder
  public :: get_current_date_time_str

contains

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

  subroutine list_files_in_folder( foldername, list_of_filenames)

    ! In/output variables:
    character(len=*),                            intent(in   ) :: foldername
    character(len=*), dimension(:), allocatable, intent(inout) :: list_of_filenames

    ! Local variables:
    integer             :: funit, n_files, ios, i, ierr
    character(len=1024) :: str

    if (par%primary) then

      ! Create a small text file listing all the files in the directory
      call system('ls ' // trim( foldername) // ' > ' // trim( foldername) // '/list_of_files.txt')

      ! Count how many there are
      open( newunit = funit, file = trim( foldername) // '/list_of_files.txt', action = 'read')
      n_files = 0
      do while (.true.)
        read( funit, fmt = '(a)', iostat = ios) str
        if (ios /= 0) exit
        if (str /= 'list_of_files.txt') n_files = n_files+1
      end do
      close( funit)

    end if
    call MPI_BCAST( n_files, 1, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    allocate( list_of_filenames( n_files))

    if (par%primary) then

      ! List all the files
      open( newunit = funit, file = trim( foldername) // '/list_of_files.txt', action = 'read')
      i = 0
      do while (.true.)
        read( funit, fmt = '(a)', iostat = ios) str
        if (ios /= 0) exit
        if (str /= 'list_of_files.txt') then
          i = i+1
          list_of_filenames( i) = trim( str)
        end if
      end do
      close( funit)

      ! Delete list_of_files.txt
      call system('rm -f ' // trim( foldername) // '/list_of_files.txt')

    end if

    do i = 1, n_files
      call MPI_BCAST( list_of_filenames( i), len( list_of_filenames( i)), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)
    end do

  end subroutine list_files_in_folder

  function get_current_date_time_str() result( datetime_str)
    character(len=8)               :: date_str
    character(len=10)              :: time_str
    character(len=5)               :: zone_str
    character(len=19)              :: datetime_str
    integer, dimension(8)          :: values

    ! Get date and time
    call date_and_time(date=date_str, time=time_str, zone=zone_str, values=values)
    ! Using internal write to format nicely
    write(datetime_str, '(I4.4, "-", I2.2, "-", I2.2, "T", I2.2, ":", I2.2, ":", I2.2)') &
         values(1), values(2), values(3), values(5), values(6), values(7)

  end function get_current_date_time_str

end module basic_model_utilities
