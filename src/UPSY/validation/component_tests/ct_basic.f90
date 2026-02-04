module ct_basic

  use mpi_basic, only: par, sync
  use UPSY_main, only: UPSY
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine

  implicit none

  private

  public :: create_component_tests_output_folder

contains

  !> Create the component test output folder
  subroutine create_component_tests_output_folder( output_dir)

    ! In/output variables:
    character(len=*), intent(in) :: output_dir

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_component_tests_output_folder'
    logical                        :: ex
    character(len=1024)            :: cwd

    ! Add routine to path
    call init_routine( routine_name)

    ! Create the directory
    if (par%primary) then

      ! Create existing folder if necessary
      inquire( file = trim( output_dir) // '/.', exist = ex)
      if (.not. ex) then
        call system('mkdir ' // trim( output_dir))
      end if

      ! Tell the user where it is
      call getcwd( cwd)
      write(0,'(A)') ''
      write(0,'(A)') ' Output directory: ' // UPSY%stru%colour_string( trim(cwd)//'/'//trim( output_dir), 'light blue')
      write(0,'(A)') ''

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_component_tests_output_folder

end module ct_basic