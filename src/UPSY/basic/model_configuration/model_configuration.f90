module model_configuration

  ! The different parameters that control the simulation.

  use precisions, only: dp
  use mpi_f08, only: MPI_BCAST, MPI_COMM_WORLD, MPI_CHAR, MPI_LOGICAL
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, &
    warning, crash, colour_string
  use basic_model_utilities, only: get_git_commit_hash, git_commit_hash, &
    check_for_uncommitted_changes, has_uncommitted_changes, generate_procedural_output_dir_name
  use model_configuration_type_and_namelist, only: type_config, copy_config_variables_to_struct, &
    read_config_file

  implicit none

  private

  public :: C, initialise_model_configuration, initialise_model_configuration_unit_tests

  ! The main config structure
  type(type_config) :: C

contains

  subroutine initialise_model_configuration( config_filename)

    ! In/output variables:
    character(len=*), intent(in) :: config_filename

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_model_configuration'

    ! Add routine to path
    call init_routine( routine_name)

    ! Figure out which git commit of the model we're running
    call initialise_model_configuration_git_commit

    ! Initialise main config parameters
    call initialise_model_configuration_config( config_filename)

    ! Set up the output directory
    call initialise_model_configuration_output_dir

    ! Copy the config file to the output directory
    if (par%primary) then
      call system('cp ' // config_filename    // ' ' // trim( C%output_dir))
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_model_configuration

  subroutine initialise_model_configuration_git_commit
    ! Figure out which git commit of the model we're running

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_model_configuration_git_commit'
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) call get_git_commit_hash( git_commit_hash)
    call mpi_bcast( git_commit_hash, len( git_commit_hash), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    if (par%primary) write(0,'(A)') ''
    if (par%primary) write(0,'(A)') ' Running UFEMISM from git commit ' // colour_string( trim( git_commit_hash), 'pink')

    if (par%primary) call check_for_uncommitted_changes
    call mpi_bcast( has_uncommitted_changes, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    if (par%primary .and. has_uncommitted_changes) then
      write(0,'(A)') colour_string( ' WARNING: You have uncommitted changes; the current simulation might not be reproducible!', 'yellow')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_model_configuration_git_commit

  subroutine initialise_model_configuration_config( config_filename)
    ! Initialise main config parameters

    ! In/output variables:
    character(len=*), intent(in) :: config_filename

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_model_configuration_config'
    character(len=1024)            :: output_dir_procedural
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) write(0,'(A)') ''
    if (par%primary) write(0,'(A)') ' Running UFEMISM with settings from configuration file: ' // colour_string( TRIM( config_filename), 'light blue')

    ! Initialise the main config structure from the config file
    call initialise_config_from_file( config_filename)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_model_configuration_config

  subroutine initialise_model_configuration_output_dir
    ! Set up the output directory

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_model_configuration_output_dir'
    character(len=1024)            :: output_dir_procedural
    integer                        :: ierr
    logical                        :: ex

    ! Add routine to path
    call init_routine( routine_name)

    ! First get the name of the output directory (either procedural, or provided in the config file)
    C%output_dir = ' '

    if (C%create_procedural_output_dir) then
      ! Automatically create an output directory with a procedural name (e.g. results_20210720_001/)

      if (par%primary) then
        call generate_procedural_output_dir_name( output_dir_procedural)
        C%output_dir( 1: len_trim( output_dir_procedural)+1) = trim( output_dir_procedural) // '/'
      end if
      call MPI_BCAST( C%output_dir, 256, MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    else
      ! Use the provided name (return an error if this directory already exists)

      C%output_dir = trim( C%fixed_output_dir) // '/'

    end if

    ! Create the directory
    if (par%primary) then

      ! Safety
      inquire( file = trim( C%output_dir) // '/.', exist = ex)
      if (ex) then
        call crash('output directory ' // trim( C%output_dir) // ' already exists!')
      end if

      ! Create output directory
      call system('mkdir ' // trim( C%output_dir))

      ! Tell the user where it is
      write(0,'(A)') ''
      write(0,'(A)') ' Output directory: ' // colour_string( trim( C%output_dir), 'light blue')

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_model_configuration_output_dir

  subroutine initialise_model_configuration_unit_tests

    ! Local variables:
    character(len=256), parameter :: routine_name = 'initialise_model_configuration_unit_tests'
    integer                       :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Figure out which git commit of the model we're running
    if (par%primary) call get_git_commit_hash( git_commit_hash)
    call mpi_bcast( git_commit_hash, len( git_commit_hash), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    if (par%primary) write(0,'(A)') ''
    if (par%primary) write(0,'(A)') ' Running UFEMISM from git commit ' // colour_string( trim( git_commit_hash), 'pink')

    if (par%primary) call check_for_uncommitted_changes
    call mpi_bcast( has_uncommitted_changes, 1, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    if (par%primary .and. has_uncommitted_changes) then
      write(0,'(A)') colour_string( ' WARNING: You have uncommitted changes; the current simulation might not be reproducible!', 'yellow')
    end if

    ! Copy values from the XXX_config variables to the config structure
    call copy_config_variables_to_struct( C)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_model_configuration_unit_tests

  subroutine initialise_config_from_file( config_filename)

    ! In/output variables:
    character(len=*), intent(in) :: config_filename

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_config_from_file'
    integer                        :: i

    ! Add routine to path
    call init_routine( routine_name)

    ! Let each of the processors read the config file in turns so there's no access conflicts
    do i = 0, par%n-1

      if (i == par%i) then

        ! Read the external file, use a Fortran NAMELIST to overwrite the default
        ! values of the XXX_config variables
        call read_config_file( config_filename)

        ! Copy values from the XXX_config variables to the config structure
        call copy_config_variables_to_struct( C)

      end if

      ! Make sure only one process at a time reads from / writes to disk
      call sync

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_config_from_file

end module model_configuration
