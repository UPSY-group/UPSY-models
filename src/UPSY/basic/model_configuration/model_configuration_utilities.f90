module model_configuration_utilities

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, &
    crash, warning, capitalise_string, remove_leading_spaces

  implicit none

  private

  public :: check_config_file_validity

contains

  subroutine check_config_file_validity( config_filename, namelist_filename)

    ! In/output variables:
    character(len=*), intent(in) :: config_filename
    character(len=*), intent(in) :: namelist_filename

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'check_config_file_validity'
    logical                        :: ex, all_are_valid, all_are_present

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    inquire( file = config_filename, exist = ex)
    if (.not. ex) call crash('config file ' // trim( config_filename) // ' could not be found!')
    inquire( file = namelist_filename, exist = ex)
    if (.not. ex) call crash('namelist file ' // trim( namelist_filename) // ' could not be found!')

    ! Check if all the variables appearing in the config file "config_filename" are valid
    call check_if_all_config_variables_are_valid( config_filename, namelist_filename, all_are_valid)

    ! Check if all the expected config variables appear in the config file "config_filename"
    !CALL check_if_all_expected_config_variables_are_present( config_filename, namelist_filename, all_are_present)
    all_are_present = .true.

    ! If not all is well, crash
    if (.not. (all_are_valid .and. all_are_present)) call crash('config file "' // TRIM( config_filename) // '" is invalid!')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_config_file_validity

  subroutine check_if_all_config_variables_are_valid( config_filename, namelist_filename, all_are_valid)
    ! Check if all the variables appearing in the config file "config_filename" are valid
    !
    ! Do this by reading one line at a time of the config file, determining the name of the variable
    ! declared in that line, and checking if that variable also exists in the namelist file.
    !
    ! The namelist file is created earlier by writing the namelist to a text file.

    ! In/output variables:
    character(len=*), intent(in   ) :: config_filename
    character(len=*), intent(in   ) :: namelist_filename
    logical,          intent(  out) :: all_are_valid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'check_if_all_config_variables_are_valid'
    integer                        :: config_unit, namelist_unit
    integer                        :: ios
    logical                        :: found_end_of_file_config, found_end_of_file_namelist
    character(len=1024)            :: single_line_config      , single_line_namelist
    integer                        :: line_counter_config     , line_counter_namelist
    logical                        :: found_match

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the config file
    open( newunit = config_unit, file = config_filename, iostat = ios)
    if (ios /= 0) call crash('couldnt open config file "' // trim( config_filename) // '"!')

    ! Read one line at a time of the config file, determine the name of the variable
    ! declared in that line, and check if that variable also exists in the namelist file

    found_end_of_file_config = .false.
    line_counter_config      = 0
    all_are_valid            = .true.

    do while (.not. found_end_of_file_config)

      line_counter_config = line_counter_config + 1

      ! Read a single line from the config file
      read( unit = config_unit, fmt = '(A)', iostat = ios) single_line_config

      ! If we've reached the end of the file before finding the
      ! terminating forward slash, this config file is not valid.
      if (ios < 0) call crash('config file "' // trim( config_filename) // &
        '" is not terminated with a forward slash!')

      ! Remove all leading spaces
      call remove_leading_spaces( single_line_config)

      ! The variable name is the part of the string left of the first (, =, or space.
      single_line_config = single_line_config( 1: scan( single_line_config, '( =')-1)

      ! Get config variable in all caps for case-insensitive comparison
      call capitalise_string( single_line_config)

      ! The forward slash at the end terminates the config file
      if (single_line_config == '/') then
        found_end_of_file_config = .true.
      end if

      ! Disregard empty lines, commented lines, and the header line
      if (single_line_config == '' .or. &
          single_line_config( 1:1) == '&' .or. &
          single_line_config( 1:1) == '!') then
        cycle
      end if

      ! Open the namelist file
      open( newunit = namelist_unit, file = namelist_filename)
      if (ios /= 0) call crash('couldnt open namelist file "' // trim( namelist_filename) // '"!')

      ! Read all variables from the namelist file and check
      ! if any of them match the current config variable

      found_end_of_file_namelist = .false.
      line_counter_namelist      = 0
      found_match                = .false.

      do while ((.not. found_end_of_file_namelist) .and. (.not. found_match))

        line_counter_namelist = line_counter_namelist + 1

        ! Read a single line from the namelist file
        read( unit = namelist_unit, fmt = '(A)', iostat = ios) single_line_namelist

        ! If we've reached the end of the file before finding the terminating forward slash, this namelist file is not valid.
        if (ios < 0) call crash('namelist file "' // trim( namelist_filename) // &
          '" is not terminated with a forward slash!')

        ! Remove all leading spaces
        call remove_leading_spaces( single_line_namelist)

        ! The variable name is the part of the string left of the first (, =, or space.
        single_line_namelist = single_line_namelist( 1: scan( single_line_namelist, '( =')-1)

        ! Get namelist variable in all caps for case-insensitive comparison
        call capitalise_string( single_line_namelist)

        ! The forward slash at the end terminates the config file
        if (single_line_namelist == '/') found_end_of_file_namelist = .true.

        ! Disregard empty lines, commented lines, and the header line
        if (single_line_namelist == '' .or. &
            single_line_namelist( 1:1) == '&' .or. &
            single_line_namelist( 1:1) == '!') then
          cycle
        end if

        ! Check if this namelist variable matches the config variable
        if (single_line_namelist == single_line_config) found_match = .true.

      end do

      ! If no matching variable was found in the namelist file, print an error
      if (.not. found_match) then
        all_are_valid = .false.
        call warning('invalid config variable "' // trim( single_line_config) // &
        '" in file "' // trim( config_filename) // '", line {int_01}', &
        int_01 = line_counter_config)
      end if

      ! Close the namelist file
      close( unit = namelist_unit)

    end do

    ! Close the config file
    close( unit = config_unit)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_if_all_config_variables_are_valid

  subroutine check_if_all_expected_config_variables_are_present( config_filename, namelist_filename, all_are_present)
    ! Check if all the expected config variables appear in the config file "config_filename"
    !
    ! Do this by reading one line at a time of the namelist file, determining the name of the variable
    ! declared in that line, and checking if that variable also exists in the config file.
    !
    ! The namelist file is created earlier by writing the namelist to a text file.

    ! In/output variables:
    character(len=*), intent(in   ) :: config_filename
    character(len=*), intent(in   ) :: namelist_filename
    logical,          intent(  out) :: all_are_present

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'check_if_all_expected_config_variables_are_present'
    integer                        :: config_unit, namelist_unit
    integer                        :: ios
    logical                        :: found_end_of_file_config, found_end_of_file_namelist
    character(len=1024)            :: single_line_config      , single_line_namelist
    integer                        :: line_counter_config     , line_counter_namelist
    logical                        :: found_match

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the namelist file
    open( newunit = namelist_unit, file = namelist_filename, iostat = ios)
    if (ios /= 0) call crash('couldnt open namelist file "' // trim( namelist_filename) // '"!')

    ! Read one line at a time of the nanmelist file, determine the name of the variable
    ! declared in that line, and check if that variable also exists in the config file

    found_end_of_file_namelist = .false.
    line_counter_namelist      = 0
    all_are_present            = .true.

    do while (.not. found_end_of_file_namelist)

      line_counter_namelist = line_counter_namelist + 1

      ! Read a single line from the config file
      read( unit = namelist_unit, fmt = '(A)', iostat = ios) single_line_namelist

      ! If we've reached the end of the file before finding the
      ! terminating forward slash, this namelist file is not valid.
      if (ios < 0) call crash('namelist file "' // trim( namelist_filename) // &
        '" is not terminated with a forward slash!')

      ! Remove all leading spaces
      call remove_leading_spaces( single_line_namelist)

      ! The variable name is the part of the string left of the first (, =, or space.
      single_line_namelist = single_line_namelist( 1: scan( single_line_namelist, '( =')-1)

      ! Get namelist variable in all caps for case-insensitive comparison
      call capitalise_string( single_line_namelist)

      ! The forward slash at the end terminates the namelist file
      if (single_line_namelist == '/') found_end_of_file_namelist = .true.

      ! Disregard empty lines, commented lines, the header line, and the final line
      if (single_line_namelist == '' .or. &
          single_line_namelist( 1:1) == '&' .or. &
          single_line_namelist( 1:1) == '!') then
        cycle
      end if

      ! Open the config file
      open( newunit = config_unit, file = config_filename)
      if (ios /= 0) call crash('couldnt open config file "' // trim( config_filename) // '"!')

      ! Read all variables from the config file and check if any of them match the current namelist variable

      found_end_of_file_config = .false.
      line_counter_config      = 0
      found_match              = .false.

      do while ((.not. found_end_of_file_config) .and. (.not. found_match))

        line_counter_config = line_counter_config + 1

        ! Read a single line from the config file
        read( unit = config_unit, fmt = '(A)', iostat = ios) single_line_config

        ! If we've reached the end of the file before finding the
        ! terminating forward slash, this config file is not valid.
        if (ios < 0) call crash('config file "' // trim( config_filename) // &
          '" is not terminated with a forward slash!')

        ! Remove all leading spaces
        call remove_leading_spaces( single_line_config)

        ! The variable name is the part of the string left of the first (, =, or space.
        single_line_config = single_line_config( 1: scan( single_line_config, '( =')-1)

        ! Get config variable in all caps for case-insensitive comparison
        call capitalise_string( single_line_config)

        ! The forward slash at the end terminates the config file
        if (single_line_config == '/') found_end_of_file_config = .true.

        ! Disregard empty lines, commented lines, and the header line
        if (single_line_config == '' .or. &
            single_line_config( 1:1) == '&' .or. &
            single_line_config( 1:1) == '!') then
          cycle
        end if

        ! Check if this namelist variable matches the config variable
        if (single_line_config == single_line_namelist) found_match = .true.

      end do

      ! If no matching variable was found in the config file, print an error
      if (.not. found_match) then
        all_are_present = .false.
        call warning('couldnt find config variable "' // trim( single_line_namelist) // &
          '" in file "' // trim( config_filename) // '"')
      end if

      ! Close the namelist file
      close( unit = config_unit)

    end do

    ! Close the config file
    close( unit = namelist_unit)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_if_all_expected_config_variables_are_present

end module model_configuration_utilities
