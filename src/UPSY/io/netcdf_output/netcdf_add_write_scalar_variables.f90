module netcdf_add_write_scalar_variables
  !< Add and write data to scalar variables

  use mpi_basic, only: par
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use netcdf_basic
  use netcdf, only: NF90_MAX_VAR_DIMS, NF90_INT, NF90_DOUBLE
  use calendar, only: convert_time_to_days

  implicit none

  private

  public :: add_field_dp_0D, add_field_int_0D, write_to_field_multopt_int_0D, &
    write_to_field_multopt_dp_0D, write_time_to_file, write_cftime_to_file

contains

  subroutine add_field_dp_0D( filename, ncid, var_name, long_name, units)
    !< Add a 0-D variable to an existing NetCDF file

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_dp_0D'
    integer                        :: id_dim_time, id_var

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time)

    ! Safety
    if (id_dim_time == -1) call crash('no time dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_DOUBLE, (/ id_dim_time /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_dp_0D

  subroutine add_field_int_0D( filename, ncid, var_name, long_name, units)
    !< Add a 0-D variable to an existing NetCDF file

    ! In/output variables:
    character(len=*),           intent(in   ) :: filename
    integer,                    intent(in   ) :: ncid
    character(len=*),           intent(in   ) :: var_name
    character(len=*), optional, intent(in   ) :: long_name
    character(len=*), optional, intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_int_0D'
    integer                        :: id_dim_time, id_var

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_time  , id_dim_time)

    ! Safety
    if (id_dim_time == -1) call crash('no time dimension could be found in file "' // trim( filename) // '"!')

    ! Create variable
    call create_variable( filename, ncid, var_name, NF90_INT, (/ id_dim_time /), id_var)

    ! Add attributes
    if (present( long_name)) call add_attribute_char( filename, ncid, id_var, 'long_name', long_name)
    if (present( units    )) call add_attribute_char( filename, ncid, id_var, 'units'    , units    )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_field_int_0D

  subroutine write_to_field_multopt_int_0D( filename, ncid, field_name_options, d)
    !< Write a 0-D data field to a NetCDF file variable

    ! Write to the last time frame of the variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: field_name_options
    integer,          intent(in   ) :: d

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'write_to_field_multopt_int_0D'
    integer                                :: id_var, id_dim_time, ti
    character(len=1024)                    :: var_name
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if the file has a time dimension and variable
    call check_time( filename, ncid)
#endif

    ! Inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Inquire file time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

    ! Check if the variable has time as a dimension
    if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, (/ d /), start = (/ ti /), count = (/ 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_int_0D

  subroutine write_to_field_multopt_dp_0D( filename, ncid, field_name_options, d)
    !< Write a 0-D data field to a NetCDF file variable

    ! Write to the last time frame of the variable

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: field_name_options
    real(dp),         intent(in   ) :: d

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'write_to_field_multopt_dp_0D'
    integer                                :: id_var, id_dim_time, ti
    character(len=1024)                    :: var_name
    integer                                :: var_type
    integer                                :: ndims_of_var
    integer, dimension( NF90_MAX_VAR_DIMS) :: dims_of_var

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire the variable
    call inquire_var_multopt( filename, ncid, field_name_options, id_var, var_name = var_name)
    if (id_var == -1) call crash('no variables for name options "' // trim( field_name_options) // '" were found in file "' // trim( filename) // '"!')

#if (DO_ASSERTIONS)
    ! Check if the file has a time dimension and variable
    call check_time( filename, ncid)
#endif

    ! Inquire variable info
    call inquire_var_info( filename, ncid, id_var, var_type = var_type, ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)

    ! Inquire file time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)

    ! Check if the variable has time as a dimension
    if (ndims_of_var /= 1) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" has {int_01} dimensions!', int_01 = ndims_of_var)
    if (.not. any( dims_of_var == id_dim_time)) call crash('variable "' // trim( var_name) // '" in file "' // trim( filename) // '" does not have time as a dimension!')

    ! Inquire length of time dimension
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

    ! Write data to the variable
    call write_var_primary( filename, ncid, id_var, (/ d /), start = (/ ti /), count = (/ 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_field_multopt_dp_0D

  subroutine write_time_to_file( filename, ncid, time)
    ! Write new time value to file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    real(dp),         intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_time_to_file'
    integer                        :: id_dim_time
    integer                        :: id_var_time
    integer                        :: nt

    ! Add routine to path
    call init_routine( routine_name)

#if (DO_ASSERTIONS)
    ! Check time dimension and variable
    call check_time( filename, ncid)
#endif

    ! Determine current length of time dimension in file
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)

    ! Inquire variable id
    call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

    ! Write time
    nt = nt + 1
    call write_var_primary( filename, ncid, id_var_time, (/ time /), start = (/ nt /), count = (/ 1 /) )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_time_to_file

  subroutine write_cftime_to_file( filename, ncid, time, with_bounds)
    ! Write new CFtime value to file

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    real(dp),         intent(in   ) :: time
    logical,          intent(in   ) :: with_bounds

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_cftime_to_file'
    integer                        :: id_dim_time
    integer                        :: id_var_time, id_var_time_bnds
    integer                        :: nt
    real(dp)                       :: days
    real(dp), dimension(2)         :: days_bnds
    real(dp), dimension(2,1)       :: days_bnds_with_time

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine current length of time dimension in file
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)

    ! Inquire variable ids
    call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

    ! Determine next time index
    nt = nt + 1

    ! Convert time to cftime
    if (with_bounds) then
      ! Days centered at 1 July of last year, days_bnds span full last year from 1 Jan to 1 Jan

      ! Convert time
      call convert_time_to_days( time, days, days_bnds)

      ! Write time
      call write_var_primary( filename, ncid, id_var_time, (/ days /), start = (/ nt /), count = (/ 1 /) )

      ! Add "pretend" time dimension to days_bnds
      days_bnds_with_time( :,1) = days_bnds

      ! Write time bounds
      call inquire_var_multopt( filename, ncid, 'time_bnds', id_var_time_bnds)
      call write_var_primary( filename, ncid, id_var_time_bnds, days_bnds_with_time, start = (/ 1, nt /), count = (/ 2, 1 /) )

    else
      ! Days centered at 1 January, no bounds

      ! Convert time
      call convert_time_to_days( time, days)

      ! Write time
      call write_var_primary( filename, ncid, id_var_time, (/ days /), start = (/ nt /), count = (/ 1 /) )

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_cftime_to_file

end module netcdf_add_write_scalar_variables
