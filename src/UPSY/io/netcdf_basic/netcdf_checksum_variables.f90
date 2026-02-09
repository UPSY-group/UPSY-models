module netcdf_checksum_variables

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use netcdf_basic_wrappers, only: inquire_dim_info, create_variable, create_scalar_variable, &
    inquire_var
  use netcdf_write_var_primary, only: write_var_primary
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_SUM, &
    MPI_MIN, MPI_MAX, MPI_COMM_WORLD

  implicit none

  private

  public :: create_checksum_variable_set, write_var_checksum_notime

  interface write_var_checksum_notime
    procedure :: write_var_checksum_dp_1D_notime
    procedure :: write_var_checksum_dp_2D_notime
  end interface write_var_checksum_notime

  interface calc_checksums
    procedure :: calc_checksums_dp_1D
    procedure :: calc_checksums_dp_2D
  end interface calc_checksums

contains

  subroutine create_checksum_variable_set( filename, ncid, var_name, var_type, dim_ids)

    ! In/output variables:
    character(len=*),      intent(in   ) :: filename
    integer,               intent(in   ) :: ncid
    character(len=*),      intent(in   ) :: var_name
    integer,               intent(in   ) :: var_type
    integer, dimension(:), intent(in   ) :: dim_ids

    ! Local variables:
    character(len=*), parameter :: routine_name = 'create_checksum_variable_set'
    logical                     :: has_time
    integer                     :: i, id_dim_time
    character(len=1024)         :: dim_name

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Check if the variable has a time dimension
    has_time = .false.
    id_dim_time = -1
    do i = 1, size( dim_ids,1)
      call inquire_dim_info( filename, ncid, dim_ids(i), dim_name = dim_name)
      if (dim_name == 'time') then
        has_time = .true.
        id_dim_time = dim_ids( i)
        exit
      end if
    end do

    if (has_time) then
      call crash('fixme')
      ! call create_checksum_variable_set_time( filename, ncid, var_name, var_type, id_dim_time)
    else
      call create_checksum_variable_set_notime( filename, ncid, var_name, var_type)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_checksum_variable_set

  subroutine create_checksum_variable_set_notime( filename, ncid, var_name, var_type)

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: var_name
    integer,          intent(in   ) :: var_type

    ! Local variables:
    character(len=*), parameter :: routine_name = 'create_checksum_variable_set_notime'
    integer                     :: id_var

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! The four variable versions: sum, sum_abs, min, and max
    call create_scalar_variable( filename, ncid, trim( var_name) // '_sum'    , var_type, id_var)
    call create_scalar_variable( filename, ncid, trim( var_name) // '_sum_abs', var_type, id_var)
    call create_scalar_variable( filename, ncid, trim( var_name) // '_min'    , var_type, id_var)
    call create_scalar_variable( filename, ncid, trim( var_name) // '_max'    , var_type, id_var)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_checksum_variable_set_notime

  subroutine write_var_checksum_dp_1D_notime( filename, ncid, var_name, d)

    ! In/output variables:
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: var_name
    real(dp), dimension(:), intent(in   ) :: d

    ! Local variables:
    character(len=*), parameter :: routine_name = 'write_var_checksum_dp_1D_notime'
    integer                     :: id_var_sum
    integer                     :: id_var_sum_abs
    integer                     :: id_var_min
    integer                     :: id_var_max
    real(dp)                    :: d_sum, d_sum_abs, d_min, d_max

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Find variable IDs
    call inquire_var( filename, ncid, trim( var_name) // '_sum'    , id_var_sum    )
    call inquire_var( filename, ncid, trim( var_name) // '_sum_abs', id_var_sum_abs)
    call inquire_var( filename, ncid, trim( var_name) // '_min'    , id_var_min    )
    call inquire_var( filename, ncid, trim( var_name) // '_max'    , id_var_max    )
    if (id_var_sum == -1 .or. id_var_sum_abs == -1 .or. id_var_min == -1 .or. id_var_max == -1) &
      call crash('couldnt find checksum variables for ' // trim( var_name) // ' in ' // trim( filename))

    ! Calculate checksums
    call calc_checksums( d, d_sum, d_sum_abs, d_min, d_max)

    ! Write to netcdf
    call write_var_primary( filename, ncid, id_var_sum    , d_sum    )
    call write_var_primary( filename, ncid, id_var_sum_abs, d_sum_abs)
    call write_var_primary( filename, ncid, id_var_min    , d_min    )
    call write_var_primary( filename, ncid, id_var_max    , d_max    )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_checksum_dp_1D_notime

  subroutine write_var_checksum_dp_2D_notime( filename, ncid, var_name, d)

    ! In/output variables:
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: var_name
    real(dp), dimension(:,:), intent(in   ) :: d

    ! Local variables:
    character(len=*), parameter :: routine_name = 'write_var_checksum_dp_1D_notime'
    integer                     :: id_var_sum
    integer                     :: id_var_sum_abs
    integer                     :: id_var_min
    integer                     :: id_var_max
    real(dp)                    :: d_sum, d_sum_abs, d_min, d_max

    ! Add routine to path
    call init_routine( routine_name, do_track_resource_use = .false.)

    ! Find variable IDs
    call inquire_var( filename, ncid, trim( var_name) // '_sum'    , id_var_sum    )
    call inquire_var( filename, ncid, trim( var_name) // '_sum_abs', id_var_sum_abs)
    call inquire_var( filename, ncid, trim( var_name) // '_min'    , id_var_min    )
    call inquire_var( filename, ncid, trim( var_name) // '_max'    , id_var_max    )
    if (id_var_sum == -1 .or. id_var_sum_abs == -1 .or. id_var_min == -1 .or. id_var_max == -1) &
      call crash('couldnt find checksum variables for ' // trim( var_name) // ' in ' // trim( filename))

    ! Calculate checksums
    call calc_checksums( d, d_sum, d_sum_abs, d_min, d_max)

    ! Write to netcdf
    call write_var_primary( filename, ncid, id_var_sum    , d_sum    )
    call write_var_primary( filename, ncid, id_var_sum_abs, d_sum_abs)
    call write_var_primary( filename, ncid, id_var_min    , d_min    )
    call write_var_primary( filename, ncid, id_var_max    , d_max    )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_var_checksum_dp_2D_notime

  subroutine calc_checksums_dp_1D( d, d_sum, d_sum_abs, d_min, d_max)
    real(dp), dimension(:), intent(in   ) :: d
    real(dp),               intent(  out) :: d_sum, d_sum_abs, d_min, d_max
    integer                               :: ierr
    d_sum     = sum( d)
    d_sum_abs = sum( abs( d))
    d_min     = minval( d)
    d_max     = maxval( d)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_sum    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_sum_abs, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_min    , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_max    , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  end subroutine calc_checksums_dp_1D

  subroutine calc_checksums_dp_2D( d, d_sum, d_sum_abs, d_min, d_max)
    real(dp), dimension(:,:), intent(in   ) :: d
    real(dp),                 intent(  out) :: d_sum, d_sum_abs, d_min, d_max
    integer                                 :: ierr
    d_sum     = sum( d)
    d_sum_abs = sum( abs( d))
    d_min     = minval( d)
    d_max     = maxval( d)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_sum    , 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_sum_abs, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_min    , 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, d_max    , 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
  end subroutine calc_checksums_dp_2D

end module netcdf_checksum_variables