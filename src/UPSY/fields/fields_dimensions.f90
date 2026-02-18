module fields_dimensions

  use precisions, only: dp
  use parameters, only: NaN
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_zeta, only: initialise_scaled_vertical_coordinate_regular, &
    initialise_scaled_vertical_coordinate_irregular_log, &
    initialise_scaled_vertical_coordinate_old_15_layer
  use netcdf_io_main
  use netcdf, only: NF90_DOUBLE
  use mpi_f08, only: MPI_BCAST, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD

  implicit none

  private

  public :: third_dimension, type_third_dimension

  type type_third_dimension_factory
  contains
    procedure :: ice_zeta    => third_dimension_factory_ice_zeta
    procedure :: month       => third_dimension_factory_month
    procedure :: ocean_depth => third_dimension_factory_ocean_depth
  end type type_third_dimension_factory

  type type_third_dimension
    character(:), allocatable           :: name
    integer                             :: n
    real(dp), dimension(:), allocatable :: val
    character(:), allocatable           :: units
  contains
    generic,   public  :: operator(==) => eq
    procedure, private :: eq => test_third_dimension_equality
    procedure, public  :: create_dim_and_var_in_netcdf
    procedure, public  :: is_dim_and_var_in_netcdf
  end type type_third_dimension

  type(type_third_dimension_factory) :: third_dimension

contains

  function third_dimension_factory_ice_zeta( self, n, choice_zeta_grid, R) result( res)

    ! In/output variables:
    class(type_third_dimension_factory), intent(in) :: self
    integer,                             intent(in) :: n
    character(len=*),                    intent(in) :: choice_zeta_grid
    real(dp), optional,                  intent(in) :: R
    type(type_third_dimension)                      :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_scaled_vertical_coordinate'
    real(dp), dimension(n  )       :: zeta
    real(dp), dimension(n-1)       :: zeta_stag

    ! Add routine to call stack
    call init_routine( routine_name)

    res%name  = 'zeta'
    res%n     = n
    res%units = '-'

    ! Calculate zeta values
    select case (choice_zeta_grid)
    case default
      call crash('unknown choice_zeta_grid "' // trim( choice_zeta_grid) // '"')
    case ('regular')
      call initialise_scaled_vertical_coordinate_regular( n, zeta, zeta_stag)
    case ('irregular_log')
      call initialise_scaled_vertical_coordinate_irregular_log( n, R, zeta, zeta_stag)
    case ('old_15_layer_zeta')
      call initialise_scaled_vertical_coordinate_old_15_layer( n, zeta, zeta_stag)
    end select

    res%val = zeta

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function third_dimension_factory_ice_zeta

  function third_dimension_factory_month( self) result( res)
    class(type_third_dimension_factory), intent(in) :: self
    type(type_third_dimension)                      :: res
    res%name = 'month'
    res%n    = 12
    res%val = [1._dp, 2._dp, 3._dp, 4._dp, 5._dp, 6._dp, &
      7._dp, 8._dp, 9._dp, 10._dp, 11._dp, 12._dp]
    res%units = '-'
  end function third_dimension_factory_month

  function third_dimension_factory_ocean_depth( self, n) result( res)
    class(type_third_dimension_factory), intent(in) :: self
    integer,                             intent(in) :: n
    type(type_third_dimension)                      :: res
    res%name = 'ocean_depth'
    res%n    = n
    ! FIXME
    allocate( res%val( n), source = NaN)
    res%units = 'm'
  end function third_dimension_factory_ocean_depth

  pure function test_third_dimension_equality( dim1, dim2) result( res)
    class(type_third_dimension), intent(in) :: dim1
    class(type_third_dimension), intent(in) :: dim2
    logical                                 :: res
    res = dim1%name == dim2%name
  end function test_third_dimension_equality

  subroutine create_dim_and_var_in_netcdf( dim, filename, ncid)

    ! In/output variables:
    class(type_third_dimension), intent(in) :: dim
    character(len=*),            intent(in) :: filename
    integer,                     intent(in) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_dim_and_var_in_netcdf'
    integer                        :: dim_length, id_dim, id_var

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Check if this dimension is already present in the file; if so, skip it
    call inquire_dim( filename, ncid, trim( dim%name), dim_length, id_dim)
    if (id_dim /= -1) then
      call finalise_routine( routine_name)
      return
    end if

    call create_dimension( filename, ncid, trim( dim%name), dim%n, id_dim)
    call create_variable( filename, ncid, trim( dim%name), NF90_DOUBLE, [id_dim], id_var)
    call write_var_primary( filename, ncid, id_var, dim%val)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_dim_and_var_in_netcdf

  function is_dim_and_var_in_netcdf( dim, filename, ncid) result( res)

    ! In/output variables:
    class(type_third_dimension), intent(in) :: dim
    character(len=*),            intent(in) :: filename
    integer,                     intent(in) :: ncid
    logical                                 :: res

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'is_dim_and_var_in_netcdf'
    integer                             :: dim_length, id_dim, id_var
    real(dp), dimension(:), allocatable :: d_from_file
    integer                             :: ierr
    integer                             :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    res = .false.

    ! Check if the corresponding dimension exists in the file, and is of the correct length
    call inquire_dim( filename, ncid, trim( dim%name), dim_length, id_dim)
    if (id_dim == -1 .or. dim_length /= dim%n) goto 888

    ! Check if the corresponding variable exists in the file
    call inquire_var( filename, ncid, trim( dim%name), id_var)
    if (id_var == -1) goto 888

    ! Read it, and check if it has the right values
    allocate( d_from_file( dim_length))
    call read_var_primary( filename, ncid, id_var, d_from_file)
    call MPI_BCAST( d_from_file, dim_length, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    do i = 1, dim_length
      if (abs( d_from_file( i) / dim%val( i) - 1._dp) > 1e-9_dp) goto 888
    end do

    res = .true.

    ! Remove routine from call stack
888 call finalise_routine( routine_name)

  end function is_dim_and_var_in_netcdf

end module fields_dimensions