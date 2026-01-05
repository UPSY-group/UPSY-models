module fields_init_field

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use fields_basic, only: &
    atype_field, atype_field_2D, atype_field_3D, &
    type_field_logical_2D, type_field_int_2D, type_field_dp_2D, &
    type_field_logical_3D, type_field_int_3D, type_field_dp_3D
  use Arakawa_grid_mod, only: type_Arakawa_grid, Arakawa_grid
  use fields_dimensions, only: type_third_dimension
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: initialise_field

  interface initialise_field
    !< Returns an allocated instance of any of the concrete type_field_TYPE_DIM types,
    !< with its metadata and grid set, and its array pointer associated
    !< with the actual data field.
    procedure :: initialise_field_logical_2D
    procedure :: initialise_field_logical_3D
    procedure :: initialise_field_int_2D
    procedure :: initialise_field_int_3D
    procedure :: initialise_field_dp_2D
    procedure :: initialise_field_dp_3D
  end interface initialise_field

contains

  subroutine initialise_field_logical_2D( field, d, w, field_grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    class(atype_field), allocatable,            intent(  out) :: field
    logical, dimension(:), contiguous, pointer, intent(in   ) :: d
    type(MPI_WIN),                              intent(in   ) :: w
    class(*), target,                           intent(in   ) :: field_grid
    type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
    character(len=*),                           intent(in   ) :: name
    character(len=*),                           intent(in   ) :: long_name
    character(len=*),                           intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_field_logical_2D'
    integer                        :: lb, ub

    ! Add routine to call stack
    call init_routine( routine_name)

    allocate( type_field_logical_2D :: field)

    call initialise_field_meta_grid( field, field_grid, &
      field_Arakawa_grid, name, long_name, units)

    select type(f => field)
    class default
      call crash('programming error?')
    type is (type_field_logical_2D)
      lb = lbound( d,1)
      ub = ubound( d,1)
      f%d( lb: ub) => d
      f%w = w
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_field_logical_2D

  subroutine initialise_field_logical_3D( field, d, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units)

    ! In/output variables:
    class(atype_field), allocatable,              intent(  out) :: field
    logical, dimension(:,:), contiguous, pointer, intent(in   ) :: d
    type(MPI_WIN),                                intent(in   ) :: w
    class(*), target,                             intent(in   ) :: field_grid
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    character(len=*),                             intent(in   ) :: name
    character(len=*),                             intent(in   ) :: long_name
    character(len=*),                             intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_field_logical_3D'
    integer                        :: lb1, ub1, lb2, ub2

    ! Add routine to call stack
    call init_routine( routine_name)

    allocate( type_field_logical_3D :: field)

    call initialise_field_meta_grid( field, field_grid, &
      field_Arakawa_grid, name, long_name, units)

    select type(f => field)
    class default
      call crash('programming error?')
    type is (type_field_logical_3D)
      call f%set_third_dimension( field_third_dimension)
      lb1 = lbound( d,1)
      ub1 = ubound( d,1)
      lb2 = lbound( d,2)
      ub2 = ubound( d,2)
      f%d( lb1: ub1, lb2: ub2) => d
      f%w = w
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_field_logical_3D

  subroutine initialise_field_int_2D( field, d, w, field_grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    class(atype_field), allocatable,             intent(  out) :: field
    integer, dimension(:), contiguous, pointer,  intent(in   ) :: d
    type(MPI_WIN),                               intent(in   ) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                            intent(in   ) :: name
    character(len=*),                            intent(in   ) :: long_name
    character(len=*),                            intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_field_int_2D'
    integer                        :: lb, ub

    ! Add routine to call stack
    call init_routine( routine_name)

    allocate( type_field_int_2D :: field)

    call initialise_field_meta_grid( field, field_grid, &
      field_Arakawa_grid, name, long_name, units)

    select type(f => field)
    class default
      call crash('programming error?')
    type is (type_field_int_2D)
      lb = lbound( d,1)
      ub = ubound( d,1)
      f%d( lb: ub) => d
      f%w = w
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_field_int_2D

  subroutine initialise_field_int_3D( field, d, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units)

    ! In/output variables:
    class(atype_field), allocatable,              intent(  out) :: field
    integer, dimension(:,:), contiguous, pointer, intent(in   ) :: d
    type(MPI_WIN),                                intent(in   ) :: w
    class(*), target,                             intent(in   ) :: field_grid
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    character(len=*),                             intent(in   ) :: name
    character(len=*),                             intent(in   ) :: long_name
    character(len=*),                             intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_field_int_3D'
    integer                        :: lb1, ub1, lb2, ub2

    ! Add routine to call stack
    call init_routine( routine_name)

    allocate( type_field_int_3D :: field)

    call initialise_field_meta_grid( field, field_grid, &
      field_Arakawa_grid, name, long_name, units)

    select type(f => field)
    class default
      call crash('programming error?')
    type is (type_field_int_3D)
      call f%set_third_dimension( field_third_dimension)
      lb1 = lbound( d,1)
      ub1 = ubound( d,1)
      lb2 = lbound( d,2)
      ub2 = ubound( d,2)
      f%d( lb1: ub1, lb2: ub2) => d
      f%w = w
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_field_int_3D

  subroutine initialise_field_dp_2D( field, d, w, field_grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    class(atype_field), allocatable,             intent(  out) :: field
    real(dp), dimension(:), contiguous, pointer, intent(in   ) :: d
    type(MPI_WIN),                               intent(in   ) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                            intent(in   ) :: name
    character(len=*),                            intent(in   ) :: long_name
    character(len=*),                            intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_field_dp_2D'
    integer                        :: lb, ub

    ! Add routine to call stack
    call init_routine( routine_name)

    allocate( type_field_dp_2D :: field)

    call initialise_field_meta_grid( field, field_grid, &
      field_Arakawa_grid, name, long_name, units)

    select type(f => field)
    class default
      call crash('programming error?')
    type is (type_field_dp_2D)
      lb = lbound( d,1)
      ub = ubound( d,1)
      f%d( lb: ub) => d
      f%w = w
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_field_dp_2D

  subroutine initialise_field_dp_3D( field, d, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units)

    ! In/output variables:
    class(atype_field), allocatable,               intent(  out) :: field
    real(dp), dimension(:,:), contiguous, pointer, intent(in   ) :: d
    type(MPI_WIN),                                 intent(in   ) :: w
    class(*), target,                              intent(in   ) :: field_grid
    type(type_Arakawa_grid),                       intent(in   ) :: field_Arakawa_grid
    type(type_third_dimension),                    intent(in   ) :: field_third_dimension
    character(len=*),                              intent(in   ) :: name
    character(len=*),                              intent(in   ) :: long_name
    character(len=*),                              intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_field_dp_3D'
    integer                        :: lb1, ub1, lb2, ub2

    ! Add routine to call stack
    call init_routine( routine_name)

    allocate( type_field_dp_3D :: field)

    call initialise_field_meta_grid( field, field_grid, &
      field_Arakawa_grid, name, long_name, units)

    select type(f => field)
    class default
      call crash('programming error?')
    type is (type_field_dp_3D)
      call f%set_third_dimension( field_third_dimension)
      lb1 = lbound( d,1)
      ub1 = ubound( d,1)
      lb2 = lbound( d,2)
      ub2 = ubound( d,2)
      f%d( lb1: ub1, lb2: ub2) => d
      f%w = w
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_field_dp_3D

  subroutine initialise_field_meta_grid( field, field_grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    class(atype_field),      intent(inout) :: field
    class(*), target,        intent(in   ) :: field_grid
    type(type_Arakawa_grid), intent(in   ) :: field_Arakawa_grid
    character(len=*),        intent(in   ) :: name
    character(len=*),        intent(in   ) :: long_name
    character(len=*),        intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_field_meta_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Metadata
    call field%set_name     ( name)
    call field%set_long_name( long_name)
    call field%set_units    ( units)

    ! Grid
    call field%set_grid        ( field_grid)
    call field%set_Arakawa_grid( field_Arakawa_grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_field_meta_grid

end module fields_init_field