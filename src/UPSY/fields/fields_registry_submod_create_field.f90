submodule( fields_registry) fields_registry_submod_create_field

contains

  subroutine create_field_logical_2D( flds_reg, d, w, field_grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    class(type_fields_registry),                 intent(inout) :: flds_reg
    logical, dimension(:), contiguous, pointer,  intent(inout) :: d
    type(MPI_WIN),                               intent(inout) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                            intent(in   ) :: name
    character(len=*),                            intent(in   ) :: long_name
    character(len=*),                            intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_logical_2D'
    integer                         :: lb, ub, n
    class(atype_field), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call determine_bounds_2D( field_grid, field_Arakawa_grid, lb, ub, n)
    call allocate_dist_shared( d, w, n)
    d( lb: ub) => d
    call initialise_field( field, d, w, field_grid, field_Arakawa_grid, name, long_name, units)
    call flds_reg%add( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_logical_2D

  subroutine create_field_logical_3D( flds_reg, d, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units)

    ! In/output variables:
    class(type_fields_registry),                  intent(inout) :: flds_reg
    logical, dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                intent(inout) :: w
    class(*), target,                             intent(in   ) :: field_grid
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    character(len=*),                             intent(in   ) :: name
    character(len=*),                             intent(in   ) :: long_name
    character(len=*),                             intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_logical_3D'
    integer                         :: lb1, ub1, n1, lb2, ub2, n2
    class(atype_field), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call determine_bounds_2D( field_grid, field_Arakawa_grid, lb1, ub1, n1)
    call determine_bounds_3D( field_third_dimension, lb2, ub2, n2)
    call allocate_dist_shared( d, w, n1, n2)
    d( lb1: ub1, lb2: ub2) => d
    call initialise_field( field, d, w, field_grid, field_Arakawa_grid, &
      field_third_dimension, name, long_name, units)
    call flds_reg%add( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_logical_3D

  subroutine create_field_int_2D( flds_reg, d, w, field_grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    class(type_fields_registry),                 intent(inout) :: flds_reg
    integer, dimension(:), contiguous, pointer,  intent(inout) :: d
    type(MPI_WIN),                               intent(inout) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                            intent(in   ) :: name
    character(len=*),                            intent(in   ) :: long_name
    character(len=*),                            intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_int_2D'
    integer                         :: lb, ub, n
    class(atype_field), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call determine_bounds_2D( field_grid, field_Arakawa_grid, lb, ub, n)
    call allocate_dist_shared( d, w, n)
    d( lb: ub) => d
    call initialise_field( field, d, w, field_grid, field_Arakawa_grid, name, long_name, units)
    call flds_reg%add( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_int_2D

  subroutine create_field_int_3D( flds_reg, d, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units)

    ! In/output variables:
    class(type_fields_registry),                  intent(inout) :: flds_reg
    integer, dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                intent(inout) :: w
    class(*), target,                             intent(in   ) :: field_grid
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    character(len=*),                             intent(in   ) :: name
    character(len=*),                             intent(in   ) :: long_name
    character(len=*),                             intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_int_3D'
    integer                         :: lb1, ub1, n1, lb2, ub2, n2
    class(atype_field), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call determine_bounds_2D( field_grid, field_Arakawa_grid, lb1, ub1, n1)
    call determine_bounds_3D( field_third_dimension, lb2, ub2, n2)
    call allocate_dist_shared( d, w, n1, n2)
    d( lb1: ub1, lb2: ub2) => d
    call initialise_field( field, d, w, field_grid, field_Arakawa_grid, &
      field_third_dimension, name, long_name, units)
    call flds_reg%add( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_int_3D

  subroutine create_field_dp_2D( flds_reg, d, w, field_grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    class(type_fields_registry),                 intent(inout) :: flds_reg
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                               intent(inout) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                            intent(in   ) :: name
    character(len=*),                            intent(in   ) :: long_name
    character(len=*),                            intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_dp_2D'
    integer                         :: lb, ub, n
    class(atype_field), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call determine_bounds_2D( field_grid, field_Arakawa_grid, lb, ub, n)
    call allocate_dist_shared( d, w, n)
    d( lb: ub) => d
    call initialise_field( field, d, w, field_grid, field_Arakawa_grid, name, long_name, units)
    call flds_reg%add( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_dp_2D

  subroutine create_field_dp_3D( flds_reg, d, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units)

    ! In/output variables:
    class(type_fields_registry),                   intent(inout) :: flds_reg
    real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                 intent(inout) :: w
    class(*), target,                              intent(in   ) :: field_grid
    type(type_Arakawa_grid),                       intent(in   ) :: field_Arakawa_grid
    type(type_third_dimension),                    intent(in   ) :: field_third_dimension
    character(len=*),                              intent(in   ) :: name
    character(len=*),                              intent(in   ) :: long_name
    character(len=*),                              intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_dp_3D'
    integer                         :: lb1, ub1, n1, lb2, ub2, n2
    class(atype_field), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call determine_bounds_2D( field_grid, field_Arakawa_grid, lb1, ub1, n1)
    call determine_bounds_3D( field_third_dimension, lb2, ub2, n2)
    call allocate_dist_shared( d, w, n1, n2)
    d( lb1: ub1, lb2: ub2) => d
    call initialise_field( field, d, w, field_grid, field_Arakawa_grid, &
      field_third_dimension, name, long_name, units)
    call flds_reg%add( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_dp_3D

  subroutine determine_bounds_2D( field_grid, field_Arakawa_grid, lb, ub, n)

    ! In/output variables:
    class(*), target,        intent(in   ) :: field_grid
    type(type_Arakawa_grid), intent(in   ) :: field_Arakawa_grid
    integer,                 intent(  out) :: lb, ub, n

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'determine_bounds_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (g => field_grid)
    class default
      call crash('invalid field_grid type')

    class is (type_grid)

      if (field_Arakawa_grid == Arakawa_grid%a()) then
        lb = g%n1
        ub = g%n2
      else
        call crash('staggered x/y-grids not supported')
      end if

    class is (type_mesh)

      if (field_Arakawa_grid == Arakawa_grid%a()) then
        lb = g%vi1
        ub = g%vi2
      elseif (field_Arakawa_grid == Arakawa_grid%b()) then
        lb = g%ti1
        ub = g%ti2
      elseif (field_Arakawa_grid == Arakawa_grid%c()) then
        lb = g%ei1
        ub = g%ei2
      else
        call crash('invalid Arakawa grid for type_mesh')
      end if

    end select

    n = ub + 1 - lb

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine determine_bounds_2D

  subroutine determine_bounds_3D( field_third_dimension, lb, ub, n)

    ! In/output variables:
    type(type_third_dimension), intent(in   ) :: field_third_dimension
    integer,                    intent(  out) :: lb, ub, n

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'determine_bounds_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    lb = 1
    ub = field_third_dimension%n
    n  = field_third_dimension%n

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine determine_bounds_3D

end submodule fields_registry_submod_create_field