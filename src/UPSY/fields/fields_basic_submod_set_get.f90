submodule (fields_basic) fields_basic_submod_set_get

contains

  ! Metadata

  subroutine set_name( field, name)
    class(atype_field), intent(inout) :: field
    character(len=*),   intent(in   ) :: name
    field%name_val = name
  end subroutine set_name

  subroutine set_long_name( field, long_name)
    class(atype_field), intent(inout) :: field
    character(len=*),   intent(in   ) :: long_name
    field%long_name_val = long_name
  end subroutine set_long_name

  subroutine set_units( field, units)
    class(atype_field), intent(inout) :: field
    character(len=*),   intent(in   ) :: units
    field%units_val = units
  end subroutine set_units

  function get_name( field) result( name)
    class(atype_field), intent(in) :: field
    character(:), allocatable      :: name
    name = field%name_val
  end function get_name

  function get_long_name( field) result( long_name)
    class(atype_field), intent(in) :: field
    character(:), allocatable      :: long_name
    long_name = field%long_name_val
  end function get_long_name

  function get_units( field) result( units)
    class(atype_field), intent(in) :: field
    character(:), allocatable      :: units
    units = field%units_val
  end function get_units

  function is_name( field, name) result( res)
    class(atype_field), intent(in) :: field
    character(len=*),   intent(in) :: name
    logical                        :: res
    res = field%name_val == name
  end function is_name

  function is_long_name( field, long_name) result( res)
    class(atype_field), intent(in) :: field
    character(len=*),   intent(in) :: long_name
    logical                        :: res
    res = field%long_name_val == long_name
  end function is_long_name

  function is_units( field, units) result( res)
    class(atype_field), intent(in) :: field
    character(len=*),   intent(in) :: units
    logical                        :: res
    res = field%units_val == units
  end function is_units

  ! Grid

  subroutine set_grid( field, grid)

    ! In/output variables:
    class(atype_field), intent(inout) :: field
    class(*), target,   intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (p => grid)
    class default
      call crash('invalid grid class')
    class is (type_grid)
      field%grid_val => grid
    class is (type_mesh)
      field%grid_val => grid
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_grid

  subroutine set_Arakawa_grid( field, field_Arakawa_grid)
    class(atype_field),      intent(inout) :: field
    type(type_Arakawa_grid), intent(in   ) :: field_Arakawa_grid
    field%Arakawa_grid_val = field_Arakawa_grid
  end subroutine set_Arakawa_grid

  function get_grid( field) result( grid)

    ! In/output variables:
    class(atype_field), intent(in) :: field
    class(*), pointer              :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'get_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (p => field%grid_val)
    class default
      call crash('invalid grid class')
    class is (type_grid)
      grid => field%grid_val
    class is (type_mesh)
      grid => field%grid_val
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function get_grid

  function get_Arakawa_grid( field) result( field_Arakawa_grid)
    class(atype_field), intent(in) :: field
    type(type_Arakawa_grid)        :: field_Arakawa_grid
    field_Arakawa_grid = field%Arakawa_grid_val
  end function get_Arakawa_grid

  function is_grid( field, grid) result( res)

    ! In/output variables:
    class(atype_field), intent(in) :: field
    class(*),           intent(in) :: grid
    logical                        :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_grid'
    character(len=1024)            :: name1, name2

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (p => field%grid_val)
    class default
      call crash('invalid field%grid class')
    class is (type_grid)
      name1 = p%name
    class is (type_mesh)
      name1 = p%name
    end select

    select type (p => grid)
    class default
      call crash('invalid grid class')
    class is (type_grid)
      name2 = p%name
    class is (type_mesh)
      name2 = p%name
    end select

    res = name1 == name2

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_grid

  function is_Arakawa_grid( field, field_Arakawa_grid) result( res)
    class(atype_field),      intent(in) :: field
    type(type_Arakawa_grid), intent(in) :: field_Arakawa_grid
    logical                             :: res
    res = field_Arakawa_grid == field%Arakawa_grid_val
  end function is_Arakawa_grid

  subroutine set_third_dimension( field, field_third_dimension)
    class(atype_field_3D),      intent(inout) :: field
    type(type_third_dimension), intent(in   ) :: field_third_dimension
    field%third_dimension_val = field_third_dimension
  end subroutine set_third_dimension

  function get_third_dimension( field) result( field_third_dimension)
    class(atype_field_3D), intent(in) :: field
    type(type_third_dimension)        :: field_third_dimension
    field_third_dimension = field%third_dimension_val
  end function get_third_dimension

  function is_third_dimension( field, field_third_dimension) result( res)

    class(atype_field),         intent(in) :: field
    type(type_third_dimension), intent(in) :: field_third_dimension
    logical                                :: res

    res = .false.

    select type (field)
    class default
      call crash('invalid field type')
    class is (atype_field_2D)
    class is (atype_field_3D)
      res = field_third_dimension == field%third_dimension_val
    end select

  end function is_third_dimension

end submodule fields_basic_submod_set_get