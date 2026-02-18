submodule (fields_basic) fields_basic_submod_set_get

contains

  ! Metadata

  subroutine set_name( self, name)
    class(atype_field), intent(inout) :: self
    character(len=*),   intent(in   ) :: name
    self%name_val = name
  end subroutine set_name

  subroutine set_long_name( self, long_name)
    class(atype_field), intent(inout) :: self
    character(len=*),   intent(in   ) :: long_name
    self%long_name_val = long_name
  end subroutine set_long_name

  subroutine set_units( self, units)
    class(atype_field), intent(inout) :: self
    character(len=*),   intent(in   ) :: units
    self%units_val = units
  end subroutine set_units

  subroutine set_remap_method( self, remap_method)
    class(atype_field), intent(inout) :: self
    character(len=*),   intent(in   ) :: remap_method
    self%remap_method_val = remap_method
  end subroutine set_remap_method

  function get_name( self) result( name)
    class(atype_field), intent(in) :: self
    character(:), allocatable      :: name
    name = self%name_val
  end function get_name

  function get_long_name( self) result( long_name)
    class(atype_field), intent(in) :: self
    character(:), allocatable      :: long_name
    long_name = self%long_name_val
  end function get_long_name

  function get_units( self) result( units)
    class(atype_field), intent(in) :: self
    character(:), allocatable      :: units
    units = self%units_val
  end function get_units

  function get_remap_method( self) result( remap_method)
    class(atype_field), intent(in) :: self
    character(:), allocatable      :: remap_method
    remap_method = self%remap_method_val
  end function get_remap_method

  function is_name( self, name) result( res)
    class(atype_field), intent(in) :: self
    character(len=*),   intent(in) :: name
    logical                        :: res
    res = self%name_val == name
  end function is_name

  function is_long_name( self, long_name) result( res)
    class(atype_field), intent(in) :: self
    character(len=*),   intent(in) :: long_name
    logical                        :: res
    res = self%long_name_val == long_name
  end function is_long_name

  function is_units( self, units) result( res)
    class(atype_field), intent(in) :: self
    character(len=*),   intent(in) :: units
    logical                        :: res
    res = self%units_val == units
  end function is_units

  function is_remap_method( self, remap_method) result( res)
    class(atype_field), intent(in) :: self
    character(len=*),   intent(in) :: remap_method
    logical                        :: res
    res = self%remap_method_val == remap_method
  end function is_remap_method

  ! Grid

  subroutine set_grid( self, grid)

    ! In/output variables:
    class(atype_field), intent(inout) :: self
    class(*), target,   intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (p => grid)
    class default
      call crash('invalid grid class')
    class is (type_grid)
      self%grid_val => grid
    class is (type_mesh)
      self%grid_val => grid
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_grid

  subroutine set_Arakawa_grid( self, field_Arakawa_grid)
    class(atype_field),      intent(inout) :: self
    type(type_Arakawa_grid), intent(in   ) :: field_Arakawa_grid
    self%Arakawa_grid_val = field_Arakawa_grid
  end subroutine set_Arakawa_grid

  subroutine set_pai( self, field_pai)
    class(atype_field),      intent(inout) :: self
    type(type_par_arr_info), intent(in   ) :: field_pai
    self%pai_val = field_pai
  end subroutine set_pai

  subroutine set_third_dimension( self, field_third_dimension)
    class(atype_field_3D),      intent(inout) :: self
    type(type_third_dimension), intent(in   ) :: field_third_dimension
    self%third_dimension_val = field_third_dimension
  end subroutine set_third_dimension

  function get_grid( self) result( grid)

    ! In/output variables:
    class(atype_field), intent(in) :: self
    class(*), pointer              :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'get_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (p => self%grid_val)
    class default
      call crash('invalid grid class')
    class is (type_grid)
      grid => self%grid_val
    class is (type_mesh)
      grid => self%grid_val
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function get_grid

  function get_Arakawa_grid( self) result( field_Arakawa_grid)
    class(atype_field), intent(in) :: self
    type(type_Arakawa_grid)        :: field_Arakawa_grid
    field_Arakawa_grid = self%Arakawa_grid_val
  end function get_Arakawa_grid

  function get_pai( self) result( field_pai)
    class(atype_field), intent(in) :: self
    type(type_par_arr_info)        :: field_pai
    field_pai = self%pai_val
  end function get_pai

  function get_third_dimension( self) result( field_third_dimension)
    class(atype_field_3D), intent(in) :: self
    type(type_third_dimension)        :: field_third_dimension
    field_third_dimension = self%third_dimension_val
  end function get_third_dimension

  function is_grid( self, grid) result( res)

    ! In/output variables:
    class(atype_field), intent(in) :: self
    class(*),           intent(in) :: grid
    logical                        :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_grid'
    character(len=1024)            :: name1, name2

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (p => self%grid_val)
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

  function is_Arakawa_grid( self, field_Arakawa_grid) result( res)
    class(atype_field),      intent(in) :: self
    type(type_Arakawa_grid), intent(in) :: field_Arakawa_grid
    logical                             :: res
    res = field_Arakawa_grid == self%Arakawa_grid_val
  end function is_Arakawa_grid

  function is_pai( self, field_pai) result( res)
    class(atype_field),      intent(in) :: self
    type(type_par_arr_info), intent(in) :: field_pai
    logical                             :: res
    res = field_pai == self%pai_val
  end function is_pai

  function is_third_dimension( self, field_third_dimension) result( res)

    class(atype_field),         intent(in) :: self
    type(type_third_dimension), intent(in) :: field_third_dimension
    logical                                :: res

    res = .false.

    select type (f => self)
    class default
      call crash('invalid field type')
    class is (atype_field_2D)
    class is (atype_field_3D)
      res = field_third_dimension == f%third_dimension()
    end select

  end function is_third_dimension

end submodule fields_basic_submod_set_get