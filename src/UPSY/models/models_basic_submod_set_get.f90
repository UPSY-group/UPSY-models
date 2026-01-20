submodule( models_basic) models_basic_submod_set_get

contains

  ! Metadata

  subroutine set_name( self, name)
    class(atype_model), intent(inout) :: self
    character(len=*),   intent(in   ) :: name
    self%name_val = name
  end subroutine set_name

  function get_name( self) result( name)
    class(atype_model), intent(in) :: self
    character(:), allocatable      :: name
    name = self%name_val
  end function get_name

  function is_name( self, name) result( res)
    class(atype_model), intent(in) :: self
    character(len=*),   intent(in) :: name
    logical                        :: res
    res = self%name_val == name
  end function is_name

end submodule models_basic_submod_set_get
