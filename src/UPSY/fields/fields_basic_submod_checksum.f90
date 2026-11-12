submodule (fields_basic) fields_basic_submod_checksum

contains

  ! Logical
  ! =======

  subroutine checksum_field( self)

    ! In/output variables:
    class(atype_field), intent(in) :: self

    ! Downcast field and grid
    select type (f => self)
    class default
      call crash('invalid field type')
    class is (type_field_logical_2D)
      call checksum_field_logical_2D( f)
    class is (type_field_logical_3D)
      call checksum_field_logical_3D( f)
    class is (type_field_int_2D)
      call checksum_field_int_2D( f)
    class is (type_field_int_3D)
      call checksum_field_int_3D( f)
    class is (type_field_dp_2D)
      call checksum_field_dp_2D( f)
    class is (type_field_dp_3D)
      call checksum_field_dp_3D( f)
    end select

  end subroutine checksum_field

  subroutine checksum_field_logical_2D( field)

    ! In/output variables:
    type(type_field_logical_2D), intent(in) :: field

    ! Downcast grid
    select type (g => field%grid())
    class default
      call crash('invalid field%grid type')

    class is (type_mesh)

      if (field%Arakawa_grid() == Arakawa_grid%a()) then
        call checksum( g%pai_V, field%d_nih, field%name())
      elseif (field%Arakawa_grid() == Arakawa_grid%b()) then
        call checksum( g%pai_Tri, field%d_nih, field%name())
      elseif (field%Arakawa_grid() == Arakawa_grid%c()) then
        call checksum( g%pai_E, field%d_nih, field%name())
      else
        call crash('invalid Arakawa grid')
      end if

    class is (type_xy_grid)

      if (field%Arakawa_grid() == Arakawa_grid%a()) then
        call checksum( g%pai, field%d_nih, field%name())
      else
        call crash('invalid Arakawa grid')
      end if

    end select

  end subroutine checksum_field_logical_2D

  subroutine checksum_field_logical_3D( field)

    ! In/output variables:
    type(type_field_logical_3D), intent(in) :: field

    ! Downcast grid
    select type (g => field%grid())
    class default
      call crash('invalid field%grid type')

    class is (type_mesh)

      if (field%Arakawa_grid() == Arakawa_grid%a()) then
        call checksum( g%pai_V, field%d_nih, field%name())
      elseif (field%Arakawa_grid() == Arakawa_grid%b()) then
        call checksum( g%pai_Tri, field%d_nih, field%name())
      elseif (field%Arakawa_grid() == Arakawa_grid%c()) then
        call checksum( g%pai_E, field%d_nih, field%name())
      else
        call crash('invalid Arakawa grid')
      end if

    class is (type_xy_grid)

      if (field%Arakawa_grid() == Arakawa_grid%a()) then
        call checksum( g%pai, field%d_nih, field%name())
      else
        call crash('invalid Arakawa grid')
      end if

    end select

  end subroutine checksum_field_logical_3D

  subroutine checksum_field_int_2D( field)

    ! In/output variables:
    type(type_field_int_2D), intent(in) :: field

    ! Downcast grid
    select type (g => field%grid())
    class default
      call crash('invalid field%grid type')

    class is (type_mesh)

      if (field%Arakawa_grid() == Arakawa_grid%a()) then
        call checksum( g%pai_V, field%d_nih, field%name())
      elseif (field%Arakawa_grid() == Arakawa_grid%b()) then
        call checksum( g%pai_Tri, field%d_nih, field%name())
      elseif (field%Arakawa_grid() == Arakawa_grid%c()) then
        call checksum( g%pai_E, field%d_nih, field%name())
      else
        call crash('invalid Arakawa grid')
      end if

    class is (type_xy_grid)

      if (field%Arakawa_grid() == Arakawa_grid%a()) then
        call checksum( g%pai, field%d_nih, field%name())
      else
        call crash('invalid Arakawa grid')
      end if

    end select

  end subroutine checksum_field_int_2D

  subroutine checksum_field_int_3D( field)

    ! In/output variables:
    type(type_field_int_3D), intent(in) :: field

    ! Downcast grid
    select type (g => field%grid())
    class default
      call crash('invalid field%grid type')

    class is (type_mesh)

      if (field%Arakawa_grid() == Arakawa_grid%a()) then
        call checksum( g%pai_V, field%d_nih, field%name())
      elseif (field%Arakawa_grid() == Arakawa_grid%b()) then
        call checksum( g%pai_Tri, field%d_nih, field%name())
      elseif (field%Arakawa_grid() == Arakawa_grid%c()) then
        call checksum( g%pai_E, field%d_nih, field%name())
      else
        call crash('invalid Arakawa grid')
      end if

    class is (type_xy_grid)

      if (field%Arakawa_grid() == Arakawa_grid%a()) then
        call checksum( g%pai, field%d_nih, field%name())
      else
        call crash('invalid Arakawa grid')
      end if

    end select

  end subroutine checksum_field_int_3D

  subroutine checksum_field_dp_2D( field)

    ! In/output variables:
    type(type_field_dp_2D), intent(in) :: field

    ! Downcast grid
    select type (g => field%grid())
    class default
      call crash('invalid field%grid type')

    class is (type_mesh)

      if (field%Arakawa_grid() == Arakawa_grid%a()) then
        call checksum( g%pai_V, field%d_nih, field%name())
      elseif (field%Arakawa_grid() == Arakawa_grid%b()) then
        call checksum( g%pai_Tri, field%d_nih, field%name())
      elseif (field%Arakawa_grid() == Arakawa_grid%c()) then
        call checksum( g%pai_E, field%d_nih, field%name())
      else
        call crash('invalid Arakawa grid')
      end if

    class is (type_xy_grid)

      if (field%Arakawa_grid() == Arakawa_grid%a()) then
        call checksum( g%pai, field%d_nih, field%name())
      else
        call crash('invalid Arakawa grid')
      end if

    end select

  end subroutine checksum_field_dp_2D

  subroutine checksum_field_dp_3D( field)

    ! In/output variables:
    type(type_field_dp_3D), intent(in) :: field

    ! Downcast grid
    select type (g => field%grid())
    class default
      call crash('invalid field%grid type')

    class is (type_mesh)

      if (field%Arakawa_grid() == Arakawa_grid%a()) then
        call checksum( g%pai_V, field%d_nih, field%name())
      elseif (field%Arakawa_grid() == Arakawa_grid%b()) then
        call checksum( g%pai_Tri, field%d_nih, field%name())
      elseif (field%Arakawa_grid() == Arakawa_grid%c()) then
        call checksum( g%pai_E, field%d_nih, field%name())
      else
        call crash('invalid Arakawa grid')
      end if

    class is (type_xy_grid)

      if (field%Arakawa_grid() == Arakawa_grid%a()) then
        call checksum( g%pai, field%d_nih, field%name())
      else
        call crash('invalid Arakawa grid')
      end if

    end select

  end subroutine checksum_field_dp_3D

end submodule fields_basic_submod_checksum