submodule (fields_basic) fields_basic_submod_eq

  use tests_main, only: test_eqv, test_eq

contains

  function test_field_equality( field1, field2) result( res)

    ! In/output variables:
    class(atype_field), intent(in) :: field1, field2
    logical                        :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_field_equality'

    ! Add routine to call stack
    call init_routine( routine_name)

    res = .true. .and. &
      field1%name()      == field2%name()            .and. &
      field1%long_name() == field2%long_name()       .and. &
      field1%units()     == field2%units()           .and. &
      field1%is_grid        ( field2%grid())         .and. &
      field1%is_Arakawa_grid( field2%Arakawa_grid())

    if (res .eqv. .false.) then
      call finalise_routine( routine_name)
      return
    end if

    select type (f1 => field1)
    class default
      call crash('invalid field1 type')
    class is (type_field_logical_2D)
      res = test_field_equality_logical_2D( f1, field2)
    class is (type_field_int_2D)
      res = test_field_equality_int_2D( f1, field2)
    class is (type_field_dp_2D)
      res = test_field_equality_dp_2D( f1, field2)
    class is (type_field_logical_3D)
      res = test_field_equality_logical_3D( f1, field2)
    class is (type_field_int_3D)
      res = test_field_equality_int_3D( f1, field2)
    class is (type_field_dp_3D)
      res = test_field_equality_dp_3D( f1, field2)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function test_field_equality

  function test_field_equality_logical_2D( field1, field2) result( res)

    ! In/output variables:
    class(type_field_logical_2D), intent(in) :: field1
    class(atype_field),           intent(in) :: field2
    logical                                  :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_field_equality_logical_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f2 => field2)
    class default
      res = .false.
    class is (type_field_logical_2D)
      res = test_eqv( field1%d_nih, f2%d_nih)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function test_field_equality_logical_2D

  function test_field_equality_int_2D( field1, field2) result( res)

    ! In/output variables:
    class(type_field_int_2D), intent(in) :: field1
    class(atype_field),       intent(in) :: field2
    logical                              :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_field_equality_int_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f2 => field2)
    class default
      res = .false.
    class is (type_field_int_2D)
      res = test_eq( field1%d_nih, f2%d_nih)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function test_field_equality_int_2D

  function test_field_equality_dp_2D( field1, field2) result( res)

    ! In/output variables:
    class(type_field_dp_2D), intent(in) :: field1
    class(atype_field),      intent(in) :: field2
    logical                             :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_field_equality_dp_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f2 => field2)
    class default
      res = .false.
    class is (type_field_dp_2D)
      res = test_eq( field1%d_nih, f2%d_nih)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function test_field_equality_dp_2D

  function test_field_equality_logical_3D( field1, field2) result( res)

    ! In/output variables:
    class(type_field_logical_3D), intent(in) :: field1
    class(atype_field),           intent(in) :: field2
    logical                                  :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_field_equality_logical_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f2 => field2)
    class default
      res = .false.
    class is (type_field_logical_3D)
      res = test_eqv( field1%d_nih, f2%d_nih) .and. &
        field1%third_dimension() == f2%third_dimension()
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function test_field_equality_logical_3D

  function test_field_equality_int_3D( field1, field2) result( res)

    ! In/output variables:
    class(type_field_int_3D), intent(in) :: field1
    class(atype_field),       intent(in) :: field2
    logical                              :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_field_equality_int_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f2 => field2)
    class default
      res = .false.
    class is (type_field_int_3D)
      res = test_eq( field1%d_nih, f2%d_nih) .and. &
        field1%third_dimension() == f2%third_dimension()
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function test_field_equality_int_3D

  function test_field_equality_dp_3D( field1, field2) result( res)

    ! In/output variables:
    class(type_field_dp_3D), intent(in) :: field1
    class(atype_field),      intent(in) :: field2
    logical                              :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_field_equality_dp_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f2 => field2)
    class default
      res = .false.
    class is (type_field_dp_3D)
      res = test_eq( field1%d_nih, f2%d_nih) .and. &
        field1%third_dimension() == f2%third_dimension()
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function test_field_equality_dp_3D

end submodule fields_basic_submod_eq