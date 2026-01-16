submodule( models_basic) models_basic_submod_create_field

contains

  subroutine create_field_logical_2D( model, d_nih, w, field_grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    class(atype_model),                          intent(inout) :: model
    logical, dimension(:), contiguous, pointer,  intent(inout) :: d_nih
    type(MPI_WIN),                               intent(inout) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                            intent(in   ) :: name
    character(len=*),                            intent(in   ) :: long_name
    character(len=*),                            intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_logical_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call model%flds_reg%create_field( d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_logical_2D

  subroutine create_field_logical_3D( model, d_nih, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units)

    ! In/output variables:
    class(atype_model),                           intent(inout) :: model
    logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    type(MPI_WIN),                                intent(inout) :: w
    class(*), target,                             intent(in   ) :: field_grid
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    character(len=*),                             intent(in   ) :: name
    character(len=*),                             intent(in   ) :: long_name
    character(len=*),                             intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_logical_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call model%flds_reg%create_field( d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_logical_3D

  subroutine create_field_int_2D( model, d_nih, w, field_grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    class(atype_model),                          intent(inout) :: model
    integer, dimension(:), contiguous, pointer,  intent(inout) :: d_nih
    type(MPI_WIN),                               intent(inout) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                            intent(in   ) :: name
    character(len=*),                            intent(in   ) :: long_name
    character(len=*),                            intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_int_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call model%flds_reg%create_field( d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_int_2D

  subroutine create_field_int_3D( model, d_nih, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units)

    ! In/output variables:
    class(atype_model),                           intent(inout) :: model
    integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    type(MPI_WIN),                                intent(inout) :: w
    class(*), target,                             intent(in   ) :: field_grid
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    character(len=*),                             intent(in   ) :: name
    character(len=*),                             intent(in   ) :: long_name
    character(len=*),                             intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_int_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call model%flds_reg%create_field( d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_int_3D

  subroutine create_field_dp_2D( model, d_nih, w, field_grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    class(atype_model),                          intent(inout) :: model
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
    type(MPI_WIN),                               intent(inout) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                            intent(in   ) :: name
    character(len=*),                            intent(in   ) :: long_name
    character(len=*),                            intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_dp_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call model%flds_reg%create_field( d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_dp_2D

  subroutine create_field_dp_3D( model, d_nih, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units)

    ! In/output variables:
    class(atype_model),                            intent(inout) :: model
    real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    type(MPI_WIN),                                 intent(inout) :: w
    class(*), target,                              intent(in   ) :: field_grid
    type(type_Arakawa_grid),                       intent(in   ) :: field_Arakawa_grid
    type(type_third_dimension),                    intent(in   ) :: field_third_dimension
    character(len=*),                              intent(in   ) :: name
    character(len=*),                              intent(in   ) :: long_name
    character(len=*),                              intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_dp_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call model%flds_reg%create_field( d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_dp_3D

end submodule models_basic_submod_create_field
