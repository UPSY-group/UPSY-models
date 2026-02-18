submodule( models_basic) models_basic_submod_create_field

contains

  subroutine create_field_logical_2D( self, d_nih, w, field_grid, &
    field_Arakawa_grid, name, long_name, units, remap_method)

    ! In/output variables:
    class(atype_model),                          intent(inout) :: self
    logical, dimension(:), contiguous, pointer,  intent(inout) :: d_nih
    type(MPI_WIN),                               intent(inout) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                  optional, intent(in   ) :: name
    character(len=*),                  optional, intent(in   ) :: long_name
    character(len=*),                  optional, intent(in   ) :: units
    character(len=*),                  optional, intent(in   ) :: remap_method

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'atype_model_create_field_logical_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_nonoptional_optionals( name, long_name, units)

    call self%flds_reg%create_field( d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units, remap_method)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_logical_2D

  subroutine create_field_logical_3D( self, d_nih, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)

    ! In/output variables:
    class(atype_model),                           intent(inout) :: self
    logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    type(MPI_WIN),                                intent(inout) :: w
    class(*), target,                             intent(in   ) :: field_grid
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    character(len=*),                   optional, intent(in   ) :: name
    character(len=*),                   optional, intent(in   ) :: long_name
    character(len=*),                   optional, intent(in   ) :: units
    character(len=*),                   optional, intent(in   ) :: remap_method

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'atype_model_create_field_logical_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_nonoptional_optionals( name, long_name, units)

    call self%flds_reg%create_field( d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_logical_3D

  subroutine create_field_int_2D( self, d_nih, w, field_grid, &
    field_Arakawa_grid, name, long_name, units, remap_method)

    ! In/output variables:
    class(atype_model),                          intent(inout) :: self
    integer, dimension(:), contiguous, pointer,  intent(inout) :: d_nih
    type(MPI_WIN),                               intent(inout) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                  optional, intent(in   ) :: name
    character(len=*),                  optional, intent(in   ) :: long_name
    character(len=*),                  optional, intent(in   ) :: units
    character(len=*),                  optional, intent(in   ) :: remap_method

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'atype_model_create_field_int_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_nonoptional_optionals( name, long_name, units)

    call self%flds_reg%create_field( d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units, remap_method)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_int_2D

  subroutine create_field_int_3D( self, d_nih, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)

    ! In/output variables:
    class(atype_model),                           intent(inout) :: self
    integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    type(MPI_WIN),                                intent(inout) :: w
    class(*), target,                             intent(in   ) :: field_grid
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    character(len=*),                   optional, intent(in   ) :: name
    character(len=*),                   optional, intent(in   ) :: long_name
    character(len=*),                   optional, intent(in   ) :: units
    character(len=*),                   optional, intent(in   ) :: remap_method

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'atype_model_create_field_int_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_nonoptional_optionals( name, long_name, units)

    call self%flds_reg%create_field( d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_int_3D

  subroutine create_field_dp_2D( self, d_nih, w, field_grid, &
    field_Arakawa_grid, name, long_name, units, remap_method)

    ! In/output variables:
    class(atype_model),                          intent(inout) :: self
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
    type(MPI_WIN),                               intent(inout) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                  optional, intent(in   ) :: name
    character(len=*),                  optional, intent(in   ) :: long_name
    character(len=*),                  optional, intent(in   ) :: units
    character(len=*),                  optional, intent(in   ) :: remap_method

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'atype_model_create_field_dp_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_nonoptional_optionals( name, long_name, units)

    call self%flds_reg%create_field( d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units, remap_method)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_dp_2D

  subroutine create_field_dp_3D( self, d_nih, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)

    ! In/output variables:
    class(atype_model),                            intent(inout) :: self
    real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    type(MPI_WIN),                                 intent(inout) :: w
    class(*), target,                              intent(in   ) :: field_grid
    type(type_Arakawa_grid),                       intent(in   ) :: field_Arakawa_grid
    type(type_third_dimension),                    intent(in   ) :: field_third_dimension
    character(len=*),                    optional, intent(in   ) :: name
    character(len=*),                    optional, intent(in   ) :: long_name
    character(len=*),                    optional, intent(in   ) :: units
    character(len=*),                    optional, intent(in   ) :: remap_method

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'atype_model_create_field_dp_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_nonoptional_optionals( name, long_name, units)

    call self%flds_reg%create_field( d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_dp_3D

  subroutine check_nonoptional_optionals( name, long_name, units)
    character(len=*),  optional, intent(in   ) :: name
    character(len=*),  optional, intent(in   ) :: long_name
    character(len=*),  optional, intent(in   ) :: units
    if (.not. present( name        )) call crash('missing input argument "name"')
    if (.not. present( long_name   )) call crash('missing input argument "long_name"')
    if (.not. present( units       )) call crash('missing input argument "units"')
  end subroutine check_nonoptional_optionals

end submodule models_basic_submod_create_field
