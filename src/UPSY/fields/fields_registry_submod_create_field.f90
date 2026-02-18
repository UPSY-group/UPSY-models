submodule( fields_registry) fields_registry_submod_create_field

contains

  subroutine create_field_logical_2D( self, d_nih, w, field_grid, &
    field_Arakawa_grid, name, long_name, units, remap_method)

    ! In/output variables:
    class(type_fields_registry),                 intent(inout) :: self
    logical, dimension(:), contiguous, pointer,  intent(inout) :: d_nih
    type(MPI_WIN),                               intent(inout) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                  optional, intent(in   ) :: name
    character(len=*),                  optional, intent(in   ) :: long_name
    character(len=*),                  optional, intent(in   ) :: units
    character(len=*),                  optional, intent(in   ) :: remap_method

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_logical_2D'
    character(:), allocatable       :: remap_method_
    integer                         :: lb, ub, n
    class(atype_field), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    if (present( remap_method)) then
      remap_method_ = remap_method
    else
      remap_method_ = '2nd_order_conservative'
    end if

    call check_nonoptional_optionals( name, long_name, units)
    call determine_bounds_2D( field_grid, field_Arakawa_grid, lb, ub, n)
    call allocate_dist_shared( d_nih, w, n)
    d_nih( lb: ub) => d_nih
    call initialise_field( field, d_nih, w, field_grid, field_Arakawa_grid, &
      name, long_name, units, remap_method_)
    call self%add( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_logical_2D

  subroutine create_field_logical_3D( self, d_nih, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)

    ! In/output variables:
    class(type_fields_registry),                  intent(inout) :: self
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
    character(len=1024), parameter  :: routine_name = 'create_field_logical_3D'
    character(:), allocatable       :: remap_method_
    integer                         :: lb1, ub1, n1, lb2, ub2, n2
    class(atype_field), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    if (present( remap_method)) then
      remap_method_ = remap_method
    else
      remap_method_ = '2nd_order_conservative'
    end if

    call check_nonoptional_optionals( name, long_name, units)
    call determine_bounds_2D( field_grid, field_Arakawa_grid, lb1, ub1, n1)
    call determine_bounds_3D( field_third_dimension, lb2, ub2, n2)
    call allocate_dist_shared( d_nih, w, n1, n2)
    d_nih( lb1: ub1, lb2: ub2) => d_nih
    call initialise_field( field, d_nih, w, field_grid, field_Arakawa_grid, &
      field_third_dimension, name, long_name, units, remap_method_)
    call self%add( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_logical_3D

  subroutine create_field_int_2D( self, d_nih, w, field_grid, &
    field_Arakawa_grid, name, long_name, units, remap_method)

    ! In/output variables:
    class(type_fields_registry),                 intent(inout) :: self
    integer, dimension(:), contiguous, pointer,  intent(inout) :: d_nih
    type(MPI_WIN),                               intent(inout) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                  optional, intent(in   ) :: name
    character(len=*),                  optional, intent(in   ) :: long_name
    character(len=*),                  optional, intent(in   ) :: units
    character(len=*),                  optional, intent(in   ) :: remap_method

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_int_2D'
    character(:), allocatable       :: remap_method_
    integer                         :: lb, ub, n
    class(atype_field), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    if (present( remap_method)) then
      remap_method_ = remap_method
    else
      remap_method_ = '2nd_order_conservative'
    end if

    call check_nonoptional_optionals( name, long_name, units)
    call determine_bounds_2D( field_grid, field_Arakawa_grid, lb, ub, n)
    call allocate_dist_shared( d_nih, w, n)
    d_nih( lb: ub) => d_nih
    call initialise_field( field, d_nih, w, field_grid, field_Arakawa_grid, &
      name, long_name, units, remap_method_)
    call self%add( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_int_2D

  subroutine create_field_int_3D( self, d_nih, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)

    ! In/output variables:
    class(type_fields_registry),                  intent(inout) :: self
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
    character(len=1024), parameter  :: routine_name = 'create_field_int_3D'
    character(:), allocatable       :: remap_method_
    integer                         :: lb1, ub1, n1, lb2, ub2, n2
    class(atype_field), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    if (present( remap_method)) then
      remap_method_ = remap_method
    else
      remap_method_ = '2nd_order_conservative'
    end if

    call check_nonoptional_optionals( name, long_name, units)
    call determine_bounds_2D( field_grid, field_Arakawa_grid, lb1, ub1, n1)
    call determine_bounds_3D( field_third_dimension, lb2, ub2, n2)
    call allocate_dist_shared( d_nih, w, n1, n2)
    d_nih( lb1: ub1, lb2: ub2) => d_nih
    call initialise_field( field, d_nih, w, field_grid, field_Arakawa_grid, &
      field_third_dimension, name, long_name, units, remap_method_)
    call self%add( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_int_3D

  subroutine create_field_dp_2D( self, d_nih, w, field_grid, &
    field_Arakawa_grid, name, long_name, units, remap_method)

    ! In/output variables:
    class(type_fields_registry),                 intent(inout) :: self
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
    type(MPI_WIN),                               intent(inout) :: w
    class(*), target,                            intent(in   ) :: field_grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                  optional, intent(in   ) :: name
    character(len=*),                  optional, intent(in   ) :: long_name
    character(len=*),                  optional, intent(in   ) :: units
    character(len=*),                  optional, intent(in   ) :: remap_method

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_field_dp_2D'
    character(:), allocatable       :: remap_method_
    integer                         :: lb, ub, n
    class(atype_field), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    if (present( remap_method)) then
      remap_method_ = remap_method
    else
      remap_method_ = '2nd_order_conservative'
    end if

    call check_nonoptional_optionals( name, long_name, units)
    call determine_bounds_2D( field_grid, field_Arakawa_grid, lb, ub, n)
    call allocate_dist_shared( d_nih, w, n)
    d_nih( lb: ub) => d_nih
    call initialise_field( field, d_nih, w, field_grid, field_Arakawa_grid, &
      name, long_name, units, remap_method_)
    call self%add( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_dp_2D

  subroutine create_field_dp_3D( self, d_nih, w, field_grid, &
    field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)

    ! In/output variables:
    class(type_fields_registry),                   intent(inout) :: self
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
    character(len=1024), parameter  :: routine_name = 'create_field_dp_3D'
    character(:), allocatable       :: remap_method_
    integer                         :: lb1, ub1, n1, lb2, ub2, n2
    class(atype_field), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    if (present( remap_method)) then
      remap_method_ = remap_method
    else
      remap_method_ = '2nd_order_conservative'
    end if

    call check_nonoptional_optionals( name, long_name, units)
    call determine_bounds_2D( field_grid, field_Arakawa_grid, lb1, ub1, n1)
    call determine_bounds_3D( field_third_dimension, lb2, ub2, n2)
    call allocate_dist_shared( d_nih, w, n1, n2)
    d_nih( lb1: ub1, lb2: ub2) => d_nih
    call initialise_field( field, d_nih, w, field_grid, field_Arakawa_grid, &
      field_third_dimension, name, long_name, units, remap_method_)
    call self%add( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_dp_3D

  ! Utilities

  subroutine check_nonoptional_optionals( name, long_name, units)
    character(len=*),  optional, intent(in   ) :: name
    character(len=*),  optional, intent(in   ) :: long_name
    character(len=*),  optional, intent(in   ) :: units
    if (.not. present( name        )) call crash('missing input argument "name"')
    if (.not. present( long_name   )) call crash('missing input argument "long_name"')
    if (.not. present( units       )) call crash('missing input argument "units"')
  end subroutine check_nonoptional_optionals

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
        lb = g%pai%i1_nih
        ub = g%pai%i2_nih
      else
        call crash('staggered x/y-grids not supported')
      end if

    class is (type_mesh)

      if (field_Arakawa_grid == Arakawa_grid%a()) then
        lb = g%pai_V%i1_nih
        ub = g%pai_V%i2_nih
      elseif (field_Arakawa_grid == Arakawa_grid%b()) then
        lb = g%pai_Tri%i1_nih
        ub = g%pai_Tri%i2_nih
      elseif (field_Arakawa_grid == Arakawa_grid%c()) then
        lb = g%pai_E%i1_nih
        ub = g%pai_E%i2_nih
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