submodule (fields_basic) fields_basic_submod_reallocate

contains

  subroutine reallocate( self, mesh_new)

    ! In/output variables:
    class(atype_field),      intent(inout) :: self
    type(type_mesh), target, intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_field_reallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f => self)
    class default
      call crash('invalid field class')
    class is (type_field_logical_2D)
      call reallocate_field_logical_2D( f, mesh_new)
    class is (type_field_logical_3D)
      call reallocate_field_logical_3D( f, mesh_new)
    class is (type_field_int_2D)
      call reallocate_field_int_2D( f, mesh_new)
    class is (type_field_int_3D)
      call reallocate_field_int_3D( f, mesh_new)
    class is (type_field_dp_2D)
      call reallocate_field_dp_2D( f, mesh_new)
    class is (type_field_dp_3D)
      call reallocate_field_dp_3D( f, mesh_new)
    end select

    ! Update the grid and parallel array info of the field
    call self%set_grid( mesh_new)

    if (self%is_Arakawa_grid( Arakawa_grid%a())) then
      call self%set_pai( mesh_new%pai_V)
    elseif (self%is_Arakawa_grid( Arakawa_grid%b())) then
      call self%set_pai( mesh_new%pai_Tri)
    elseif (self%is_Arakawa_grid( Arakawa_grid%c())) then
      call self%set_pai( mesh_new%pai_E)
    else
      call crash('invalid Arakawa grid')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate

  subroutine reallocate_field_logical_2D( field, mesh_new)

    ! In/output variables:
    type(type_field_logical_2D), intent(inout) :: field
    type(type_mesh), target,     intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_field_logical_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    if (field%is_Arakawa_grid( Arakawa_grid%a())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_V%n_nih)
      field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih) => field%d_nih
      field%d_nih( mesh_new%pai_V%i1: mesh_new%pai_V%i2) = .false.
    elseif (field%is_Arakawa_grid( Arakawa_grid%b())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_Tri%n_nih)
      field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih) => field%d_nih
      field%d_nih( mesh_new%pai_Tri%i1: mesh_new%pai_Tri%i2) = .false.
    elseif (field%is_Arakawa_grid( Arakawa_grid%c())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_E%n_nih)
      field%d_nih( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih) => field%d_nih
      field%d_nih( mesh_new%pai_E%i1: mesh_new%pai_E%i2) = .false.
    else
      call crash('invalid Arakawa grid')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_field_logical_2D

  subroutine reallocate_field_logical_3D( field, mesh_new)

    ! In/output variables:
    type(type_field_logical_3D), intent(inout) :: field
    type(type_mesh), target,     intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_field_logical_3D'
    type(type_third_dimension)     :: field_third_dimension
    integer                        :: nz

    ! Add routine to call stack
    call init_routine( routine_name)

    field_third_dimension = field%third_dimension()
    nz = field_third_dimension%n

    if (field%is_Arakawa_grid( Arakawa_grid%a())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_V%n_nih, nz)
      field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:nz) => field%d_nih
      field%d_nih( mesh_new%pai_V%i1: mesh_new%pai_V%i2,:) = .false.
    elseif (field%is_Arakawa_grid( Arakawa_grid%b())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_Tri%n_nih, nz)
      field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih, 1:nz) => field%d_nih
      field%d_nih( mesh_new%pai_Tri%i1: mesh_new%pai_Tri%i2,:) = .false.
    elseif (field%is_Arakawa_grid( Arakawa_grid%c())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_E%n_nih, nz)
      field%d_nih( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih, 1:nz) => field%d_nih
      field%d_nih( mesh_new%pai_E%i1: mesh_new%pai_E%i2,:) = .false.
    else
      call crash('invalid Arakawa grid')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_field_logical_3D

  subroutine reallocate_field_int_2D( field, mesh_new)

    ! In/output variables:
    type(type_field_int_2D), intent(inout) :: field
    type(type_mesh), target, intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_field_int_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    if (field%is_Arakawa_grid( Arakawa_grid%a())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_V%n_nih)
      field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih) => field%d_nih
      field%d_nih( mesh_new%pai_V%i1: mesh_new%pai_V%i2) = -42
    elseif (field%is_Arakawa_grid( Arakawa_grid%b())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_Tri%n_nih)
      field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih) => field%d_nih
      field%d_nih( mesh_new%pai_Tri%i1: mesh_new%pai_Tri%i2) = -42
    elseif (field%is_Arakawa_grid( Arakawa_grid%c())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_E%n_nih)
      field%d_nih( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih) => field%d_nih
      field%d_nih( mesh_new%pai_E%i1: mesh_new%pai_E%i2) = -42
    else
      call crash('invalid Arakawa grid')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_field_int_2D

  subroutine reallocate_field_int_3D( field, mesh_new)

    ! In/output variables:
    type(type_field_int_3D), intent(inout) :: field
    type(type_mesh), target, intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_field_int_3D'
    type(type_third_dimension)     :: field_third_dimension
    integer                        :: nz

    ! Add routine to call stack
    call init_routine( routine_name)

    field_third_dimension = field%third_dimension()
    nz = field_third_dimension%n

    if (field%is_Arakawa_grid( Arakawa_grid%a())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_V%n_nih, nz)
      field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:nz) => field%d_nih
      field%d_nih( mesh_new%pai_V%i1: mesh_new%pai_V%i2,:) = -42
    elseif (field%is_Arakawa_grid( Arakawa_grid%b())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_Tri%n_nih, nz)
      field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih, 1:nz) => field%d_nih
      field%d_nih( mesh_new%pai_Tri%i1: mesh_new%pai_Tri%i2,:) = -42
    elseif (field%is_Arakawa_grid( Arakawa_grid%c())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_E%n_nih, nz)
      field%d_nih( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih, 1:nz) => field%d_nih
      field%d_nih( mesh_new%pai_E%i1: mesh_new%pai_E%i2,:) = -42
    else
      call crash('invalid Arakawa grid')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_field_int_3D

  subroutine reallocate_field_dp_2D( field, mesh_new)

    ! In/output variables:
    type(type_field_dp_2D), intent(inout) :: field
    type(type_mesh), target, intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_field_dp_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    if (field%is_Arakawa_grid( Arakawa_grid%a())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_V%n_nih)
      field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih) => field%d_nih
      field%d_nih( mesh_new%pai_V%i1: mesh_new%pai_V%i2) = NaN
    elseif (field%is_Arakawa_grid( Arakawa_grid%b())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_Tri%n_nih)
      field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih) => field%d_nih
      field%d_nih( mesh_new%pai_Tri%i1: mesh_new%pai_Tri%i2) = NaN
    elseif (field%is_Arakawa_grid( Arakawa_grid%c())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_E%n_nih)
      field%d_nih( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih) => field%d_nih
      field%d_nih( mesh_new%pai_E%i1: mesh_new%pai_E%i2) = NaN
    else
      call crash('invalid Arakawa grid')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_field_dp_2D

  subroutine reallocate_field_dp_3D( field, mesh_new)

    ! In/output variables:
    type(type_field_dp_3D), intent(inout) :: field
    type(type_mesh), target, intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_field_dp_3D'
    type(type_third_dimension)     :: field_third_dimension
    integer                        :: nz

    ! Add routine to call stack
    call init_routine( routine_name)

    field_third_dimension = field%third_dimension()
    nz = field_third_dimension%n

    if (field%is_Arakawa_grid( Arakawa_grid%a())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_V%n_nih, nz)
      field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:nz) => field%d_nih
      field%d_nih( mesh_new%pai_V%i1: mesh_new%pai_V%i2,:) = NaN
    elseif (field%is_Arakawa_grid( Arakawa_grid%b())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_Tri%n_nih, nz)
      field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih, 1:nz) => field%d_nih
      field%d_nih( mesh_new%pai_Tri%i1: mesh_new%pai_Tri%i2,:) = NaN
    elseif (field%is_Arakawa_grid( Arakawa_grid%c())) then
      call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_E%n_nih, nz)
      field%d_nih( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih, 1:nz) => field%d_nih
      field%d_nih( mesh_new%pai_E%i1: mesh_new%pai_E%i2,:) = NaN
    else
      call crash('invalid Arakawa grid')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_field_dp_3D

end submodule fields_basic_submod_reallocate