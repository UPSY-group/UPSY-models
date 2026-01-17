submodule (fields_basic) fields_basic_submod_remap

contains

  subroutine remap( self, mesh_new)

    ! In/output variables:
    class(atype_field),      intent(inout) :: self
    type(type_mesh), target, intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_field_remap'
    type(type_mesh), pointer       :: mesh_old

    ! Add routine to call stack
    call init_routine( routine_name)

    if (self%remap_method() == 'reallocate') then
      call self%reallocate( mesh_new)
      call finalise_routine( routine_name)
      return
    end if

    select type (g => self%grid())
    class default
      call crash('remapping only defined for mesh-based fields')
    class is (type_mesh)

      select type (f => self)
      class default
        call crash('invalid field class')

      class is (type_field_dp_2D)

        if (f%Arakawa_grid() == Arakawa_grid%a()) then
          call remap_field_dp_2D_a( f, g, mesh_new)
        elseif (f%Arakawa_grid() == Arakawa_grid%b()) then
          call remap_field_dp_2D_b( f, g, mesh_new)
        else
          call crash('invalid Arakawa grid')
        end if

      class is (type_field_dp_3D)

        if (f%Arakawa_grid() == Arakawa_grid%a()) then
          call remap_field_dp_3D_a( f, g, mesh_new)
        elseif (f%Arakawa_grid() == Arakawa_grid%b()) then
          call remap_field_dp_3D_b( f, g, mesh_new)
        else
          call crash('invalid Arakawa grid')
        end if

      end select

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap

  subroutine remap_field_dp_2D_a( field, mesh_old, mesh_new)

    ! In/output variables:
    type(type_field_dp_2D),  intent(inout) :: field
    type(type_mesh), target, intent(in   ) :: mesh_old, mesh_new

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'remap_field_dp_2D_a'
    real(dp), dimension(:), allocatable :: d_old
    real(dp), dimension(:), pointer     :: d_new

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Copy data on old mesh to d_old
    allocate( d_old( mesh_old%vi1: mesh_old%vi2), source = &
        field%d_nih( mesh_old%vi1: mesh_old%vi2))

    ! Reallocate array to dimension of new mesh
    call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_V%n_nih)
    field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih) => field%d_nih

    ! Associate pointer with the process-local part of the new array
    ! (since the remapping routines have not yet been adapted to work on hybrid memory)
    d_new( mesh_new%vi1: mesh_new%vi2) => field%d_nih( mesh_new%vi1: mesh_new%vi2)

    ! Perform the actual remapping operation
    call map_from_mesh_to_mesh_2D( mesh_old, mesh_new, trim( C%output_dir), &
      d_old, d_new, method = field%remap_method())

    ! Clean up after yourself
    deallocate( d_old)
    nullify( d_new)

    ! Update the grid of the field
    call field%set_grid( mesh_new)
    call field%set_pai( mesh_new%pai_V)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_field_dp_2D_a

  subroutine remap_field_dp_2D_b( field, mesh_old, mesh_new)

    ! In/output variables:
    type(type_field_dp_2D),  intent(inout) :: field
    type(type_mesh), target, intent(in   ) :: mesh_old, mesh_new

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'remap_field_dp_2D_b'
    real(dp), dimension(:), allocatable :: d_old
    real(dp), dimension(:), pointer     :: d_new

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Copy data on old mesh to d_old
    allocate( d_old( mesh_old%ti1: mesh_old%ti2), source = &
        field%d_nih( mesh_old%ti1: mesh_old%ti2))

    ! Reallocate array to dimension of new mesh
    call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_Tri%n_nih)
    field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih) => field%d_nih

    ! Associate pointer with the process-local part of the new array
    ! (since the remapping routines have not yet been adapted to work on hybrid memory)
    d_new( mesh_new%ti1: mesh_new%ti2) => field%d_nih( mesh_new%ti1: mesh_new%ti2)

    ! Perform the actual remapping operation
    call map_from_mesh_tri_to_mesh_tri_2D( mesh_old, mesh_new, trim( C%output_dir), &
      d_old, d_new, method = field%remap_method())

    ! Clean up after yourself
    deallocate( d_old)
    nullify( d_new)

    ! Update the grid of the field
    call field%set_grid( mesh_new)
    call field%set_pai( mesh_new%pai_Tri)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_field_dp_2D_b

  subroutine remap_field_dp_3D_a( field, mesh_old, mesh_new)

    ! In/output variables:
    type(type_field_dp_3D),  intent(inout) :: field
    type(type_mesh), target, intent(in   ) :: mesh_old, mesh_new

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'remap_field_dp_3D_a'
    type(type_third_dimension)                    :: field_third_dimension
    integer                                       :: nz, k
    real(dp), dimension(:,:), allocatable, target :: d_old
    real(dp), dimension(:  ), pointer             :: d_old_k, d_new_k

    ! Add routine to call stack
    call init_routine( routine_name)

    field_third_dimension = field%third_dimension()
    nz = field_third_dimension%n

    ! Copy data on old mesh to d_old
    allocate( d_old( mesh_old%vi1: mesh_old%vi2, 1:nz), source = &
        field%d_nih( mesh_old%vi1: mesh_old%vi2, 1:nz))

    ! Reallocate array to dimension of new mesh
    call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_V%n_nih, nz)
    field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:nz) => field%d_nih

    do k = 1, nz

      d_old_k => d_old( :,k)

      ! Associate pointer with the process-local part of the new array
      ! (since the remapping routines have not yet been adapted to work on hybrid memory)
      d_new_k => field%d_nih( mesh_new%vi1: mesh_new%vi2, k)

      ! Perform the actual remapping operation
      call map_from_mesh_to_mesh_2D( mesh_old, mesh_new, trim( C%output_dir), &
        d_old_k, d_new_k, method = field%remap_method())

    end do

    ! Clean up after yourself
    deallocate( d_old)
    nullify( d_old_k)
    nullify( d_new_k)

    ! Update the grid of the field
    call field%set_grid( mesh_new)
    call field%set_pai( mesh_new%pai_V)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_field_dp_3D_a

  subroutine remap_field_dp_3D_b( field, mesh_old, mesh_new)

    ! In/output variables:
    type(type_field_dp_3D),  intent(inout) :: field
    type(type_mesh), target, intent(in   ) :: mesh_old, mesh_new

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'remap_field_dp_3D_b'
    type(type_third_dimension)                    :: field_third_dimension
    integer                                       :: nz, k
    real(dp), dimension(:,:), allocatable, target :: d_old
    real(dp), dimension(:  ), pointer             :: d_old_k, d_new_k

    ! Add routine to call stack
    call init_routine( routine_name)

    field_third_dimension = field%third_dimension()
    nz = field_third_dimension%n

    ! Copy data on old mesh to d_old
    allocate( d_old( mesh_old%ti1: mesh_old%ti2, 1:nz), source = &
        field%d_nih( mesh_old%ti1: mesh_old%ti2, 1:nz))

    ! Reallocate array to dimension of new mesh
    call reallocate_dist_shared( field%d_nih, field%w, mesh_new%pai_Tri%n_nih, nz)
    field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih, 1:nz) => field%d_nih

    do k = 1, nz

      d_old_k => d_old( :,k)

      ! Associate pointer with the process-local part of the new array
      ! (since the remapping routines have not yet been adapted to work on hybrid memory)
      d_new_k => field%d_nih( mesh_new%ti1: mesh_new%ti2, k)

      ! Perform the actual remapping operation
      call map_from_mesh_tri_to_mesh_tri_2D( mesh_old, mesh_new, trim( C%output_dir), &
        d_old_k, d_new_k, method = field%remap_method())

    end do

    ! Clean up after yourself
    deallocate( d_old)
    nullify( d_old_k)
    nullify( d_new_k)

    ! Update the grid of the field
    call field%set_grid( mesh_new)
    call field%set_pai( mesh_new%pai_Tri)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_field_dp_3D_b

end submodule fields_basic_submod_remap