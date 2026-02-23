submodule (fields_basic) fields_basic_submod_reallocate

contains

  ! Logical
  ! =======

  subroutine reallocate_logical_2D( self, mesh_new, d_nih)

    ! In/output variables:
    class(atype_field),                         intent(inout) :: self
    type(type_mesh),                            intent(in   ) :: mesh_new
    logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_field_reallocate_logical_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast field and grid
    select type (f => self)
    class default
      call crash('invalid field type')
    class is (type_field_logical_2D)

      select type (g => self%grid())
      class default
        call crash('invalid field%grid type')
      class is (type_mesh)

        if (f%Arakawa_grid() == Arakawa_grid%a()) then
          call reallocate_logical_2D_a( f, mesh_new, d_nih)
        elseif (f%Arakawa_grid() == Arakawa_grid%b()) then
          call reallocate_logical_2D_b( f, mesh_new, d_nih)
        elseif (f%Arakawa_grid() == Arakawa_grid%c()) then
          call reallocate_logical_2D_c( f, mesh_new, d_nih)
        else
          call crash('invalid Arakawa grid')
        end if

      end select

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_logical_2D

  subroutine reallocate_logical_2D_a( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_logical_2D),                intent(inout) :: field
    type(type_mesh),                            intent(in   ) :: mesh_new
    logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_logical_2D_a'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_V%n_nih)
    d_nih      ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih) => d_nih
    field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_V)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_logical_2D_a

  subroutine reallocate_logical_2D_b( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_logical_2D),                intent(inout) :: field
    type(type_mesh),                            intent(in   ) :: mesh_new
    logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_logical_2D_b'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_Tri%n_nih)
    d_nih      ( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih) => d_nih
    field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_Tri)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_logical_2D_b

  subroutine reallocate_logical_2D_c( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_logical_2D),                intent(inout) :: field
    type(type_mesh),                            intent(in   ) :: mesh_new
    logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_logical_2D_c'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_E%n_nih)
    d_nih      ( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih) => d_nih
    field%d_nih( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_E)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_logical_2D_c

  subroutine reallocate_logical_3D( self, mesh_new, d_nih)

    ! In/output variables:
    class(atype_field),                           intent(inout) :: self
    type(type_mesh),                              intent(in   ) :: mesh_new
    logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_field_reallocate_logical_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast field and grid
    select type (f => self)
    class default
      call crash('invalid field type')
    class is (type_field_logical_3D)

      select type (g => self%grid())
      class default
        call crash('invalid field%grid type')
      class is (type_mesh)

        if (f%Arakawa_grid() == Arakawa_grid%a()) then
          call reallocate_logical_3D_a( f, mesh_new, d_nih)
        elseif (f%Arakawa_grid() == Arakawa_grid%b()) then
          call reallocate_logical_3D_b( f, mesh_new, d_nih)
        elseif (f%Arakawa_grid() == Arakawa_grid%c()) then
          call reallocate_logical_3D_c( f, mesh_new, d_nih)
        else
          call crash('invalid Arakawa grid')
        end if

      end select

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_logical_3D

  subroutine reallocate_logical_3D_a( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_logical_3D),                  intent(inout) :: field
    type(type_mesh),                              intent(in   ) :: mesh_new
    logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_logical_3D_a'
    integer                        :: nz

    ! Add routine to call stack
    call init_routine( routine_name)

    nz = size( d_nih,2)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_V%n_nih, nz)
    d_nih      ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:nz) => d_nih
    field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:nz) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_V)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_logical_3D_a

  subroutine reallocate_logical_3D_b( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_logical_3D),                  intent(inout) :: field
    type(type_mesh),                              intent(in   ) :: mesh_new
    logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_logical_3D_b'
    integer                        :: nz

    ! Add routine to call stack
    call init_routine( routine_name)

    nz = size( d_nih,2)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_Tri%n_nih, nz)
    d_nih      ( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih, 1:nz) => d_nih
    field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih, 1:nz) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_Tri)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_logical_3D_b

  subroutine reallocate_logical_3D_c( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_logical_3D),                  intent(inout) :: field
    type(type_mesh),                              intent(in   ) :: mesh_new
    logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_logical_3D_c'
    integer                        :: nz

    ! Add routine to call stack
    call init_routine( routine_name)

    nz = size( d_nih,2)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_E%n_nih, nz)
    d_nih      ( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih, 1:nz) => d_nih
    field%d_nih( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih, 1:nz) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_E)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_logical_3D_c

  ! Integer
  ! =======

  subroutine reallocate_int_2D( self, mesh_new, d_nih)

    ! In/output variables:
    class(atype_field),                         intent(inout) :: self
    type(type_mesh),                            intent(in   ) :: mesh_new
    integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_field_reallocate_int_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast field and grid
    select type (f => self)
    class default
      call crash('invalid field type')
    class is (type_field_int_2D)

      select type (g => self%grid())
      class default
        call crash('invalid field%grid type')
      class is (type_mesh)

        if (f%Arakawa_grid() == Arakawa_grid%a()) then
          call reallocate_int_2D_a( f, mesh_new, d_nih)
        elseif (f%Arakawa_grid() == Arakawa_grid%b()) then
          call reallocate_int_2D_b( f, mesh_new, d_nih)
        elseif (f%Arakawa_grid() == Arakawa_grid%c()) then
          call reallocate_int_2D_c( f, mesh_new, d_nih)
        else
          call crash('invalid Arakawa grid')
        end if

      end select

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_int_2D

  subroutine reallocate_int_2D_a( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_int_2D),                    intent(inout) :: field
    type(type_mesh),                            intent(in   ) :: mesh_new
    integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_int_2D_a'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_V%n_nih)
    d_nih      ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih) => d_nih
    field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_V)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_int_2D_a

  subroutine reallocate_int_2D_b( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_int_2D),                   intent(inout) :: field
    type(type_mesh),                            intent(in   ) :: mesh_new
    integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_int_2D_b'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_Tri%n_nih)
    d_nih      ( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih) => d_nih
    field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_Tri)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_int_2D_b

  subroutine reallocate_int_2D_c( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_int_2D),                   intent(inout) :: field
    type(type_mesh),                            intent(in   ) :: mesh_new
    integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_int_2D_c'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_E%n_nih)
    d_nih      ( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih) => d_nih
    field%d_nih( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_E)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_int_2D_c

  subroutine reallocate_int_3D( self, mesh_new, d_nih)

    ! In/output variables:
    class(atype_field),                           intent(inout) :: self
    type(type_mesh),                              intent(in   ) :: mesh_new
    integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_field_reallocate_int_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast field and grid
    select type (f => self)
    class default
      call crash('invalid field type')
    class is (type_field_int_3D)

      select type (g => self%grid())
      class default
        call crash('invalid field%grid type')
      class is (type_mesh)

        if (f%Arakawa_grid() == Arakawa_grid%a()) then
          call reallocate_int_3D_a( f, mesh_new, d_nih)
        elseif (f%Arakawa_grid() == Arakawa_grid%b()) then
          call reallocate_int_3D_b( f, mesh_new, d_nih)
        elseif (f%Arakawa_grid() == Arakawa_grid%c()) then
          call reallocate_int_3D_c( f, mesh_new, d_nih)
        else
          call crash('invalid Arakawa grid')
        end if

      end select

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_int_3D

  subroutine reallocate_int_3D_a( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_int_3D),                     intent(inout) :: field
    type(type_mesh),                              intent(in   ) :: mesh_new
    integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_int_3D_a'
    integer                        :: nz

    ! Add routine to call stack
    call init_routine( routine_name)

    nz = size( d_nih,2)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_V%n_nih, nz)
    d_nih      ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:nz) => d_nih
    field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:nz) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_V)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_int_3D_a

  subroutine reallocate_int_3D_b( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_int_3D),                     intent(inout) :: field
    type(type_mesh),                              intent(in   ) :: mesh_new
    integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_int_3D_b'
    integer                        :: nz

    ! Add routine to call stack
    call init_routine( routine_name)

    nz = size( d_nih,2)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_Tri%n_nih, nz)
    d_nih      ( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih, 1:nz) => d_nih
    field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih, 1:nz) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_Tri)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_int_3D_b

  subroutine reallocate_int_3D_c( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_int_3D),                     intent(inout) :: field
    type(type_mesh),                              intent(in   ) :: mesh_new
    integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_int_3D_c'
    integer                        :: nz

    ! Add routine to call stack
    call init_routine( routine_name)

    nz = size( d_nih,2)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_E%n_nih, nz)
    d_nih      ( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih, 1:nz) => d_nih
    field%d_nih( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih, 1:nz) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_E)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_int_3D_c

  ! Double precision
  ! ================

  subroutine reallocate_dp_2D( self, mesh_new, d_nih)

    ! In/output variables:
    class(atype_field),                          intent(inout) :: self
    type(type_mesh),                             intent(in   ) :: mesh_new
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_field_reallocate_dp_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast field and grid
    select type (f => self)
    class default
      call crash('invalid field type')
    class is (type_field_dp_2D)

      select type (g => self%grid())
      class default
        call crash('invalid field%grid type')
      class is (type_mesh)

        if (f%Arakawa_grid() == Arakawa_grid%a()) then
          call reallocate_dp_2D_a( f, mesh_new, d_nih)
        elseif (f%Arakawa_grid() == Arakawa_grid%b()) then
          call reallocate_dp_2D_b( f, mesh_new, d_nih)
        elseif (f%Arakawa_grid() == Arakawa_grid%c()) then
          call reallocate_dp_2D_c( f, mesh_new, d_nih)
        else
          call crash('invalid Arakawa grid')
        end if

      end select

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_dp_2D

  subroutine reallocate_dp_2D_a( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_dp_2D),                      intent(inout) :: field
    type(type_mesh),                             intent(in   ) :: mesh_new
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_dp_2D_a'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_V%n_nih)
    d_nih      ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih) => d_nih
    field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_V)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_dp_2D_a

  subroutine reallocate_dp_2D_b( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_dp_2D),                      intent(inout) :: field
    type(type_mesh),                             intent(in   ) :: mesh_new
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_dp_2D_b'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_Tri%n_nih)
    d_nih      ( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih) => d_nih
    field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_Tri)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_dp_2D_b

  subroutine reallocate_dp_2D_c( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_dp_2D),                      intent(inout) :: field
    type(type_mesh),                             intent(in   ) :: mesh_new
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_dp_2D_c'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_E%n_nih)
    d_nih      ( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih) => d_nih
    field%d_nih( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_E)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_dp_2D_c

  subroutine reallocate_dp_3D( self, mesh_new, d_nih)

    ! In/output variables:
    class(atype_field),                            intent(inout) :: self
    type(type_mesh),                               intent(in   ) :: mesh_new
    real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_field_reallocate_dp_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast field and grid
    select type (f => self)
    class default
      call crash('invalid field type')
    class is (type_field_dp_3D)

      select type (g => self%grid())
      class default
        call crash('invalid field%grid type')
      class is (type_mesh)

        if (f%Arakawa_grid() == Arakawa_grid%a()) then
          call reallocate_dp_3D_a( f, mesh_new, d_nih)
        elseif (f%Arakawa_grid() == Arakawa_grid%b()) then
          call reallocate_dp_3D_b( f, mesh_new, d_nih)
        elseif (f%Arakawa_grid() == Arakawa_grid%c()) then
          call reallocate_dp_3D_c( f, mesh_new, d_nih)
        else
          call crash('invalid Arakawa grid')
        end if

      end select

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_dp_3D

  subroutine reallocate_dp_3D_a( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_dp_3D),                        intent(inout) :: field
    type(type_mesh),                               intent(in   ) :: mesh_new
    real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_dp_3D_a'
    integer                        :: nz

    ! Add routine to call stack
    call init_routine( routine_name)

    nz = size( d_nih,2)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_V%n_nih, nz)
    d_nih      ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:nz) => d_nih
    field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:nz) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_V)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_dp_3D_a

  subroutine reallocate_dp_3D_b( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_dp_3D),                        intent(inout) :: field
    type(type_mesh),                               intent(in   ) :: mesh_new
    real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_dp_3D_b'
    integer                        :: nz

    ! Add routine to call stack
    call init_routine( routine_name)

    nz = size( d_nih,2)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_Tri%n_nih, nz)
    d_nih      ( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih, 1:nz) => d_nih
    field%d_nih( mesh_new%pai_Tri%i1_nih: mesh_new%pai_Tri%i2_nih, 1:nz) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_Tri)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_dp_3D_b

  subroutine reallocate_dp_3D_c( field, mesh_new, d_nih)

    ! In/output variables:
    type(type_field_dp_3D),                        intent(inout) :: field
    type(type_mesh),                               intent(in   ) :: mesh_new
    real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_dp_3D_c'
    integer                        :: nz

    ! Add routine to call stack
    call init_routine( routine_name)

    nz = size( d_nih,2)

    ! Reallocate hybrid memory and update bounds and pointers
    call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_E%n_nih, nz)
    d_nih      ( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih, 1:nz) => d_nih
    field%d_nih( mesh_new%pai_E%i1_nih: mesh_new%pai_E%i2_nih, 1:nz) => d_nih

    ! Update field grid and parallel array info
    call field%set_grid( mesh_new)
    call field%set_pai ( mesh_new%pai_E)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_dp_3D_c

end submodule fields_basic_submod_reallocate