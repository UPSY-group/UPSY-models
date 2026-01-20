submodule( SMB_basic) SMB_basic_submod_allocate

contains

  function SMB_model_context_allocate( self, mesh) result( context)
    !< Return an instance of the SMB model allocation context type
    class(atype_SMB_model),    intent(in) :: self
    type(type_mesh), target,   intent(in) :: mesh
    type(type_SMB_model_context_allocate) :: context
    context%mesh => mesh
  end function SMB_model_context_allocate

  subroutine allocate_SMB_model_common_abs( self, context)
    !< Create an allocated but uninitialised instance of the SMB model

    ! In/output variables:
    class(atype_SMB_model),     intent(inout) :: self
    class(atype_model_context), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SMB_model_common_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (ct => context)
    class default
      call crash('context should be of type type_SMB_model_context_allocate')
    class is (type_SMB_model_context_allocate)

      ! Common part
      call allocate_SMB_model_common( self, ct%mesh)

      ! Specific part
      call self%allocate_SMB_model( ct)

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_common_abs

  subroutine allocate_SMB_model_common( self, mesh)

    ! In/output variables:
    class(atype_SMB_model), intent(inout) :: self
    type(type_mesh),        intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SMB_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Set model metadata and mesh
    call self%set_name('abstract_SMB_model')
    call self%set_grid( mesh)

    ! Create model fields
    ! ===================

    call self%create_field( self%SMB, self%wSMB, &
      mesh, Arakawa_grid%a(), &
      name      = 'SMB', &
      long_name = 'surface mass balance', &
      units     = 'm yr^-1', &
      remap_method = '2nd_order_conservative')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_common

end submodule SMB_basic_submod_allocate