submodule( SMB_basic) SMB_basic_submod_remap

contains

  function SMB_model_context_remap( self, mesh_new) result( context)
    !< Return an instance of the SMB model remapping context type
    class(atype_SMB_model),  intent(in) :: self
    type(type_mesh), target, intent(in) :: mesh_new
    type(type_SMB_model_context_remap)  :: context
    context%mesh_new => mesh_new
  end function SMB_model_context_remap

  subroutine remap_SMB_model_common_abs( self, context)
    !< Remap an instance of the SMB model

    ! In/output variables:
    class(atype_SMB_model),     intent(inout) :: self
    class(atype_model_context), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_SMB_model_common_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (ct => context)
    class default
      call crash('context should be of type type_SMB_model_context_remap')
    class is (type_SMB_model_context_remap)

      ! Common part
      call remap_SMB_model_common( self, ct%mesh_new)

      ! Specific part
      call self%remap_SMB_model( ct)

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_common_abs

  subroutine remap_SMB_model_common( self, mesh_new)

    ! In/output variables:
    class(atype_SMB_model), intent(inout) :: self
    type(type_mesh),        intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_SMB_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_common

end submodule SMB_basic_submod_remap