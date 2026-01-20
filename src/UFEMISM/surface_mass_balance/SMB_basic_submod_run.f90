submodule( SMB_basic) SMB_basic_submod_run

contains

  function SMB_model_context_run( self) result( context)
    !< Return an instance of the SMB model run context type
    class(atype_SMB_model), intent(in) :: self
    type(type_SMB_model_context_run)   :: context
  end function SMB_model_context_run

  subroutine run_SMB_model_common_abs( self, context)
    !< Run an instance of the SMB model

    ! In/output variables:
    class(atype_SMB_model),     intent(inout) :: self
    class(atype_model_context), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_common_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (ct => context)
    class default
      call crash('context should be of type type_SMB_model_context_run')
    class is (type_SMB_model_context_run)

      ! Common part
      call run_SMB_model_common( self)

      ! Specific part
      call self%run_SMB_model( ct)

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_common_abs

  subroutine run_SMB_model_common( self)

    ! In/output variables:
    class(atype_SMB_model), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_common

end submodule SMB_basic_submod_run