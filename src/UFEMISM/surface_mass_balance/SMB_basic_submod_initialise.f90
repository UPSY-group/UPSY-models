submodule( SMB_basic) SMB_basic_submod_initialise

contains

  function SMB_model_context_initialise( self, region_name) result( context)
    !< Return an instance of the SMB model initialisation context type
    class(atype_SMB_model),      intent(in) :: self
    character(len=3),            intent(in) :: region_name
    type(type_SMB_model_context_initialise) :: context
    context%region_name = region_name
  end function SMB_model_context_initialise

  subroutine initialise_SMB_model_common_abs( self, context)
    !< Initialise an instance of the SMB model

    ! In/output variables:
    class(atype_SMB_model),     intent(inout) :: self
    class(atype_model_context), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_common_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (ct => context)
    class default
      call crash('context should be of type type_SMB_model_context_initialise')
    class is (type_SMB_model_context_initialise)

      ! Common part
      call initialise_SMB_model_common( self)

      ! Specific part
      call self%initialise_SMB_model( ct)

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_common_abs

  subroutine initialise_SMB_model_common( self)

    ! In/output variables:
    class(atype_SMB_model), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_common

end submodule SMB_basic_submod_initialise