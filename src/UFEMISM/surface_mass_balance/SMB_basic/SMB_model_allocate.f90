submodule(SMB_model_basic) SMB_model_allocate

contains

  subroutine deallocate_model( self)

    ! In/output variables:
    class(atype_SMB_model), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_SMB_model_deallocate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_SMB_model
    call deallocate_model_common( self)

    ! Part specific to the model classes inheriting from atype_SMB_model
    call self%deallocate_SMB_model

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_model

  subroutine deallocate_model_common( self)

    ! In/output variables:
    class(atype_SMB_model), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'deallocate_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    nullify( self%SMB)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_model_common

end submodule SMB_model_allocate
