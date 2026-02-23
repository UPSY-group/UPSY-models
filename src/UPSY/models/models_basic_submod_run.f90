submodule(models_basic) models_basic_submod_run

contains

  subroutine run( self, context)
    !< Run an instance of a model

    ! In/output variables:
    class(atype_model),                     intent(inout) :: self
    class(atype_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_model_run'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_model
    call run_common( self, context)

    ! Part specific to the model classes inheriting from atype_model
    call self%run_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run

  subroutine run_common( self, context)

    ! In/output variables:
    class(atype_model),                     intent(inout) :: self
    class(atype_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_common

end submodule models_basic_submod_run