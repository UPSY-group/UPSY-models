submodule(models_basic) models_basic_submod_initialise

contains

  subroutine initialise( self, context)
    !< Initialise an instance of a model

    ! In/output variables:
    class(atype_model),                            intent(inout) :: self
    class(atype_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_model_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_model
    call initialise_common( self, context)

    ! Part specific to the model classes inheriting from atype_model
    call self%initialise_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise

  subroutine initialise_common( self, context)

    ! In/output variables:
    class(atype_model),                            intent(inout) :: self
    class(atype_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_common

end submodule models_basic_submod_initialise