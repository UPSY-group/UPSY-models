submodule(models_basic) models_basic_submod_remap

contains

  subroutine remap( self, context)
    !< Remap an instance of a model

    ! In/output variables:
    class(atype_model),                       intent(inout) :: self
    class(atype_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_model_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_model
    call remap_common( self, context)

    ! Part specific to the model classes inheriting from atype_model
    call self%remap_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap

  subroutine remap_common( self, context)

    ! In/output variables:
    class(atype_model),                       intent(inout) :: self
    class(atype_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    self%mesh => context%mesh_new

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_common

end submodule models_basic_submod_remap