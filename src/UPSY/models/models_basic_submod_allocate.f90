submodule(models_basic) models_basic_submod_allocate

contains

  subroutine allocate( self, context)
    !< Allocate an instance of a model

    ! In/output variables:
    class(atype_model),                          intent(inout) :: self
    class(atype_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_model_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_model
    call allocate_common( self, context)

    ! Part specific to the model classes inheriting from atype_model
    call self%allocate_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate

  subroutine allocate_common( self, context)

    ! In/output variables:
    class(atype_model),                          intent(inout) :: self
    class(atype_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Set model metadata and mesh
    call self%set_name       ( context%name)
    call self%set_region_name( context%region_name)
    self%mesh => context%mesh

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_common

  subroutine deallocate( self)
    !< Deallocate an instance of a model

    ! In/output variables:
    class(atype_model), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_model_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_model
    call deallocate_common( self)

    ! Part specific to the model classes inheriting from atype_model
    call self%deallocate_model

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate

  subroutine deallocate_common( self)

    ! In/output variables:
    class(atype_model), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'deallocate_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Set model metadata and mesh
    call self%set_name       ( 'empty_model')
    call self%set_region_name( '!!!')
    nullify( self%mesh)
    call self%flds_reg%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_common

end submodule models_basic_submod_allocate