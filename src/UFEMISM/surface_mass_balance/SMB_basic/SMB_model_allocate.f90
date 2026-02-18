submodule(SMB_model_basic) SMB_model_allocate

contains

  function ct_allocate( name, region_name, mesh) result( context)
    !< Create a contect object for SMB_model%allocate
    character(len=*),           intent(in) :: name
    character(len=*),           intent(in) :: region_name
    type(type_mesh), target,    intent(in) :: mesh
    type(type_SMB_model_context_allocate) :: context
    context%name        =  name
    context%region_name =  region_name
    context%mesh        => mesh
  end function ct_allocate

  subroutine allocate_model_abs( self, context)

    ! In/output variables:
    class(atype_SMB_model),                      intent(inout) :: self
    class(atype_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_SMB_model_allocate_model_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast context class
    select type (ct => context)
    class default
      call crash('invalid context class; should be atype_SMB_model_context_allocate')
    class is (type_SMB_model_context_allocate)
      call allocate_model( self, ct)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_model_abs

  subroutine allocate_model( self, context)

    ! In/output variables:
    class(atype_SMB_model),                        intent(inout) :: self
    type(type_SMB_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_SMB_model_allocate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_SMB_model
    call allocate_model_common( self, context)

    ! Part specific to the model classes inheriting from atype_SMB_model
    call self%allocate_SMB_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_model

  subroutine allocate_model_common( self, context)

    ! In/output variables:
    class(atype_SMB_model),                        intent(inout) :: self
    type(type_SMB_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%create_field( self%SMB, self%wSMB, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'SMB', &
      long_name = 'surface mass balance', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_model_common

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