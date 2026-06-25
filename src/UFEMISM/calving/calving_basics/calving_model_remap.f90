submodule(calving_model_basic) calving_model_remap

contains

  function ct_remap( mesh_new) result( context)
    !< Create a contect object for calving_model%remap
    type(type_mesh),                   target, intent(in) :: mesh_new
    type(type_calving_model_context_remap)                :: context
    context%mesh_new    => mesh_new
  end function ct_remap

  subroutine remap_model_abs( self, context)

    ! In/output variables:
    class(atype_calving_model),               intent(inout) :: self
    class(atype_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_calving_model_remap_model_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast context class
    select type (ct => context)
    class default
      call crash('invalid context class; should be atype_calving_model_context_remap')
    class is (type_calving_model_context_remap)
      call remap_model( self, ct)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_model_abs

  subroutine remap_model( self, context)

    ! In/output variables:
    class(atype_calving_model),                     intent(inout) :: self
    type(type_calving_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_calving_model_remap_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_calving_model
    call remap_model_common( self, context%mesh_new)

    ! Part specific to the model classes inheriting from atype_calving_model
    call self%remap_calving_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_model

  subroutine remap_model_common( self, mesh_new)

    ! In/output variables:
    class(atype_calving_model), intent(inout) :: self
    type(type_mesh),        intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%remap_field( mesh_new, 'Hi_calved', self%Hi_calved)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_model_common

end submodule calving_model_remap