submodule( models_demo) models_demo_submod_remap

contains

  function demo_model_context_remap( self, mesh_new) result( context)
    !< Return an instance of the model remap context type
    class(type_demo_model),  intent(in   ) :: self
    type(type_mesh), target, intent(in   ) :: mesh_new
    type(type_demo_model_context_remap)   :: context
    context%mesh_new => mesh_new
  end function demo_model_context_remap

  subroutine remap_demo_model_abs( self, context)

    ! In/output variables:
    class(type_demo_model),     intent(inout) :: self
    class(atype_model_context), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_demo_model_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (ct => context)
    class default
      call crash('context should be of type type_demo_model_context_remap')
    class is (type_demo_model_context_remap)
      call remap_demo_model( self, self%mesh, ct%mesh_new)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_demo_model_abs

  subroutine remap_demo_model( self, mesh_old, mesh_new)

    ! In/output variables:
    class(type_demo_model), intent(inout) :: self
    type(type_mesh),        intent(in   ) :: mesh_old, mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_demo_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! call self%remap_field( )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_demo_model

end submodule models_demo_submod_remap