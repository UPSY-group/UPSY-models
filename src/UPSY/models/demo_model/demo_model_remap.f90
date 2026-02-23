submodule(demo_model_basic) demo_model_remap

contains

  function ct_remap( mesh_new) result( context)
    !< Create a contect object for demo_model%remap
    type(type_mesh), target, intent(in) :: mesh_new
    type(type_demo_model_context_remap) :: context
    context%mesh_new => mesh_new
  end function ct_remap

  subroutine remap_model_abs( self, context)

    ! In/output variables:
    class(atype_demo_model),                  intent(inout) :: self
    class(atype_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_demo_model_remap_model_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast context class
    select type (ct => context)
    class default
      call crash('invalid context class; should be atype_demo_model_context_remap')
    class is (type_demo_model_context_remap)
      call remap_model( self, ct)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_model_abs

  subroutine remap_model( self, context)

    ! In/output variables:
    class(atype_demo_model),                     intent(inout) :: self
    type(type_demo_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_demo_model_remap_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_demo_model
    call remap_model_common( self, self%s, context%mesh_new)

    ! Part specific to the model classes inheriting from atype_demo_model
    call self%remap_demo_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_model

  subroutine remap_model_common( self, demo, mesh_new)

    ! In/output variables:
    class(atype_demo_model),     intent(inout) :: self
    type(type_demo_model_state), intent(inout) :: demo
    type(type_mesh),             intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%remap_field( mesh_new, 'H'       , demo%H       )
    call self%remap_field( mesh_new, 'u_3D'    , demo%u_3D    )
    call self%remap_field( mesh_new, 'v_3D'    , demo%v_3D    )
    call self%remap_field( mesh_new, 'mask_ice', demo%mask_ice)
    call self%remap_field( mesh_new, 'T2m'     , demo%T2m     )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_model_common

end submodule demo_model_remap