submodule( models_demo) models_demo_submod_create

contains

  function demo_model_context_create( self, mesh) result( context)
    !< Return an instance of the demo model creation context type
    class(type_demo_model),  intent(in   ) :: self
    type(type_mesh), target, intent(in   ) :: mesh
    type(type_demo_model_context_create)  :: context
    context%mesh => mesh
  end function demo_model_context_create

  subroutine create_demo_model_abs( self, context)
    !< Create an allocated but uninitialised instance of the demo model

    ! In/output variables:
    class(type_demo_model),     intent(inout) :: self
    class(atype_model_context), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_demo_model_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (ct => context)
    class default
      call crash('context should be of type type_demo_model_context_create')
    class is (type_demo_model_context_create)
      call create_demo_model( self, ct%mesh)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_demo_model_abs

  subroutine create_demo_model( self, mesh)

    ! In/output variables:
    class(type_demo_model), intent(inout) :: self
    type(type_mesh),        intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_demo_model'
    integer, parameter             :: nz = 10

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Set model metadata and mesh
    call self%set_name('demo_model')
    call self%set_grid( mesh)

    ! Create model fields
    ! ===================

    call self%create_field( self%H, self%wH, &
      mesh, Arakawa_grid%a(), &
      name      = 'H', &
      long_name = 'ice thickness', &
      units     = 'm', &
      remap_method = '2nd_order_conservative')

    call self%create_field( self%u_3D, self%wu_3D, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = 'u_3D', &
      long_name = 'depth-dependent horizontal ice velocity in x-direction', &
      units     = 'm yr^-1', &
      remap_method = '2nd_order_conservative')

    call self%create_field( self%v_3D, self%wv_3D, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = 'v_3D', &
      long_name = 'depth-dependent horizontal ice velocity in y-direction', &
      units     = 'm yr^-1', &
      remap_method = '2nd_order_conservative')

    call self%create_field( self%mask_ice, self%wmask_ice, &
      mesh, Arakawa_grid%a(), &
      name      = 'mask_ice', &
      long_name = 'ice mask', &
      units     = '-', &
      remap_method = 'reallocate')

    call self%create_field( self%T2m, self%wT2m, &
      mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m', &
      long_name = 'Monthly 2-m air temperature', &
      units     = 'K', &
      remap_method = '2nd_order_conservative')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_demo_model

end submodule models_demo_submod_create