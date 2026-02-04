submodule(demo_model_basic) demo_model_allocate

contains

  function ct_allocate( name, region_name, mesh, nz) result( context)
    !< Create a contect object for demo_model%allocate
    character(len=*),           intent(in) :: name
    character(len=*),           intent(in) :: region_name
    type(type_mesh), target,    intent(in) :: mesh
    integer,                    intent(in) :: nz
    type(type_demo_model_context_allocate) :: context
    context%name        =  name
    context%region_name =  region_name
    context%mesh        => mesh
    context%nz          =  nz
  end function ct_allocate

  subroutine allocate_model_abs( self, context)

    ! In/output variables:
    class(atype_demo_model),                     intent(inout) :: self
    class(atype_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_demo_model_allocate_model_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast context class
    select type (ct => context)
    class default
      call crash('invalid context class; should be atype_demo_model_context_allocate')
    class is (type_demo_model_context_allocate)
      call allocate_model( self, ct)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_model_abs

  subroutine allocate_model( self, context)

    ! In/output variables:
    class(atype_demo_model),                        intent(inout) :: self
    type(type_demo_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_demo_model_allocate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_demo_model
    call allocate_model_common( self%mesh, self, self%s, context)

    ! Part specific to the model classes inheriting from atype_demo_model
    call self%allocate_demo_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_model

  subroutine allocate_model_common( mesh, self, demo, context)

    ! In/output variables:
    type(type_mesh),                                intent(in   ) :: mesh
    class(atype_demo_model),                        intent(inout) :: self
    type(type_demo_model_state),                    intent(inout) :: demo
    type(type_demo_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%create_field( demo%H, demo%wH, &
      mesh, Arakawa_grid%a(), &
      name      = 'H', &
      long_name = 'ice thickness', &
      units     = 'm', &
      remap_method = '2nd_order_conservative')

    call self%create_field( demo%u_3D, demo%wu_3D, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( context%nz, 'regular'), &
      name      = 'u_3D', &
      long_name = 'depth-dependent horizontal ice velocity in x-direction', &
      units     = 'm yr^-1', &
      remap_method = '2nd_order_conservative')

    call self%create_field( demo%v_3D, demo%wv_3D, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( context%nz, 'regular'), &
      name      = 'v_3D', &
      long_name = 'depth-dependent horizontal ice velocity in y-direction', &
      units     = 'm yr^-1', &
      remap_method = '2nd_order_conservative')

    call self%create_field( demo%mask_ice, demo%wmask_ice, &
      mesh, Arakawa_grid%a(), &
      name      = 'mask_ice', &
      long_name = 'ice mask', &
      units     = '-', &
      remap_method = 'reallocate')

    call self%create_field( demo%T2m, demo%wT2m, &
      mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m', &
      long_name = 'Monthly 2-m air temperature', &
      units     = 'K', &
      remap_method = '2nd_order_conservative')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_model_common

  subroutine deallocate_model( self)

    ! In/output variables:
    class(atype_demo_model), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_demo_model_deallocate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_demo_model
    call deallocate_model_common( self%s)

    ! Part specific to the model classes inheriting from atype_demo_model
    call self%deallocate_demo_model

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_model

  subroutine deallocate_model_common( demo)

    ! In/output variables:
    type(type_demo_model_state), intent(inout) :: demo

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'deallocate_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    nullify( demo%H)
    nullify( demo%u_3D)
    nullify( demo%v_3D)
    nullify( demo%mask_ice)
    nullify( demo%T2m)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_model_common

end submodule demo_model_allocate