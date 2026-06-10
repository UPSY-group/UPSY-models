submodule(climate_model_basic) climate_model_allocate

contains

  function ct_allocate( name, region_name, mesh) result( context)
    !< Create a contect object for climate_model%allocate
    character(len=*),              intent(in) :: name
    character(len=*),              intent(in) :: region_name
    type(type_mesh), target,       intent(in) :: mesh
    type(type_climate_model_context_allocate) :: context
    context%name        =  name
    context%region_name =  region_name
    context%mesh        => mesh
  end function ct_allocate

  subroutine allocate_model_abs( self, context)

    ! In/output variables:
    class(atype_climate_model),                  intent(inout) :: self
    class(atype_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=*), parameter :: routine_name = 'atype_climate_model_allocate_model_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast context class
    select type (ct => context)
    class default
      call crash('invalid context class; should be atype_climate_model_context_allocate')
    class is (type_climate_model_context_allocate)
      call allocate_model( self, ct)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_model_abs

  subroutine allocate_model( self, context)

    ! In/output variables:
    class(atype_climate_model),                        intent(inout) :: self
    type(type_climate_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=*), parameter :: routine_name = 'atype_climate_model_allocate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_climate_model
    call allocate_model_common( self, context)

    ! Part specific to the model classes inheriting from atype_climate_model
    call self%allocate_climate_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_model

  subroutine allocate_model_common( self, context)

    ! In/output variables:
    class(atype_climate_model),                        intent(inout) :: self
    type(type_climate_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Set generic names
    call self%set_name('climate')
    call self%set_region_name(context%region_name)

    ! Allocate generic fields
    call self%create_field( self%T2m, self%wT2m, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m', &
      long_name = 'Monthly mean 2-m air temperature', &
      units     = 'K', &
      remap_method = 'reallocate')

    call self%create_field( self%Precip, self%wPrecip, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Precip', &
      long_name = 'Monthly total precipitation', &
      units     = 'm.w.e.', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_model_common

  subroutine deallocate_model( self)

    ! In/output variables:
    class(atype_climate_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'atype_climate_model_deallocate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_climate_model
    call deallocate_model_common( self)

    ! Part specific to the model classes inheriting from atype_climate_model
    call self%deallocate_climate_model

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_model

  subroutine deallocate_model_common( self)

    ! In/output variables:
    class(atype_climate_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    nullify( self%T2m)
    nullify( self%Precip)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_model_common

end submodule climate_model_allocate
