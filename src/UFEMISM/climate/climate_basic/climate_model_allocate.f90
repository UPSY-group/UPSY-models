submodule(climate_model_basic) climate_model_allocate

contains

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
