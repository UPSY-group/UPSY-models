submodule(climate_model_basic) climate_model_initialise

contains

  function ct_initialise( refgeo_init, refgeo_PD) result( context)
    !< Create a contect object for climate_model%initialise
    type(type_reference_geometry), target, intent(in) :: refgeo_init, refgeo_PD
    type(type_climate_model_context_initialise)  :: context
    context%refgeo_init => refgeo_init
    context%refgeo_PD   => refgeo_PD
  end function ct_initialise

  subroutine initialise_model_abs( self, context)

    ! In/output variables:
    class(atype_climate_model),                    intent(inout) :: self
    class(atype_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=*), parameter :: routine_name = 'atype_climate_model_initialise_model_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Set time of next calculation to start time
    self%t_next = C%start_time_of_run

    ! Downcast context class
    select type (ct => context)
    class default
      call crash('invalid context class; should be atype_climate_model_context_initialise')
    class is (type_climate_model_context_initialise)
      call initialise_model( self, ct)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_model_abs

  subroutine initialise_model( self, context)

    ! In/output variables:
    class(atype_climate_model),                          intent(inout) :: self
    type(type_climate_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=*), parameter :: routine_name = 'atype_climate_model_initialise_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_climate_model
    call initialise_model_common( self)

    ! Part specific to the model classes inheriting from atype_climate_model
    call self%initialise_climate_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_model

  subroutine initialise_model_common( self)

    ! In/output variables:
    class(atype_climate_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_model_common

end submodule climate_model_initialise