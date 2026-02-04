submodule(SMB_model_basic) SMB_model_run

contains

  function ct_run( time, ice, climate, grid_smooth) result( context)
    !< Create a contect object for SMB_model%run
    real(dp),                          intent(in) :: time
    type(type_ice_model),     pointer, intent(in) :: ice
    type(type_climate_model), pointer, intent(in) :: climate
    type(type_grid),          pointer, intent(in) :: grid_smooth
    type(type_SMB_model_context_run)              :: context
    context%time        =  time
    context%ice         => ice
    context%climate     => climate
    context%grid_smooth => grid_smooth
  end function ct_run

  subroutine run_model_abs( self, context)

    ! In/output variables:
    class(atype_SMB_model),                 intent(inout) :: self
    class(atype_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_SMB_model_run_model_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast context class
    select type (ct => context)
    class default
      call crash('invalid context class; should be atype_SMB_model_context_run')
    class is (type_SMB_model_context_run)
      call run_model( self, ct)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_model_abs

  subroutine run_model( self, context)

    ! In/output variables:
    class(atype_SMB_model),                   intent(inout) :: self
    type(type_SMB_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_SMB_model_run_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Check if we need to calculate a new SMB
    if (C%do_asynchronous_SMB) then
      ! Asynchronous coupling: do not calculate a new SMB in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next SMB time step
      if (context%time == self%t_next) then
        ! Go on to calculate a new SMB
        self%t_next = context%time + C%dt_SMB
      elseif (context%time > self%t_next) then
        ! This should not be possible
        call crash('overshot the SMB time step')
      else
        ! It is not yet time to calculate a new SMB
        call finalise_routine( routine_name)
        return
      end if

    else
      ! Synchronous coupling: calculate a new SMB in every model loop
      self%t_next = context%time + C%dt_SMB
    end if

    ! Part common to all models of atype_SMB_model
    call run_model_common( self, context)

    ! Part specific to the model classes inheriting from atype_SMB_model
    call self%run_SMB_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_model

  subroutine run_model_common( self, context)

    ! In/output variables:
    class(atype_SMB_model),                   intent(inout) :: self
    type(type_SMB_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_model_common

end submodule SMB_model_run