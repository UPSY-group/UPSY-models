submodule(ocean_model_basic) ocean_model_run

contains

  function ct_run( time) result( context)
    !< Create a contect object for ocean_model%run
    real(dp), intent(in) :: time
    type(type_ocean_model_context_run) :: context
    context%time =  time
  end function ct_run

  subroutine run_model_abs( self, context)

    ! In/output variables:
    class(atype_ocean_model),               intent(inout) :: self
    class(atype_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_ocean_model_run_model_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast context class
    select type (ct => context)
    class default
      call crash('invalid context class; should be atype_ocean_model_context_run')
    class is (type_ocean_model_context_run)
      call run_model( self, ct)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_model_abs

  subroutine run_model( self, context)

    ! In/output variables:
    class(atype_ocean_model),                   intent(inout) :: self
    type(type_ocean_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_ocean_model_run_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Check if we need to calculate a new ocean
    if (C%do_asynchronous_ocean) then
      ! Asynchronous coupling: do not calculate a new ocean in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next ocean time step
      if (context%time == self%t_next) then
        ! Go on to calculate a new ocean
        self%t_next = context%time + C%dt_ocean
      elseif (context%time > self%t_next) then
        ! This should not be possible
        call crash('overshot the ocean time step')
      else
        ! It is not yet time to calculate a new ocean
        call finalise_routine( routine_name)
        return
      end if

    else
      ! Synchronous coupling: calculate a new ocean in every model loop
      self%t_next = context%time + C%dt_ocean
    end if

    ! Part common to all models of atype_ocean_model
    call run_model_common( self, context)

    ! Part specific to the model classes inheriting from atype_ocean_model
    call self%run_ocean_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_model

  subroutine run_model_common( self, context)

    ! In/output variables:
    class(atype_ocean_model),                   intent(inout) :: self
    type(type_ocean_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_model_common

end submodule ocean_model_run