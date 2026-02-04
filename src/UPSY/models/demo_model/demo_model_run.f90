submodule(demo_model_basic) demo_model_run

contains

  function ct_run( H_new, dH) result( context)
    !< Create a contect object for demo_model%run
    real(dp),              intent(in) :: H_new
    real(dp),              intent(in) :: dH
    type(type_demo_model_context_run) :: context
    context%H_new =  H_new
    context%dH    =  dH
  end function ct_run

  subroutine run_model_abs( self, context)

    ! In/output variables:
    class(atype_demo_model),                intent(inout) :: self
    class(atype_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_demo_model_run_model_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast context class
    select type (ct => context)
    class default
      call crash('invalid context class; should be atype_demo_model_context_run')
    class is (type_demo_model_context_run)
      call run_model( self, ct)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_model_abs

  subroutine run_model( self, context)

    ! In/output variables:
    class(atype_demo_model),                   intent(inout) :: self
    type(type_demo_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_demo_model_run_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_demo_model
    call run_model_common( self, context)

    ! Part specific to the model classes inheriting from atype_demo_model
    call self%run_demo_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_model

  subroutine run_model_common( self, context)

    ! In/output variables:
    class(atype_demo_model),                   intent(inout) :: self
    type(type_demo_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_model_common

end submodule demo_model_run