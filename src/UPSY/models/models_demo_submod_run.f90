submodule( models_demo) models_demo_submod_run

contains

  function demo_model_context_run( self, c, d) result( context)
    !< Return an instance of the model run context type
    class(type_demo_model), intent(in   ) :: self
    integer,                intent(in   ) :: c, d
    type(type_demo_model_context_run)     :: context
    context%c = c
    context%d = d
  end function demo_model_context_run

  subroutine run_demo_model_abs( self, context)
    !< Run the demo model

    ! In/output variables:
    class(type_demo_model),     intent(inout) :: self
    class(atype_model_context), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_demo_model_abs'
    type(type_mesh), pointer       :: mesh

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (ct => context)
    class default
      call crash('context should be of type type_demo_model_context_run')
    class is (type_demo_model_context_run)
      call run_demo_model( self, self%mesh, ct%c, ct%d)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_demo_model_abs

  subroutine run_demo_model( self, mesh, c, d)

    ! In/output variables:
    class(type_demo_model), intent(inout) :: self
    type(type_mesh),        intent(in   ) :: mesh
    integer,                intent(in   ) :: c, d

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_demo_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Physics! Maths! Science!

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_demo_model

end submodule models_demo_submod_run