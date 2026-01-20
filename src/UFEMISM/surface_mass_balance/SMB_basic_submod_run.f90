submodule( SMB_basic) SMB_basic_submod_run

contains

  function SMB_model_context_run( self, &
    ice, climate, grid_smooth, time, region_name) result( context)
    !< Return an instance of the SMB model run context type

    ! In/output variables:
    class(atype_SMB_model),           intent(in) :: self
    type(type_ice_model),     target, intent(in) :: ice
    type(type_climate_model), target, intent(in) :: climate
    type(type_grid),          target, intent(in) :: grid_smooth
    real(dp),                         intent(in) :: time
    character(len=3),                 intent(in) :: region_name
    type(type_SMB_model_context_run)             :: context

    context%ice         => ice
    context%climate     => climate
    context%grid_smooth => grid_smooth
    context%time        =  time
    context%region_name =  region_name

  end function SMB_model_context_run

  subroutine run_SMB_model_common_abs( self, context)
    !< Run an instance of the SMB model

    ! In/output variables:
    class(atype_SMB_model),     intent(inout) :: self
    class(atype_model_context), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_common_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (ct => context)
    class default
      call crash('context should be of type type_SMB_model_context_run')
    class is (type_SMB_model_context_run)

      ! Common part
      call run_SMB_model_common( self)

      ! Specific part
      call self%run_SMB_model( ct)

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_common_abs

  subroutine run_SMB_model_common( self)

    ! In/output variables:
    class(atype_SMB_model), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_common

end submodule SMB_basic_submod_run