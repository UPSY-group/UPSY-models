submodule(demo_model_basic) demo_model_initialise

contains

  function ct_initialise( H0, till_friction_angle_uniform, beta_sq_uniform) result( context)
    !< Create a contect object for demo_model%initialise
    real(dp),                     intent(in) :: H0
    real(dp),                     intent(in) :: till_friction_angle_uniform
    real(dp),                     intent(in) :: beta_sq_uniform
    type(type_demo_model_context_initialise) :: context
    context%H0                          =  H0
    context%till_friction_angle_uniform = till_friction_angle_uniform
    context%beta_sq_uniform             = beta_sq_uniform
  end function ct_initialise

  subroutine initialise_model_abs( self, context)

    ! In/output variables:
    class(atype_demo_model),                       intent(inout) :: self
    class(atype_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_demo_model_initialise_model_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast context class
    select type (ct => context)
    class default
      call crash('invalid context class; should be atype_demo_model_context_initialise')
    class is (type_demo_model_context_initialise)
      call initialise_model( self, ct)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_model_abs

  subroutine initialise_model( self, context)

    ! In/output variables:
    class(atype_demo_model),                          intent(inout) :: self
    type(type_demo_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_demo_model_initialise_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Part common to all models of atype_demo_model
    call initialise_model_common( self%mesh, self%s, context%H0)

    ! Part specific to the model classes inheriting from atype_demo_model
    call self%initialise_demo_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_model

  subroutine initialise_model_common( mesh, demo, H0)

    ! In/output variables:
    type(type_mesh),              intent(in   ) :: mesh
    class(type_demo_model_state), intent(inout) :: demo
    real(dp),                     intent(in   ) :: H0

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_model_common'
    integer                        :: vi,ti,k,m
    real(dp)                       :: x,y,cx,cy

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise fields with simple analytical functions

    cx = mesh%xmax - mesh%xmin
    cy = mesh%ymax - mesh%ymin

    do vi = mesh%vi1, mesh%vi2

      x = mesh%V( vi,1)
      y = mesh%V( vi,2)

      demo%H( vi) = max( H0, cos( x * pi / cx) * cos( y * pi / cy) - 0.2_dp)
      demo%mask_ice( vi) = demo%H( vi) > 0._dp

      do m = 1, 12
        demo%T2m( vi,m) = hypot( x,y) + sin( real(m,dp) * 2._dp * pi / 12._dp)
      end do

    end do

    do ti = mesh%ti1, mesh%ti2

      x = mesh%Trigc( ti,1)
      y = mesh%Trigc( ti,2)

      do k = 1, size( demo%u_3D,2)
        demo%u_3D( ti,k) = max( 0._dp, cos( x * pi / cx) * cos( y * pi / cy) - 0.2_dp) + real( k,dp)
        demo%v_3D( ti,k) =             sin( x * pi / cx) * sin( y * pi / cy)           + real( k,dp)
      end do

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_model_common

end submodule demo_model_initialise