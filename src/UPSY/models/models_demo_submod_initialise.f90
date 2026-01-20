submodule( models_demo) models_demo_submod_initialise

contains

  function demo_model_context_initialise( self, a, b) result( context)
    !< Return an instance of the model initialisation context type
    class(type_demo_model),    intent(in   ) :: self
    integer,                   intent(in   ) :: a, b
    type(type_demo_model_context_initialise) :: context
    context%a = a
    context%b = b
  end function demo_model_context_initialise

  subroutine initialise_demo_model_abs( self, context)

    ! In/output variables:
    class(type_demo_model),     intent(inout) :: self
    class(atype_model_context), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_demo_model_abs'
    type(type_mesh), pointer       :: mesh

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (ct => context)
    class default
      call crash('context should be of type type_demo_model_context_initialise')
    class is (type_demo_model_context_initialise)
      call initialise_demo_model( self, self%mesh, ct%a, ct%b)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_demo_model_abs

  subroutine initialise_demo_model( self, mesh, a, b)

    ! In/output variables:
    class(type_demo_model), intent(inout) :: self
    type(type_mesh),        intent(in   ) :: mesh
    integer,                intent(in   ) :: a, b

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_demo_model'
    integer                        :: nz,vi,ti,k,m
    real(dp)                       :: x,y,cx,cy

    ! Add routine to call stack
    call init_routine( routine_name)

    nz = size( self%u_3D,2)

    ! Initialise fields with simple analytical functions
    cx = mesh%xmax - mesh%xmin
    cy = mesh%ymax - mesh%ymin

    do vi = mesh%vi1, mesh%vi2

      x = mesh%V( vi,1)
      y = mesh%V( vi,2)

      self%H( vi) = max( 0._dp, cos( x * pi / cx) * cos( y * pi / cy) - 0.2_dp)
      self%mask_ice( vi) = self%H( vi) > 0._dp

      do m = 1, 12
        self%T2m( vi,m) = hypot( x,y) + sin( real(m,dp) * 2._dp * pi / 12._dp)
      end do

    end do

    do ti = mesh%ti1, mesh%ti2

      x = mesh%Trigc( ti,1)
      y = mesh%Trigc( ti,2)

      do k = 1, nz
        self%u_3D( ti,k) = max( 0._dp, cos( x * pi / cx) * cos( y * pi / cy) - 0.2_dp) + real( k,dp)
        self%v_3D( ti,k) =             sin( x * pi / cx) * sin( y * pi / cy)           + real( k,dp)
      end do

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_demo_model

end submodule models_demo_submod_initialise