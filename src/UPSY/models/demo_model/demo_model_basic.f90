module demo_model_basic

  use precisions, only: dp
  use parameters, only: pi, NaN
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use demo_model_data, only: atype_demo_model_data
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_demo_model

  type, abstract, extends(atype_demo_model_data) :: atype_demo_model
    !< Stuff that is common to all demo models
    !<
    !< (except for the variables that we want other models to
    !< be able to access, which are already defined in atype_demo_model_data)

    contains

      ! Type-bound procedures that apply to all demo models
      procedure, public :: allocate_demo_model
      procedure, public :: deallocate_demo_model
      procedure, public :: initialise_demo_model
      procedure, public :: run_demo_model
      procedure, public :: remap_demo_model

      ! Deferred procedures that must be defined by each individual demo model
      procedure(demo_model_allocate_ifc),   deferred :: allocate
      procedure(demo_model_deallocate_ifc), deferred :: deallocate
      procedure(demo_model_initialise_ifc), deferred :: initialise
      procedure(demo_model_run_ifc),        deferred :: run
      procedure(demo_model_remap_ifc),      deferred :: remap

  end type atype_demo_model

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine demo_model_allocate_ifc( self, region_name, mesh, nz)
      import atype_demo_model, type_mesh
      class(atype_demo_model), intent(inout) :: self
      character(len=*),        intent(in   ) :: region_name
      type(type_mesh), target, intent(in   ) :: mesh
      integer,                 intent(in   ) :: nz
    end subroutine demo_model_allocate_ifc

    subroutine demo_model_deallocate_ifc( self)
      import atype_demo_model
      class(atype_demo_model), intent(inout) :: self
    end subroutine demo_model_deallocate_ifc

    subroutine demo_model_initialise_ifc( self, H0, till_friction_angle_uniform, beta_sq_uniform)
      import atype_demo_model, dp
      class(atype_demo_model), intent(inout) :: self
      real(dp),                intent(in   ) :: H0
      real(dp),                intent(in   ) :: till_friction_angle_uniform
      real(dp),                intent(in   ) :: beta_sq_uniform
    end subroutine demo_model_initialise_ifc

    subroutine demo_model_run_ifc( self, H_new, dH)
      import atype_demo_model, dp
      class(atype_demo_model), intent(inout) :: self
      real(dp),                intent(in   ) :: H_new
      real(dp),                intent(in   ) :: dH
    end subroutine demo_model_run_ifc

    subroutine demo_model_remap_ifc( self, mesh_new)
      import atype_demo_model, type_mesh
      class(atype_demo_model), intent(inout) :: self
      type(type_mesh), target, intent(in   ) :: mesh_new
    end subroutine demo_model_remap_ifc

  end interface

contains

  subroutine allocate_demo_model( self, name, region_name, mesh, nz)
    !< Allocate stuff that is common to all demo models
    !< (call this from your demo model-specific allocate routine)

    ! In/output variables:
    class(atype_demo_model), intent(inout) :: self
    character(len=*),        intent(in   ) :: name
    character(len=*),        intent(in   ) :: region_name
    type(type_mesh), target, intent(in   ) :: mesh
    integer,                 intent(in   ) :: nz

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_demo_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate stuff that is common to all models
    call self%allocate_model( name, region_name, mesh)

    ! Allocate stuff that is specific to demo models

    call self%create_field( self%H, self%wH, &
      mesh, Arakawa_grid%a(), &
      name      = 'H', &
      long_name = 'ice thickness', &
      units     = 'm', &
      remap_method = '2nd_order_conservative')

    call self%create_field( self%u_3D, self%wu_3D, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = 'u_3D', &
      long_name = 'depth-dependent horizontal ice velocity in x-direction', &
      units     = 'm yr^-1', &
      remap_method = '2nd_order_conservative')

    call self%create_field( self%v_3D, self%wv_3D, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = 'v_3D', &
      long_name = 'depth-dependent horizontal ice velocity in y-direction', &
      units     = 'm yr^-1', &
      remap_method = '2nd_order_conservative')

    call self%create_field( self%mask_ice, self%wmask_ice, &
      mesh, Arakawa_grid%a(), &
      name      = 'mask_ice', &
      long_name = 'ice mask', &
      units     = '-', &
      remap_method = 'reallocate')

    call self%create_field( self%T2m, self%wT2m, &
      mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m', &
      long_name = 'Monthly 2-m air temperature', &
      units     = 'K', &
      remap_method = '2nd_order_conservative')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_demo_model

  subroutine deallocate_demo_model( self)
    !< Deallocate stuff that is common to all demo models
    !< (call this from your demo model-specific deallocate routine)

    ! In/output variables:
    class(atype_demo_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_demo_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate stuff that is common to all models
    call self%deallocate_model()

    ! Deallocate stuff that is specific to demo models

    nullify( self%H)
    nullify( self%u_3D)
    nullify( self%v_3D)
    nullify( self%mask_ice)
    nullify( self%T2m)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_demo_model

  subroutine initialise_demo_model( self, H0)
    !< Initialise stuff that is common to all demo models
    !< (call this from your demo model-specific initialise routine)

    ! In/output variables:
    class(atype_demo_model), intent(inout) :: self
    real(dp),                intent(in   ) :: H0

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_demo_model'
    integer                     :: vi,ti,k,m
    real(dp)                    :: x,y,cx,cy

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise stuff that is common to all models
    call self%initialise_model()

    ! Initialise stuff that is specific to demo models

    cx = self%mesh%xmax - self%mesh%xmin
    cy = self%mesh%ymax - self%mesh%ymin

    do vi = self%mesh%vi1, self%mesh%vi2

      x = self%mesh%V( vi,1)
      y = self%mesh%V( vi,2)

      self%H( vi) = max( H0, cos( x * pi / cx) * cos( y * pi / cy) - 0.2_dp)
      self%mask_ice( vi) = self%H( vi) > 0._dp

      do m = 1, 12
        self%T2m( vi,m) = hypot( x,y) + sin( real(m,dp) * 2._dp * pi / 12._dp)
      end do

    end do

    do ti = self%mesh%ti1, self%mesh%ti2

      x = self%mesh%Trigc( ti,1)
      y = self%mesh%Trigc( ti,2)

      do k = 1, size( self%u_3D,2)
        self%u_3D( ti,k) = max( 0._dp, cos( x * pi / cx) * cos( y * pi / cy) - 0.2_dp) + real( k,dp)
        self%v_3D( ti,k) =             sin( x * pi / cx) * sin( y * pi / cy)           + real( k,dp)
      end do

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_demo_model

  subroutine run_demo_model( self)
    !< Run stuff that is common to all demo models
    !< (call this from your demo model-specific run routine)

    ! In/output variables:
    class(atype_demo_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_demo_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run stuff that is common to all models
    call self%run_model()

    ! Run stuff that is specific to demo models

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_demo_model

  subroutine remap_demo_model( self, mesh_new)
    !< Remap stuff that is common to all demo models
    !< (call this from your demo model-specific remap routine)

    ! In/output variables:
    class(atype_demo_model), intent(inout) :: self
    type(type_mesh), target, intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_demo_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap stuff that is common to all models
    call self%remap_model( mesh_new)

    ! Remap stuff that is specific to demo models

    call self%remap_field( mesh_new, 'H'       , self%H       )
    call self%remap_field( mesh_new, 'u_3D'    , self%u_3D    )
    call self%remap_field( mesh_new, 'v_3D'    , self%v_3D    )
    call self%remap_field( mesh_new, 'mask_ice', self%mask_ice)
    call self%remap_field( mesh_new, 'T2m'     , self%T2m     )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_demo_model

end module demo_model_basic
