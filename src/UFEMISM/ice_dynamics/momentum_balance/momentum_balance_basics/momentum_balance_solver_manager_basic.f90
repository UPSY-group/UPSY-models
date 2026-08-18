module momentum_balance_solver_manager_basic

  ! A manager class that can contain multiple instances of different
  ! momentum balance solvers. This is needed for e.g. the SIA/SSA,
  ! and for the hybrid DIVA/BPA, both of which need to run code from
  ! two separate solvers.

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use momentum_balance_solver_manager_data, only: atype_momentum_balance_solver_manager_data
  use mesh_types, only: type_mesh
  use parameters, only: NaN
  use reallocate_mod, only: reallocate_bounds
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use model_configuration, only: C
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_data, only: atype_ice_velocity_model_data

  implicit none

  private

  public :: atype_momentum_balance_solver_manager

  type, abstract, extends(atype_momentum_balance_solver_manager_data) :: atype_momentum_balance_solver_manager

    contains

      ! Procedures for model memory management and operation
      procedure, public :: allocate   => momentum_balance_solver_manager_allocate
      procedure, public :: deallocate => momentum_balance_solver_manager_deallocate
      procedure, public :: initialise => momentum_balance_solver_manager_initialise
      procedure, public :: run        => momentum_balance_solver_manager_run
      procedure, public :: remap      => momentum_balance_solver_manager_remap

      ! Deferred procedures that must be overridden by each individual momentum balance implementation
      procedure(momentum_balance_solver_manager_allocate_ifc),   deferred :: allocate_momentum_balance_solver_manager
      procedure(momentum_balance_solver_manager_deallocate_ifc), deferred :: deallocate_momentum_balance_solver_manager
      procedure(momentum_balance_solver_manager_initialise_ifc), deferred :: initialise_momentum_balance_solver_manager
      procedure(momentum_balance_solver_manager_run_ifc),        deferred :: run_momentum_balance_solver_manager
      procedure(momentum_balance_solver_manager_remap_ifc),      deferred :: remap_momentum_balance_solver_manager

      procedure, public                                                 :: get_model_name
      procedure(get_momentum_balance_solver_manager_name_ifc), deferred :: get_momentum_balance_solver_manager_name

  end type atype_momentum_balance_solver_manager

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine momentum_balance_solver_manager_allocate_ifc( self, region_name, mesh)
      import atype_momentum_balance_solver_manager, type_mesh
      class(atype_momentum_balance_solver_manager),  intent(inout) :: self
      character(len=*),                              intent(in   ) :: region_name
      type(type_mesh), target,                       intent(in   ) :: mesh
    end subroutine momentum_balance_solver_manager_allocate_ifc

    subroutine momentum_balance_solver_manager_deallocate_ifc( self)
      import atype_momentum_balance_solver_manager
      class(atype_momentum_balance_solver_manager), intent(inout) :: self
    end subroutine momentum_balance_solver_manager_deallocate_ifc

    subroutine momentum_balance_solver_manager_initialise_ifc( self)
      import atype_momentum_balance_solver_manager
      class(atype_momentum_balance_solver_manager), intent(inout) :: self
    end subroutine momentum_balance_solver_manager_initialise_ifc

    subroutine momentum_balance_solver_manager_run_ifc( self, ice, geom, vel)
      import atype_momentum_balance_solver_manager, atype_ice_model_data, atype_ice_geometry_model_data, &
        atype_ice_velocity_model_data
      class(atype_momentum_balance_solver_manager), intent(inout) :: self
      class(atype_ice_model_data),                  intent(inout) :: ice
      class(atype_ice_geometry_model_data),         intent(in   ) :: geom
      class(atype_ice_velocity_model_data),         intent(inout) :: vel
    end subroutine momentum_balance_solver_manager_run_ifc

    subroutine momentum_balance_solver_manager_remap_ifc( self, mesh_new)
      import atype_momentum_balance_solver_manager, type_mesh
      class(atype_momentum_balance_solver_manager), intent(inout) :: self
      type(type_mesh), target,                      intent(in   ) :: mesh_new
    end subroutine momentum_balance_solver_manager_remap_ifc

    function get_momentum_balance_solver_manager_name_ifc( self) result( momentum_balance_solver_manager_name)
      import atype_momentum_balance_solver_manager
      class(atype_momentum_balance_solver_manager), intent(in) :: self
      character(len=:), allocatable :: momentum_balance_solver_manager_name
    end function get_momentum_balance_solver_manager_name_ifc

  end interface

contains

  subroutine momentum_balance_solver_manager_allocate( self, region_name, mesh)

    ! In/output variables:
    class(atype_momentum_balance_solver_manager), intent(inout) :: self
    character(len=*),              intent(in   ) :: region_name
    type(type_mesh), target,       intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_manager_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all models
    call self%allocate_model( region_name, mesh)

    ! Allocate all the stuff that is common to all momentum balance solver managers


    ! Allocate stuff that is specific to each individual momentum balance solver manager
    call self%allocate_momentum_balance_solver_manager( region_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_manager_allocate

  subroutine momentum_balance_solver_manager_deallocate( self)

    ! In/output variables:
    class(atype_momentum_balance_solver_manager), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_manager_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate stuff that is common to all models
    call self%deallocate_model()

    ! Deallocate all the stuff that is common to all momentum balance solver managers


    ! Deallocate stuff that is specific to each individual momentum balance solver manager
    call self%deallocate_momentum_balance_solver_manager()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_manager_deallocate

  subroutine momentum_balance_solver_manager_initialise( self)

    ! In/output variables:
    class(atype_momentum_balance_solver_manager), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_manager_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise stuff that is common to all models
    call self%initialise_model()

    ! Initialise all the stuff that is common to all momentum balance solver managers


    ! Initialise stuff that is specific to each individual momentum balance solver manager
    call self%initialise_momentum_balance_solver_manager()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_manager_initialise

  subroutine momentum_balance_solver_manager_run( self, ice, geom, vel)

    ! In/output variables:
    class(atype_momentum_balance_solver_manager),        intent(inout) :: self
    class(atype_ice_model_data),          intent(inout) :: ice
    class(atype_ice_geometry_model_data), intent(in   ) :: geom
    class(atype_ice_velocity_model_data), intent(inout) :: vel

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_manager_run'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run stuff that is common to all models
    call self%run_model()

    ! Run all the stuff that is common to all momentum balance solver managers

    ! Run stuff that is specific to each individual momentum balance solver manager
    call self%run_momentum_balance_solver_manager( ice, geom, vel)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_manager_run

  subroutine momentum_balance_solver_manager_remap( self, mesh_new)

    ! In/output variables:
    class(atype_momentum_balance_solver_manager), intent(inout) :: self
    type(type_mesh),                              intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_manager_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap stuff that is common to all models
    call self%remap_model( mesh_new)

    ! Remap all the stuff that is common to all momentum balance solver managers


    ! Remap stuff that is specific to each individual momentum balance solver manager
    call self%remap_momentum_balance_solver_manager( mesh_new)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_manager_remap

  function get_model_name( self) result( model_name)
    class(atype_momentum_balance_solver_manager), intent(in) :: self
    character(len=:), allocatable      :: model_name
    model_name = 'momentum_balance_solver_manager_' // self%get_momentum_balance_solver_manager_name()
  end function get_model_name

end module momentum_balance_solver_manager_basic