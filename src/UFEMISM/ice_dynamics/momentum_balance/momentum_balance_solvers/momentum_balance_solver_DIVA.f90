module momentum_balance_solver_DIVA

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use momentum_balance_solver_plain_DIVA, only: type_momentum_balance_solver_plain_DIVA
  use bed_roughness_model_types, only: type_bed_roughness_model

  implicit none

  private

  public :: type_momentum_balance_solver_DIVA

  type, extends(atype_momentum_balance_solver) :: type_momentum_balance_solver_DIVA

      type(type_momentum_balance_solver_plain_DIVA) :: solver

    contains

      ! Procedures for model memory management and operation
      procedure, public :: allocate_momentum_balance_solver   => momentum_balance_solver_DIVA_allocate
      procedure, public :: deallocate_momentum_balance_solver => momentum_balance_solver_DIVA_deallocate
      procedure, public :: initialise_momentum_balance_solver => momentum_balance_solver_DIVA_initialise
      procedure, public :: run_momentum_balance_solver        => momentum_balance_solver_DIVA_run
      procedure, public :: remap_momentum_balance_solver      => momentum_balance_solver_DIVA_remap

      procedure, public :: get_momentum_balance_solver_name

  end type type_momentum_balance_solver_DIVA

contains

  subroutine momentum_balance_solver_DIVA_allocate( self, region_name, mesh)

    ! In/output variables:
    class(type_momentum_balance_solver_DIVA), intent(inout) :: self
    character(len=*),                        intent(in   ) :: region_name
    type(type_mesh), target,                 intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_DIVA_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is specific to the DIVA momentum balance solver
    call self%solver%allocate( region_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_DIVA_allocate

  subroutine momentum_balance_solver_DIVA_deallocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_DIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_DIVA_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to the DIVA momentum balance solver
    call self%solver%deallocate()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_DIVA_deallocate

  subroutine momentum_balance_solver_DIVA_initialise( self)

    ! In/output variables:
    class(type_momentum_balance_solver_DIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_DIVA_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise all the stuff that is specific to the DIVA momentum balance solver
    call self%solver%initialise()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_DIVA_initialise

  subroutine momentum_balance_solver_DIVA_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)

    ! In/output variables:
    class(type_momentum_balance_solver_DIVA), intent(inout) :: self
    class(atype_ice_model_data),             intent(inout) :: ice
    class(atype_ice_geometry_model_data),    intent(in   ) :: geom
    type(type_bed_roughness_model),          intent(in   ) :: bed_roughness
    integer,  dimension(:), optional,        intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:), optional,        intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(:), optional,        intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_momentum_balance_solver_DIVA'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is specific to the DIVA momentum balance solver
    call self%solver%run( ice, geom, bed_roughness, &
      BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)

    self%n_visc_its = self%solver%n_visc_its
    self%n_Axb_its  = self%solver%n_Axb_its

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_DIVA_run

  subroutine momentum_balance_solver_DIVA_remap( self, mesh_old, mesh_new)

    ! In/output variables:
    class(type_momentum_balance_solver_DIVA), intent(inout) :: self
    type(type_mesh),                         intent(in   ) :: mesh_old
    type(type_mesh), target,                 intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_DIVA_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap all the stuff that is specific to the DIVA momentum balance solver
    call self%solver%remap( mesh_old, mesh_new)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_DIVA_remap

  function get_momentum_balance_solver_name( self) result( model_name)
    class(type_momentum_balance_solver_DIVA), intent(in) :: self
    character(len=:), allocatable :: model_name
    model_name = 'DIVA'
  end function get_momentum_balance_solver_name

end module momentum_balance_solver_DIVA
