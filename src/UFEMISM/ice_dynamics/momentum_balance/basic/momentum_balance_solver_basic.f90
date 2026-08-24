module momentum_balance_solver_basic

  use UPSY_main, only: UPSY
  use precisions, only: dp
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use momentum_balance_solver_data, only: atype_momentum_balance_solver_data
  use mesh_types, only: type_mesh
  use parameters, only: NaN
  use reallocate_mod, only: reallocate_bounds
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use model_configuration, only: C
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_basic, only: atype_ice_velocity_model
  use bed_roughness_model_types, only: type_bed_roughness_model

  implicit none

  private

  public :: atype_momentum_balance_solver

  type, abstract, extends(atype_momentum_balance_solver_data) :: atype_momentum_balance_solver

    contains

      ! Procedures for model memory management and operation
      procedure, public :: allocate   => momentum_balance_solver_allocate
      procedure, public :: deallocate => momentum_balance_solver_deallocate
      procedure, public :: initialise => momentum_balance_solver_initialise
      procedure, public :: run        => momentum_balance_solver_run
      procedure, public :: remap      => momentum_balance_solver_remap

      ! Deferred procedures that must be overridden by each individual momentum balance implementation
      procedure(momentum_balance_solver_allocate_ifc),   deferred :: allocate_momentum_balance_solver
      procedure(momentum_balance_solver_deallocate_ifc), deferred :: deallocate_momentum_balance_solver
      procedure(momentum_balance_solver_initialise_ifc), deferred :: initialise_momentum_balance_solver
      procedure(momentum_balance_solver_run_ifc),        deferred :: run_momentum_balance_solver
      procedure(momentum_balance_solver_set_vel_ifc),    deferred :: set_velocities_to_solver_results
      procedure(momentum_balance_solver_remap_ifc),      deferred :: remap_momentum_balance_solver

      procedure, public                                         :: get_model_name
      procedure(get_momentum_balance_solver_name_ifc), deferred :: get_momentum_balance_solver_name

      procedure(create_restart_file_old_ifc),   deferred :: create_restart_file_old
      procedure(write_to_restart_file_old_ifc), deferred :: write_to_restart_file_old

  end type atype_momentum_balance_solver

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine momentum_balance_solver_allocate_ifc( self)
      import atype_momentum_balance_solver
      class(atype_momentum_balance_solver), intent(inout) :: self
    end subroutine momentum_balance_solver_allocate_ifc

    subroutine momentum_balance_solver_deallocate_ifc( self)
      import atype_momentum_balance_solver
      class(atype_momentum_balance_solver), intent(inout) :: self
    end subroutine momentum_balance_solver_deallocate_ifc

    subroutine momentum_balance_solver_initialise_ifc( self)
      import atype_momentum_balance_solver
      class(atype_momentum_balance_solver), intent(inout) :: self
    end subroutine momentum_balance_solver_initialise_ifc

    subroutine momentum_balance_solver_run_ifc( self, ice, geom, bed_roughness, &
      BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
      import atype_momentum_balance_solver, atype_ice_model_data, atype_ice_geometry_model_data, &
        type_bed_roughness_model, dp
      class(atype_momentum_balance_solver), intent(inout) :: self
      class(atype_ice_model_data),          intent(inout) :: ice
      class(atype_ice_geometry_model_data), intent(in   ) :: geom
      type(type_bed_roughness_model),       intent(in   ) :: bed_roughness
      integer,  dimension(self%mesh%ti1:self%mesh%ti2),                optional, intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2),                optional, intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2),                optional, intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
      integer,  dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), optional, intent(in   ) :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), optional, intent(in   ) :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
      real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), optional, intent(in   ) :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction
    end subroutine momentum_balance_solver_run_ifc

    subroutine momentum_balance_solver_set_vel_ifc( self, ice, vel)
      import atype_momentum_balance_solver, atype_ice_model_data, atype_ice_velocity_model
      class(atype_momentum_balance_solver), intent(in   ) :: self
      class(atype_ice_model_data),          intent(inout) :: ice
      class(atype_ice_velocity_model),      intent(inout) :: vel
    end subroutine momentum_balance_solver_set_vel_ifc

    subroutine momentum_balance_solver_remap_ifc( self, mesh_old, mesh_new)
      import atype_momentum_balance_solver, type_mesh
      class(atype_momentum_balance_solver), intent(inout) :: self
      type(type_mesh),                      intent(in   ) :: mesh_old
      type(type_mesh), target,              intent(in   ) :: mesh_new
    end subroutine momentum_balance_solver_remap_ifc

    function get_momentum_balance_solver_name_ifc( self) result( momentum_balance_solver_name)
      import atype_momentum_balance_solver
      class(atype_momentum_balance_solver), intent(in) :: self
      character(len=:), allocatable :: momentum_balance_solver_name
    end function get_momentum_balance_solver_name_ifc

    subroutine create_restart_file_old_ifc( self)
      import atype_momentum_balance_solver
      class(atype_momentum_balance_solver), intent(inout) :: self
    end subroutine create_restart_file_old_ifc

    subroutine write_to_restart_file_old_ifc( self, time)
      import atype_momentum_balance_solver, dp
      class(atype_momentum_balance_solver), intent(in   ) :: self
      real(dp),                             intent(in   ) :: time
    end subroutine write_to_restart_file_old_ifc

  end interface

contains

  recursive subroutine momentum_balance_solver_allocate( self, region_name, mesh)

    ! In/output variables:
    class(atype_momentum_balance_solver), intent(inout) :: self
    character(len=*),                     intent(in   ) :: region_name
    type(type_mesh), target,              intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all models
    call self%allocate_model( region_name, mesh)

    ! Allocate stuff that is common to all basic momentum balance solvers


    ! Allocate stuff that is specific to each individual basic momentum balance solver
    call self%allocate_momentum_balance_solver()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_allocate

  recursive subroutine momentum_balance_solver_deallocate( self)

    ! In/output variables:
    class(atype_momentum_balance_solver), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate stuff that is common to all models
    call self%deallocate_model()

    ! Deallocate stuff that is common to all basic momentum balance solvers


    ! Deallocate stuff that is specific to each individual basic momentum balance solver
    call self%deallocate_momentum_balance_solver()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_deallocate

  recursive subroutine momentum_balance_solver_initialise( self)

    ! In/output variables:
    class(atype_momentum_balance_solver), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(*,"(A)") '   Initialising ' // &
      UPSY%stru%colour_string( trim( C%choice_stress_balance_approximation),'light blue') // ' solver...'

    ! Initialise stuff that is common to all models
    call self%initialise_model()

    ! Initialise stuff that is common to all basic momentum balance solvers
    self%PETSc_rtol   = C%stress_balance_PETSc_rtol
    self%PETSc_abstol = C%stress_balance_PETSc_abstol
    self%n_visc_its   = 0
    self%n_Axb_its    = 0

    ! Initialise stuff that is specific to each individual basic momentum balance solver
    call self%initialise_momentum_balance_solver()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_initialise

  recursive subroutine momentum_balance_solver_run( self, ice, geom, bed_roughness, vel, BMB, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)

    ! In/output variables:
    class(atype_momentum_balance_solver),             intent(inout) :: self
    class(atype_ice_model_data),                      intent(inout) :: ice
    class(atype_ice_geometry_model_data),             intent(in   ) :: geom
    type(type_bed_roughness_model),                   intent(in   ) :: bed_roughness
    class(atype_ice_velocity_model),                  intent(inout) :: vel
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: BMB
    integer,  dimension(self%mesh%ti1:self%mesh%ti2),                optional, intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2),                optional, intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2),                optional, intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    integer,  dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), optional, intent(in   ) :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), optional, intent(in   ) :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), optional, intent(in   ) :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_run'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run stuff that is common to all models
    call self%run_model()

    ! Run stuff that is common to all basic momentum balance solvers

    ! Run stuff that is specific to each individual basic momentum balance solver
    call self%run_momentum_balance_solver( ice, geom, bed_roughness, &
      BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)

    ! Update the velocity model
    call self%set_velocities_to_solver_results( ice, vel)
    call vel%calc_secondary_velocities( ice, geom, BMB)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_run

  recursive subroutine momentum_balance_solver_remap( self, mesh_old, mesh_new, ice, geom, vel, BMB)

    ! In/output variables:
    class(atype_momentum_balance_solver),           intent(inout) :: self
    type(type_mesh),                                intent(in   ) :: mesh_old
    type(type_mesh), target,                        intent(in   ) :: mesh_new
    class(atype_ice_model_data),                    intent(inout) :: ice
    class(atype_ice_geometry_model_data),           intent(in   ) :: geom
    class(atype_ice_velocity_model),                intent(inout) :: vel
    real(dp), dimension(mesh_new%vi1:mesh_new%vi2), intent(in   ) :: BMB

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap stuff that is common to all models
    call self%remap_model( mesh_new)

    ! Remap stuff that is common to all basic momentum balance solvers


    ! Remap stuff that is specific to each individual basic momentum balance solver
    call self%remap_momentum_balance_solver( mesh_old, mesh_new)

    ! Update the velocity model
    call self%set_velocities_to_solver_results( ice, vel)
    call vel%calc_secondary_velocities( ice, geom, BMB)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_remap

  function get_model_name( self) result( model_name)
    class(atype_momentum_balance_solver), intent(in) :: self
    character(len=:), allocatable :: model_name
    model_name = 'momentum_balance_solver_' // self%get_momentum_balance_solver_name()
  end function get_model_name

end module momentum_balance_solver_basic