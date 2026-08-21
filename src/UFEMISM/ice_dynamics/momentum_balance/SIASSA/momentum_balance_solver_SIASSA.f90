module momentum_balance_solver_SIASSA

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_basic, only: atype_ice_velocity_model
  use parameters, only: grav, ice_density, NaN
  use reallocate_mod, only: reallocate_bounds
  use constitutive_equation, only: calc_ice_rheology_Glen
  use mesh_disc_apply_operators, only: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, map_a_b_3D, ddx_a_a_2D, ddy_a_a_2D
  use mesh_zeta, only: integrate_from_zeta_is_one_to_zeta_is_zetap
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use momentum_balance_solver_SIA, only: type_momentum_balance_solver_SIA
  use momentum_balance_solver_SSA, only: type_momentum_balance_solver_SSA
  use bed_roughness_model_types, only: type_bed_roughness_model

  implicit none

  private

  public :: type_momentum_balance_solver_SIASSA

  type, extends(atype_momentum_balance_solver) :: type_momentum_balance_solver_SIASSA

      type(type_momentum_balance_solver_SIA) :: SIA
      type(type_momentum_balance_solver_SSA) :: SSA

    contains

      ! Procedures for model memory management and operation
      procedure, public :: allocate_momentum_balance_solver   => momentum_balance_solver_SIASSA_allocate
      procedure, public :: deallocate_momentum_balance_solver => momentum_balance_solver_SIASSA_deallocate
      procedure, public :: initialise_momentum_balance_solver => momentum_balance_solver_SIASSA_initialise
      procedure, public :: run_momentum_balance_solver        => momentum_balance_solver_SIASSA_run
      procedure, public :: set_velocities_to_solver_results   => momentum_balance_solver_SIASSA_set_velocities
      procedure, public :: remap_momentum_balance_solver      => momentum_balance_solver_SIASSA_remap

      procedure, public :: get_momentum_balance_solver_name

  end type type_momentum_balance_solver_SIASSA

contains

  subroutine momentum_balance_solver_SIASSA_allocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SIASSA), intent(inout) :: self
    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SIASSA_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is specific to the SIA/SSA momentum balance solver
    call self%SIA%allocate( self%region_name(), self%mesh)
    call self%SSA%allocate( self%region_name(), self%mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIASSA_allocate

  subroutine momentum_balance_solver_SIASSA_deallocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SIASSA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SIASSA_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to the SIA/SSA momentum balance solver
    call self%SIA%deallocate()
    call self%SSA%deallocate()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIASSA_deallocate

  subroutine momentum_balance_solver_SIASSA_initialise( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SIASSA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SIASSA_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise all the stuff that is specific to the SIA/SSA momentum balance solver
    call self%SIA%initialise()
    call self%SSA%initialise()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIASSA_initialise

  subroutine momentum_balance_solver_SIASSA_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)

    ! In/output variables:
    class(type_momentum_balance_solver_SIASSA), intent(inout) :: self
    class(atype_ice_model_data),                intent(inout) :: ice
    class(atype_ice_geometry_model_data),       intent(in   ) :: geom
    type(type_bed_roughness_model),             intent(in   ) :: bed_roughness
    integer,  dimension(:  ), optional,         intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:  ), optional,         intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(:  ), optional,         intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    integer,  dimension(:,:), optional,         intent(in   ) :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:,:), optional,         intent(in   ) :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
    real(dp), dimension(:,:), optional,         intent(in   ) :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_momentum_balance_solver_SIASSA'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is specific to the SIA/SSA momentum balance solver
    call self%SIA%run_momentum_balance_solver( ice, geom, bed_roughness, &
      BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    call self%SSA%run_momentum_balance_solver( ice, geom, bed_roughness, &
      BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)

    self%n_visc_its = self%SSA%n_visc_its
    self%n_Axb_its  = self%SSA%n_Axb_its

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIASSA_run

  subroutine momentum_balance_solver_SIASSA_set_velocities( self, ice, vel)

    ! In/output variables:
    class(type_momentum_balance_solver_SIASSA), intent(in   ) :: self
    class(atype_ice_model_data),                intent(inout) :: ice
    class(atype_ice_velocity_model),            intent(inout) :: vel

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SIASSA_set_velocities'
    integer                     :: ti, vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Velocities
    do ti = self%mesh%ti1, self%mesh%ti2
      vel%u_3D_b( ti,:) = self%SIA%u_3D_b( ti,:) + self%SSA%u_vav_b( ti)
      vel%v_3D_b( ti,:) = self%SIA%v_3D_b( ti,:) + self%SSA%v_vav_b( ti)
    end do

    ! Strain rates
    do vi = self%mesh%vi1, self%mesh%vi2
      vel%du_dx_3D( vi,:) = self%SSA%du_dx_a( vi  )
      vel%du_dy_3D( vi,:) = self%SSA%du_dy_a( vi  )
      vel%du_dz_3D( vi,:) = self%SIA%du_dz_3D  ( vi,:)
      vel%dv_dx_3D( vi,:) = self%SSA%dv_dx_a( vi  )
      vel%dv_dy_3D( vi,:) = self%SSA%dv_dy_a( vi  )
      vel%dv_dz_3D( vi,:) = self%SIA%dv_dz_3D  ( vi,:)
    end do

    ! In the SIA/SSA, all gradients of w are neglected
    vel%dw_dx_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp
    vel%dw_dy_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp
    ! vel%dw_dz_3D = 0._dp ! Because we now always calculate dw/dz in calc_vertical_velocities

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIASSA_set_velocities

  subroutine momentum_balance_solver_SIASSA_remap( self, mesh_old, mesh_new)

    ! In/output variables:
    class(type_momentum_balance_solver_SIASSA), intent(inout) :: self
    type(type_mesh),                            intent(in   ) :: mesh_old
    type(type_mesh), target,                    intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SIASSA_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap all the stuff that is specific to the SIA/SSA momentum balance solver
    call self%SIA%remap( mesh_old, mesh_new)
    call self%SSA%remap( mesh_old, mesh_new)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIASSA_remap

  function get_momentum_balance_solver_name( self) result( model_name)
    class(type_momentum_balance_solver_SIASSA), intent(in) :: self
    character(len=:), allocatable :: model_name
    model_name = 'SIASSA'
  end function get_momentum_balance_solver_name

end module momentum_balance_solver_SIASSA
