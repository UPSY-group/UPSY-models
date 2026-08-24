module momentum_balance_solver_hybrid_DIVA_BPA

  ! Routines for calculating ice velocities using the hybrid DIVA/BPA

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_BCAST, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, &
    MPI_INTEGER, MPI_LOGICAL, MPI_LOR, MPI_MAX, MPI_MIN, MPI_SUM
  use precisions, only: dp
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: warning, crash, init_routine, finalise_routine
  use model_configuration, only: C
  use parameters
  use petsc_basic, only: solve_matrix_equation_CSR_PETSc
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph_pair
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_basic, only: atype_ice_velocity_model
  use reallocate_mod, only: reallocate_bounds
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D, map_from_mesh_to_mesh_with_reallocation_3D
  use mesh_disc_apply_operators, only: map_a_b_2D, map_b_a_2D, map_b_a_3D, map_a_b_3D
  use mesh_disc_calc_matrix_operators_3D, only: calc_3D_matrix_operators_mesh
  use mesh_ROI_polygons
  use plane_geometry, only: is_in_polygon, is_in_polygons
  use mpi_distributed_memory, only: gather_to_all
  use zeta_gradients, only: calc_zeta_gradients
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use grid_basic, only: type_grid, calc_grid_mask_as_polygons
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary
  use netcdf_io_main
  use bed_roughness_model_types, only: type_bed_roughness_model
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use momentum_balance_solver_DIVA, only: type_momentum_balance_solver_DIVA
  use momentum_balance_solver_BPA, only: type_momentum_balance_solver_BPA

  implicit none

  private

  public :: type_momentum_balance_solver_hybrid_DIVA_BPA

  type, extends(atype_momentum_balance_solver) :: type_momentum_balance_solver_hybrid_DIVA_BPA

      ! Solution
      real(dp), dimension(:  ), allocatable :: u_vav_b                     ! Vertically averaged horizontal ice velocity [m yr^-1]
      real(dp), dimension(:  ), allocatable :: v_vav_b
      real(dp), dimension(:,:), allocatable :: u_bk                        ! 3-D horizontal ice velocity [m yr^-1]
      real(dp), dimension(:,:), allocatable :: v_bk

      ! DIVA and BPA solvers
      class(type_momentum_balance_solver_DIVA), allocatable :: DIVA        ! Depth-Integrated Viscosity Approximation
      class(type_momentum_balance_solver_BPA),  allocatable :: BPA         ! Blatter-Pattyn Approximation

      ! Solving masks
      logical,  dimension(:  ), allocatable :: mask_DIVA_b                 ! T: solve the DIVA here, F: otherwise
      logical,  dimension(:  ), allocatable :: mask_BPA_b                  ! T: solve the BPA  here, F: otherwise
      logical,  dimension(:  ), allocatable :: mask_3D_from_DIVA_b         ! T: calculate 3-D velocities from the vertically averaged DIVA solution here, F: otherwise
      logical,  dimension(:  ), allocatable :: mask_vav_from_BPA_b         ! T: calculate vertically averaged velocities from the 3-D BPA  solution here, F: otherwise

      ! Intermediate data fields
      real(dp), dimension(:,:), allocatable :: u_bk_prev                   ! Previous velocity solution
      real(dp), dimension(:,:), allocatable :: v_bk_prev

      ! Restart file
      character(len=256)                    :: restart_filename

    contains

      ! Procedures for model memory management and operation
      procedure, public :: get_momentum_balance_solver_name
      procedure, public :: allocate_momentum_balance_solver   => momentum_balance_solver_hybrid_DIVA_BPA_allocate
      procedure, public :: deallocate_momentum_balance_solver => momentum_balance_solver_hybrid_DIVA_BPA_deallocate
      procedure, public :: initialise_momentum_balance_solver => momentum_balance_solver_hybrid_DIVA_BPA_initialise
      procedure, public :: run_momentum_balance_solver        => momentum_balance_solver_hybrid_DIVA_BPA_run
      procedure, public :: set_velocities_to_solver_results   => momentum_balance_solver_hybrid_DIVA_BPA_set_velocities
      procedure, public :: remap_momentum_balance_solver      => momentum_balance_solver_hybrid_DIVA_BPA_remap
      procedure, public :: create_restart_file_old            => create_restart_file_old_hybrid_DIVA_BPA
      procedure, public :: write_to_restart_file_old          => write_to_restart_file_old_hybrid_DIVA_BPA

      procedure, public :: calc_hybrid_solver_masks_basic
      procedure, public :: calc_hybrid_solver_masks_basic_ROI
      procedure, public :: calc_hybrid_solver_masks_basic_read_from_file

      procedure, public :: solve_hybrid_DIVA_BPA_linearised
      procedure, public :: calc_masked_DIVA_stiffness_matrix_and_load_vector
      procedure, public :: calc_masked_BPA_stiffness_matrix_and_load_vector
      procedure, public :: calc_hybrid_solver_masks_transition
      procedure, public :: calc_hybrid_solver_translation_tables

      procedure, public :: apply_velocity_limits
      procedure, public :: relax_viscosity_iterations
      procedure, public :: calc_visc_iter_UV_resid

  end type type_momentum_balance_solver_hybrid_DIVA_BPA

contains

! == Main routines

  subroutine momentum_balance_solver_hybrid_DIVA_BPA_allocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_hybrid_DIVA_BPA_allocate'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    allocate( self%u_vav_b( self%mesh%ti1:self%mesh%ti2               ), source = 0._dp)
    allocate( self%v_vav_b( self%mesh%ti1:self%mesh%ti2               ), source = 0._dp)
    allocate( self%u_bk   ( self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), source = 0._dp)
    allocate( self%v_bk   ( self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), source = 0._dp)

    ! Separate DIVA/BPA solvers
    call self%DIVA%allocate( self%region_name(), self%mesh)
    call self%BPA %allocate( self%region_name(), self%mesh)

    ! Solver masks
    allocate( self%mask_DIVA_b        ( self%mesh%ti1:self%mesh%ti2), source = .false.)
    allocate( self%mask_BPA_b         ( self%mesh%ti1:self%mesh%ti2), source = .false.)
    allocate( self%mask_3D_from_DIVA_b( self%mesh%ti1:self%mesh%ti2), source = .false.)
    allocate( self%mask_vav_from_BPA_b( self%mesh%ti1:self%mesh%ti2), source = .false.)

    ! Intermediate data fields
    allocate( self%u_bk_prev( self%mesh%nTri,self%mesh%nz), source = 0._dp)
    allocate( self%v_bk_prev( self%mesh%nTri,self%mesh%nz), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_hybrid_DIVA_BPA_allocate

  subroutine momentum_balance_solver_hybrid_DIVA_BPA_deallocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_hybrid_DIVA_BPA_deallocate'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    deallocate( self%u_vav_b)
    deallocate( self%v_vav_b)
    deallocate( self%u_bk)
    deallocate( self%v_bk)

    ! Separate DIVA/BPA solvers
    call self%DIVA%deallocate()
    call self%BPA %deallocate()

    ! Solver masks
    deallocate( self%mask_DIVA_b)
    deallocate( self%mask_BPA_b)
    deallocate( self%mask_3D_from_DIVA_b)
    deallocate( self%mask_vav_from_BPA_b)

    ! Intermediate data fields
    deallocate( self%u_bk_prev)
    deallocate( self%v_bk_prev)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_hybrid_DIVA_BPA_deallocate

  subroutine momentum_balance_solver_hybrid_DIVA_BPA_initialise( self)

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'momentum_balance_solver_hybrid_DIVA_BPA_initialise'
    character(len=:), allocatable :: choice_initial_velocity

    ! Add routine to path
    call init_routine( routine_name)

    call crash('FIXME')

    ! Determine the choice of initial velocities for this model region
    select case (self%region_name())
    case default
      call crash('unknown model region "' // self%region_name() // '"!')
    case ('NAM')
      choice_initial_velocity  = C%choice_initial_velocity_NAM
    case ('EAS')
      choice_initial_velocity  = C%choice_initial_velocity_EAS
    case ('GRL')
      choice_initial_velocity  = C%choice_initial_velocity_GRL
    case ('ANT')
      choice_initial_velocity  = C%choice_initial_velocity_ANT
    end select

    ! Initialise velocities according to the specified method
    select case (choice_initial_velocity)
    case default
      call crash('unknown choice_initial_velocity "' // trim( choice_initial_velocity) // '"!')
    case ('zero')
      self%u_vav_b( self%mesh%ti1:self%mesh%ti2  ) = 0._dp
      self%v_vav_b( self%mesh%ti1:self%mesh%ti2  ) = 0._dp
      self%u_bk   ( self%mesh%ti1:self%mesh%ti2,:) = 0._dp
      self%v_bk   ( self%mesh%ti1:self%mesh%ti2,:) = 0._dp
    case ('read_from_file')
      call crash('restarting ice velocities not yet possible for the hybrid DIVA/BPA!')
    end select

    ! Set tolerances for PETSc matrix solver for the linearised hybrid DIVA/BPA
    self%PETSc_rtol   = C%stress_balance_PETSc_rtol
    self%PETSc_abstol = C%stress_balance_PETSc_abstol

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_hybrid_DIVA_BPA_initialise

  subroutine momentum_balance_solver_hybrid_DIVA_BPA_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate ice velocities by solving the hybrid DIVA/BPA

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self
    class(atype_ice_model_data),                         intent(inout) :: ice
    class(atype_ice_geometry_model_data),                intent(in   ) :: geom
    type(type_bed_roughness_model),                      intent(in   ) :: bed_roughness
    integer,  dimension(self%mesh%ti1:self%mesh%ti2),                optional, intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2),                optional, intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2),                optional, intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    integer,  dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), optional, intent(in   ) :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), optional, intent(in   ) :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), optional, intent(in   ) :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'solve_hybrid_DIVA_BPA'
    logical                             :: grounded_ice_exists
    integer                             :: ierr
    integer,  dimension(:), allocatable :: BC_prescr_mask_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_u_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_v_b_applied
    integer                             :: viscosity_iteration_i
    logical                             :: has_converged
    real(dp)                            :: resid_UV, resid_UV_prev
    real(dp)                            :: uv_min, uv_max
    real(dp)                            :: visc_it_relax_applied
    real(dp)                            :: Glens_flow_law_epsilon_sq_0_applied
    integer                             :: nit_diverg_consec
    integer                             :: n_Axb_its_visc_it

    ! Add routine to path
    call init_routine( routine_name)

    ! if there is no grounded ice, no need (in fact, no way) to solve the velocities
    grounded_ice_exists = ANY( geom%mask_grounded_ice)
    call MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.not. grounded_ice_exists) then
      self%u_vav_b( self%mesh%ti1:self%mesh%ti2  ) = 0._dp
      self%v_vav_b( self%mesh%ti1:self%mesh%ti2  ) = 0._dp
      self%u_bk   ( self%mesh%ti1:self%mesh%ti2,:) = 0._dp
      self%v_bk   ( self%mesh%ti1:self%mesh%ti2,:) = 0._dp
      call finalise_routine( routine_name)
      return
    end if

    ! Handle the optional prescribed u,v boundary conditions
    allocate( BC_prescr_mask_b_applied( self%mesh%ti1:self%mesh%ti2))
    allocate( BC_prescr_u_b_applied(    self%mesh%ti1:self%mesh%ti2))
    allocate( BC_prescr_v_b_applied(    self%mesh%ti1:self%mesh%ti2))
    if (present( BC_prescr_mask_b) .or. present( BC_prescr_u_b) .or. present( BC_prescr_v_b)) then
      ! Safety
      if (.not. (present( BC_prescr_mask_b) .and. present( BC_prescr_u_b) .and. present( BC_prescr_v_b))) then
        call crash('need to provide prescribed u,v fields and mask!')
      end if
      BC_prescr_mask_b_applied( self%mesh%ti1:self%mesh%ti2) = BC_prescr_mask_b( self%mesh%ti1:self%mesh%ti2)
      BC_prescr_u_b_applied   ( self%mesh%ti1:self%mesh%ti2) = BC_prescr_u_b   ( self%mesh%ti1:self%mesh%ti2)
      BC_prescr_v_b_applied   ( self%mesh%ti1:self%mesh%ti2) = BC_prescr_v_b   ( self%mesh%ti1:self%mesh%ti2)
    else
      BC_prescr_mask_b_applied( self%mesh%ti1:self%mesh%ti2) = 0
      BC_prescr_u_b_applied   ( self%mesh%ti1:self%mesh%ti2) = 0._dp
      BC_prescr_v_b_applied   ( self%mesh%ti1:self%mesh%ti2) = 0._dp
    end if

    ! Calculate zeta gradients
    call calc_zeta_gradients( self%mesh, ice, geom)

    ! Calculate 3-D matrix operators for the current ice geometry
    call calc_3D_matrix_operators_mesh( self%mesh, &
      ice%dzeta_dx_ak, ice%dzeta_dy_ak, ice%dzeta_dx_bk, ice%dzeta_dy_bk, &
      ice%dzeta_dz_bk, ice%dzeta_dz_bks, &
      ice%d2zeta_dx2_bk, ice%d2zeta_dxdy_bk, ice%d2zeta_dy2_bk)

    ! Calculate the driving stress
    call self%DIVA%calc_driving_stress( geom)
    call self%BPA %calc_driving_stress( geom)

    ! Calculate the solving masks for the hybrid solver
    call self%calc_hybrid_solver_masks_basic()

    ! Adaptive relaxation parameter for the viscosity iteration
    resid_UV                            = 1E9_dp
    nit_diverg_consec                   = 0
    visc_it_relax_applied               = C%visc_it_relax
    Glens_flow_law_epsilon_sq_0_applied = C%Glens_flow_law_epsilon_sq_0

    ! Initialise stability info
    self%n_visc_its = 0
    self%n_Axb_its  = 0

    ! The viscosity iteration
    viscosity_iteration_i = 0
    has_converged         = .false.
    viscosity_iteration: do while (.not. has_converged)
      viscosity_iteration_i = viscosity_iteration_i + 1

      ! == Calculate secondary terms in the DIVA
      ! ========================================

      ! Calculate the horizontal strain rates for the current velocity solution
      call self%DIVA%calc_horizontal_strain_rates()

      ! Calculate the vertical shear strain rates
      call self%DIVA%calc_vertical_shear_strain_rates()

      ! Calculate the effective viscosity for the current velocity solution
      call self%DIVA%calc_effective_viscosity( ice, geom, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the F-integrals
      call self%DIVA%calc_F_integrals( geom)

      ! Calculate the "effective" friction coefficient (turning the SSA into the DIVA)
      call self%DIVA%calc_effective_basal_friction_coefficient( ice, geom, bed_roughness)

      ! == Calculate secondary terms in the BPA
      ! =======================================

      ! Calculate the strain rates for the current velocity solution
      call self%BPA%calc_strain_rates()

      ! Calculate the effective viscosity for the current velocity solution
      call self%BPA%calc_effective_viscosity( ice, geom, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the basal friction coefficient betab for the current velocity solution
      call self%BPA%calc_applied_basal_friction_coefficient( ice, geom, bed_roughness)

      ! == Solve the linearised hybrid DIVA/BPA
      ! =======================================

      ! Solve the linearised hybrid DIVA/BPA to calculate a new velocity solution
      call self%solve_hybrid_DIVA_BPA_linearised( ice, n_Axb_its_visc_it, &
        BC_prescr_mask_b_applied, BC_prescr_u_b_applied, BC_prescr_v_b_applied)

      ! Update stability info
      self%n_Axb_its = self%n_Axb_its + n_Axb_its_visc_it

      ! Copy results to the DIVA and BPA structures
      self%DIVA%u_vav_b( self%mesh%ti1:self%mesh%ti2  ) = self%u_vav_b( self%mesh%ti1:self%mesh%ti2  )
      self%DIVA%v_vav_b( self%mesh%ti1:self%mesh%ti2  ) = self%v_vav_b( self%mesh%ti1:self%mesh%ti2  )
      self%DIVA%u_3D_b ( self%mesh%ti1:self%mesh%ti2,:) = self%u_bk   ( self%mesh%ti1:self%mesh%ti2,:)
      self%DIVA%v_3D_b ( self%mesh%ti1:self%mesh%ti2,:) = self%v_bk   ( self%mesh%ti1:self%mesh%ti2,:)
      self%BPA%u_bk    ( self%mesh%ti1:self%mesh%ti2,:) = self%u_bk   ( self%mesh%ti1:self%mesh%ti2,:)
      self%BPA%v_bk    ( self%mesh%ti1:self%mesh%ti2,:) = self%v_bk   ( self%mesh%ti1:self%mesh%ti2,:)

      ! == Calculate more secondary terms in the DIVA
      ! =============================================

      ! Calculate basal velocities
      call self%DIVA%calc_basal_velocities()

      ! Calculate basal shear stress
      call self%DIVA%calc_basal_shear_stress()

      ! == Improve stability and check for convergence
      ! ==============================================

      ! Limit velocities for improved stability
      call self%apply_velocity_limits()

      ! Reduce the change between velocity solutions
      call self%relax_viscosity_iterations( visc_it_relax_applied)

      ! Calculate the L2-norm of the two consecutive velocity solutions
      resid_UV_prev = resid_UV
      call self%calc_visc_iter_UV_resid( resid_UV)

      ! if the viscosity iteration diverges, lower the relaxation parameter
      if (resid_UV > resid_UV_prev) then
        nit_diverg_consec = nit_diverg_consec + 1
      else
        nit_diverg_consec = 0
      end if
      if (nit_diverg_consec > 2) then
        nit_diverg_consec = 0
        visc_it_relax_applied               = visc_it_relax_applied               * 0.9_dp
        Glens_flow_law_epsilon_sq_0_applied = Glens_flow_law_epsilon_sq_0_applied * 1.2_dp
      end if
      if (visc_it_relax_applied <= 0.05_dp .or. Glens_flow_law_epsilon_sq_0_applied >= 1E-5_dp) then
        if (visc_it_relax_applied < 0.05_dp) then
          call crash('viscosity iteration still diverges even with very low relaxation factor!')
        elseif (Glens_flow_law_epsilon_sq_0_applied > 1E-5_dp) then
          call crash('viscosity iteration still diverges even with very high effective strain rate regularisation!')
        end if
      end if

      ! DENK DROM
     uv_min = MINVAL( self%u_bk( self%mesh%ti1:self%mesh%ti2,:))
     uv_max = MAXVAL( self%u_bk( self%mesh%ti1:self%mesh%ti2,:))
     call MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
    !  if (par%primary) WRITE(0,*) '    hybrid DIVA/BPA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], resid = ', resid_UV

      ! if the viscosity iteration has converged, or has reached the maximum allowed number of iterations, stop it.
      has_converged = .false.
      if (resid_UV < C%visc_it_norm_dUV_tol) then
        has_converged = .true.
      end if

      ! if we've reached the maximum allowed number of iterations without converging, throw a warning
      if (viscosity_iteration_i > C%visc_it_nit) then
        if (par%primary) call warning('viscosity iteration failed to converge within {int_01} iterations!', int_01 = C%visc_it_nit)
        exit viscosity_iteration
      end if

    end do viscosity_iteration

    ! Stability info
    self%n_visc_its = viscosity_iteration_i

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_hybrid_DIVA_BPA_run

  subroutine momentum_balance_solver_hybrid_DIVA_BPA_set_velocities( self, ice, vel)

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(in   ) :: self
    class(atype_ice_model_data),                         intent(inout) :: ice
    class(atype_ice_velocity_model),                     intent(inout) :: vel

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_hybrid_DIVA_BPA_set_velocities'
    integer                     :: ti, vi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Velocities
    do ti = self%mesh%ti1, self%mesh%ti2
      vel%u_3D_b( ti,:) = self%u_bk( ti,:)
      vel%v_3D_b( ti,:) = self%v_bk( ti,:)
    end do

    ! Strain rates
    do vi = self%mesh%vi1, self%mesh%vi2
      vel%du_dx_3D( vi,:) = self%BPA%du_dx_ak( vi,:)
      vel%du_dy_3D( vi,:) = self%BPA%du_dy_ak( vi,:)
      vel%du_dz_3D( vi,:) = self%BPA%du_dz_ak( vi,:)
      vel%dv_dx_3D( vi,:) = self%BPA%dv_dx_ak( vi,:)
      vel%dv_dy_3D( vi,:) = self%BPA%dv_dy_ak( vi,:)
      vel%dv_dz_3D( vi,:) = self%BPA%dv_dz_ak( vi,:)
    end do

    ! In the hybrid DIVA/BPA, gradients of w are neglected
    vel%dw_dx_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp
    vel%dw_dy_3D( self%mesh%vi1:self%mesh%vi2,:) = 0._dp
    ! vel%dw_dz_3D = 0._dp ! Because we now always calculate dw/dz in calc_vertical_velocities

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_hybrid_DIVA_BPA_set_velocities

  subroutine momentum_balance_solver_hybrid_DIVA_BPA_remap( self, mesh_old, mesh_new)

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self
    type(type_mesh),                                     intent(in   ) :: mesh_old
    type(type_mesh), target,                             intent(in   ) :: mesh_new

    ! Local variables:
    character(*), parameter               :: routine_name = 'remap_hybrid_DIVA_BPA_solver'
    real(dp), dimension(:  ), allocatable :: u_vav_a
    real(dp), dimension(:  ), allocatable :: v_vav_a
    real(dp), dimension(:,:), allocatable :: u_ak
    real(dp), dimension(:,:), allocatable :: v_ak

    ! Add routine to path
    call init_routine( routine_name)

    call crash('FIXME')

    ! ! Remap the fields that are re-used during the viscosity iteration
    ! ! ================================================================

    ! ! allocate memory for velocities on the a-grid (vertices)
    ! allocate( u_vav_a ( mesh_old%vi1: mesh_old%vi2             ))
    ! allocate( v_vav_a ( mesh_old%vi1: mesh_old%vi2             ))
    ! allocate( u_ak    ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    ! allocate( v_ak    ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))

    ! ! Map data from the triangles of the old mesh to the vertices of the old mesh
    ! call map_b_a_2D( mesh_old, self%u_vav_b, u_vav_a)
    ! call map_b_a_2D( mesh_old, self%v_vav_b, v_vav_a)
    ! call map_b_a_3D( mesh_old, self%u_bk   , u_ak   )
    ! call map_b_a_3D( mesh_old, self%v_bk   , v_ak   )

    ! ! Remap data from the vertices of the old mesh to the vertices of the new mesh
    ! call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, u_vav_a, '2nd_order_conservative')
    ! call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, v_vav_a, '2nd_order_conservative')
    ! call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, u_ak   , '2nd_order_conservative')
    ! call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, v_ak   , '2nd_order_conservative')

    ! ! reallocate memory for the data on the triangles
    ! call reallocate_bounds( self%u_vav_b, mesh_new%ti1, mesh_new%ti2             )
    ! call reallocate_bounds( self%v_vav_b, mesh_new%ti1, mesh_new%ti2             )
    ! call reallocate_bounds( self%u_bk   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    ! call reallocate_bounds( self%v_bk   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)

    ! ! Map data from the vertices of the new mesh to the triangles of the new mesh
    ! call map_a_b_2D( mesh_new, u_vav_a, self%u_vav_b)
    ! call map_a_b_2D( mesh_new, v_vav_a, self%v_vav_b)
    ! call map_a_b_3D( mesh_new, u_ak   , self%u_bk   )
    ! call map_a_b_3D( mesh_new, v_ak   , self%v_bk   )

    ! ! Remap data of the separate DIVA and BPA solvers
    ! ! ===============================================

    ! call self%DIVA%remap( mesh_old, mesh_new)
    ! call self%BPA %remap( mesh_old, mesh_new)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_hybrid_DIVA_BPA_remap

  function get_momentum_balance_solver_name( self) result( model_name)
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(in) :: self
    character(len=:), allocatable :: model_name
    model_name = 'hybrid_DIVA_BPA'
  end function get_momentum_balance_solver_name

! == Basic masks and translation tables for the hybrid solver

  subroutine calc_hybrid_solver_masks_basic( self)
    !< Calculate the solving masks for the hybrid solver

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'calc_hybrid_solver_masks_basic'
    character(len=:), allocatable :: choice_hybrid_DIVA_BPA_mask

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    self%mask_BPA_b  = .false.
    self%mask_DIVA_b = .false.

    ! Determine filename for this model region
    select case (self%region_name())
      case default
        call crash('unknown region name "' // trim( self%region_name()) // '"!')
      case ('NAM')
        choice_hybrid_DIVA_BPA_mask = C%choice_hybrid_DIVA_BPA_mask_NAM
      case ('EAS')
        choice_hybrid_DIVA_BPA_mask = C%choice_hybrid_DIVA_BPA_mask_EAS
      case ('GRL')
        choice_hybrid_DIVA_BPA_mask = C%choice_hybrid_DIVA_BPA_mask_GRL
      case ('ANT')
        choice_hybrid_DIVA_BPA_mask = C%choice_hybrid_DIVA_BPA_mask_ANT
    end select

    select case (choice_hybrid_DIVA_BPA_mask)
      case default
        call crash('unknown choice_hybrid_DIVA_BPA_mask "' // trim( choice_hybrid_DIVA_BPA_mask) // '"!')
      case ('ROI')
        call self%calc_hybrid_solver_masks_basic_ROI()
      case ('read_from_file')
        call self%calc_hybrid_solver_masks_basic_read_from_file()
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_hybrid_solver_masks_basic

  subroutine calc_hybrid_solver_masks_basic_ROI( self)
    !< Calculate the solving masks for the hybrid solver

    ! Solve the BPA only in the specified regions of interest

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'calc_hybrid_solver_masks_basic_ROI'
    character(len=:), allocatable         :: all_names_ROI, name_ROI
    integer                               :: i
    real(dp), dimension(:,:), allocatable :: poly_ROI
    integer                               :: ti
    real(dp), dimension(2)                :: p

    ! Add routine to path
    call init_routine( routine_name)

    ! Go over all listed regions of interest
    all_names_ROI = C%choice_regions_of_interest

    do while (.true.)

      ! Get the first region of interest from the list
      i = index( all_names_ROI, '||')
      if (i == 0) then
        ! There is only one left in the list
        name_ROI = trim( all_names_ROI)
        all_names_ROI = ''
      else
        ! Get the first first one from the list and remove it
        name_ROI = all_names_ROI( 1:i-1)
        all_names_ROI = all_names_ROI( i+2:LEN_trim( all_names_ROI))
      end if


      select case (self%region_name())
        case default
          call crash('unknown region name "' // self%region_name() // '"!')
        case ('NAM')
          ! North america

          select case (name_ROI)
            case default
              call crash('unknown region of interest "' // trim( name_ROI) // '"!')
            case ('')
              ! Don't need to do anything
              exit
            case ('PineIsland')
              ! Don't need to do anything
              exit
            case ('Thwaites')
              ! Don't need to do anything
              exit
          end select

        case ('EAS')
          ! Eurasia

          select case (name_ROI)
            case default
              call crash('unknown region of interest "' // trim( name_ROI) // '"!')
            case ('')
              ! Don't need to do anything
              exit
            case ('PineIsland')
              ! Don't need to do anything
              exit
            case ('Thwaites')
              ! Don't need to do anything
              exit
          end select

        case ('GRL')
          ! Greenland

          select case (name_ROI)
            case default
              call crash('unknown region of interest "' // trim( name_ROI) // '"!')
            case ('')
              ! Don't need to do anything
              exit
            case ('PineIsland')
              ! Don't need to do anything
              exit
            case ('Thwaites')
              ! Don't need to do anything
              exit
          end select

        case ('ANT')

          select case (name_ROI)
            case default
              call crash('unknown region of interest "' // trim( name_ROI) // '"!')
            case ('')
              ! Don't need to do anything
              exit
            case ('PineIsland')
              call calc_polygon_Pine_Island_Glacier( poly_ROI)
            case ('Thwaites')
              call calc_polygon_Thwaites_Glacier( poly_ROI)
          end select

      end select

      ! Find all triangles that lie within this region of interest
      do ti = self%mesh%ti1, self%mesh%ti2
        p = self%mesh%TriGC( ti,:)
        if (is_in_polygon( poly_ROI, p)) then
          self%mask_BPA_b(  ti) = .true.
          self%mask_DIVA_b( ti) = .false.
        else
          self%mask_BPA_b(  ti) = .false.
          self%mask_DIVA_b( ti) = .true.
        end if
      end DO

      ! Clean up after yourself
      deallocate( poly_ROI)

      ! if no names are left, we are finished
      if (all_names_ROI == '') exit

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_hybrid_solver_masks_basic_ROI

  subroutine calc_hybrid_solver_masks_basic_read_from_file( self)
    !< Calculate the solving masks for the hybrid solver

    ! Read the mask that determines where to solve the DIVA
    ! and where to solve the BPA from an external NetCDF file

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'calc_hybrid_solver_masks_basic_read_from_file'
    character(len=:), allocatable         :: filename_hybrid_DIVA_BPA_mask
    integer                               :: ncid
    type(type_grid)                       :: grid
    integer,  dimension(:  ), allocatable :: mask_int_grid_vec_partial
    integer,  dimension(:,:), allocatable :: mask_int_grid
    logical,  dimension(:,:), allocatable :: mask_grid
    integer                               :: i,j,ierr
    real(dp), dimension(:,:), allocatable :: poly_mult
    integer                               :: ti
    real(dp), dimension(2)                :: p

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine filename for this model region
    select case (self%region_name())
      case default
        call crash('unknown region name "' // trim( self%region_name()) // '"!')
      case ('NAM')
        filename_hybrid_DIVA_BPA_mask = C%filename_hybrid_DIVA_BPA_mask_NAM
      case ('EAS')
        filename_hybrid_DIVA_BPA_mask = C%filename_hybrid_DIVA_BPA_mask_EAS
      case ('GRL')
        filename_hybrid_DIVA_BPA_mask = C%filename_hybrid_DIVA_BPA_mask_GRL
      case ('ANT')
        filename_hybrid_DIVA_BPA_mask = C%filename_hybrid_DIVA_BPA_mask_ANT
    end select

    ! Read grid from file
    call open_existing_netcdf_file_for_reading( filename_hybrid_DIVA_BPA_mask, ncid)
    call setup_xy_grid_from_file( filename_hybrid_DIVA_BPA_mask, ncid, grid)
    call close_netcdf_file( ncid)

    ! Read gridded mask from file
    allocate( mask_int_grid_vec_partial( grid%n1: grid%n2))
    call read_field_from_xy_file_int_2D( filename_hybrid_DIVA_BPA_mask, 'mask_BPA', mask_int_grid_vec_partial)

    ! Gather partial gridded data to the primary and broadcast the total field to all processes
    allocate( mask_int_grid( grid%nx, grid%ny))
    call gather_gridded_data_to_primary( grid, mask_int_grid_vec_partial, mask_int_grid)
    call MPI_BCAST( mask_int_grid(:,:), grid%nx * grid%ny, MPI_integer, 0, MPI_COMM_WORLD, ierr)

    ! Calculate logical mask (assumes data from file is integer 0 for FALSE and integer 1 for true)
    allocate( mask_grid( grid%nx, grid%ny), source = .false.)
    do i = 1, grid%nx
    do j = 1, grid%ny
      if (mask_int_grid( i,j) == 1) mask_grid( i,j) = .true.
    end DO
    end DO

    ! Calculate contour from gridded mask
    call calc_grid_mask_as_polygons( grid, mask_grid, poly_mult)

    ! Determine BPA solving masks on the mesh
    do ti = self%mesh%ti1, self%mesh%ti2
      p = self%mesh%TriGC( ti,:)
      if (is_in_polygons( poly_mult, p)) then
        self%mask_BPA_b(  ti) = .true.
        self%mask_DIVA_b( ti) = .false.
      else
        self%mask_BPA_b(  ti) = .false.
        self%mask_DIVA_b( ti) = .true.
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_hybrid_solver_masks_basic_read_from_file

! == Assemble and solve the linearised hybrid DIVA/BPA

  subroutine solve_hybrid_DIVA_BPA_linearised( self, ice, n_Axb_its, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    !< Solve the linearised hybrid DIVA/BPA

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self
    class(atype_ice_model_data),                         intent(in   ) :: ice
    integer,                                             intent(  out) :: n_Axb_its             ! Number of iterations used in the iterative solver
    integer,  dimension(self%mesh%ti1:self%mesh%ti2),    intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2),    intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2),    intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter             :: routine_name = 'solve_hybrid_DIVA_BPA_linearised'
    type(type_CSR_matrix_dp)                :: A_DIVA, A_BPA
    real(dp), dimension(:    ), allocatable :: b_DIVA, b_BPA
    integer,  dimension(:,:  ), allocatable :: tiuv2nh
    integer,  dimension(:,:,:), allocatable :: tikuv2nh
    integer,  dimension(:,:  ), allocatable :: nh2tiuv_tikuv
    integer                                 :: neq,i1,i2
    integer                                 :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    type(type_CSR_matrix_dp)                :: A_combi
    real(dp), dimension(:    ), allocatable :: b_combi
    real(dp), dimension(:    ), allocatable :: uv_combi
    integer                                 :: neq_loc
    integer                                 :: row_nh,ti,k,uv,row_tiuv,row_tikuv,kk1,kk2,kk,col_tiuv,col_tikuv,tin,kn,uvn,col_nh
    real(dp)                                :: val, dzeta
    integer                                 :: nhu, nhv

    ! Add routine to path
    call init_routine( routine_name)

    ! Store the previous solution
    self%u_bk_prev( self%mesh%ti1:self%mesh%ti2,:) = self%u_bk( self%mesh%ti1:self%mesh%ti2,:)
    self%v_bk_prev( self%mesh%ti1:self%mesh%ti2,:) = self%v_bk( self%mesh%ti1:self%mesh%ti2,:)

    ! Calculate the stiffness matrix and load vector for the DIVA and the BPA
    call self%calc_masked_DIVA_stiffness_matrix_and_load_vector(      BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, A_DIVA, b_DIVA)
    call self%calc_masked_BPA_stiffness_matrix_and_load_vector ( ice, BC_prescr_mask_b                              , A_BPA , b_BPA )

    ! Calculate the "transition" solver masks
    call self%calc_hybrid_solver_masks_transition( A_DIVA, A_BPA)

    ! Calculate combined DIVA/BPA translation tables
    call self%calc_hybrid_solver_translation_tables( tiuv2nh, tikuv2nh, nh2tiuv_tikuv, neq, i1, i2)
    neq_loc = i2 + 1 - i1

    ! == Construct combined stiffness matrix and load vector
    ! ======================================================

    ! Initialise the stiffness matrix using the native UFEMISM CSR-matrix format

    ! Matrix size
    ncols           = neq      ! from
    ncols_loc       = neq_loc
    nrows           = neq      ! to
    nrows_loc       = neq_loc
    nnz_est_proc    = ceiling( 1.1_dp * real( A_DIVA%nnz + A_BPA%nnz, dp))

    call A_combi%allocate( nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for the load vector and the solution
    allocate( b_combi(  i1:i2))
    allocate( uv_combi( i1:i2))

    do row_nh = i1, i2

      if (nh2tiuv_tikuv( row_nh,1) == 1) then
        ! This row represents a vertically averaged velocity

        ! This row in the combined matrix corresponds to this triangle and vertically averaged velocity component
        ti = nh2tiuv_tikuv( row_nh,2)
        uv = nh2tiuv_tikuv( row_nh,4)

        if (self%mask_DIVA_b( ti)) then
          ! Copy the corresponding row from the DIVA stiffness matrix

          ! This triangle and vertically averaged velocity component correspond to this row in the DIVA matrix
          row_tiuv = self%mesh%tiuv2n( ti,uv)

          ! This row in the DIVA matrix contains these columns
          kk1 = A_DIVA%ptr( row_tiuv)
          kk2 = A_DIVA%ptr( row_tiuv+1)-1

          ! Loop over the columns of this row of the DIVA matrix
          do kk = kk1, kk2
            ! This column index and coefficient of this entry in the DIVA matrix
            col_tiuv = A_DIVA%ind( kk)
            val      = A_DIVA%val( kk)
            ! This column in the DIVA matrix corresponds to this neighbouring triangle and vertically averaged velocity component
            tin = self%mesh%n2tiuv( col_tiuv,1)
            uvn = self%mesh%n2tiuv( col_tiuv,2)
            ! This neighbouring triangle and vertically averaged velocity component corresponds to this column in the combined matrix
            col_nh = tiuv2nh( tin,uvn)
            ! Add the coefficient from the DIVA matrix to the combined matrix
            call A_combi%add_entry( row_nh, col_nh, val)
          end do ! do kk = kk1, kk2

          ! Copy the DIVA load vector
          b_combi( row_nh) = b_DIVA( row_tiuv)

          ! Take the previous velocity solution as the initial guess
          uv_combi( row_nh) = self%DIVA%u_vav_b( ti)

        elseif (self%mask_vav_from_BPA_b( ti)) then
          ! Define the vertically averaged velocities here from the 3-D BPA velocities
          !
          ! -u_vav + SUM_k [ u_3D( k) * dzeta( k)] = 0

          ! Add the coefficient of -1 for the vertically averaged velocity to the combined matrix
          call A_combi%add_entry( row_nh, row_nh, -1._dp)

          ! Loop over the vertical column
          do k = 1, self%mesh%nz

            ! Calculate the weight dzeta for the vertical average
            if     (k == 1) then
              dzeta = self%mesh%zeta_stag( 1)
            elseif (k == self%mesh%nz) then
              dzeta = 1._dp - self%mesh%zeta_stag( self%mesh%nz-1)
            else
              dzeta = self%mesh%zeta_stag( k) - self%mesh%zeta_stag( k-1)
            end if

            ! The 3-D velocity for this layer in this triangle corresponds to this column in the combined matrix
            col_nh = tikuv2nh( ti,k,uv)

            ! Add the coefficient to the combined matrix
            call A_combi%add_entry( row_nh, col_nh, dzeta)

          end do ! do k = 1, mesh%nz

          ! The load vector is zero in this case
          b_combi( row_nh) = 0._dp

          ! Take the previous velocity solution as the initial guess
          uv_combi( row_nh) = self%DIVA%u_vav_b( ti)

        else
          call crash('mask inconsistency; expected vertically averaged velocities, but both mask_DIVA_b and mask_vav_from_BPA_b are false!')
        end if

      elseif (nh2tiuv_tikuv( row_nh,1) == 2) then
        ! This row represents a 3-D velocity

        ! This row in the combined matrix corresponds to this triangle, layer, and 3-D velocity component
        ti = nh2tiuv_tikuv( row_nh,2)
        k  = nh2tiuv_tikuv( row_nh,3)
        uv = nh2tiuv_tikuv( row_nh,4)

        if     (self%mask_BPA_b( ti)) then
          ! Copy the corresponding row from the BPA stiffness matrix

          ! This triangle, layer, and 3-D velocity component correspond to this row in the BPA matrix
          row_tikuv = self%mesh%tikuv2n( ti,k,uv)

          ! This row in the BPA matrix contains these columns
          kk1 = A_BPA%ptr( row_tikuv)
          kk2 = A_BPA%ptr( row_tikuv+1)-1

          ! Loop over the columns of this row of the BPA matrix
          do kk = kk1, kk2
            ! This column index and coefficient of this entry in the BPA matrix
            col_tikuv = A_BPA%ind( kk)
            val       = A_BPA%val( kk)
            ! This column in the BPA matrix corresponds to this neighbouring triangle, layer, and 3-D velocity component
            tin = self%mesh%n2tikuv( col_tikuv,1)
            kn  = self%mesh%n2tikuv( col_tikuv,2)
            uvn = self%mesh%n2tikuv( col_tikuv,3)
            ! This neighbouring triangle, layer, and 3-D velocity component corresponds to this column in the combined matrix
            col_nh = tikuv2nh( tin,kn,uvn)
            ! Add the coefficient from the DIVA matrix to the combined matrix
            call A_combi%add_entry( row_nh, col_nh, val)
          end do ! do kk = kk1, kk2

          ! Copy the BPA load vector
          b_combi( row_nh) = b_BPA( row_tikuv)

          ! Take the previous velocity solution as the initial guess
          uv_combi( row_nh) = self%BPA%u_bk( ti,k)

        elseif (self%mask_3D_from_DIVA_b( ti)) then
          ! Define the 3-D velocities here from the vertically averaged DIVA velocities

          if (C%choice_sliding_law == 'no_sliding') then
            ! Exception for the case of no sliding
            !
            ! According to Lipscomb et al., 2019, Eq. 29, and text between Eqs. 33 and 34:
            !
            !   [1] u( z) = tau_bx * F1( z)
            !
            ! Also, according to the text just above Eq. 33:
            !
            !   [2] tau_bx = u_vav * beta_eff
            !
            ! Substituting [2] into [1] yields:
            !
            !   [3] u( z) = u_vav * beta_eff * F1( z)
            !
            ! This can be rearranged to read:
            !
            !   [4] u_vav * beta_eff * F1( z) - u( z) = 0

            ! u_vav term
            col_nh = tiuv2nh( ti,uv)
            val = self%DIVA%beta_eff_b( ti) * self%DIVA%F1_3D_b( ti,k)
            call A_combi%add_entry( row_nh, col_nh, val)

            ! u( z) term
            val = -1._dp
            call A_combi%add_entry( row_nh, row_nh, val)

            ! The load vector is zero in this case
            b_combi( row_nh) = 0._dp

            ! Take the previous velocity solution as the initial guess
            uv_combi( row_nh) = self%BPA%u_bk( ti,k)

          else ! if (C%choice_sliding_law == 'no_sliding') then
            ! The case default of finite sliding
            !
            ! According to Lipscomb et al., 2019, Eq. 29:
            !
            !   [1] u( z) = u_b * (1 + beta * F1( z))
            !
            ! Also, according to Eq. 32:
            !
            !   [2] u_b = u_vav / (1 + beta * F2( z=s))
            !
            ! Substituting [2] into [1] yields:
            !
            !   [3] u( z) = u_vav * (1 + beta * F1( z)) / (1 + beta * F2( z=s))
            !
            ! This can be rearranged to read:
            !
            !   [4] u_vav * (1 + beta * F1( z)) / (1 + beta * F2( z=s)) - u( z) = 0

            ! u_vav term
            col_nh = tiuv2nh( ti,uv)
            val = (1._dp + self%DIVA%basal_friction_coefficient_b( ti) * self%DIVA%F1_3D_b( ti,k)) / &
                  (1._dp + self%DIVA%basal_friction_coefficient_b( ti) * self%DIVA%F2_3D_b( ti,1))
            call A_combi%add_entry( row_nh, col_nh, val)

            ! u( z) term
            val = -1._dp
            call A_combi%add_entry( row_nh, row_nh, val)

            ! The load vector is zero in this case
            b_combi( row_nh) = 0._dp

            ! Take the previous velocity solution as the initial guess
            uv_combi( row_nh) = self%BPA%u_bk( ti,k)

          end if ! if (C%choice_sliding_law == 'no_sliding') then

        else
          call crash('mask inconsistency; expected 3-D velocities, but both mask_BPA_b and mask_3D_from_DIVA_b are false!')
        end if

      else
        call crash('nh2tiuv_tikuv( row_nh,1) = {int_01}, should be only 1 or 2!', int_01 = nh2tiuv_tikuv( row_nh,1))
      end if ! if     (nh2tiuv_tikuv( row_nh,1) == 1) then

    end do ! do row_nh = i1, i2

    call A_combi%finalise()

    ! == Solve the matrix equation
    ! ============================

    ! use PETSc to solve the matrix equation
    call solve_matrix_equation_CSR_PETSc( A_combi, b_combi, uv_combi, &
      self%PETSc_rtol, self%PETSc_abstol, n_Axb_its, &
      PETSc_KSPtype = C%stress_balance_PETSc_KSPtype, PETSc_PCtype = C%stress_balance_PETSc_PCtype)

    ! Get velocities back from the combined vector
    do ti = self%mesh%ti1, self%mesh%ti2

      if (self%mask_DIVA_b( ti)) then
        ! The DIVA was solved here

        ! Get vertically averaged DIVA velocities back from the combined vector
        nhu = tiuv2nh( ti,1)
        nhv = tiuv2nh( ti,2)
        self%u_vav_b( ti) = uv_combi( nhu)
        self%v_vav_b( ti) = uv_combi( nhv)

        ! Calculate 3-D velocities from the vertically averaged DIVA velocities
        if (C%choice_sliding_law == 'no_sliding') then
          ! Lipscomb et al., 2019, Eq. 29, and text between Eqs. 33 and 34
          do k = 1, self%mesh%nz
            self%u_bk( ti,k) = self%u_vav_b( ti) * self%DIVA%beta_eff_b( ti) * self%DIVA%F1_3D_b( ti,k)
            self%v_bk( ti,k) = self%v_vav_b( ti) * self%DIVA%beta_eff_b( ti) * self%DIVA%F1_3D_b( ti,k)
          end DO
        else ! if (C%choice_sliding_law == 'no_sliding') then
          ! Lipscomb et al., 2019, Eq. 29
          do k = 1, self%mesh%nz
            self%u_bk( ti,k) = self%u_vav_b( ti) * &
              (1._dp + self%DIVA%basal_friction_coefficient_b( ti) * self%DIVA%F1_3D_b( ti,k)) / &
              (1._dp + self%DIVA%basal_friction_coefficient_b( ti) * self%DIVA%F2_3D_b( ti,1))
            self%v_bk( ti,k) = self%v_vav_b( ti) * &
              (1._dp + self%DIVA%basal_friction_coefficient_b( ti) * self%DIVA%F1_3D_b( ti,k)) / &
              (1._dp + self%DIVA%basal_friction_coefficient_b( ti) * self%DIVA%F2_3D_b( ti,1))
          end DO
        end if ! if (C%choice_sliding_law == 'no_sliding') then

      elseif (self%mask_BPA_b( ti)) then
        ! The BPA was solved here

        ! Get 3-D BPA velocities back from the combined vector
        do k = 1, self%mesh%nz
          nhu = tikuv2nh( ti,k,1)
          nhv = tikuv2nh( ti,k,2)
          self%u_bk( ti,k) = uv_combi( nhu)
          self%v_bk( ti,k) = uv_combi( nhv)
        end DO

        ! Calculate vertically averaged velocities from the 3-D BPA velocities
        self%u_vav_b( ti) = 0._dp
        self%v_vav_b( ti) = 0._dp

        do k = 1, self%mesh%nz
          if     (k == 1) then
            dzeta = self%mesh%zeta_stag( 1)
          elseif (k == self%mesh%nz) then
            dzeta = 1._dp - self%mesh%zeta_stag( self%mesh%nz-1)
          else
            dzeta = self%mesh%zeta_stag( k) - self%mesh%zeta_stag( k-1)
          end if
          self%u_vav_b( ti) = self%u_vav_b( ti) + dzeta * self%u_bk( ti,k)
          self%v_vav_b( ti) = self%v_vav_b( ti) + dzeta * self%v_bk( ti,k)
        end do ! do k = 1, mesh%nz

      else
        ! Safety
        call crash('neither the DIVA nor the BPA was apparently solved here!')
      end if ! if (hybrid%mask_DIVA_b( ti)) then
    end do ! do row_nh = i1, i2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_hybrid_DIVA_BPA_linearised

  subroutine calc_masked_DIVA_stiffness_matrix_and_load_vector( self, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, A_DIVA, b_DIVA)
    !< Calculate the stiffness matrix for the masked DIVA

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self
    integer,  dimension(self%mesh%ti1:self%mesh%ti2),    intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2),    intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2),    intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    type(type_CSR_matrix_dp),                            intent(  out) :: A_DIVA
    real(dp), dimension(:), allocatable,                 intent(  out) :: b_DIVA

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'calc_masked_DIVA_stiffness_matrix_and_load_vector'
    integer                       :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    integer                       :: row_tiuv,ti,uv
    character(len=:), allocatable :: choice_BC_u, choice_BC_v

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
    ! =============================================================================

    ! Matrix size
    ncols           = self%mesh%nTri     * 2      ! from
    ncols_loc       = self%mesh%nTri_loc * 2
    nrows           = self%mesh%nTri     * 2      ! to
    nrows_loc       = self%mesh%nTri_loc * 2
    nnz_est_proc    = self%mesh%M2_ddx_b_b%nnz * 4

    call A_DIVA%allocate( nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for the load vector
    allocate( b_DIVA( self%mesh%ti1*2-1: self%mesh%ti2*2))

    ! == Construct the stiffness matrix for the linearised DIVA
    ! ========================================================

    do row_tiuv = A_DIVA%i1, A_DIVA%i2

      ti = self%mesh%n2tiuv( row_tiuv,1)
      uv = self%mesh%n2tiuv( row_tiuv,2)

      if (BC_prescr_mask_b( ti) == 1) then
        ! Dirichlet boundary condition; velocities are prescribed for this triangle

        ! Stiffness matrix: diagonal element set to 1
        call A_DIVA%add_entry( row_tiuv, row_tiuv, 1._dp)

        ! Load vector: prescribed velocity
        if     (uv == 1) then
          b_DIVA( row_tiuv) = BC_prescr_u_b( ti)
        elseif (uv == 2) then
          b_DIVA( row_tiuv) = BC_prescr_v_b( ti)
        else
          call crash('uv can only be 1 or 2!')
        end if

      elseif (.not. self%mask_DIVA_b( ti)) then
        ! The BPA is solved here, not the DIVA

        call A_DIVA%add_empty_row( row_tiuv)

      elseif (self%mesh%TriBI( ti) > 0) then
        ! Domain border: apply boundary conditions

        select case (self%mesh%TriBI( ti))
        case default
          call crash('invalid TriBI value at triangle {int_01}', int_01 = ti)
        case (1,2)
          ! Northern domain border
          choice_BC_u = C%BC_u_north
          choice_BC_v = C%BC_v_north
        case (3,4)
          ! Eastern domain border
          choice_BC_u = C%BC_u_east
          choice_BC_v = C%BC_v_east
        case (5,6)
          ! Southern domain border
          choice_BC_u = C%BC_u_south
          choice_BC_v = C%BC_v_south
        case (7,8)
          ! Western domain border
          choice_BC_u = C%BC_u_west
          choice_BC_v = C%BC_v_west
        end select

        call self%DIVA%calc_SSA_DIVA_stiffness_matrix_row_BC( A_DIVA, b_DIVA, row_tiuv, choice_BC_u, choice_BC_v)

      else
        ! No boundary conditions apply; solve the DIVA

        if (C%do_include_SSADIVA_crossterms) then
          ! Calculate matrix coefficients for the full DIVA
          call self%DIVA%calc_SSA_DIVA_stiffness_matrix_row_free( self%DIVA%beta_eff_b, A_DIVA, b_DIVA, row_tiuv)
        else
          ! Calculate matrix coefficients for the DIVA sans the gradients of the effective viscosity (the "cross-terms")
          call self%DIVA%calc_SSA_DIVA_sans_stiffness_matrix_row_free( self%DIVA%beta_eff_b, A_DIVA, b_DIVA, row_tiuv)
        end if

      end if

    end do

    call A_DIVA%finalise()

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_masked_DIVA_stiffness_matrix_and_load_vector

  subroutine calc_masked_BPA_stiffness_matrix_and_load_vector( self, ice, &
    BC_prescr_mask_b, A_BPA, b_BPA)
    !< Calculate the stiffness matrix for the masked BPA

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self
    class(atype_ice_model_data),                         intent(in   ) :: ice
    integer,  dimension(self%mesh%ti1:self%mesh%ti2),    intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    type(type_CSR_matrix_dp),                            intent(  out) :: A_BPA
    real(dp), dimension(:), allocatable,                 intent(  out) :: b_BPA

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_masked_BPA_stiffness_matrix_and_load_vector'
    integer                     :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    integer                     :: row_tikuv,ti,k,uv

    ! Add routine to path
    call init_routine( routine_name)

    ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
    ! =============================================================================

    ! Matrix size
    ncols           = self%mesh%nTri     * self%mesh%nz * 2      ! from
    ncols_loc       = self%mesh%nTri_loc * self%mesh%nz * 2
    nrows           = self%mesh%nTri     * self%mesh%nz * 2      ! to
    nrows_loc       = self%mesh%nTri_loc * self%mesh%nz * 2
    nnz_est_proc    = self%mesh%M2_ddx_bk_bk%nnz * 4

    call A_BPA%allocate( nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! allocate memory for the load vector
    allocate( b_BPA( A_BPA%i1:A_BPA%i2))

    ! == Construct the stiffness matrix for the linearised BPA
    ! ========================================================

    do row_tikuv = A_BPA%i1, A_BPA%i2

      ti = self%mesh%n2tikuv( row_tikuv,1)
      k  = self%mesh%n2tikuv( row_tikuv,2)
      uv = self%mesh%n2tikuv( row_tikuv,3)

      if (BC_prescr_mask_b( ti) == 1 .or. .not. self%mask_BPA_b( ti)) then
        ! The DIVA is solved here, not the BPA

        call A_BPA%add_empty_row( row_tikuv)

      elseif (self%mesh%TriBI( ti) == 1 .or. self%mesh%TriBI( ti) == 2) then
        ! Northern domain border

        call self%BPA%calc_BPA_stiffness_matrix_row_BC_north( A_BPA, b_BPA, row_tikuv)

      elseif (self%mesh%TriBI( ti) == 3 .or. self%mesh%TriBI( ti) == 4) then
        ! Eastern domain border

        call self%BPA%calc_BPA_stiffness_matrix_row_BC_east( A_BPA, b_BPA, row_tikuv)

      elseif (self%mesh%TriBI( ti) == 5 .or. self%mesh%TriBI( ti) == 6) then
        ! Southern domain border

        call self%BPA%calc_BPA_stiffness_matrix_row_BC_south( A_BPA, b_BPA, row_tikuv)

      elseif (self%mesh%TriBI( ti) == 7 .or. self%mesh%TriBI( ti) == 8) then
        ! Western domain border

        call self%BPA%calc_BPA_stiffness_matrix_row_BC_west( A_BPA, b_BPA, row_tikuv)

      elseif (k == 1) then
        ! Ice surface

        call self%BPA%calc_BPA_stiffness_matrix_row_BC_surf( ice, A_BPA, b_BPA, row_tikuv)

      elseif (k == self%mesh%nz) then
        ! Ice base

        call self%BPA%calc_BPA_stiffness_matrix_row_BC_base( ice, A_BPA, b_BPA, row_tikuv)

      else
        ! No boundary conditions apply; solve the BPA

        call self%BPA%calc_BPA_stiffness_matrix_row_free( A_BPA, b_BPA, row_tikuv)

      end if

    end do

    call A_BPA%finalise()

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_masked_BPA_stiffness_matrix_and_load_vector

  subroutine calc_hybrid_solver_masks_transition( self, A_DIVA, A_BPA)
    !< Calculate the "transition" solver masks

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self
    type(type_CSR_matrix_dp),                            intent(in   ) :: A_DIVA
    type(type_CSR_matrix_dp),                            intent(in   ) :: A_BPA

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'calc_hybrid_solver_masks_transition'
    logical, dimension(self%mesh%nTri) :: mask_DIVA_halo_b, mask_BPA_halo_b
    integer                            :: row_tiuv,ti,kk1,kk2,kk,col_tjuv,tj
    integer                            :: row_tikuv,col_tjkuv
    integer                            :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Mark all triangles where the DIVA solver needs vertically averaged velocities
    mask_DIVA_halo_b = .false.
    do row_tiuv = A_DIVA%i1, A_DIVA%i2

      ti = self%mesh%n2tiuv( row_tiuv,1)

      kk1 = A_DIVA%ptr( row_tiuv)
      kk2 = A_DIVA%ptr( row_tiuv+1) - 1

      if (kk2 >= kk1) mask_DIVA_halo_b( ti) = .true.

      do kk = kk1, kk2

        col_tjuv = A_DIVA%ind( kk)
        tj = self%mesh%n2tiuv( col_tjuv,1)

        mask_DIVA_halo_b( tj) = .true.

      end do

    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, mask_DIVA_halo_b, self%mesh%nTri, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)

    ! Mark all triangles where the BPA solver needs vertically averaged velocities
    mask_BPA_halo_b = .false.
    do row_tikuv = A_BPA%i1, A_BPA%i2

      ti = self%mesh%n2tikuv( row_tikuv,1)

      kk1 = A_BPA%ptr( row_tikuv)
      kk2 = A_BPA%ptr( row_tikuv+1) - 1

      if (kk2 >= kk1) mask_BPA_halo_b( ti) = .true.

      do kk = kk1, kk2

        col_tjkuv = A_BPA%ind( kk)
        tj = self%mesh%n2tikuv( col_tjkuv,1)

        mask_BPA_halo_b( tj) = .true.

      end do

    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, mask_BPA_halo_b, self%mesh%nTri, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)

    ! Mark all triangles where the DIVA is solved, but a nearby BPA triangle needs 3-D velocities
    self%mask_3D_from_DIVA_b = .false.
    do ti = self%mesh%ti1, self%mesh%ti2
      if (self%mask_DIVA_b( ti) .and. mask_BPA_halo_b( ti)) then
        self%mask_3D_from_DIVA_b( ti) = .true.
      end if
    end do

    ! Mark all triangles where the BPA is solved, but a nearby DIVA triangle needs vertically averaged velocities
    self%mask_vav_from_BPA_b = .false.
    do ti = self%mesh%ti1, self%mesh%ti2
      if (self%mask_BPA_b( ti) .and. mask_DIVA_halo_b( ti)) then
        self%mask_vav_from_BPA_b( ti) = .true.
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_hybrid_solver_masks_transition

  subroutine calc_hybrid_solver_translation_tables( self, &
    tiuv2nh, tikuv2nh, nh2tiuv_tikuv, neq, i1, i2)
    !< Calculate combined DIVA/BPA translation tables

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self
    integer,  dimension(:,:  ), allocatable,             intent(  out) :: tiuv2nh
    integer,  dimension(:,:,:), allocatable,             intent(  out) :: tikuv2nh
    integer,  dimension(:,:  ), allocatable,             intent(  out) :: nh2tiuv_tikuv
    integer,                                             intent(  out) :: neq,i1,i2

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'calc_hybrid_solver_translation_tables'
    logical, dimension(self%mesh%nTri) :: mask_DIVA_b_tot
    logical, dimension(self%mesh%nTri) :: mask_BPA_b_tot
    logical, dimension(self%mesh%nTri) :: mask_3D_from_DIVA_b_tot
    logical, dimension(self%mesh%nTri) :: mask_vav_from_BPA_b_tot
    integer                            :: ti,k,uv

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather global masks
    call gather_to_all( self%mask_DIVA_b        , mask_DIVA_b_tot        )
    call gather_to_all( self%mask_BPA_b         , mask_BPA_b_tot         )
    call gather_to_all( self%mask_3D_from_DIVA_b, mask_3D_from_DIVA_b_tot)
    call gather_to_all( self%mask_vav_from_BPA_b, mask_vav_from_BPA_b_tot)

    ! allocate memory
    allocate( tiuv2nh      ( self%mesh%nTri               ,  2   ))
    allocate( tikuv2nh     ( self%mesh%nTri,  self%mesh%nz,  2   ))
    allocate( nh2tiuv_tikuv( self%mesh%nTri * self%mesh%nz * 2, 4))

    neq = 0
    i1  = 0
    i2  = 0

    do ti = 1, self%mesh%nTri

      if (ti == self%mesh%ti1) then
        i1 = neq + 1
      end if

      if (mask_DIVA_b_tot( ti) .or. mask_vav_from_BPA_b_tot( ti)) then
        ! Vertically averaged velocities must be defined here
        do uv = 1, 2
          neq = neq + 1
          tiuv2nh( ti,uv) = neq
          nh2tiuv_tikuv( neq,:) = [1,ti,0,uv]
        end do ! do uv = 1, 2
      end if

      if (mask_BPA_b_tot( ti) .or. mask_3D_from_DIVA_b_tot( ti)) then
        ! 3-D velocities must be defined here
        do k = 1, self%mesh%nz
        do uv = 1, 2
          neq = neq + 1
          tikuv2nh( ti,k,uv) = neq
          nh2tiuv_tikuv( neq,:) = [2,ti,k,uv]
        end DO
        end DO
      end if

      if (ti == self%mesh%ti2) then
        i2 = neq
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_hybrid_solver_translation_tables

! == Some useful tools for improving numerical stability of the viscosity iteration

  subroutine apply_velocity_limits( self)
    ! Limit velocities for improved stability

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'apply_velocity_limits'
    integer                     :: ti,k
    real(dp)                    :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    do ti = self%mesh%ti1, self%mesh%ti2
    do k  = 1, self%mesh%nz

      ! Calculate absolute speed
      uabs = SQRT( self%u_bk( ti,k)**2 + self%v_bk( ti,k)**2)

      ! Reduce velocities if neceBPAry
      if (uabs > C%vel_max) then
        self%u_bk( ti,k) = self%u_bk( ti,k) * C%vel_max / uabs
        self%v_bk( ti,k) = self%v_bk( ti,k) * C%vel_max / uabs
      end if

    end DO
    end DO

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_velocity_limits

  subroutine relax_viscosity_iterations( self, visc_it_relax)
    ! Reduce the change between velocity solutions

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self
    real(dp),                                            intent(in   ) :: visc_it_relax

    ! Local variables:
    character(len=*), parameter :: routine_name = 'relax_viscosity_iterations'
    integer                     :: ti

    ! Add routine to path
    call init_routine( routine_name)

    do ti = self%mesh%ti1, self%mesh%ti2
      self%u_bk( ti,:) = (visc_it_relax * self%u_bk( ti,:)) + ((1._dp - visc_it_relax) * self%u_bk_prev( ti,:))
      self%v_bk( ti,:) = (visc_it_relax * self%v_bk( ti,:)) + ((1._dp - visc_it_relax) * self%v_bk_prev( ti,:))
    end DO

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine relax_viscosity_iterations

  subroutine calc_visc_iter_UV_resid( self, resid_UV)
    ! Calculate the L2-norm of the two consecutive velocity solutions

    ! In/output variables:
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self
    real(dp),                                            intent(  out) :: resid_UV

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_visc_iter_UV_resid'
    integer                     :: ierr
    integer                     :: ti,k
    real(dp)                    :: res1, res2

    ! Add routine to path
    call init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    do ti = self%mesh%ti1, self%mesh%ti2
    do k  = 1, self%mesh%nz

      res1 = res1 + (self%u_bk( ti,k) - self%u_bk_prev( ti,k))**2
      res1 = res1 + (self%v_bk( ti,k) - self%v_bk_prev( ti,k))**2

      res2 = res2 + (self%u_bk( ti,k) + self%u_bk_prev( ti,k))**2
      res2 = res2 + (self%v_bk( ti,k) + self%v_bk_prev( ti,k))**2

    end DO
    end DO

    ! Combine results from all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate residual
    resid_UV = 2._dp * res1 / MAX( res2, 1E-8_dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_visc_iter_UV_resid

! == Initialisation

  subroutine create_restart_file_old_hybrid_DIVA_BPA( self)
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(inout) :: self
  end subroutine create_restart_file_old_hybrid_DIVA_BPA

  subroutine write_to_restart_file_old_hybrid_DIVA_BPA( self, time)
    class(type_momentum_balance_solver_hybrid_DIVA_BPA), intent(in   ) :: self
    real(dp),                                            intent(in   ) :: time
  end subroutine write_to_restart_file_old_hybrid_DIVA_BPA

end module momentum_balance_solver_hybrid_DIVA_BPA
