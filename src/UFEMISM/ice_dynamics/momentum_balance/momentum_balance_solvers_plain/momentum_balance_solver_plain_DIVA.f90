module momentum_balance_solver_plain_DIVA

  ! Routines for calculating ice velocities using the Depth-Integrated Viscosity Approximation (DIVA)

  use mpi_basic, only: par
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use netcdf_io_main
  use mesh_disc_apply_operators, only: map_a_b_2D, map_a_b_3D, map_b_a_2D, map_b_a_3D
  use reallocate_mod, only: reallocate_bounds, reallocate_clean
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D, &
    map_from_mesh_to_mesh_with_reallocation_3D
  use bed_roughness_model_types, only: type_bed_roughness_model
  use mpi_distributed_memory, only: gather_to_all
  use momentum_balance_solver_plain_SSADIVA, only: atype_momentum_balance_solver_plain_SSADIVA
  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, &
    MPI_LOR, MPI_LOGICAL, MPI_MIN, MPI_MAX
  use checksum_mod, only: checksum
  use sliding_laws, only: calc_basal_friction_coefficient
  use mesh_zeta, only: integrate_from_zeta_is_one_to_zeta_is_zetap, vertical_average
  use constitutive_equation, only: calc_ice_rheology_Glen, calc_effective_viscosity_Glen_3D_uv_only
  use mesh_disc_apply_operators, only: map_a_b_2D, map_a_b_3D, ddx_a_b_2D, ddy_a_b_2D, &
    map_b_a_2D, map_b_a_3D, ddx_b_a_2D, ddy_b_a_2D

  implicit none

  private

  public :: type_momentum_balance_solver_plain_DIVA

  type, extends(atype_momentum_balance_solver_plain_SSADIVA) :: type_momentum_balance_solver_plain_DIVA

    ! Solution
    real(dp), dimension(:  ), allocatable :: u_base_b                    ! [m yr^-1] 2-D horizontal ice velocity at the ice base
    real(dp), dimension(:  ), allocatable :: v_base_b
    real(dp), dimension(:,:), allocatable :: u_3D_b                      ! [m yr^-1] 3-D horizontal ice velocity
    real(dp), dimension(:,:), allocatable :: v_3D_b

    ! Intermediate data fields
    real(dp), dimension(:,:), allocatable :: du_dz_3D_a                  ! [yr^-1] 3-D vertical shear strain rates
    real(dp), dimension(:,:), allocatable :: dv_dz_3D_a
    real(dp), dimension(:,:), allocatable :: eta_3D_a                    ! Effective viscosity
    real(dp), dimension(:,:), allocatable :: eta_3D_b
    real(dp), dimension(:,:), allocatable :: F1_3D_a                     ! F-integrals
    real(dp), dimension(:,:), allocatable :: F2_3D_a
    real(dp), dimension(:,:), allocatable :: F1_3D_b
    real(dp), dimension(:,:), allocatable :: F2_3D_b
    real(dp), dimension(:  ), allocatable :: beta_eff_a                  ! "Effective" friction coefficient (turning the SSA into the DIVA)
    real(dp), dimension(:  ), allocatable :: beta_eff_b
    real(dp), dimension(:  ), allocatable :: tau_bx_b                    ! Basal shear stress
    real(dp), dimension(:  ), allocatable :: tau_by_b

    contains

      ! Procedures for model memory management and operation
      procedure, public :: get_momentum_balance_solver_plain_name
      procedure, public :: allocate_momentum_balance_solver_plain   => momentum_balance_solver_plain_DIVA_allocate
      procedure, public :: deallocate_momentum_balance_solver_plain => momentum_balance_solver_plain_DIVA_deallocate
      procedure, public :: initialise_momentum_balance_solver_plain => momentum_balance_solver_plain_DIVA_initialise
      procedure, public :: run_momentum_balance_solver_plain        => momentum_balance_solver_plain_DIVA_run
      procedure, public :: remap_momentum_balance_solver_plain      => momentum_balance_solver_plain_DIVA_remap

      procedure, public :: create_restart_file_DIVA
      procedure, public :: write_to_restart_file_DIVA

      procedure, private :: initialise_DIVA_velocities_from_file
      procedure, private :: calc_vertical_shear_strain_rates
      procedure, private :: calc_effective_viscosity
      procedure, private :: calc_F_integrals
      procedure, private :: calc_effective_basal_friction_coefficient
      procedure, private :: calc_basal_velocities
      procedure, private :: calc_basal_shear_stress
      procedure, private :: calc_3D_velocities

  end type type_momentum_balance_solver_plain_DIVA

contains

  ! == Main routines

  subroutine momentum_balance_solver_plain_DIVA_allocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_plain_DIVA_allocate'

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate variables that are shared between the SSA and DIVA solvers
    call self%allocate_shared_SSA_DIVA_variables()

    ! Allocate variables that are specific to the DIVA solver

    ! Solution
    allocate( self%u_base_b(                     self%mesh%ti1:self%mesh%ti2               ), source = 0._dp)
    allocate( self%v_base_b(                     self%mesh%ti1:self%mesh%ti2               ), source = 0._dp)
    allocate( self%u_3D_b(                       self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), source = 0._dp)
    allocate( self%v_3D_b(                       self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), source = 0._dp)

    ! Intermediate data fields
    allocate( self%du_dz_3D_a(                   self%mesh%vi1:self%mesh%vi2,1:self%mesh%nz), source = 0._dp)
    allocate( self%dv_dz_3D_a(                   self%mesh%vi1:self%mesh%vi2,1:self%mesh%nz), source = 0._dp)
    allocate( self%eta_3D_a(                     self%mesh%vi1:self%mesh%vi2,1:self%mesh%nz), source = 0._dp)
    allocate( self%eta_3D_b(                     self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), source = 0._dp)
    allocate( self%F1_3D_a(                      self%mesh%vi1:self%mesh%vi2,1:self%mesh%nz), source = 0._dp)
    allocate( self%F2_3D_a(                      self%mesh%vi1:self%mesh%vi2,1:self%mesh%nz), source = 0._dp)
    allocate( self%F1_3D_b(                      self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), source = 0._dp)
    allocate( self%F2_3D_b(                      self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz), source = 0._dp)
    allocate( self%beta_eff_a(                   self%mesh%vi1:self%mesh%vi2               ), source = 0._dp)
    allocate( self%beta_eff_b(                   self%mesh%ti1:self%mesh%ti2               ), source = 0._dp)
    allocate( self%tau_bx_b(                     self%mesh%ti1:self%mesh%ti2               ), source = 0._dp)
    allocate( self%tau_by_b(                     self%mesh%ti1:self%mesh%ti2               ), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_plain_DIVA_allocate

  subroutine momentum_balance_solver_plain_DIVA_deallocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_plain_DIVA_deallocate'

    ! Add routine to path
    call init_routine( routine_name)

    ! Deallocate variables that are shared between the SSA and DIVA solvers
    call self%deallocate_shared_SSA_DIVA_variables()

    ! Deallocate variables that are specific to the DIVA solver

    ! Solution
    deallocate( self%u_base_b)
    deallocate( self%v_base_b)
    deallocate( self%u_3D_b)
    deallocate( self%v_3D_b)

    ! Intermediate data fields
    deallocate( self%du_dz_3D_a)
    deallocate( self%dv_dz_3D_a)
    deallocate( self%eta_3D_a)
    deallocate( self%eta_3D_b)
    deallocate( self%F1_3D_a)
    deallocate( self%F2_3D_a)
    deallocate( self%F1_3D_b)
    deallocate( self%F2_3D_b)
    deallocate( self%beta_eff_a)
    deallocate( self%beta_eff_b)
    deallocate( self%tau_bx_b)
    deallocate( self%tau_by_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_plain_DIVA_deallocate

  subroutine momentum_balance_solver_plain_DIVA_initialise( self)

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'momentum_balance_solver_plain_DIVA_initialise'
    character(len=:), allocatable :: choice_initial_velocity

    ! Add routine to path
    call init_routine( routine_name)

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
      self%u_vav_b( self%mesh%ti1:self%mesh%ti2) = 0._dp
      self%v_vav_b( self%mesh%ti1:self%mesh%ti2) = 0._dp
    case ('read_from_file')
      call self%initialise_DIVA_velocities_from_file()
    end select

    ! Set tolerances for PETSc matrix solver for the linearised DIVA
    self%PETSc_rtol   = C%stress_balance_PETSc_rtol
    self%PETSc_abstol = C%stress_balance_PETSc_abstol

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_plain_DIVA_initialise

  subroutine initialise_DIVA_velocities_from_file( self)
    !< Initialise the velocities for the DIVA solver from an external NetCDF file

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'initialise_DIVA_velocities_from_file'
    character(len=:), allocatable                    :: filename
    real(dp)                                         :: timeframe
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2) :: u_vav_b_prev_loc
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2) :: v_vav_b_prev_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine the filename and timeframe to read for this model region
    select case (self%region_name())
    case default
      call crash('unknown model region "' // self%region_name() // '"!')
    case ('NAM')
      filename  = C%filename_initial_velocity_NAM
      timeframe = C%timeframe_initial_velocity_NAM
    case ('EAS')
      filename  = C%filename_initial_velocity_EAS
      timeframe = C%timeframe_initial_velocity_EAS
    case ('GRL')
      filename  = C%filename_initial_velocity_GRL
      timeframe = C%timeframe_initial_velocity_GRL
    case ('ANT')
      filename  = C%filename_initial_velocity_ANT
      timeframe = C%timeframe_initial_velocity_ANT
    end select

    ! Exception for when we want to flexible read the last output file of a previous UFEMISM simulation
    if (index( filename,'_LAST.nc') > 1) then
      call find_last_output_file( filename)
      call find_last_timeframe(   filename, timeframe)
    end if

    ! Write to terminal
    if (par%primary) write(0,*) '   Initialising DIVA velocities from file "' // &
      UPSY%stru%colour_string( trim( filename),'light blue') // '"...'

    ! Solution
    call read_field_from_mesh_file_dp_2D_b( filename, 'u_vav_b'                     , self%u_vav_b                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'v_vav_b'                     , self%v_vav_b                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'u_base_b'                    , self%u_base_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'v_base_b'                    , self%v_base_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D_b( filename, 'u_3D_b'                      , self%u_3D_b                      , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D_b( filename, 'v_3D_b'                      , self%v_3D_b                      , time_to_read = timeframe)

    ! Intermediate data fields
    call read_field_from_mesh_file_dp_2D  ( filename, 'du_dx_a'                     , self%du_dx_a                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D  ( filename, 'du_dy_a'                     , self%du_dy_a                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D  ( filename, 'dv_dx_a'                     , self%dv_dx_a                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D  ( filename, 'dv_dy_a'                     , self%dv_dy_a                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D  ( filename, 'du_dz_3D_a'                  , self%du_dz_3D_a                  , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D  ( filename, 'dv_dz_3D_a'                  , self%dv_dz_3D_a                  , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D  ( filename, 'eta_3D_a'                    , self%eta_3D_a                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D_b( filename, 'eta_3D_b'                    , self%eta_3D_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D  ( filename, 'eta_vav_a'                   , self%eta_vav_a                   , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D  ( filename, 'N_a'                         , self%N_a                         , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'N_b'                         , self%N_b                         , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'dN_dx_b'                     , self%dN_dx_b                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'dN_dy_b'                     , self%dN_dy_b                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D  ( filename, 'F1_3D_a'                     , self%F1_3D_a                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D  ( filename, 'F2_3D_a'                     , self%F2_3D_a                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D_b( filename, 'F1_3D_b'                     , self%F1_3D_b                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D_b( filename, 'F2_3D_b'                     , self%F2_3D_b                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'basal_friction_coefficient_b', self%basal_friction_coefficient_b, time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D  ( filename, 'beta_eff_a'                  , self%beta_eff_a                  , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'beta_eff_b'                  , self%beta_eff_b                  , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'tau_bx_b'                    , self%tau_bx_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'tau_by_b'                    , self%tau_by_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'tau_dx_b'                    , self%tau_dx_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'tau_dy_b'                    , self%tau_dy_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'u_vav_b_prev'                , u_vav_b_prev_loc                 , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'v_vav_b_prev'                , v_vav_b_prev_loc                 , time_to_read = timeframe)

    call gather_to_all( u_vav_b_prev_loc, self%u_vav_b_prev)
    call gather_to_all( v_vav_b_prev_loc, self%v_vav_b_prev)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_DIVA_velocities_from_file

  subroutine momentum_balance_solver_plain_DIVA_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate ice velocities by solving the Depth-Integrated Viscosity Approximation

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self
    class(atype_ice_model_data),                    intent(inout) :: ice
    class(atype_ice_geometry_model_data),           intent(in   ) :: geom
    type(type_bed_roughness_model),                 intent(in   ) :: bed_roughness
    integer,  dimension(:  ), optional,             intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:  ), optional,             intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(:  ), optional,             intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    integer,  dimension(:,:), optional,             intent(in   ) :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:,:), optional,             intent(in   ) :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
    real(dp), dimension(:,:), optional,             intent(in   ) :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'momentum_balance_solver_plain_DIVA_run'
    logical                             :: grounded_ice_exists
    integer                             :: ierr
    integer,  dimension(:), allocatable :: BC_prescr_mask_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_u_b_applied
    real(dp), dimension(:), allocatable :: BC_prescr_v_b_applied
    integer                             :: viscosity_iteration_i
    logical                             :: has_converged
    real(dp)                            :: L2_uv, L2_uv_prev
    real(dp)                            :: uv_min, uv_max
    real(dp)                            :: visc_it_relax_applied
    real(dp)                            :: Glens_flow_law_epsilon_sq_0_applied
    integer                             :: nit_diverg_consec
    integer                             :: n_Axb_its_visc_it

    ! Add routine to path
    call init_routine( routine_name)

    ! if there is no grounded ice, no need (in fact, no way) to solve the DIVA
    grounded_ice_exists = any( geom%mask_grounded_ice)
    call MPI_ALLREDUCE( MPI_IN_PLACE, grounded_ice_exists, 1, MPI_logical, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (.not. grounded_ice_exists) then
      self%u_vav_b ( self%mesh%ti1:self%mesh%ti2  ) = 0._dp
      self%v_vav_b ( self%mesh%ti1:self%mesh%ti2  ) = 0._dp
      self%u_base_b( self%mesh%ti1:self%mesh%ti2  ) = 0._dp
      self%v_base_b( self%mesh%ti1:self%mesh%ti2  ) = 0._dp
      self%u_3D_b  ( self%mesh%ti1:self%mesh%ti2,:) = 0._dp
      self%v_3D_b  ( self%mesh%ti1:self%mesh%ti2,:) = 0._dp
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
      BC_prescr_mask_b_applied = BC_prescr_mask_b
      BC_prescr_u_b_applied    = BC_prescr_u_b
      BC_prescr_v_b_applied    = BC_prescr_v_b
    else
      BC_prescr_mask_b_applied = 0
      BC_prescr_u_b_applied    = 0._dp
      BC_prescr_v_b_applied    = 0._dp
    end if

    ! Calculate the driving stress
    call self%calc_driving_stress( geom)

    ! Adaptive relaxation parameter for the viscosity iteration
    L2_uv                               = 1E9_dp
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

      ! Calculate the horizontal strain rates for the current velocity solution
      call self%calc_horizontal_strain_rates()

      ! Calculate the vertical shear strain rates
      call self%calc_vertical_shear_strain_rates()

      ! Calculate the effective viscosity for the current velocity solution
      call self%calc_effective_viscosity( ice, geom, Glens_flow_law_epsilon_sq_0_applied)

      ! Calculate the F-integrals (Lipscomb et al. (2019), Eq. 30)
      call self%calc_F_integrals( geom)

      ! Calculate the "effective" friction coefficient (turning the SSA into the DIVA)
      call self%calc_effective_basal_friction_coefficient( ice, geom, bed_roughness)

      ! Solve the linearised DIVA to calculate a new velocity solution
      call self%solve_SSA_DIVA_linearised( n_Axb_its_visc_it, &
        BC_prescr_mask_b_applied, BC_prescr_u_b_applied, BC_prescr_v_b_applied)

      ! Update stability info
      self%n_Axb_its = self%n_Axb_its + n_Axb_its_visc_it

      ! Limit velocities for improved stability
      call self%apply_velocity_limits()

      ! Reduce the change between velocity solutions
      call self%relax_viscosity_iterations( visc_it_relax_applied)

      ! Calculate basal velocities
      call self%calc_basal_velocities ()

      ! Calculate basal shear stress
      call self%calc_basal_shear_stress()

      ! Calculate the L2-norm of the two consecutive velocity solutions
      L2_uv_prev = L2_uv
      call self%calc_L2_norm_uv( L2_uv)

      ! if the viscosity iteration diverges, lower the relaxation parameter
      if (L2_uv > L2_uv_prev) then
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
      uv_min = minval( self%u_vav_b)
      uv_max = maxval( self%u_vav_b)
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_min, 1, MPI_doUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
      call MPI_ALLREDUCE( MPI_IN_PLACE, uv_max, 1, MPI_doUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)
      ! if (par%primary) WRITE(0,*) '    DIVA - viscosity iteration ', viscosity_iteration_i, ', u = [', uv_min, ' - ', uv_max, '], L2_uv = ', L2_uv

      ! if the viscosity iteration has converged, or has reached the maximum allowed number of iterations, stop it.
      has_converged = .false.
      if (L2_uv < C%visc_it_norm_dUV_tol) then
        has_converged = .true.
      end if

      ! if we've reached the maximum allowed number of iterations without converging, throw a warning
      if (viscosity_iteration_i > C%visc_it_nit) then
        if (par%primary) call warning('viscosity iteration failed to converge within {int_01} iterations!', int_01 = C%visc_it_nit)
        exit viscosity_iteration
      end if

    end do viscosity_iteration

    ! Calculate 3-D ice velocities
    call self%calc_3D_velocities()

    ! Stability info
    self%n_visc_its = viscosity_iteration_i

    call checksum( self%mesh%pai_Tri, self%u_vav_b, 'DIVA%u_vav_b')
    call checksum( self%mesh%pai_Tri, self%v_vav_b, 'DIVA%v_vav_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_plain_DIVA_run

  subroutine momentum_balance_solver_plain_DIVA_remap( self, mesh_old, mesh_new)

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self
    type(type_mesh),                                intent(in   ) :: mesh_old
    type(type_mesh), target,                        intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'momentum_balance_solver_plain_DIVA_remap'
    real(dp), dimension(:  ), allocatable :: tau_bx_a
    real(dp), dimension(:  ), allocatable :: tau_by_a
    real(dp), dimension(:,:), allocatable :: eta_3D_a
    real(dp), dimension(:,:), allocatable :: u_3D_a
    real(dp), dimension(:,:), allocatable :: v_3D_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap variables that are shared between the SSA and DIVA solvers
    call self%remap_shared_SSA_DIVA_variables( mesh_old, mesh_new)

    ! Remap variables that are specific to the DIVA solver

    ! Remap the fields that are re-used during the viscosity iteration
    ! ================================================================

    ! allocate memory for velocities on the a-grid (vertices)
    allocate( tau_bx_a( mesh_old%vi1: mesh_old%vi2             ))
    allocate( tau_by_a( mesh_old%vi1: mesh_old%vi2             ))
    allocate( eta_3D_a( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    allocate( u_3D_a  ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    allocate( v_3D_a  ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))

    ! Map data from the triangles of the old mesh to the vertices of the old mesh
    call map_b_a_2D( mesh_old, self%tau_bx_b, tau_bx_a)
    call map_b_a_2D( mesh_old, self%tau_by_b, tau_by_a)
    call map_b_a_3D( mesh_old, self%eta_3D_b, eta_3D_a)
    call map_b_a_3D( mesh_old, self%u_3D_b  , u_3D_a  )
    call map_b_a_3D( mesh_old, self%v_3D_b  , v_3D_a  )

    ! Remap data from the vertices of the old mesh to the vertices of the new mesh
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, tau_bx_a, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, tau_by_a, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, eta_3D_a, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, u_3D_a  , '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, v_3D_a  , '2nd_order_conservative')

    ! reallocate memory for the data on the triangles
    call reallocate_bounds( self%tau_bx_b , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( self%tau_by_b , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( self%eta_3D_b , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( self%u_3D_b   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( self%v_3D_b   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)

    ! Map data from the vertices of the new mesh to the triangles of the new mesh
    call map_a_b_2D( mesh_new, tau_bx_a, self%tau_bx_b)
    call map_a_b_2D( mesh_new, tau_by_a, self%tau_by_b)
    call map_a_b_3D( mesh_new, eta_3D_a, self%eta_3D_b)
    call map_a_b_3D( mesh_new, u_3D_a  , self%u_3D_b  )
    call map_a_b_3D( mesh_new, v_3D_a  , self%v_3D_b  )

    ! reallocate everything else
    ! ==========================

    call reallocate_bounds( self%u_base_b  , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( self%v_base_b  , mesh_new%ti1, mesh_new%ti2             )
   !call reallocate_bounds( DIVA%u_3D_b    , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
   !call reallocate_bounds( DIVA%v_3D_b    , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( self%du_dz_3D_a, mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    call reallocate_bounds( self%dv_dz_3D_a, mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    call reallocate_bounds( self%eta_3D_a  , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
   !call reallocate_bounds( DIVA%eta_3D_b  , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( self%F1_3D_a   , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    call reallocate_bounds( self%F2_3D_a   , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    call reallocate_bounds( self%F1_3D_b   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( self%F2_3D_b   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( self%beta_eff_a, mesh_new%vi1, mesh_new%vi2             )
    call reallocate_bounds( self%beta_eff_b, mesh_new%ti1, mesh_new%ti2             )
   !call reallocate_bounds( DIVA%tau_bx_b  , mesh_new%ti1, mesh_new%ti2             )
   !call reallocate_bounds( DIVA%tau_by_b  , mesh_new%ti1, mesh_new%ti2             )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_plain_DIVA_remap

  function get_momentum_balance_solver_plain_name( self) result( model_name)
    class(type_momentum_balance_solver_plain_DIVA), intent(in) :: self
    character(len=:), allocatable :: model_name
    model_name = 'DIVA'
  end function get_momentum_balance_solver_plain_name

  ! == Calculate several intermediate terms in the SSA

  subroutine calc_vertical_shear_strain_rates( self)
    ! Calculate the vertical shear strain rates

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter                                     :: routine_name = 'calc_vertical_shear_strain_rates'
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz) :: du_dz_3D_b
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2,1:self%mesh%nz) :: dv_dz_3D_b
    integer                                                         :: ti,k

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate (parameterised) vertical shear strain rates on the b-grid (Lipscomb et al., 2019, Eq. 36)
    do ti = self%mesh%ti1, self%mesh%ti2
    do k = 1, self%mesh%nz
      du_dz_3D_b( ti,k) = self%tau_bx_b( ti) * self%mesh%zeta( k) / max( C%visc_eff_min, self%eta_3D_b( ti,k))
      dv_dz_3D_b( ti,k) = self%tau_by_b( ti) * self%mesh%zeta( k) / max( C%visc_eff_min, self%eta_3D_b( ti,k))
    end do
    end do

    ! Map vertical shear strain rates from the b-grid to the a-grid
    call map_b_a_3D( self%mesh, du_dz_3D_b, self%du_dz_3D_a)
    call map_b_a_3D( self%mesh, dv_dz_3D_b, self%dv_dz_3D_a)

    call checksum( self%mesh%pai_V, self%du_dz_3D_a, 'DIVA%du_dz_3D_a')
    call checksum( self%mesh%pai_V, self%dv_dz_3D_a, 'DIVA%dv_dz_3D_a')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_vertical_shear_strain_rates

  subroutine calc_effective_viscosity( self, ice, geom, Glens_flow_law_epsilon_sq_0_applied)

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self
    class(atype_ice_model_data),                    intent(inout) :: ice
    class(atype_ice_geometry_model_data),           intent(in   ) :: geom
    real(dp),                                       intent(in   ) :: Glens_flow_law_epsilon_sq_0_applied

    ! Local variables:
    character(len=*), parameter       :: routine_name = 'calc_effective_viscosity'
    integer                           :: vi,k
    real(dp)                          :: A_min, eta_max
    real(dp), dimension(self%mesh%nz) :: prof

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate maximum allowed effective viscosity, for stability
    A_min = 1E-18_dp
    eta_max = 0.5_dp * A_min**(-1._dp / C%Glens_flow_law_exponent) * &
      (Glens_flow_law_epsilon_sq_0_applied)**((1._dp - C%Glens_flow_law_exponent)/(2._dp*C%Glens_flow_law_exponent))

    ! Calculate the effective viscosity eta
    if (C%choice_flow_law == 'Glen') then
      ! Calculate the effective viscosity eta according to Glen's flow law

      ! Calculate flow factors
      call calc_ice_rheology_Glen( self%mesh, ice, geom)

      ! Calculate effective viscosity
      do vi = self%mesh%vi1, self%mesh%vi2
      do k  = 1, self%mesh%nz
        self%eta_3D_a( vi,k) = calc_effective_viscosity_Glen_3D_uv_only( &
          Glens_flow_law_epsilon_sq_0_applied, &
          self%du_dx_a( vi), self%du_dy_a( vi), self%du_dz_3D_a( vi,k), &
          self%dv_dx_a( vi), self%dv_dy_a( vi), self%dv_dz_3D_a( vi,k), ice%A_flow( vi,k))
      end do
      end do

    else
      call crash('unknown choice_flow_law "' // TRIM( C%choice_flow_law) // '"!')
    end if

    ! Safety
    self%eta_3D_a = min( max( self%eta_3D_a, C%visc_eff_min), eta_max)

    ! Map effective viscosity to the b-grid
    call map_a_b_3D( self%mesh, self%eta_3D_a, self%eta_3D_b)

    ! Calculate vertically averaged effective viscosity on the a-grid
    do vi = self%mesh%vi1, self%mesh%vi2
      prof = self%eta_3D_a( vi,:)
      self%eta_vav_a( vi) = vertical_average( self%mesh%zeta, prof)
    end do

    ! Calculate the product term N = eta * H on the a-grid
    do vi = self%mesh%vi1, self%mesh%vi2
      self%N_a( vi) = self%eta_vav_a( vi) * max( 0.1, geom%Hi( vi))
    end do

    ! Calculate the product term N and its gradients on the b-grid
    call map_a_b_2D( self%mesh, self%N_a, self%N_b    )
    call ddx_a_b_2D( self%mesh, self%N_a, self%dN_dx_b)
    call ddy_a_b_2D( self%mesh, self%N_a, self%dN_dy_b)

    call checksum( self%mesh%pai_V  , self%eta_3D_a , 'DIVA%eta_3D_a')
    call checksum( self%mesh%pai_Tri, self%eta_3D_b , 'DIVA%eta_3D_b')
    call checksum( self%mesh%pai_V  , self%eta_vav_a, 'DIVA%eta_vav_a')
    call checksum( self%mesh%pai_V  , self%N_a      , 'DIVA%N_a')
    call checksum( self%mesh%pai_Tri, self%N_b      , 'DIVA%N_b')
    call checksum( self%mesh%pai_Tri, self%dN_dx_b  , 'DIVA%dN_dx_b')
    call checksum( self%mesh%pai_Tri, self%dN_dy_b  , 'DIVA%dN_dy_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_viscosity

  subroutine calc_F_integrals( self, geom)
    !< Calculate the F-integrals on the a-grid (Lipscomb et al. (2019), Eq. 30)

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self
    class(atype_ice_geometry_model_data),           intent(in   ) :: geom

    ! Local variables:
    character(len=*), parameter       :: routine_name = 'calc_F_integrals'
    integer                           :: vi,k
    real(dp), dimension(self%mesh%nz) :: prof

    ! Add routine to path
    call init_routine( routine_name)

    do vi = self%mesh%vi1, self%mesh%vi2

      ! F1
      do k = 1, self%mesh%nz
        prof( k) = (self%mesh%zeta( k)    / self%eta_3D_a( vi,k))
      end do
      self%F1_3D_a( vi,:) = -max( 0.1_dp, geom%Hi( vi)) * integrate_from_zeta_is_one_to_zeta_is_zetap( self%mesh%zeta, prof)

      ! F2
      do k = 1, self%mesh%nz
        prof( k) = (self%mesh%zeta( k)**2 / self%eta_3D_a( vi,k))
      end do
      self%F2_3D_a( vi,:) = -max( 0.1_dp, geom%Hi( vi)) * integrate_from_zeta_is_one_to_zeta_is_zetap( self%mesh%zeta, prof)

    end do

    ! Map F-integrals from the a-grid to the b-grid
    call map_a_b_3D( self%mesh, self%F1_3D_a, self%F1_3D_b)
    call map_a_b_3D( self%mesh, self%F2_3D_a, self%F2_3D_b)

    call checksum( self%mesh%pai_Tri, self%F1_3D_b, 'DIVA%F1_3D_b')
    call checksum( self%mesh%pai_Tri, self%F2_3D_b, 'DIVA%F2_3D_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_F_integrals

  subroutine calc_effective_basal_friction_coefficient( self, ice, geom, bed_roughness)
    !< Calculate the "effective" friction coefficient (turning the SSA into the DIVA)

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self
    class(atype_ice_model_data),                    intent(inout) :: ice
    class(atype_ice_geometry_model_data),           intent(in   ) :: geom
    type(type_bed_roughness_model),                 intent(in   ) :: bed_roughness

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'calc_effective_basal_friction_coefficient'
    integer                                          :: vi,ti
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: u_base_a
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: v_base_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the basal friction coefficient beta_b for the current velocity solution
    ! This is where the sliding law is called!
    call map_b_a_2D( self%mesh, self%u_base_b, u_base_a)
    call map_b_a_2D( self%mesh, self%v_base_b, v_base_a)
    call calc_basal_friction_coefficient( self%mesh, geom, bed_roughness, u_base_a, v_base_a, &
      ice%effective_pressure, ice%till_yield_stress, ice%basal_friction_coefficient)

    ! Calculate beta_eff on the a-grid
    if (C%choice_sliding_law == 'no_sliding') then
      ! Exception for the case of no sliding (Lipscomb et al., 2019, Eq. 35)

      do vi = self%mesh%vi1, self%mesh%vi2
        self%beta_eff_a( vi) = 1._dp / self%F2_3D_a( vi,1)
      end do

    else ! if (C%choice_sliding_law == 'no_sliding') then
      ! Lipscomb et al., 2019, Eq. 33

      do vi = self%mesh%vi1, self%mesh%vi2
        self%beta_eff_a( vi) = ice%basal_friction_coefficient( vi) / (1._dp + ice%basal_friction_coefficient( vi) * self%F2_3D_a( vi,1))
      end do

    end if ! if (C%choice_sliding_law == 'no_sliding') then

    ! Map basal friction coefficient beta_b and effective basal friction coefficient beta_eff to the b-grid
    call map_a_b_2D( self%mesh, ice%basal_friction_coefficient, self%basal_friction_coefficient_b)
    call map_a_b_2D( self%mesh, self%beta_eff_a               , self%beta_eff_b                  )

    ! Apply the sub-grid grounded fraction, and limit the friction coefficient to improve stability
    if (C%do_GL_subgrid_friction) then
      ! On the b-grid
      do ti = self%mesh%ti1, self%mesh%ti2
        self%beta_eff_b( ti) = self%beta_eff_b( ti) * geom%fraction_gr_b( ti)**C%subgrid_friction_exponent_on_B_grid
      end do
    end if

    call checksum( self%mesh%pai_V  , ice%basal_friction_coefficient   , 'ice%basal_friction_coefficient')
    call checksum( self%mesh%pai_Tri, self%basal_friction_coefficient_b, 'DIVA%basal_friction_coefficient_b')
    call checksum( self%mesh%pai_V  , self%beta_eff_a                  , 'DIVA%beta_eff_a')
    call checksum( self%mesh%pai_Tri, self%beta_eff_b                  , 'DIVA%beta_eff_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_basal_friction_coefficient

  subroutine calc_basal_velocities( self)
    !< Calculate basal velocities (Lipscomb et al., 2019, Eq. 32)

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_basal_shear_stress'
    integer                     :: ti

    ! Add routine to path
    call init_routine( routine_name)

    if (C%choice_sliding_law == 'no_sliding') then
      ! Exception for the case of no sliding

      self%u_base_b( self%mesh%ti1:self%mesh%ti2) = 0._dp
      self%v_base_b( self%mesh%ti1:self%mesh%ti2) = 0._dp

    else

      ! Calculate basal velocities (Lipscomb et al., 2019, Eq. 32)
      do ti = self%mesh%ti1, self%mesh%ti2
        self%u_base_b( ti) = self%u_vav_b( ti) / (1._dp + self%basal_friction_coefficient_b( ti) * self%F2_3D_b( ti,1))
        self%v_base_b( ti) = self%v_vav_b( ti) / (1._dp + self%basal_friction_coefficient_b( ti) * self%F2_3D_b( ti,1))
      end do

    end if

    call checksum( self%mesh%pai_Tri, self%u_base_b, 'DIVA%u_base_b')
    call checksum( self%mesh%pai_Tri, self%v_base_b, 'DIVA%v_base_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_basal_velocities

  subroutine calc_basal_shear_stress( self)
    !< Calculate the basal shear stress (Lipscomb et al., 2019, just above Eq. 33)

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_basal_shear_stress'
    integer                     :: ti

    ! Add routine to path
    call init_routine( routine_name)

    do ti = self%mesh%ti1, self%mesh%ti2
      ! Lipscomb et al., 2019, just above Eq. 33
      self%tau_bx_b( ti) = self%u_vav_b( ti) * self%beta_eff_b( ti)
      self%tau_by_b( ti) = self%v_vav_b( ti) * self%beta_eff_b( ti)
    end do

    call checksum( self%mesh%pai_Tri, self%tau_bx_b, 'DIVA%tau_bx_b')
    call checksum( self%mesh%pai_Tri, self%tau_by_b, 'DIVA%tau_by_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_basal_shear_stress

  subroutine calc_3D_velocities( self)
    !< Calculate 3D velocities (Lipscomb et al., 2019, Eq. 29)

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_3D_velocities'
    integer                     :: ti,k

    ! Add routine to path
    call init_routine( routine_name)

    if (C%choice_sliding_law == 'no_sliding') then
      ! Exception for the case of no sliding

      do ti = self%mesh%ti1, self%mesh%ti2
      do k = 1, self%mesh%nz
        ! Lipscomb et al., 2019, Eq. 29, and text between Eqs. 33 and 34
        self%u_3D_b( ti,k) = self%tau_bx_b( ti) * self%F1_3D_b( ti,k)
        self%v_3D_b( ti,k) = self%tau_by_b( ti) * self%F1_3D_b( ti,k)
      end do
      end do

    else

      do ti = self%mesh%ti1, self%mesh%ti2
      do k = 1, self%mesh%nz
        ! Lipscomb et al., 2019, Eq. 29
        self%u_3D_b( ti,k) = self%u_base_b( ti) * (1._dp + self%basal_friction_coefficient_b( ti) * self%F1_3D_b( ti,k))
        self%v_3D_b( ti,k) = self%v_base_b( ti) * (1._dp + self%basal_friction_coefficient_b( ti) * self%F1_3D_b( ti,k))
      end do
      end do

    end if

    call checksum( self%mesh%pai_Tri, self%u_3D_b, 'DIVA%u_3D_b')
    call checksum( self%mesh%pai_Tri, self%v_3D_b, 'DIVA%v_3D_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_3D_velocities

  ! == Restart NetCDF files

  subroutine write_to_restart_file_DIVA( self, time)
    ! Write to the restart NetCDF file for the DIVA solver

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(in   ) :: self
    real(dp),                                       intent(in   ) :: time

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'write_to_restart_file_DIVA'
    integer                                          :: ncid
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2) :: u_vav_b_prev_loc
    real(dp), dimension(self%mesh%ti1:self%mesh%ti2) :: v_vav_b_prev_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to DIVA restart file "' // &
      UPSY%stru%colour_string( trim( self%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( self%restart_filename, ncid)

    ! Write the time to the file
    call write_time_to_file( self%restart_filename, ncid, time)

    ! Solution
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'u_vav_b'                     , self%u_vav_b)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'v_vav_b'                     , self%v_vav_b)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'u_base_b'                    , self%u_base_b)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'v_base_b'                    , self%v_base_b)
    call write_to_field_multopt_mesh_dp_3D_b( self%mesh, self%restart_filename, ncid, 'u_3D_b'                      , self%u_3D_b)
    call write_to_field_multopt_mesh_dp_3D_b( self%mesh, self%restart_filename, ncid, 'v_3D_b'                      , self%v_3D_b)

    ! Intermediate data fields
    call write_to_field_multopt_mesh_dp_2D  ( self%mesh, self%restart_filename, ncid, 'du_dx_a'                     , self%du_dx_a)
    call write_to_field_multopt_mesh_dp_2D  ( self%mesh, self%restart_filename, ncid, 'du_dy_a'                     , self%du_dy_a)
    call write_to_field_multopt_mesh_dp_2D  ( self%mesh, self%restart_filename, ncid, 'dv_dx_a'                     , self%dv_dx_a)
    call write_to_field_multopt_mesh_dp_2D  ( self%mesh, self%restart_filename, ncid, 'dv_dy_a'                     , self%dv_dy_a)
    call write_to_field_multopt_mesh_dp_3D  ( self%mesh, self%restart_filename, ncid, 'du_dz_3D_a'                  , self%du_dz_3D_a)
    call write_to_field_multopt_mesh_dp_3D  ( self%mesh, self%restart_filename, ncid, 'dv_dz_3D_a'                  , self%dv_dz_3D_a)
    call write_to_field_multopt_mesh_dp_3D  ( self%mesh, self%restart_filename, ncid, 'eta_3D_a'                    , self%eta_3D_a)
    call write_to_field_multopt_mesh_dp_3D_b( self%mesh, self%restart_filename, ncid, 'eta_3D_b'                    , self%eta_3D_b)
    call write_to_field_multopt_mesh_dp_2D  ( self%mesh, self%restart_filename, ncid, 'eta_vav_a'                   , self%eta_vav_a)
    call write_to_field_multopt_mesh_dp_2D  ( self%mesh, self%restart_filename, ncid, 'N_a'                         , self%N_a)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'N_b'                         , self%N_b)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'dN_dx_b'                     , self%dN_dx_b)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'dN_dy_b'                     , self%dN_dy_b)
    call write_to_field_multopt_mesh_dp_3D  ( self%mesh, self%restart_filename, ncid, 'F1_3D_a'                     , self%F1_3D_a)
    call write_to_field_multopt_mesh_dp_3D  ( self%mesh, self%restart_filename, ncid, 'F2_3D_a'                     , self%F2_3D_a)
    call write_to_field_multopt_mesh_dp_3D_b( self%mesh, self%restart_filename, ncid, 'F1_3D_b'                     , self%F1_3D_b)
    call write_to_field_multopt_mesh_dp_3D_b( self%mesh, self%restart_filename, ncid, 'F2_3D_b'                     , self%F2_3D_b)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'basal_friction_coefficient_b', self%basal_friction_coefficient_b)
    call write_to_field_multopt_mesh_dp_2D  ( self%mesh, self%restart_filename, ncid, 'beta_eff_a'                  , self%beta_eff_a)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'beta_eff_b'                  , self%beta_eff_b)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'tau_bx_b'                    , self%tau_bx_b)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'tau_by_b'                    , self%tau_by_b)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'tau_dx_b'                    , self%tau_dx_b)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'tau_dy_b'                    , self%tau_dy_b)
    u_vav_b_prev_loc = self%u_vav_b_prev( self%mesh%ti1:self%mesh%ti2)
    v_vav_b_prev_loc = self%v_vav_b_prev( self%mesh%ti1:self%mesh%ti2)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'u_vav_b_prev'                , u_vav_b_prev_loc)
    call write_to_field_multopt_mesh_dp_2D_b( self%mesh, self%restart_filename, ncid, 'v_vav_b_prev'                , v_vav_b_prev_loc)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file_DIVA

  subroutine create_restart_file_DIVA( self)
    ! Create a restart NetCDF file for the DIVA solver
    ! Includes generation of the procedural filename (e.g. "restart_DIVA_00001.nc")

    ! In/output variables:
    class(type_momentum_balance_solver_plain_DIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'create_restart_file_DIVA'
    character(len=:), allocatable :: filename_base
    integer                       :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'restart_ice_velocity_DIVA'
    call generate_filename_XXXXXdotnc( filename_base, self%restart_filename)

    ! Print to terminal
    if (par%primary) WRITE(0,'(A)') '   Creating DIVA restart file "' // &
      UPSY%stru%colour_string( TRIM( self%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( self%restart_filename, ncid)

    ! Set up the mesh in the file
    call setup_mesh_in_netcdf_file( self%restart_filename, ncid, self%mesh)

    ! Add a time dimension to the file
    call add_time_dimension_to_file( self%restart_filename, ncid)

    ! Add a zeta dimension to the file
    call add_zeta_dimension_to_file( self%restart_filename, ncid, self%mesh%zeta)

    ! Solution
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'u_vav_b'                     , long_name = 'Vertically averaged horizontal ice velocity in the x-direction', units = 'm/yr')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'v_vav_b'                     , long_name = 'Vertically averaged horizontal ice velocity in the y-direction', units = 'm/yr')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'u_base_b'                    , long_name = 'Basal ice velocity in the x-direction', units = 'm/yr')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'v_base_b'                    , long_name = 'Basal ice velocity in the y-direction', units = 'm/yr')
    call add_field_mesh_dp_3D_b( self%restart_filename, ncid, 'u_3D_b'                      , long_name = '3-D horizontal ice velocity in the x-direction', units = 'm/yr')
    call add_field_mesh_dp_3D_b( self%restart_filename, ncid, 'v_3D_b'                      , long_name = '3-D horizontal ice velocity in the y-direction', units = 'm/yr')

    ! Intermediate data fields
    call add_field_mesh_dp_2D  ( self%restart_filename, ncid, 'du_dx_a'                     , long_name = 'Vertically averaged xx strain rate', units = 'yr^-1')
    call add_field_mesh_dp_2D  ( self%restart_filename, ncid, 'du_dy_a'                     , long_name = 'Vertically averaged xy strain rate', units = 'yr^-1')
    call add_field_mesh_dp_2D  ( self%restart_filename, ncid, 'dv_dx_a'                     , long_name = 'Vertically averaged yx strain rate', units = 'yr^-1')
    call add_field_mesh_dp_2D  ( self%restart_filename, ncid, 'dv_dy_a'                     , long_name = 'Vertically averaged yy strain rate', units = 'yr^-1')
    call add_field_mesh_dp_3D  ( self%restart_filename, ncid, 'du_dz_3D_a'                  , long_name = '3-D xz strain rate', units = 'yr^-1')
    call add_field_mesh_dp_3D  ( self%restart_filename, ncid, 'dv_dz_3D_a'                  , long_name = '3-D yz strain rate', units = 'yr^-1')
    call add_field_mesh_dp_3D  ( self%restart_filename, ncid, 'eta_3D_a'                    , long_name = '3-D effective viscosity on a-grid')
    call add_field_mesh_dp_3D_b( self%restart_filename, ncid, 'eta_3D_b'                    , long_name = '3-D effective viscosity on b-grid')
    call add_field_mesh_dp_2D  ( self%restart_filename, ncid, 'eta_vav_a'                   , long_name = 'Vertically averaged effective viscosity on a-grid')
    call add_field_mesh_dp_2D  ( self%restart_filename, ncid, 'N_a'                         , long_name = 'Product term N = eta * H on a-grid')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'N_b'                         , long_name = 'Product term N = eta * H on b-grid')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'dN_dx_b'                     , long_name = 'Gradient of N in x-direction on b-grid')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'dN_dy_b'                     , long_name = 'Gradient of N in y-direction on b-grid')
    call add_field_mesh_dp_3D  ( self%restart_filename, ncid, 'F1_3D_a'                     , long_name = '3-D F_1 integral on a-grid')
    call add_field_mesh_dp_3D  ( self%restart_filename, ncid, 'F2_3D_a'                     , long_name = '3-D F_2 integral on a-grid')
    call add_field_mesh_dp_3D_b( self%restart_filename, ncid, 'F1_3D_b'                     , long_name = '3-D F_1 integral on b-grid')
    call add_field_mesh_dp_3D_b( self%restart_filename, ncid, 'F2_3D_b'                     , long_name = '3-D F_2 integral on b-grid')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'basal_friction_coefficient_b', long_name = 'Basal friction coefficient on b-grid')
    call add_field_mesh_dp_2D  ( self%restart_filename, ncid, 'beta_eff_a'                  , long_name = 'Beta_eff on a-grid')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'beta_eff_b'                  , long_name = 'Beta_eff on b-grid')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'tau_bx_b'                    , long_name = 'Basal shear stress in the x-direction', units = 'Pa')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'tau_by_b'                    , long_name = 'Basal shear stress in the y-direction', units = 'Pa')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'tau_dx_b'                    , long_name = 'Driving stress in the x-direction', units = 'Pa')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'tau_dy_b'                    , long_name = 'Driving stress in the y-direction', units = 'Pa')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'u_vav_b_prev'                , long_name = 'Previous iteration of u_b', units = 'm yr^-1')
    call add_field_mesh_dp_2D_b( self%restart_filename, ncid, 'v_vav_b_prev'                , long_name = 'Previous iteration of v_b', units = 'm yr^-1')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_DIVA

end module momentum_balance_solver_plain_DIVA