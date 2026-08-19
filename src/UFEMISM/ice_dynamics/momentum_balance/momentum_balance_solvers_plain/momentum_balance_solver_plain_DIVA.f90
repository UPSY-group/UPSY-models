module momentum_balance_solver_plain_DIVA

  ! Routines for calculating ice velocities using the Depth-Integrated Viscosity Approximation (DIVA)

  use mpi_basic, only: par
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data, type_ice_velocity_solver_DIVA
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use netcdf_io_main
  use mesh_disc_apply_operators, only: map_a_b_2D, map_a_b_3D, map_b_a_2D, map_b_a_3D
  use reallocate_mod, only: reallocate_bounds, reallocate_clean
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D, &
    map_from_mesh_to_mesh_with_reallocation_3D
  use bed_roughness_model_types, only: type_bed_roughness_model
  use DIVA_solver_infinite_slab, only: solve_DIVA_infinite_slab
  use DIVA_solver_ocean_pressure, only: solve_DIVA_ocean_pressure
  use mpi_distributed_memory, only: gather_to_all

  implicit none

contains

  ! == Main routines

  subroutine initialise_DIVA_solver( mesh, DIVA, region_name)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(  out) :: DIVA
    character(len=3),                    intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_DIVA_solver'
    character(len=256)             :: choice_initial_velocity

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    call allocate_DIVA_solver( mesh, DIVA)

    ! Determine the choice of initial velocities for this model region
    select case (region_name)
    case default
      call crash('unknown model region "' // region_name // '"!')
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
      DIVA%u_vav_b = 0._dp
      DIVA%v_vav_b = 0._dp
    case ('read_from_file')
      call initialise_DIVA_velocities_from_file( mesh, DIVA, region_name)
    end select

    ! Set tolerances for PETSc matrix solver for the linearised DIVA
    DIVA%PETSc_rtol   = C%stress_balance_PETSc_rtol
    DIVA%PETSc_abstol = C%stress_balance_PETSc_abstol

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_DIVA_solver

  subroutine solve_DIVA( mesh, ice, geom, bed_roughness, DIVA, n_visc_its, n_Axb_its, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    !< Calculate ice velocities by solving the Depth-Integrated Viscosity Approximation

    ! In/output variables:
    type(type_mesh),                      intent(in   ) :: mesh
    class(atype_ice_model_data),          intent(inout) :: ice
    class(atype_ice_geometry_model_data), intent(in   ) :: geom
    type(type_bed_roughness_model),       intent(in   ) :: bed_roughness
    type(type_ice_velocity_solver_DIVA),  intent(inout) :: DIVA
    integer,                              intent(  out) :: n_visc_its               ! Number of non-linear viscosity iterations
    integer,                              intent(  out) :: n_Axb_its                ! Number of iterations in iterative solver for linearised momentum balance
    integer,  dimension(:), optional,     intent(in   ) :: BC_prescr_mask_b         ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:), optional,     intent(in   ) :: BC_prescr_u_b            ! Prescribed velocities in the x-direction
    real(dp), dimension(:), optional,     intent(in   ) :: BC_prescr_v_b            ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'solve_DIVA'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%BC_ice_front)
    case default
      call crash('unknown BC_ice_front "' // trim( C%BC_ice_front) // '"')
    case ('infinite_slab')
      call solve_DIVA_infinite_slab( mesh, ice, geom, bed_roughness, DIVA, n_visc_its, n_Axb_its, &
        BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    case ('ocean_pressure')
      call solve_DIVA_ocean_pressure( mesh, ice, geom, bed_roughness, DIVA, n_visc_its, n_Axb_its, &
        BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_DIVA

  subroutine remap_DIVA_solver( mesh_old, mesh_new, DIVA)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh_old
    type(type_mesh),                     intent(in   ) :: mesh_new
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'remap_DIVA_solver'
    real(dp), dimension(:  ), allocatable :: u_vav_a
    real(dp), dimension(:  ), allocatable :: v_vav_a
    real(dp), dimension(:  ), allocatable :: tau_bx_a
    real(dp), dimension(:  ), allocatable :: tau_by_a
    real(dp), dimension(:,:), allocatable :: eta_3D_a
    real(dp), dimension(:,:), allocatable :: u_3D_a
    real(dp), dimension(:,:), allocatable :: v_3D_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap the fields that are re-used during the viscosity iteration
    ! ================================================================

    ! allocate memory for velocities on the a-grid (vertices)
    allocate( u_vav_a ( mesh_old%vi1: mesh_old%vi2             ))
    allocate( v_vav_a ( mesh_old%vi1: mesh_old%vi2             ))
    allocate( tau_bx_a( mesh_old%vi1: mesh_old%vi2             ))
    allocate( tau_by_a( mesh_old%vi1: mesh_old%vi2             ))
    allocate( eta_3D_a( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    allocate( u_3D_a  ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))
    allocate( v_3D_a  ( mesh_old%vi1: mesh_old%vi2, mesh_old%nz))

    ! Map data from the triangles of the old mesh to the vertices of the old mesh
    call map_b_a_2D( mesh_old, DIVA%u_vav_b , u_vav_a )
    call map_b_a_2D( mesh_old, DIVA%v_vav_b , v_vav_a )
    call map_b_a_2D( mesh_old, DIVA%tau_bx_b, tau_bx_a)
    call map_b_a_2D( mesh_old, DIVA%tau_by_b, tau_by_a)
    call map_b_a_3D( mesh_old, DIVA%eta_3D_b, eta_3D_a)
    call map_b_a_3D( mesh_old, DIVA%u_3D_b  , u_3D_a  )
    call map_b_a_3D( mesh_old, DIVA%v_3D_b  , v_3D_a  )

    ! Remap data from the vertices of the old mesh to the vertices of the new mesh
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, u_vav_a , '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, v_vav_a , '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, tau_bx_a, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, tau_by_a, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, eta_3D_a, '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, u_3D_a  , '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, C%output_dir, v_3D_a  , '2nd_order_conservative')

    ! reallocate memory for the data on the triangles
    call reallocate_bounds( DIVA%u_vav_b  , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%v_vav_b  , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%tau_bx_b , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%tau_by_b , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%eta_3D_b , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( DIVA%u_3D_b   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( DIVA%v_3D_b   , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)

    ! Map data from the vertices of the new mesh to the triangles of the new mesh
    call map_a_b_2D( mesh_new, u_vav_a , DIVA%u_vav_b )
    call map_a_b_2D( mesh_new, v_vav_a , DIVA%v_vav_b )
    call map_a_b_2D( mesh_new, tau_bx_a, DIVA%tau_bx_b)
    call map_a_b_2D( mesh_new, tau_by_a, DIVA%tau_by_b)
    call map_a_b_3D( mesh_new, eta_3D_a, DIVA%eta_3D_b)
    call map_a_b_3D( mesh_new, u_3D_a  , DIVA%u_3D_b  )
    call map_a_b_3D( mesh_new, v_3D_a  , DIVA%v_3D_b  )

    ! reallocate everything else
    ! ==========================

    ! call reallocate_bounds( DIVA%u_vav_b                     , mesh_new%ti1, mesh_new%ti2             )           ! [m yr^-1] 2-D vertically averaged horizontal ice velocity
    ! call reallocate_bounds( DIVA%v_vav_b                     , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%u_base_b                    , mesh_new%ti1, mesh_new%ti2             )           ! [m yr^-1] 2-D horizontal ice velocity at the ice base
    call reallocate_bounds( DIVA%v_base_b                    , mesh_new%ti1, mesh_new%ti2             )
    ! call reallocate_bounds( DIVA%u_3D_b                      , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)           ! [m yr^-1] 3-D horizontal ice velocity
    ! call reallocate_bounds( DIVA%v_3D_b                      , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( DIVA%du_dx_a                     , mesh_new%vi1, mesh_new%vi2             )           ! [yr^-1] 2-D horizontal strain rates
    call reallocate_bounds( DIVA%du_dy_a                     , mesh_new%vi1, mesh_new%vi2             )
    call reallocate_bounds( DIVA%dv_dx_a                     , mesh_new%vi1, mesh_new%vi2             )
    call reallocate_bounds( DIVA%dv_dy_a                     , mesh_new%vi1, mesh_new%vi2             )
    call reallocate_bounds( DIVA%du_dz_3D_a                  , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)           ! [yr^-1] 3-D vertical shear strain rates
    call reallocate_bounds( DIVA%dv_dz_3D_a                  , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    call reallocate_bounds( DIVA%eta_3D_a                    , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)           ! Effective viscosity
    ! call reallocate_bounds( DIVA%eta_3D_b                    , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( DIVA%eta_vav_a                   , mesh_new%vi1, mesh_new%vi2             )
    call reallocate_bounds( DIVA%N_a                         , mesh_new%vi1, mesh_new%vi2             )           ! Product term N = eta * H
    call reallocate_bounds( DIVA%N_b                         , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%dN_dx_b                     , mesh_new%ti1, mesh_new%ti2             )           ! Gradients of N
    call reallocate_bounds( DIVA%dN_dy_b                     , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%F1_3D_a                     , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)           ! F-integrals
    call reallocate_bounds( DIVA%F2_3D_a                     , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    call reallocate_bounds( DIVA%F1_3D_b                     , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( DIVA%F2_3D_b                     , mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( DIVA%basal_friction_coefficient_b, mesh_new%ti1, mesh_new%ti2             )           ! Basal friction coefficient (basal_shear_stress = u * basal_friction_coefficient)
    call reallocate_bounds( DIVA%beta_eff_a                  , mesh_new%vi1, mesh_new%vi2             )           ! "Effective" friction coefficient (turning the SSA into the DIVA)
    call reallocate_bounds( DIVA%beta_eff_b                  , mesh_new%ti1, mesh_new%ti2             )
    ! call reallocate_bounds( DIVA%tau_bx_b                    , mesh_new%ti1, mesh_new%ti2             )           ! Basal shear stress
    ! call reallocate_bounds( DIVA%tau_by_b                    , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_bounds( DIVA%tau_dx_b                    , mesh_new%ti1, mesh_new%ti2             )           ! Driving stress
    call reallocate_bounds( DIVA%tau_dy_b                    , mesh_new%ti1, mesh_new%ti2             )
    call reallocate_clean ( DIVA%u_b_prev                    , mesh_new%nTri                          )           ! Velocity solution from previous viscosity iteration
    call reallocate_clean ( DIVA%v_b_prev                    , mesh_new%nTri                          )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_DIVA_solver

  ! == Initialisation

  subroutine initialise_DIVA_velocities_from_file( mesh, DIVA, region_name)
    !< Initialise the velocities for the DIVA solver from an external NetCDF file

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA
    character(len=3),                    intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_DIVA_velocities_from_file'
    real(dp)                       :: dummy1
    character(len=1024)            :: filename
    real(dp)                       :: timeframe
    real(dp), dimension(mesh%ti1:mesh%ti2) :: u_b_prev_loc
    real(dp), dimension(mesh%ti1:mesh%ti2) :: v_b_prev_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! To prevent compiler warnings
    dummy1 = mesh%xmin

    ! Determine the filename and timeframe to read for this model region
    select case (region_name)
    case default
      call crash('unknown model region "' // region_name // '"!')
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
    call read_field_from_mesh_file_dp_2D_b( filename, 'u_vav_b'                     , DIVA%u_vav_b                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'v_vav_b'                     , DIVA%v_vav_b                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'u_base_b'                    , DIVA%u_base_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'v_base_b'                    , DIVA%v_base_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D_b( filename, 'u_3D_b'                      , DIVA%u_3D_b                      , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D_b( filename, 'v_3D_b'                      , DIVA%v_3D_b                      , time_to_read = timeframe)

    ! Intermediate data fields
    call read_field_from_mesh_file_dp_2D  ( filename, 'du_dx_a'                     , DIVA%du_dx_a                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D  ( filename, 'du_dy_a'                     , DIVA%du_dy_a                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D  ( filename, 'dv_dx_a'                     , DIVA%dv_dx_a                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D  ( filename, 'dv_dy_a'                     , DIVA%dv_dy_a                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D  ( filename, 'du_dz_3D_a'                  , DIVA%du_dz_3D_a                  , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D  ( filename, 'dv_dz_3D_a'                  , DIVA%dv_dz_3D_a                  , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D  ( filename, 'eta_3D_a'                    , DIVA%eta_3D_a                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D_b( filename, 'eta_3D_b'                    , DIVA%eta_3D_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D  ( filename, 'eta_vav_a'                   , DIVA%eta_vav_a                   , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D  ( filename, 'N_a'                         , DIVA%N_a                         , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'N_b'                         , DIVA%N_b                         , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'dN_dx_b'                     , DIVA%dN_dx_b                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'dN_dy_b'                     , DIVA%dN_dy_b                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D  ( filename, 'F1_3D_a'                     , DIVA%F1_3D_a                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D  ( filename, 'F2_3D_a'                     , DIVA%F2_3D_a                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D_b( filename, 'F1_3D_b'                     , DIVA%F1_3D_b                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_3D_b( filename, 'F2_3D_b'                     , DIVA%F2_3D_b                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'basal_friction_coefficient_b', DIVA%basal_friction_coefficient_b, time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D  ( filename, 'beta_eff_a'                  , DIVA%beta_eff_a                  , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'beta_eff_b'                  , DIVA%beta_eff_b                  , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'tau_bx_b'                    , DIVA%tau_bx_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'tau_by_b'                    , DIVA%tau_by_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'tau_dx_b'                    , DIVA%tau_dx_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'tau_dy_b'                    , DIVA%tau_dy_b                    , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'u_b_prev'                    , u_b_prev_loc                     , time_to_read = timeframe)
    call read_field_from_mesh_file_dp_2D_b( filename, 'v_b_prev'                    , v_b_prev_loc                     , time_to_read = timeframe)

    call gather_to_all( u_b_prev_loc, DIVA%u_b_prev)
    call gather_to_all( v_b_prev_loc, DIVA%v_b_prev)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_DIVA_velocities_from_file

  subroutine allocate_DIVA_solver( mesh, DIVA)
    ! allocate memory the DIVA solver

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(  out) :: DIVA

    ! Local variables:
    character(len=1024), parameter                     :: routine_name = 'allocate_DIVA_solver'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    allocate( DIVA%u_vav_b(  mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%v_vav_b(  mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%u_base_b( mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%v_base_b( mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%u_3D_b(   mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( DIVA%v_3D_b(   mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)

    ! Intermediate data fields
    allocate( DIVA%du_dx_a(                      mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%du_dy_a(                      mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%dv_dx_a(                      mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%dv_dy_a(                      mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%du_dz_3D_a(                   mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( DIVA%dv_dz_3D_a(                   mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( DIVA%eta_3D_a(                     mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( DIVA%eta_3D_b(                     mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( DIVA%eta_vav_a(                    mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%N_a(                          mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%N_b(                          mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%dN_dx_b(                      mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%dN_dy_b(                      mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%F1_3D_a(                      mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( DIVA%F2_3D_a(                      mesh%vi1:mesh%vi2,mesh%nz), source = 0._dp)
    allocate( DIVA%F1_3D_b(                      mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( DIVA%F2_3D_b(                      mesh%ti1:mesh%ti2,mesh%nz), source = 0._dp)
    allocate( DIVA%basal_friction_coefficient_b( mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%beta_eff_a(                   mesh%vi1:mesh%vi2        ), source = 0._dp)
    allocate( DIVA%beta_eff_b(                   mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%tau_bx_b(                     mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%tau_by_b(                     mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%tau_dx_b(                     mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%tau_dy_b(                     mesh%ti1:mesh%ti2        ), source = 0._dp)
    allocate( DIVA%u_b_prev(                     mesh%nTri                ), source = 0._dp)
    allocate( DIVA%v_b_prev(                     mesh%nTri                ), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_DIVA_solver

  ! == Restart NetCDF files

  subroutine write_to_restart_file_DIVA( mesh, DIVA, time)
    ! Write to the restart NetCDF file for the DIVA solver

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(in   ) :: DIVA
    real(dp),                            intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_restart_file_DIVA'
    integer                        :: ncid
    real(dp), dimension(mesh%ti1:mesh%ti2) :: u_b_prev_loc
    real(dp), dimension(mesh%ti1:mesh%ti2) :: v_b_prev_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to DIVA restart file "' // &
      UPSY%stru%colour_string( trim( DIVA%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( DIVA%restart_filename, ncid)

    ! Write the time to the file
    call write_time_to_file( DIVA%restart_filename, ncid, time)

    ! Solution
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'u_vav_b'                     , DIVA%u_vav_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'v_vav_b'                     , DIVA%v_vav_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'u_base_b'                    , DIVA%u_base_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'v_base_b'                    , DIVA%v_base_b)
    call write_to_field_multopt_mesh_dp_3D_b( mesh, DIVA%restart_filename, ncid, 'u_3D_b'                      , DIVA%u_3D_b)
    call write_to_field_multopt_mesh_dp_3D_b( mesh, DIVA%restart_filename, ncid, 'v_3D_b'                      , DIVA%v_3D_b)

    ! Intermediate data fields
    call write_to_field_multopt_mesh_dp_2D  ( mesh, DIVA%restart_filename, ncid, 'du_dx_a'                     , DIVA%du_dx_a)
    call write_to_field_multopt_mesh_dp_2D  ( mesh, DIVA%restart_filename, ncid, 'du_dy_a'                     , DIVA%du_dy_a)
    call write_to_field_multopt_mesh_dp_2D  ( mesh, DIVA%restart_filename, ncid, 'dv_dx_a'                     , DIVA%dv_dx_a)
    call write_to_field_multopt_mesh_dp_2D  ( mesh, DIVA%restart_filename, ncid, 'dv_dy_a'                     , DIVA%dv_dy_a)
    call write_to_field_multopt_mesh_dp_3D  ( mesh, DIVA%restart_filename, ncid, 'du_dz_3D_a'                  , DIVA%du_dz_3D_a)
    call write_to_field_multopt_mesh_dp_3D  ( mesh, DIVA%restart_filename, ncid, 'dv_dz_3D_a'                  , DIVA%dv_dz_3D_a)
    call write_to_field_multopt_mesh_dp_3D  ( mesh, DIVA%restart_filename, ncid, 'eta_3D_a'                    , DIVA%eta_3D_a)
    call write_to_field_multopt_mesh_dp_3D_b( mesh, DIVA%restart_filename, ncid, 'eta_3D_b'                    , DIVA%eta_3D_b)
    call write_to_field_multopt_mesh_dp_2D  ( mesh, DIVA%restart_filename, ncid, 'eta_vav_a'                   , DIVA%eta_vav_a)
    call write_to_field_multopt_mesh_dp_2D  ( mesh, DIVA%restart_filename, ncid, 'N_a'                         , DIVA%N_a)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'N_b'                         , DIVA%N_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'dN_dx_b'                     , DIVA%dN_dx_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'dN_dy_b'                     , DIVA%dN_dy_b)
    call write_to_field_multopt_mesh_dp_3D  ( mesh, DIVA%restart_filename, ncid, 'F1_3D_a'                     , DIVA%F1_3D_a)
    call write_to_field_multopt_mesh_dp_3D  ( mesh, DIVA%restart_filename, ncid, 'F2_3D_a'                     , DIVA%F2_3D_a)
    call write_to_field_multopt_mesh_dp_3D_b( mesh, DIVA%restart_filename, ncid, 'F1_3D_b'                     , DIVA%F1_3D_b)
    call write_to_field_multopt_mesh_dp_3D_b( mesh, DIVA%restart_filename, ncid, 'F2_3D_b'                     , DIVA%F2_3D_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'basal_friction_coefficient_b', DIVA%basal_friction_coefficient_b)
    call write_to_field_multopt_mesh_dp_2D  ( mesh, DIVA%restart_filename, ncid, 'beta_eff_a'                  , DIVA%beta_eff_a)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'beta_eff_b'                  , DIVA%beta_eff_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'tau_bx_b'                    , DIVA%tau_bx_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'tau_by_b'                    , DIVA%tau_by_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'tau_dx_b'                    , DIVA%tau_dx_b)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'tau_dy_b'                    , DIVA%tau_dy_b)
    u_b_prev_loc = DIVA%u_b_prev( mesh%ti1:mesh%ti2)
    v_b_prev_loc = DIVA%v_b_prev( mesh%ti1:mesh%ti2)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'u_b_prev'                    , u_b_prev_loc)
    call write_to_field_multopt_mesh_dp_2D_b( mesh, DIVA%restart_filename, ncid, 'v_b_prev'                    , v_b_prev_loc)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file_DIVA

  subroutine create_restart_file_DIVA( mesh, DIVA)
    ! Create a restart NetCDF file for the DIVA solver
    ! Includes generation of the procedural filename (e.g. "restart_DIVA_00001.nc")

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_restart_file_DIVA'
    character(len=1024)            :: filename_base
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'restart_ice_velocity_DIVA'
    call generate_filename_XXXXXdotnc( filename_base, DIVA%restart_filename)

    ! Print to terminal
    if (par%primary) WRITE(0,'(A)') '   Creating DIVA restart file "' // &
      UPSY%stru%colour_string( TRIM( DIVA%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( DIVA%restart_filename, ncid)

    ! Set up the mesh in the file
    call setup_mesh_in_netcdf_file( DIVA%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    call add_time_dimension_to_file( DIVA%restart_filename, ncid)

    ! Add a zeta dimension to the file
    call add_zeta_dimension_to_file( DIVA%restart_filename, ncid, mesh%zeta)

    ! Solution
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'u_vav_b'                     , long_name = 'Vertically averaged horizontal ice velocity in the x-direction', units = 'm/yr')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'v_vav_b'                     , long_name = 'Vertically averaged horizontal ice velocity in the y-direction', units = 'm/yr')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'u_base_b'                    , long_name = 'Basal ice velocity in the x-direction', units = 'm/yr')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'v_base_b'                    , long_name = 'Basal ice velocity in the y-direction', units = 'm/yr')
    call add_field_mesh_dp_3D_b( DIVA%restart_filename, ncid, 'u_3D_b'                      , long_name = '3-D horizontal ice velocity in the x-direction', units = 'm/yr')
    call add_field_mesh_dp_3D_b( DIVA%restart_filename, ncid, 'v_3D_b'                      , long_name = '3-D horizontal ice velocity in the y-direction', units = 'm/yr')

    ! Intermediate data fields
    call add_field_mesh_dp_2D  ( DIVA%restart_filename, ncid, 'du_dx_a'                     , long_name = 'Vertically averaged xx strain rate', units = 'yr^-1')
    call add_field_mesh_dp_2D  ( DIVA%restart_filename, ncid, 'du_dy_a'                     , long_name = 'Vertically averaged xy strain rate', units = 'yr^-1')
    call add_field_mesh_dp_2D  ( DIVA%restart_filename, ncid, 'dv_dx_a'                     , long_name = 'Vertically averaged yx strain rate', units = 'yr^-1')
    call add_field_mesh_dp_2D  ( DIVA%restart_filename, ncid, 'dv_dy_a'                     , long_name = 'Vertically averaged yy strain rate', units = 'yr^-1')
    call add_field_mesh_dp_3D  ( DIVA%restart_filename, ncid, 'du_dz_3D_a'                  , long_name = '3-D xz strain rate', units = 'yr^-1')
    call add_field_mesh_dp_3D  ( DIVA%restart_filename, ncid, 'dv_dz_3D_a'                  , long_name = '3-D yz strain rate', units = 'yr^-1')
    call add_field_mesh_dp_3D  ( DIVA%restart_filename, ncid, 'eta_3D_a'                    , long_name = '3-D effective viscosity on a-grid')
    call add_field_mesh_dp_3D_b( DIVA%restart_filename, ncid, 'eta_3D_b'                    , long_name = '3-D effective viscosity on b-grid')
    call add_field_mesh_dp_2D  ( DIVA%restart_filename, ncid, 'eta_vav_a'                   , long_name = 'Vertically averaged effective viscosity on a-grid')
    call add_field_mesh_dp_2D  ( DIVA%restart_filename, ncid, 'N_a'                         , long_name = 'Product term N = eta * H on a-grid')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'N_b'                         , long_name = 'Product term N = eta * H on b-grid')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'dN_dx_b'                     , long_name = 'Gradient of N in x-direction on b-grid')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'dN_dy_b'                     , long_name = 'Gradient of N in y-direction on b-grid')
    call add_field_mesh_dp_3D  ( DIVA%restart_filename, ncid, 'F1_3D_a'                     , long_name = '3-D F_1 integral on a-grid')
    call add_field_mesh_dp_3D  ( DIVA%restart_filename, ncid, 'F2_3D_a'                     , long_name = '3-D F_2 integral on a-grid')
    call add_field_mesh_dp_3D_b( DIVA%restart_filename, ncid, 'F1_3D_b'                     , long_name = '3-D F_1 integral on b-grid')
    call add_field_mesh_dp_3D_b( DIVA%restart_filename, ncid, 'F2_3D_b'                     , long_name = '3-D F_2 integral on b-grid')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'basal_friction_coefficient_b', long_name = 'Basal friction coefficient on b-grid')
    call add_field_mesh_dp_2D  ( DIVA%restart_filename, ncid, 'beta_eff_a'                  , long_name = 'Beta_eff on a-grid')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'beta_eff_b'                  , long_name = 'Beta_eff on b-grid')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'tau_bx_b'                    , long_name = 'Basal shear stress in the x-direction', units = 'Pa')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'tau_by_b'                    , long_name = 'Basal shear stress in the y-direction', units = 'Pa')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'tau_dx_b'                    , long_name = 'Driving stress in the x-direction', units = 'Pa')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'tau_dy_b'                    , long_name = 'Driving stress in the y-direction', units = 'Pa')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'u_b_prev'                    , long_name = 'Previous iteration of u_b', units = 'm yr^-1')
    call add_field_mesh_dp_2D_b( DIVA%restart_filename, ncid, 'v_b_prev'                    , long_name = 'Previous iteration of v_b', units = 'm yr^-1')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_DIVA

end module momentum_balance_solver_plain_DIVA
