module ice_velocities_main

  !< Contains all the routines needed to solve for conservation of momentum
  !< and calculate instantaneous ice velocities for the current modelled ice-sheet geometry.

  use mpi_basic, only: par
  use precisions, only: dp
  use UPSY_main, only: UPSY
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use parameters, only: ice_density, seawater_density, pi
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data, &
    type_ice_velocity_solver_DIVA, type_ice_velocity_solver_BPA, type_ice_velocity_solver_hybrid
  use momentum_balance_solver_plain_DIVA, only: initialise_DIVA_solver, solve_DIVA, remap_DIVA_solver, &
    create_restart_file_DIVA, write_to_restart_file_DIVA
  use BPA_main, only: initialise_BPA_solver, solve_BPA, remap_BPA_solver, &
    create_restart_file_BPA, write_to_restart_file_BPA
  use hybrid_DIVA_BPA_main, only: initialise_hybrid_DIVA_BPA_solver, solve_hybrid_DIVA_BPA, &
    remap_hybrid_DIVA_BPA_solver
  use mesh_disc_apply_operators, only: ddx_a_a_2D, ddy_a_a_2D, map_b_a_2D, map_b_a_3D
  use mpi_distributed_memory, only: gather_to_all
  use mesh_zeta, only: vertical_average
  use map_velocities_to_c_grid
  use vertical_velocities, only: calc_vertical_velocities
  use bed_roughness_model_types, only: type_bed_roughness_model
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_data, only: atype_ice_velocity_model_data
  use ice_velocity_model_basic, only: atype_ice_velocity_model
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use momentum_balance_solver_SIA, only: type_momentum_balance_solver_SIA
  use momentum_balance_solver_SSA, only: type_momentum_balance_solver_SSA
  use momentum_balance_solver_SIASSA, only: type_momentum_balance_solver_SIASSA

  implicit none

contains

  ! == The main routines, to be called from the ice dynamics module

  subroutine initialise_velocity_solver( vel, ice, momentum_balance_solver, region_name)
    !< Initialise the velocity solver for the chosen Stokes approximation

    ! In/output variables:
    class(atype_ice_velocity_model),      intent(inout) :: vel
    class(atype_ice_model_data),          intent(inout) :: ice
    class(atype_momentum_balance_solver), intent(inout) :: momentum_balance_solver
    character(len=3),                     intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_velocity_solver'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) write(*,"(A)") '   Initialising ' // &
      UPSY%stru%colour_string( trim( C%choice_stress_balance_approximation),'light blue') // ' solver...'

    select case (C%choice_stress_balance_approximation)
      case default
        call crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
      case ('none')
        ! No need to do anything
      case ('SIA')
        call momentum_balance_solver%initialise()
      case ('SSA')
        call momentum_balance_solver%initialise()
      case ('SIA/SSA')
        call momentum_balance_solver%initialise()
      case ('DIVA')
        call initialise_DIVA_solver           ( vel%mesh, ice%DIVA  , region_name)
      case ('BPA')
        call initialise_BPA_solver            ( vel%mesh, ice%BPA   , region_name)
      case ('hybrid DIVA/BPA')
        call initialise_hybrid_DIVA_BPA_solver( vel%mesh, ice%hybrid, region_name)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_velocity_solver

  subroutine solve_stress_balance( mesh, ice, geom, vel, momentum_balance_solver, bed_roughness, BMB, region_name, n_visc_its, n_Axb_its, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    !< Calculate all ice velocities based on the chosen stress balance approximation

    ! In/output variables:
    type(type_mesh),                        intent(inout) :: mesh
    class(atype_ice_model_data),            intent(inout) :: ice
    class(atype_ice_geometry_model_data),   intent(in   ) :: geom
    class(atype_ice_velocity_model),        intent(inout) :: vel
    class(atype_momentum_balance_solver),   intent(inout) :: momentum_balance_solver
    type(type_bed_roughness_model),         intent(in   ) :: bed_roughness
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: BMB
    character(len=3),                       intent(in   ) :: region_name
    integer,                                intent(out)   :: n_visc_its            ! Number of non-linear viscosity iterations
    integer,                                intent(out)   :: n_Axb_its             ! Number of iterations in iterative solver for linearised momentum balance
    ! Prescribed velocities for the SSA/DIVA
    integer,  dimension(:  ), optional,     intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:  ), optional,     intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(:  ), optional,     intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    ! Prescribed velocities for the BPA
    integer,  dimension(:,:), optional,     intent(in   ) :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:,:), optional,     intent(in   ) :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
    real(dp), dimension(:,:), optional,     intent(in   ) :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter :: routine_name = 'solve_stress_balance'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_stress_balance_approximation)

      case default
        call crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')

      case ('none')

        vel%u_3D_b( mesh%ti1:mesh%ti2,:) = 0._dp
        vel%v_3D_b( mesh%ti1:mesh%ti2,:) = 0._dp

        n_visc_its = 0
        n_Axb_its  = 0

      case ('SIA')
        ! Calculate velocities according to the Shallow Ice Approximation

        call momentum_balance_solver%run( ice, geom, bed_roughness, &
          BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)

        select type (SIA => momentum_balance_solver)
        class default
          call crash('invalid momentum_balance_solver class')
        class is (type_momentum_balance_solver_SIA)
          call set_ice_velocities_to_SIA_results( mesh, vel, SIA)
        end select

        n_visc_its = 0
        n_Axb_its  = 0

      case ('SSA')
        ! Calculate velocities according to the Shallow Shelf Approximation

        call momentum_balance_solver%run( ice, geom, bed_roughness, &
          BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
        n_visc_its = momentum_balance_solver%n_visc_its
        n_Axb_its  = momentum_balance_solver%n_Axb_its

        select type (SSA => momentum_balance_solver)
        class default
          call crash('invalid momentum_balance_solver class')
        class is (type_momentum_balance_solver_SSA)
          call set_ice_velocities_to_SSA_results( mesh, vel, SSA)
        end select

      case ('SIA/SSA')
        ! Calculate velocities according to the hybrid SIA/SSA

        call momentum_balance_solver%run( ice, geom, bed_roughness, &
          BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
        n_visc_its = momentum_balance_solver%n_visc_its
        n_Axb_its  = momentum_balance_solver%n_Axb_its

        select type (SIASSA => momentum_balance_solver)
        class default
          call crash('invalid momentum_balance_solver class')
        class is (type_momentum_balance_solver_SIASSA)
          call set_ice_velocities_to_SIASSA_results( mesh, vel, SIASSA)
        end select

      case ('DIVA')
        ! Calculate velocities according to the Depth-Integrated Viscosity Approximation

        call solve_DIVA( mesh, ice, geom, bed_roughness, ice%DIVA, &
          n_visc_its, n_Axb_its, &
          BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
        call set_ice_velocities_to_DIVA_results( mesh, ice, vel, ice%DIVA)

      case ('BPA')
        ! Calculate velocities according to the Blatter-Pattyn Approximation

        call solve_BPA( mesh, ice, geom, bed_roughness, ice%BPA, &
          n_visc_its, n_Axb_its, &
          BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
        call set_ice_velocities_to_BPA_results( mesh, vel, ice%BPA)

      case ('hybrid DIVA/BPA')
        ! Calculate velocities according to the hybrid DIVA/BPA

        call solve_hybrid_DIVA_BPA( mesh, ice, geom, bed_roughness, ice%hybrid, region_name, &
          n_visc_its, n_Axb_its, &
          BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
        call set_ice_velocities_to_hybrid_DIVA_BPA_results( mesh, vel, ice%hybrid)

    end select

    ! Calculate all secondary ice velocities (surface, base, vertical average)
    call calc_secondary_velocities( mesh, vel)

    ! Calculate vertical velocities
    call calc_vertical_velocities( vel, mesh, ice, geom, BMB)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_stress_balance

  subroutine calc_secondary_velocities( mesh, vel)
    !< Calculate all secondary ice velocities (surface, base, vertical average)

    ! In/output variables:
    type(type_mesh),                      intent(in   ) :: mesh
    class(atype_ice_velocity_model_data), intent(inout) :: vel

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'calc_secondary_velocities'
    integer                           :: vi,ti
    real(dp), dimension(mesh%nz)      :: u_prof, v_prof

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2

      ! Surface
      vel%u_surf_b(    ti) = vel%u_3D_b( ti,1)
      vel%v_surf_b(    ti) = vel%v_3D_b( ti,1)
      vel%uabs_surf_b( ti) = SQRT( vel%u_surf_b( ti)**2 + vel%v_surf_b( ti)**2)

      ! Base
      vel%u_base_b(    ti) = vel%u_3D_b( ti,C%nz)
      vel%v_base_b(    ti) = vel%v_3D_b( ti,C%nz)
      vel%uabs_base_b( ti) = SQRT( vel%u_base_b( ti)**2 + vel%v_base_b( ti)**2)

      ! Vertical average
      u_prof = vel%u_3D_b( ti,:)
      v_prof = vel%v_3D_b( ti,:)
      vel%u_vav_b( ti) = vertical_average( mesh%zeta, u_prof)
      vel%v_vav_b( ti) = vertical_average( mesh%zeta, v_prof)
      vel%uabs_vav_b( ti) = SQRT( vel%u_vav_b( ti)**2 + vel%v_vav_b( ti)**2)

    end do

    ! == Calculate velocities on the a-grid (needed to calculate the vertical velocity w, and for writing to output)

    ! 3-D
    call map_b_a_3D( mesh, vel%u_3D_b  , vel%u_3D  )
    call map_b_a_3D( mesh, vel%v_3D_b  , vel%v_3D  )

    ! Surface
    call map_b_a_2D( mesh, vel%u_surf_b, vel%u_surf)
    call map_b_a_2D( mesh, vel%v_surf_b, vel%v_surf)

    ! Base
    call map_b_a_2D( mesh, vel%u_base_b, vel%u_base)
    call map_b_a_2D( mesh, vel%v_base_b, vel%v_base)

    ! Vertical average
    call map_b_a_2D( mesh, vel%u_vav_b , vel%u_vav )
    call map_b_a_2D( mesh, vel%v_vav_b , vel%v_vav )

    ! Absolute
    do vi = mesh%vi1, mesh%vi2
      vel%uabs_surf( vi) = sqrt( vel%u_surf( vi)**2 + vel%v_surf( vi)**2)
      vel%uabs_base( vi) = sqrt( vel%u_base( vi)**2 + vel%v_base( vi)**2)
      vel%uabs_vav(  vi) = sqrt( vel%u_vav(  vi)**2 + vel%v_vav(  vi)**2)
    end do

    call calc_u_vav_perp( mesh, vel)

    ! Slide/shear ratio
    do vi = mesh%vi1, mesh%vi2
      vel%R_shear( vi) = (vel%uabs_base( vi) + 0.1_dp) / (vel%uabs_surf( vi) + 0.1_dp)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_secondary_velocities

  subroutine calc_u_vav_perp( mesh, vel)
    !< Calculate the vertically averaged ice velocity component
    !< perpendicular to the shared Voronoi cell boundaries

    ! In/output variables:
    type(type_mesh),                      intent(in   ) :: mesh
    class(atype_ice_velocity_model_data), intent(inout) :: vel

    ! Local variables:
    character(len=*), parameter            :: routine_name = 'calc_u_vav_perp'
    real(dp), dimension(mesh%ei1:mesh%ei2) :: u_vav_c, v_vav_c
    real(dp), dimension(mesh%nE)           :: u_vav_c_tot, v_vav_c_tot
    integer                                :: vi, ci, ei, vj

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate vertically averaged ice velocities on the edges
    call map_velocities_from_b_to_c_2D( mesh, vel%u_vav_b, vel%v_vav_b, u_vav_c, v_vav_c)
    call gather_to_all( u_vav_c, u_vav_c_tot)
    call gather_to_all( v_vav_c, v_vav_c_tot)

    do vi = mesh%vi1, mesh%vi2
      do ci = 1, mesh%nC( vi)

        ! Connection ci from vertex vi leads through edge ei to vertex vj
        ei = mesh%VE( vi,ci)

        ! Calculate vertically averaged ice velocity component perpendicular to this shared Voronoi cell boundary section
        vel%u_vav_perp( vi, ci) = &
          u_vav_c_tot( ei) * mesh%D_x( vi, ci)/mesh%D( vi, ci) + &
          v_vav_c_tot( ei) * mesh%D_y( vi, ci)/mesh%D( vi, ci)

      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_u_vav_perp

  subroutine remap_velocity_solver( vel, momentum_balance_solver, mesh_old, mesh_new, ice)
    !< Remap the velocity solver for the chosen stress balance approximation

    ! In/output variables:
    class(atype_ice_velocity_model),      intent(inout) :: vel
    class(atype_momentum_balance_solver), intent(inout) :: momentum_balance_solver
    type(type_mesh),                      intent(in   ) :: mesh_old
    type(type_mesh),                      intent(in   ) :: mesh_new
    class(atype_ice_model_data),          intent(inout) :: ice

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_velocity_solver'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_stress_balance_approximation)

      case default
        call crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')

      case ('none')
      ! No need to do anything

      case ('SIA')

        call momentum_balance_solver%remap( mesh_old, mesh_new)

        select type (SIA => momentum_balance_solver)
        class default
          call crash('invalid momentum_balance_solver class')
        class is (type_momentum_balance_solver_SIA)
          call set_ice_velocities_to_SIA_results( mesh_new, vel, SIA)
        end select

      case ('SSA')

        call momentum_balance_solver%remap( mesh_old, mesh_new)

        select type (SSA => momentum_balance_solver)
        class default
          call crash('invalid momentum_balance_solver class')
        class is (type_momentum_balance_solver_SSA)
          call set_ice_velocities_to_SSA_results( mesh_new, vel, SSA)
        end select

      case ('SIA/SSA')

        call momentum_balance_solver%remap( mesh_old, mesh_new)

        select type (SIASSA => momentum_balance_solver)
        class default
          call crash('invalid momentum_balance_solver class')
        class is (type_momentum_balance_solver_SIASSA)
          call set_ice_velocities_to_SIASSA_results( mesh_new, vel, SIASSA)
        end select

      case ('DIVA')

        call remap_DIVA_solver( mesh_old, mesh_new, ice%DIVA)
        call set_ice_velocities_to_DIVA_results( mesh_new, ice, vel, ice%DIVA)

      case ('BPA')

        call remap_BPA_solver(  mesh_old, mesh_new, ice%BPA)
        call set_ice_velocities_to_BPA_results( mesh_new, vel, ice%BPA)

      case ('hybrid DIVA/BPA')

        call remap_hybrid_DIVA_BPA_solver(  mesh_old, mesh_new, ice%hybrid)
        call set_ice_velocities_to_hybrid_DIVA_BPA_results( mesh_new, vel, ice%hybrid)

    end select

    ! Calculate all secondary ice velocities (surface, base, vertical average)
    call calc_secondary_velocities( mesh_new, vel)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_velocity_solver

  ! == Set applied ice model velocities to stress balance results

  subroutine set_ice_velocities_to_SIA_results( mesh, vel, SIA)
    !< Set applied ice model velocities to SIA results

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    class(atype_ice_velocity_model_data),    intent(inout) :: vel
    class(type_momentum_balance_solver_SIA), intent(inout) :: SIA

    ! Local variables:
    character(len=*), parameter :: routine_name = 'set_ice_velocities_to_SIA_results'
    integer                     :: ti,vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Velocities
    do ti = mesh%ti1, mesh%ti2
      vel%u_3D_b( ti,:) = SIA%solver%u_3D_b( ti,:)
      vel%v_3D_b( ti,:) = SIA%solver%v_3D_b( ti,:)
    end do

    ! Strain rates
    do vi = mesh%vi1, mesh%vi2
      vel%du_dz_3D( vi,:) = SIA%solver%du_dz_3D( vi,:)
      vel%dv_dz_3D( vi,:) = SIA%solver%dv_dz_3D( vi,:)
    end do

    ! In the SIA, horizontal gradients of u,v, and all gradients of w, are neglected
    vel%du_dx_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    vel%du_dy_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    vel%dv_dx_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    vel%dv_dy_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    vel%dw_dx_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    vel%dw_dy_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    ! vel%dw_dz_3D = 0._dp ! Because we now always calculate dw/dz in calc_vertical_velocities

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_ice_velocities_to_SIA_results

  subroutine set_ice_velocities_to_SSA_results( mesh, vel, SSA)
    !< Set applied ice model velocities to SSA results

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    class(atype_ice_velocity_model_data),    intent(inout) :: vel
    class(type_momentum_balance_solver_SSA), intent(inout) :: SSA

    ! Local variables:
    character(len=*), parameter :: routine_name = 'set_ice_velocities_to_SSA_results'
    integer                     :: ti,vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Velocities
    do ti = mesh%ti1, mesh%ti2
      vel%u_3D_b( ti,:) = SSA%solver%u_b( ti)
      vel%v_3D_b( ti,:) = SSA%solver%v_b( ti)
    end do

    ! Strain rates
    do vi = mesh%vi1, mesh%vi2
      vel%du_dx_3D( vi,:) = SSA%solver%du_dx_a( vi)
      vel%du_dy_3D( vi,:) = SSA%solver%du_dy_a( vi)
      vel%dv_dx_3D( vi,:) = SSA%solver%dv_dx_a( vi)
      vel%dv_dy_3D( vi,:) = SSA%solver%dv_dy_a( vi)
    end do

    ! In the SSA, vertical gradients of u,v, and all gradients of w, are neglected
    vel%du_dz_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    vel%dv_dz_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    vel%dw_dx_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    vel%dw_dy_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    ! vel%dw_dz_3D = 0._dp ! Because we now always calculate dw/dz in calc_vertical_velocities

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_ice_velocities_to_SSA_results

  subroutine set_ice_velocities_to_SIASSA_results( mesh, vel, SIASSA)
    !< Set applied ice model velocities to SIA/SSA results

    ! In/output variables:
    type(type_mesh),                            intent(in   ) :: mesh
    class(atype_ice_velocity_model_data),       intent(inout) :: vel
    class(type_momentum_balance_solver_SIASSA), intent(inout) :: SIASSA

    ! Local variables:
    character(len=*), parameter :: routine_name = 'set_ice_velocities_to_SIASSA_results'
    integer                     :: ti,vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Velocities
    do ti = mesh%ti1, mesh%ti2
      vel%u_3D_b( ti,:) = SIASSA%solver_SIA%u_3D_b( ti,:) + SIASSA%solver_SSA%u_b( ti)
      vel%v_3D_b( ti,:) = SIASSA%solver_SIA%v_3D_b( ti,:) + SIASSA%solver_SSA%v_b( ti)
    end do

    ! Strain rates
    do vi = mesh%vi1, mesh%vi2
      vel%du_dx_3D( vi,:) = SIASSA%solver_SSA%du_dx_a( vi  )
      vel%du_dy_3D( vi,:) = SIASSA%solver_SSA%du_dy_a( vi  )
      vel%du_dz_3D( vi,:) = SIASSA%solver_SIA%du_dz_3D  ( vi,:)
      vel%dv_dx_3D( vi,:) = SIASSA%solver_SSA%dv_dx_a( vi  )
      vel%dv_dy_3D( vi,:) = SIASSA%solver_SSA%dv_dy_a( vi  )
      vel%dv_dz_3D( vi,:) = SIASSA%solver_SIA%dv_dz_3D  ( vi,:)
    end do

    ! In the SIA/SSA, all gradients of w are neglected
    vel%dw_dx_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    vel%dw_dy_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    ! vel%dw_dz_3D = 0._dp ! Because we now always calculate dw/dz in calc_vertical_velocities

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_ice_velocities_to_SIASSA_results

  subroutine set_ice_velocities_to_DIVA_results( mesh, ice, vel, DIVA)
    !< Set applied ice model velocities to DIVA results

    ! In/output variables:
    type(type_mesh),                      intent(in   ) :: mesh
    class(atype_ice_model_data),          intent(inout) :: ice
    class(atype_ice_velocity_model_data), intent(inout) :: vel
    type(type_ice_velocity_solver_DIVA),  intent(in   ) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_ice_velocities_to_DIVA_results'
    integer                        :: ti,vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Velocities
    do ti = mesh%ti1, mesh%ti2
      vel%u_3D_b( ti,:) = DIVA%u_3D_b( ti,:)
      vel%v_3D_b( ti,:) = DIVA%v_3D_b( ti,:)
    end do

    ! Strain rates
    do vi = mesh%vi1, mesh%vi2
      vel%du_dx_3D( vi,:) = DIVA%du_dx_a(    vi  )
      vel%du_dy_3D( vi,:) = DIVA%du_dy_a(    vi  )
      vel%du_dz_3D( vi,:) = DIVA%du_dz_3D_a( vi,:)
      vel%dv_dx_3D( vi,:) = DIVA%dv_dx_a(    vi  )
      vel%dv_dy_3D( vi,:) = DIVA%dv_dy_a(    vi  )
      vel%dv_dz_3D( vi,:) = DIVA%dv_dz_3D_a( vi,:)
    end do

    ! In the DIVA, gradients of w are neglected
    vel%dw_dx_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    vel%dw_dy_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    ! vel%dw_dz_3D = 0._dp ! Because we now always calculate dw/dz in calc_vertical_velocities

    ! Stresses
    do ti = mesh%ti1, mesh%ti2
      ice%basal_shear_stress( ti) = hypot( ice%DIVA%tau_bx_b( ti), ice%DIVA%tau_by_b( ti))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_ice_velocities_to_DIVA_results

  subroutine set_ice_velocities_to_BPA_results( mesh, vel, BPA)
    ! Set applied ice model velocities and strain rates to BPA results

    ! In/output variables:
    type(type_mesh),                      intent(in   ) :: mesh
    class(atype_ice_velocity_model_data), intent(inout) :: vel
    type(type_ice_velocity_solver_BPA),   intent(in   ) :: BPA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_ice_velocities_to_BPA_results'
    integer                        :: ti,vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Velocities
    do ti = mesh%ti1, mesh%ti2
      vel%u_3D_b( ti,:) = BPA%u_bk( ti,:)
      vel%v_3D_b( ti,:) = BPA%v_bk( ti,:)
    end do

    ! Strain rates
    do vi = mesh%vi1, mesh%vi2
      vel%du_dx_3D( vi,:) = BPA%du_dx_ak( vi,:)
      vel%du_dy_3D( vi,:) = BPA%du_dy_ak( vi,:)
      vel%du_dz_3D( vi,:) = BPA%du_dz_ak( vi,:)
      vel%dv_dx_3D( vi,:) = BPA%dv_dx_ak( vi,:)
      vel%dv_dy_3D( vi,:) = BPA%dv_dy_ak( vi,:)
      vel%dv_dz_3D( vi,:) = BPA%dv_dz_ak( vi,:)
    end do

    ! In the BPA, gradients of w are neglected
    vel%dw_dx_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    vel%dw_dy_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    ! vel%dw_dz_3D = 0._dp ! Because we now always calculate dw/dz in calc_vertical_velocities

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_ice_velocities_to_BPA_results

  subroutine set_ice_velocities_to_hybrid_DIVA_BPA_results( mesh, vel, hybrid)
    !< Set applied ice model velocities and strain rates to hybrid DIVA/BPA results

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    class(atype_ice_velocity_model_data),  intent(inout) :: vel
    type(type_ice_velocity_solver_hybrid), intent(in   ) :: hybrid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_ice_velocities_to_hybrid_DIVA_BPA_results'
    integer                        :: ti,vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Velocities
    do ti = mesh%ti1, mesh%ti2
      vel%u_3D_b( ti,:) = hybrid%u_bk( ti,:)
      vel%v_3D_b( ti,:) = hybrid%v_bk( ti,:)
    end do

    ! Strain rates
    do vi = mesh%vi1, mesh%vi2
      vel%du_dx_3D( vi,:) = hybrid%BPA%du_dx_ak( vi,:)
      vel%du_dy_3D( vi,:) = hybrid%BPA%du_dy_ak( vi,:)
      vel%du_dz_3D( vi,:) = hybrid%BPA%du_dz_ak( vi,:)
      vel%dv_dx_3D( vi,:) = hybrid%BPA%dv_dx_ak( vi,:)
      vel%dv_dy_3D( vi,:) = hybrid%BPA%dv_dy_ak( vi,:)
      vel%dv_dz_3D( vi,:) = hybrid%BPA%dv_dz_ak( vi,:)
    end do

    ! In the hybrid DIVA/BPA, gradients of w are neglected
    vel%dw_dx_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    vel%dw_dy_3D( mesh%vi1:mesh%vi2,:) = 0._dp
    ! vel%dw_dz_3D = 0._dp ! Because we now always calculate dw/dz in calc_vertical_velocities

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_ice_velocities_to_hybrid_DIVA_BPA_results

  ! == Restart NetCDF files

  subroutine write_to_restart_file_ice_velocity( mesh, ice, momentum_balance_solver, time)
    !< Write to the restart NetCDF file for the ice velocity solver

    ! In/output variables:
    type(type_mesh),                      intent(in   ) :: mesh
    class(atype_ice_model_data),          intent(in   ) :: ice
    class(atype_momentum_balance_solver), intent(in   ) :: momentum_balance_solver
    real(dp),                             intent(in   ) :: time

    ! Local variables:
    character(len=*), parameter :: routine_name = 'write_to_restart_file_ice_velocity'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_stress_balance_approximation)
      case default
        call crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
      case ('none')
        ! No need to do anything
      case ('SIA')
        ! The SIA doesn't have a restart file
      case ('SSA')

        select type (SSA => momentum_balance_solver)
        class default
          call crash('invalid momentum_balance_solver class')
        class is (type_momentum_balance_solver_SSA)
          call SSA%solver%write_to_restart_file_SSA( time)
        end select

      case ('SIA/SSA')

        select type (SIASSA => momentum_balance_solver)
        class default
          call crash('invalid momentum_balance_solver class')
        class is (type_momentum_balance_solver_SIASSA)
          call SIASSA%solver_SSA%write_to_restart_file_SSA( time)
        end select

      case ('DIVA')
        call write_to_restart_file_DIVA( mesh, ice%DIVA, time)
      case ('BPA')
        call write_to_restart_file_BPA( mesh, ice%BPA, time)
      case ('hybrid DIVA/BPA')
        call warning('the hybrid DIVA/BPA does not have a restart file yet!')
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file_ice_velocity

  subroutine create_restart_file_ice_velocity( mesh, ice, momentum_balance_solver)
    !< Create a restart NetCDF file for the ice velocity solver

    ! In/output variables:
    type(type_mesh),                      intent(in   ) :: mesh
    class(atype_ice_model_data),          intent(inout) :: ice
    class(atype_momentum_balance_solver), intent(inout) :: momentum_balance_solver

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_restart_file_ice_velocity'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_stress_balance_approximation)
    case default
      call crash('unknown choice_stress_balance_approximation "' // TRIM( C%choice_stress_balance_approximation) // '"!')
    case ('none')
      ! No need to do anything
    case ('SIA')
      ! The SIA doesn't have a restart file
    case ('SSA')

      select type (SSA => momentum_balance_solver)
      class default
        call crash('invalid momentum_balance_solver class')
      class is (type_momentum_balance_solver_SSA)
        call SSA%solver%create_restart_file_SSA()
      end select

    case ('SIA/SSA')

      select type (SIASSA => momentum_balance_solver)
      class default
        call crash('invalid momentum_balance_solver class')
      class is (type_momentum_balance_solver_SIASSA)
        call SIASSA%solver_SSA%create_restart_file_SSA()
      end select

    case ('DIVA')
      call create_restart_file_DIVA( mesh, ice%DIVA)
    case ('BPA')
      call create_restart_file_BPA( mesh, ice%BPA)
    case ('hybrid DIVA/BPA')
      call warning('the hybrid DIVA/BPA does not have a restart file yet!')
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_ice_velocity

end module ice_velocities_main
