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
  use ice_model_data, only: atype_ice_model_data
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
  use momentum_balance_solver_DIVA, only: type_momentum_balance_solver_DIVA
  use momentum_balance_solver_BPA, only: type_momentum_balance_solver_BPA
  use momentum_balance_solver_hybrid_DIVA_BPA, only: type_momentum_balance_solver_hybrid_DIVA_BPA

  implicit none

contains

  ! == The main routines, to be called from the ice dynamics module

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

    call momentum_balance_solver%run( ice, geom, bed_roughness, &
      BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)
    call momentum_balance_solver%set_velocities_to_solver_results( ice, vel)
    call calc_secondary_velocities( mesh, vel)
    call calc_vertical_velocities( vel, mesh, ice, geom, BMB)

    n_visc_its = momentum_balance_solver%n_visc_its
    n_Axb_its  = momentum_balance_solver%n_Axb_its

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

    call momentum_balance_solver%remap( mesh_old, mesh_new)
    call momentum_balance_solver%set_velocities_to_solver_results( ice, vel)
    call calc_secondary_velocities( mesh_new, vel)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_velocity_solver

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
          call SSA%write_to_restart_file_SSA( time)
        end select

      case ('SIA/SSA')

        select type (SIASSA => momentum_balance_solver)
        class default
          call crash('invalid momentum_balance_solver class')
        class is (type_momentum_balance_solver_SIASSA)
          call SIASSA%SSA%write_to_restart_file_SSA( time)
        end select

      case ('DIVA')

        select type (DIVA => momentum_balance_solver)
        class default
          call crash('invalid momentum_balance_solver class')
        class is (type_momentum_balance_solver_DIVA)
          call DIVA%write_to_restart_file_DIVA( time)
        end select

      case ('BPA')

        select type (BPA => momentum_balance_solver)
        class default
          call crash('invalid momentum_balance_solver class')
        class is (type_momentum_balance_solver_BPA)
          call BPA%write_to_restart_file_BPA( time)
        end select

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
        call SSA%create_restart_file_SSA()
      end select

    case ('SIA/SSA')

      select type (SIASSA => momentum_balance_solver)
      class default
        call crash('invalid momentum_balance_solver class')
      class is (type_momentum_balance_solver_SIASSA)
        call SIASSA%SSA%create_restart_file_SSA()
      end select

    case ('DIVA')

      select type (DIVA => momentum_balance_solver)
      class default
        call crash('invalid momentum_balance_solver class')
      class is (type_momentum_balance_solver_DIVA)
        call DIVA%create_restart_file_DIVA()
      end select

    case ('BPA')

        select type (BPA => momentum_balance_solver)
        class default
          call crash('invalid momentum_balance_solver class')
        class is (type_momentum_balance_solver_BPA)
          call BPA%create_restart_file_BPA()
        end select

    case ('hybrid DIVA/BPA')
      call warning('the hybrid DIVA/BPA does not have a restart file yet!')
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_ice_velocity

end module ice_velocities_main
