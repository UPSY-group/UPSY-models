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

  subroutine remap_velocity_solver( vel, momentum_balance_solver, mesh_old, mesh_new, ice, geom, BMB)
    !< Remap the velocity solver for the chosen stress balance approximation

    ! In/output variables:
    class(atype_ice_velocity_model),                intent(inout) :: vel
    class(atype_momentum_balance_solver),           intent(inout) :: momentum_balance_solver
    type(type_mesh),                                intent(in   ) :: mesh_old
    type(type_mesh),                                intent(in   ) :: mesh_new
    class(atype_ice_model_data),                    intent(inout) :: ice
    class(atype_ice_geometry_model_data),           intent(in   ) :: geom
    real(dp), dimension(mesh_new%vi1:mesh_new%vi2), intent(in   ) :: BMB

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_velocity_solver'

    ! Add routine to path
    call init_routine( routine_name)

    call momentum_balance_solver%remap( mesh_old, mesh_new, ice, geom, vel, BMB)

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
