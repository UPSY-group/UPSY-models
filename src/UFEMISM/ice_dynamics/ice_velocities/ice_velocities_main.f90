module ice_velocities_main

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use momentum_balance_solver_SIA, only: type_momentum_balance_solver_SIA
  use momentum_balance_solver_SSA, only: type_momentum_balance_solver_SSA
  use momentum_balance_solver_SIASSA, only: type_momentum_balance_solver_SIASSA
  use momentum_balance_solver_DIVA, only: type_momentum_balance_solver_DIVA
  use momentum_balance_solver_BPA, only: type_momentum_balance_solver_BPA
  use momentum_balance_solver_hybrid_DIVA_BPA, only: type_momentum_balance_solver_hybrid_DIVA_BPA

  implicit none

contains

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
