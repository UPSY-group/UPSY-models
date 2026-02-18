module basal_hydrology_main

  ! Contains all the different basal hydrology models.

  use mpi_basic, only: par
  use precisions, only: dp
  use erf_mod, only: error_function
  use call_stack_and_comp_time_tracking, only: crash, init_routine, finalise_routine
  use model_configuration, only: C
  use parameters, only: grav, ice_density, pi, seawater_density
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model

  implicit none

contains

  subroutine run_basal_hydrology_model( mesh, ice)

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_basal_hydrology_model'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate pore water pressure using the chosen basal hydrology model
    ! ====================================================================

    select case (C%choice_basal_hydrology_model)
    case default
      call crash('unknown choice_basal_hydrology_model "' // trim( C%choice_basal_hydrology_model) // '"')
    case ('none')
      call calc_pore_water_pressure_none( mesh, ice)
      ! Calculate overburden and effective pressure
      do vi = mesh%vi1, mesh%vi2
        ice%overburden_pressure( vi) = ice_density * grav * ice%Hi_eff( vi)
        ice%effective_pressure(  vi) = max( 0._dp, ice%overburden_pressure( vi) - ice%pore_water_pressure( vi))
      end do
    case ('Martin2011')
      call calc_pore_water_pressure_Martin2011( mesh, ice)
      ! Calculate overburden and effective pressure
      do vi = mesh%vi1, mesh%vi2
        ice%overburden_pressure( vi) = ice_density * grav * ice%Hi_eff( vi)
        ice%effective_pressure(  vi) = max( 0._dp, ice%overburden_pressure( vi) - ice%pore_water_pressure( vi))
      end do
    case ('Leguy2014')
      call calc_pore_water_pressure_none( mesh, ice)
      call calc_effective_pressure_Leguy2014(mesh, ice)
    case ('error_function_Martin2011')
      call calc_pore_water_pressure_Martin2011( mesh, ice)   ! we need that for the maximum inland effective pressure
      call calc_effective_pressure_error_function_M11(mesh, ice)
    case ('error_function_constant')
      call calc_effective_pressure_error_function_constant(mesh, ice)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_basal_hydrology_model

  subroutine calc_pore_water_pressure_none( mesh, ice)

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_pore_water_pressure_none'
    integer                        :: vi
    real(dp)                       :: weight_gr

    ! Add routine to path
    call init_routine( routine_name)

    ! Compute pore water pressure based on the pore water fraction as
    ! the fraction of the overburden pressure supported by basal water
    do vi = mesh%vi1, mesh%vi2
      ice%pore_water_pressure( vi) = ice%pore_water_fraction(vi) * ice_density * grav * ice%Hi_eff( vi)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_pore_water_pressure_none

  subroutine calc_pore_water_pressure_Martin2011( mesh, ice)
    ! Calculate pore water pressure according to the parameterisation from Martin et al. (2011)

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_pore_water_pressure_Martin2011'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      ! Pore water pressure scaling factor (Martin et al., 2011, Eq. 12)
      ice%pore_water_fraction( vi) = min( 1._dp, max( 0._dp, &
        1._dp - (ice%Hb( vi) - ice%SL( vi) - C%Martin2011_hydro_Hb_min) / (C%Martin2011_hydro_Hb_max - C%Martin2011_hydro_Hb_min) ))

      ! Pore water pressure (Martin et al., 2011, Eq. 11)
      ice%pore_water_pressure( vi) = 0.96_dp * ice_density * grav * ice%Hi_eff( vi) * ice%pore_water_fraction( vi)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_pore_water_pressure_Martin2011

  subroutine calc_effective_pressure_error_function_M11( mesh, ice)
    ! Calculate pore water pressure as an error function

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_effective_pressure_error_function_M11'
    integer                        :: vi
    real(dp)                       :: effective_pressure_max

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
        ice%overburden_pressure( vi) = ice_density * grav * ice%Hi_eff( vi)
        effective_pressure_max = max( 0._dp, ice%overburden_pressure( vi) - ice%pore_water_pressure( vi))

        if (effective_pressure_max == 0._dp) then
          ice%effective_pressure( vi)  = 0.0_dp
        else
          ice%effective_pressure(  vi) = error_function(ice%overburden_pressure( vi)*sqrt(pi)/2._dp/effective_pressure_max)*effective_pressure_max
        end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_pressure_error_function_M11


  subroutine calc_effective_pressure_Leguy2014(mesh, ice)
    ! Calculate effective pressure based on Leguy et al. (2014)

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_effective_pressure_Leguy2014'
    integer                        :: vi
    real(dp)                       :: Hi_f

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      ! check if ice is grounded first
      if (ice%mask_grounded_ice( vi) .eqv. .TRUE.) then
         ! prevent division by zero
        if (ice%Hi_eff( vi) == 0._dp) then
          ice%effective_pressure( vi) = 0.0_dp
        else
          ice%overburden_pressure( vi) = ice_density * grav * ice%Hi_eff( vi)
          Hi_f = max(0._dp, - seawater_density/ice_density * ice%Hb( vi))
          ! if (Hi_f == 0._dp) then
          !   ice%effective_pressure( vi) = 0.0_dp
          ! else
          ice%effective_pressure( vi) = ice%overburden_pressure( vi) * ((1 - Hi_f/ice%Hi_eff( vi)) ** C%Leguy2014_hydro_connect_exponent)
          ! end if
        end if
      else
        ice%effective_pressure( vi) = 0.0_dp
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)
  end subroutine calc_effective_pressure_Leguy2014

  subroutine calc_effective_pressure_error_function_constant( mesh, ice)
    ! Calculate pore water pressure as an error function

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_effective_pressure_error_function_constant'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
        ice%overburden_pressure( vi) = ice_density * grav * ice%Hi_eff( vi)
        ice%effective_pressure(  vi) = error_function(ice%overburden_pressure( vi)*sqrt(pi)/2._dp/C%error_function_max_effective_pressure)*C%error_function_max_effective_pressure
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_pressure_error_function_constant

end module basal_hydrology_main
