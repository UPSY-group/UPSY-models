module LMB_inverted

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use model_configuration, only: C

  implicit none

  private

  public :: calc_inverted_LMB

contains

  subroutine calc_inverted_LMB( mesh, divQ, fraction_margin, SMB, BMB, dHi_dt_target, LMB, region_name)
    ! In/output variables:
    type(type_mesh),                        intent(in   )           :: mesh                  ! [-]       The model mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: divQ                  ! [m yr^-1] Horizontal ice flux divergence
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: fraction_margin       ! [0-1]     Sub-grid ice-filled fraction
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: SMB                   ! [m yr^-1] Surface mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: BMB                   ! [m yr^-1] Basal   mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   )           :: dHi_dt_target         ! [m yr^-1] Target ice thickness rate of change
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout)           :: LMB                   ! [m yr^-1] Lateral mass balance
    character(len=3),                       intent(in   )           :: region_name

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'calc_inverted_LMB'
    character(len=256)              :: choice_LMB_model
    integer                         :: vi

    ! Determine whether inversion is needed
    select case (region_name)
      case ('NAM')
        choice_LMB_model = C%choice_LMB_model_NAM
      case ('EAS')
        choice_LMB_model = C%choice_LMB_model_EAS
      case ('GRL')
        choice_LMB_model = C%choice_LMB_model_GRL
      case ('ANT')
        choice_LMB_model = C%choice_LMB_model_ANT
      case default
        call crash('unknown region_name "' // region_name // '"')
    end select

    if (choice_LMB_model == 'inverted') then
      ! Compute LMB based on the current mass balance terms
      ! This routine should be called from the conversion_of_mass scheme,
      ! after divQ is defined, and before dHi_dt is computed
      ! It ensures that dHi_dt = 0 on margin cells

      do vi = mesh%vi1, mesh%vi2

        if (fraction_margin( vi) > 0.0_dp .and. fraction_margin( vi) < 1.0_dp) then

          LMB( vi) = min(0._dp, divQ( vi) & 
            - fraction_margin( vi) * ( &
              SMB( vi) + BMB( vi) - dHi_dt_target( vi)))

        else

          LMB( vi) = 0.0_dp

        end if

      end do

    end if

  end subroutine calc_inverted_LMB

end module LMB_inverted
