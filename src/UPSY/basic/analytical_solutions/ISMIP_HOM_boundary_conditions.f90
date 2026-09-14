module ISMIP_HOM_boundary_conditions

  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking , only:  crash

  implicit none

  private

  public :: u_BC_ISMIP_HOM

contains

  function u_BC_ISMIP_HOM() result( u)
    !< Since periodic BCs are an absolute pain in the rear, instead
    !< we just prescribe the value shown at [normalized x] = 0.5
    !< in the figures from Pattyn et al. (2008)

    real(dp) :: u

    u = 0._dp

    select case (C%choice_refgeo_init_idealised)
    case default
      call crash('invalid choice_refgeo_init_idealised ' // trim( C%choice_refgeo_init_idealised))

    case ('ISMIP-HOM_A')

      if     (C%refgeo_idealised_ISMIP_HOM_L == 160e3_dp) then
        u = 20._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 80e3_dp) then
        u = 22._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 40e3_dp) then
        u = 27._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 20e3_dp) then
        u = 25._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 10e3_dp) then
        u = 22._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 5e3_dp) then
        u = 14._dp
      else
        call crash('invalid value for refgeo_idealised_ISMIP_HOM_L')
      end if

    case ('ISMIP-HOM_B')

      if     (C%refgeo_idealised_ISMIP_HOM_L == 160e3_dp) then
        u = 21._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 80e3_dp) then
        u = 25._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 40e3_dp) then
        u = 30._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 20e3_dp) then
        u = 27._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 10e3_dp) then
        u = 20._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 5e3_dp) then
        u = 9._dp
      else
        call crash('invalid value for refgeo_idealised_ISMIP_HOM_L')
      end if

    case ('ISMIP-HOM_C')

      if     (C%refgeo_idealised_ISMIP_HOM_L == 160e3_dp) then
        u = 20._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 80e3_dp) then
        u = 18._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 40e3_dp) then
        u = 17._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 20e3_dp) then
        u = 16.5_dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 10e3_dp) then
        u = 15.5_dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 5e3_dp) then
        u = 12._dp
      else
        call crash('invalid value for refgeo_idealised_ISMIP_HOM_L')
      end if

    case ('ISMIP-HOM_D')

      if     (C%refgeo_idealised_ISMIP_HOM_L == 160e3_dp) then
        u = 20._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 80e3_dp) then
        u = 20._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 40e3_dp) then
        u = 20._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 20e3_dp) then
        u = 18._dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 10e3_dp) then
        u = 16.5_dp
      elseif (C%refgeo_idealised_ISMIP_HOM_L == 5e3_dp) then
        u = 13._dp
      else
        call crash('invalid value for refgeo_idealised_ISMIP_HOM_L')
      end if

    end select

  end function u_BC_ISMIP_HOM

end module ISMIP_HOM_boundary_conditions