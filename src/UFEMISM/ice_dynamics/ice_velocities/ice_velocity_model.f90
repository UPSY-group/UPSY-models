module ice_velocity_model

  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use ice_velocity_model_basic, only: atype_ice_velocity_model, type_ice_velocity_model_clean

  implicit none

  private

  public :: atype_ice_velocity_model, create_ice_velocity_model

contains

  subroutine create_ice_velocity_model( vel, choice_stress_balance_approximation)
    !< Allocate a concrete implementation of an ice_velocity_model

    ! In/output variables:
    class(atype_ice_velocity_model), allocatable, intent(inout) :: vel
    character(len=*),                             intent(in   ) :: choice_stress_balance_approximation

    ! Local variables:
    character(len=*), parameter :: routine_name = 'create_ice_velocity_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    select case (choice_stress_balance_approximation)
    case default
      call crash('unknown choice_stress_balance_approximation "' // trim( choice_stress_balance_approximation) // '"!')
    case ('none')
      allocate( type_ice_velocity_model_clean :: vel)
    case ('SIA')
      allocate( type_ice_velocity_model_clean :: vel)
    case ('SSA')
      allocate( type_ice_velocity_model_clean :: vel)
    case ('SIA/SSA')
      allocate( type_ice_velocity_model_clean :: vel)
    case ('DIVA')
      allocate( type_ice_velocity_model_clean :: vel)
    case ('BPA')
      allocate( type_ice_velocity_model_clean :: vel)
    case ('hybrid DIVA/BPA')
      allocate( type_ice_velocity_model_clean :: vel)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_ice_velocity_model

end module ice_velocity_model