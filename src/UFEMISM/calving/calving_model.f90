module calving_model

  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use calving_model_basic, only: atype_calving_model
  use calving_position, only: type_calving_model_position
  ! use calving_rate, only: type_calving_model_rate
  use model_configuration, only: C

  implicit none

  private

  public :: atype_calving_model, create_calving_model

contains

  subroutine create_calving_model( calving)
    !< Allocate a concrete implementation of a calving_model

    ! In/output variables:
    class(atype_calving_model), allocatable, intent(inout) :: calving
    ! character(len=*),                        intent(in   ) :: choice_calving_model

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_calving_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    select case (trim(C%choice_calving_law))
    case default
      call crash('invalid choice_calving_law "' // trim( C%choice_calving_law) // '"')
    case ('threshold_thickness_front')
      allocate( type_calving_model_position :: calving)
    case ('threshold_thickness_all')
      allocate( type_calving_model_position :: calving)
    ! case ('rate_law')
    !   allocate( type_calving_model_rate :: calving)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_calving_model

end module calving_model