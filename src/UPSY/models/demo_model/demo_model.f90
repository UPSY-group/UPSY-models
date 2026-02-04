module demo_model

  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use demo_model_basic, only: atype_demo_model
  use demo_model_a, only: type_demo_model_a
  use demo_model_b, only: type_demo_model_b

  implicit none

  private

  public :: atype_demo_model, create_demo_model

contains

  subroutine create_demo_model( demo, choice_demo_model)
    !< Allocate a concrete implementation of a demo_model

    ! In/output variables:
    class(atype_demo_model), allocatable, intent(inout) :: demo
    character(len=*),                     intent(in   ) :: choice_demo_model

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_demo_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    select case (choice_demo_model)
    case default
      call crash('invalid choice_demo_model "' // trim( choice_demo_model) // '"')
    case ('demo_a')
      allocate( type_demo_model_a :: demo)
    case ('demo_b')
      allocate( type_demo_model_b :: demo)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_demo_model

end module demo_model