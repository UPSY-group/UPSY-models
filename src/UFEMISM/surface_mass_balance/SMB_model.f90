module SMB_model

  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use SMB_model_basic, only: atype_SMB_model
  use SMB_idealised, only: type_SMB_model_idealised
  use SMB_prescribed, only: type_SMB_model_prescribed
  use SMB_reconstructed, only: type_SMB_model_reconstructed
  use SMB_snapshot_plus_anomalies, only: type_SMB_model_snp_p_anml
  use SMB_IMAU_ITM, only: type_SMB_model_IMAU_ITM

  implicit none

  private

  public :: atype_SMB_model, create_SMB_model

contains

  subroutine create_SMB_model( SMB, choice_SMB_model)
    !< Allocate a concrete implementation of a SMB_model

    ! In/output variables:
    class(atype_SMB_model), allocatable, intent(inout) :: SMB
    character(len=*),                    intent(in   ) :: choice_SMB_model

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_SMB_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    select case (choice_SMB_model)
    case default
      call crash('invalid choice_SMB_model "' // trim( choice_SMB_model) // '"')
    case ('idealised')
      allocate( type_SMB_model_idealised :: SMB)
    case ('prescribed')
      allocate( type_SMB_model_prescribed :: SMB)
    case ('reconstructed')
      allocate( type_SMB_model_reconstructed :: SMB)
    case ('snapshot_plus_anomalies')
      allocate( type_SMB_model_snp_p_anml :: SMB)
    case ('IMAU-ITM')
      allocate( type_SMB_model_IMAU_ITM :: SMB)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_SMB_model

end module SMB_model