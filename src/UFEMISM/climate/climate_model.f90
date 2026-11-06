module climate_model

	use call_stack_and_comp_time_tracking,    only: init_routine, finalise_routine, crash
    use climate_model_basic,                  only: atype_climate_model
    use climate_idealised,                    only: type_climate_model_idealised
    use climate_snapshot,                     only: type_climate_model_snapshot
    use climate_snapshot_lapse,               only: type_climate_model_snapshot_lapse
    ! use climate_snapshot_plus_anomalies,    only: type_climate_model_snapshot_plus_anomalies
    ! use climate_snapshot_plus_transient_deltaT, only: type_climate_snapshot_plus_transient_deltaT
    use climate_snapshot_plus_uniform_deltaT, only: type_climate_model_snapshot_plus_uniform_deltaT
    use ISMIP7_climate,                       only: type_climate_model_ISMIP7
	use TIPMIP_climate,                       only: type_climate_model_TIPMIP

	implicit none

	private

	public :: atype_climate_model, create_climate_model

  contains

  subroutine create_climate_model(climate, choice_climate_model)

	class(atype_climate_model), allocatable, intent(inout) :: climate
	character(len=*), intent(in) :: choice_climate_model
	
	character(len=1024), parameter :: routine_name = 'create_climate_model'
	
	! Add routine to call stack
	call init_routine( routine_name)
	
	select case (choice_climate_model)
	    case default
	      call crash('invalid choice_climate_model "' // trim( choice_climate_model) // '"')
	    case ('idealised')
	      allocate( type_climate_model_idealised :: climate)
	    case ('snapshot')
	      allocate( type_climate_model_snapshot :: climate)
	    case ('snapshot_lapse')
	      allocate( type_climate_model_snapshot_lapse :: climate)
		case ('ISMIP7')
	      allocate( type_climate_model_ISMIP7 :: climate)
		case ('TIPMIP')
	      allocate( type_climate_model_TIPMIP :: climate)    
	    ! case ('snapshot_plus_anomalies')
	    !   allocate( type_climate_model_snapshot_p_anml :: climate)
	    ! case ('snapshot_plus_transient_deltaT')
	    !   allocate( type_climate_model_snapshot_trans_dT :: climate)
	    case ('snapshot_plus_uniform_deltaT')
	      allocate( type_climate_model_snapshot_plus_uniform_deltaT :: climate)
	    end select
	
	    ! Remove routine from call stack
	    call finalise_routine( routine_name)
	
	end subroutine create_climate_model

end module climate_model