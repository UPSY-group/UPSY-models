module climate_retreat_mask

! ===== Preamble =====
! ====================

  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, sync
  use call_stack_and_comp_time_tracking                      , only: crash, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C
  use parameters
  use mesh_types                                             , only: type_mesh
  use ice_model_types                                        , only: type_ice_model
  use climate_model_types                                    , only: type_climate_model
  use reallocate_mod                                         , only: reallocate_bounds
  use netcdf_io_main
  use assertions_basic, only: assert

 implicit none

  private

  public :: run_climate_retreat_mask
  public :: update_ISMIP_style_future_timeframes
  public :: initialise_climate_retreat_mask

contains

! == ISMIP-style Antarctica future retreat mask
! ======================================================================

  subroutine run_climate_retreat_mask( mesh, climate, time, ice)

    ! In/output variables:
    type(type_mesh),                    intent(in)    :: mesh
    type(type_climate_model),           intent(inout) :: climate
    real(dp),                           intent(in)    :: time
    type(type_ice_model),               intent(in)    :: ice

    ! Local variables:
    character(len=256), parameter                           :: routine_name = 'run_climate_retreat_mask'
    integer                                                 :: vi
    real(dp)                                                :: wt0, wt1

    ! Add routine to path
    call init_routine( routine_name)

    if (.not. C%retreat_mask_without_time) then

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    if (time < climate%ISMIP_style%shelf_collapse_mask_t0 .OR. time > climate%ISMIP_style%shelf_collapse_mask_t1) then

      ! Find and read the two global time frames
      call update_ISMIP_style_future_timeframes( mesh, climate, time)

    end if ! IF (time >= climate_matrix%SMB_direct%t0 .AND. time <= climate_matrix%SMB_direct%t1) THEN

    ! Interpolate the two timeframes in time
    wt0 = (climate%ISMIP_style%shelf_collapse_mask_t1 - time) / (climate%ISMIP_style%shelf_collapse_mask_t1 - climate%ISMIP_style%shelf_collapse_mask_t0)
    wt1 = 1._dp - wt0

    do vi = mesh%vi1, mesh%vi2

      climate%ISMIP_style%shelf_collapse_mask( vi) = (wt0 * climate%ISMIP_style%shelf_collapse_mask0( vi)) + &
                                                        (wt1 * climate%ISMIP_style%shelf_collapse_mask1( vi))

    end do

    else ! if (C%retreat_mask_without_time)
      if (time == C%start_time_of_run) then
        ! If the retreat mask is time-independent, just read it once and keep it constant
        call read_field_from_file_2D( C%ISMIP_future_shelf_collapse_forcing_filename, 'mask', mesh, C%output_dir, climate%ISMIP_style%shelf_collapse_mask)
      else
        ! do nothing, no need to read the mask again
      end if ! if (time == C%start_time_of_run)
    end if
    
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_climate_retreat_mask

  subroutine update_ISMIP_style_future_timeframes( mesh, climate, time)
    ! Update the two timeframes of the ISMIP-style forcing data

    ! In/output variables:
    type(type_mesh),                          intent(in   ) :: mesh
    type(type_climate_model),                 intent(inout) :: climate
    real(dp),                                 intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                           :: routine_name = 'update_ISMIP_style_future_timeframes'
    real(dp)                                                :: time0, time1

    ! Add routine to path
    call init_routine( routine_name)

    time0= real( floor( time, dp), dp)
    time1= real( floor( time, dp), dp) + 10._dp
    
    ! Read timeframes from file
    call read_field_from_file_2D( C%ISMIP_future_shelf_collapse_forcing_filename, 'mask', mesh, C%output_dir, climate%ISMIP_style%shelf_collapse_mask0, time_to_read = time0)
    call read_field_from_file_2D( C%ISMIP_future_shelf_collapse_forcing_filename, 'mask', mesh, C%output_dir, climate%ISMIP_style%shelf_collapse_mask1, time_to_read = time1)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_ISMIP_style_future_timeframes

  subroutine initialise_climate_retreat_mask( mesh, climate)
    ! Use the ISMIP-style forcing

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_climate_model),           intent(inout) :: climate

    ! Local variables:
    character(len=256), parameter                           :: routine_name = 'initialise_climate_retreat_mask'

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate memory
    allocate( climate%ISMIP_style%shelf_collapse_mask0( mesh%vi1:mesh%vi2))
    allocate( climate%ISMIP_style%shelf_collapse_mask1( mesh%vi1:mesh%vi2))
    allocate( climate%ISMIP_style%shelf_collapse_mask( mesh%vi1:mesh%vi2))

    ! Give impossible values to timeframes, so that the first call to run_climate_model_ISMIP_style
    ! is guaranteed to first read two new timeframes from the NetCDF file
    climate%ISMIP_style%shelf_collapse_mask_t0 = C%start_time_of_run - 100._dp
    climate%ISMIP_style%shelf_collapse_mask_t1 = C%start_time_of_run - 90._dp

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_climate_retreat_mask

end module climate_retreat_mask