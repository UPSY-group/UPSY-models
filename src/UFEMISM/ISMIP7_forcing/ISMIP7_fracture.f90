module ISMIP7_fracture

  ! The ISMIP7 protocol provides a single NetCDF file containing a mask
  ! for every year (1950-2299). The mask is integer-valued (as NetCDF
  ! doesn't support logicals), with 1 indicating fracture, implying
  ! that no (floating?) ice is allowed to exist there anymore.

  ! This setting is controlled by the config parameter 'do_apply_ISMIP7_fracture_mask_config'

  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash

  implicit none

  private

contains

end module ISMIP7_fracture
