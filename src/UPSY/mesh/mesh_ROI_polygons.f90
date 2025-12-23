module mesh_ROI_polygons

  ! Pre-defined regions of interest in Greenland and Antarctica

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine

  implicit none

  private

  public :: calc_polygon_Pine_Island_Glacier
  public :: calc_polygon_Thwaites_Glacier
  public :: calc_polygon_Amery_ice_shelf
  public :: calc_polygon_Riiser_Larsen_ice_shelf
  public :: calc_polygon_Siple_Coast
  public :: calc_polygon_Larsen_ice_shelf
  public :: calc_polygon_Transantarctic_Mountains
  public :: calc_polygon_DotsonCrosson_ice_shelf
  public :: calc_polygon_Patagonia
  public :: calc_polygon_Narsarsuaq
  public :: calc_polygon_Nuuk
  public :: calc_polygon_Jakobshavn
  public :: calc_polygon_NGIS
  public :: calc_polygon_Qaanaaq
  public :: calc_polygon_CalvMIP_quarter
  public :: calc_polygon_Franka_WAIS
  public :: calc_polygon_Dotson_channel
  public :: calc_polygon_Mulock_glacier
  public :: calc_polygon_Byrd_glacier
  public :: calc_polygon_Nimrod_glacier
  public :: calc_polygon_Beardmore_glacier
  public :: calc_polygon_Shackleton_glacier
  public :: calc_polygon_Amundsen_glacier
  public :: calc_polygon_Scott_glacier
  public :: calc_polygon_Mercer_glacier
  public :: calc_polygon_Wilkes_basins
  public :: calc_polygon_Institute_basin

contains

subroutine calc_polygon_Pine_Island_Glacier( poly)
  ! Return a polygon enveloping the Pine Island Glacier catchment basin
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Pine_Island_Glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 42,2))

  poly(  1,:) = [ -1.64e6_dp, -3.4e5_dp]
  poly(  2,:) = [-1.60e6_dp, -3.5e5_dp]
  poly(  3,:) = [-1.55e6_dp, -3.4e5_dp]
  poly(  4,:) = [-1.50e6_dp, -3.2e5_dp]
  poly(  5,:) = [-1.45e6_dp, -2.9e5_dp]
  poly(  6,:) = [-1.40e6_dp, -2.5e5_dp]
  poly(  7,:) = [-1.37e6_dp, -2.0e5_dp]
  poly(  8,:) = [-1.34e6_dp, -1.7e5_dp]
  poly(  9,:) = [-1.30e6_dp, -1.6e5_dp]
  poly( 10,:) = [-1.26e6_dp, -1.6e5_dp]
  poly( 11,:) = [-1.22e6_dp, -1.7e5_dp]
  poly( 12,:) = [-1.18e6_dp, -1.75e5_dp]
  poly( 13,:) = [-1.14e6_dp, -1.75e5_dp]
  poly( 14,:) = [-1.11e6_dp, -1.72e5_dp]
  poly( 15,:) = [-1.09e6_dp, -1.6e5_dp]
  poly( 16,:) = [-1.085e6_dp, -1.4e5_dp]
  poly( 17,:) = [-1.09e6_dp, -1.2e5_dp]
  poly( 18,:) = [-1.1e6_dp, -1.0e5_dp]
  poly( 19,:) = [-1.13e6_dp, -0.7e5_dp]
  poly( 20,:) = [-1.17e6_dp, -0.4e5_dp]
  poly( 21,:) = [-1.21e6_dp, -0.2e5_dp]
  poly( 22,:) = [-1.26e6_dp, -0.0e5_dp]
  poly( 23,:) = [-1.32e6_dp, 0.1e5_dp]
  poly( 24,:) = [-1.45e6_dp, 0.1e5_dp]
  poly( 25,:) = [-1.48e6_dp, 0.15e5_dp]
  poly( 26,:) = [-1.51e6_dp, 0.35e5_dp]
  poly( 27,:) = [-1.53e6_dp, 0.75e5_dp]
  poly( 28,:) = [-1.55e6_dp, 0.95e5_dp]
  poly( 29,:) = [-1.58e6_dp, 0.1e6_dp]
  poly( 30,:) = [-1.62e6_dp, 0.11e6_dp]
  poly( 31,:) = [-1.65e6_dp, 0.12e6_dp]
  poly( 32,:) = [-1.67e6_dp, 0.10e6_dp]
  poly( 33,:) = [-1.69e6_dp, 0.9e5_dp]
  poly( 34,:) = [-1.71e6_dp, 0.5e5_dp]
  poly( 35,:) = [-1.74e6_dp, 0.1e5_dp]
  poly( 36,:) = [-1.75e6_dp, -0.5e5_dp]
  poly( 37,:) = [-1.75e6_dp, -0.15e6_dp]
  poly( 38,:) = [-1.71e6_dp, -0.19e6_dp]
  poly( 39,:) = [-1.66e6_dp, -0.2e6_dp]
  poly( 40,:) = [-1.64e6_dp, -0.21e6_dp]
  poly( 41,:) = [-1.63e6_dp, -0.23e6_dp]
  poly( 42,:) = [-1.63e6_dp, -0.29e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Pine_Island_Glacier

subroutine calc_polygon_Thwaites_Glacier( poly)
  ! Return a polygon enveloping the Pine Island Glacier catchment basin
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Thwaites_Glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 44,2))

  poly(  1,:) = [-1.6e6_dp, -5.4e5_dp]
  poly(  2,:) = [-1.55e6_dp, -5.4e5_dp]
  poly(  3,:) = [-1.50e6_dp, -5.5e5_dp]
  poly(  4,:) = [-1.45e6_dp, -5.6e5_dp]
  poly(  5,:) = [-1.40e6_dp, -5.65e5_dp]
  poly(  6,:) = [-1.37e6_dp, -5.75e5_dp]
  poly(  7,:) = [-1.35e6_dp, -6e5_dp]
  poly(  8,:) = [-1.35e6_dp, -6.5e5_dp]
  poly(  9,:) = [-1.34e6_dp, -6.9e5_dp]
  poly( 10,:) = [-1.32e6_dp, -7.3e5_dp]
  poly( 11,:) = [-1.29e6_dp, -7.6e5_dp]
  poly( 12,:) = [-1.25e6_dp, -7.8e5_dp]
  poly( 13,:) = [-1.22e6_dp, -7.8e5_dp]
  poly( 14,:) = [-1.20e6_dp, -7.6e5_dp]
  poly( 15,:) = [-1.18e6_dp, -7.4e5_dp]
  poly( 16,:) = [-1.15e6_dp, -6.9e5_dp]
  poly( 17,:) = [-1.14e6_dp, -6.4e5_dp]
  poly( 18,:) = [-1.14e6_dp, -5.9e5_dp]
  poly( 19,:) = [-1.11e6_dp, -5.6e5_dp]
  poly( 20,:) = [-1.08e6_dp, -5.5e5_dp]
  poly( 21,:) = [-1.04e6_dp, -5.4e5_dp]
  poly( 22,:) = [-1.01e6_dp, -5.2e5_dp]
  poly( 23,:) = [-0.99e6_dp, -5.0e5_dp]
  poly( 24,:) = [-0.99e6_dp, -4.6e5_dp]
  poly( 25,:) = [-1.02e6_dp, -4.4e5_dp]
  poly( 26,:) = [-1.04e6_dp, -4.2e5_dp]
  poly( 27,:) = [-1.06e6_dp, -3.9e5_dp]
  poly( 28,:) = [-1.07e6_dp, -3.5e5_dp]
  poly( 29,:) = [-1.07e6_dp, -3.2e5_dp]
  poly( 30,:) = [-1.09e6_dp, -2.8e5_dp]
  poly( 31,:) = [-1.12e6_dp, -2.5e5_dp]
  poly( 32,:) = [-1.15e6_dp, -2.2e5_dp]
  poly( 33,:) = [-1.18e6_dp, -1.9e5_dp]
  poly( 34,:) = [-1.22e6_dp, -1.7e5_dp]
  poly( 35,:) = [-1.26e6_dp, -1.6e5_dp]
  poly( 36,:) = [-1.30e6_dp, -1.6e5_dp]
  poly( 37,:) = [-1.34e6_dp, -1.7e5_dp]
  poly( 38,:) = [-1.37e6_dp, -2.0e5_dp]
  poly( 39,:) = [-1.40e6_dp, -2.5e5_dp]
  poly( 40,:) = [-1.45e6_dp, -2.9e5_dp]
  poly( 41,:) = [-1.50e6_dp, -3.2e5_dp]
  poly( 42,:) = [-1.55e6_dp, -3.4e5_dp]
  poly( 43,:) = [-1.60e6_dp, -3.5e5_dp]
  poly( 44,:) = [-1.64e6_dp, -3.4e5_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Thwaites_Glacier

subroutine calc_polygon_Amery_ice_shelf( poly)
  ! Return a polygon enveloping the Amery ice shelf catchment basin
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Amery_ice_shelf'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 11,2))

  poly(  1,:) = [2.2798e6_dp, 0.8624e6_dp]
  poly(  2,:) = [1.9637e6_dp, 0.9955e6_dp]
  poly(  3,:) = [1.4229e6_dp, 0.9234e6_dp]
  poly(  4,:) = [1.3480e6_dp, 0.7792e6_dp]
  poly(  5,:) = [1.2981e6_dp, 0.6711e6_dp]
  poly(  6,:) = [1.4340e6_dp, 0.4353e6_dp]
  poly(  7,:) = [1.6337e6_dp, 0.4742e6_dp]
  poly(  8,:) = [1.8056e6_dp, 0.5019e6_dp]
  poly(  9,:) = [1.8777e6_dp, 0.4215e6_dp]
  poly( 10,:) = [2.1079e6_dp, 0.4520e6_dp]
  poly( 11,:) = [2.3075e6_dp, 0.6711e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Amery_ice_shelf

subroutine calc_polygon_Riiser_Larsen_ice_shelf( poly)
  ! Return a polygon enveloping the Riiser-Larsen ice shelf catchment basin
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Riiser_Larsen_ice_shelf'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 31,2))

  poly(  1,:) = [-0.6469e6_dp, 1.6448e6_dp]
  poly(  2,:) = [-0.6507e6_dp, 1.7370e6_dp]
  poly(  3,:) = [-0.6411e6_dp, 1.8005e6_dp]
  poly(  4,:) = [-0.5989e6_dp, 1.8370e6_dp]
  poly(  5,:) = [-0.5508e6_dp, 1.8639e6_dp]
  poly(  6,:) = [-0.5104e6_dp, 1.9081e6_dp]
  poly(  7,:) = [-0.4758e6_dp, 1.9331e6_dp]
  poly(  8,:) = [-0.4451e6_dp, 1.9542e6_dp]
  poly(  9,:) = [-0.4393e6_dp, 1.9946e6_dp]
  poly( 10,:) = [-0.3336e6_dp, 1.9696e6_dp]
  poly( 11,:) = [-0.3048e6_dp, 1.9292e6_dp]
  poly( 12,:) = [-0.2644e6_dp, 1.9081e6_dp]
  poly( 13,:) = [-0.2029e6_dp, 1.8927e6_dp]
  poly( 14,:) = [-0.1741e6_dp, 1.8716e6_dp]
  poly( 15,:) = [-0.1644e6_dp, 1.8351e6_dp]
  poly( 16,:) = [-0.1414e6_dp, 1.8043e6_dp]
  poly( 17,:) = [-0.1222e6_dp, 1.7659e6_dp]
  poly( 18,:) = [-0.1202e6_dp, 1.7313e6_dp]
  poly( 19,:) = [-0.1318e6_dp, 1.6928e6_dp]
  poly( 20,:) = [-0.1644e6_dp, 1.6640e6_dp]
  poly( 21,:) = [-0.2125e6_dp, 1.6275e6_dp]
  poly( 22,:) = [-0.2394e6_dp, 1.5948e6_dp]
  poly( 23,:) = [-0.2663e6_dp, 1.5833e6_dp]
  poly( 24,:) = [-0.3259e6_dp, 1.5813e6_dp]
  poly( 25,:) = [-0.3778e6_dp, 1.5717e6_dp]
  poly( 26,:) = [-0.4201e6_dp, 1.5640e6_dp]
  poly( 27,:) = [-0.4528e6_dp, 1.5640e6_dp]
  poly( 28,:) = [-0.4931e6_dp, 1.5660e6_dp]
  poly( 29,:) = [-0.5354e6_dp, 1.5698e6_dp]
  poly( 30,:) = [-0.5758e6_dp, 1.5871e6_dp]
  poly( 31,:) = [-0.6142e6_dp, 1.6102e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Riiser_Larsen_ice_shelf

subroutine calc_polygon_Siple_Coast( poly)
  ! Return a polygon enveloping the Siple Coast area
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Siple_Coast'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 21,2))

  poly(  1,:) = [-0.6394e6_dp, -0.2184e6_dp]
  poly(  2,:) = [-0.6852e6_dp, -0.3498e6_dp]
  poly(  3,:) = [-0.7219e6_dp, -0.3101e6_dp]
  poly(  4,:) = [-0.8165e6_dp, -0.2979e6_dp]
  poly(  5,:) = [-0.8288e6_dp, -0.3681e6_dp]
  poly(  6,:) = [-0.7402e6_dp, -0.4567e6_dp]
  poly(  7,:) = [-1.0059e6_dp, -0.3803e6_dp]
  poly(  8,:) = [-1.0029e6_dp, -0.4689e6_dp]
  poly(  9,:) = [-0.9326e6_dp, -0.5514e6_dp]
  poly( 10,:) = [-0.8440e6_dp, -0.6125e6_dp]
  poly( 11,:) = [-1.0609e6_dp, -0.6033e6_dp]
  poly( 12,:) = [-0.8807e6_dp, -0.6980e6_dp]
  poly( 13,:) = [-1.0273e6_dp, -0.7652e6_dp]
  poly( 14,:) = [-1.0609e6_dp, -0.9210e6_dp]
  poly( 15,:) = [-0.9876e6_dp, -1.0737e6_dp]
  poly( 16,:) = [-0.7463e6_dp, -1.0004e6_dp]
  poly( 17,:) = [-0.6363e6_dp, -1.0981e6_dp]
  poly( 18,:) = [-0.5019e6_dp, -1.1287e6_dp]
  poly( 19,:) = [ 0.0051e6_dp, -0.8355e6_dp]
  poly( 20,:) = [-0.0132e6_dp, -0.2887e6_dp]
  poly( 21,:) = [-0.3034e6_dp, -0.1573e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Siple_Coast

subroutine calc_polygon_Larsen_ice_shelf( poly)
  ! Return a polygon enveloping the Larsen C ice shelf area
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Larsen_ice_shelf'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 5,2))

  poly(  1,:) = [-2.3962e6_dp, 0.8370e6_dp]
  poly(  2,:) = [-1.9819e6_dp, 0.8482e6_dp]
  poly(  3,:) = [-1.8363e6_dp, 1.0721e6_dp]
  poly(  4,:) = [-2.4857e6_dp, 1.6880e6_dp]
  poly(  5,:) = [-2.6985e6_dp, 1.3968e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Larsen_ice_shelf

subroutine calc_polygon_Transantarctic_Mountains( poly)
  ! Return a polygon enveloping the Transantarctic Mountains
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Transantarctic_Mountains'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 12,2))

  poly(  1,:) = [ 0.2911e6_dp, -1.3464e6_dp]
  poly(  2,:) = [ 0.5487e6_dp, -1.2233e6_dp]
  poly(  3,:) = [ 0.6158e6_dp, -1.1225e6_dp]
  poly(  4,:) = [ 0.5934e6_dp, -0.7978e6_dp]
  poly(  5,:) = [ 0.3695e6_dp, -0.5067e6_dp]
  poly(  6,:) = [ 0.1680e6_dp, -0.3387e6_dp]
  poly(  7,:) = [-0.1792e6_dp, -0.1708e6_dp]
  poly(  8,:) = [-0.4143e6_dp, -0.1484e6_dp]
  poly(  9,:) = [-0.4591e6_dp, -0.3947e6_dp]
  poly( 10,:) = [ 0.0672e6_dp, -0.7306e6_dp]
  poly( 11,:) = [ 0.2127e6_dp, -0.8762e6_dp]
  poly( 12,:) = [ 0.3359e6_dp, -1.0217e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Transantarctic_Mountains

subroutine calc_polygon_DotsonCrosson_ice_shelf( poly)
  ! Return a polygon enveloping the Dotson-Crosson ice shelf area
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_DotsonCrosson_ice_shelf'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 14,2))

  poly(  1,:) = [-1.5260e6_dp, -0.5303e6_dp]
  poly(  2,:) = [-1.4997e6_dp, -0.5339e6_dp]
  poly(  3,:) = [-1.4156e6_dp, -0.5703e6_dp]
  poly(  4,:) = [-1.3637e6_dp, -0.6060e6_dp]
  poly(  5,:) = [-1.4103e6_dp, -0.6627e6_dp]
  poly(  6,:) = [-1.3691e6_dp, -0.7253e6_dp]
  poly(  7,:) = [-1.4210e6_dp, -0.7212e6_dp]
  poly(  8,:) = [-1.4789e6_dp, -0.7021e6_dp]
  poly(  9,:) = [-1.5176e6_dp, -0.6949e6_dp]
  poly( 10,:) = [-1.5689e6_dp, -0.7074e6_dp]
  poly( 11,:) = [-1.6011e6_dp, -0.6955e6_dp]
  poly( 12,:) = [-1.6148e6_dp, -0.6013e6_dp]
  poly( 13,:) = [-1.5862e6_dp, -0.5488e6_dp]
  poly( 14,:) = [-1.5457e6_dp, -0.5219e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_DotsonCrosson_ice_shelf

subroutine calc_polygon_Patagonia( poly)
  ! Return a polygon enveloping the region where the former
  ! Patagonian ice sheet peaked during the last glacial maximum
  !
  ! (based on manual analysis of the PATICE reconstruction,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do mesh refinement)
  !
  ! This assumes a very particular stereographic projection for
  ! the model domain, which as of now reads:
  !
  ! lambda_M_ANT_config    = 289.0
  ! phi_M_ANT_config       = -47.0
  ! beta_stereo_ANT_config = 71.0
  ! xmin_ANT_config        = -400000.0
  ! xmax_ANT_config        =  400000.0
  ! ymin_ANT_config        = -1110000.0
  ! ymax_ANT_config        =  1110000.0

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Patagonia'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 130,2))

  poly(  1,:) = [-0.5409e5_dp,  0.8469e5_dp]
  poly(  2,:) = [-0.1447e5_dp,  0.7109e5_dp]
  poly(  3,:) = [-0.4643e5_dp,  0.4652e5_dp]
  poly(  4,:) = [-0.8643e5_dp,  0.3775e5_dp]
  poly(  5,:) = [-0.6965e5_dp, -0.0131e5_dp]
  poly(  6,:) = [-0.3570e5_dp, -0.2868e5_dp]
  poly(  7,:) = [-0.5385e5_dp, -0.6522e5_dp]
  poly(  8,:) = [-0.7107e5_dp, -1.0253e5_dp]
  poly(  9,:) = [-0.8151e5_dp, -1.4277e5_dp]
  poly( 10,:) = [-1.0274e5_dp, -1.8160e5_dp]
  poly( 11,:) = [-0.8547e5_dp, -2.1939e5_dp]
  poly( 12,:) = [-1.2204e5_dp, -2.3569e5_dp]
  poly( 13,:) = [-1.0413e5_dp, -2.7274e5_dp]
  poly( 14,:) = [-0.6790e5_dp, -2.9327e5_dp]
  poly( 15,:) = [-1.0808e5_dp, -2.9403e5_dp]
  poly( 16,:) = [-1.2326e5_dp, -3.3331e5_dp]
  poly( 17,:) = [-0.8064e5_dp, -3.3984e5_dp]
  poly( 18,:) = [-0.8601e5_dp, -3.8041e5_dp]
  poly( 19,:) = [-0.6702e5_dp, -4.1567e5_dp]
  poly( 20,:) = [-1.0990e5_dp, -4.1746e5_dp]
  poly( 21,:) = [-1.0944e5_dp, -4.5829e5_dp]
  poly( 22,:) = [-1.0065e5_dp, -5.0013e5_dp]
  poly( 23,:) = [-0.7625e5_dp, -5.3305e5_dp]
  poly( 24,:) = [-0.6571e5_dp, -5.7201e5_dp]
  poly( 25,:) = [-0.2665e5_dp, -5.9098e5_dp]
  poly( 26,:) = [-0.6042e5_dp, -6.1652e5_dp]
  poly( 27,:) = [-0.1109e5_dp, -6.1846e5_dp]
  poly( 28,:) = [-0.2399e5_dp, -6.6280e5_dp]
  poly( 29,:) = [-0.3598e5_dp, -7.0219e5_dp]
  poly( 30,:) = [-0.0160e5_dp, -6.7598e5_dp]
  poly( 31,:) = [ 0.1932e5_dp, -6.3562e5_dp]
  poly( 32,:) = [ 0.3721e5_dp, -6.7754e5_dp]
  poly( 33,:) = [ 0.7725e5_dp, -6.8812e5_dp]
  poly( 34,:) = [ 0.5899e5_dp, -7.2985e5_dp]
  poly( 35,:) = [ 0.6939e5_dp, -7.6885e5_dp]
  poly( 36,:) = [ 1.0711e5_dp, -7.9010e5_dp]
  poly( 37,:) = [ 1.4723e5_dp, -8.0512e5_dp]
  poly( 38,:) = [ 1.8838e5_dp, -8.2100e5_dp]
  poly( 39,:) = [ 2.3419e5_dp, -8.1880e5_dp]
  poly( 40,:) = [ 1.9120e5_dp, -8.2741e5_dp]
  poly( 41,:) = [ 2.2133e5_dp, -8.5626e5_dp]
  poly( 42,:) = [ 1.7064e5_dp, -8.5981e5_dp]
  poly( 43,:) = [ 1.3405e5_dp, -8.7757e5_dp]
  poly( 44,:) = [ 1.5893e5_dp, -9.2305e5_dp]
  poly( 45,:) = [ 1.2073e5_dp, -8.9282e5_dp]
  poly( 46,:) = [ 0.8931e5_dp, -9.2499e5_dp]
  poly( 47,:) = [ 0.4329e5_dp, -8.9716e5_dp]
  poly( 48,:) = [-0.0783e5_dp, -8.7450e5_dp]
  poly( 49,:) = [-0.5789e5_dp, -8.3810e5_dp]
  poly( 50,:) = [-1.3746e5_dp, -7.7720e5_dp]
  poly( 51,:) = [-1.7175e5_dp, -7.2582e5_dp]
  poly( 52,:) = [-2.0364e5_dp, -6.9912e5_dp]
  poly( 53,:) = [-2.3851e5_dp, -6.5499e5_dp]
  poly( 54,:) = [-2.7530e5_dp, -5.3553e5_dp]
  poly( 55,:) = [-2.8841e5_dp, -4.9232e5_dp]
  poly( 56,:) = [-2.9286e5_dp, -4.4711e5_dp]
  poly( 57,:) = [-3.0745e5_dp, -3.3945e5_dp]
  poly( 58,:) = [-3.2386e5_dp, -2.9319e5_dp]
  poly( 59,:) = [-3.2507e5_dp, -2.5020e5_dp]
  poly( 60,:) = [-3.3118e5_dp, -2.0832e5_dp]
  poly( 61,:) = [-3.2731e5_dp, -1.4841e5_dp]
  poly( 62,:) = [-3.1580e5_dp, -0.9411e5_dp]
  poly( 63,:) = [-2.7306e5_dp, -0.6904e5_dp]
  poly( 64,:) = [-2.3183e5_dp, -0.1091e5_dp]
  poly( 65,:) = [-2.2407e5_dp,  0.3389e5_dp]
  poly( 66,:) = [-1.9999e5_dp,  0.6980e5_dp]
  poly( 67,:) = [-2.0187e5_dp,  1.1338e5_dp]
  poly( 68,:) = [-2.0304e5_dp,  1.5620e5_dp]
  poly( 69,:) = [-1.8445e5_dp,  2.0050e5_dp]
  poly( 70,:) = [-1.8952e5_dp,  2.5343e5_dp]
  poly( 71,:) = [-2.0799e5_dp,  2.1778e5_dp]
  poly( 72,:) = [-2.0690e5_dp,  1.7765e5_dp]
  poly( 73,:) = [-2.1957e5_dp,  1.3968e5_dp]
  poly( 74,:) = [-2.1956e5_dp,  0.9968e5_dp]
  poly( 75,:) = [-2.3816e5_dp,  0.5980e5_dp]
  poly( 76,:) = [-2.6377e5_dp,  0.2888e5_dp]
  poly( 77,:) = [-3.0360e5_dp,  0.3382e5_dp]
  poly( 78,:) = [-3.0505e5_dp,  0.7408e5_dp]
  poly( 79,:) = [-3.0559e5_dp,  1.1415e5_dp]
  poly( 80,:) = [-2.7731e5_dp,  1.4247e5_dp]
  poly( 81,:) = [-2.6206e5_dp,  1.7961e5_dp]
  poly( 82,:) = [-2.6578e5_dp,  2.1960e5_dp]
  poly( 83,:) = [-2.6453e5_dp,  3.0531e5_dp]
  poly( 84,:) = [-2.3148e5_dp,  3.5499e5_dp]
  poly( 85,:) = [-2.3681e5_dp,  4.0124e5_dp]
  poly( 86,:) = [-2.2282e5_dp,  4.5178e5_dp]
  poly( 87,:) = [-2.1428e5_dp,  4.9543e5_dp]
  poly( 88,:) = [-1.9766e5_dp,  5.3229e5_dp]
  poly( 89,:) = [-1.8472e5_dp,  5.7164e5_dp]
  poly( 90,:) = [-1.5207e5_dp,  5.9587e5_dp]
  poly( 91,:) = [-1.1078e5_dp,  5.7979e5_dp]
  poly( 92,:) = [-1.1706e5_dp,  6.2322e5_dp]
  poly( 93,:) = [-1.5765e5_dp,  6.1397e5_dp]
  poly( 94,:) = [-1.1696e5_dp,  6.3771e5_dp]
  poly( 95,:) = [-1.3889e5_dp,  6.7127e5_dp]
  poly( 96,:) = [-0.9800e5_dp,  6.8870e5_dp]
  poly( 97,:) = [-1.2072e5_dp,  7.2285e5_dp]
  poly( 98,:) = [-1.0288e5_dp,  7.5925e5_dp]
  poly( 99,:) = [-1.0287e5_dp,  7.9975e5_dp]
  poly(100,:) = [-0.7675e5_dp,  8.3210e5_dp]
  poly(101,:) = [-0.6032e5_dp,  8.6911e5_dp]
  poly(102,:) = [-0.6162e5_dp,  9.1558e5_dp]
  poly(103,:) = [-0.5947e5_dp,  9.5816e5_dp]
  poly(104,:) = [ 0.0682e5_dp,  9.7742e5_dp]
  poly(105,:) = [ 0.1766e5_dp,  9.3773e5_dp]
  poly(106,:) = [-0.2248e5_dp,  9.5232e5_dp]
  poly(107,:) = [-0.2538e5_dp,  9.1150e5_dp]
  poly(108,:) = [-0.0421e5_dp,  8.7645e5_dp]
  poly(109,:) = [-0.2659e5_dp,  8.4246e5_dp]
  poly(110,:) = [-0.2557e5_dp,  8.0233e5_dp]
  poly(111,:) = [-0.2657e5_dp,  7.6161e5_dp]
  poly(112,:) = [-0.3503e5_dp,  7.2132e5_dp]
  poly(113,:) = [-0.1780e5_dp,  6.8387e5_dp]
  poly(114,:) = [-0.3278e5_dp,  6.4600e5_dp]
  poly(115,:) = [-0.1889e5_dp,  6.0735e5_dp]
  poly(116,:) = [-0.3383e5_dp,  5.6766e5_dp]
  poly(117,:) = [-0.3034e5_dp,  5.2726e5_dp]
  poly(118,:) = [-0.4553e5_dp,  4.8891e5_dp]
  poly(119,:) = [-0.4680e5_dp,  4.4792e5_dp]
  poly(120,:) = [-0.1740e5_dp,  4.2013e5_dp]
  poly(121,:) = [-0.1911e5_dp,  3.7940e5_dp]
  poly(122,:) = [-0.3552e5_dp,  3.4222e5_dp]
  poly(123,:) = [-0.1063e5_dp,  3.0976e5_dp]
  poly(124,:) = [-0.4183e5_dp,  2.8087e5_dp]
  poly(125,:) = [-0.7996e5_dp,  2.6792e5_dp]
  poly(126,:) = [-0.4217e5_dp,  2.5454e5_dp]
  poly(127,:) = [-0.7106e5_dp,  2.2639e5_dp]
  poly(128,:) = [-0.7213e5_dp,  1.8527e5_dp]
  poly(129,:) = [-0.8023e5_dp,  1.4559e5_dp]
  poly(130,:) = [-0.6977e5_dp,  1.0662e5_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Patagonia

subroutine calc_polygon_Narsarsuaq( poly)
  ! Return a polygon enveloping the Narsarsuaq area in Southern Greenland
  !
  ! (based on manual analysis of the Rignot velocity data,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Narsarsuaq'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 22,2))

  poly(  1,:) = [ 0.0182e6, -3.0978e6]
  poly(  2,:) = [ 0.0544e6, -3.1334e6]
  poly(  3,:) = [ 0.0550e6, -3.1469e6]
  poly(  4,:) = [ 0.0495e6, -3.1579e6]
  poly(  5,:) = [ 0.0556e6, -3.1634e6]
  poly(  6,:) = [ 0.0495e6, -3.1720e6]
  poly(  7,:) = [ 0.0354e6, -3.1781e6]
  poly(  8,:) = [ 0.0434e6, -3.2008e6]
  poly(  9,:) = [ 0.0403e6, -3.2162e6]
  poly( 10,:) = [ 0.0219e6, -3.2107e6]
  poly( 11,:) = [ 0.0035e6, -3.2174e6]
  poly( 12,:) = [-0.0131e6, -3.2217e6]
  poly( 13,:) = [-0.0247e6, -3.2254e6]
  poly( 14,:) = [-0.0775e6, -3.2015e6]
  poly( 15,:) = [-0.1075e6, -3.1518e6]
  poly( 16,:) = [-0.1088e6, -3.1285e6]
  poly( 17,:) = [-0.0990e6, -3.1064e6]
  poly( 18,:) = [-0.0830e6, -3.0953e6]
  poly( 19,:) = [-0.0511e6, -3.0800e6]
  poly( 20,:) = [-0.0321e6, -3.0708e6]
  poly( 21,:) = [-0.0180e6, -3.0555e6]
  poly( 22,:) = [ 0.0059e6, -3.0555e6]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Narsarsuaq

subroutine calc_polygon_Nuuk( poly)
  ! Return a polygon enveloping the Nuuk area in Southwest Greenland
  !
  ! (based on manual analysis of a nice figure I once saw,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Nuuk'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 8,2))

  poly(  1,:) = [ -0.3411e6, -2.7256e6]
  poly(  2,:) = [ -0.2326e6, -2.6803e6]
  poly(  3,:) = [ -0.1396e6, -2.6743e6]
  poly(  4,:) = [ -0.0955e6, -2.7781e6]
  poly(  5,:) = [ -0.1193e6, -2.9044e6]
  poly(  6,:) = [ -0.2219e6, -2.9271e6]
  poly(  7,:) = [ -0.3184e6, -2.9199e6]
  poly(  8,:) = [ -0.3578e6, -2.8210e6]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Nuuk

subroutine calc_polygon_Jakobshavn( poly)
  ! Return a polygon enveloping the Jakobshavn area in West Greenland
  !
  ! (based on manual analysis of a nice figure I once saw,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Jakobshavn'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 7,2))

  poly(  1,:) = [ -0.3003e6, -2.2592e6]
  poly(  2,:) = [ -0.2908e6, -2.1644e6]
  poly(  3,:) = [ -0.2114e6, -2.1430e6]
  poly(  4,:) = [ -0.1046e6, -2.1679e6]
  poly(  5,:) = [ -0.0833e6, -2.2901e6]
  poly(  6,:) = [ -0.1212e6, -2.3778e6]
  poly(  7,:) = [ -0.2624e6, -2.3731e6]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Jakobshavn

subroutine calc_polygon_NGIS( poly)
  ! Return a polygon enveloping the Northern Greenland ice stream
  !
  ! (based on manual analysis of a nice figure I once saw,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_NGIS'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 20,2))

  poly(  1,:) = [0.5064e6_dp, -0.9936e6_dp]
  poly(  2,:) = [0.4467e6_dp, -1.0000e6_dp]
  poly(  3,:) = [0.3999e6_dp, -1.0469e6_dp]
  poly(  4,:) = [0.2805e6_dp, -1.0959e6_dp]
  poly(  5,:) = [0.2699e6_dp, -1.1322e6_dp]
  poly(  6,:) = [0.3338e6_dp, -1.1556e6_dp]
  poly(  7,:) = [0.3658e6_dp, -1.1471e6_dp]
  poly(  8,:) = [0.3295e6_dp, -1.2068e6_dp]
  poly(  9,:) = [0.2869e6_dp, -1.3261e6_dp]
  poly( 10,:) = [0.2208e6_dp, -1.4625e6_dp]
  poly( 11,:) = [0.2017e6_dp, -1.6308e6_dp]
  poly( 12,:) = [0.3295e6_dp, -1.5072e6_dp]
  poly( 13,:) = [0.5022e6_dp, -1.3453e6_dp]
  poly( 14,:) = [0.5362e6_dp, -1.3900e6_dp]
  poly( 15,:) = [0.5703e6_dp, -1.3794e6_dp]
  poly( 16,:) = [0.5938e6_dp, -1.3261e6_dp]
  poly( 17,:) = [0.5320e6_dp, -1.2004e6_dp]
  poly( 18,:) = [0.5576e6_dp, -1.1812e6_dp]
  poly( 19,:) = [0.5533e6_dp, -1.1045e6_dp]
  poly( 20,:) = [0.5405e6_dp, -1.0469e6_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_NGIS

subroutine calc_polygon_Qaanaaq( poly)
  ! Return a polygon enveloping the Qaanaaq area in Northwest Greenland
  !
  ! (based on manual analysis of a nice figure I once saw,
  ! not meant for basin-integrated SMB stuff or such, but accurate
  ! enough to do region-specific mesh refinement)

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Qaanaaq'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 10,2))

  poly(  1,:) = [ -0.5994e6, -1.1429e6]
  poly(  2,:) = [ -0.4913e6, -1.1981e6]
  poly(  3,:) = [ -0.4315e6, -1.2119e6]
  poly(  4,:) = [ -0.2681e6, -1.3132e6]
  poly(  5,:) = [ -0.4062e6, -1.3500e6]
  poly(  6,:) = [ -0.4683e6, -1.3523e6]
  poly(  7,:) = [ -0.6201e6, -1.3270e6]
  poly(  8,:) = [ -0.6270e6, -1.2579e6]
  poly(  9,:) = [ -0.5764e6, -1.2395e6]
  poly( 10,:) = [ -0.6063e6, -1.1751e6]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Qaanaaq

subroutine calc_polygon_CalvMIP_quarter( poly)
  ! Return a polygon enveloping one of the radially
  ! symmetrical quarters of the domain

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_CalvMIP_quarter'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 4,2))

  poly(  1,:) = [  0.0_dp,    0.0_dp]
  poly(  2,:) = [  0.0_dp,   8.e5_dp]
  poly(  3,:) = [ 8.e5_dp,   8.e5_dp]
  poly(  4,:) = [ 8.e5_dp,     0._dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_CalvMIP_quarter

subroutine calc_polygon_Franka_WAIS( poly)
  ! Return a polygon enveloping PIG, THW, and CD drainage basins

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Franka_WAIS'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 49,2))

  poly(  1,:) = [ -1.067649e6_dp, -2.389460e5_dp]
  poly(  2,:) = [ -1.083626e6_dp, -2.749513e5_dp]
  poly(  3,:) = [ -1.099971e6_dp, -3.261952e5_dp]
  poly(  4,:) = [ -1.086908e6_dp, -3.794976e5_dp]
  poly(  5,:) = [ -1.086874e6_dp, -4.297385e5_dp]
  poly(  6,:) = [ -1.106102e6_dp, -4.962918e5_dp]
  poly(  7,:) = [ -1.127440e6_dp, -5.594330e5_dp]
  poly(  8,:) = [ -1.150831e6_dp, -6.302047e5_dp]
  poly(  9,:) = [ -1.169688e6_dp, -6.953044e5_dp]
  poly( 10,:) = [ -1.199417e6_dp, -7.510359e5_dp]
  poly( 11,:) = [ -1.218135e6_dp, -8.169631e5_dp]
  poly( 12,:) = [ -1.262344e6_dp, -7.980657e5_dp]
  poly( 13,:) = [ -1.306049e6_dp, -7.573541e5_dp]
  poly( 14,:) = [ -1.369227e6_dp, -7.392547e5_dp]
  poly( 15,:) = [ -1.431547e6_dp, -7.236201e5_dp]
  poly( 16,:) = [ -1.494554e6_dp, -7.078608e5_dp]
  poly( 17,:) = [ -1.557047e6_dp, -7.071623e5_dp]
  poly( 18,:) = [ -1.595950e6_dp, -6.934250e5_dp]
  poly( 19,:) = [ -1.733956e6_dp, -3.668159e5_dp]
  poly( 20,:) = [ -1.717921e6_dp, -2.903907e5_dp]
  poly( 21,:) = [ -1.725877e6_dp, -2.215324e5_dp]
  poly( 22,:) = [ -1.738146e6_dp, -1.729808e5_dp]
  poly( 23,:) = [ -1.781543e6_dp, -1.265068e5_dp]
  poly( 24,:) = [ -1.756807e6_dp, -8.943395e4_dp]
  poly( 25,:) = [ -1.746682e6_dp, -4.915999e4_dp]
  poly( 26,:) = [ -1.755599e6_dp, -1.686297e4_dp]
  poly( 27,:) = [ -1.741802e6_dp, 1.058750e4_dp]
  poly( 28,:) = [ -1.712282e6_dp, 3.618599e4_dp]
  poly( 29,:) = [ -1.683568e6_dp, 8.043884e4_dp]
  poly( 30,:) = [ -1.650550e6_dp, 9.680165e4_dp]
  poly( 31,:) = [ -1.603154e6_dp, 8.566953e4_dp]
  poly( 32,:) = [ -1.542980e6_dp, 8.095809e4_dp]
  poly( 33,:) = [ -1.509540e6_dp, 4.239154e4_dp]
  poly( 34,:) = [ -1.472178e6_dp, -1.924253e3_dp]
  poly( 35,:) = [ -1.433291e6_dp, -3.426486e4_dp]
  poly( 36,:) = [ -1.410628e6_dp, -4.678015e4_dp]
  poly( 37,:) = [ -1.389430e6_dp, -3.021530e4_dp]
  poly( 38,:) = [ -1.352290e6_dp, -1.080577e4_dp]
  poly( 39,:) = [ -1.317960e6_dp, 1.190183e4_dp]
  poly( 40,:) = [ -1.272197e6_dp, 3.837878e4_dp]
  poly( 41,:) = [ -1.257608e6_dp, 6.008563e4_dp]
  poly( 42,:) = [ -1.219491e6_dp, 1.980493e4_dp]
  poly( 43,:) = [ -1.179582e6_dp, 1.178385e4_dp]
  poly( 44,:) = [ -1.132424e6_dp, -1.729169e4_dp]
  poly( 45,:) = [ -1.093011e6_dp, -4.302789e4_dp]
  poly( 46,:) = [ -1.042725e6_dp, -8.342796e4_dp]
  poly( 47,:) = [ -1.000538e6_dp, -1.366113e5_dp]
  poly( 48,:) = [ -1.009121e6_dp, -1.815426e5_dp]
  poly( 49,:) = [ -1.048750e6_dp, -2.290763e5_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Franka_WAIS

subroutine calc_polygon_Dotson_channel( poly)
  ! Return a polygon of the Dotson channel

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=256), parameter :: routine_name = 'calc_polygon_Dotson_channel'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 5,2))

  poly(  1,:) = [ -1564290_dp, -673200_dp]
  poly(  2,:) = [ -1555650_dp, -673200_dp]
  poly(  3,:) = [ -1541070_dp, -664560_dp]
  poly(  4,:) = [ -1545660_dp, -656730_dp]
  poly(  5,:) = [ -1565910_dp, -662400_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Dotson_channel

subroutine calc_polygon_Mulock_glacier( poly)
  ! Return a polygon enveloping Mulock glacier in the Transantarctic Mountains

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Mulock_glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 16,2))

  poly( 1,:) = [3.90e5,  -1.105e6]
  poly( 2,:) = [4.05e5,  -1.115e6]
  poly( 3,:) = [4.20e5,  -1.120e6]
  poly( 4,:) = [4.30e5,  -1.120e6]
  poly( 5,:) = [4.40e5,  -1.115e6]
  poly( 6,:) = [4.60e5,  -1.110e6]
  poly( 7,:) = [4.80e5,  -1.105e6]
  poly( 8,:) = [5.00e5,  -1.130e6]
  poly( 9,:) = [4.80e5,  -1.155e6]
  poly(10,:) = [4.60e5,  -1.150e6]
  poly(11,:) = [4.40e5,  -1.150e6]
  poly(12,:) = [4.30e5,  -1.140e6]
  poly(13,:) = [4.20e5,  -1.135e6]
  poly(14,:) = [4.05e5,  -1.130e6]
  poly(15,:) = [3.90e5,  -1.130e6]
  poly(16,:) = [3.70e5,  -1.115e6]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Mulock_glacier

subroutine calc_polygon_Byrd_glacier( poly)
  ! Return a polygon enveloping Byrd glacier in the Transantarctic Mountains

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Byrd_glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 12,2))

  poly( 1,:) = [3.45e5,  -9.90e5]
  poly( 2,:) = [3.70e5,  -9.75e5]
  poly( 3,:) = [3.90e5,  -9.50e5]
  poly( 4,:) = [4.05e5,  -9.20e5]
  poly( 5,:) = [4.10e5,  -9.00e5]
  poly( 6,:) = [4.30e5,  -9.00e5]
  poly( 7,:) = [4.30e5,  -9.20e5]
  poly( 8,:) = [4.15e5,  -9.50e5]
  poly( 9,:) = [3.92e5,  -9.80e5]
  poly(10,:) = [3.67e5,  -10.1e5]
  poly(11,:) = [3.50e5,  -10.35e5]
  poly(12,:) = [3.25e5,  -10.2e5]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Byrd_glacier

subroutine calc_polygon_Nimrod_glacier( poly)
  ! Return a polygon enveloping Nimrod glacier in the Transantarctic Mountains

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Nimrod_glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 24,2))

  poly( 1,:) = [2.20e5,  -8.00e5]
  poly( 2,:) = [2.60e5,  -7.72e5]
  poly( 3,:) = [2.70e5,  -7.50e5]
  poly( 4,:) = [2.78e5,  -7.30e5]
  poly( 5,:) = [2.75e5,  -7.15e5]
  poly( 6,:) = [2.45e5,  -6.50e5]
  poly( 7,:) = [2.45e5,  -6.10e5]
  poly( 8,:) = [2.67e5,  -6.20e5]
  poly( 9,:) = [2.67e5,  -6.60e5]
  poly(10,:) = [2.75e5,  -6.90e5]
  poly(11,:) = [2.85e5,  -7.03e5]
  poly(12,:) = [2.92e5,  -7.20e5]
  poly(13,:) = [3.10e5,  -7.00e5]
  poly(14,:) = [3.05e5,  -6.60e5]
  poly(15,:) = [3.30e5,  -6.60e5]
  poly(16,:) = [3.30e5,  -7.20e5]
  poly(17,:) = [3.42e5,  -7.18e5]
  poly(18,:) = [3.42e5,  -7.25e5]
  poly(19,:) = [3.30e5,  -7.28e5]
  poly(20,:) = [3.20e5,  -7.60e5]
  poly(21,:) = [2.76e5,  -7.77e5]
  poly(22,:) = [2.55e5,  -8.05e5]
  poly(23,:) = [2.35e5,  -8.10e5]
  poly(24,:) = [2.20e5,  -8.15e5]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Nimrod_glacier

subroutine calc_polygon_Beardmore_glacier( poly)
  ! Return a polygon enveloping Beardmore glacier in the Transantarctic Mountains

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Beardmore_glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 22,2))

  poly( 1,:) = [0.80e5,  -6.85e5]
  poly( 2,:) = [0.98e5,  -6.73e5]
  poly( 3,:) = [0.98e5,  -6.65e5]
  poly( 4,:) = [0.88e5,  -6.50e5]
  poly( 5,:) = [0.88e5,  -6.35e5]
  poly( 6,:) = [0.93e5,  -6.25e5]
  poly( 7,:) = [0.85e5,  -6.10e5]
  poly( 8,:) = [0.85e5,  -5.95e5]
  poly( 9,:) = [1.25e5,  -5.15e5]
  poly(10,:) = [1.38e5,  -5.02e5]
  poly(11,:) = [1.40e5,  -4.60e5]
  poly(12,:) = [1.70e5,  -4.60e5]
  poly(13,:) = [1.70e5,  -4.97e5]
  poly(14,:) = [1.50e5,  -5.15e5]
  poly(15,:) = [1.40e5,  -5.40e5]
  poly(16,:) = [1.35e5,  -5.75e5]
  poly(17,:) = [1.15e5,  -5.95e5]
  poly(18,:) = [1.10e5,  -6.05e5]
  poly(19,:) = [1.10e5,  -6.50e5]
  poly(20,:) = [1.15e5,  -6.65e5]
  poly(21,:) = [1.10e5,  -6.90e5]
  poly(22,:) = [0.90e5,  -7.00e5]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Beardmore_glacier

subroutine calc_polygon_Shackleton_glacier( poly)
  ! Return a polygon enveloping Shackleton glacier in the Transantarctic Mountains

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Shackleton_glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 12,2))

  poly( 1,:) = [-0.40e5, -6.05e5]
  poly( 2,:) = [-0.30e5, -5.90e5]
  poly( 3,:) = [-0.30e5, -5.25e5]
  poly( 4,:) = [-0.15e5, -5.10e5]
  poly( 5,:) = [-0.00e5, -4.95e5]
  poly( 6,:) = [-0.00e5, -4.70e5]
  poly( 7,:) = [-0.15e5, -4.85e5]
  poly( 8,:) = [-0.25e5, -5.00e5]
  poly( 9,:) = [-0.25e5, -5.10e5]
  poly(10,:) = [-0.40e5, -5.15e5]
  poly(11,:) = [-0.40e5, -5.80e5]
  poly(12,:) = [-0.55e5, -6.05e5]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Shackleton_glacier

subroutine calc_polygon_Amundsen_glacier( poly)
  ! Return a polygon enveloping Amundsen glacier in the Transantarctic Mountains

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Amundsen_glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 10,2))

  poly( 1,:) = [-1.60e5, -4.95e5]
  poly( 2,:) = [-1.65e5, -4.65e5]
  poly( 3,:) = [-1.50e5, -4.25e5]
  poly( 4,:) = [-1.40e5, -4.05e5]
  poly( 5,:) = [-1.20e5, -4.00e5]
  poly( 6,:) = [-1.15e5, -3.85e5]
  poly( 7,:) = [-1.40e5, -3.95e5]
  poly( 8,:) = [-1.55e5, -4.15e5]
  poly( 9,:) = [-1.75e5, -4.50e5]
  poly(10,:) = [-1.85e5, -4.80e5]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Amundsen_glacier

subroutine calc_polygon_Scott_glacier( poly)
  ! Return a polygon enveloping Scott glacier in the Transantarctic Mountains

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Scott_glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 14,2))

  poly( 1,:) = [-1.95e5, -4.70e5]
  poly( 2,:) = [-2.15e5, -4.45e5]
  poly( 3,:) = [-2.05e5, -4.15e5]
  poly( 4,:) = [-1.95e5, -3.90e5]
  poly( 5,:) = [-1.95e5, -3.50e5]
  poly( 6,:) = [-1.55e5, -3.05e5]
  poly( 7,:) = [-1.55e5, -2.80e5]
  poly( 8,:) = [-1.65e5, -3.05e5]
  poly( 9,:) = [-1.90e5, -3.25e5]
  poly(10,:) = [-2.05e5, -3.45e5]
  poly(11,:) = [-2.10e5, -3.90e5]
  poly(12,:) = [-2.15e5, -4.45e5]
  poly(13,:) = [-2.35e5, -4.40e5]
  poly(14,:) = [-2.25e5, -4.75e5]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Scott_glacier

subroutine calc_polygon_Mercer_glacier( poly)
  ! Return a polygon enveloping Mercer glacier in the Transantarctic Mountains

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Mercer_glacier'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 9,2))

  poly( 1,:) = [-3.45e5, -3.65e5]
  poly( 2,:) = [-3.30e5, -3.30e5]
  poly( 3,:) = [-3.20e5, -2.80e5]
  poly( 4,:) = [-3.10e5, -2.50e5]
  poly( 5,:) = [-2.85e5, -2.20e5]
  poly( 6,:) = [-3.10e5, -2.20e5]
  poly( 7,:) = [-3.35e5, -2.65e5]
  poly( 8,:) = [-3.40e5, -3.20e5]
  poly( 9,:) = [-3.60e5, -3.45e5]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Mercer_glacier

subroutine calc_polygon_Wilkes_basins( poly)
  ! Return a polygon enveloping basins 14, 15, and 16 of Wilkes basin, as requested for WilkesMIP

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Wilkes_basins'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 180,2))

  poly(  1,:) = [8.55E+05_dp, -2.12E+06_dp]
  poly(  2,:) = [8.19E+05_dp, -2.12E+06_dp]
  poly(  3,:) = [8.01E+05_dp, -2.10E+06_dp]
  poly(  4,:) = [7.99E+05_dp, -2.11E+06_dp]
  poly(  5,:) = [7.73E+05_dp, -2.11E+06_dp]
  poly(  6,:) = [7.44E+05_dp, -2.08E+06_dp]
  poly(  7,:) = [7.63E+05_dp, -2.07E+06_dp]
  poly(  8,:) = [7.06E+05_dp, -2.05E+06_dp]
  poly(  9,:) = [6.70E+05_dp, -1.98E+06_dp]
  poly( 10,:) = [6.49E+05_dp, -1.99E+06_dp]
  poly( 11,:) = [6.26E+05_dp, -1.98E+06_dp]
  poly( 12,:) = [6.63E+05_dp, -2.04E+06_dp]
  poly( 13,:) = [6.55E+05_dp, -2.06E+06_dp]
  poly( 14,:) = [6.08E+05_dp, -2.03E+06_dp]
  poly( 15,:) = [5.90E+05_dp, -2.04E+06_dp]
  poly( 16,:) = [5.84E+05_dp, -2.02E+06_dp]
  poly( 17,:) = [5.84E+05_dp, -2.06E+06_dp]
  poly( 18,:) = [5.35E+05_dp, -2.07E+06_dp]
  poly( 19,:) = [5.15E+05_dp, -2.05E+06_dp]
  poly( 20,:) = [4.99E+05_dp, -2.07E+06_dp]
  poly( 21,:) = [4.82E+05_dp, -2.06E+06_dp]
  poly( 22,:) = [4.93E+05_dp, -2.05E+06_dp]
  poly( 23,:) = [4.42E+05_dp, -2.05E+06_dp]
  poly( 24,:) = [4.64E+05_dp, -2.03E+06_dp]
  poly( 25,:) = [4.35E+05_dp, -2.03E+06_dp]
  poly( 26,:) = [4.38E+05_dp, -2.01E+06_dp]
  poly( 27,:) = [4.05E+05_dp, -2.02E+06_dp]
  poly( 28,:) = [3.42E+05_dp, -1.98E+06_dp]
  poly( 29,:) = [3.47E+05_dp, -2.02E+06_dp]
  poly( 30,:) = [3.20E+05_dp, -1.98E+06_dp]
  poly( 31,:) = [3.13E+05_dp, -1.96E+06_dp]
  poly( 32,:) = [3.24E+05_dp, -1.95E+06_dp]
  poly( 33,:) = [3.36E+05_dp, -1.96E+06_dp]
  poly( 34,:) = [3.46E+05_dp, -1.92E+06_dp]
  poly( 35,:) = [3.38E+05_dp, -1.89E+06_dp]
  poly( 36,:) = [3.27E+05_dp, -1.91E+06_dp]
  poly( 37,:) = [3.23E+05_dp, -1.87E+06_dp]
  poly( 38,:) = [3.91E+05_dp, -1.89E+06_dp]
  poly( 39,:) = [3.30E+05_dp, -1.86E+06_dp]
  poly( 40,:) = [3.42E+05_dp, -1.80E+06_dp]
  poly( 41,:) = [3.75E+05_dp, -1.83E+06_dp]
  poly( 42,:) = [3.71E+05_dp, -1.80E+06_dp]
  poly( 43,:) = [4.08E+05_dp, -1.80E+06_dp]
  poly( 44,:) = [4.41E+05_dp, -1.82E+06_dp]
  poly( 45,:) = [4.17E+05_dp, -1.77E+06_dp]
  poly( 46,:) = [4.08E+05_dp, -1.79E+06_dp]
  poly( 47,:) = [3.80E+05_dp, -1.78E+06_dp]
  poly( 48,:) = [4.08E+05_dp, -1.75E+06_dp]
  poly( 49,:) = [4.34E+05_dp, -1.75E+06_dp]
  poly( 50,:) = [4.41E+05_dp, -1.74E+06_dp]
  poly( 51,:) = [4.21E+05_dp, -1.72E+06_dp]
  poly( 52,:) = [4.41E+05_dp, -1.72E+06_dp]
  poly( 53,:) = [4.36E+05_dp, -1.70E+06_dp]
  poly( 54,:) = [4.51E+05_dp, -1.71E+06_dp]
  poly( 55,:) = [4.61E+05_dp, -1.76E+06_dp]
  poly( 56,:) = [4.87E+05_dp, -1.76E+06_dp]
  poly( 57,:) = [4.43E+05_dp, -1.69E+06_dp]
  poly( 58,:) = [4.47E+05_dp, -1.68E+06_dp]
  poly( 59,:) = [4.65E+05_dp, -1.69E+06_dp]
  poly( 60,:) = [4.72E+05_dp, -1.69E+06_dp]
  poly( 61,:) = [4.22E+05_dp, -1.62E+06_dp]
  poly( 62,:) = [4.57E+05_dp, -1.63E+06_dp]
  poly( 63,:) = [4.70E+05_dp, -1.61E+06_dp]
  poly( 64,:) = [4.88E+05_dp, -1.64E+06_dp]
  poly( 65,:) = [5.02E+05_dp, -1.58E+06_dp]
  poly( 66,:) = [4.83E+05_dp, -1.54E+06_dp]
  poly( 67,:) = [5.07E+05_dp, -1.54E+06_dp]
  poly( 68,:) = [5.31E+05_dp, -1.51E+06_dp]
  poly( 69,:) = [4.87E+05_dp, -1.52E+06_dp]
  poly( 70,:) = [4.85E+05_dp, -1.50E+06_dp]
  poly( 71,:) = [4.55E+05_dp, -1.48E+06_dp]
  poly( 72,:) = [4.64E+05_dp, -1.43E+06_dp]
  poly( 73,:) = [4.45E+05_dp, -1.44E+06_dp]
  poly( 74,:) = [4.35E+05_dp, -1.41E+06_dp]
  poly( 75,:) = [4.40E+05_dp, -1.39E+06_dp]
  poly( 76,:) = [4.22E+05_dp, -1.38E+06_dp]
  poly( 77,:) = [4.33E+05_dp, -1.36E+06_dp]
  poly( 78,:) = [4.11E+05_dp, -1.35E+06_dp]
  poly( 79,:) = [3.82E+05_dp, -1.32E+06_dp]
  poly( 80,:) = [3.81E+05_dp, -1.28E+06_dp]
  poly( 81,:) = [3.61E+05_dp, -1.28E+06_dp]
  poly( 82,:) = [3.71E+05_dp, -1.27E+06_dp]
  poly( 83,:) = [3.52E+05_dp, -1.27E+06_dp]
  poly( 84,:) = [3.57E+05_dp, -1.23E+06_dp]
  poly( 85,:) = [3.36E+05_dp, -1.23E+06_dp]
  poly( 86,:) = [3.28E+05_dp, -1.26E+06_dp]
  poly( 87,:) = [3.24E+05_dp, -1.22E+06_dp]
  poly( 88,:) = [3.57E+05_dp, -1.20E+06_dp]
  poly( 89,:) = [3.74E+05_dp, -1.20E+06_dp]
  poly( 90,:) = [3.92E+05_dp, -1.23E+06_dp]
  poly( 91,:) = [5.41E+05_dp, -1.25E+06_dp]
  poly( 92,:) = [5.62E+05_dp, -1.28E+06_dp]
  poly( 93,:) = [6.05E+05_dp, -1.25E+06_dp]
  poly( 94,:) = [7.12E+05_dp, -1.24E+06_dp]
  poly( 95,:) = [7.37E+05_dp, -1.22E+06_dp]
  poly( 96,:) = [7.58E+05_dp, -1.23E+06_dp]
  poly( 97,:) = [7.94E+05_dp, -1.18E+06_dp]
  poly( 98,:) = [8.26E+05_dp, -1.19E+06_dp]
  poly( 99,:) = [8.40E+05_dp, -1.15E+06_dp]
  poly(100,:) = [8.85E+05_dp, -1.15E+06_dp]
  poly(101,:) = [9.27E+05_dp, -1.17E+06_dp]
  poly(102,:) = [9.64E+05_dp, -1.13E+06_dp]
  poly(103,:) = [9.89E+05_dp, -1.15E+06_dp]
  poly(104,:) = [9.89E+05_dp, -1.13E+06_dp]
  poly(105,:) = [1.08E+06_dp, -1.13E+06_dp]
  poly(106,:) = [1.11E+06_dp, -1.10E+06_dp]
  poly(107,:) = [1.22E+06_dp, -1.17E+06_dp]
  poly(108,:) = [1.27E+06_dp, -1.18E+06_dp]
  poly(109,:) = [1.32E+06_dp, -1.15E+06_dp]
  poly(110,:) = [1.37E+06_dp, -1.18E+06_dp]
  poly(111,:) = [1.43E+06_dp, -1.15E+06_dp]
  poly(112,:) = [1.47E+06_dp, -1.22E+06_dp]
  poly(113,:) = [1.45E+06_dp, -1.31E+06_dp]
  poly(114,:) = [1.50E+06_dp, -1.34E+06_dp]
  poly(115,:) = [1.50E+06_dp, -1.37E+06_dp]
  poly(116,:) = [1.55E+06_dp, -1.44E+06_dp]
  poly(117,:) = [1.63E+06_dp, -1.60E+06_dp]
  poly(118,:) = [1.79E+06_dp, -1.67E+06_dp]
  poly(119,:) = [1.88E+06_dp, -1.65E+06_dp]
  poly(120,:) = [2.00E+06_dp, -1.68E+06_dp]
  poly(121,:) = [1.99E+06_dp, -1.72E+06_dp]
  poly(122,:) = [1.95E+06_dp, -1.74E+06_dp]
  poly(123,:) = [1.91E+06_dp, -1.81E+06_dp]
  poly(124,:) = [1.86E+06_dp, -1.83E+06_dp]
  poly(125,:) = [1.85E+06_dp, -1.81E+06_dp]
  poly(126,:) = [1.84E+06_dp, -1.85E+06_dp]
  poly(127,:) = [1.85E+06_dp, -1.87E+06_dp]
  poly(128,:) = [1.82E+06_dp, -1.86E+06_dp]
  poly(129,:) = [1.81E+06_dp, -1.89E+06_dp]
  poly(130,:) = [1.77E+06_dp, -1.88E+06_dp]
  poly(131,:) = [1.79E+06_dp, -1.90E+06_dp]
  poly(132,:) = [1.75E+06_dp, -1.93E+06_dp]
  poly(133,:) = [1.72E+06_dp, -1.92E+06_dp]
  poly(134,:) = [1.71E+06_dp, -1.94E+06_dp]
  poly(135,:) = [1.67E+06_dp, -1.96E+06_dp]
  poly(136,:) = [1.65E+06_dp, -1.95E+06_dp]
  poly(137,:) = [1.58E+06_dp, -2.01E+06_dp]
  poly(138,:) = [1.54E+06_dp, -2.01E+06_dp]
  poly(139,:) = [1.52E+06_dp, -2.05E+06_dp]
  poly(140,:) = [1.50E+06_dp, -2.05E+06_dp]
  poly(141,:) = [1.49E+06_dp, -2.03E+06_dp]
  poly(142,:) = [1.47E+06_dp, -2.04E+06_dp]
  poly(143,:) = [1.47E+06_dp, -2.06E+06_dp]
  poly(144,:) = [1.46E+06_dp, -2.06E+06_dp]
  poly(145,:) = [1.46E+06_dp, -2.01E+06_dp]
  poly(146,:) = [1.43E+06_dp, -1.96E+06_dp]
  poly(147,:) = [1.41E+06_dp, -2.03E+06_dp]
  poly(148,:) = [1.38E+06_dp, -2.02E+06_dp]
  poly(149,:) = [1.38E+06_dp, -2.05E+06_dp]
  poly(150,:) = [1.35E+06_dp, -2.03E+06_dp]
  poly(151,:) = [1.35E+06_dp, -2.04E+06_dp]
  poly(152,:) = [1.34E+06_dp, -2.02E+06_dp]
  poly(153,:) = [1.33E+06_dp, -2.03E+06_dp]
  poly(154,:) = [1.32E+06_dp, -2.00E+06_dp]
  poly(155,:) = [1.31E+06_dp, -2.03E+06_dp]
  poly(156,:) = [1.31E+06_dp, -2.00E+06_dp]
  poly(157,:) = [1.30E+06_dp, -2.00E+06_dp]
  poly(158,:) = [1.29E+06_dp, -1.99E+06_dp]
  poly(159,:) = [1.21E+06_dp, -2.05E+06_dp]
  poly(160,:) = [1.19E+06_dp, -2.04E+06_dp]
  poly(161,:) = [1.20E+06_dp, -2.05E+06_dp]
  poly(162,:) = [1.17E+06_dp, -2.06E+06_dp]
  poly(163,:) = [1.16E+06_dp, -2.08E+06_dp]
  poly(164,:) = [1.13E+06_dp, -2.04E+06_dp]
  poly(165,:) = [1.10E+06_dp, -2.04E+06_dp]
  poly(166,:) = [1.11E+06_dp, -2.06E+06_dp]
  poly(167,:) = [1.10E+06_dp, -2.07E+06_dp]
  poly(168,:) = [1.05E+06_dp, -2.07E+06_dp]
  poly(169,:) = [1.03E+06_dp, -2.09E+06_dp]
  poly(170,:) = [1.05E+06_dp, -2.14E+06_dp]
  poly(171,:) = [1.01E+06_dp, -2.14E+06_dp]
  poly(172,:) = [1.00E+06_dp, -2.10E+06_dp]
  poly(173,:) = [9.95E+05_dp, -2.12E+06_dp]
  poly(174,:) = [9.51E+05_dp, -2.09E+06_dp]
  poly(175,:) = [9.51E+05_dp, -2.11E+06_dp]
  poly(176,:) = [9.17E+05_dp, -2.09E+06_dp]
  poly(177,:) = [8.89E+05_dp, -2.11E+06_dp]
  poly(178,:) = [8.67E+05_dp, -2.08E+06_dp]
  poly(179,:) = [8.72E+05_dp, -2.11E+06_dp]
  poly(180,:) = [8.55E+05_dp, -2.12E+06_dp]

  
  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Wilkes_basins

subroutine calc_polygon_Institute_basin( poly)
  ! Return a polygon enveloping Scott glacier in the Transantarctic Mountains

  ! In/output variables:
  real(dp), dimension(:,:), allocatable, intent(out) :: poly

  ! Local variables:
  character(len=1024), parameter :: routine_name = 'calc_polygon_Institute_basin'

  ! Add routine to path
  call init_routine( routine_name)

  allocate( poly( 499,2))

  poly(  1,:) = [-1.60E+06_dp, 7.53E+05_dp]
  poly(  2,:) = [-1.59E+06_dp, 7.56E+05_dp]
  poly(  3,:) = [-1.58E+06_dp, 7.51E+05_dp]
  poly(  4,:) = [-1.56E+06_dp, 7.53E+05_dp]
  poly(  5,:) = [-1.55E+06_dp, 7.60E+05_dp]
  poly(  6,:) = [-1.53E+06_dp, 7.65E+05_dp]
  poly(  7,:) = [-1.52E+06_dp, 7.70E+05_dp]
  poly(  8,:) = [-1.51E+06_dp, 7.71E+05_dp]
  poly(  9,:) = [-1.50E+06_dp, 7.68E+05_dp]
  poly( 10,:) = [-1.49E+06_dp, 7.75E+05_dp]
  poly( 11,:) = [-1.48E+06_dp, 7.72E+05_dp]
  poly( 12,:) = [-1.48E+06_dp, 7.73E+05_dp]
  poly( 13,:) = [-1.47E+06_dp, 7.67E+05_dp]
  poly( 14,:) = [-1.46E+06_dp, 7.83E+05_dp]
  poly( 15,:) = [-1.45E+06_dp, 7.71E+05_dp]
  poly( 16,:) = [-1.45E+06_dp, 7.58E+05_dp]
  poly( 17,:) = [-1.46E+06_dp, 7.57E+05_dp]
  poly( 18,:) = [-1.48E+06_dp, 7.63E+05_dp]
  poly( 19,:) = [-1.48E+06_dp, 7.62E+05_dp]
  poly( 20,:) = [-1.49E+06_dp, 7.63E+05_dp]
  poly( 21,:) = [-1.49E+06_dp, 7.58E+05_dp]
  poly( 22,:) = [-1.50E+06_dp, 7.55E+05_dp]
  poly( 23,:) = [-1.49E+06_dp, 7.51E+05_dp]
  poly( 24,:) = [-1.48E+06_dp, 7.51E+05_dp]
  poly( 25,:) = [-1.47E+06_dp, 7.48E+05_dp]
  poly( 26,:) = [-1.47E+06_dp, 7.46E+05_dp]
  poly( 27,:) = [-1.47E+06_dp, 7.41E+05_dp]
  poly( 28,:) = [-1.48E+06_dp, 7.37E+05_dp]
  poly( 29,:) = [-1.48E+06_dp, 7.34E+05_dp]
  poly( 30,:) = [-1.47E+06_dp, 7.33E+05_dp]
  poly( 31,:) = [-1.47E+06_dp, 7.18E+05_dp]
  poly( 32,:) = [-1.47E+06_dp, 7.16E+05_dp]
  poly( 33,:) = [-1.47E+06_dp, 7.17E+05_dp]
  poly( 34,:) = [-1.47E+06_dp, 7.22E+05_dp]
  poly( 35,:) = [-1.46E+06_dp, 7.34E+05_dp]
  poly( 36,:) = [-1.46E+06_dp, 7.36E+05_dp]
  poly( 37,:) = [-1.46E+06_dp, 7.37E+05_dp]
  poly( 38,:) = [-1.45E+06_dp, 7.35E+05_dp]
  poly( 39,:) = [-1.45E+06_dp, 7.36E+05_dp]
  poly( 40,:) = [-1.45E+06_dp, 7.36E+05_dp]
  poly( 41,:) = [-1.45E+06_dp, 7.27E+05_dp]
  poly( 42,:) = [-1.45E+06_dp, 7.23E+05_dp]
  poly( 43,:) = [-1.45E+06_dp, 7.15E+05_dp]
  poly( 44,:) = [-1.45E+06_dp, 7.05E+05_dp]
  poly( 45,:) = [-1.45E+06_dp, 6.94E+05_dp]
  poly( 46,:) = [-1.44E+06_dp, 6.90E+05_dp]
  poly( 47,:) = [-1.44E+06_dp, 6.97E+05_dp]
  poly( 48,:) = [-1.44E+06_dp, 6.99E+05_dp]
  poly( 49,:) = [-1.44E+06_dp, 7.06E+05_dp]
  poly( 50,:) = [-1.43E+06_dp, 7.20E+05_dp]
  poly( 51,:) = [-1.43E+06_dp, 7.24E+05_dp]
  poly( 52,:) = [-1.42E+06_dp, 7.21E+05_dp]
  poly( 53,:) = [-1.41E+06_dp, 6.87E+05_dp]
  poly( 54,:) = [-1.41E+06_dp, 6.59E+05_dp]
  poly( 55,:) = [-1.41E+06_dp, 6.44E+05_dp]
  poly( 56,:) = [-1.40E+06_dp, 6.37E+05_dp]
  poly( 57,:) = [-1.40E+06_dp, 6.15E+05_dp]
  poly( 58,:) = [-1.40E+06_dp, 6.03E+05_dp]
  poly( 59,:) = [-1.40E+06_dp, 5.86E+05_dp]
  poly( 60,:) = [-1.40E+06_dp, 5.65E+05_dp]
  poly( 61,:) = [-1.40E+06_dp, 5.60E+05_dp]
  poly( 62,:) = [-1.40E+06_dp, 5.61E+05_dp]
  poly( 63,:) = [-1.40E+06_dp, 5.56E+05_dp]
  poly( 64,:) = [-1.40E+06_dp, 5.53E+05_dp]
  poly( 65,:) = [-1.40E+06_dp, 5.54E+05_dp]
  poly( 66,:) = [-1.39E+06_dp, 5.60E+05_dp]
  poly( 67,:) = [-1.39E+06_dp, 5.47E+05_dp]
  poly( 68,:) = [-1.40E+06_dp, 5.34E+05_dp]
  poly( 69,:) = [-1.40E+06_dp, 5.32E+05_dp]
  poly( 70,:) = [-1.40E+06_dp, 5.14E+05_dp]
  poly( 71,:) = [-1.40E+06_dp, 5.13E+05_dp]
  poly( 72,:) = [-1.39E+06_dp, 5.03E+05_dp]
  poly( 73,:) = [-1.39E+06_dp, 5.00E+05_dp]
  poly( 74,:) = [-1.39E+06_dp, 5.02E+05_dp]
  poly( 75,:) = [-1.38E+06_dp, 5.03E+05_dp]
  poly( 76,:) = [-1.38E+06_dp, 5.01E+05_dp]
  poly( 77,:) = [-1.38E+06_dp, 4.96E+05_dp]
  poly( 78,:) = [-1.37E+06_dp, 4.97E+05_dp]
  poly( 79,:) = [-1.37E+06_dp, 4.92E+05_dp]
  poly( 80,:) = [-1.37E+06_dp, 4.92E+05_dp]
  poly( 81,:) = [-1.37E+06_dp, 4.84E+05_dp]
  poly( 82,:) = [-1.37E+06_dp, 4.83E+05_dp]
  poly( 83,:) = [-1.36E+06_dp, 4.86E+05_dp]
  poly( 84,:) = [-1.36E+06_dp, 4.82E+05_dp]
  poly( 85,:) = [-1.37E+06_dp, 4.54E+05_dp]
  poly( 86,:) = [-1.38E+06_dp, 4.48E+05_dp]
  poly( 87,:) = [-1.38E+06_dp, 4.50E+05_dp]
  poly( 88,:) = [-1.38E+06_dp, 4.46E+05_dp]
  poly( 89,:) = [-1.38E+06_dp, 4.41E+05_dp]
  poly( 90,:) = [-1.39E+06_dp, 4.38E+05_dp]
  poly( 91,:) = [-1.38E+06_dp, 4.36E+05_dp]
  poly( 92,:) = [-1.39E+06_dp, 4.35E+05_dp]
  poly( 93,:) = [-1.39E+06_dp, 4.32E+05_dp]
  poly( 94,:) = [-1.39E+06_dp, 4.27E+05_dp]
  poly( 95,:) = [-1.39E+06_dp, 4.24E+05_dp]
  poly( 96,:) = [-1.39E+06_dp, 4.22E+05_dp]
  poly( 97,:) = [-1.39E+06_dp, 4.23E+05_dp]
  poly( 98,:) = [-1.39E+06_dp, 4.18E+05_dp]
  poly( 99,:) = [-1.40E+06_dp, 4.05E+05_dp]
  poly(100,:) = [-1.39E+06_dp, 4.03E+05_dp]
  poly(101,:) = [-1.40E+06_dp, 3.99E+05_dp]
  poly(102,:) = [-1.40E+06_dp, 3.93E+05_dp]
  poly(103,:) = [-1.42E+06_dp, 3.80E+05_dp]
  poly(104,:) = [-1.42E+06_dp, 3.77E+05_dp]
  poly(105,:) = [-1.42E+06_dp, 3.73E+05_dp]
  poly(106,:) = [-1.42E+06_dp, 3.72E+05_dp]
  poly(107,:) = [-1.42E+06_dp, 3.70E+05_dp]
  poly(108,:) = [-1.42E+06_dp, 3.64E+05_dp]
  poly(109,:) = [-1.42E+06_dp, 3.63E+05_dp]
  poly(110,:) = [-1.41E+06_dp, 3.64E+05_dp]
  poly(111,:) = [-1.41E+06_dp, 3.62E+05_dp]
  poly(112,:) = [-1.41E+06_dp, 3.59E+05_dp]
  poly(113,:) = [-1.41E+06_dp, 3.56E+05_dp]
  poly(114,:) = [-1.42E+06_dp, 3.55E+05_dp]
  poly(115,:) = [-1.42E+06_dp, 3.56E+05_dp]
  poly(116,:) = [-1.43E+06_dp, 3.54E+05_dp]
  poly(117,:) = [-1.43E+06_dp, 3.50E+05_dp]
  poly(118,:) = [-1.43E+06_dp, 3.45E+05_dp]
  poly(119,:) = [-1.43E+06_dp, 3.39E+05_dp]
  poly(120,:) = [-1.45E+06_dp, 3.37E+05_dp]
  poly(121,:) = [-1.47E+06_dp, 3.28E+05_dp]
  poly(122,:) = [-1.48E+06_dp, 3.35E+05_dp]
  poly(123,:) = [-1.50E+06_dp, 3.33E+05_dp]
  poly(124,:) = [-1.51E+06_dp, 3.27E+05_dp]
  poly(125,:) = [-1.51E+06_dp, 3.20E+05_dp]
  poly(126,:) = [-1.51E+06_dp, 3.16E+05_dp]
  poly(127,:) = [-1.50E+06_dp, 3.16E+05_dp]
  poly(128,:) = [-1.49E+06_dp, 3.19E+05_dp]
  poly(129,:) = [-1.49E+06_dp, 3.18E+05_dp]
  poly(130,:) = [-1.49E+06_dp, 3.11E+05_dp]
  poly(131,:) = [-1.48E+06_dp, 3.07E+05_dp]
  poly(132,:) = [-1.48E+06_dp, 3.01E+05_dp]
  poly(133,:) = [-1.47E+06_dp, 2.98E+05_dp]
  poly(134,:) = [-1.46E+06_dp, 3.06E+05_dp]
  poly(135,:) = [-1.46E+06_dp, 3.07E+05_dp]
  poly(136,:) = [-1.45E+06_dp, 3.12E+05_dp]
  poly(137,:) = [-1.45E+06_dp, 3.11E+05_dp]
  poly(138,:) = [-1.44E+06_dp, 3.15E+05_dp]
  poly(139,:) = [-1.44E+06_dp, 3.13E+05_dp]
  poly(140,:) = [-1.42E+06_dp, 3.14E+05_dp]
  poly(141,:) = [-1.41E+06_dp, 3.12E+05_dp]
  poly(142,:) = [-1.38E+06_dp, 3.22E+05_dp]
  poly(143,:) = [-1.38E+06_dp, 3.22E+05_dp]
  poly(144,:) = [-1.39E+06_dp, 3.20E+05_dp]
  poly(145,:) = [-1.39E+06_dp, 3.16E+05_dp]
  poly(146,:) = [-1.39E+06_dp, 3.14E+05_dp]
  poly(147,:) = [-1.39E+06_dp, 3.14E+05_dp]
  poly(148,:) = [-1.37E+06_dp, 3.19E+05_dp]
  poly(149,:) = [-1.37E+06_dp, 3.25E+05_dp]
  poly(150,:) = [-1.36E+06_dp, 3.28E+05_dp]
  poly(151,:) = [-1.34E+06_dp, 3.29E+05_dp]
  poly(152,:) = [-1.33E+06_dp, 3.25E+05_dp]
  poly(153,:) = [-1.32E+06_dp, 3.26E+05_dp]
  poly(154,:) = [-1.32E+06_dp, 3.31E+05_dp]
  poly(155,:) = [-1.33E+06_dp, 3.34E+05_dp]
  poly(156,:) = [-1.33E+06_dp, 3.36E+05_dp]
  poly(157,:) = [-1.32E+06_dp, 3.43E+05_dp]
  poly(158,:) = [-1.32E+06_dp, 3.53E+05_dp]
  poly(159,:) = [-1.31E+06_dp, 3.78E+05_dp]
  poly(160,:) = [-1.29E+06_dp, 4.01E+05_dp]
  poly(161,:) = [-1.29E+06_dp, 4.00E+05_dp]
  poly(162,:) = [-1.29E+06_dp, 4.01E+05_dp]
  poly(163,:) = [-1.28E+06_dp, 3.97E+05_dp]
  poly(164,:) = [-1.26E+06_dp, 3.72E+05_dp]
  poly(165,:) = [-1.25E+06_dp, 3.72E+05_dp]
  poly(166,:) = [-1.25E+06_dp, 3.69E+05_dp]
  poly(167,:) = [-1.25E+06_dp, 3.65E+05_dp]
  poly(168,:) = [-1.25E+06_dp, 3.63E+05_dp]
  poly(169,:) = [-1.24E+06_dp, 3.52E+05_dp]
  poly(170,:) = [-1.25E+06_dp, 3.36E+05_dp]
  poly(171,:) = [-1.24E+06_dp, 3.31E+05_dp]
  poly(172,:) = [-1.25E+06_dp, 3.31E+05_dp]
  poly(173,:) = [-1.26E+06_dp, 3.23E+05_dp]
  poly(174,:) = [-1.28E+06_dp, 3.00E+05_dp]
  poly(175,:) = [-1.28E+06_dp, 2.88E+05_dp]
  poly(176,:) = [-1.29E+06_dp, 2.82E+05_dp]
  poly(177,:) = [-1.29E+06_dp, 2.71E+05_dp]
  poly(178,:) = [-1.29E+06_dp, 2.69E+05_dp]
  poly(179,:) = [-1.29E+06_dp, 2.60E+05_dp]
  poly(180,:) = [-1.30E+06_dp, 2.56E+05_dp]
  poly(181,:) = [-1.31E+06_dp, 2.52E+05_dp]
  poly(182,:) = [-1.31E+06_dp, 2.45E+05_dp]
  poly(183,:) = [-1.31E+06_dp, 2.28E+05_dp]
  poly(184,:) = [-1.32E+06_dp, 2.18E+05_dp]
  poly(185,:) = [-1.32E+06_dp, 2.16E+05_dp]
  poly(186,:) = [-1.32E+06_dp, 2.15E+05_dp]
  poly(187,:) = [-1.31E+06_dp, 2.08E+05_dp]
  poly(188,:) = [-1.30E+06_dp, 2.10E+05_dp]
  poly(189,:) = [-1.30E+06_dp, 2.07E+05_dp]
  poly(190,:) = [-1.31E+06_dp, 2.04E+05_dp]
  poly(191,:) = [-1.31E+06_dp, 1.98E+05_dp]
  poly(192,:) = [-1.31E+06_dp, 1.95E+05_dp]
  poly(193,:) = [-1.31E+06_dp, 1.95E+05_dp]
  poly(194,:) = [-1.29E+06_dp, 2.09E+05_dp]
  poly(195,:) = [-1.25E+06_dp, 2.64E+05_dp]
  poly(196,:) = [-1.24E+06_dp, 2.71E+05_dp]
  poly(197,:) = [-1.24E+06_dp, 2.74E+05_dp]
  poly(198,:) = [-1.23E+06_dp, 2.76E+05_dp]
  poly(199,:) = [-1.22E+06_dp, 2.74E+05_dp]
  poly(200,:) = [-1.21E+06_dp, 2.75E+05_dp]
  poly(201,:) = [-1.21E+06_dp, 2.74E+05_dp]
  poly(202,:) = [-1.20E+06_dp, 2.63E+05_dp]
  poly(203,:) = [-1.20E+06_dp, 2.62E+05_dp]
  poly(204,:) = [-1.19E+06_dp, 2.40E+05_dp]
  poly(205,:) = [-1.19E+06_dp, 2.23E+05_dp]
  poly(206,:) = [-1.19E+06_dp, 2.15E+05_dp]
  poly(207,:) = [-1.20E+06_dp, 2.06E+05_dp]
  poly(208,:) = [-1.20E+06_dp, 2.02E+05_dp]
  poly(209,:) = [-1.22E+06_dp, 1.86E+05_dp]
  poly(210,:) = [-1.24E+06_dp, 1.74E+05_dp]
  poly(211,:) = [-1.27E+06_dp, 1.59E+05_dp]
  poly(212,:) = [-1.27E+06_dp, 1.57E+05_dp]
  poly(213,:) = [-1.26E+06_dp, 1.54E+05_dp]
  poly(214,:) = [-1.25E+06_dp, 1.53E+05_dp]
  poly(215,:) = [-1.25E+06_dp, 1.55E+05_dp]
  poly(216,:) = [-1.25E+06_dp, 1.51E+05_dp]
  poly(217,:) = [-1.26E+06_dp, 1.47E+05_dp]
  poly(218,:) = [-1.26E+06_dp, 1.41E+05_dp]
  poly(219,:) = [-1.26E+06_dp, 1.37E+05_dp]
  poly(220,:) = [-1.25E+06_dp, 1.35E+05_dp]
  poly(221,:) = [-1.24E+06_dp, 1.36E+05_dp]
  poly(222,:) = [-1.23E+06_dp, 1.40E+05_dp]
  poly(223,:) = [-1.22E+06_dp, 1.40E+05_dp]
  poly(224,:) = [-1.22E+06_dp, 1.43E+05_dp]
  poly(225,:) = [-1.21E+06_dp, 1.46E+05_dp]
  poly(226,:) = [-1.20E+06_dp, 1.51E+05_dp]
  poly(227,:) = [-1.20E+06_dp, 1.53E+05_dp]
  poly(228,:) = [-1.20E+06_dp, 1.53E+05_dp]
  poly(229,:) = [-1.19E+06_dp, 1.57E+05_dp]
  poly(230,:) = [-1.18E+06_dp, 1.68E+05_dp]
  poly(231,:) = [-1.17E+06_dp, 1.76E+05_dp]
  poly(232,:) = [-1.16E+06_dp, 1.82E+05_dp]
  poly(233,:) = [-1.16E+06_dp, 1.77E+05_dp]
  poly(234,:) = [-1.14E+06_dp, 1.76E+05_dp]
  poly(235,:) = [-1.13E+06_dp, 1.81E+05_dp]
  poly(236,:) = [-1.12E+06_dp, 1.83E+05_dp]
  poly(237,:) = [-1.12E+06_dp, 1.87E+05_dp]
  poly(238,:) = [-1.13E+06_dp, 1.91E+05_dp]
  poly(239,:) = [-1.14E+06_dp, 1.90E+05_dp]
  poly(240,:) = [-1.15E+06_dp, 1.94E+05_dp]
  poly(241,:) = [-1.16E+06_dp, 2.00E+05_dp]
  poly(242,:) = [-1.15E+06_dp, 2.18E+05_dp]
  poly(243,:) = [-1.14E+06_dp, 2.58E+05_dp]
  poly(244,:) = [-1.13E+06_dp, 2.69E+05_dp]
  poly(245,:) = [-1.13E+06_dp, 2.71E+05_dp]
  poly(246,:) = [-1.12E+06_dp, 2.75E+05_dp]
  poly(247,:) = [-1.11E+06_dp, 2.75E+05_dp]
  poly(248,:) = [-1.08E+06_dp, 2.64E+05_dp]
  poly(249,:) = [-1.07E+06_dp, 2.54E+05_dp]
  poly(250,:) = [-1.06E+06_dp, 2.39E+05_dp]
  poly(251,:) = [-1.06E+06_dp, 2.22E+05_dp]
  poly(252,:) = [-1.07E+06_dp, 2.02E+05_dp]
  poly(253,:) = [-1.07E+06_dp, 1.96E+05_dp]
  poly(254,:) = [-1.07E+06_dp, 1.95E+05_dp]
  poly(255,:) = [-1.06E+06_dp, 2.03E+05_dp]
  poly(256,:) = [-1.06E+06_dp, 2.05E+05_dp]
  poly(257,:) = [-1.05E+06_dp, 2.17E+05_dp]
  poly(258,:) = [-1.05E+06_dp, 2.28E+05_dp]
  poly(259,:) = [-1.04E+06_dp, 2.38E+05_dp]
  poly(260,:) = [-1.05E+06_dp, 2.55E+05_dp]
  poly(261,:) = [-1.04E+06_dp, 2.51E+05_dp]
  poly(262,:) = [-1.04E+06_dp, 2.53E+05_dp]
  poly(263,:) = [-1.03E+06_dp, 2.54E+05_dp]
  poly(264,:) = [-1.01E+06_dp, 2.58E+05_dp]
  poly(265,:) = [-1.00E+06_dp, 2.58E+05_dp]
  poly(266,:) = [-1.00E+06_dp, 2.57E+05_dp]
  poly(267,:) = [-1.00E+06_dp, 2.54E+05_dp]
  poly(268,:) = [-1.02E+06_dp, 2.52E+05_dp]
  poly(269,:) = [-1.03E+06_dp, 2.48E+05_dp]
  poly(270,:) = [-1.02E+06_dp, 2.42E+05_dp]
  poly(271,:) = [-1.01E+06_dp, 2.45E+05_dp]
  poly(272,:) = [-9.97E+05_dp, 2.50E+05_dp]
  poly(273,:) = [-9.87E+05_dp, 2.49E+05_dp]
  poly(274,:) = [-9.81E+05_dp, 2.50E+05_dp]
  poly(275,:) = [-9.71E+05_dp, 2.59E+05_dp]
  poly(276,:) = [-9.71E+05_dp, 2.51E+05_dp]
  poly(277,:) = [-9.69E+05_dp, 2.45E+05_dp]
  poly(278,:) = [-9.66E+05_dp, 2.44E+05_dp]
  poly(279,:) = [-9.61E+05_dp, 2.47E+05_dp]
  poly(280,:) = [-9.57E+05_dp, 2.56E+05_dp]
  poly(281,:) = [-9.56E+05_dp, 2.59E+05_dp]
  poly(282,:) = [-9.60E+05_dp, 2.64E+05_dp]
  poly(283,:) = [-9.60E+05_dp, 2.69E+05_dp]
  poly(284,:) = [-9.62E+05_dp, 2.71E+05_dp]
  poly(285,:) = [-9.41E+05_dp, 2.85E+05_dp]
  poly(286,:) = [-9.37E+05_dp, 2.91E+05_dp]
  poly(287,:) = [-9.38E+05_dp, 2.94E+05_dp]
  poly(288,:) = [-9.49E+05_dp, 2.97E+05_dp]
  poly(289,:) = [-9.58E+05_dp, 3.03E+05_dp]
  poly(290,:) = [-9.56E+05_dp, 3.19E+05_dp]
  poly(291,:) = [-9.60E+05_dp, 3.26E+05_dp]
  poly(292,:) = [-9.65E+05_dp, 3.30E+05_dp]
  poly(293,:) = [-9.60E+05_dp, 3.38E+05_dp]
  poly(294,:) = [-9.56E+05_dp, 3.38E+05_dp]
  poly(295,:) = [-9.44E+05_dp, 3.34E+05_dp]
  poly(296,:) = [-9.33E+05_dp, 3.32E+05_dp]
  poly(297,:) = [-9.26E+05_dp, 3.38E+05_dp]
  poly(298,:) = [-9.25E+05_dp, 3.41E+05_dp]
  poly(299,:) = [-9.27E+05_dp, 3.45E+05_dp]
  poly(300,:) = [-9.25E+05_dp, 3.46E+05_dp]
  poly(301,:) = [-9.27E+05_dp, 3.47E+05_dp]
  poly(302,:) = [-9.15E+05_dp, 3.52E+05_dp]
  poly(303,:) = [-9.09E+05_dp, 3.59E+05_dp]
  poly(304,:) = [-9.02E+05_dp, 3.60E+05_dp]
  poly(305,:) = [-8.70E+05_dp, 3.76E+05_dp]
  poly(306,:) = [-8.49E+05_dp, 3.90E+05_dp]
  poly(307,:) = [-8.42E+05_dp, 3.98E+05_dp]
  poly(308,:) = [-8.24E+05_dp, 4.07E+05_dp]
  poly(309,:) = [-8.18E+05_dp, 4.21E+05_dp]
  poly(310,:) = [-8.13E+05_dp, 4.26E+05_dp]
  poly(311,:) = [-8.10E+05_dp, 4.27E+05_dp]
  poly(312,:) = [-8.01E+05_dp, 4.24E+05_dp]
  poly(313,:) = [-7.99E+05_dp, 4.17E+05_dp]
  poly(314,:) = [-8.07E+05_dp, 4.00E+05_dp]
  poly(315,:) = [-8.10E+05_dp, 4.01E+05_dp]
  poly(316,:) = [-8.18E+05_dp, 3.88E+05_dp]
  poly(317,:) = [-8.19E+05_dp, 3.84E+05_dp]
  poly(318,:) = [-8.19E+05_dp, 3.79E+05_dp]
  poly(319,:) = [-8.18E+05_dp, 3.77E+05_dp]
  poly(320,:) = [-8.09E+05_dp, 3.73E+05_dp]
  poly(321,:) = [-8.06E+05_dp, 3.79E+05_dp]
  poly(322,:) = [-8.06E+05_dp, 3.98E+05_dp]
  poly(323,:) = [-7.94E+05_dp, 3.95E+05_dp]
  poly(324,:) = [-7.94E+05_dp, 3.80E+05_dp]
  poly(325,:) = [-7.99E+05_dp, 3.69E+05_dp]
  poly(326,:) = [-8.08E+05_dp, 3.62E+05_dp]
  poly(327,:) = [-8.07E+05_dp, 3.56E+05_dp]
  poly(328,:) = [-8.00E+05_dp, 3.53E+05_dp]
  poly(329,:) = [-7.90E+05_dp, 3.61E+05_dp]
  poly(330,:) = [-7.86E+05_dp, 3.69E+05_dp]
  poly(331,:) = [-7.89E+05_dp, 3.58E+05_dp]
  poly(332,:) = [-7.89E+05_dp, 3.54E+05_dp]
  poly(333,:) = [-7.78E+05_dp, 3.43E+05_dp]
  poly(334,:) = [-7.71E+05_dp, 3.42E+05_dp]
  poly(335,:) = [-7.65E+05_dp, 3.45E+05_dp]
  poly(336,:) = [-7.64E+05_dp, 3.51E+05_dp]
  poly(337,:) = [-7.62E+05_dp, 3.52E+05_dp]
  poly(338,:) = [-7.54E+05_dp, 3.50E+05_dp]
  poly(339,:) = [-7.49E+05_dp, 3.56E+05_dp]
  poly(340,:) = [-7.46E+05_dp, 3.74E+05_dp]
  poly(341,:) = [-7.48E+05_dp, 3.78E+05_dp]
  poly(342,:) = [-7.48E+05_dp, 3.82E+05_dp]
  poly(343,:) = [-7.46E+05_dp, 3.90E+05_dp]
  poly(344,:) = [-7.43E+05_dp, 3.92E+05_dp]
  poly(345,:) = [-7.43E+05_dp, 4.14E+05_dp]
  poly(346,:) = [-7.39E+05_dp, 4.17E+05_dp]
  poly(347,:) = [-7.32E+05_dp, 4.07E+05_dp]
  poly(348,:) = [-7.28E+05_dp, 3.91E+05_dp]
  poly(349,:) = [-7.25E+05_dp, 3.86E+05_dp]
  poly(350,:) = [-7.23E+05_dp, 3.78E+05_dp]
  poly(351,:) = [-7.25E+05_dp, 3.73E+05_dp]
  poly(352,:) = [-7.30E+05_dp, 3.70E+05_dp]
  poly(353,:) = [-7.30E+05_dp, 3.67E+05_dp]
  poly(354,:) = [-7.28E+05_dp, 3.62E+05_dp]
  poly(355,:) = [-7.18E+05_dp, 3.59E+05_dp]
  poly(356,:) = [-7.11E+05_dp, 3.60E+05_dp]
  poly(357,:) = [-7.06E+05_dp, 3.55E+05_dp]
  poly(358,:) = [-7.02E+05_dp, 3.56E+05_dp]
  poly(359,:) = [-6.99E+05_dp, 3.61E+05_dp]
  poly(360,:) = [-6.95E+05_dp, 3.59E+05_dp]
  poly(361,:) = [-6.89E+05_dp, 3.51E+05_dp]
  poly(362,:) = [-6.82E+05_dp, 3.50E+05_dp]
  poly(363,:) = [-6.80E+05_dp, 3.55E+05_dp]
  poly(364,:) = [-6.73E+05_dp, 3.48E+05_dp]
  poly(365,:) = [-6.71E+05_dp, 3.58E+05_dp]
  poly(366,:) = [-6.65E+05_dp, 3.54E+05_dp]
  poly(367,:) = [-6.64E+05_dp, 3.59E+05_dp]
  poly(368,:) = [-6.71E+05_dp, 3.66E+05_dp]
  poly(369,:) = [-6.69E+05_dp, 3.69E+05_dp]
  poly(370,:) = [-6.66E+05_dp, 3.69E+05_dp]
  poly(371,:) = [-6.48E+05_dp, 3.56E+05_dp]
  poly(372,:) = [-6.35E+05_dp, 3.39E+05_dp]
  poly(373,:) = [-6.33E+05_dp, 3.42E+05_dp]
  poly(374,:) = [-6.27E+05_dp, 3.43E+05_dp]
  poly(375,:) = [-6.14E+05_dp, 3.37E+05_dp]
  poly(376,:) = [-6.03E+05_dp, 3.21E+05_dp]
  poly(377,:) = [-5.86E+05_dp, 3.04E+05_dp]
  poly(378,:) = [-5.65E+05_dp, 2.87E+05_dp]
  poly(379,:) = [-5.55E+05_dp, 2.84E+05_dp]
  poly(380,:) = [-5.42E+05_dp, 2.71E+05_dp]
  poly(381,:) = [-5.35E+05_dp, 2.67E+05_dp]
  poly(382,:) = [-5.36E+05_dp, 2.61E+05_dp]
  poly(383,:) = [-5.29E+05_dp, 2.60E+05_dp]
  poly(384,:) = [-5.23E+05_dp, 2.63E+05_dp]
  poly(385,:) = [-5.19E+05_dp, 2.54E+05_dp]
  poly(386,:) = [-5.18E+05_dp, 2.45E+05_dp]
  poly(387,:) = [-5.14E+05_dp, 2.39E+05_dp]
  poly(388,:) = [-5.10E+05_dp, 2.36E+05_dp]
  poly(389,:) = [-5.02E+05_dp, 2.41E+05_dp]
  poly(390,:) = [-4.96E+05_dp, 2.42E+05_dp]
  poly(391,:) = [-4.74E+05_dp, 2.32E+05_dp]
  poly(392,:) = [-4.75E+05_dp, 2.24E+05_dp]
  poly(393,:) = [-4.73E+05_dp, 2.21E+05_dp]
  poly(394,:) = [-4.68E+05_dp, 2.20E+05_dp]
  poly(395,:) = [-4.40E+05_dp, 2.00E+05_dp]
  poly(396,:) = [-4.06E+05_dp, 1.96E+05_dp]
  poly(397,:) = [-3.54E+05_dp, 1.79E+05_dp]
  poly(398,:) = [-3.07E+05_dp, 1.54E+05_dp]
  poly(399,:) = [-2.70E+05_dp, 1.39E+05_dp]
  poly(400,:) = [-1.93E+05_dp, 9.35E+04_dp]
  poly(401,:) = [-1.74E+05_dp, 8.04E+04_dp]
  poly(402,:) = [-1.68E+05_dp, 6.26E+04_dp]
  poly(403,:) = [-1.64E+05_dp, 2.60E+04_dp]
  poly(404,:) = [-1.67E+05_dp, -4.10E+04_dp]
  poly(405,:) = [-1.67E+05_dp, -7.76E+04_dp]
  poly(406,:) = [-1.61E+05_dp, -1.01E+05_dp]
  poly(407,:) = [-1.41E+05_dp, -1.20E+05_dp]
  poly(408,:) = [-1.73E+05_dp, -1.16E+05_dp]
  poly(409,:) = [-2.59E+05_dp, -1.09E+05_dp]
  poly(410,:) = [-2.97E+05_dp, -1.16E+05_dp]
  poly(411,:) = [-3.22E+05_dp, -1.18E+05_dp]
  poly(412,:) = [-3.66E+05_dp, -1.17E+05_dp]
  poly(413,:) = [-3.79E+05_dp, -1.13E+05_dp]
  poly(414,:) = [-4.43E+05_dp, -1.09E+05_dp]
  poly(415,:) = [-4.69E+05_dp, -1.25E+05_dp]
  poly(416,:) = [-4.83E+05_dp, -1.37E+05_dp]
  poly(417,:) = [-5.15E+05_dp, -1.38E+05_dp]
  poly(418,:) = [-5.35E+05_dp, -1.33E+05_dp]
  poly(419,:) = [-5.69E+05_dp, -9.87E+04_dp]
  poly(420,:) = [-6.04E+05_dp, -9.95E+04_dp]
  poly(421,:) = [-6.22E+05_dp, -1.02E+05_dp]
  poly(422,:) = [-6.51E+05_dp, -1.14E+05_dp]
  poly(423,:) = [-6.85E+05_dp, -1.38E+05_dp]
  poly(424,:) = [-7.22E+05_dp, -1.57E+05_dp]
  poly(425,:) = [-7.45E+05_dp, -1.65E+05_dp]
  poly(426,:) = [-7.70E+05_dp, -1.77E+05_dp]
  poly(427,:) = [-7.80E+05_dp, -1.87E+05_dp]
  poly(428,:) = [-7.88E+05_dp, -1.91E+05_dp]
  poly(429,:) = [-7.96E+05_dp, -1.86E+05_dp]
  poly(430,:) = [-8.28E+05_dp, -1.75E+05_dp]
  poly(431,:) = [-8.67E+05_dp, -1.67E+05_dp]
  poly(432,:) = [-9.48E+05_dp, -1.57E+05_dp]
  poly(433,:) = [-9.90E+05_dp, -1.41E+05_dp]
  poly(434,:) = [-1.00E+06_dp, -1.34E+05_dp]
  poly(435,:) = [-1.01E+06_dp, -1.23E+05_dp]
  poly(436,:) = [-1.02E+06_dp, -1.12E+05_dp]
  poly(437,:) = [-1.04E+06_dp, -8.67E+04_dp]
  poly(438,:) = [-1.11E+06_dp, -3.44E+04_dp]
  poly(439,:) = [-1.17E+06_dp, 2.81E+03_dp]
  poly(440,:) = [-1.19E+06_dp, 2.31E+04_dp]
  poly(441,:) = [-1.22E+06_dp, 1.82E+04_dp]
  poly(442,:) = [-1.23E+06_dp, 1.28E+04_dp]
  poly(443,:) = [-1.25E+06_dp, 5.96E+03_dp]
  poly(444,:) = [-1.26E+06_dp, -3.17E+03_dp]
  poly(445,:) = [-1.27E+06_dp, -2.60E+04_dp]
  poly(446,:) = [-1.30E+06_dp, -4.51E+04_dp]
  poly(447,:) = [-1.31E+06_dp, -5.04E+04_dp]
  poly(448,:) = [-1.34E+06_dp, -5.49E+04_dp]
  poly(449,:) = [-1.35E+06_dp, -5.62E+04_dp]
  poly(450,:) = [-1.36E+06_dp, -5.51E+04_dp]
  poly(451,:) = [-1.39E+06_dp, -4.64E+04_dp]
  poly(452,:) = [-1.44E+06_dp, -2.70E+04_dp]
  poly(453,:) = [-1.46E+06_dp, -1.59E+04_dp]
  poly(454,:) = [-1.50E+06_dp, 3.04E+04_dp]
  poly(455,:) = [-1.54E+06_dp, 7.98E+04_dp]
  poly(456,:) = [-1.55E+06_dp, 8.34E+04_dp]
  poly(457,:) = [-1.59E+06_dp, 7.81E+04_dp]
  poly(458,:) = [-1.63E+06_dp, 8.57E+04_dp]
  poly(459,:) = [-1.65E+06_dp, 9.67E+04_dp]
  poly(460,:) = [-1.66E+06_dp, 1.03E+05_dp]
  poly(461,:) = [-1.66E+06_dp, 1.12E+05_dp]
  poly(462,:) = [-1.66E+06_dp, 1.31E+05_dp]
  poly(463,:) = [-1.67E+06_dp, 2.44E+05_dp]
  poly(464,:) = [-1.67E+06_dp, 2.54E+05_dp]
  poly(465,:) = [-1.68E+06_dp, 2.68E+05_dp]
  poly(466,:) = [-1.68E+06_dp, 3.07E+05_dp]
  poly(467,:) = [-1.68E+06_dp, 3.37E+05_dp]
  poly(468,:) = [-1.68E+06_dp, 3.60E+05_dp]
  poly(469,:) = [-1.66E+06_dp, 3.81E+05_dp]
  poly(470,:) = [-1.66E+06_dp, 3.88E+05_dp]
  poly(471,:) = [-1.66E+06_dp, 4.03E+05_dp]
  poly(472,:) = [-1.65E+06_dp, 4.46E+05_dp]
  poly(473,:) = [-1.64E+06_dp, 4.71E+05_dp]
  poly(474,:) = [-1.63E+06_dp, 4.83E+05_dp]
  poly(475,:) = [-1.62E+06_dp, 4.89E+05_dp]
  poly(476,:) = [-1.60E+06_dp, 4.97E+05_dp]
  poly(477,:) = [-1.60E+06_dp, 5.00E+05_dp]
  poly(478,:) = [-1.61E+06_dp, 5.18E+05_dp]
  poly(479,:) = [-1.61E+06_dp, 5.23E+05_dp]
  poly(480,:) = [-1.60E+06_dp, 5.31E+05_dp]
  poly(481,:) = [-1.59E+06_dp, 5.36E+05_dp]
  poly(482,:) = [-1.59E+06_dp, 5.40E+05_dp]
  poly(483,:) = [-1.59E+06_dp, 5.45E+05_dp]
  poly(484,:) = [-1.59E+06_dp, 5.51E+05_dp]
  poly(485,:) = [-1.59E+06_dp, 5.58E+05_dp]
  poly(486,:) = [-1.59E+06_dp, 5.65E+05_dp]
  poly(487,:) = [-1.57E+06_dp, 5.74E+05_dp]
  poly(488,:) = [-1.56E+06_dp, 5.90E+05_dp]
  poly(489,:) = [-1.56E+06_dp, 6.07E+05_dp]
  poly(490,:) = [-1.57E+06_dp, 6.43E+05_dp]
  poly(491,:) = [-1.58E+06_dp, 6.68E+05_dp]
  poly(492,:) = [-1.59E+06_dp, 6.88E+05_dp]
  poly(493,:) = [-1.59E+06_dp, 6.94E+05_dp]
  poly(494,:) = [-1.58E+06_dp, 7.01E+05_dp]
  poly(495,:) = [-1.58E+06_dp, 7.04E+05_dp]
  poly(496,:) = [-1.59E+06_dp, 7.13E+05_dp]
  poly(497,:) = [-1.60E+06_dp, 7.34E+05_dp]
  poly(498,:) = [-1.59E+06_dp, 7.48E+05_dp]
  poly(499,:) = [-1.60E+06_dp, 7.53E+05_dp]

  ! Finalise routine path
  call finalise_routine( routine_name)

end subroutine calc_polygon_Institute_basin



end module mesh_ROI_polygons
