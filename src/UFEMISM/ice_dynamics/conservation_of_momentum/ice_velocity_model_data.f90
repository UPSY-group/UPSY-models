module ice_velocity_model_data

  use models_basic, only: atype_model
  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_ice_velocity_model_data

  type, abstract, extends(atype_model) :: atype_ice_velocity_model_data

    ! 3-D
    real(dp), dimension(:,:), allocatable :: u_3D                        ! [m yr^-1]   3-D ice velocity in the x-direction
    real(dp), dimension(:,:), allocatable :: v_3D                        ! [m yr^-1]   3-D ice velocity in the y-direction
    real(dp), dimension(:,:), allocatable :: u_3D_b                      ! [m yr^-1]   3-D ice velocity in the x-direction on the b-grid
    real(dp), dimension(:,:), allocatable :: v_3D_b                      ! [m yr^-1]   3-D ice velocity in the y-direction on the b-grid
    real(dp), dimension(:,:), allocatable :: w_3D                        ! [m yr^-1]   3-D ice velocity in the z-direction

    ! Vertically averaged
    real(dp), dimension(:  ), allocatable :: u_vav                       ! [m yr^-1]   Vertically averaged ice velocity in the x-direction
    real(dp), dimension(:  ), allocatable :: v_vav                       ! [m yr^-1]   Vertically averaged ice velocity in the y-direction
    real(dp), dimension(:  ), allocatable :: u_vav_b                     ! [m yr^-1]   Vertically averaged ice velocity in the x-direction on the b-grid
    real(dp), dimension(:  ), allocatable :: v_vav_b                     ! [m yr^-1]   Vertically averaged ice velocity in the y-direction on the b-grid
    real(dp), dimension(:  ), allocatable :: uabs_vav                    ! [m yr^-1]   Vertically averaged horizontal ice speed
    real(dp), dimension(:  ), allocatable :: uabs_vav_b                  ! [m yr^-1]   Vertically averaged horizontal ice speed on the b-grid

    ! Surface
    real(dp), dimension(:  ), allocatable :: u_surf                      ! [m yr^-1]   Surface ice velocity in the x-direction
    real(dp), dimension(:  ), allocatable :: v_surf                      ! [m yr^-1]   Surface ice velocity in the y-direction
    real(dp), dimension(:  ), allocatable :: u_surf_b                    ! [m yr^-1]   Surface ice velocity in the x-direction on the b-grid
    real(dp), dimension(:  ), allocatable :: v_surf_b                    ! [m yr^-1]   Surface ice velocity in the x-direction on the b-grid
    real(dp), dimension(:  ), allocatable :: w_surf                      ! [m yr^-1]   Surface ice velocity in the z-direction
    real(dp), dimension(:  ), allocatable :: uabs_surf                   ! [m yr^-1]   Surface ice speed
    real(dp), dimension(:  ), allocatable :: uabs_surf_b                 ! [m yr^-1]   Surface ice speed on the b-grid

    ! Base
    real(dp), dimension(:  ), allocatable :: u_base                      ! [m yr^-1]   Basal ice velocity in the x-direction
    real(dp), dimension(:  ), allocatable :: v_base                      ! [m yr^-1]   Basal ice velocity in the y-direction
    real(dp), dimension(:  ), allocatable :: u_base_b                    ! [m yr^-1]   Basal ice velocity in the x-direction on the b-grid
    real(dp), dimension(:  ), allocatable :: v_base_b                    ! [m yr^-1]   Basal ice velocity in the x-direction on the b-grid
    real(dp), dimension(:  ), allocatable :: w_base                      ! [m yr^-1]   Basal ice velocity in the z-direction
    real(dp), dimension(:  ), allocatable :: uabs_base                   ! [m yr^-1]   Basal ice speed
    real(dp), dimension(:  ), allocatable :: uabs_base_b                 ! [m yr^-1]   Basal ice speed on the b-grid

    ! Strain rates
    real(dp), dimension(:,:), allocatable :: du_dx_3D                    ! [yr^-1]     3-D xx strain rate
    real(dp), dimension(:,:), allocatable :: du_dy_3D                    ! [yr^-1]     3-D xy strain rate
    real(dp), dimension(:,:), allocatable :: du_dz_3D                    ! [yr^-1]     3-D xz strain rate
    real(dp), dimension(:,:), allocatable :: dv_dx_3D                    ! [yr^-1]     3-D yx strain rate
    real(dp), dimension(:,:), allocatable :: dv_dy_3D                    ! [yr^-1]     3-D yy strain rate
    real(dp), dimension(:,:), allocatable :: dv_dz_3D                    ! [yr^-1]     3-D yz strain rate
    real(dp), dimension(:,:), allocatable :: dw_dx_3D                    ! [yr^-1]     3-D zx strain rate
    real(dp), dimension(:,:), allocatable :: dw_dy_3D                    ! [yr^-1]     3-D zy strain rate
    real(dp), dimension(:,:), allocatable :: dw_dz_3D                    ! [yr^-1]     3-D zz strain rate

  end type atype_ice_velocity_model_data

end module ice_velocity_model_data