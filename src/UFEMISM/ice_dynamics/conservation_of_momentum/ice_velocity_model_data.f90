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

  end type atype_ice_velocity_model_data

end module ice_velocity_model_data