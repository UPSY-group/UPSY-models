module ocean_model_common

  use models_basic, only: atype_model
  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_ocean_model_common

  type, abstract, extends(atype_model) :: type_ocean_model_common

    ! Main data fields
    real(dp), dimension(:,:), contiguous, pointer :: T => null()   !< [degrees Celsius] Temperature
    real(dp), dimension(:,:), contiguous, pointer :: S => null()   !< [PSU]             Salinity
    type(MPI_WIN) :: wT, wS

    ! Secondary data fields
    real(dp), dimension(:), contiguous, pointer :: T_draft          => null()   !< [degrees Celsius] Temperature at ice base
    real(dp), dimension(:), contiguous, pointer :: T_freezing_point => null()   !< [degrees Celsius] Pressure freezing point of water
    type(MPI_WIN) :: wT_draft, wT_freezing_point

  end type type_ocean_model_common

end module ocean_model_common