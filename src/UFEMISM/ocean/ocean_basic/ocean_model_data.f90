module ocean_model_data

  use models_basic, only: atype_model
  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_ocean_model_data

  type, abstract, extends(atype_model) :: atype_ocean_model_data
    !< Variables that are common to all ocean models,
    !< which we want other models to be able to access

    ! Main data fields
    real(dp), dimension(:,:), contiguous, pointer :: T => null()   !< [degrees Celsius] Temperature
    real(dp), dimension(:,:), contiguous, pointer :: S => null()   !< [PSU]             Salinity
    type(MPI_WIN) :: wT, wS

    ! Secondary data fields
    real(dp), dimension(:), contiguous, pointer :: T_draft          => null()   !< [degrees Celsius] Temperature at ice base
    real(dp), dimension(:), contiguous, pointer :: T_freezing_point => null()   !< [degrees Celsius] Pressure freezing point of water
    type(MPI_WIN) :: wT_draft, wT_freezing_point

  end type atype_ocean_model_data

end module ocean_model_data