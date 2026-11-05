module climate_model_data

  use models_basic, only: atype_model
  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_climate_model_data

  type, abstract, extends(atype_model) :: atype_climate_model_data
    !< Variables that are common to all climate models,
    !< which we want other models to be able to access

    real(dp), dimension(:,:), contiguous, pointer :: T2m    => null()
    real(dp), dimension(:,:), contiguous, pointer :: Precip => null()
    real(dp), dimension(:,:), contiguous, pointer :: Q_TOA  => null()
    type(MPI_WIN) :: wT2m, wPrecip, wQ_TOA

  end type atype_climate_model_data

end module climate_model_data
