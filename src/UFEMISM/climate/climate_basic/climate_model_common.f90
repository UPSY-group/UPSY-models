module climate_model_common

  use models_basic, only: atype_model
  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_climate_model_common

  type, abstract, extends(atype_model) :: type_climate_model_common

    real(dp), dimension(:,:), contiguous, pointer :: T2m    => null()
    real(dp), dimension(:,:), contiguous, pointer :: Precip => null()
    type(MPI_WIN) :: wT2m, wPrecip

  end type type_climate_model_common

end module climate_model_common