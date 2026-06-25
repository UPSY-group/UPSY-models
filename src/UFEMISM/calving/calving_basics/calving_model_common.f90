module calving_model_common

  use models_basic, only: atype_model
  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_calving_model_common

  type, abstract, extends(atype_model) :: type_calving_model_common

    real(dp), dimension(:), contiguous, pointer :: Hi_calved => null()
    type(MPI_WIN) :: wHi_calved

  end type type_calving_model_common

end module calving_model_common