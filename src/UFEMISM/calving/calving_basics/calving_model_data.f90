module calving_model_data

  use models_basic, only: atype_model
  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_calving_model_data

  type, abstract, extends(atype_model) :: atype_calving_model_data
    !< Variables that are common to all calving models,
    !< which we want other models to be able to access

    real(dp), dimension(:), contiguous, pointer :: Hi_calved => null()
    type(MPI_WIN) :: wHi_calved

  end type atype_calving_model_data

end module calving_model_data