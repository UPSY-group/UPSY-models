module SMB_model_common

  use models_basic, only: atype_model
  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_SMB_model_common

  type, abstract, extends(atype_model) :: type_SMB_model_common

    real(dp), dimension(:), contiguous, pointer :: SMB => null()
    type(MPI_WIN) :: wSMB

  end type type_SMB_model_common

end module SMB_model_common