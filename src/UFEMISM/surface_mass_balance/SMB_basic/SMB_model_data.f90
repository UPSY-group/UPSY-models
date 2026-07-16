module SMB_model_data

  use models_basic, only: atype_model
  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_SMB_model_data

  type, abstract, extends(atype_model) :: atype_SMB_model_data
    !< Variables that are common to all SMB models,
    !< which we want other models to be able to access

    real(dp), dimension(:), contiguous, pointer :: SMB => null()   !< [m.i.e. yr^-1] Yearly total SMB
    type(MPI_WIN) :: wSMB

  end type atype_SMB_model_data

end module SMB_model_data