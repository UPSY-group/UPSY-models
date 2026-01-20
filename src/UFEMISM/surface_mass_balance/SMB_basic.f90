module SMB_basic
  !< The basic surface mass balance model

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mpi_f08, only: MPI_WIN
  use UPSY_main, only: atype_model

  implicit none

  private

  public :: atype_SMB_model

  type, abstract, extends(atype_model) :: atype_SMB_model
    !< The basic surface mass balance model

      ! Main data fields
      real(dp), dimension(:), contiguous, pointer :: SMB   !< [m] Net annual  SMB
      type(MPI_WIN) :: wSMB

      ! Timestepping
      real(dp) :: t_next

    contains

      private

  end type atype_SMB_model

contains

end module SMB_basic