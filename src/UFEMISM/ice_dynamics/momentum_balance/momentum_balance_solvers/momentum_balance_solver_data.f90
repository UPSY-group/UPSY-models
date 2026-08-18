module momentum_balance_solver_data

  ! A wrapper class that can contain multiple instances of different
  ! momentum balance solvers. This is needed for e.g. the SIA/SSA,
  ! and for the hybrid DIVA/BPA, both of which need to run code from
  ! two separate solvers.

  use models_basic, only: atype_model
  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_momentum_balance_solver_data

  type, abstract, extends(atype_model) :: atype_momentum_balance_solver_data

  end type atype_momentum_balance_solver_data

end module momentum_balance_solver_data