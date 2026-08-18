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

    integer  :: n_visc_its            ! Number of non-linear viscosity/friction iterations
    integer  :: n_Axb_its             ! Number of iterations needed by PETSc to solve the linearised momentum balance

  end type atype_momentum_balance_solver_data

end module momentum_balance_solver_data