module momentum_balance_solver_plain_data

  use models_basic, only: atype_model
  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_momentum_balance_solver_plain_data

  type, abstract, extends(atype_model) :: atype_momentum_balance_solver_plain_data

    ! Parameters for solving the linearised momentum balance with PETSc
    real(dp) :: PETSc_rtol
    real(dp) :: PETSc_abstol
    integer  :: n_visc_its            ! Number of non-linear viscosity/friction iterations
    integer  :: n_Axb_its             ! Number of iterations needed by PETSc to solve the linearised momentum balance

  end type atype_momentum_balance_solver_plain_data

end module momentum_balance_solver_plain_data