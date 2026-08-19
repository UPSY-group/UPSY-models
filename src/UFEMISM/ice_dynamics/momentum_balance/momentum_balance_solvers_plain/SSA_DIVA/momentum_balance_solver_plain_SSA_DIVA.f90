module momentum_balance_solver_plain_SSA_DIVA

  use momentum_balance_solver_plain_basic, only: atype_momentum_balance_solver_plain

  implicit none

  private

  public :: atype_momentum_balance_solver_plain_SSA_DIVA

  type, abstract, extends(atype_momentum_balance_solver_plain) :: atype_momentum_balance_solver_plain_SSA_DIVA

    contains

  end type atype_momentum_balance_solver_plain_SSA_DIVA

contains

end module momentum_balance_solver_plain_SSA_DIVA
