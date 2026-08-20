module momentum_balance_solver_plain_SSADIVA

  ! Intermediate class containing data and functions that are
  ! shared between the SSA and DIVA momentum balance solvers

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use momentum_balance_solver_plain_basic, only: atype_momentum_balance_solver_plain

  implicit none

  private

  public :: atype_momentum_balance_solver_plain_SSADIVA

  type, abstract, extends(atype_momentum_balance_solver_plain) :: atype_momentum_balance_solver_plain_SSADIVA

    contains

  end type atype_momentum_balance_solver_plain_SSADIVA

contains

end module momentum_balance_solver_plain_SSADIVA
