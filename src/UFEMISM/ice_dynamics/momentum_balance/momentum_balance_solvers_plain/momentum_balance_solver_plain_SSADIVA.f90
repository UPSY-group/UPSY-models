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

      ! Solution
      real(dp), dimension(:), allocatable :: u_vav_b                     ! [m yr^-1] 2-D horizontal ice velocity
      real(dp), dimension(:), allocatable :: v_vav_b

      ! Intermediate data fields
      real(dp), dimension(:), allocatable :: du_dx_a                     ! [yr^-1] 2-D horizontal strain rates
      real(dp), dimension(:), allocatable :: du_dy_a
      real(dp), dimension(:), allocatable :: dv_dx_a
      real(dp), dimension(:), allocatable :: dv_dy_a
      real(dp), dimension(:), allocatable :: eta_vav_a                   ! Effective viscosity
      real(dp), dimension(:), allocatable :: N_a                         ! Product term N = eta * H
      real(dp), dimension(:), allocatable :: N_b
      real(dp), dimension(:), allocatable :: dN_dx_b                     ! Gradients of N
      real(dp), dimension(:), allocatable :: dN_dy_b
      real(dp), dimension(:), allocatable :: basal_friction_coefficient_b! Basal friction coefficient (tau_b = u * beta_b)
      real(dp), dimension(:), allocatable :: tau_dx_b                    ! Driving stress
      real(dp), dimension(:), allocatable :: tau_dy_b
      real(dp), dimension(:), allocatable :: u_vav_b_prev                ! Velocity solution from previous viscosity iteration
      real(dp), dimension(:), allocatable :: v_vav_b_prev

      ! Restart file
      character(len=256)                  :: restart_filename

    contains

  end type atype_momentum_balance_solver_plain_SSADIVA

contains

end module momentum_balance_solver_plain_SSADIVA
