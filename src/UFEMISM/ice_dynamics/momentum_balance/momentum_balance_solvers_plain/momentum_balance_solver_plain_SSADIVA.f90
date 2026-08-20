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

      procedure, public :: allocate_shared_SSA_DIVA_variables
      procedure, public :: deallocate_shared_SSA_DIVA_variables

  end type atype_momentum_balance_solver_plain_SSADIVA

contains

  subroutine allocate_shared_SSA_DIVA_variables( self)

    ! In/output variables:
    class(atype_momentum_balance_solver_plain_SSADIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_shared_SSA_DIVA_variables'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    allocate( self%u_vav_b(                      self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%v_vav_b(                      self%mesh%ti1:self%mesh%ti2), source = 0._dp)

    ! Intermediate data fields
    allocate( self%du_dx_a(                      self%mesh%vi1:self%mesh%vi2), source = 0._dp)
    allocate( self%du_dy_a(                      self%mesh%vi1:self%mesh%vi2), source = 0._dp)
    allocate( self%dv_dx_a(                      self%mesh%vi1:self%mesh%vi2), source = 0._dp)
    allocate( self%dv_dy_a(                      self%mesh%vi1:self%mesh%vi2), source = 0._dp)
    allocate( self%eta_vav_a(                    self%mesh%vi1:self%mesh%vi2), source = 0._dp)
    allocate( self%N_a(                          self%mesh%vi1:self%mesh%vi2), source = 0._dp)
    allocate( self%N_b(                          self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%dN_dx_b(                      self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%dN_dy_b(                      self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%basal_friction_coefficient_b( self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%tau_dx_b(                     self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%tau_dy_b(                     self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%u_vav_b_prev(                 self%mesh%nTri             ), source = 0._dp)
    allocate( self%v_vav_b_prev(                 self%mesh%nTri             ), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_shared_SSA_DIVA_variables

  subroutine deallocate_shared_SSA_DIVA_variables( self)

    ! In/output variables:
    class(atype_momentum_balance_solver_plain_SSADIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_shared_SSA_DIVA_variables'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    deallocate( self%u_vav_b)
    deallocate( self%v_vav_b)

    ! Intermediate data fields
    deallocate( self%du_dx_a)
    deallocate( self%du_dy_a)
    deallocate( self%dv_dx_a)
    deallocate( self%dv_dy_a)
    deallocate( self%eta_vav_a)
    deallocate( self%N_a)
    deallocate( self%N_b)
    deallocate( self%dN_dx_b)
    deallocate( self%dN_dy_b)
    deallocate( self%basal_friction_coefficient_b)
    deallocate( self%tau_dx_b)
    deallocate( self%tau_dy_b)
    deallocate( self%u_vav_b_prev)
    deallocate( self%v_vav_b_prev)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_shared_SSA_DIVA_variables

end module momentum_balance_solver_plain_SSADIVA
