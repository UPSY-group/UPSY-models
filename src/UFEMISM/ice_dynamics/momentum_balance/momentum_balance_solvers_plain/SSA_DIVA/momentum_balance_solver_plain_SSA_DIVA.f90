module momentum_balance_solver_plain_SSA_DIVA

  ! Intermediate class containing all the data and functions
  ! that are shared between the SSA and DIVA solvers

  use precisions, only: dp
  use momentum_balance_solver_plain_basic, only: atype_momentum_balance_solver_plain
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use mesh_disc_apply_operators, only: map_b_a_2D, map_a_b_2D
  use remapping_main, only: map_from_mesh_to_mesh_with_reallocation_2D
  use reallocate_mod, only: reallocate_bounds, reallocate_clean
  use model_configuration, only: C

  implicit none

  private

  public :: atype_momentum_balance_solver_plain_SSA_DIVA

  type, abstract, extends(atype_momentum_balance_solver_plain) :: atype_momentum_balance_solver_plain_SSA_DIVA

      ! Solution
      real(dp), dimension(:  ), allocatable :: u_vav_b                        ! [m yr^-1] 2-D horizontal ice velocity
      real(dp), dimension(:  ), allocatable :: v_vav_b

      ! Intermediate data fields
      real(dp), dimension(:  ), allocatable :: tau_dx_b                       ! Driving stress
      real(dp), dimension(:  ), allocatable :: tau_dy_b
      real(dp), dimension(:  ), allocatable :: du_dx_a                        ! [yr^-1] 2-D horizontal strain rates
      real(dp), dimension(:  ), allocatable :: du_dy_a
      real(dp), dimension(:  ), allocatable :: dv_dx_a
      real(dp), dimension(:  ), allocatable :: dv_dy_a
      real(dp), dimension(:  ), allocatable :: eta_vav_a                      ! Effective viscosity
      real(dp), dimension(:  ), allocatable :: N_a                            ! Product term N = eta * H
      real(dp), dimension(:  ), allocatable :: N_b
      real(dp), dimension(:  ), allocatable :: dN_dx_b                        ! Gradients of N
      real(dp), dimension(:  ), allocatable :: dN_dy_b
      real(dp), dimension(:  ), allocatable :: basal_friction_coefficient_b   ! Basal friction coefficient (tau_b = u * beta_b)
      real(dp), dimension(:  ), allocatable :: u_vav_b_prev                   ! Velocity solution from previous viscosity iteration
      real(dp), dimension(:  ), allocatable :: v_vav_b_prev

      ! Restart file
      character(len=:), allocatable :: restart_filename

    contains

      procedure, public :: allocate_shared_SSA_DIVA_fields
      procedure, public :: deallocate_shared_SSA_DIVA_fields
      procedure, public :: remap_shared_SSA_DIVA_fields

  end type atype_momentum_balance_solver_plain_SSA_DIVA

contains

  subroutine allocate_shared_SSA_DIVA_fields( self)

    ! In/output variables:
    class(atype_momentum_balance_solver_plain_SSA_DIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_shared_SSA_DIVA_fields'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    allocate( self%u_vav_b(                      self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%v_vav_b(                      self%mesh%ti1:self%mesh%ti2), source = 0._dp)

    ! Intermediate data fields
    allocate( self%tau_dx_b(                     self%mesh%ti1:self%mesh%ti2), source = 0._dp)
    allocate( self%tau_dy_b(                     self%mesh%ti1:self%mesh%ti2), source = 0._dp)
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
    allocate( self%u_vav_b_prev(                 self%mesh%nTri             ), source = 0._dp)
    allocate( self%v_vav_b_prev(                 self%mesh%nTri             ), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_shared_SSA_DIVA_fields

  subroutine deallocate_shared_SSA_DIVA_fields( self)

    ! In/output variables:
    class(atype_momentum_balance_solver_plain_SSA_DIVA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_shared_SSA_DIVA_fields'

    ! Add routine to path
    call init_routine( routine_name)

    ! Solution
    deallocate( self%u_vav_b)
    deallocate( self%v_vav_b)

    ! Intermediate data fields
    deallocate( self%tau_dx_b)
    deallocate( self%tau_dy_b)
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
    deallocate( self%u_vav_b_prev)
    deallocate( self%v_vav_b_prev)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_shared_SSA_DIVA_fields

  subroutine remap_shared_SSA_DIVA_fields( self, mesh_old, mesh_new)

    ! In/output variables:
    class(atype_momentum_balance_solver_plain_SSA_DIVA), intent(inout) :: self
    type(type_mesh),                                     intent(in   ) :: mesh_old
    type(type_mesh), target,                             intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'remap_shared_SSA_DIVA_fields'
    real(dp), dimension(:  ), allocatable :: u_vav_a
    real(dp), dimension(:  ), allocatable :: v_vav_a

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap the fields that are re-used during the viscosity iteration
    ! ================================================================

    ! allocate memory for velocities on the a-grid (vertices)
    allocate( u_vav_a( mesh_old%vi1: mesh_old%vi2))
    allocate( v_vav_a( mesh_old%vi1: mesh_old%vi2))

    ! Map data from the triangles of the old mesh to the vertices of the old mesh
    call map_b_a_2D( mesh_old, self%u_vav_b, u_vav_a)
    call map_b_a_2D( mesh_old, self%v_vav_b, v_vav_a)

    ! Remap data from the vertices of the old mesh to the vertices of the new mesh
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, u_vav_a , '2nd_order_conservative')
    call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, C%output_dir, v_vav_a , '2nd_order_conservative')

    ! reallocate memory for the data on the triangles
    call reallocate_bounds( self%u_vav_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%v_vav_b, mesh_new%ti1, mesh_new%ti2)

    ! Map data from the vertices of the new mesh to the triangles of the new mesh
    call map_a_b_2D( mesh_new, u_vav_a, self%u_vav_b)
    call map_a_b_2D( mesh_new, v_vav_a, self%v_vav_b)

    ! reallocate everything else
    ! ==========================

    call reallocate_bounds( self%du_dx_a                     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%du_dy_a                     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%dv_dx_a                     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%dv_dy_a                     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%eta_vav_a                   , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%N_a                         , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%N_b                         , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%dN_dx_b                     , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%dN_dy_b                     , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%basal_friction_coefficient_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%tau_dx_b                    , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%tau_dy_b                    , mesh_new%ti1, mesh_new%ti2)
    call reallocate_clean ( self%u_vav_b_prev                , mesh_new%nTri             )
    call reallocate_clean ( self%v_vav_b_prev                , mesh_new%nTri             )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_shared_SSA_DIVA_fields

end module momentum_balance_solver_plain_SSA_DIVA
