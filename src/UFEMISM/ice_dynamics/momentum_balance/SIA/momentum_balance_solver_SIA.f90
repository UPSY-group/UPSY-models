module momentum_balance_solver_SIA

  !< Routines for solving the Shallow Ice Approximation (SIA) to the momentum balance

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_data, only: atype_ice_velocity_model_data
  use parameters, only: grav, ice_density, NaN
  use reallocate_mod, only: reallocate_bounds
  use constitutive_equation, only: calc_ice_rheology_Glen
  use mesh_disc_apply_operators, only: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, map_a_b_3D, ddx_a_a_2D, ddy_a_a_2D
  use mesh_zeta, only: integrate_from_zeta_is_one_to_zeta_is_zetap
  use momentum_balance_solver_basic, only: atype_momentum_balance_solver
  use mpi_f08, only: MPI_WIN
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use bed_roughness_model_types, only: type_bed_roughness_model

  implicit none

  private

  public :: type_momentum_balance_solver_SIA

  type, extends(atype_momentum_balance_solver) :: type_momentum_balance_solver_SIA

      real(dp), dimension(:,:), contiguous, pointer :: u_3D_b    => null()   ! [m yr^-1] 3-D ice velocity in the x-direction
      real(dp), dimension(:,:), contiguous, pointer :: v_3D_b    => null()   ! [m yr^-1] 3-D ice velocity in the y-direction
      real(dp), dimension(:,:), contiguous, pointer :: D_3D_b    => null()   ! [m yr^-1] Diffusivity
      real(dp), dimension(:,:), contiguous, pointer :: du_dz_3D  => null()   ! [yr^-1]   3-D xz strain rate
      real(dp), dimension(:,:), contiguous, pointer :: dv_dz_3D  => null()   ! [yr^-1]   3-D yz strain rate
      type(MPI_WIN) :: wu_3D_b, wv_3D_b, wD_3D_b, wdu_dz_3D, wdv_dz_3D

    contains

      ! Procedures for model memory management and operation
      procedure, public :: allocate_momentum_balance_solver   => momentum_balance_solver_SIA_allocate
      procedure, public :: deallocate_momentum_balance_solver => momentum_balance_solver_SIA_deallocate
      procedure, public :: initialise_momentum_balance_solver => momentum_balance_solver_SIA_initialise
      procedure, public :: run_momentum_balance_solver        => momentum_balance_solver_SIA_run
      procedure, public :: remap_momentum_balance_solver      => momentum_balance_solver_SIA_remap

      procedure, public :: get_momentum_balance_solver_name

  end type type_momentum_balance_solver_SIA

contains

  subroutine momentum_balance_solver_SIA_allocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SIA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SIA_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is specific to the SIA momentum balance solver
    call self%create_field( self%u_3D_b, self%wu_3D_b, &
      self%mesh, Arakawa_grid%b(), third_dimension%ice_zeta( C%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'u_3D_b', &
      long_name = '3-D ice velocity in the x-direction', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%v_3D_b, self%wv_3D_b, &
      self%mesh, Arakawa_grid%b(), third_dimension%ice_zeta( C%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'v_3D_b', &
      long_name = '3-D ice velocity in the y-direction', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%D_3D_b, self%wD_3D_b, &
      self%mesh, Arakawa_grid%b(), third_dimension%ice_zeta( C%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'D_3D_b', &
      long_name = '3-D SIA ice diffusivity on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%du_dz_3D, self%wdu_dz_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( C%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'du_dz_3D', &
      long_name = '3-D xz strain rate', &
      units     = 'yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%dv_dz_3D, self%wdv_dz_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( C%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'dv_dz_3D', &
      long_name = '3-D yz strain rate', &
      units     = 'yr^-1', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIA_allocate

  subroutine momentum_balance_solver_SIA_deallocate( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SIA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SIA_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to momentum balance solver SIA
    nullify( self%u_3D_b)
    nullify( self%v_3D_b)
    nullify( self%D_3D_b)
    nullify( self%du_dz_3D)
    nullify( self%dv_dz_3D)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIA_deallocate

  subroutine momentum_balance_solver_SIA_initialise( self)

    ! In/output variables:
    class(type_momentum_balance_solver_SIA), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SIA_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise all the stuff that is specific to momentum balance solver SIA


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIA_initialise

  subroutine momentum_balance_solver_SIA_run( self, ice, geom, bed_roughness, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b, BC_prescr_mask_bk, BC_prescr_u_bk, BC_prescr_v_bk)

    ! In/output variables:
    class(type_momentum_balance_solver_SIA), intent(inout) :: self
    class(atype_ice_model_data),                   intent(inout) :: ice
    class(atype_ice_geometry_model_data),          intent(in   ) :: geom
    type(type_bed_roughness_model),                intent(in   ) :: bed_roughness
    integer,  dimension(:  ), optional,            intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:  ), optional,            intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(:  ), optional,            intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction
    integer,  dimension(:,:), optional,            intent(in   ) :: BC_prescr_mask_bk     ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:,:), optional,            intent(in   ) :: BC_prescr_u_bk        ! Prescribed velocities in the x-direction
    real(dp), dimension(:,:), optional,            intent(in   ) :: BC_prescr_v_bk        ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'run_momentum_balance_solver_SIA'
    real(dp), dimension(:  ), allocatable :: Hi_b, Hs_b, dHs_dx, dHs_dy, dHs_dx_b, dHs_dy_b
    real(dp), dimension(:,:), allocatable :: A_flow_b
    integer                               :: vi,ti,k
    real(dp)                              :: abs_grad_Hs
    real(dp), dimension(self%mesh%nz)     :: z, int_A_hminzetan

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is specific to momentum balance solver SIA

    ! Safety
    if (.not. C%choice_flow_law == 'Glen') then
      call crash('the analytical solution to the SIA is only valid when using Glens flow law!')
    end if

    ! Allocate memory
    allocate( Hi_b(     self%mesh%ti1:self%mesh%ti2                ))
    allocate( Hs_b(     self%mesh%ti1:self%mesh%ti2                ))
    allocate( dHs_dx(   self%mesh%vi1:self%mesh%vi2                ))
    allocate( dHs_dy(   self%mesh%vi1:self%mesh%vi2                ))
    allocate( dHs_dx_b( self%mesh%ti1:self%mesh%ti2                ))
    allocate( dHs_dy_b( self%mesh%ti1:self%mesh%ti2                ))
    allocate( A_flow_b( self%mesh%ti1:self%mesh%ti2, 1:self%mesh%nz))

    ! Calculate flow factors
    call calc_ice_rheology_Glen( self%mesh, ice, geom)

    ! Calculate ice thickness, surface elevation, surface slopes, and ice flow factor on the b-grid
    call map_a_b_2D( self%mesh, geom%Hi   , Hi_b    )
    call map_a_b_2D( self%mesh, geom%Hs   , Hs_b    )
    call ddx_a_a_2D( self%mesh, geom%Hs   , dHs_dx  )
    call ddy_a_a_2D( self%mesh, geom%Hs   , dHs_dy  )
    call ddx_a_b_2D( self%mesh, geom%Hs   , dHs_dx_b)
    call ddy_a_b_2D( self%mesh, geom%Hs   , dHs_dy_b)
    call map_a_b_3D( self%mesh, ice%A_flow, A_flow_b)

    ! Calculate velocities and strain rates according to the analytical solution of the SIA:
    ! (see also Bueler and Brown, 2009, Eqs. 12-13)
    !
    !   D( z) = -2 (rho g)^n (abs(grad H))^(n-1) int_b_z( A(T*) (h - zeta)^n ) dzeta
    !   u( z) = dh/dx D( z)
    !   v( z) = dh/dy D( z)
    !
    !   du/dz( z) = -2 (rho g)^n (abs(grad H))^(n-1) A(T*) (h - z)^n dh/dx
    !   dv/dz( z) = -2 (rho g)^n (abs(grad H))^(n-1) A(T*) (h - z)^n dh/dy

    ! Calculate velocities
    do ti = self%mesh%ti1, self%mesh%ti2

      ! Calculate the integral from b to z of (A_flow * (h - zeta)^n) dzeta
      z = Hs_b( ti) - self%mesh%zeta * Hi_b( ti)
      int_A_hminzetan = integrate_from_zeta_is_one_to_zeta_is_zetap( z, A_flow_b( ti,:) * (Hs_b( ti) - z)**C%Glens_flow_law_exponent)

      ! Calculate the diffusivity term
      abs_grad_Hs = sqrt( dHs_dx_b( ti)**2 + dHs_dy_b( ti)**2)
      self%D_3D_b( ti,:) = -2._dp * (ice_density * grav)**C%Glens_flow_law_exponent * &
        abs_grad_Hs**(C%Glens_flow_law_exponent - 1._dp) * int_A_hminzetan

      ! Safety
      self%D_3D_b( ti,:) = max( -C%SIA_maximum_diffusivity, self%D_3D_b( ti,:))

      ! Calculate the velocities
      self%u_3D_b( ti,:) = self%D_3D_b( ti,:) * dHs_dx_b( ti)
      self%v_3D_b( ti,:) = self%D_3D_b( ti,:) * dHs_dy_b( ti)

    end do

    ! Calculate vertical shear strain rates (needed later to calculate strain heating in thermodynamics)
    do vi = self%mesh%vi1, self%mesh%vi2

      abs_grad_Hs = sqrt( dHs_dx( vi)**2 + dHs_dy( vi)**2)
      z = geom%Hs( vi) - self%mesh%zeta * geom%Hi( vi)

      do k = 1, self%mesh%nz
        self%du_dz_3D( vi,k) = -2._dp * (ice_density * grav)**C%Glens_flow_law_exponent * &
          abs_grad_Hs**(C%Glens_flow_law_exponent - 1._dp) * &
          ice%A_flow( vi,k) * (geom%Hs( vi) - z( k))**C%Glens_flow_law_exponent * dHs_dx( vi)
        self%dv_dz_3D( vi,k) = -2._dp * (ice_density * grav)**C%Glens_flow_law_exponent * &
          abs_grad_Hs**(C%Glens_flow_law_exponent - 1._dp) * &
          ice%A_flow( vi,k) * (geom%Hs( vi) - z( k))**C%Glens_flow_law_exponent * dHs_dy( vi)
      end do

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIA_run

  subroutine momentum_balance_solver_SIA_remap( self, mesh_old, mesh_new)

    ! In/output variables:
    class(type_momentum_balance_solver_SIA), intent(inout) :: self
    type(type_mesh),                               intent(in   ) :: mesh_old
    type(type_mesh), target,                       intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'momentum_balance_solver_SIA_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap all the stuff that is specific to momentum balance solver SIA
    call self%remap_field( mesh_new, 'u_3D_b'  , self%u_3D_b  )
    call self%remap_field( mesh_new, 'v_3D_b'  , self%v_3D_b  )
    call self%remap_field( mesh_new, 'D_3D_b'  , self%D_3D_b  )
    call self%remap_field( mesh_new, 'du_dz_3D', self%du_dz_3D)
    call self%remap_field( mesh_new, 'dv_dz_3D', self%dv_dz_3D)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine momentum_balance_solver_SIA_remap

  function get_momentum_balance_solver_name( self) result( model_name)
    class(type_momentum_balance_solver_SIA), intent(in) :: self
    character(len=:), allocatable                  :: model_name
    model_name = 'SIA'
  end function get_momentum_balance_solver_name

end module momentum_balance_solver_SIA
