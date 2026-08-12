module ice_velocity_model_basic

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use ice_velocity_model_data, only: atype_ice_velocity_model_data
  use mesh_types, only: type_mesh
  use parameters, only: NaN
  use reallocate_mod, only: reallocate_bounds
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use model_configuration, only: C

  implicit none

  private

  public :: atype_ice_velocity_model

  type, abstract, extends(atype_ice_velocity_model_data) :: atype_ice_velocity_model

    contains

      ! Procedures for model memory management and operation
      procedure, public :: allocate   => ice_velocity_model_allocate
      procedure, public :: deallocate => ice_velocity_model_deallocate
      procedure, public :: initialise => ice_velocity_model_initialise
      procedure, public :: run        => ice_velocity_model_run
      procedure, public :: remap      => ice_velocity_model_remap

      ! Deferred procedures that must be overridden by each individual ice velocity model implementation
      procedure(ice_velocity_model_allocate_ifc),   deferred :: allocate_ice_velocity_model
      procedure(ice_velocity_model_deallocate_ifc), deferred :: deallocate_ice_velocity_model
      procedure(ice_velocity_model_initialise_ifc), deferred :: initialise_ice_velocity_model
      procedure(ice_velocity_model_run_ifc),        deferred :: run_ice_velocity_model
      procedure(ice_velocity_model_remap_ifc),      deferred :: remap_ice_velocity_model

      procedure, public                                    :: get_model_name
      procedure(get_ice_velocity_model_name_ifc), deferred :: get_ice_velocity_model_name

  end type atype_ice_velocity_model

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine ice_velocity_model_allocate_ifc( self)
      import atype_ice_velocity_model
      class(atype_ice_velocity_model),  intent(inout) :: self
    end subroutine ice_velocity_model_allocate_ifc

    subroutine ice_velocity_model_deallocate_ifc( self)
      import atype_ice_velocity_model
      class(atype_ice_velocity_model), intent(inout) :: self
    end subroutine ice_velocity_model_deallocate_ifc

    subroutine ice_velocity_model_initialise_ifc( self)
      import atype_ice_velocity_model
      class(atype_ice_velocity_model), intent(inout) :: self
    end subroutine ice_velocity_model_initialise_ifc

    subroutine ice_velocity_model_run_ifc( self)
      import atype_ice_velocity_model
      class(atype_ice_velocity_model), intent(inout) :: self
    end subroutine ice_velocity_model_run_ifc

    subroutine ice_velocity_model_remap_ifc( self, mesh_new)
      import atype_ice_velocity_model, type_mesh
      class(atype_ice_velocity_model), intent(inout) :: self
      type(type_mesh), target,         intent(in   ) :: mesh_new
    end subroutine ice_velocity_model_remap_ifc

    function get_ice_velocity_model_name_ifc( self) result( ice_velocity_model_name)
      import atype_ice_velocity_model
      class(atype_ice_velocity_model), intent(in) :: self
      character(len=:), allocatable :: ice_velocity_model_name
    end function get_ice_velocity_model_name_ifc

  end interface

contains

  subroutine ice_velocity_model_allocate( self, region_name, mesh)

    ! In/output variables:
    class(atype_ice_velocity_model), intent(inout) :: self
    character(len=*),                intent(in   ) :: region_name
    type(type_mesh), target,         intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'ice_velocity_model_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all models
    call self%allocate_model( region_name, mesh)

    ! Allocate all the stuff that is common to all ice velocity models

    ! 3-D
    allocate( self%u_3D  ( mesh%vi1:mesh%vi2, 1:mesh%nz), source = NaN)
    allocate( self%v_3D  ( mesh%vi1:mesh%vi2, 1:mesh%nz), source = NaN)
    allocate( self%u_3D_b( mesh%ti1:mesh%ti2, 1:mesh%nz), source = NaN)
    allocate( self%v_3D_b( mesh%ti1:mesh%ti2, 1:mesh%nz), source = NaN)
    allocate( self%w_3D  ( mesh%vi1:mesh%vi2, 1:mesh%nz), source = NaN)

    ! Vertically averaged
    allocate( self%u_vav     ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%v_vav     ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%u_vav_b   ( mesh%ti1:mesh%ti2), source = NaN)
    allocate( self%v_vav_b   ( mesh%ti1:mesh%ti2), source = NaN)
    allocate( self%uabs_vav  ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%uabs_vav_b( mesh%ti1:mesh%ti2), source = NaN)
    allocate( self%u_vav_perp( mesh%vi1:mesh%vi2, 1:mesh%nC_mem), source = 0._dp)

    ! Surface
    allocate( self%u_surf     ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%v_surf     ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%u_surf_b   ( mesh%ti1:mesh%ti2), source = NaN)
    allocate( self%v_surf_b   ( mesh%ti1:mesh%ti2), source = NaN)
    allocate( self%w_surf     ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%uabs_surf  ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%uabs_surf_b( mesh%ti1:mesh%ti2), source = NaN)

    ! Base
    call self%create_field( self%u_base, self%wu_base, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'u_base', &
      long_name = 'Basal ice velocity in the x-direction', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%v_base, self%wv_base, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'v_base', &
      long_name = 'Basal ice velocity in the y-direction', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%u_base_b, self%wu_base_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'u_base_b', &
      long_name = 'Basal ice velocity in the x-direction on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%v_base_b, self%wv_base_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'v_base_b', &
      long_name = 'Basal ice velocity in the y-direction on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%w_base, self%ww_base, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'w_base', &
      long_name = 'Basal ice velocity in the z-direction', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%uabs_base, self%wuabs_base, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'uabs_base', &
      long_name = 'Basal ice speed', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%uabs_base_b, self%wuabs_base_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'uabs_base_b', &
      long_name = 'Basal ice speed on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    ! Strain rates
    call self%create_field( self%du_dx_3D, self%wdu_dx_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( mesh%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'du_dx_3D', &
      long_name = '3-D xx strain rate', &
      units     = '', &
      remap_method = 'reallocate')

    call self%create_field( self%du_dy_3D, self%wdu_dy_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( mesh%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'du_dy_3D', &
      long_name = '3-D xy strain rate', &
      units     = '', &
      remap_method = 'reallocate')

    call self%create_field( self%du_dz_3D, self%wdu_dz_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( mesh%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'du_dz_3D', &
      long_name = '3-D xz strain rate', &
      units     = '', &
      remap_method = 'reallocate')

    call self%create_field( self%dv_dx_3D, self%wdv_dx_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( mesh%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'dv_dx_3D', &
      long_name = '3-D yx strain rate', &
      units     = '', &
      remap_method = 'reallocate')

    call self%create_field( self%dv_dy_3D, self%wdv_dy_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( mesh%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'dv_dy_3D', &
      long_name = '3-D yy strain rate', &
      units     = '', &
      remap_method = 'reallocate')

    call self%create_field( self%dv_dz_3D, self%wdv_dz_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( mesh%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'dv_dz_3D', &
      long_name = '3-D yz strain rate', &
      units     = '', &
      remap_method = 'reallocate')

    call self%create_field( self%dw_dx_3D, self%wdw_dx_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( mesh%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'dw_dx_3D', &
      long_name = '3-D zx strain rate', &
      units     = '', &
      remap_method = 'reallocate')

    call self%create_field( self%dw_dy_3D, self%wdw_dy_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( mesh%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'dw_dy_3D', &
      long_name = '3-D zy strain rate', &
      units     = '', &
      remap_method = 'reallocate')

    call self%create_field( self%dw_dz_3D, self%wdw_dz_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( mesh%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'dw_dz_3D', &
      long_name = '3-D zz strain rate', &
      units     = '', &
      remap_method = 'reallocate')

    ! Flow regime
    call self%create_field( self%R_shear, self%wR_shear, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'R_shear', &
      long_name = 'Slide/shear ratio', &
      units     = '', &
      remap_method = 'reallocate')

    ! Allocate stuff that is specific to each individual ice velocity model implementation
    call self%allocate_ice_velocity_model()

    ! Remove routine from call stack
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 10)

  end subroutine ice_velocity_model_allocate

  subroutine ice_velocity_model_deallocate( self)

    ! In/output variables:
    class(atype_ice_velocity_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'ice_velocity_model_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate stuff that is common to all models
    call self%deallocate_model()

    ! Deallocate stuff that is common to all ice velocity models

    ! 3-D
    deallocate( self%u_3D  )
    deallocate( self%v_3D  )
    deallocate( self%u_3D_b)
    deallocate( self%v_3D_b)
    deallocate( self%w_3D  )

    ! Vertically averaged
    deallocate( self%u_vav     )
    deallocate( self%v_vav     )
    deallocate( self%u_vav_b   )
    deallocate( self%v_vav_b   )
    deallocate( self%uabs_vav  )
    deallocate( self%uabs_vav_b)
    deallocate( self%u_vav_perp)

    ! Surface
    deallocate( self%u_surf     )
    deallocate( self%v_surf     )
    deallocate( self%u_surf_b   )
    deallocate( self%v_surf_b   )
    deallocate( self%w_surf     )
    deallocate( self%uabs_surf  )
    deallocate( self%uabs_surf_b)

    ! Base
    nullify( self%u_base     )
    nullify( self%v_base     )
    nullify( self%u_base_b   )
    nullify( self%v_base_b   )
    nullify( self%w_base     )
    nullify( self%uabs_base  )
    nullify( self%uabs_base_b)

    ! Strain rates
    nullify( self%du_dx_3D)
    nullify( self%du_dy_3D)
    nullify( self%du_dz_3D)
    nullify( self%dv_dx_3D)
    nullify( self%dv_dy_3D)
    nullify( self%dv_dz_3D)
    nullify( self%dw_dx_3D)
    nullify( self%dw_dy_3D)
    nullify( self%dw_dz_3D)

    ! Flow regime
    nullify( self%R_shear)

    ! Deallocate stuff that is specific to each individual ice velocity model implementation
    call self%deallocate_ice_velocity_model()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ice_velocity_model_deallocate

  subroutine ice_velocity_model_initialise( self)

    ! In/output variables:
    class(atype_ice_velocity_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'ice_velocity_model_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise stuff that is common to all models
    call self%initialise_model()

    ! Initialise stuff that is common to all ice_velocity models

    ! Initialise stuff that is specific to each individual ice_velocity model implementation
    call self%initialise_ice_velocity_model()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ice_velocity_model_initialise

  subroutine ice_velocity_model_run( self)

    ! In/output variables:
    class(atype_ice_velocity_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'ice_velocity_model_run'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run stuff that is common to all models
    call self%run_model()

    ! Run stuff that is common to all ice_velocity models

    ! Run stuff that is specific to each individual ice_velocity model implementation
    call self%run_ice_velocity_model()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ice_velocity_model_run

  subroutine ice_velocity_model_remap( self, mesh_new)

    ! In/output variables:
    class(atype_ice_velocity_model), intent(inout) :: self
    type(type_mesh),                 intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'ice_velocity_model_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap stuff that is common to all models
    call self%remap_model( mesh_new)

    ! Remap stuff that is common to all ice_velocity models

    ! 3-D
    call reallocate_bounds( self%u_3D  , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    call reallocate_bounds( self%v_3D  , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)
    call reallocate_bounds( self%u_3D_b, mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( self%v_3D_b, mesh_new%ti1, mesh_new%ti2, mesh_new%nz)
    call reallocate_bounds( self%w_3D  , mesh_new%vi1, mesh_new%vi2, mesh_new%nz)

    ! Vertically averaged
    call reallocate_bounds( self%u_vav     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%v_vav     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%u_vav_b   , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%v_vav_b   , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%uabs_vav  , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%uabs_vav_b, mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%u_vav_perp, mesh_new%vi1, mesh_new%vi2, mesh_new%nC_mem)

    ! Surface
    call reallocate_bounds( self%u_surf     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%v_surf     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%u_surf_b   , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%v_surf_b   , mesh_new%ti1, mesh_new%ti2)
    call reallocate_bounds( self%w_surf     , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%uabs_surf  , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( self%uabs_surf_b, mesh_new%ti1, mesh_new%ti2)

    ! Base
    call self%remap_field( mesh_new, 'u_base'     , self%u_base     )
    call self%remap_field( mesh_new, 'v_base'     , self%v_base     )
    call self%remap_field( mesh_new, 'u_base_b'   , self%u_base_b   )
    call self%remap_field( mesh_new, 'v_base_b'   , self%v_base_b   )
    call self%remap_field( mesh_new, 'w_base'     , self%w_base     )
    call self%remap_field( mesh_new, 'uabs_base'  , self%uabs_base  )
    call self%remap_field( mesh_new, 'uabs_base_b', self%uabs_base_b)

    ! Strain rates
    call self%remap_field( mesh_new, 'du_dx_3D', self%du_dx_3D)
    call self%remap_field( mesh_new, 'du_dy_3D', self%du_dy_3D)
    call self%remap_field( mesh_new, 'du_dz_3D', self%du_dz_3D)
    call self%remap_field( mesh_new, 'dv_dx_3D', self%dv_dx_3D)
    call self%remap_field( mesh_new, 'dv_dy_3D', self%dv_dy_3D)
    call self%remap_field( mesh_new, 'dv_dz_3D', self%dv_dz_3D)
    call self%remap_field( mesh_new, 'dw_dx_3D', self%dw_dx_3D)
    call self%remap_field( mesh_new, 'dw_dy_3D', self%dw_dy_3D)
    call self%remap_field( mesh_new, 'dw_dz_3D', self%dw_dz_3D)

    ! Flow regime
    call self%remap_field( mesh_new, 'R_shear', self%R_shear)

    ! Remap stuff that is specific to each individual ice_velocity model implementation
    call self%remap_ice_velocity_model( mesh_new)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ice_velocity_model_remap

  function get_model_name( self) result( model_name)
    class(atype_ice_velocity_model), intent(in) :: self
    character(len=:), allocatable      :: model_name
    model_name = 'ice_velocity_' // self%get_ice_velocity_model_name()
  end function get_model_name

end module ice_velocity_model_basic