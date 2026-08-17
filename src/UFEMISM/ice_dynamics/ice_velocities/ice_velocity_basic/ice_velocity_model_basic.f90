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
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data

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

    subroutine ice_velocity_model_run_ifc( self, ice, geom)
      import atype_ice_velocity_model, atype_ice_model_data, atype_ice_geometry_model_data
      class(atype_ice_velocity_model),      intent(inout) :: self
      class(atype_ice_model_data),          intent(inout) :: ice
      class(atype_ice_geometry_model_data), intent(in   ) :: geom
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
    call self%create_field( self%u_3D, self%wu_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( C%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'u_3D', &
      long_name = '3-D ice velocity in the x-direction', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%v_3D, self%wv_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( C%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'v_3D', &
      long_name = '3-D ice velocity in the y-direction', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%u_3D_b, self%wu_3D_b, &
      self%mesh, Arakawa_grid%b(), third_dimension%ice_zeta( C%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'u_3D_b', &
      long_name = '3-D ice velocity in the x-direction on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%v_3D_b, self%wv_3D_b, &
      self%mesh, Arakawa_grid%b(), third_dimension%ice_zeta( C%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'v_3D_b', &
      long_name = '3-D ice velocity in the y-direction on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%w_3D, self%ww_3D, &
      self%mesh, Arakawa_grid%a(), third_dimension%ice_zeta( C%nz, C%choice_zeta_grid, C%zeta_irregular_log_R), &
      name      = 'w_3D', &
      long_name = '3-D ice velocity in the z-direction', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    ! Vertically averaged
    call self%create_field( self%u_vav, self%wu_vav, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'u_vav', &
      long_name = 'Vertically averaged ice velocity in the x-direction', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%v_vav, self%wv_vav, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'v_vav', &
      long_name = 'Vertically averaged ice velocity in the y-direction', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%u_vav_b, self%wu_vav_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'u_vav_b', &
      long_name = 'Vertically averaged ice velocity in the x-direction on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%v_vav_b, self%wv_vav_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'v_vav_b', &
      long_name = 'Vertically averaged ice velocity in the y-direction on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%uabs_vav, self%wuabs_vav, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'uabs_vav', &
      long_name = 'Vertically averaged ice speed', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%uabs_vav_b, self%wuabs_vav_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'uabs_vav_b', &
      long_name = 'Vertically averaged ice speed on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%u_vav_perp, self%wu_vav_perp, &
      self%mesh, Arakawa_grid%b(), third_dimension%vertex_connectivity( mesh%nC_mem), &
      name      = 'u_vav_perp', &
      long_name = 'Vertically averaged ice velocity perpendicular to Voronoi cell boundaries', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    ! DENK DROM
    self%u_vav     ( self%mesh%vi1:self%mesh%vi2) = 0._dp
    self%v_vav     ( self%mesh%vi1:self%mesh%vi2) = 0._dp
    self%u_vav_b   ( self%mesh%ti1:self%mesh%ti2) = 0._dp
    self%v_vav_b   ( self%mesh%ti1:self%mesh%ti2) = 0._dp
    self%uabs_vav  ( self%mesh%vi1:self%mesh%vi2) = 0._dp
    self%uabs_vav_b( self%mesh%ti1:self%mesh%ti2) = 0._dp
    self%u_vav_perp( self%mesh%ti1:self%mesh%ti2,:) = 0._dp

    ! Surface
    call self%create_field( self%u_surf, self%wu_surf, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'u_surf', &
      long_name = 'Surface ice velocity in the x-direction', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%v_surf, self%wv_surf, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'v_surf', &
      long_name = 'Surface ice velocity in the y-direction', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%u_surf_b, self%wu_surf_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'u_surf_b', &
      long_name = 'Surface ice velocity in the x-direction on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%v_surf_b, self%wv_surf_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'v_surf_b', &
      long_name = 'Surface ice velocity in the y-direction on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%w_surf, self%ww_surf, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'w_surf', &
      long_name = 'Surface ice velocity in the z-direction', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%uabs_surf, self%wuabs_surf, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'uabs_surf', &
      long_name = 'Surface ice speed', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%uabs_surf_b, self%wuabs_surf_b, &
      self%mesh, Arakawa_grid%b(), &
      name      = 'uabs_surf_b', &
      long_name = 'Surface ice speed on the triangles', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

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
    nullify( self%u_3D  )
    nullify( self%v_3D  )
    nullify( self%u_3D_b)
    nullify( self%v_3D_b)
    nullify( self%w_3D  )

    ! Vertically averaged
    nullify( self%u_vav     )
    nullify( self%v_vav     )
    nullify( self%u_vav_b   )
    nullify( self%v_vav_b   )
    nullify( self%uabs_vav  )
    nullify( self%uabs_vav_b)
    nullify( self%u_vav_perp)

    ! Surface
    nullify( self%u_surf     )
    nullify( self%v_surf     )
    nullify( self%u_surf_b   )
    nullify( self%v_surf_b   )
    nullify( self%w_surf     )
    nullify( self%uabs_surf  )
    nullify( self%uabs_surf_b)

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

  subroutine ice_velocity_model_run( self, ice, geom)

    ! In/output variables:
    class(atype_ice_velocity_model),      intent(inout) :: self
    class(atype_ice_model_data),          intent(inout) :: ice
    class(atype_ice_geometry_model_data), intent(in   ) :: geom

    ! Local variables:
    character(len=*), parameter :: routine_name = 'ice_velocity_model_run'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run stuff that is common to all models
    call self%run_model()

    ! Run stuff that is common to all ice_velocity models

    ! Run stuff that is specific to each individual ice_velocity model implementation
    call self%run_ice_velocity_model( ice, geom)

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
    call self%remap_field( mesh_new, 'u_3D'  , self%u_3D  )
    call self%remap_field( mesh_new, 'v_3D'  , self%v_3D  )
    call self%remap_field( mesh_new, 'u_3D_b', self%u_3D_b)
    call self%remap_field( mesh_new, 'v_3D_b', self%v_3D_b)
    call self%remap_field( mesh_new, 'w_3D'  , self%w_3D  )

    ! Vertically averaged
    call self%remap_field( mesh_new, 'u_vav'     , self%u_vav     )
    call self%remap_field( mesh_new, 'v_vav'     , self%v_vav     )
    call self%remap_field( mesh_new, 'u_vav_b'   , self%u_vav_b   )
    call self%remap_field( mesh_new, 'v_vav_b'   , self%v_vav_b   )
    call self%remap_field( mesh_new, 'uabs_vav'  , self%uabs_vav  )
    call self%remap_field( mesh_new, 'uabs_vav_b', self%uabs_vav_b)
    call self%remap_field( mesh_new, 'u_vav_perp', self%u_vav_perp)

    ! Surface
    call self%remap_field( mesh_new, 'u_surf'     , self%u_surf     )
    call self%remap_field( mesh_new, 'v_surf'     , self%v_surf     )
    call self%remap_field( mesh_new, 'u_surf_b'   , self%u_surf_b   )
    call self%remap_field( mesh_new, 'v_surf_b'   , self%v_surf_b   )
    call self%remap_field( mesh_new, 'w_surf'     , self%w_surf     )
    call self%remap_field( mesh_new, 'uabs_surf'  , self%uabs_surf  )
    call self%remap_field( mesh_new, 'uabs_surf_b', self%uabs_surf_b)

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