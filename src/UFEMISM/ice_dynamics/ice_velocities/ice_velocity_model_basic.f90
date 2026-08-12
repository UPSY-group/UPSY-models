module ice_velocity_model_basic

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use ice_velocity_model_data, only: atype_ice_velocity_model_data
  use mesh_types, only: type_mesh
  use parameters, only: NaN

  implicit none

  private

  public :: atype_ice_velocity_model, type_ice_velocity_model_clean

  type, abstract, extends(atype_ice_velocity_model_data) :: atype_ice_velocity_model

    contains

      procedure, public :: allocate   => allocate_ice_velocity_model
      procedure, public :: deallocate => deallocate_ice_velocity_model
      ! procedure, public :: remap      => remap_ice_velocity_model

      procedure, public :: get_model_name

  end type atype_ice_velocity_model

  ! Temporary concrete extension, to be used until all velocity models have been converted
  ! to extensions of the abstract base type
  type, extends(atype_ice_velocity_model) :: type_ice_velocity_model_clean
  end type type_ice_velocity_model_clean

  ! Interfaces for procedures defined in submodules
  interface

    ! module subroutine remap_ice_velocity_model( self, mesh_old, mesh_new)
    !   class(type_ice_velocity_model),       intent(inout) :: self
    !   type(type_mesh),                      intent(in   ) :: mesh_old
    !   type(type_mesh),                      intent(in   ) :: mesh_new
    ! end subroutine remap_ice_velocity_model

  end interface

contains

  subroutine allocate_ice_velocity_model( self, region_name, mesh)

    ! In/output variables:
    class(atype_ice_velocity_model), intent(inout) :: self
    character(len=*),                intent(in   ) :: region_name
    type(type_mesh), target,         intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_ice_velocity_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all models
    call self%allocate_model( region_name, mesh)

    ! Allocate all the stuff that is specific to the ice_velocity model

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
    allocate( self%u_base     ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%v_base     ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%u_base_b   ( mesh%ti1:mesh%ti2), source = NaN)
    allocate( self%v_base_b   ( mesh%ti1:mesh%ti2), source = NaN)
    allocate( self%w_base     ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%uabs_base  ( mesh%vi1:mesh%vi2), source = NaN)
    allocate( self%uabs_base_b( mesh%ti1:mesh%ti2), source = NaN)

    ! Strain rates
    allocate( self%du_dx_3D( mesh%vi1:mesh%vi2, 1:mesh%nz), source = NaN)
    allocate( self%du_dy_3D( mesh%vi1:mesh%vi2, 1:mesh%nz), source = NaN)
    allocate( self%du_dz_3D( mesh%vi1:mesh%vi2, 1:mesh%nz), source = NaN)
    allocate( self%dv_dx_3D( mesh%vi1:mesh%vi2, 1:mesh%nz), source = NaN)
    allocate( self%dv_dy_3D( mesh%vi1:mesh%vi2, 1:mesh%nz), source = NaN)
    allocate( self%dv_dz_3D( mesh%vi1:mesh%vi2, 1:mesh%nz), source = NaN)
    allocate( self%dw_dx_3D( mesh%vi1:mesh%vi2, 1:mesh%nz), source = NaN)
    allocate( self%dw_dy_3D( mesh%vi1:mesh%vi2, 1:mesh%nz), source = NaN)
    allocate( self%dw_dz_3D( mesh%vi1:mesh%vi2, 1:mesh%nz), source = NaN)

    ! Flow regime
    allocate( self%R_shear( mesh%vi1:mesh%vi2), source = NaN)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_ice_velocity_model

  subroutine deallocate_ice_velocity_model( self)

    ! In/output variables:
    class(atype_ice_velocity_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_ice_velocity_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate stuff that is common to all models
    call self%deallocate_model()

    ! Deallocate stuff that is specific to the ice_velocity model

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
    deallocate( self%u_base     )
    deallocate( self%v_base     )
    deallocate( self%u_base_b   )
    deallocate( self%v_base_b   )
    deallocate( self%w_base     )
    deallocate( self%uabs_base  )
    deallocate( self%uabs_base_b)

    ! Strain rates
    deallocate( self%du_dx_3D)
    deallocate( self%du_dy_3D)
    deallocate( self%du_dz_3D)
    deallocate( self%dv_dx_3D)
    deallocate( self%dv_dy_3D)
    deallocate( self%dv_dz_3D)
    deallocate( self%dw_dx_3D)
    deallocate( self%dw_dy_3D)
    deallocate( self%dw_dz_3D)

    ! Flow regime
    deallocate( self%R_shear)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_ice_velocity_model

  function get_model_name( self) result( model_name)
    class(atype_ice_velocity_model), intent(in) :: self
    character(len=:), allocatable              :: model_name
    model_name = 'ice_velocity'
  end function get_model_name

end module ice_velocity_model_basic