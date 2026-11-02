module ugrid_basic

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use ugrid_data, only: atype_ugrid_data
  use basic_grid_mod, only: atype_grid
  use Arakawa_grid_mod, only: type_Arakawa_grid, Arakawa_grid
  use mesh_basic, only: type_mesh
  use xy_grid_basic, only: type_xy_grid
  use lonlat_grid_basic, only: type_lonlat_grid

  implicit none

  private

  public :: type_ugrid

  type, extends(atype_ugrid_data) :: type_ugrid
    contains

      procedure, public  :: build_from_grid
      procedure, private :: build_from_mesh
      procedure, private :: build_from_xy_grid
      procedure, private :: build_from_lonlat_grid

  end type type_ugrid

  ! Interfaces for procedures defined in submodules
  interface

    module subroutine build_from_mesh( self, mesh, mesh_Arakawa_grid)
      class(type_ugrid),       intent(inout) :: self
      type(type_mesh),         intent(in   ) :: mesh
      type(type_Arakawa_grid), intent(in   ) :: mesh_Arakawa_grid
    end subroutine build_from_mesh

    module subroutine build_from_xy_grid( self, xy_grid, xy_grid_Arakawa_grid)
      class(type_ugrid),       intent(inout) :: self
      type(type_xy_grid),      intent(in   ) :: xy_grid
      type(type_Arakawa_grid), intent(in   ) :: xy_grid_Arakawa_grid
    end subroutine build_from_xy_grid

    module subroutine build_from_lonlat_grid( self, lonlat_grid, lonlat_grid_Arakawa_grid)
      class(type_ugrid),       intent(inout) :: self
      type(type_lonlat_grid),  intent(in   ) :: lonlat_grid
      type(type_Arakawa_grid), intent(in   ) :: lonlat_grid_Arakawa_grid
    end subroutine build_from_lonlat_grid

  end interface

contains

  subroutine build_from_grid( self, grid, grid_Arakawa_grid)
    ! Build a ugrid object from an existing mesh, x/y-grid, or lon/lat-grid object

    ! In/output variables:
    class(type_ugrid),       intent(inout) :: self
    class(atype_grid),       intent(in   ) :: grid
    type(type_Arakawa_grid), intent(in   ) :: grid_Arakawa_grid

    ! Local variables:
    character(len=*), parameter :: routine_name = 'build_from_grid'

    ! Add routine to path
    call init_routine( routine_name)

    select type (g => grid)
      class default
        call crash('unsupported grid type')
      class is (type_mesh)
        call self%build_from_mesh( g, grid_Arakawa_grid)
      class is (type_xy_grid)
        call self%build_from_xy_grid( g, grid_Arakawa_grid)
      class is (type_lonlat_grid)
        call self%build_from_lonlat_grid( g, grid_Arakawa_grid)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine build_from_grid

end module ugrid_basic
