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
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared
  use mpi_basic, only: par

  implicit none

  private

  public :: type_ugrid

  type, extends(atype_ugrid_data) :: type_ugrid
    contains

      procedure, public  :: build_from_grid

      final :: finalise

  end type type_ugrid

  ! Interfaces for procedures defined in submodules
  interface

    module subroutine build_from_mesh( ugrid, mesh, mesh_Arakawa_grid)
      type(type_ugrid),        intent(inout) :: ugrid
      type(type_mesh),         intent(in   ) :: mesh
      type(type_Arakawa_grid), intent(in   ) :: mesh_Arakawa_grid
    end subroutine build_from_mesh

    module subroutine build_from_xy_grid( ugrid, xy_grid, xy_grid_Arakawa_grid)
      type(type_ugrid),        intent(inout) :: ugrid
      type(type_xy_grid),      intent(in   ) :: xy_grid
      type(type_Arakawa_grid), intent(in   ) :: xy_grid_Arakawa_grid
    end subroutine build_from_xy_grid

    module subroutine build_from_lonlat_grid( ugrid, lonlat_grid, lonlat_grid_Arakawa_grid)
      type(type_ugrid),        intent(inout) :: ugrid
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
        call build_from_mesh( self, g, grid_Arakawa_grid)
      class is (type_xy_grid)
        call build_from_xy_grid( self, g, grid_Arakawa_grid)
      class is (type_lonlat_grid)
        call build_from_lonlat_grid( self, g, grid_Arakawa_grid)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine build_from_grid

  subroutine finalise( ugrid)
    !< Deallocate hybrid distributed/shared components (as those are
    !< not captured by automatic finalisation)

    type(type_ugrid) :: ugrid

    if (associated( ugrid%node_coordinates        )) call deallocate_dist_shared( ugrid%node_coordinates        , ugrid%wnode_coordinates        )
    if (associated( ugrid%face_coordinates        )) call deallocate_dist_shared( ugrid%face_coordinates        , ugrid%wface_coordinates        )
    if (associated( ugrid%face_node_connectivity_n)) call deallocate_dist_shared( ugrid%face_node_connectivity_n, ugrid%wface_node_connectivity_n)
    if (associated( ugrid%face_node_connectivity  )) call deallocate_dist_shared( ugrid%face_node_connectivity  , ugrid%wface_node_connectivity  )
    if (associated( ugrid%face_face_connectivity_n)) call deallocate_dist_shared( ugrid%face_face_connectivity_n, ugrid%wface_face_connectivity_n)
    if (associated( ugrid%face_face_connectivity  )) call deallocate_dist_shared( ugrid%face_face_connectivity  , ugrid%wface_face_connectivity  )

  end subroutine finalise

end module ugrid_basic
