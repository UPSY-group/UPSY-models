submodule( ugrid_basic) build_from_mesh_mod

contains

  subroutine build_from_mesh( self, mesh, mesh_Arakawa_grid)

    ! In/output variables:
    class(type_ugrid),       intent(inout) :: self
    type(type_mesh),         intent(in   ) :: mesh
    type(type_Arakawa_grid), intent(in   ) :: mesh_Arakawa_grid

    ! Local variables:
    character(len=*), parameter :: routine_name = 'build_from_mesh'

    ! Add routine to path
    call init_routine( routine_name)

    ! Build from the mesh a-grid (vertices), b-grid (triangles), or c-grid (edges)
    if     (mesh_Arakawa_grid == Arakawa_grid%a()) then
      call build_from_mesh_vertices( self, mesh)
    elseif (mesh_Arakawa_grid == Arakawa_grid%b()) then
      call crash('fixme!')
    elseif (mesh_Arakawa_grid == Arakawa_grid%c()) then
      call crash('fixme!')
    else
      call crash('unsupported Arakawa grid')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine build_from_mesh

  subroutine build_from_mesh_vertices( self, mesh)

    ! In/output variables:
    class(type_ugrid),       intent(inout) :: self
    type(type_mesh),         intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'build_from_mesh_vertices'

    ! Add routine to path
    call init_routine( routine_name)

    call crash('fixme!')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine build_from_mesh_vertices

end submodule build_from_mesh_mod
