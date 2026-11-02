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

  subroutine build_from_mesh_vertices( ugrid, mesh)
    ! Build a ugrid from the mesh vertices (i.e. the Voronoi mesh)

    ! In/output variables:
    type(type_ugrid), intent(inout) :: ugrid
    type(type_mesh),  intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'build_from_mesh_vertices'
    integer                     :: face_node_connectivity_n_max
    integer                     :: face_face_connectivity_n_max
    integer                     :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Number of nodes and faces in the ugrid
    ugrid%nNodes = mesh%nVor
    ugrid%nFaces = mesh%nV

    ! x/y-coordinates of nodes and faces
    call allocate_dist_shared( ugrid%node_coordinates, ugrid%wnode_coordinates, [1, ugrid%nNodes], [1, 2])
    call allocate_dist_shared( ugrid%face_coordinates, ugrid%wface_coordinates, [1, ugrid%nFaces], [1, 2])
    if (par%node_primary) then
      ugrid%node_coordinates = mesh%Vor
      ugrid%face_coordinates = mesh%V
    end if

    ! Face-node connectivity
    ugrid%face_node_connectivity_fill_value = -999
    face_node_connectivity_n_max = maxval( mesh%nVVor)

    call allocate_dist_shared( ugrid%face_node_connectivity_n, ugrid%wface_node_connectivity_n, &
      [1, ugrid%nFaces])
    call allocate_dist_shared( ugrid%face_node_connectivity, ugrid%wface_node_connectivity, &
      [1, ugrid%nFaces], [1, face_node_connectivity_n_max])

    if (par%node_primary) then
      ugrid%face_node_connectivity_n = mesh%nVVor
      ugrid%face_node_connectivity = ugrid%face_node_connectivity_fill_value
      do vi = 1, mesh%nV
        ugrid%face_node_connectivity( vi, 1: mesh%nVVor( vi)) = mesh%VVor( vi, 1: mesh%nVVor( vi))
      end do
    end if

    ! Face-face connectivity
    ugrid%face_face_connectivity_fill_value = -999
    face_face_connectivity_n_max = maxval( mesh%nC)

    call allocate_dist_shared( ugrid%face_face_connectivity_n, ugrid%wface_face_connectivity_n, &
      [1, ugrid%nFaces])
    call allocate_dist_shared( ugrid%face_face_connectivity, ugrid%wface_face_connectivity, &
      [1, ugrid%nFaces], [1, face_face_connectivity_n_max])

    if (par%node_primary) then
      ugrid%face_face_connectivity_n = mesh%nC
      ugrid%face_face_connectivity = ugrid%face_face_connectivity_fill_value
      do vi = 1, mesh%nV
        ugrid%face_face_connectivity( vi, 1: mesh%nC( vi)) = mesh%C( vi, 1: mesh%nC( vi))
      end do
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine build_from_mesh_vertices

end submodule build_from_mesh_mod
