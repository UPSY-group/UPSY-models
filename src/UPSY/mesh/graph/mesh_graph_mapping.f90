module mesh_graph_mapping

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph
  use mpi_f08, only: MPI_WIN
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use deallocate_dist_shared_mod, only: deallocate_dist_shared
  use mpi_distributed_memory, only: gather_to_all
  use mpi_distributed_shared_memory, only: gather_dist_shared_to_all
  use mpi_basic, only: par, sync
  use parameters, only: NaN

  implicit none

  private

  public :: map_mesh_vertices_to_graph, map_mesh_triangles_to_graph, map_mesh_edges_to_graph, &
    map_graph_to_mesh_vertices, map_graph_to_mesh_triangles, map_graph_to_mesh_edges

  interface map_mesh_vertices_to_graph
    procedure :: map_mesh_vertices_to_graph_logical_2D
    procedure :: map_mesh_vertices_to_graph_logical_3D
    procedure :: map_mesh_vertices_to_graph_int_2D
    procedure :: map_mesh_vertices_to_graph_int_3D
    procedure :: map_mesh_vertices_to_graph_dp_2D
    procedure :: map_mesh_vertices_to_graph_dp_3D
  end interface map_mesh_vertices_to_graph

  interface map_mesh_triangles_to_graph
    procedure :: map_mesh_triangles_to_graph_logical_2D
    procedure :: map_mesh_triangles_to_graph_logical_3D
    procedure :: map_mesh_triangles_to_graph_int_2D
    procedure :: map_mesh_triangles_to_graph_int_3D
    procedure :: map_mesh_triangles_to_graph_dp_2D
    procedure :: map_mesh_triangles_to_graph_dp_3D
  end interface map_mesh_triangles_to_graph

  interface map_mesh_edges_to_graph
    procedure :: map_mesh_edges_to_graph_logical_2D
    procedure :: map_mesh_edges_to_graph_logical_3D
    procedure :: map_mesh_edges_to_graph_int_2D
    procedure :: map_mesh_edges_to_graph_int_3D
    procedure :: map_mesh_edges_to_graph_dp_2D
    procedure :: map_mesh_edges_to_graph_dp_3D
  end interface map_mesh_edges_to_graph

  interface map_graph_to_mesh_vertices
    procedure :: map_graph_to_mesh_vertices_logical_2D
    procedure :: map_graph_to_mesh_vertices_logical_3D
    procedure :: map_graph_to_mesh_vertices_int_2D
    procedure :: map_graph_to_mesh_vertices_int_3D
    procedure :: map_graph_to_mesh_vertices_dp_2D
    procedure :: map_graph_to_mesh_vertices_dp_3D
  end interface map_graph_to_mesh_vertices

  interface map_graph_to_mesh_triangles
    procedure :: map_graph_to_mesh_triangles_logical_2D
    procedure :: map_graph_to_mesh_triangles_logical_3D
    procedure :: map_graph_to_mesh_triangles_int_2D
    procedure :: map_graph_to_mesh_triangles_int_3D
    procedure :: map_graph_to_mesh_triangles_dp_2D
    procedure :: map_graph_to_mesh_triangles_dp_3D
  end interface map_graph_to_mesh_triangles

  interface map_graph_to_mesh_edges
    procedure :: map_graph_to_mesh_edges_logical_2D
    procedure :: map_graph_to_mesh_edges_logical_3D
    procedure :: map_graph_to_mesh_edges_int_2D
    procedure :: map_graph_to_mesh_edges_int_3D
    procedure :: map_graph_to_mesh_edges_dp_2D
    procedure :: map_graph_to_mesh_edges_dp_3D
  end interface map_graph_to_mesh_edges

contains

  ! Mesh vertices to graph
  ! ======================

  subroutine map_mesh_vertices_to_graph_logical_2D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                                       intent(in   ) :: mesh
    logical, dimension(:),                                 intent(in   ) :: d_mesh
    type(type_graph),                                      intent(in   ) :: graph
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_vertices_to_graph_logical_2D'
    logical, dimension(:), pointer :: d_mesh_tot => null()
    type(MPI_WIN)                  :: wd_mesh_tot
    integer                        :: ni, vi

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_mesh_tot, wd_mesh_tot, mesh%nV)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_V%n_loc) then
      call gather_to_all( d_mesh, d_mesh_tot)
    elseif (size( d_mesh,1) == mesh%pai_V%n_nih) then
      call gather_dist_shared_to_all( mesh%pai_V, d_mesh, d_mesh_tot)
    else
      call crash('invalid size for d_mesh')
    end if

    do ni = graph%pai%i1, graph%pai%i2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        d_graph_nih( ni) = d_mesh_tot( vi)
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_tot, wd_mesh_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_vertices_to_graph_logical_2D

  subroutine map_mesh_vertices_to_graph_logical_3D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    logical, dimension(:,:), target, intent(in   ) :: d_mesh
    type(type_graph),                intent(in   ) :: graph
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_vertices_to_graph_logical_3D'
    integer                        :: k
    logical, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_mesh_vertices_to_graph_logical_2D( mesh, d_mesh_k, graph, d_graph_nih_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_vertices_to_graph_logical_3D

  subroutine map_mesh_vertices_to_graph_int_2D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                                       intent(in   ) :: mesh
    integer, dimension(:),                                 intent(in   ) :: d_mesh
    type(type_graph),                                      intent(in   ) :: graph
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_vertices_to_graph_int_2D'
    integer, dimension(:), pointer :: d_mesh_tot => null()
    type(MPI_WIN)                  :: wd_mesh_tot
    integer                        :: ni, vi

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_mesh_tot, wd_mesh_tot, mesh%nV)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_V%n_loc) then
      call gather_to_all( d_mesh, d_mesh_tot)
    elseif (size( d_mesh,1) == mesh%pai_V%n_nih) then
      call gather_dist_shared_to_all( mesh%pai_V, d_mesh, d_mesh_tot)
    else
      call crash('invalid size for d_mesh')
    end if

    do ni = graph%pai%i1, graph%pai%i2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        d_graph_nih( ni) = d_mesh_tot( vi)
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_tot, wd_mesh_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_vertices_to_graph_int_2D

  subroutine map_mesh_vertices_to_graph_int_3D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    integer, dimension(:,:), target, intent(in   ) :: d_mesh
    type(type_graph),                intent(in   ) :: graph
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_vertices_to_graph_int_3D'
    integer                        :: k
    integer, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_mesh_vertices_to_graph_int_2D( mesh, d_mesh_k, graph, d_graph_nih_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_vertices_to_graph_int_3D

  subroutine map_mesh_vertices_to_graph_dp_2D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                                        intent(in   ) :: mesh
    real(dp), dimension(:),                                 intent(in   ) :: d_mesh
    type(type_graph),                                       intent(in   ) :: graph
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_mesh_vertices_to_graph_dp_2D'
    real(dp), dimension(:), pointer :: d_mesh_tot => null()
    type(MPI_WIN)                   :: wd_mesh_tot
    integer                         :: ni, vi

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_mesh_tot, wd_mesh_tot, mesh%nV)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_V%n_loc) then
      call gather_to_all( d_mesh, d_mesh_tot)
    elseif (size( d_mesh,1) == mesh%pai_V%n_nih) then
      call gather_dist_shared_to_all( mesh%pai_V, d_mesh, d_mesh_tot)
    else
      call crash('invalid size for d_mesh')
    end if

    do ni = graph%pai%i1, graph%pai%i2
      vi = graph%ni2vi( ni)
      if (vi > 0) then
        d_graph_nih( ni) = d_mesh_tot( vi)
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_tot, wd_mesh_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_vertices_to_graph_dp_2D

  subroutine map_mesh_vertices_to_graph_dp_3D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                  intent(in   ) :: mesh
    real(dp), dimension(:,:), target, intent(in   ) :: d_mesh
    type(type_graph),                 intent(in   ) :: graph
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),      target, intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_mesh_vertices_to_graph_dp_3D'
    integer                         :: k
    real(dp), dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_mesh_vertices_to_graph_dp_2D( mesh, d_mesh_k, graph, d_graph_nih_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_vertices_to_graph_dp_3D

  ! Mesh triangles to graph
  ! =======================

  subroutine map_mesh_triangles_to_graph_logical_2D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                                       intent(in   ) :: mesh
    logical, dimension(:),                                 intent(in   ) :: d_mesh
    type(type_graph),                                      intent(in   ) :: graph
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_triangles_to_graph_logical_2D'
    logical, dimension(:), pointer :: d_mesh_tot => null()
    type(MPI_WIN)                  :: wd_mesh_tot
    integer                        :: ni, ti

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_mesh_tot, wd_mesh_tot, mesh%nTri)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_Tri%n_loc) then
      call gather_to_all( d_mesh, d_mesh_tot)
    elseif (size( d_mesh,1) == mesh%pai_Tri%n_nih) then
      call gather_dist_shared_to_all( mesh%pai_Tri, d_mesh, d_mesh_tot)
    else
      call crash('invalid size for d_mesh')
    end if

    do ni = graph%pai%i1, graph%pai%i2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        d_graph_nih( ni) = d_mesh_tot( ti)
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_tot, wd_mesh_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_triangles_to_graph_logical_2D

  subroutine map_mesh_triangles_to_graph_logical_3D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    logical, dimension(:,:), target, intent(in   ) :: d_mesh
    type(type_graph),                intent(in   ) :: graph
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_triangles_to_graph_logical_3D'
    integer                        :: k
    logical, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_mesh_triangles_to_graph_logical_2D( mesh, d_mesh_k, graph, d_graph_nih_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_triangles_to_graph_logical_3D

  subroutine map_mesh_triangles_to_graph_int_2D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                                       intent(in   ) :: mesh
    integer, dimension(:),                                 intent(in   ) :: d_mesh
    type(type_graph),                                      intent(in   ) :: graph
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_triangles_to_graph_int_2D'
    integer, dimension(:), pointer :: d_mesh_tot => null()
    type(MPI_WIN)                  :: wd_mesh_tot
    integer                        :: ni, ti

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_mesh_tot, wd_mesh_tot, mesh%nTri)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_Tri%n_loc) then
      call gather_to_all( d_mesh, d_mesh_tot)
    elseif (size( d_mesh,1) == mesh%pai_Tri%n_nih) then
      call gather_dist_shared_to_all( mesh%pai_Tri, d_mesh, d_mesh_tot)
    else
      call crash('invalid size for d_mesh')
    end if

    do ni = graph%pai%i1, graph%pai%i2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        d_graph_nih( ni) = d_mesh_tot( ti)
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_tot, wd_mesh_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_triangles_to_graph_int_2D

  subroutine map_mesh_triangles_to_graph_int_3D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    integer, dimension(:,:), target, intent(in   ) :: d_mesh
    type(type_graph),                intent(in   ) :: graph
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_triangles_to_graph_int_3D'
    integer                        :: k
    integer, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_mesh_triangles_to_graph_int_2D( mesh, d_mesh_k, graph, d_graph_nih_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_triangles_to_graph_int_3D

  subroutine map_mesh_triangles_to_graph_dp_2D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                                        intent(in   ) :: mesh
    real(dp), dimension(:),                                 intent(in   ) :: d_mesh
    type(type_graph),                                       intent(in   ) :: graph
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_mesh_triangles_to_graph_dp_2D'
    real(dp), dimension(:), pointer :: d_mesh_tot => null()
    type(MPI_WIN)                   :: wd_mesh_tot
    integer                         :: ni, ti

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_mesh_tot, wd_mesh_tot, mesh%nTri)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_Tri%n_loc) then
      call gather_to_all( d_mesh, d_mesh_tot)
    elseif (size( d_mesh,1) == mesh%pai_Tri%n_nih) then
      call gather_dist_shared_to_all( mesh%pai_Tri, d_mesh, d_mesh_tot)
    else
      call crash('invalid size for d_mesh')
    end if

    do ni = graph%pai%i1, graph%pai%i2
      ti = graph%ni2ti( ni)
      if (ti > 0) then
        d_graph_nih( ni) = d_mesh_tot( ti)
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_tot, wd_mesh_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_triangles_to_graph_dp_2D

  subroutine map_mesh_triangles_to_graph_dp_3D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                  intent(in   ) :: mesh
    real(dp), dimension(:,:), target, intent(in   ) :: d_mesh
    type(type_graph),                 intent(in   ) :: graph
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),      target, intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_mesh_triangles_to_graph_dp_3D'
    integer                         :: k
    real(dp), dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_mesh_triangles_to_graph_dp_2D( mesh, d_mesh_k, graph, d_graph_nih_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_triangles_to_graph_dp_3D

  ! Mesh edges to graph
  ! =======================

  subroutine map_mesh_edges_to_graph_logical_2D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                                       intent(in   ) :: mesh
    logical, dimension(:),                                 intent(in   ) :: d_mesh
    type(type_graph),                                      intent(in   ) :: graph
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_edges_to_graph_logical_2D'
    logical, dimension(:), pointer :: d_mesh_tot => null()
    type(MPI_WIN)                  :: wd_mesh_tot
    integer                        :: ni, ei

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_mesh_tot, wd_mesh_tot, mesh%nE)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_E%n_loc) then
      call gather_to_all( d_mesh, d_mesh_tot)
    elseif (size( d_mesh,1) == mesh%pai_E%n_nih) then
      call gather_dist_shared_to_all( mesh%pai_E, d_mesh, d_mesh_tot)
    else
      call crash('invalid size for d_mesh')
    end if

    do ni = graph%pai%i1, graph%pai%i2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        d_graph_nih( ni) = d_mesh_tot( ei)
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_tot, wd_mesh_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_edges_to_graph_logical_2D

  subroutine map_mesh_edges_to_graph_logical_3D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    logical, dimension(:,:), target, intent(in   ) :: d_mesh
    type(type_graph),                intent(in   ) :: graph
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_edges_to_graph_logical_3D'
    integer                        :: k
    logical, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_mesh_edges_to_graph_logical_2D( mesh, d_mesh_k, graph, d_graph_nih_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_edges_to_graph_logical_3D

  subroutine map_mesh_edges_to_graph_int_2D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                                       intent(in   ) :: mesh
    integer, dimension(:),                                 intent(in   ) :: d_mesh
    type(type_graph),                                      intent(in   ) :: graph
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_edges_to_graph_int_2D'
    integer, dimension(:), pointer :: d_mesh_tot => null()
    type(MPI_WIN)                  :: wd_mesh_tot
    integer                        :: ni, ei

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_mesh_tot, wd_mesh_tot, mesh%nE)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_E%n_loc) then
      call gather_to_all( d_mesh, d_mesh_tot)
    elseif (size( d_mesh,1) == mesh%pai_E%n_nih) then
      call gather_dist_shared_to_all( mesh%pai_E, d_mesh, d_mesh_tot)
    else
      call crash('invalid size for d_mesh')
    end if

    do ni = graph%pai%i1, graph%pai%i2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        d_graph_nih( ni) = d_mesh_tot( ei)
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_tot, wd_mesh_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_edges_to_graph_int_2D

  subroutine map_mesh_edges_to_graph_int_3D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    integer, dimension(:,:), target, intent(in   ) :: d_mesh
    type(type_graph),                intent(in   ) :: graph
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_edges_to_graph_int_3D'
    integer                        :: k
    integer, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_mesh_edges_to_graph_int_2D( mesh, d_mesh_k, graph, d_graph_nih_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_edges_to_graph_int_3D

  subroutine map_mesh_edges_to_graph_dp_2D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                                        intent(in   ) :: mesh
    real(dp), dimension(:),                                 intent(in   ) :: d_mesh
    type(type_graph),                                       intent(in   ) :: graph
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_mesh_edges_to_graph_dp_2D'
    real(dp), dimension(:), pointer :: d_mesh_tot => null()
    type(MPI_WIN)                   :: wd_mesh_tot
    integer                         :: ni, ei

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_mesh_tot, wd_mesh_tot, mesh%nE)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_E%n_loc) then
      call gather_to_all( d_mesh, d_mesh_tot)
    elseif (size( d_mesh,1) == mesh%pai_E%n_nih) then
      call gather_dist_shared_to_all( mesh%pai_E, d_mesh, d_mesh_tot)
    else
      call crash('invalid size for d_mesh')
    end if

    do ni = graph%pai%i1, graph%pai%i2
      ei = graph%ni2ei( ni)
      if (ei > 0) then
        d_graph_nih( ni) = d_mesh_tot( ei)
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_tot, wd_mesh_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_edges_to_graph_dp_2D

  subroutine map_mesh_edges_to_graph_dp_3D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                  intent(in   ) :: mesh
    real(dp), dimension(:,:), target, intent(in   ) :: d_mesh
    type(type_graph),                 intent(in   ) :: graph
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),      target, intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_mesh_edges_to_graph_dp_3D'
    integer                         :: k
    real(dp), dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_mesh_edges_to_graph_dp_2D( mesh, d_mesh_k, graph, d_graph_nih_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_edges_to_graph_dp_3D

  ! Graph to mesh vertices
  ! ======================

  subroutine map_graph_to_mesh_vertices_logical_2D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                                      intent(in   ) :: graph
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(in   ) :: d_graph_nih
    type(type_mesh),                                       intent(in   ) :: mesh
    logical, dimension(:), target,                         intent(  out) :: d_mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_graph_to_mesh_vertices_logical_2D'
    logical, dimension(:), pointer :: d_graph_tot => null()
    type(MPI_WIN)                  :: wd_graph_tot
    logical, dimension(:), pointer :: d_mesh_nih, d_mesh_loc
    integer                        :: ni, vi

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_graph_tot, wd_graph_tot, graph%n)
    call gather_dist_shared_to_all( graph%pai, d_graph_nih, d_graph_tot)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_V%n_loc) then
      d_mesh_loc( mesh%vi1: mesh%vi2) => d_mesh
    elseif (size( d_mesh,1) == mesh%pai_V%n_nih) then
      d_mesh_nih( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh
      d_mesh_loc( mesh%vi1: mesh%vi2) => d_mesh_nih( mesh%vi1: mesh%vi2)
    else
      call crash('invalid size for d_mesh')
    end if

    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        d_mesh_loc( vi) = d_graph_tot( ni)
      else
        d_mesh_loc( vi) = .false.
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_tot, wd_graph_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_vertices_logical_2D

  subroutine map_graph_to_mesh_vertices_logical_3D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                intent(in   ) :: graph
    logical, dimension(:,:), target, intent(  out) :: d_mesh
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(in   ) :: d_graph_nih
    type(type_mesh),                 intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_graph_to_mesh_vertices_logical_3D'
    integer                        :: k
    logical, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_graph_to_mesh_vertices_logical_2D( graph, d_graph_nih_k, mesh, d_mesh_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_vertices_logical_3D

  subroutine map_graph_to_mesh_vertices_int_2D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                                      intent(in   ) :: graph
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(in   ) :: d_graph_nih
    type(type_mesh),                                       intent(in   ) :: mesh
    integer, dimension(:), target,                         intent(  out) :: d_mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_graph_to_mesh_vertices_int_2D'
    integer, dimension(:), pointer :: d_graph_tot => null()
    type(MPI_WIN)                  :: wd_graph_tot
    integer, dimension(:), pointer :: d_mesh_nih, d_mesh_loc
    integer                        :: ni, vi

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_graph_tot, wd_graph_tot, graph%n)
    call gather_dist_shared_to_all( graph%pai, d_graph_nih, d_graph_tot)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_V%n_loc) then
      d_mesh_loc( mesh%vi1: mesh%vi2) => d_mesh
    elseif (size( d_mesh,1) == mesh%pai_V%n_nih) then
      d_mesh_nih( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh
      d_mesh_loc( mesh%vi1: mesh%vi2) => d_mesh_nih( mesh%vi1: mesh%vi2)
    else
      call crash('invalid size for d_mesh')
    end if

    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        d_mesh_loc( vi) = d_graph_tot( ni)
      else
        d_mesh_loc( vi) = 0
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_tot, wd_graph_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_vertices_int_2D

  subroutine map_graph_to_mesh_vertices_int_3D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                intent(in   ) :: graph
    integer, dimension(:,:), target, intent(  out) :: d_mesh
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(in   ) :: d_graph_nih
    type(type_mesh),                 intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_graph_to_mesh_vertices_int_3D'
    integer                        :: k
    integer, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_graph_to_mesh_vertices_int_2D( graph, d_graph_nih_k, mesh, d_mesh_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_vertices_int_3D

  subroutine map_graph_to_mesh_vertices_dp_2D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                                       intent(in   ) :: graph
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(in   ) :: d_graph_nih
    type(type_mesh),                                        intent(in   ) :: mesh
    real(dp), dimension(:), target,                         intent(  out) :: d_mesh

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_graph_to_mesh_vertices_dp_2D'
    real(dp), dimension(:), pointer :: d_graph_tot => null()
    type(MPI_WIN)                   :: wd_graph_tot
    real(dp), dimension(:), pointer :: d_mesh_nih, d_mesh_loc
    integer                         :: ni, vi

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_graph_tot, wd_graph_tot, graph%n)
    call gather_dist_shared_to_all( graph%pai, d_graph_nih, d_graph_tot)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_V%n_loc) then
      d_mesh_loc( mesh%vi1: mesh%vi2) => d_mesh
    elseif (size( d_mesh,1) == mesh%pai_V%n_nih) then
      d_mesh_nih( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh
      d_mesh_loc( mesh%vi1: mesh%vi2) => d_mesh_nih( mesh%vi1: mesh%vi2)
    else
      call crash('invalid size for d_mesh')
    end if

    do vi = mesh%vi1, mesh%vi2
      ni = graph%vi2ni( vi)
      if (ni > 0) then
        d_mesh_loc( vi) = d_graph_tot( ni)
      else
        d_mesh_loc( vi) = NaN
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_tot, wd_graph_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_vertices_dp_2D

  subroutine map_graph_to_mesh_vertices_dp_3D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                 intent(in   ) :: graph
    real(dp), dimension(:,:), target, intent(  out) :: d_mesh
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),      target, intent(in   ) :: d_graph_nih
    type(type_mesh),                  intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_graph_to_mesh_vertices_dp_3D'
    integer                         :: k
    real(dp), dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_graph_to_mesh_vertices_dp_2D( graph, d_graph_nih_k, mesh, d_mesh_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_vertices_dp_3D

  ! Graph to mesh triangles
  ! =======================

  subroutine map_graph_to_mesh_triangles_logical_2D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                                      intent(in   ) :: graph
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(in   ) :: d_graph_nih
    type(type_mesh),                                       intent(in   ) :: mesh
    logical, dimension(:), target,                         intent(  out) :: d_mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_graph_to_mesh_triangles_logical_2D'
    logical, dimension(:), pointer :: d_graph_tot => null()
    type(MPI_WIN)                  :: wd_graph_tot
    logical, dimension(:), pointer :: d_mesh_nih, d_mesh_loc
    integer                        :: ni, ti

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_graph_tot, wd_graph_tot, graph%n)
    call gather_dist_shared_to_all( graph%pai, d_graph_nih, d_graph_tot)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_Tri%n_loc) then
      d_mesh_loc( mesh%ti1: mesh%ti2) => d_mesh
    elseif (size( d_mesh,1) == mesh%pai_Tri%n_nih) then
      d_mesh_nih( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh
      d_mesh_loc( mesh%ti1: mesh%ti2) => d_mesh_nih( mesh%ti1: mesh%ti2)
    else
      call crash('invalid size for d_mesh')
    end if

    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        d_mesh_loc( ti) = d_graph_tot( ni)
      else
        d_mesh_loc( ti) = .false.
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_tot, wd_graph_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_triangles_logical_2D

  subroutine map_graph_to_mesh_triangles_logical_3D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                intent(in   ) :: graph
    logical, dimension(:,:), target, intent(  out) :: d_mesh
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(in   ) :: d_graph_nih
    type(type_mesh),                 intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_graph_to_mesh_triangles_logical_3D'
    integer                        :: k
    logical, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_graph_to_mesh_triangles_logical_2D( graph, d_graph_nih_k, mesh, d_mesh_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_triangles_logical_3D

  subroutine map_graph_to_mesh_triangles_int_2D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                                      intent(in   ) :: graph
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(in   ) :: d_graph_nih
    type(type_mesh),                                       intent(in   ) :: mesh
    integer, dimension(:), target,                         intent(  out) :: d_mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_graph_to_mesh_triangles_int_2D'
    integer, dimension(:), pointer :: d_graph_tot => null()
    type(MPI_WIN)                  :: wd_graph_tot
    integer, dimension(:), pointer :: d_mesh_nih, d_mesh_loc
    integer                        :: ni, ti

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_graph_tot, wd_graph_tot, graph%n)
    call gather_dist_shared_to_all( graph%pai, d_graph_nih, d_graph_tot)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_Tri%n_loc) then
      d_mesh_loc( mesh%ti1: mesh%ti2) => d_mesh
    elseif (size( d_mesh,1) == mesh%pai_Tri%n_nih) then
      d_mesh_nih( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh
      d_mesh_loc( mesh%ti1: mesh%ti2) => d_mesh_nih( mesh%ti1: mesh%ti2)
    else
      call crash('invalid size for d_mesh')
    end if

    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        d_mesh_loc( ti) = d_graph_tot( ni)
      else
        d_mesh_loc( ti) = 0
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_tot, wd_graph_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_triangles_int_2D

  subroutine map_graph_to_mesh_triangles_int_3D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                intent(in   ) :: graph
    integer, dimension(:,:), target, intent(  out) :: d_mesh
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(in   ) :: d_graph_nih
    type(type_mesh),                 intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_graph_to_mesh_triangles_int_3D'
    integer                        :: k
    integer, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_graph_to_mesh_triangles_int_2D( graph, d_graph_nih_k, mesh, d_mesh_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_triangles_int_3D

  subroutine map_graph_to_mesh_triangles_dp_2D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                                       intent(in   ) :: graph
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(in   ) :: d_graph_nih
    type(type_mesh),                                        intent(in   ) :: mesh
    real(dp), dimension(:), target,                         intent(  out) :: d_mesh

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_graph_to_mesh_triangles_dp_2D'
    real(dp), dimension(:), pointer :: d_graph_tot => null()
    type(MPI_WIN)                   :: wd_graph_tot
    real(dp), dimension(:), pointer :: d_mesh_nih, d_mesh_loc
    integer                         :: ni, ti

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_graph_tot, wd_graph_tot, graph%n)
    call gather_dist_shared_to_all( graph%pai, d_graph_nih, d_graph_tot)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_Tri%n_loc) then
      d_mesh_loc( mesh%ti1: mesh%ti2) => d_mesh
    elseif (size( d_mesh,1) == mesh%pai_Tri%n_nih) then
      d_mesh_nih( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh
      d_mesh_loc( mesh%ti1: mesh%ti2) => d_mesh_nih( mesh%ti1: mesh%ti2)
    else
      call crash('invalid size for d_mesh')
    end if

    do ti = mesh%ti1, mesh%ti2
      ni = graph%ti2ni( ti)
      if (ni > 0) then
        d_mesh_loc( ti) = d_graph_tot( ni)
      else
        d_mesh_loc( ti) = NaN
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_tot, wd_graph_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_triangles_dp_2D

  subroutine map_graph_to_mesh_triangles_dp_3D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                 intent(in   ) :: graph
    real(dp), dimension(:,:), target, intent(  out) :: d_mesh
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),      target, intent(in   ) :: d_graph_nih
    type(type_mesh),                  intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_graph_to_mesh_triangles_dp_3D'
    integer                         :: k
    real(dp), dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_graph_to_mesh_triangles_dp_2D( graph, d_graph_nih_k, mesh, d_mesh_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_triangles_dp_3D

  ! Graph to mesh edges
  ! =======================

  subroutine map_graph_to_mesh_edges_logical_2D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                                      intent(in   ) :: graph
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(in   ) :: d_graph_nih
    type(type_mesh),                                       intent(in   ) :: mesh
    logical, dimension(:), target,                         intent(  out) :: d_mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_graph_to_mesh_edges_logical_2D'
    logical, dimension(:), pointer :: d_graph_tot => null()
    type(MPI_WIN)                  :: wd_graph_tot
    logical, dimension(:), pointer :: d_mesh_nih, d_mesh_loc
    integer                        :: ni, ei

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_graph_tot, wd_graph_tot, graph%n)
    call gather_dist_shared_to_all( graph%pai, d_graph_nih, d_graph_tot)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_E%n_loc) then
      d_mesh_loc( mesh%ei1: mesh%ei2) => d_mesh
    elseif (size( d_mesh,1) == mesh%pai_E%n_nih) then
      d_mesh_nih( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih) => d_mesh
      d_mesh_loc( mesh%ei1: mesh%ei2) => d_mesh_nih( mesh%ei1: mesh%ei2)
    else
      call crash('invalid size for d_mesh')
    end if

    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        d_mesh_loc( ei) = d_graph_tot( ni)
      else
        d_mesh_loc( ei) = .false.
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_tot, wd_graph_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_edges_logical_2D

  subroutine map_graph_to_mesh_edges_logical_3D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                intent(in   ) :: graph
    logical, dimension(:,:), target, intent(  out) :: d_mesh
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(in   ) :: d_graph_nih
    type(type_mesh),                 intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_graph_to_mesh_edges_logical_3D'
    integer                        :: k
    logical, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_graph_to_mesh_edges_logical_2D( graph, d_graph_nih_k, mesh, d_mesh_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_edges_logical_3D

  subroutine map_graph_to_mesh_edges_int_2D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                                      intent(in   ) :: graph
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(in   ) :: d_graph_nih
    type(type_mesh),                                       intent(in   ) :: mesh
    integer, dimension(:), target,                         intent(  out) :: d_mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_graph_to_mesh_edges_int_2D'
    integer, dimension(:), pointer :: d_graph_tot => null()
    type(MPI_WIN)                  :: wd_graph_tot
    integer, dimension(:), pointer :: d_mesh_nih, d_mesh_loc
    integer                        :: ni, ei

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_graph_tot, wd_graph_tot, graph%n)
    call gather_dist_shared_to_all( graph%pai, d_graph_nih, d_graph_tot)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_E%n_loc) then
      d_mesh_loc( mesh%ei1: mesh%ei2) => d_mesh
    elseif (size( d_mesh,1) == mesh%pai_E%n_nih) then
      d_mesh_nih( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih) => d_mesh
      d_mesh_loc( mesh%ei1: mesh%ei2) => d_mesh_nih( mesh%ei1: mesh%ei2)
    else
      call crash('invalid size for d_mesh')
    end if

    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        d_mesh_loc( ei) = d_graph_tot( ni)
      else
        d_mesh_loc( ei) = 0
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_tot, wd_graph_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_edges_int_2D

  subroutine map_graph_to_mesh_edges_int_3D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                intent(in   ) :: graph
    integer, dimension(:,:), target, intent(  out) :: d_mesh
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(in   ) :: d_graph_nih
    type(type_mesh),                 intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_graph_to_mesh_edges_int_3D'
    integer                        :: k
    integer, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_graph_to_mesh_edges_int_2D( graph, d_graph_nih_k, mesh, d_mesh_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_edges_int_3D

  subroutine map_graph_to_mesh_edges_dp_2D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                                       intent(in   ) :: graph
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(in   ) :: d_graph_nih
    type(type_mesh),                                        intent(in   ) :: mesh
    real(dp), dimension(:), target,                         intent(  out) :: d_mesh

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_graph_to_mesh_edges_dp_2D'
    real(dp), dimension(:), pointer :: d_graph_tot => null()
    type(MPI_WIN)                   :: wd_graph_tot
    real(dp), dimension(:), pointer :: d_mesh_nih, d_mesh_loc
    integer                         :: ni, ei

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_graph_tot, wd_graph_tot, graph%n)
    call gather_dist_shared_to_all( graph%pai, d_graph_nih, d_graph_tot)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_E%n_loc) then
      d_mesh_loc( mesh%ei1: mesh%ei2) => d_mesh
    elseif (size( d_mesh,1) == mesh%pai_E%n_nih) then
      d_mesh_nih( mesh%pai_E%i1_nih: mesh%pai_E%i2_nih) => d_mesh
      d_mesh_loc( mesh%ei1: mesh%ei2) => d_mesh_nih( mesh%ei1: mesh%ei2)
    else
      call crash('invalid size for d_mesh')
    end if

    do ei = mesh%ei1, mesh%ei2
      ni = graph%ei2ni( ei)
      if (ni > 0) then
        d_mesh_loc( ei) = d_graph_tot( ni)
      else
        d_mesh_loc( ei) = NaN
      end if
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_tot, wd_graph_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_edges_dp_2D

  subroutine map_graph_to_mesh_edges_dp_3D( graph, d_graph_nih, mesh, d_mesh)

    ! In/output variables:
    type(type_graph),                 intent(in   ) :: graph
    real(dp), dimension(:,:), target, intent(  out) :: d_mesh
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),      target, intent(in   ) :: d_graph_nih
    type(type_mesh),                  intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_graph_to_mesh_edges_dp_3D'
    integer                         :: k
    real(dp), dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_graph_to_mesh_edges_dp_2D( graph, d_graph_nih_k, mesh, d_mesh_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_graph_to_mesh_edges_dp_3D

end module mesh_graph_mapping