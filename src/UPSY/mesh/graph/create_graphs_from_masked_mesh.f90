module create_graphs_from_masked_mesh

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph
  use graph_memory, only: allocate_graph_primary, crop_graph_primary
  use mpi_distributed_memory, only: gather_to_all
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared
  use plane_geometry, only: projection_of_p_on_qr
  use assertions_basic, only: assert
  use graph_contiguous_domains, only: enforce_contiguous_process_domains_graph
  use graph_parallelisation, only: setup_graph_parallelisation
  use tests_main, only: test_graph_is_self_consistent
  use parameters, only: NaN
  use mpi_basic, only: par
  use plane_geometry, only: cross2

  implicit none

  private

  public :: create_graph_from_masked_mesh_a, create_graph_from_masked_mesh_b

contains

  subroutine create_graph_from_masked_mesh_a( mesh, mask_a, nz, graph)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_a
    integer,                               intent(in   ) :: nz
    type(type_graph),                      intent(  out) :: graph

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_graph_from_masked_mesh_a'
    logical, dimension(1:mesh%nV  ) :: mask_a_tot
    integer                         :: ni, vi, ci, vj, nj

    ! Add routine to path
    call init_routine( routine_name)

    call gather_to_all( mask_a, mask_a_tot)

    ! Allocate graph and set metadata
    call allocate_graph_primary( graph, mesh%nV, mesh%nV, mesh%nTri, mesh%nE, maxval( mesh%nC))

    graph%parent_mesh_name = trim( mesh%name) // '_vertices'
    graph%xmin             = mesh%xmin
    graph%xmax             = mesh%xmax
    graph%ymin             = mesh%ymin
    graph%ymax             = mesh%ymax

    ! Create vertex-to-node mapping
    graph%n = 0
    do vi = 1, mesh%nV
      if (mask_a_tot( vi)) then
        graph%n = graph%n + 1
        graph%vi2ni( vi     ) = graph%n
        graph%ni2vi( graph%n) = vi
      end if
    end do

    ! Mark regular and border nodes
    graph%is_border = .false.
    do ni = 1, graph%n
      vi = graph%ni2vi( ni)
      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        if (.not. mask_a_tot( vj)) graph%is_border( ni) = .true.
      end do
    end do

    ! Calculate node coordinates and connectivity
    do ni = 1, graph%n

      vi = graph%ni2vi( ni)

      graph%V( ni,:) = mesh%V( vi,:)

      graph%nC( ni) = 0
      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        if (mask_a_tot( vj)) then
          nj = graph%vi2ni( vj)
          graph%nC( ni) = graph%nC( ni) + 1
          graph%C( ni, graph%nC( ni)) = nj
        end if
      end do

    end do

    ! Ensure border node connectivity is sorted counter-clockwise, without passing
    ! outside the graph interior (i.e. the same as for mesh border vertices)
    call sort_border_nodes_connectivity( graph)

    ! Finalisation
    call crop_graph_primary( graph)
    call enforce_contiguous_process_domains_graph( graph)
    call setup_graph_parallelisation( graph, nz)

#if (DO_ASSERTIONS)
    call assert(test_graph_is_self_consistent( mesh, graph), 'inconsistent graph connectivity')
#endif

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 4)

  end subroutine create_graph_from_masked_mesh_a

  subroutine create_graph_from_masked_mesh_b( mesh, mask_a, nz, graph)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_a
    integer,                               intent(in   ) :: nz
    type(type_graph),                      intent(  out) :: graph

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_graph_from_masked_mesh_b'
    logical, dimension(1:mesh%nV)   :: mask_a_tot, mask_border_a_tot
    logical, dimension(1:mesh%nTri) :: mask_b_tot, mask_border_b_tot
    integer                         :: ti, n, tj, nj, ei
    integer                         :: ni_reg, ni_ghost
    integer                         :: ni, vi, vj, ci, ej, nj1, nj2
    real(dp), dimension(2)          :: r, s, outward_normal_vector

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate b-grid mask
    call calc_extra_masks( mesh, mask_a, &
      mask_a_tot, mask_border_a_tot, mask_b_tot, mask_border_b_tot)

    ! Allocate graph and set metadata
    call allocate_graph_primary( graph, mesh%nTri*3, mesh%nV, mesh%nTri, mesh%nE, 3)

    graph%parent_mesh_name = trim( mesh%name) // '_triangles'
    graph%xmin             = mesh%xmin
    graph%xmax             = mesh%xmax
    graph%ymin             = mesh%ymin
    graph%ymax             = mesh%ymax

    ! Create triangle-to-node mapping
    graph%n = 0
    do ti = 1, mesh%nTri
      if (mask_b_tot( ti)) then
        graph%n = graph%n + 1
        graph%ti2ni( ti     ) = graph%n
        graph%ni2ti( graph%n) = ti
      end if
    end do

    ! Mark regular and border nodes
    graph%is_border = .false.
    do ni = 1, graph%n
      ti = graph%ni2ti( ni)
      if (mask_border_b_tot( ti)) graph%is_border( ni) = .true.
    end do

    ! Calculate node coordinates and connectivity
    do ni = 1, graph%n

      ti = graph%ni2ti( ni)

      graph%V( ni,:) = mesh%Trigc( ti,:)

      graph%nC( ni) = 0
      do n = 1, 3
        tj = mesh%TriC( ti,n)
        if (tj == 0) cycle
        if (mask_b_tot( tj)) then
          nj = graph%ti2ni( tj)
          graph%nC( ni) = graph%nC( ni) + 1
          graph%C( ni, graph%nC( ni)) = nj
        end if
      end do

    end do

    ! Ensure border node connectivity is sorted counter-clockwise, without passing
    ! outside the graph interior (i.e. the same as for mesh border vertices)
    call sort_border_nodes_connectivity( graph)
    call calc_outward_unit_normal_vectors_at_border_triangles( graph, mesh, &
      mask_a_tot, mask_border_a_tot, mask_b_tot, mask_border_b_tot)

    ! Finalisation
    call crop_graph_primary( graph)
    call enforce_contiguous_process_domains_graph( graph)
    call setup_graph_parallelisation( graph, nz)

#if (DO_ASSERTIONS)
    call assert(test_graph_is_self_consistent( mesh, graph), 'inconsistent graph connectivity')
#endif

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 4)

  end subroutine create_graph_from_masked_mesh_b

  subroutine calc_extra_masks( mesh, mask_a, &
      mask_a_tot, mask_border_a_tot, mask_b_tot, mask_border_b_tot)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_a
    logical, dimension(mesh%nV),           intent(  out) :: mask_a_tot, mask_border_a_tot
    logical, dimension(mesh%nTri),         intent(  out) :: mask_b_tot, mask_border_b_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_extra_masks'
    integer                        :: vi, ci, vj, ti, n_a, n

    ! Add routine to path
    call init_routine( routine_name)

    ! a-grid mask
    call gather_to_all( mask_a, mask_a_tot)

    ! a-grid border mask
    mask_border_a_tot = .false.
    do vi = 1, mesh%nV
      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        if (.not. mask_a_tot( vj)) mask_border_a_tot( vi) = .true.
      end do
    end do

    ! b-grid mask
    mask_b_tot = .false.
    do ti = 1, mesh%nTri
      n_a = 0
      do n = 1, 3
        vi = mesh%Tri( ti,n)
        if (mask_a_tot( vi)) n_a = n_a + 1
      end do
      if (n_a == 3) then
        mask_b_tot( ti) = .true.
      end if
    end do

    ! b-grid border mask
    mask_border_b_tot = .false.
    do ti = 1, mesh%nTri
      n_a = 0
      do n = 1, 3
        vi = mesh%Tri( ti,n)
        if (mask_border_a_tot( vi)) n_a = n_a + 1
      end do
      if (n_a >= 2) then
        mask_border_b_tot( ti) = .true.
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_extra_masks

  subroutine sort_border_nodes_connectivity( graph)
    ! Ensure border node connectivity is sorted counter-clockwise, without passing
    ! outside the graph interior (i.e. the same as for mesh border vertices)

    ! In/output variables:
    type(type_graph), intent(inout) :: graph

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'sort_border_nodes_connectivity'
    integer                        :: ni

    ! Add routine to path
    call init_routine( routine_name)

    do ni = 1, graph%n
      if (.not. graph%is_border( ni)) cycle

      if (graph%nC( ni) == 2) then
        call sort_border_node_connectivity_nCeq2( graph, ni)
      else
        call sort_border_node_connectivity_nCgt2( graph, ni)
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine sort_border_nodes_connectivity

  subroutine sort_border_node_connectivity_nCeq2( graph, ni)

    ! In/output variables:
    type(type_graph), intent(inout) :: graph
    integer         , intent(in   ) :: ni

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'sort_border_node_connectivity_nCeq2'
    integer                        :: nj, nk
    real(dp), dimension(2)         :: dj, dk
    real(dp)                       :: z

    ! Add routine to path
    call init_routine( routine_name)

    nj = graph%C( ni,1)
    nk = graph%C( ni,2)

    dj = graph%V( nj,:) - graph%V( ni,:)
    dk = graph%V( nk,:) - graph%V( ni,:)

    z = cross2( dj,dk)

    if (z < 0._dp) then
      ! They are sorted the wrong way; switch them
      graph%C( ni,1:2) = [nk,nj]
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine sort_border_node_connectivity_nCeq2

  subroutine sort_border_node_connectivity_nCgt2( graph, ni)

    ! In/output variables:
    type(type_graph), intent(inout) :: graph
    integer         , intent(in   ) :: ni

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'sort_border_node_connectivity_nCgt2'
    integer                        :: vj_first, vj_last, cj, vk
    logical                        :: first_and_last_are_connected

    ! Add routine to path
    call init_routine( routine_name)

    do while (.true.)

      vj_first = graph%C( ni, 1)
      vj_last  = graph%C( ni, graph%nC( ni))

      first_and_last_are_connected = .false.
      do cj = 1, graph%nC( vj_first)
        vk = graph%C( vj_first, cj)
        if (vk == vj_last) first_and_last_are_connected = .true.
      end do

      if (.not. first_and_last_are_connected) then
        ! They are now sorted
        exit
      end if

      ! Cycle by one element
      graph%C( ni, 1:graph%nC( ni)) = [graph%C( ni, graph%nC( ni)), graph%C( ni, 1:graph%nC( ni)-1)]

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine sort_border_node_connectivity_nCgt2

  subroutine calc_outward_unit_normal_vectors_at_border_triangles( graph, mesh, &
    mask_a_tot, mask_border_a_tot, mask_b_tot, mask_border_b_tot)

    ! In/output variables:
    type(type_graph),              intent(inout) :: graph
    type(type_mesh),               intent(in   ) :: mesh
    logical, dimension(mesh%nV),   intent(in   ) :: mask_a_tot, mask_border_a_tot
    logical, dimension(mesh%nTri), intent(in   ) :: mask_b_tot, mask_border_b_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_outward_unit_normal_vectors_at_border_triangles'
    integer                        :: ni, ti, n1, n2, vi, vj
    real(dp), dimension(2)         :: outward_normal_vector

    ! Add routine to path
    call init_routine( routine_name)

    do ni = 1, graph%n
      if (.not. graph%is_border( ni)) cycle

      ti = graph%ni2ti( ni)

      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1

        vi = mesh%Tri( ti, n1)
        vj = mesh%Tri( ti, n2)

        if (mask_border_a_tot( vi) .and. mask_border_a_tot( vj)) exit

      end do

      outward_normal_vector = [mesh%V( vj,2) - mesh%V( vi,2), mesh%V( vi,1) - mesh%V( vj,1)]

      graph%border_nhat( ni,:) = outward_normal_vector / norm2( outward_normal_vector)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_outward_unit_normal_vectors_at_border_triangles

end module create_graphs_from_masked_mesh