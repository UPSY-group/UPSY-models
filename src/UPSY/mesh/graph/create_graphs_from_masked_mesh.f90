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

    ! Create vertex-to-regular-node mapping
    graph%n = 0
    do vi = 1, mesh%nV
      if (mask_a_tot( vi)) then
        graph%n = graph%n + 1
        graph%vi2ni( vi     ) = graph%n
        graph%ni2vi( graph%n) = vi
      end if
    end do

    ! Construct regular nodes and their interconnectivity
    do ni = 1, graph%n

      ! This regular graph node corresponds to a masked vertex
      vi = graph%ni2vi( ni)

      graph%V         ( ni,:) = mesh%V( vi,:)
      graph%is_ghost  ( ni  ) = .false.
      graph%ghost_nhat( ni,:) = NaN

      ! Add connections to other regular nodes
      graph%nC( ni) = 0
      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        if (mask_a_tot( vj)) then
          ! This connection points to another regular node
          nj = graph%vi2ni( vj)
          graph%nC( ni) = graph%nC( ni) + 1
          graph%C( ni, graph%nC( ni)) = nj
        end if
      end do

    end do

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
    logical, dimension(1:mesh%nTri) :: mask_b_tot
    integer                         :: ti, n, tj, nj, ei
    integer                         :: ni_reg, ni_ghost
    integer                         :: ni, vi, vj, ci, ej, nj1, nj2
    real(dp), dimension(2)          :: r, s, outward_normal_vector

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate b-grid mask
    call calc_mask_b( mesh, mask_a, mask_b_tot)

    ! Allocate graph and set metadata
    call allocate_graph_primary( graph, mesh%nTri*3, mesh%nV, mesh%nTri, mesh%nE, 3)

    graph%parent_mesh_name = trim( mesh%name) // '_triangles'
    graph%xmin             = mesh%xmin
    graph%xmax             = mesh%xmax
    graph%ymin             = mesh%ymin
    graph%ymax             = mesh%ymax

    ! Create triangle-to-regular-node mapping
    graph%n = 0
    do ti = 1, mesh%nTri
      if (mask_b_tot( ti)) then
        graph%n = graph%n + 1
        graph%ti2ni( ti     ) = graph%n
        graph%ni2ti( graph%n) = ti
      end if
    end do

    ! Construct regular nodes and their interconnectivity
    do ni = 1, graph%n

      ! This regular graph node corresponds to a masked triangle
      ti = graph%ni2ti( ni)

      graph%V         ( ni,:) = mesh%TriGC( ti,:)
      graph%is_ghost  ( ni  ) = .false.
      graph%ghost_nhat( ni,:) = NaN

      ! Add connections to other regular nodes
      graph%nC( ni) = 0
      do n = 1, 3
        tj = mesh%TriC( ti,n)
        nj = graph%ti2ni( tj)
        if (nj > 0) then
          ! This connection points to another regular node
          graph%nC( ni) = graph%nC( ni) + 1
          graph%C( ni, graph%nC( ni)) = nj
        end if
      end do

    end do

    ! Add ghost nodes
    do ti = 1, mesh%nTri
      if (mask_b_tot( ti)) then

        ! Regular node
        ni_reg = graph%ti2ni( ti)

        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (.not. mask_b_tot( tj)) then
            ! This connection points from a mask to an unmasked triangle;
            ! add a ghost node here

            graph%n = graph%n + 1
            ni_ghost = graph%n

            ei = mesh%TriE( ti,n)
            graph%ei2ni( ei      ) = ni_ghost
            graph%ni2ei( ni_ghost) = ei

            graph%V ( ni_ghost,:) = mesh%E( ei,:)
            graph%nC( ni_ghost  ) = 1
            graph%C ( ni_ghost,1) = ni_reg

            ! Add reverse connection
            graph%nC( ni_reg) = graph%nC( ni_reg) + 1
            graph%C( ni_reg, graph%nC( ni_reg)) = ni_ghost

            graph%is_ghost( ni_ghost) = .true.

          end if
        end do

      end if
    end do

    ! Connect ghost nodes to each other along the front
    do ni = 1, graph%n
      if (graph%is_ghost( ni)) then

        ei = graph%ni2ei( ni)
        if (mask_b_tot( mesh%ETri( ei,1))) then
          ti = mesh%ETri( ei,1)
        elseif (mask_b_tot( mesh%ETri( ei,2))) then
          ti = mesh%ETri( ei,2)
        else
          call crash('whaa')
        end if
        vi = mesh%EV( ei,1)
        vj = mesh%EV( ei,2)

        nj1 = 0
        do ci = 1, mesh%nC( vi)
          ej = mesh%VE( vi,ci)
          nj = graph%ei2ni( ej)
          if (nj > 0 .and. nj /= ni) then
            nj1 = nj
            exit
          end if
        end do

        nj2 = 0
        do ci = 1, mesh%nC( vj)
          ej = mesh%VE( vj,ci)
          nj = graph%ei2ni( ej)
          if (nj > 0 .and. nj /= ni) then
            nj2 = nj
            exit
          end if
        end do

        if (nj1 == 0 .or. nj2 == 0) call crash('couldnt find both nearby ghost nodes')

        ! Add connections, maintaining counter-clockwise ordering
        graph%nC( ni) = 3
        graph%C( ni,2) = graph%C( ni,1)
        if (mesh%ETri( ei,1) == ti) then
          graph%C( ni, 1) = nj2
          graph%C( ni, 3) = nj1
        elseif (mesh%ETri( ei,2) == ti) then
          graph%C( ni, 1) = nj1
          graph%C( ni, 3) = nj2
        else
          call crash('inconsistency in mesh%ETri')
        end if

      end if
    end do

    ! Calculate outward unit normal vectors for ghost nodes
    do ni = 1, graph%n
      if (graph%is_ghost( ni)) then

        nj1 = graph%C( ni,1)
        nj2 = graph%C( ni,3)

        r = graph%V( nj1,:)
        s = graph%V( nj2,:)

        outward_normal_vector = [r(2) - s(2), s(1) - r(1)]
        graph%ghost_nhat( ni,:) = outward_normal_vector / norm2( outward_normal_vector)

      end if
    end do

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

  subroutine calc_mask_b( mesh, mask_a, mask_b_tot)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_a
    logical, dimension(mesh%nTri),         intent(  out) :: mask_b_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_mask_b'
    logical, dimension(mesh%nV  )  :: mask_a_tot
    integer                        :: ti, n_a, n, vi

    ! Add routine to path
    call init_routine( routine_name)

    call gather_to_all( mask_a, mask_a_tot)

    ! b mask
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

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_mask_b

end module create_graphs_from_masked_mesh