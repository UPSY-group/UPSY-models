module graph_pair_creation

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph_pair
  use ice_model_types, only: type_ice_model
  use create_graphs_from_masked_mesh, only: create_graph_from_masked_mesh_a, &
    create_graph_from_masked_mesh_b
  use graph_operators, only: calc_graph_matrix_operators_2nd_order, &
    calc_graph_matrix_operators_1st_order, &
    calc_graph_a_to_graph_b_matrix_operators, calc_graph_b_to_graph_a_matrix_operators
  use graph_memory, only: deallocate_graph
  use CSR_matrix_basics, only: deallocate_matrix_CSR_dist

  implicit none

  private

  public :: create_ice_only_graph_pair, deallocate_graph_pair

contains

  subroutine create_ice_only_graph_pair( mesh, ice, graphs)

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    type(type_ice_model),  intent(in   ) :: ice
    type(type_graph_pair), intent(  out) :: graphs

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'create_ice_only_graph_pair'
    logical, dimension(mesh%vi1:mesh%vi2) :: mask_ice_a
    integer                               :: vi, ni, uv, n

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the ice mask
    do vi = mesh%vi1, mesh%vi2
      mask_ice_a( vi) = ice%Hi( vi) > 0._dp
    end do

    ! Create graphs from the masked vertices and triangles
    call create_graph_from_masked_mesh_a( mesh, mask_ice_a, mesh%nz, graphs%graph_a)
    call create_graph_from_masked_mesh_b( mesh, mask_ice_a, mesh%nz, graphs%graph_b)

    ! Calculate matrix operators
    call calc_graph_matrix_operators_2nd_order( graphs%graph_b, &
      graphs%M2_ddx_b_b, &
      graphs%M2_ddy_b_b, &
      graphs%M2_d2dx2_b_b, &
      graphs%M2_d2dxdy_b_b, &
      graphs%M2_d2dy2_b_b)

    call calc_graph_matrix_operators_1st_order( graphs%graph_b, &
      graphs%M_ddx_b_b, graphs%M_ddy_b_b)

    call calc_graph_a_to_graph_b_matrix_operators( mesh, graphs%graph_a, graphs%graph_b, &
      graphs%M_map_a_b, graphs%M_ddx_a_b, graphs%M_ddy_a_b)

    call calc_graph_b_to_graph_a_matrix_operators( mesh, graphs%graph_b, graphs%graph_a, &
      graphs%M_map_b_a, graphs%M_ddx_b_a, graphs%M_ddy_b_a)

    ! Calculate b-graph-node-to-matrix-row translation tables
    allocate( graphs%biuv2n( graphs%graph_b%n   ,  2), source = 0)
    allocate( graphs%n2biuv( graphs%graph_b%n * 2, 2), source = 0)

    n = 0
    do ni = 1, graphs%graph_b%n
      do uv = 1, 2
        n = n+1
        graphs%biuv2n( ni,uv) = n
        graphs%n2biuv( n,1) = ni
        graphs%n2biuv( n,2) = uv
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 8)

  end subroutine create_ice_only_graph_pair

  subroutine deallocate_graph_pair( graphs)

    ! In/output variables:
    type(type_graph_pair), intent(inout) :: graphs

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'deallocate_graph_pair'

    ! Add routine to path
    call init_routine( routine_name)

    ! Create graphs from the masked vertices and triangles
    call deallocate_graph( graphs%graph_a)
    call deallocate_graph( graphs%graph_b)

    ! Calculate matrix operators
    call deallocate_matrix_CSR_dist( graphs%M2_ddx_b_b   )
    call deallocate_matrix_CSR_dist( graphs%M2_ddy_b_b   )
    call deallocate_matrix_CSR_dist( graphs%M2_d2dx2_b_b )
    call deallocate_matrix_CSR_dist( graphs%M2_d2dxdy_b_b)
    call deallocate_matrix_CSR_dist( graphs%M2_d2dy2_b_b )

    call deallocate_matrix_CSR_dist( graphs%M_ddx_b_b)
    call deallocate_matrix_CSR_dist( graphs%M_ddy_b_b)

    call deallocate_matrix_CSR_dist( graphs%M_map_a_b)
    call deallocate_matrix_CSR_dist( graphs%M_ddx_a_b)
    call deallocate_matrix_CSR_dist( graphs%M_ddy_a_b)

    call deallocate_matrix_CSR_dist( graphs%M_map_b_a)
    call deallocate_matrix_CSR_dist( graphs%M_ddx_b_a)
    call deallocate_matrix_CSR_dist( graphs%M_ddy_b_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_graph_pair

end module graph_pair_creation
