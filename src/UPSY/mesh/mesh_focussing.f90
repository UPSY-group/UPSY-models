module mesh_focussing

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use polyline_types, only: type_polyline
  use remapping_types, only: type_single_row_mapping_matrices
  use line_tracing_Voronoi, only: trace_line_Vor
  use delete_vertices, only: delete_vertex
  use split_triangles, only: split_triangle
  use mesh_memory, only: duplicate_mesh_primary, extend_mesh_primary, crop_mesh_primary
  use mesh_edges, only: construct_mesh_edges
  use mesh_secondary, only: calc_all_secondary_mesh_data

  implicit none

  private

  public :: delete_vertices_along_polyline, focus_mesh_on_polyline

contains

  subroutine focus_mesh_on_polyline( mesh, ll, mesh_focused, &
    vi_new2vi_old, vi_old2vi_new, vi_ll2vi, vi2vi_ll)
    !< Focus a mesh on a polyline

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_polyline),                intent(in   ) :: ll
    type(type_mesh),                    intent(  out) :: mesh_focused
    integer, dimension(:), allocatable, intent(  out) :: vi_new2vi_old, vi_old2vi_new
    integer, dimension(:), allocatable, intent(  out) :: vi_ll2vi, vi2vi_ll

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'focus_mesh_on_polyline'
    integer, dimension(:), allocatable :: vi_new2vi_old_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Duplicate the mesh
    call duplicate_mesh_primary( mesh, mesh_focused)
    call construct_mesh_edges( mesh_focused)

    ! Focus the mesh
    call delete_vertices_along_polyline( mesh_focused, ll, vi_new2vi_old, vi_old2vi_new)
    call add_polyline_vertices_to_mesh ( mesh_focused, ll, vi_ll2vi, vi2vi_ll)
    call calc_all_secondary_mesh_data( mesh_focused, mesh%lambda_M, mesh%phi_M, mesh%beta_stereo)

    ! Update translation table
    vi_new2vi_old_temp = vi_new2vi_old
    deallocate( vi_new2vi_old)
    allocate( vi_new2vi_old( mesh%nV), source = 0)
    vi_new2vi_old( 1: size( vi_new2vi_old_temp)) = vi_new2vi_old_temp
    deallocate( vi_new2vi_old_temp)

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 12)

  end subroutine focus_mesh_on_polyline

  subroutine delete_vertices_along_polyline( mesh, ll, vi_new2vi_old, vi_old2vi_new)
    !< Delete all vertices whose Voronoi cells are crossed by the polyline

    ! In/output variables:
    type(type_mesh),                    intent(inout) :: mesh
    type(type_polyline),                intent(in   ) :: ll
    integer, dimension(:), allocatable, intent(  out) :: vi_new2vi_old, vi_old2vi_new

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'delete_vertices_along_polyline'
    integer, dimension(:), allocatable :: vi_list
    integer                            :: vi, i, vi_kill, ii
    integer, dimension(:), allocatable :: vi_new2vi_prev, vi_prev2vi_new

    ! Add routine to path
    call init_routine( routine_name)

    call list_vertices_crossed_by_polyline( mesh, ll, vi_list)

    allocate( vi_new2vi_old( mesh%nV))
    allocate( vi_old2vi_new( mesh%nV))
    do vi = 1, mesh%nV
      vi_new2vi_old( vi) = vi
      vi_old2vi_new( vi) = vi
    end do

    do i = 1, size( vi_list)

      ! Delete vertex
      vi_kill = vi_list( i)
      call delete_vertex( mesh, vi_kill, vi_new2vi_prev, vi_prev2vi_new)

      ! Update indices of remaining to-be-deleted vertices
      do ii = i+1, size( vi_list)
        vi_list( ii) = vi_prev2vi_new( vi_list( ii))
      end do

      ! Update vertex translation tables
      call delete_vertices_along_polyline_update_vertex_translation_tables( mesh, &
        vi_new2vi_prev, vi_prev2vi_new, vi_new2vi_old, vi_old2vi_new)

      ! Clean up after yourself (not sure if memory is automatically deallocated when arrays are allocated in another procedure...)
      deallocate( vi_new2vi_prev, vi_prev2vi_new)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertices_along_polyline

  subroutine list_vertices_crossed_by_polyline( mesh, ll, vi_list)
    !< List all vertices whose Voronoi cells are crossed by the polyline

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_polyline),                intent(in   ) :: ll
    integer, dimension(:), allocatable, intent(  out) :: vi_list

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'delete_vertices_along_polyline'
    integer                                :: vi_hint
    integer                                :: i_end, i1, i2
    real(dp), dimension(2)                 :: p,q
    type(type_single_row_mapping_matrices) :: single_row

    ! Add routine to path
    call init_routine( routine_name)

    if (ll%is_closed) then
      ! Include the line segment from the last to the first vertex
      i_end = ll%n
    else
      ! Don't
      i_end = ll%n-1
    end if

    ! Allocate memory for single row results
    single_row%n_max = mesh%nV
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    vi_hint = 1
    do i1 = 1, i_end
      i2 = i1 + 1
      if (i2 == ll%n+1) i2 = 1

      p = ll%p( i1,:)
      q = ll%p( i2,:)

      call trace_line_Vor( mesh, p, q, single_row, .true., vi_hint)

    end do

    vi_list = single_row%index_left( 1: single_row%n)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine list_vertices_crossed_by_polyline

  subroutine delete_vertices_along_polyline_update_vertex_translation_tables( mesh, &
    vi_new2vi_prev, vi_prev2vi_new, vi_new2vi_old, vi_old2vi_new)

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    integer, dimension(:), allocatable, intent(inout) :: vi_new2vi_prev, vi_prev2vi_new, vi_new2vi_old, vi_old2vi_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertices_along_polyline'
    integer, dimension(:), allocatable :: vi_old2vi_prev, vi_prev2vi_old
    integer :: vi_new, vi_prev, vi_old

    ! Add routine to path
    call init_routine( routine_name)

    ! Save old lists
    vi_old2vi_prev = vi_old2vi_new
    vi_prev2vi_old = vi_new2vi_old

    deallocate( vi_old2vi_new, vi_new2vi_old)
    allocate( vi_old2vi_new( size( vi_old2vi_prev)), source = 0)
    allocate( vi_new2vi_old( mesh%nV), source = 0)

    ! Update lists
    vi_old2vi_new = 0
    do vi_new = 1, mesh%nV
      vi_prev = vi_new2vi_prev( vi_new)
      vi_old  = vi_prev2vi_old( vi_prev)
      vi_old2vi_new( vi_old) = vi_new
      vi_new2vi_old( vi_new) = vi_old
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertices_along_polyline_update_vertex_translation_tables

  subroutine add_polyline_vertices_to_mesh( mesh, ll, vi_ll2vi, vi2vi_ll)
    !< Add the vertices of the polyline to the mesh

    ! In/output variables:
    type(type_mesh),                    intent(inout) :: mesh
    type(type_polyline),                intent(in   ) :: ll
    integer, dimension(:), allocatable, intent(  out) :: vi_ll2vi, vi2vi_ll

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_polyline_vertices_to_mesh'
    integer                        :: i, ti_in, vi
    real(dp), dimension(2)         :: p

    ! Add routine to path
    call init_routine( routine_name)

    ! Make sure there's enough memory for the new vertices
    call extend_mesh_primary( mesh, mesh%nV + ll%n, mesh%nTri + ll%n*3)

    ti_in = 1
    do i = 1, ll%n
      p = ll%p( i,:)
      call split_triangle( mesh, ti_in, p)
    end do

    call crop_mesh_primary( mesh)

    ! Translation tables
    allocate( vi_ll2vi( ll%n   ), source = 0)
    allocate( vi2vi_ll( mesh%nV), source = 0)

    do i = 1, ll%n
      vi = mesh%nV - ll%n + i
      vi_ll2vi( i ) = vi
      vi2vi_ll( vi) = i
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_polyline_vertices_to_mesh

end module mesh_focussing
