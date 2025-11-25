module ut_mesh_delete_vertices

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mpi_basic, only: par
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary, crop_mesh_primary, deallocate_mesh
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use parameters, only: pi
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_edges, only: construct_mesh_edges
  use delete_vertices, only: delete_vertex
  use polyline_types, only: type_polyline
  use mesh_memory, only: duplicate_mesh_primary
  use mesh_focussing, only: delete_vertices_along_polyline, focus_mesh_on_polyline

  use netcdf_io_main
  use control_resources_and_error_messaging, only: insert_val_into_string_dp

  implicit none

  private

  public :: test_delete_vertices

contains

  subroutine test_delete_vertices( test_name_parent)
    ! Test the vertex deletion subroutines

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_delete_vertices'
    character(len=1024), parameter     :: test_name_local = 'delete_vertices'
    character(len=1024)                :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_delete_single_vertex          ( test_name)
    call test_delete_vertices_along_polyline( test_name)
    call test_focus_mesh_on_polyline        ( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_delete_vertices

  subroutine test_delete_single_vertex( test_name_parent)
    ! Test the delete_vertex subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_delete_vertex'
    character(len=1024), parameter     :: test_name_local = 'delete_vertex'
    character(len=1024)                :: test_name
    real(dp), parameter                :: xmin = 0._dp
    real(dp), parameter                :: xmax = 1._dp
    real(dp), parameter                :: ymin = 0._dp
    real(dp), parameter                :: ymax = 1._dp
    real(dp)                           :: alpha_min
    real(dp)                           :: res_max
    type(type_mesh)                    :: mesh, mesh2
    integer                            :: vi, ci, vj, vi_kill
    logical                            :: is_border_vertex, has_border_neighbours
    integer, dimension(:), allocatable :: vi_new2vi_old, vi_old2vi_new
    logical                            :: verified_vi_new2vi_old, verified_vi_old2vi_new
    integer                            :: vi_new, vi_old

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call allocate_mesh_primary( mesh, 'test_mesh', 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 20._dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call crop_mesh_primary( mesh)

    ! Find a nice vertex to kill
    ! (must not be a border vertex, nor be adjacent to any border vertices)
    vi_kill = 0
    do vi = 1, mesh%nV
      is_border_vertex = mesh%VBI( vi) > 0
      has_border_neighbours = .false.
      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        if (mesh%VBI( vj) > 0) then
          has_border_neighbours = .true.
          exit
        end if
      end do
      if (.not. is_border_vertex .and. .not. has_border_neighbours) then
        vi_kill = vi
        exit
      end if
    end do

    call duplicate_mesh_primary( mesh, mesh2)
    call delete_vertex( mesh2, vi_kill, vi_new2vi_old, vi_old2vi_new)

    verified_vi_new2vi_old = .true.
    do vi_new = 1, mesh2%nV
      vi_old = vi_new2vi_old( vi_new)
      verified_vi_new2vi_old = verified_vi_new2vi_old .and. &
        norm2( mesh2%V( vi_new,:) - mesh%V( vi_old,:)) < mesh%tol_dist
    end do

    verified_vi_old2vi_new = .true.
    do vi_old = 1, mesh%nV
      vi_new = vi_old2vi_new( vi_old)
      if (vi_new > 0) then
        verified_vi_old2vi_new = verified_vi_old2vi_new .and. &
          norm2( mesh2%V( vi_new,:) - mesh%V( vi_old,:)) < mesh%tol_dist
      end if
    end do

    call unit_test( test_mesh_is_self_consistent( mesh2), trim( test_name) // '/mesh_self_consistency')
    call unit_test( mesh2%nV   == mesh%nV   - 1, trim( test_name) // '/nV')
    call unit_test( mesh2%nTri == mesh%nTri - 2, trim( test_name) // '/nTri')
    call unit_test( verified_vi_new2vi_old, trim( test_name) // '/vi_new2vi_old')
    call unit_test( verified_vi_old2vi_new, trim( test_name) // '/vi_old2vi_new')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_delete_single_vertex

  subroutine test_delete_vertices_along_polyline( test_name_parent)
    ! Test the delete_vertices_along_polyline subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_delete_vertices_along_polyline'
    character(len=1024), parameter     :: test_name_local = 'delete_vertices_along_polyline'
    character(len=1024)                :: test_name
    real(dp), parameter                :: xmin = -1._dp
    real(dp), parameter                :: xmax =  1._dp
    real(dp), parameter                :: ymin = -1._dp
    real(dp), parameter                :: ymax =  1._dp
    real(dp)                           :: alpha_min
    real(dp)                           :: res_max
    type(type_mesh)                    :: mesh, mesh2
    type(type_polyline)                :: ll
    integer                            :: i
    real(dp)                           :: r, theta, x, y
    integer, dimension(:), allocatable :: vi_new2vi_old, vi_old2vi_new
    logical                            :: verified_vi_new2vi_old, verified_vi_old2vi_new
    integer                            :: vi_new, vi_old

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call allocate_mesh_primary( mesh, 'test_mesh', 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 50._dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call crop_mesh_primary( mesh)

    ! Define a nice circular polyline
    ll%is_closed = .true.
    ll%n         = 100
    allocate( ll%p( ll%n,2), source = 0._dp)

    r = mesh%xmax * 0.6_dp
    do i = 1, ll%n
      theta = real( i,dp) / real( ll%n,dp) * 2._dp * pi
      x = r * cos( theta)
      y = r * sin( theta)
      ll%p( i,:) = [x,y]
    end do

    ! Duplicate the mesh
    call duplicate_mesh_primary( mesh, mesh2)

    ! Delete vertices along the polyline
    call construct_mesh_edges( mesh2)
    call delete_vertices_along_polyline( mesh2, ll, vi_new2vi_old, vi_old2vi_new)

    verified_vi_new2vi_old = .true.
    do vi_new = 1, mesh2%nV
      vi_old = vi_new2vi_old( vi_new)
      verified_vi_new2vi_old = verified_vi_new2vi_old .and. &
        norm2( mesh2%V( vi_new,:) - mesh%V( vi_old,:)) < mesh%tol_dist
    end do

    verified_vi_old2vi_new = .true.
    do vi_old = 1, mesh%nV
      vi_new = vi_old2vi_new( vi_old)
      if (vi_new > 0) then
        verified_vi_old2vi_new = verified_vi_old2vi_new .and. &
          norm2( mesh2%V( vi_new,:) - mesh%V( vi_old,:)) < mesh%tol_dist
      end if
    end do

    call unit_test( test_mesh_is_self_consistent( mesh2), trim( test_name) // '/mesh_self_consistency')
    call unit_test( mesh2%nV < mesh%nV, trim( test_name) // '/nV')
    call unit_test( verified_vi_new2vi_old, trim( test_name) // '/vi_new2vi_old')
    call unit_test( verified_vi_old2vi_new, trim( test_name) // '/vi_old2vi_new')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_delete_vertices_along_polyline

  subroutine test_focus_mesh_on_polyline( test_name_parent)
    ! Test the focus_mesh_on_polyline subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_focus_mesh_on_polyline'
    character(len=1024), parameter     :: test_name_local = 'focus_mesh_on_polyline'
    character(len=1024)                :: test_name
    real(dp), parameter                :: xmin = -1._dp
    real(dp), parameter                :: xmax =  1._dp
    real(dp), parameter                :: ymin = -1._dp
    real(dp), parameter                :: ymax =  1._dp
    real(dp)                           :: alpha_min
    real(dp)                           :: res_max
    type(type_mesh)                    :: mesh, mesh_focused
    type(type_polyline)                :: ll
    integer                            :: i, vi
    real(dp)                           :: r, theta, x, y
    integer, dimension(:), allocatable :: vi_new2vi_old, vi_old2vi_new
    integer, dimension(:), allocatable :: vi_ll2vi, vi2vi_ll
    logical                            :: verified_vi_new2vi_old, verified_vi_old2vi_new
    logical                            :: verified_vi_ll2vi, verified_vi2vi_ll
    integer                            :: vi_new, vi_old

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call allocate_mesh_primary( mesh, 'test_mesh', 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 50._dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call crop_mesh_primary( mesh)

    ! Define a nice circular polyline
    ll%is_closed = .true.
    ll%n         = 100
    allocate( ll%p( ll%n,2), source = 0._dp)

    r = mesh%xmax * 0.6_dp
    do i = 1, ll%n
      theta = real( i,dp) / real( ll%n,dp) * 2._dp * pi
      x = r * cos( theta)
      y = r * sin( theta)
      ll%p( i,:) = [x,y]
    end do

    ! Focus mesh on the polyline
    call focus_mesh_on_polyline( mesh, ll, mesh_focused, &
      vi_new2vi_old, vi_old2vi_new, vi_ll2vi, vi2vi_ll)

    verified_vi_new2vi_old = .true.
    do vi_new = 1, mesh_focused%nV
      vi_old = vi_new2vi_old( vi_new)
      if (vi_old > 0) then
        verified_vi_new2vi_old = verified_vi_new2vi_old .and. &
          norm2( mesh_focused%V( vi_new,:) - mesh%V( vi_old,:)) < mesh%tol_dist
      end if
    end do

    verified_vi_old2vi_new = .true.
    do vi_old = 1, mesh%nV
      vi_new = vi_old2vi_new( vi_old)
      if (vi_new > 0) then
        verified_vi_old2vi_new = verified_vi_old2vi_new .and. &
          norm2( mesh_focused%V( vi_new,:) - mesh%V( vi_old,:)) < mesh%tol_dist
      end if
    end do

    verified_vi2vi_ll = .true.
    do vi = 1, mesh_focused%nV
      i = vi2vi_ll( vi)
      if (i > 0) then
        verified_vi2vi_ll = verified_vi2vi_ll .and. &
          norm2( mesh_focused%V( vi,:) - ll%p( i,:)) < mesh%tol_dist
      end if
    end do

    verified_vi_ll2vi = .true.
    do i = 1, ll%n
      vi = vi_ll2vi( i)
      verified_vi_ll2vi = verified_vi_ll2vi .and. &
        norm2( mesh_focused%V( vi,:) - ll%p( i,:)) < mesh%tol_dist
    end do

    call unit_test( test_mesh_is_self_consistent( mesh_focused), trim( test_name) // '/mesh_self_consistency')
    call unit_test( verified_vi_new2vi_old, trim( test_name) // '/vi_new2vi_old')
    call unit_test( verified_vi_old2vi_new, trim( test_name) // '/vi_old2vi_new')
    call unit_test( verified_vi2vi_ll, trim( test_name) // '/vi2vi_ll')
    call unit_test( verified_vi_ll2vi, trim( test_name) // '/vi_ll2vi')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_focus_mesh_on_polyline

end module ut_mesh_delete_vertices