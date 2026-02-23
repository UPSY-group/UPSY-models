module graph_operators

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use graph_types, only: type_graph
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_empty_row_CSR_dist, &
    add_entry_CSR_dist, finalise_matrix_CSR_dist
  use shape_functions, only: calc_shape_functions_2D_reg_2nd_order, &
    calc_shape_functions_2D_stag_1st_order, calc_shape_functions_2D_reg_1st_order
  use mesh_types, only: type_mesh
  use tests_main

  implicit none

  private

  public :: calc_graph_matrix_operators_2nd_order, calc_graph_matrix_operators_1st_order, &
    calc_graph_a_to_graph_b_matrix_operators, calc_graph_b_to_graph_a_matrix_operators

contains

  subroutine calc_graph_matrix_operators_2nd_order( graph, M2_ddx, M2_ddy, M2_d2dx2, M2_d2dxdy, M2_d2dy2)
    !< Calculate 2nd-order accurate matrix operators on a graph

    ! In/output variables:
    type(type_graph),                intent(in   ) :: graph
    type(type_sparse_matrix_CSR_dp), intent(  out) :: M2_ddx, M2_ddy, M2_d2dx2, M2_d2dxdy, M2_d2dy2

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'calc_graph_matrix_operators_2nd_order'
    integer                             :: ncols, ncols_loc, nrows, nrows_loc
    integer                             :: nnz_per_row_est, nnz_est_proc
    integer                             :: ni
    real(dp)                            :: x, y
    integer                             :: nj
    integer                             :: n_neighbours_min
    integer                             :: n_neighbours_max
    integer,  dimension(graph%n)        :: map, stack, list
    integer                             :: stackN, listN
    integer                             :: i
    integer                             :: n_c
    integer,  dimension(:), allocatable :: i_c
    real(dp), dimension(:), allocatable :: x_c, y_c
    real(dp)                            :: Nfx_i, Nfy_i, Nfxx_i, Nfxy_i, Nfyy_i
    real(dp), dimension(:), allocatable :: Nfx_c, Nfy_c, Nfxx_c, Nfxy_c, Nfyy_c
    logical                             :: succeeded

    ! Add routine to path
    call init_routine( routine_name)

    n_neighbours_min = 5

    ! == Initialise the matrices using the native UPSY CSR-matrix format
    ! ==================================================================

    ! Matrix size
    ncols           = graph%n        ! from
    ncols_loc       = graph%n_loc
    nrows           = graph%n        ! to
    nrows_loc       = graph%n_loc
    nnz_per_row_est = graph%nC_mem+1
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( M2_ddx   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = graph%pai, pai_y = graph%pai)
    call allocate_matrix_CSR_dist( M2_ddy   , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = graph%pai, pai_y = graph%pai)
    call allocate_matrix_CSR_dist( M2_d2dx2 , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = graph%pai, pai_y = graph%pai)
    call allocate_matrix_CSR_dist( M2_d2dxdy, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = graph%pai, pai_y = graph%pai)
    call allocate_matrix_CSR_dist( M2_d2dy2 , nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = graph%pai, pai_y = graph%pai)

    ! Calculate shape functions and fill them into the matrices
    ! =========================================================

    ! allocate memory for map, stack, and shape functions
    n_neighbours_max = 32
    allocate( i_c(    n_neighbours_max))
    allocate( x_c(    n_neighbours_max))
    allocate( y_c(    n_neighbours_max))
    allocate( Nfx_c(  n_neighbours_max))
    allocate( Nfy_c(  n_neighbours_max))
    allocate( Nfxx_c( n_neighbours_max))
    allocate( Nfxy_c( n_neighbours_max))
    allocate( Nfyy_c( n_neighbours_max))

    map    = 0
    stack  = 0
    stackN = 0

    do ni = graph%ni1, graph%ni2

      ! Calculate shape functions for this node
      x = graph%V( ni,1)
      y = graph%V( ni,2)

      ! Initialise local neighbourhood with node ni
      map    = 0
      stackN = 1
      stack( 1) = ni
      map( ni) = 1
      listN = 0

      ! Calculate shape functions; if this fails, add more neighbours until it succeeds
      succeeded = .false.
      do while (.not. succeeded)

        call extend_group_single_iteration_graph( graph, map, stack, stackN, list, listN)

        ! Remove ni from local neighbourhood
        n_c = listN - 1
        i_c( 1:n_c) = list( 2: listN)

        if (n_c < n_neighbours_min) cycle

        ! List coordinates of neighbourhood nodes
        do i = 1, n_c
          nj = i_c( i)
          x_c( i) = graph%V( nj,1)
          y_c( i) = graph%V( nj,2)
        end do

        ! Calculate shape functions
        call calc_shape_functions_2D_reg_2nd_order( x, y, n_neighbours_max, n_c, x_c, y_c, &
          Nfx_i, Nfy_i, Nfxx_i, Nfxy_i, Nfyy_i, &
          Nfx_c, Nfy_c, Nfxx_c, Nfxy_c, Nfyy_c, succeeded)

      end do

      ! Fill them into the matrices

      ! Diagonal elements: shape functions for the home element
      call add_entry_CSR_dist( M2_ddx   , ni, ni, Nfx_i )
      call add_entry_CSR_dist( M2_ddy   , ni, ni, Nfy_i )
      call add_entry_CSR_dist( M2_d2dx2 , ni, ni, Nfxx_i)
      call add_entry_CSR_dist( M2_d2dxdy, ni, ni, Nfxy_i)
      call add_entry_CSR_dist( M2_d2dy2 , ni, ni, Nfyy_i)

      ! Off-diagonal elements: shape functions for the neighbours
      do i = 1, n_c
        nj = i_c( i)
        call add_entry_CSR_dist( M2_ddx   , ni, nj, Nfx_c(  i))
        call add_entry_CSR_dist( M2_ddy   , ni, nj, Nfy_c(  i))
        call add_entry_CSR_dist( M2_d2dx2 , ni, nj, Nfxx_c( i))
        call add_entry_CSR_dist( M2_d2dxdy, ni, nj, Nfxy_c( i))
        call add_entry_CSR_dist( M2_d2dy2 , ni, nj, Nfyy_c( i))
      end do

    end do

    ! Crop matrix memory
    call finalise_matrix_CSR_dist( M2_ddx   )
    call finalise_matrix_CSR_dist( M2_ddy   )
    call finalise_matrix_CSR_dist( M2_d2dx2 )
    call finalise_matrix_CSR_dist( M2_d2dxdy)
    call finalise_matrix_CSR_dist( M2_d2dy2 )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_graph_matrix_operators_2nd_order

  subroutine calc_graph_matrix_operators_1st_order( graph, M_ddx, M_ddy)
    !< Calculate 1st-order accurate matrix operators on a graph

    ! In/output variables:
    type(type_graph),                intent(in   ) :: graph
    type(type_sparse_matrix_CSR_dp), intent(  out) :: M_ddx, M_ddy

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'calc_graph_matrix_operators_1st_order'
    integer                             :: ncols, ncols_loc, nrows, nrows_loc
    integer                             :: nnz_per_row_est, nnz_est_proc
    integer                             :: ni
    real(dp)                            :: x, y
    integer                             :: nj
    integer                             :: n_neighbours_min
    integer                             :: n_neighbours_max
    integer,  dimension(graph%n)        :: map, stack, list
    integer                             :: stackN, listN
    integer                             :: i
    integer                             :: n_c
    integer,  dimension(:), allocatable :: i_c
    real(dp), dimension(:), allocatable :: x_c, y_c
    real(dp)                            :: Nfx_i, Nfy_i
    real(dp), dimension(:), allocatable :: Nfx_c, Nfy_c
    logical                             :: succeeded

    ! Add routine to path
    call init_routine( routine_name)

    n_neighbours_min = 2

    ! == Initialise the matrices using the native UPSY CSR-matrix format
    ! ==================================================================

    ! Matrix size
    ncols           = graph%n        ! from
    ncols_loc       = graph%n_loc
    nrows           = graph%n        ! to
    nrows_loc       = graph%n_loc
    nnz_per_row_est = graph%nC_mem+1
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( M_ddx, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = graph%pai, pai_y = graph%pai)
    call allocate_matrix_CSR_dist( M_ddy, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = graph%pai, pai_y = graph%pai)

    ! Calculate shape functions and fill them into the matrices
    ! =========================================================

    ! allocate memory for map, stack, and shape functions
    n_neighbours_max = 32
    allocate( i_c(   n_neighbours_max))
    allocate( x_c(   n_neighbours_max))
    allocate( y_c(   n_neighbours_max))
    allocate( Nfx_c( n_neighbours_max))
    allocate( Nfy_c( n_neighbours_max))

    map    = 0
    stack  = 0
    stackN = 0

    do ni = graph%ni1, graph%ni2

      ! Calculate shape functions for this node
      x = graph%V( ni,1)
      y = graph%V( ni,2)

      ! Initialise local neighbourhood with node ni
      map    = 0
      stackN = 1
      stack( 1) = ni
      map( ni) = 1
      listN = 0

      ! Calculate shape functions; if this fails, add more neighbours until it succeeds
      succeeded = .false.
      do while (.not. succeeded)

        call extend_group_single_iteration_graph( graph, map, stack, stackN, list, listN)

        ! Remove ni from local neighbourhood
        n_c = listN - 1
        i_c( 1:n_c) = list( 2: listN)

        if (n_c < n_neighbours_min) cycle

        ! List coordinates of neighbourhood nodes
        do i = 1, n_c
          nj = i_c( i)
          x_c( i) = graph%V( nj,1)
          y_c( i) = graph%V( nj,2)
        end do

        ! Calculate shape functions
        call calc_shape_functions_2D_reg_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, &
          Nfx_i, Nfy_i, &
          Nfx_c, Nfy_c, succeeded)

      end do

      ! Fill them into the matrices

      ! Diagonal elements: shape functions for the home element
      call add_entry_CSR_dist( M_ddx, ni, ni, Nfx_i)
      call add_entry_CSR_dist( M_ddy, ni, ni, Nfy_i)

      ! Off-diagonal elements: shape functions for the neighbours
      do i = 1, n_c
        nj = i_c( i)
        call add_entry_CSR_dist( M_ddx, ni, nj, Nfx_c( i))
        call add_entry_CSR_dist( M_ddy, ni, nj, Nfy_c( i))
      end do

    end do

    ! Crop matrix memory
    call finalise_matrix_CSR_dist( M_ddx)
    call finalise_matrix_CSR_dist( M_ddy)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_graph_matrix_operators_1st_order

  subroutine calc_graph_a_to_graph_b_matrix_operators( mesh, graph_a, graph_b, &
    M_map_a_b, M_ddx_a_b, M_ddy_a_b)
    !< Calculate 1st-order accurate matrix operators between two graphs

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_graph),                intent(in   ) :: graph_a, graph_b
    type(type_sparse_matrix_CSR_dp), intent(  out) :: M_map_a_b, M_ddx_a_b, M_ddy_a_b

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'calc_graph_a_to_graph_b_matrix_operators'
    integer                             :: ncols, ncols_loc, nrows, nrows_loc
    integer                             :: nnz_per_row_est, nnz_est_proc
    integer                             :: ni
    integer                             :: ti, via, vib, vic, nja, njb, njc
    real(dp)                            :: x, y
    integer                             :: nj
    integer                             :: n_neighbours_min
    integer                             :: n_neighbours_max
    integer                             :: i
    integer                             :: n_c
    integer,  dimension(:), allocatable :: i_c
    real(dp), dimension(:), allocatable :: x_c, y_c
    real(dp), dimension(:), allocatable :: Nf_c, Nfx_c, Nfy_c
    logical                             :: succeeded
    integer                             :: ei, vi, vj

    ! Add routine to path
    call init_routine( routine_name)

    n_neighbours_min = 3

    ! == Initialise the matrices using the native UPSY CSR-matrix format
    ! ==================================================================

    ! Matrix size
    ncols           = graph_a%n        ! from
    ncols_loc       = graph_a%n_loc
    nrows           = graph_b%n        ! to
    nrows_loc       = graph_b%n_loc
    nnz_per_row_est = graph_a%nC_mem+1
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( M_map_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = graph_a%pai, pai_y = graph_b%pai)
    call allocate_matrix_CSR_dist( M_ddx_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = graph_a%pai, pai_y = graph_b%pai)
    call allocate_matrix_CSR_dist( M_ddy_a_b, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = graph_a%pai, pai_y = graph_b%pai)

    ! Calculate shape functions and fill them into the matrices
    ! =========================================================

    ! allocate memory for map, stack, and shape functions
    n_neighbours_max = 32
    allocate( i_c(    n_neighbours_max))
    allocate( x_c(    n_neighbours_max))
    allocate( y_c(    n_neighbours_max))
    allocate( Nf_c(   n_neighbours_max))
    allocate( Nfx_c(  n_neighbours_max))
    allocate( Nfy_c(  n_neighbours_max))

    do ni = graph_b%ni1, graph_b%ni2

      ! Calculate shape functions at this graph node
      x = graph_b%V( ni,1)
      y = graph_b%V( ni,2)

      ! Set local neighbourhood to the vertices spanning triangle ti
      ti = graph_b%ni2ti( ni)
      via = mesh%Tri( ti,1)
      vib = mesh%Tri( ti,2)
      vic = mesh%Tri( ti,3)
      nja = graph_a%vi2ni( via)
      njb = graph_a%vi2ni( vib)
      njc = graph_a%vi2ni( vic)

      n_c = 3
      i_c( 1:3) = [nja, njb, njc]

      ! List coordinates of neighbourhood nodes
      do i = 1, n_c
        nj = i_c( i)
        x_c( i) = graph_a%V( nj,1)
        y_c( i) = graph_a%V( nj,2)
      end do

      ! Calculate shape functions
      call calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, &
        Nf_c, Nfx_c, Nfy_c, succeeded)

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      do i = 1, n_c
        nj = i_c( i)
        call add_entry_CSR_dist( M_map_a_b, ni, nj, Nf_c ( i))
        call add_entry_CSR_dist( M_ddx_a_b, ni, nj, Nfx_c( i))
        call add_entry_CSR_dist( M_ddy_a_b, ni, nj, Nfy_c( i))
      end do

    end do

    ! Crop matrix memory
    call finalise_matrix_CSR_dist( M_map_a_b)
    call finalise_matrix_CSR_dist( M_ddx_a_b)
    call finalise_matrix_CSR_dist( M_ddy_a_b)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_graph_a_to_graph_b_matrix_operators

  subroutine calc_graph_b_to_graph_a_matrix_operators( mesh, graph_b, graph_a, &
    M_map_b_a, M_ddx_b_a, M_ddy_b_a)
    !< Calculate 1st-order accurate matrix operators between two graphs

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_graph),                intent(in   ) :: graph_b, graph_a
    type(type_sparse_matrix_CSR_dp), intent(  out) :: M_map_b_a, M_ddx_b_a, M_ddy_b_a

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'calc_graph_b_to_graph_a_matrix_operators'
    integer                             :: ncols, ncols_loc, nrows, nrows_loc
    integer                             :: nnz_per_row_est, nnz_est_proc
    integer                             :: ni
    integer                             :: vi, iti, ti
    real(dp)                            :: x, y
    integer                             :: nj
    integer                             :: n_neighbours_min
    integer                             :: n_neighbours_max
    integer,  dimension(graph_b%n)      :: map, stack, list
    integer                             :: stackN, listN
    integer                             :: i
    integer                             :: n_c
    integer,  dimension(:), allocatable :: i_c
    real(dp), dimension(:), allocatable :: x_c, y_c
    real(dp), dimension(:), allocatable :: Nf_c, Nfx_c, Nfy_c
    logical                             :: succeeded

    ! Add routine to path
    call init_routine( routine_name)

    n_neighbours_min = 3

    ! == Initialise the matrices using the native UPSY CSR-matrix format
    ! ==================================================================

    ! Matrix size
    ncols           = graph_b%n        ! from
    ncols_loc       = graph_b%n_loc
    nrows           = graph_a%n        ! to
    nrows_loc       = graph_a%n_loc
    nnz_per_row_est = graph_a%nC_mem+1
    nnz_est_proc    = nrows_loc * nnz_per_row_est

    call allocate_matrix_CSR_dist( M_map_b_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = graph_a%pai, pai_y = graph_b%pai)
    call allocate_matrix_CSR_dist( M_ddx_b_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = graph_a%pai, pai_y = graph_b%pai)
    call allocate_matrix_CSR_dist( M_ddy_b_a, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc, &
      pai_x = graph_a%pai, pai_y = graph_b%pai)

    ! Calculate shape functions and fill them into the matrices
    ! =========================================================

    ! allocate memory for map, stack, and shape functions
    n_neighbours_max = 32
    allocate( i_c(    n_neighbours_max))
    allocate( x_c(    n_neighbours_max))
    allocate( y_c(    n_neighbours_max))
    allocate( Nf_c(   n_neighbours_max))
    allocate( Nfx_c(  n_neighbours_max))
    allocate( Nfy_c(  n_neighbours_max))

    map    = 0
    stack  = 0
    stackN = 0

    do ni = graph_a%ni1, graph_a%ni2

      ! Calculate shape functions at this graph node
      x = graph_a%V( ni,1)
      y = graph_a%V( ni,2)

      ! Initialise the local neighbourhood with the masked itriangles around this vertex
      map    = 0
      stackN = 0
      listN = 0

      vi = graph_a%ni2vi( ni)
      do iti = 1, mesh%niTri( vi)
        ti = mesh%iTri( vi,iti)
        nj = graph_b%ti2ni( ti)
        if (nj > 0) then
          stackN = stackN+1
          stack( stackN) = nj
          map( ni) = 1
        end if
      end do

      ! Calculate shape functions; if this fails, add more neighbours until it succeeds
      succeeded = .false.
      do while (.not. succeeded)

        call extend_group_single_iteration_graph( graph_b, map, stack, stackN, list, listN)

        n_c = listN
        i_c( 1:n_c) = list( 1: listN)

        if (n_c < n_neighbours_min) cycle

        ! List coordinates of neighbourhood nodes
        do i = 1, n_c
          nj = i_c( i)
          x_c( i) = graph_b%V( nj,1)
          y_c( i) = graph_b%V( nj,2)
        end do

        ! Calculate shape functions
        call calc_shape_functions_2D_stag_1st_order( x, y, n_neighbours_max, n_c, x_c, y_c, &
          Nf_c, Nfx_c, Nfy_c, succeeded)

      end do

      ! Fill them into the matrices

      ! Off-diagonal elements: shape functions for the neighbours
      do i = 1, n_c
        nj = i_c( i)
        call add_entry_CSR_dist( M_map_b_a, ni, nj, Nf_c ( i))
        call add_entry_CSR_dist( M_ddx_b_a, ni, nj, Nfx_c( i))
        call add_entry_CSR_dist( M_ddy_b_a, ni, nj, Nfy_c( i))
      end do

    end do

    ! Crop matrix memory
    call finalise_matrix_CSR_dist( M_map_b_a)
    call finalise_matrix_CSR_dist( M_ddx_b_a)
    call finalise_matrix_CSR_dist( M_ddy_b_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_graph_b_to_graph_a_matrix_operators

  subroutine extend_group_single_iteration_graph( graph, map, stack, stackN, list, listN)

    ! In/output variables:
    type(type_graph),      intent(in   ) :: graph
    integer, dimension(:), intent(inout) :: map, stack, list
    integer,               intent(inout) :: stackN, listN

    ! Local variables:
    integer, dimension(size(stack)) :: stack2
    integer                         :: stack2N
    integer                         :: i, ni, ci, nj

    stack2N = 0

    do i = 1, stackN

      ni = stack( i)

      map( ni) = 2

      listN = listN + 1
      list( listN) = ni

      do ci = 1, graph%nC( ni)
        nj = graph%C( ni,ci)
        if (map( nj) == 0) then
          map( nj) = 1
          stack2N = stack2N + 1
          stack2( stack2N) = nj
        end if
      end do

    end do

    stack( 1:stack2N) = stack2( 1: stack2N)
    stackN = stack2N

  end subroutine extend_group_single_iteration_graph

end module graph_operators