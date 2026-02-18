module ct_discretisation_mapping_derivatives_graph

  ! Test the mesh matrix operators for mapping and derivatives

  use mpi_f08, only: MPI_COMM_WORLD, MPI_BCAST, MPI_CHAR
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use parameters
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use string_module, only: strrep
  use mpi_basic, only: par, sync
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph
  use create_graphs_from_masked_mesh, only: create_graph_from_masked_mesh_a, create_graph_from_masked_mesh_b
  use netcdf_io_main
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use graph_operators, only: calc_graph_matrix_operators_2nd_order, &
    calc_graph_a_to_graph_b_matrix_operators, calc_graph_b_to_graph_a_matrix_operators
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_1D_wrapper
  use ct_discretisation_mapping_derivatives, only: test_function, test_function_linear, &
    test_function_quadratic, test_function_periodic
  use assertions_basic
  use tests_main
  use mpi_f08, only: MPI_WIN
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared
  use mesh_graph_mapping, only: map_graph_to_mesh_vertices, map_graph_to_mesh_triangles

  implicit none

  private

  public :: run_all_map_deriv_tests_graph

contains

  !> Run all mapping/derivative tests.
  subroutine run_all_map_deriv_tests_graph( foldername_discretisation, test_mesh_filenames)

    ! In/output variables:
    character(len=*),               intent(in) :: foldername_discretisation
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_all_map_deriv_tests_graph'
    character(len=1024)            :: foldername_map_deriv
    character(len=1024)            :: test_mesh_filename
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '    Running graph mapping/derivative tests...'
    if (par%primary) write(0,*) ''

    call create_graph_map_deriv_tests_output_folder( foldername_discretisation, foldername_map_deriv)

    do i = 1, size(test_mesh_filenames)
      test_mesh_filename = trim(test_mesh_filenames( i))
      call run_all_map_deriv_tests_on_graphs( foldername_map_deriv, test_mesh_filename)
    end do

    if (par%primary) write(0,*) ''

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_map_deriv_tests_graph

  !> Create the output folder for the mappoing/derivative tests
  subroutine create_graph_map_deriv_tests_output_folder( foldername_discretisation, foldername_map_deriv)

    ! In/output variables:
    character(len=*), intent(in)  :: foldername_discretisation
    character(len=*), intent(out) :: foldername_map_deriv

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_graph_map_deriv_tests_output_folder'
    logical                        :: ex
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    foldername_map_deriv = trim(foldername_discretisation) // '/mapping_derivatives_graph'

    if (par%primary) then

      ! Remove existing folder if necessary
      inquire( file = trim( foldername_map_deriv) // '/.', exist = ex)
      if (ex) then
        call system('rm -rf ' // trim( foldername_map_deriv))
      end if

      ! Create the directory
      call system('mkdir ' // trim( foldername_map_deriv))

    end if
    call MPI_BCAST( foldername_map_deriv, len(foldername_map_deriv), MPI_CHAR, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_graph_map_deriv_tests_output_folder

  !> Run all the mapping/derivative tests on a particular graph
  subroutine run_all_map_deriv_tests_on_graphs( foldername_discretisation, test_mesh_filename)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_discretisation
    character(len=*), intent(in) :: test_mesh_filename

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'run_all_map_deriv_tests_on_graphs'
    type(type_mesh)                    :: mesh
    integer                            :: ncid
    logical, dimension(:), allocatable :: mask_a
    real(dp)                           :: xmid, ymid, wx, wy, wmin
    integer                            :: nz
    integer                            :: vi
    type(type_graph)                   :: graph_a, graph_b
    procedure(test_function), pointer  :: test_function_ptr
    character(len=1024)                :: function_name

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '      Running mapping/derivative tests on graphs from mesh ', &
      UPSY%stru%colour_string(trim(test_mesh_filename( index( test_mesh_filename,'/',back=.true.)+1:&
      len_trim( test_mesh_filename))),'light blue'), '...'

    ! Set up the mesh from the file (includes calculating secondary geometry data and matrix operators)
    call open_existing_netcdf_file_for_reading( trim(test_mesh_filename), ncid)
    call setup_mesh_from_file( test_mesh_filename, ncid, mesh)
    call close_netcdf_file( ncid)

    ! Check mesh self-consistency
    call assert( test_mesh_is_self_consistent( mesh), 'mesh is not self-consistent')

    ! Define a masked set of vertices
    allocate( mask_a( mesh%vi1:mesh%vi2), source = .false.)
    xmid = (mesh%xmin + mesh%xmax) / 2._dp
    ymid = (mesh%ymin + mesh%ymax) / 2._dp
    wx = mesh%xmax - xmid
    wy = mesh%ymax - ymid
    wmin = min( wx, wy)
    do vi = mesh%vi1, mesh%vi2
      if (hypot( mesh%V( vi,1) - xmid, mesh%V( vi,2) - ymid) < wmin * 0.75_dp) mask_a( vi) = .true.
    end do

    ! Create graphs from the masked vertices and triangles
    nz = 12
    call create_graph_from_masked_mesh_a( mesh, mask_a, nz, graph_a)
    call create_graph_from_masked_mesh_b( mesh, mask_a, nz, graph_b)

    call test_2nd_order_operators_on_regular_nodes( foldername_discretisation, mesh, graph_b)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_all_map_deriv_tests_on_graphs

  !> Test 2nd-order operators on the regular nodes of the graph from the mesh triangles
  subroutine test_2nd_order_operators_on_regular_nodes( foldername_discretisation, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_discretisation
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'test_2nd_order_operators_on_regular_nodes'
    character(len=1024)               :: graph_name_disp
    type(type_sparse_matrix_CSR_dp)   :: M2_ddx, M2_ddy, M2_d2dx2, M2_d2dxdy, M2_d2dy2
    procedure(test_function), pointer :: test_function_ptr
    character(len=1024)               :: function_name

    ! Add routine to call stack
    call init_routine( routine_name)

    graph_name_disp = graph%parent_mesh_name
    graph_name_disp = strrep( graph_name_disp, '.nc"', '')
    graph_name_disp = graph_name_disp( index( graph_name_disp,'/',back=.true.)+1: &
      len_trim( graph_name_disp))

    if (par%primary) write(0,*) '       Testing 2nd-order operators on the regular nodes of graph ', &
      UPSY%stru%colour_string( trim( graph_name_disp),'light blue'), '...'

    call calc_graph_matrix_operators_2nd_order( graph, M2_ddx, M2_ddy, M2_d2dx2, M2_d2dxdy, M2_d2dy2)

    ! Run all the tests with different test functions

    test_function_ptr => test_function_linear
    function_name = 'linear'
    call test_2nd_order_operators_on_regular_nodes_with_function( foldername_discretisation, &
      mesh, graph, M2_ddx, M2_ddy, M2_d2dx2, M2_d2dxdy, M2_d2dy2, test_function_ptr, function_name)

    test_function_ptr => test_function_quadratic
    function_name = 'quadratic'
    call test_2nd_order_operators_on_regular_nodes_with_function( foldername_discretisation, &
      mesh, graph, M2_ddx, M2_ddy, M2_d2dx2, M2_d2dxdy, M2_d2dy2, test_function_ptr, function_name)

    test_function_ptr => test_function_periodic
    function_name = 'periodic'
    call test_2nd_order_operators_on_regular_nodes_with_function( foldername_discretisation, &
      mesh, graph, M2_ddx, M2_ddy, M2_d2dx2, M2_d2dxdy, M2_d2dy2, test_function_ptr, function_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_2nd_order_operators_on_regular_nodes

  subroutine test_2nd_order_operators_on_regular_nodes_with_function( foldername_discretisation, &
    mesh, graph, M2_ddx, M2_ddy, M2_d2dx2, M2_d2dxdy, M2_d2dy2, test_function_ptr, function_name)

    ! In/output variables:
    character(len=*),                intent(in) :: foldername_discretisation
    type(type_mesh),                 intent(in) :: mesh
    type(type_graph),                intent(in) :: graph
    type(type_sparse_matrix_CSR_dp), intent(in) :: M2_ddx, M2_ddy, M2_d2dx2, M2_d2dxdy, M2_d2dy2
    procedure(test_function),           pointer :: test_function_ptr
    character(len=1024),             intent(in) :: function_name

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_2nd_order_operators_on_regular_nodes_with_function'

    real(dp), dimension(:), pointer :: d_ex      => null()
    real(dp), dimension(:), pointer :: ddx_ex    => null()
    real(dp), dimension(:), pointer :: ddy_ex    => null()
    real(dp), dimension(:), pointer :: d2dx2_ex  => null()
    real(dp), dimension(:), pointer :: d2dxdy_ex => null()
    real(dp), dimension(:), pointer :: d2dy2_ex  => null()
    type(MPI_WIN)                   :: wd_ex, wddx_ex, wddy_ex, wd2dx2_ex, wd2dxdy_ex, wd2dy2_ex
    real(dp), dimension(:), pointer :: ddx_disc    => null()
    real(dp), dimension(:), pointer :: ddy_disc    => null()
    real(dp), dimension(:), pointer :: d2dx2_disc  => null()
    real(dp), dimension(:), pointer :: d2dxdy_disc => null()
    real(dp), dimension(:), pointer :: d2dy2_disc  => null()
    type(MPI_WIN)                   :: wddx_disc, wddy_disc, wd2dx2_disc, wd2dxdy_disc, wd2dy2_disc
    integer                         :: ni
    real(dp)                        :: x,y
    real(dp), dimension(6)          :: d_all
    real(dp)                        :: d, ddx, ddy, d2dx2, d2dxdy, d2dy2

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '        Running all mapping/derivative tests on function ', &
      UPSY%stru%colour_string( trim( function_name),'light blue'), '...'

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_ex       , wd_ex       , graph%pai%n_nih)
    call allocate_dist_shared( ddx_ex     , wddx_ex     , graph%pai%n_nih)
    call allocate_dist_shared( ddy_ex     , wddy_ex     , graph%pai%n_nih)
    call allocate_dist_shared( d2dx2_ex   , wd2dx2_ex   , graph%pai%n_nih)
    call allocate_dist_shared( d2dxdy_ex  , wd2dxdy_ex  , graph%pai%n_nih)
    call allocate_dist_shared( d2dy2_ex   , wd2dy2_ex   , graph%pai%n_nih)

    call allocate_dist_shared( ddx_disc   , wddx_disc   , graph%pai%n_nih)
    call allocate_dist_shared( ddy_disc   , wddy_disc   , graph%pai%n_nih)
    call allocate_dist_shared( d2dx2_disc , wd2dx2_disc , graph%pai%n_nih)
    call allocate_dist_shared( d2dxdy_disc, wd2dxdy_disc, graph%pai%n_nih)
    call allocate_dist_shared( d2dy2_disc , wd2dy2_disc , graph%pai%n_nih)

    ! Calculate exact solutions
    do ni = graph%ni1, graph%ni2

      x = graph%V( ni,1)
      y = graph%V( ni,2)

      d_all = test_function_ptr( x, y, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax)
      d      = d_all( 1)
      ddx    = d_all( 2)
      ddy    = d_all( 3)
      d2dx2  = d_all( 4)
      d2dxdy = d_all( 5)
      d2dy2  = d_all( 6)

      d_ex     ( ni) = d
      ddx_ex   ( ni) = ddx
      ddy_ex   ( ni) = ddy
      d2dx2_ex ( ni) = d2dx2
      d2dxdy_ex( ni) = d2dxdy
      d2dy2_ex ( ni) = d2dy2
    end do

    ! Calculate discretised approximations
    call multiply_CSR_matrix_with_vector_1D_wrapper( M2_ddx, graph%pai, d_ex, graph%pai, ddx_disc, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = graph%buffer1_g_nih, buffer_yy_nih = graph%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( M2_ddy, graph%pai, d_ex, graph%pai, ddy_disc, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = graph%buffer1_g_nih, buffer_yy_nih = graph%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( M2_d2dx2, graph%pai, d_ex, graph%pai, d2dx2_disc, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = graph%buffer1_g_nih, buffer_yy_nih = graph%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( M2_d2dxdy, graph%pai, d_ex, graph%pai, d2dxdy_disc, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = graph%buffer1_g_nih, buffer_yy_nih = graph%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( M2_d2dy2, graph%pai, d_ex, graph%pai, d2dy2_disc, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = graph%buffer1_g_nih, buffer_yy_nih = graph%buffer2_g_nih)

    ! Write results to NetCDF
    call test_2nd_or_op_on_reg_nod_w_func_write_to_file( &
      foldername_discretisation, graph, function_name, &
      d_ex, ddx_ex, ddy_ex, d2dx2_ex, d2dxdy_ex, d2dy2_ex, &
      ddx_disc, ddy_disc, d2dx2_disc, d2dxdy_disc, d2dy2_disc)

    ! Clean up after yourself
    call deallocate_dist_shared( d_ex       , wd_ex       )
    call deallocate_dist_shared( ddx_ex     , wddx_ex     )
    call deallocate_dist_shared( ddy_ex     , wddy_ex     )
    call deallocate_dist_shared( d2dx2_ex   , wd2dx2_ex   )
    call deallocate_dist_shared( d2dxdy_ex  , wd2dxdy_ex  )
    call deallocate_dist_shared( d2dy2_ex   , wd2dy2_ex   )

    call deallocate_dist_shared( ddx_disc   , wddx_disc   )
    call deallocate_dist_shared( ddy_disc   , wddy_disc   )
    call deallocate_dist_shared( d2dx2_disc , wd2dx2_disc )
    call deallocate_dist_shared( d2dxdy_disc, wd2dxdy_disc)
    call deallocate_dist_shared( d2dy2_disc , wd2dy2_disc )

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_2nd_order_operators_on_regular_nodes_with_function

  subroutine test_2nd_or_op_on_reg_nod_w_func_write_to_file( &
    foldername_discretisation, graph, function_name, &
    d_ex, ddx_ex, ddy_ex, d2dx2_ex, d2dxdy_ex, d2dy2_ex, &
    ddx_disc, ddy_disc, d2dx2_disc, d2dxdy_disc, d2dy2_disc)

    ! In/output variables:
    character(len=*),       intent(in) :: foldername_discretisation
    type(type_graph),       intent(in) :: graph
    character(len=1024),    intent(in) :: function_name
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), target, intent(in) :: d_ex
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), target, intent(in) :: ddx_ex
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), target, intent(in) :: ddy_ex
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), target, intent(in) :: d2dx2_ex
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), target, intent(in) :: d2dxdy_ex
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), target, intent(in) :: d2dy2_ex
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), target, intent(in) :: ddx_disc
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), target, intent(in) :: ddy_disc
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), target, intent(in) :: d2dx2_disc
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), target, intent(in) :: d2dxdy_disc
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), target, intent(in) :: d2dy2_disc

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'test_2nd_or_op_on_reg_nod_w_func_write_to_file'
    character(len=1024)             :: graph_name
    integer                         :: i
    character(len=1024)             :: filename
    real(dp), dimension(:), pointer :: d_ex_loc
    real(dp), dimension(:), pointer :: ddx_ex_loc
    real(dp), dimension(:), pointer :: ddy_ex_loc
    real(dp), dimension(:), pointer :: d2dx2_ex_loc
    real(dp), dimension(:), pointer :: d2dxdy_ex_loc
    real(dp), dimension(:), pointer :: d2dy2_ex_loc
    real(dp), dimension(:), pointer :: ddx_disc_loc
    real(dp), dimension(:), pointer :: ddy_disc_loc
    real(dp), dimension(:), pointer :: d2dx2_disc_loc
    real(dp), dimension(:), pointer :: d2dxdy_disc_loc
    real(dp), dimension(:), pointer :: d2dy2_disc_loc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Save the graph
    graph_name = graph%parent_mesh_name
    graph_name = strrep( graph_name, '.nc"', '')
    i = index( graph_name, '/', back = .true.)
    graph_name = graph_name( i+1:len_trim( graph_name))
    filename = trim( foldername_discretisation) // '/res_' // &
      trim( graph_name) // '_2nd_order_' // trim( function_name) // '_graph.nc'

    call save_graph_as_netcdf( filename, graph)

    ! Point to process-local data range
    ! (as the netcdf code cannot handle hybrid distributed/shared data yet)

    d_ex_loc        => d_ex       ( graph%pai%i1: graph%pai%i2)
    ddx_ex_loc      => ddx_ex     ( graph%pai%i1: graph%pai%i2)
    ddy_ex_loc      => ddy_ex     ( graph%pai%i1: graph%pai%i2)
    d2dx2_ex_loc    => d2dx2_ex   ( graph%pai%i1: graph%pai%i2)
    d2dxdy_ex_loc   => d2dxdy_ex  ( graph%pai%i1: graph%pai%i2)
    d2dy2_ex_loc    => d2dy2_ex   ( graph%pai%i1: graph%pai%i2)
    ddx_disc_loc    => ddx_disc   ( graph%pai%i1: graph%pai%i2)
    ddy_disc_loc    => ddy_disc   ( graph%pai%i1: graph%pai%i2)
    d2dx2_disc_loc  => d2dx2_disc ( graph%pai%i1: graph%pai%i2)
    d2dxdy_disc_loc => d2dxdy_disc( graph%pai%i1: graph%pai%i2)
    d2dy2_disc_loc  => d2dy2_disc ( graph%pai%i1: graph%pai%i2)

    ! Save data
    call save_variable_as_netcdf_dp_1D( trim( foldername_discretisation), d_ex_loc, &
      'res_' // trim( graph_name) // '_2nd_order_' // trim( function_name) // '_d_ex')
    call save_variable_as_netcdf_dp_1D( trim( foldername_discretisation), ddx_ex_loc, &
      'res_' // trim( graph_name) // '_2nd_order_' // trim( function_name) // '_ddx_ex')
    call save_variable_as_netcdf_dp_1D( trim( foldername_discretisation), ddy_ex_loc, &
      'res_' // trim( graph_name) // '_2nd_order_' // trim( function_name) // '_ddy_ex')
    call save_variable_as_netcdf_dp_1D( trim( foldername_discretisation), d2dx2_ex_loc, &
      'res_' // trim( graph_name) // '_2nd_order_' // trim( function_name) // '_d2dx2_ex')
    call save_variable_as_netcdf_dp_1D( trim( foldername_discretisation), d2dxdy_ex_loc, &
      'res_' // trim( graph_name) // '_2nd_order_' // trim( function_name) // '_d2dxdy_ex')
    call save_variable_as_netcdf_dp_1D( trim( foldername_discretisation), d2dy2_ex_loc, &
      'res_' // trim( graph_name) // '_2nd_order_' // trim( function_name) // '_d2dy2_ex')
    call save_variable_as_netcdf_dp_1D( trim( foldername_discretisation), ddx_disc_loc, &
      'res_' // trim( graph_name) // '_2nd_order_' // trim( function_name) // '_ddx_disc')
    call save_variable_as_netcdf_dp_1D( trim( foldername_discretisation), ddy_disc_loc, &
      'res_' // trim( graph_name) // '_2nd_order_' // trim( function_name) // '_ddy_disc')
    call save_variable_as_netcdf_dp_1D( trim( foldername_discretisation), d2dx2_disc_loc, &
      'res_' // trim( graph_name) // '_2nd_order_' // trim( function_name) // '_d2dx2_disc')
    call save_variable_as_netcdf_dp_1D( trim( foldername_discretisation), d2dxdy_disc_loc, &
      'res_' // trim( graph_name) // '_2nd_order_' // trim( function_name) // '_d2dxdy_disc')
    call save_variable_as_netcdf_dp_1D( trim( foldername_discretisation), d2dy2_disc_loc, &
      'res_' // trim( graph_name) // '_2nd_order_' // trim( function_name) // '_d2dy2_disc')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_2nd_or_op_on_reg_nod_w_func_write_to_file

end module ct_discretisation_mapping_derivatives_graph
