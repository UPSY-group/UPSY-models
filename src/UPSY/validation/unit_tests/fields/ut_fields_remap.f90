module ut_fields_remap

  use precisions, only: dp
  use mpi_basic, only: par
  use parameters, only: pi
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use fields_main, only: third_dimension, &
    atype_field, atype_field_2D, atype_field_3D, &
    type_field_logical_2D, type_field_int_2D, type_field_dp_2D, &
    type_field_logical_3D, type_field_int_3D, type_field_dp_3D, &
    type_fields_registry
  use ut_basic
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use grid_basic, only: setup_square_grid
  use Arakawa_grid_mod, only: Arakawa_grid
  use mpi_f08, only: MPI_WIN, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, &
    MPI_MIN, MPI_MAX, MPI_COMM_WORLD
  use mesh_memory, only: allocate_mesh_primary, crop_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_Lloyds_algorithm, only: Lloyds_algorithm_single_iteration
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_disc_calc_matrix_operators_2D, only: calc_all_matrix_operators_mesh
  use tests_main

  implicit none

  private

  public :: test_remap_field

contains

  subroutine test_remap_field( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_remap_field'
    character(len=1024), parameter :: test_name_local = 'remap_field'
    character(len=1024)            :: test_name
    real(dp)                       :: alpha_min, res_max
    real(dp), parameter            :: xmin = -1._dp
    real(dp), parameter            :: xmax =  1._dp
    real(dp), parameter            :: ymin = -1._dp
    real(dp), parameter            :: ymax =  1._dp
    type(type_mesh)                :: mesh_src, mesh_dst

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Source mesh
    ! ===========

    call allocate_mesh_primary( mesh_src, 'dummy_mesh_1', 100, 200)
    call initialise_dummy_mesh_5( mesh_src, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 23.2_dp
    call refine_mesh_uniform( mesh_src, res_max, alpha_min)
    call Lloyds_algorithm_single_iteration( mesh_src, alpha_min)
    call crop_mesh_primary( mesh_src)
    call calc_all_secondary_mesh_data( mesh_src, 0._dp, -90._dp, 71._dp)
    call calc_all_matrix_operators_mesh( mesh_src)

    ! Destination mesh
    ! ================

    call allocate_mesh_primary( mesh_dst, 'dummy_mesh_2', 100, 200)
    call initialise_dummy_mesh_5( mesh_dst, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 18._dp
    call refine_mesh_uniform( mesh_dst, res_max, alpha_min)
    call Lloyds_algorithm_single_iteration( mesh_dst, alpha_min)
    call crop_mesh_primary( mesh_dst)
    call calc_all_secondary_mesh_data( mesh_dst, 0._dp, -90._dp, 71._dp)
    call calc_all_matrix_operators_mesh( mesh_dst)

    ! Remapping tests
    ! ===============

    call test_remap_field_mesh_dp_a_2D      ( test_name, mesh_src, mesh_dst)
    call test_remap_field_mesh_dp_a_3D_zeta ( test_name, mesh_src, mesh_dst)
    call test_remap_field_mesh_dp_a_3D_month( test_name, mesh_src, mesh_dst)
    call test_remap_field_mesh_dp_a_3D_ocean( test_name, mesh_src, mesh_dst)

    call test_remap_field_mesh_dp_b_2D      ( test_name, mesh_src, mesh_dst)
    call test_remap_field_mesh_dp_b_3D_zeta ( test_name, mesh_src, mesh_dst)
    call test_remap_field_mesh_dp_b_3D_month( test_name, mesh_src, mesh_dst)
    call test_remap_field_mesh_dp_b_3D_ocean( test_name, mesh_src, mesh_dst)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remap_field



  subroutine test_remap_field_mesh_dp_a_2D( test_name_parent, mesh_src, mesh_dst)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh_src, mesh_dst

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_remap_field_mesh_dp_a_2D'
    character(len=1024), parameter                :: test_name_local = 'mesh_dp_a_2D'
    character(len=1024)                           :: test_name
    type(type_fields_registry)                    :: flds_reg
    character(len=1024)                           :: name, long_name, units
    real(dp), dimension(:  ), contiguous, pointer :: d => null()
    type(MPI_WIN)                                 :: wd
    integer                                       :: i, vi, ierr
    integer                                       :: lb_a, ub_a
    integer                                       :: lb_f, ub_f
    real(dp)                                      :: dmin, dmax

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd'
    long_name = 'd_long_name'
    units     = 'd_units'

    call flds_reg%create_field( d, wd, &
      mesh_src, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = '2nd_order_conservative')

    i = flds_reg%find( name)

    do vi = mesh_src%vi1, mesh_src%vi2
      d( vi) = test_function( mesh_src%xmin, mesh_src%xmax, mesh_src%ymin, mesh_src%ymax, &
        mesh_src%V( vi,1), mesh_src%V( vi,2))
    end do

    call flds_reg%remap_field( mesh_dst, 'd', d)

    lb_a = lbound( d,1)
    ub_a = ubound( d,1)

    lb_f = flds_reg%items(i)%p%lbound( 1)
    ub_f = flds_reg%items(i)%p%ubound( 1)

    dmin = minval( d( mesh_dst%vi1: mesh_dst%vi2))
    dmax = maxval( d( mesh_dst%vi1: mesh_dst%vi2))
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call unit_test( (&
      flds_reg%items(i)%p%name()      == name .and. &
      flds_reg%items(i)%p%long_name() == long_name .and. &
      flds_reg%items(i)%p%units()     == units .and. &
      flds_reg%items(i)%p%is_grid( mesh_dst) .and. &
      flds_reg%items(i)%p%is_Arakawa_grid( Arakawa_grid%a()) .and. &
      flds_reg%items(i)%p%is_pai( mesh_dst%pai_V) .and. &
      lb_a == mesh_dst%pai_V%i1_nih .and. &
      ub_a == mesh_dst%pai_V%i2_nih .and. &
      lb_f == mesh_dst%pai_V%i1_nih .and. &
      ub_f == mesh_dst%pai_V%i2_nih .and. &
      test_ge_le( dmin, -1.15_dp, -0.85_dp) .and. &
      test_ge_le( dmax,  0.85_dp,  1.15_dp)), &
      trim( test_name))

    ! Clean up after yourself
    call flds_reg%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remap_field_mesh_dp_a_2D

  subroutine test_remap_field_mesh_dp_a_3D_zeta( test_name_parent, mesh_src, mesh_dst)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh_src, mesh_dst

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_remap_field_mesh_dp_a_3D_zeta'
    character(len=1024), parameter                :: test_name_local = 'mesh_dp_a_3D_zeta'
    character(len=1024)                           :: test_name
    type(type_fields_registry)                    :: flds_reg
    character(len=1024)                           :: name, long_name, units
    integer, parameter                            :: nz = 10
    real(dp), dimension(:,:), contiguous, pointer :: d => null()
    type(MPI_WIN)                                 :: wd
    integer                                       :: i, vi, ierr
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f
    real(dp)                                      :: dmin, dmax

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd'
    long_name = 'd_long_name'
    units     = 'd_units'

    call flds_reg%create_field( d, wd, &
      mesh_src, Arakawa_grid%a(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = '2nd_order_conservative')

    i = flds_reg%find( name)

    do vi = mesh_src%vi1, mesh_src%vi2
      d( vi,:) = test_function( mesh_src%xmin, mesh_src%xmax, mesh_src%ymin, mesh_src%ymax, &
        mesh_src%V( vi,1), mesh_src%V( vi,2))
    end do

    call flds_reg%remap_field( mesh_dst, 'd', d)

    lb1_a = lbound( d,1)
    ub1_a = ubound( d,1)

    lb2_a = lbound( d,2)
    ub2_a = ubound( d,2)

    lb1_f = flds_reg%items(i)%p%lbound( 1)
    ub1_f = flds_reg%items(i)%p%ubound( 1)

    lb2_f = flds_reg%items(i)%p%lbound( 2)
    ub2_f = flds_reg%items(i)%p%ubound( 2)

    dmin = minval( d( mesh_dst%vi1: mesh_dst%vi2,:))
    dmax = maxval( d( mesh_dst%vi1: mesh_dst%vi2,:))
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call unit_test( (&
      flds_reg%items(i)%p%name()      == name .and. &
      flds_reg%items(i)%p%long_name() == long_name .and. &
      flds_reg%items(i)%p%units()     == units .and. &
      flds_reg%items(i)%p%is_grid( mesh_dst) .and. &
      flds_reg%items(i)%p%is_Arakawa_grid( Arakawa_grid%a()) .and. &
      flds_reg%items(i)%p%is_third_dimension( third_dimension%ice_zeta( nz,'regular')) .and. &
      flds_reg%items(i)%p%is_pai( mesh_dst%pai_V) .and. &
      lb1_a == mesh_dst%pai_V%i1_nih .and. &
      ub1_a == mesh_dst%pai_V%i2_nih .and. &
      lb1_f == mesh_dst%pai_V%i1_nih .and. &
      ub1_f == mesh_dst%pai_V%i2_nih .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      test_ge_le( dmin, -1.15_dp, -0.85_dp) .and. &
      test_ge_le( dmax,  0.85_dp,  1.15_dp)), &
      trim( test_name))

    ! Clean up after yourself
    call flds_reg%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remap_field_mesh_dp_a_3D_zeta

  subroutine test_remap_field_mesh_dp_a_3D_month( test_name_parent, mesh_src, mesh_dst)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh_src, mesh_dst

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_remap_field_mesh_dp_a_3D_month'
    character(len=1024), parameter                :: test_name_local = 'mesh_dp_a_3D_month'
    character(len=1024)                           :: test_name
    type(type_fields_registry)                    :: flds_reg
    character(len=1024)                           :: name, long_name, units
    real(dp), dimension(:,:), contiguous, pointer :: d => null()
    type(MPI_WIN)                                 :: wd
    integer                                       :: i, vi, ierr
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f
    real(dp)                                      :: dmin, dmax

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd'
    long_name = 'd_long_name'
    units     = 'd_units'

    call flds_reg%create_field( d, wd, &
      mesh_src, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = '2nd_order_conservative')

    i = flds_reg%find( name)

    do vi = mesh_src%vi1, mesh_src%vi2
      d( vi,:) = test_function( mesh_src%xmin, mesh_src%xmax, mesh_src%ymin, mesh_src%ymax, &
        mesh_src%V( vi,1), mesh_src%V( vi,2))
    end do

    call flds_reg%remap_field( mesh_dst, 'd', d)

    lb1_a = lbound( d,1)
    ub1_a = ubound( d,1)

    lb2_a = lbound( d,2)
    ub2_a = ubound( d,2)

    lb1_f = flds_reg%items(i)%p%lbound( 1)
    ub1_f = flds_reg%items(i)%p%ubound( 1)

    lb2_f = flds_reg%items(i)%p%lbound( 2)
    ub2_f = flds_reg%items(i)%p%ubound( 2)

    dmin = minval( d( mesh_dst%vi1: mesh_dst%vi2,:))
    dmax = maxval( d( mesh_dst%vi1: mesh_dst%vi2,:))
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call unit_test( (&
      flds_reg%items(i)%p%name()      == name .and. &
      flds_reg%items(i)%p%long_name() == long_name .and. &
      flds_reg%items(i)%p%units()     == units .and. &
      flds_reg%items(i)%p%is_grid( mesh_dst) .and. &
      flds_reg%items(i)%p%is_Arakawa_grid( Arakawa_grid%a()) .and. &
      flds_reg%items(i)%p%is_third_dimension( third_dimension%month()) .and. &
      flds_reg%items(i)%p%is_pai( mesh_dst%pai_V) .and. &
      lb1_a == mesh_dst%pai_V%i1_nih .and. &
      ub1_a == mesh_dst%pai_V%i2_nih .and. &
      lb1_f == mesh_dst%pai_V%i1_nih .and. &
      ub1_f == mesh_dst%pai_V%i2_nih .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      test_ge_le( dmin, -1.15_dp, -0.85_dp) .and. &
      test_ge_le( dmax,  0.85_dp,  1.15_dp)), &
      trim( test_name))

    ! Clean up after yourself
    call flds_reg%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remap_field_mesh_dp_a_3D_month

  subroutine test_remap_field_mesh_dp_a_3D_ocean( test_name_parent, mesh_src, mesh_dst)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh_src, mesh_dst

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_remap_field_mesh_dp_a_3D_ocean'
    character(len=1024), parameter                :: test_name_local = 'mesh_dp_a_3D_ocean'
    character(len=1024)                           :: test_name
    type(type_fields_registry)                    :: flds_reg
    character(len=1024)                           :: name, long_name, units
    integer, parameter                            :: nz = 15
    real(dp), dimension(:,:), contiguous, pointer :: d => null()
    type(MPI_WIN)                                 :: wd
    integer                                       :: i, vi, ierr
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f
    real(dp)                                      :: dmin, dmax

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd'
    long_name = 'd_long_name'
    units     = 'd_units'

    call flds_reg%create_field( d, wd, &
      mesh_src, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = '2nd_order_conservative')

    i = flds_reg%find( name)

    do vi = mesh_src%vi1, mesh_src%vi2
      d( vi,:) = test_function( mesh_src%xmin, mesh_src%xmax, mesh_src%ymin, mesh_src%ymax, &
        mesh_src%V( vi,1), mesh_src%V( vi,2))
    end do

    call flds_reg%remap_field( mesh_dst, 'd', d)

    lb1_a = lbound( d,1)
    ub1_a = ubound( d,1)

    lb2_a = lbound( d,2)
    ub2_a = ubound( d,2)

    lb1_f = flds_reg%items(i)%p%lbound( 1)
    ub1_f = flds_reg%items(i)%p%ubound( 1)

    lb2_f = flds_reg%items(i)%p%lbound( 2)
    ub2_f = flds_reg%items(i)%p%ubound( 2)

    dmin = minval( d( mesh_dst%vi1: mesh_dst%vi2,:))
    dmax = maxval( d( mesh_dst%vi1: mesh_dst%vi2,:))
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call unit_test( (&
      flds_reg%items(i)%p%name()      == name .and. &
      flds_reg%items(i)%p%long_name() == long_name .and. &
      flds_reg%items(i)%p%units()     == units .and. &
      flds_reg%items(i)%p%is_grid( mesh_dst) .and. &
      flds_reg%items(i)%p%is_Arakawa_grid( Arakawa_grid%a()) .and. &
      flds_reg%items(i)%p%is_third_dimension( third_dimension%ocean_depth( nz)) .and. &
      flds_reg%items(i)%p%is_pai( mesh_dst%pai_V) .and. &
      lb1_a == mesh_dst%pai_V%i1_nih .and. &
      ub1_a == mesh_dst%pai_V%i2_nih .and. &
      lb1_f == mesh_dst%pai_V%i1_nih .and. &
      ub1_f == mesh_dst%pai_V%i2_nih .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      test_ge_le( dmin, -1.15_dp, -0.85_dp) .and. &
      test_ge_le( dmax,  0.85_dp,  1.15_dp)), &
      trim( test_name))

    ! Clean up after yourself
    call flds_reg%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remap_field_mesh_dp_a_3D_ocean



  subroutine test_remap_field_mesh_dp_b_2D( test_name_parent, mesh_src, mesh_dst)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh_src, mesh_dst

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_remap_field_mesh_dp_b_2D'
    character(len=1024), parameter                :: test_name_local = 'mesh_dp_b_2D'
    character(len=1024)                           :: test_name
    type(type_fields_registry)                    :: flds_reg
    character(len=1024)                           :: name, long_name, units
    real(dp), dimension(:  ), contiguous, pointer :: d => null()
    type(MPI_WIN)                                 :: wd
    integer                                       :: i, ti, ierr
    integer                                       :: lb_a, ub_a
    integer                                       :: lb_f, ub_f
    real(dp)                                      :: dmin, dmax

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd'
    long_name = 'd_long_name'
    units     = 'd_units'

    call flds_reg%create_field( d, wd, &
      mesh_src, Arakawa_grid%b(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = '2nd_order_conservative')

    i = flds_reg%find( name)

    do ti = mesh_src%ti1, mesh_src%ti2
      d( ti) = test_function( mesh_src%xmin, mesh_src%xmax, mesh_src%ymin, mesh_src%ymax, &
        mesh_src%TriGC( ti,1), mesh_src%TriGC( ti,2))
    end do

    call flds_reg%remap_field( mesh_dst, 'd', d)

    lb_a = lbound( d,1)
    ub_a = ubound( d,1)

    lb_f = flds_reg%items(i)%p%lbound( 1)
    ub_f = flds_reg%items(i)%p%ubound( 1)

    dmin = minval( d( mesh_dst%ti1: mesh_dst%ti2))
    dmax = maxval( d( mesh_dst%ti1: mesh_dst%ti2))
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call unit_test( (&
      flds_reg%items(i)%p%name()      == name .and. &
      flds_reg%items(i)%p%long_name() == long_name .and. &
      flds_reg%items(i)%p%units()     == units .and. &
      flds_reg%items(i)%p%is_grid( mesh_dst) .and. &
      flds_reg%items(i)%p%is_Arakawa_grid( Arakawa_grid%b()) .and. &
      flds_reg%items(i)%p%is_pai( mesh_dst%pai_Tri) .and. &
      lb_a == mesh_dst%pai_Tri%i1_nih .and. &
      ub_a == mesh_dst%pai_Tri%i2_nih .and. &
      lb_f == mesh_dst%pai_Tri%i1_nih .and. &
      ub_f == mesh_dst%pai_Tri%i2_nih .and. &
      test_ge_le( dmin, -1.15_dp, -0.85_dp) .and. &
      test_ge_le( dmax,  0.85_dp,  1.15_dp)), &
      trim( test_name))

    ! Clean up after yourself
    call flds_reg%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remap_field_mesh_dp_b_2D

  subroutine test_remap_field_mesh_dp_b_3D_zeta( test_name_parent, mesh_src, mesh_dst)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh_src, mesh_dst

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_remap_field_mesh_dp_b_3D_zeta'
    character(len=1024), parameter                :: test_name_local = 'mesh_dp_b_3D_zeta'
    character(len=1024)                           :: test_name
    type(type_fields_registry)                    :: flds_reg
    character(len=1024)                           :: name, long_name, units
    integer, parameter                            :: nz = 10
    real(dp), dimension(:,:), contiguous, pointer :: d => null()
    type(MPI_WIN)                                 :: wd
    integer                                       :: i, ti, ierr
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f
    real(dp)                                      :: dmin, dmax

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd'
    long_name = 'd_long_name'
    units     = 'd_units'

    call flds_reg%create_field( d, wd, &
      mesh_src, Arakawa_grid%b(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = '2nd_order_conservative')

    i = flds_reg%find( name)

    do ti = mesh_src%ti1, mesh_src%ti2
      d( ti,:) = test_function( mesh_src%xmin, mesh_src%xmax, mesh_src%ymin, mesh_src%ymax, &
        mesh_src%TriGC( ti,1), mesh_src%TriGC( ti,2))
    end do

    call flds_reg%remap_field( mesh_dst, 'd', d)

    lb1_a = lbound( d,1)
    ub1_a = ubound( d,1)

    lb2_a = lbound( d,2)
    ub2_a = ubound( d,2)

    lb1_f = flds_reg%items(i)%p%lbound( 1)
    ub1_f = flds_reg%items(i)%p%ubound( 1)

    lb2_f = flds_reg%items(i)%p%lbound( 2)
    ub2_f = flds_reg%items(i)%p%ubound( 2)

    dmin = minval( d( mesh_dst%ti1: mesh_dst%ti2,:))
    dmax = maxval( d( mesh_dst%ti1: mesh_dst%ti2,:))
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call unit_test( (&
      flds_reg%items(i)%p%name()      == name .and. &
      flds_reg%items(i)%p%long_name() == long_name .and. &
      flds_reg%items(i)%p%units()     == units .and. &
      flds_reg%items(i)%p%is_grid( mesh_dst) .and. &
      flds_reg%items(i)%p%is_Arakawa_grid( Arakawa_grid%b()) .and. &
      flds_reg%items(i)%p%is_third_dimension( third_dimension%ice_zeta( nz,'regular')) .and. &
      flds_reg%items(i)%p%is_pai( mesh_dst%pai_Tri) .and. &
      lb1_a == mesh_dst%pai_Tri%i1_nih .and. &
      ub1_a == mesh_dst%pai_Tri%i2_nih .and. &
      lb1_f == mesh_dst%pai_Tri%i1_nih .and. &
      ub1_f == mesh_dst%pai_Tri%i2_nih .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      test_ge_le( dmin, -1.15_dp, -0.85_dp) .and. &
      test_ge_le( dmax,  0.85_dp,  1.15_dp)), &
      trim( test_name))

    ! Clean up after yourself
    call flds_reg%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remap_field_mesh_dp_b_3D_zeta

  subroutine test_remap_field_mesh_dp_b_3D_month( test_name_parent, mesh_src, mesh_dst)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh_src, mesh_dst

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_remap_field_mesh_dp_b_3D_month'
    character(len=1024), parameter                :: test_name_local = 'mesh_dp_b_3D_month'
    character(len=1024)                           :: test_name
    type(type_fields_registry)                    :: flds_reg
    character(len=1024)                           :: name, long_name, units
    real(dp), dimension(:,:), contiguous, pointer :: d => null()
    type(MPI_WIN)                                 :: wd
    integer                                       :: i, ti, ierr
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f
    real(dp)                                      :: dmin, dmax

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd'
    long_name = 'd_long_name'
    units     = 'd_units'

    call flds_reg%create_field( d, wd, &
      mesh_src, Arakawa_grid%b(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = '2nd_order_conservative')

    i = flds_reg%find( name)

    do ti = mesh_src%ti1, mesh_src%ti2
      d( ti,:) = test_function( mesh_src%xmin, mesh_src%xmax, mesh_src%ymin, mesh_src%ymax, &
        mesh_src%TriGC( ti,1), mesh_src%TriGC( ti,2))
    end do

    call flds_reg%remap_field( mesh_dst, 'd', d)

    lb1_a = lbound( d,1)
    ub1_a = ubound( d,1)

    lb2_a = lbound( d,2)
    ub2_a = ubound( d,2)

    lb1_f = flds_reg%items(i)%p%lbound( 1)
    ub1_f = flds_reg%items(i)%p%ubound( 1)

    lb2_f = flds_reg%items(i)%p%lbound( 2)
    ub2_f = flds_reg%items(i)%p%ubound( 2)

    dmin = minval( d( mesh_dst%ti1: mesh_dst%ti2,:))
    dmax = maxval( d( mesh_dst%ti1: mesh_dst%ti2,:))
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call unit_test( (&
      flds_reg%items(i)%p%name()      == name .and. &
      flds_reg%items(i)%p%long_name() == long_name .and. &
      flds_reg%items(i)%p%units()     == units .and. &
      flds_reg%items(i)%p%is_grid( mesh_dst) .and. &
      flds_reg%items(i)%p%is_Arakawa_grid( Arakawa_grid%b()) .and. &
      flds_reg%items(i)%p%is_third_dimension( third_dimension%month()) .and. &
      flds_reg%items(i)%p%is_pai( mesh_dst%pai_Tri) .and. &
      lb1_a == mesh_dst%pai_Tri%i1_nih .and. &
      ub1_a == mesh_dst%pai_Tri%i2_nih .and. &
      lb1_f == mesh_dst%pai_Tri%i1_nih .and. &
      ub1_f == mesh_dst%pai_Tri%i2_nih .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      test_ge_le( dmin, -1.15_dp, -0.85_dp) .and. &
      test_ge_le( dmax,  0.85_dp,  1.15_dp)), &
      trim( test_name))

    ! Clean up after yourself
    call flds_reg%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remap_field_mesh_dp_b_3D_month

  subroutine test_remap_field_mesh_dp_b_3D_ocean( test_name_parent, mesh_src, mesh_dst)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh_src, mesh_dst

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_remap_field_mesh_dp_b_3D_ocean'
    character(len=1024), parameter                :: test_name_local = 'mesh_dp_b_3D_ocean'
    character(len=1024)                           :: test_name
    type(type_fields_registry)                    :: flds_reg
    character(len=1024)                           :: name, long_name, units
    integer, parameter                            :: nz = 15
    real(dp), dimension(:,:), contiguous, pointer :: d => null()
    type(MPI_WIN)                                 :: wd
    integer                                       :: i, ti, ierr
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f
    real(dp)                                      :: dmin, dmax

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    name      = 'd'
    long_name = 'd_long_name'
    units     = 'd_units'

    call flds_reg%create_field( d, wd, &
      mesh_src, Arakawa_grid%b(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units, &
      remap_method = '2nd_order_conservative')

    i = flds_reg%find( name)

    do ti = mesh_src%ti1, mesh_src%ti2
      d( ti,:) = test_function( mesh_src%xmin, mesh_src%xmax, mesh_src%ymin, mesh_src%ymax, &
        mesh_src%TriGC( ti,1), mesh_src%TriGC( ti,2))
    end do

    call flds_reg%remap_field( mesh_dst, 'd', d)

    lb1_a = lbound( d,1)
    ub1_a = ubound( d,1)

    lb2_a = lbound( d,2)
    ub2_a = ubound( d,2)

    lb1_f = flds_reg%items(i)%p%lbound( 1)
    ub1_f = flds_reg%items(i)%p%ubound( 1)

    lb2_f = flds_reg%items(i)%p%lbound( 2)
    ub2_f = flds_reg%items(i)%p%ubound( 2)

    dmin = minval( d( mesh_dst%ti1: mesh_dst%ti2,:))
    dmax = maxval( d( mesh_dst%ti1: mesh_dst%ti2,:))
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, dmax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierr)

    call unit_test( (&
      flds_reg%items(i)%p%name()      == name .and. &
      flds_reg%items(i)%p%long_name() == long_name .and. &
      flds_reg%items(i)%p%units()     == units .and. &
      flds_reg%items(i)%p%is_grid( mesh_dst) .and. &
      flds_reg%items(i)%p%is_Arakawa_grid( Arakawa_grid%b()) .and. &
      flds_reg%items(i)%p%is_third_dimension( third_dimension%ocean_depth( nz)) .and. &
      flds_reg%items(i)%p%is_pai( mesh_dst%pai_Tri) .and. &
      lb1_a == mesh_dst%pai_Tri%i1_nih .and. &
      ub1_a == mesh_dst%pai_Tri%i2_nih .and. &
      lb1_f == mesh_dst%pai_Tri%i1_nih .and. &
      ub1_f == mesh_dst%pai_Tri%i2_nih .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      test_ge_le( dmin, -1.15_dp, -0.85_dp) .and. &
      test_ge_le( dmax,  0.85_dp,  1.15_dp)), &
      trim( test_name))

    ! Clean up after yourself
    call flds_reg%destroy

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_remap_field_mesh_dp_b_3D_ocean



  function test_function( xmin, xmax, ymin, ymax, x, y) result( f)
    real(dp), intent(in) :: xmin, xmax, ymin, ymax, x, y
    real(dp)             :: cx,cy,f
    cx = xmax - xmin
    cy = ymax - ymin
    f = cos( x * 4._dp * pi / cx) * cos( y * 5._dp * pi / cy)
  end function test_function

end module ut_fields_remap