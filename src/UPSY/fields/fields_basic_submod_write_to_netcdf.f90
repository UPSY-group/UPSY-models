submodule (fields_basic) fields_basic_submod_write_to_netcdf

  use parameters, only: NaN
  use gather_dist_shared_to_primary_mod, only: gather_dist_shared_to_primary
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary
  use netcdf_io_main
  use netcdf, only: NF90_INT, NF90_DOUBLE

contains

  subroutine write_to_netcdf( self, filename, ncid)
    ! NOTE: assumes the NetCDF file already exists, is open, and
    !       already contains the grid/mesh and third dimension

    ! In/output variables:
    class(atype_field), intent(in) :: self
    character(len=*),   intent(in) :: filename
    integer,            intent(in) :: ncid

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'write_to_netcdf'
    integer                            :: var_type
    integer, dimension(:), allocatable :: dim_ids
    integer                            :: id_var

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f => self)
    class default
      call crash('invalid field type')
    class is (atype_field_2D)
    class is (atype_field_3D)
      call f%third_dimension_val%create_dim_and_var_in_netcdf( filename, ncid)
    end select

    var_type = get_var_type( self)
    dim_ids  = get_dim_ids( self, filename, ncid)
    call create_variable( filename, ncid, trim( self%name()), var_type, dim_ids, id_var)

    select type (g => self%grid())
    class default
      call crash('invalid field%grid type')
    class is (type_grid)
      call write_to_netcdf_grid( self, g, filename, ncid, id_var)
    class is (type_mesh)
      call write_to_netcdf_mesh( self, g, filename, ncid, id_var)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf

  ! ===== Grid-based fields
  ! =======================

  subroutine write_to_netcdf_grid( field, grid, filename, ncid, id_var)

    ! In/output variables:
    class(atype_field), intent(in) :: field
    type(type_grid),    intent(in) :: grid
    character(len=*),   intent(in) :: filename
    integer,            intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_netcdf_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f => field)
    class default
      call crash('invalid field type')
    class is (type_field_logical_2D)
      call write_to_netcdf_grid_logical_2D( f, grid, filename, ncid, id_var)
    class is (type_field_int_2D)
      call write_to_netcdf_grid_int_2D    ( f, grid, filename, ncid, id_var)
    class is (type_field_dp_2D)
      call write_to_netcdf_grid_dp_2D     ( f, grid, filename, ncid, id_var)
    class is (type_field_logical_3D)
      call write_to_netcdf_grid_logical_3D( f, grid, filename, ncid, id_var)
    class is (type_field_int_3D)
      call write_to_netcdf_grid_int_3D    ( f, grid, filename, ncid, id_var)
    class is (type_field_dp_3D)
      call write_to_netcdf_grid_dp_3D     ( f, grid, filename, ncid, id_var)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_grid

  subroutine write_to_netcdf_grid_logical_2D( field, grid, filename, ncid, id_var)

    ! In/output variables:
    class(type_field_logical_2D), intent(in) :: field
    type(type_grid),              intent(in) :: grid
    character(len=*),             intent(in) :: filename
    integer,                      intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'write_to_netcdf_grid_logical_2D'
    logical, dimension(grid%n1:grid%n2)  :: d_grid_vec_partial
    integer, dimension(grid%n1:grid%n2)  :: d_grid_vec_partial_int
    integer, dimension(:,:), allocatable :: d_grid_int

    ! Add routine to call stack
    call init_routine( routine_name)

    d_grid_vec_partial( grid%n1: grid%n2) = field%d_nih( grid%n1: grid%n2)

    where (d_grid_vec_partial)
      d_grid_vec_partial_int = 1
    elsewhere
      d_grid_vec_partial_int = 0
    end where

    if (par%primary) then
      allocate( d_grid_int( grid%nx, grid%ny))
    else
      allocate( d_grid_int( 0,0))
    end if
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial_int, d_grid_int)

    call write_var_primary( filename, ncid, id_var, d_grid_int)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_grid_logical_2D

  subroutine write_to_netcdf_grid_int_2D( field, grid, filename, ncid, id_var)

    ! In/output variables:
    class(type_field_int_2D), intent(in) :: field
    type(type_grid),          intent(in) :: grid
    character(len=*),         intent(in) :: filename
    integer,                  intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'write_to_netcdf_grid_int_2D'
    integer, dimension(grid%n1:grid%n2)  :: d_grid_vec_partial
    integer, dimension(:,:), allocatable :: d_grid

    ! Add routine to call stack
    call init_routine( routine_name)

    d_grid_vec_partial( grid%n1: grid%n2) = field%d_nih( grid%n1: grid%n2)

    if (par%primary) then
      allocate( d_grid( grid%nx, grid%ny))
    else
      allocate( d_grid( 0,0))
    end if
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    call write_var_primary( filename, ncid, id_var, d_grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_grid_int_2D

  subroutine write_to_netcdf_grid_dp_2D( field, grid, filename, ncid, id_var)

    ! In/output variables:
    class(type_field_dp_2D), intent(in) :: field
    type(type_grid),         intent(in) :: grid
    character(len=*),        intent(in) :: filename
    integer,                 intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_netcdf_grid_dp_2D'
    real(dp), dimension(grid%n1:grid%n2)  :: d_grid_vec_partial
    real(dp), dimension(:,:), allocatable :: d_grid

    ! Add routine to call stack
    call init_routine( routine_name)

    d_grid_vec_partial( grid%n1: grid%n2) = field%d_nih( grid%n1: grid%n2)

    if (par%primary) then
      allocate( d_grid( grid%nx, grid%ny))
    else
      allocate( d_grid( 0,0))
    end if
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    call write_var_primary( filename, ncid, id_var, d_grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_grid_dp_2D

  subroutine write_to_netcdf_grid_logical_3D( field, grid, filename, ncid, id_var)

    ! In/output variables:
    class(type_field_logical_3D), intent(in) :: field
    type(type_grid),              intent(in) :: grid
    character(len=*),             intent(in) :: filename
    integer,                      intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_netcdf_grid_logical_3D'
    logical, dimension(grid%n1:grid%n2,1:field%third_dimension_val%n)  :: d_grid_vec_partial
    integer, dimension(grid%n1:grid%n2,1:field%third_dimension_val%n)  :: d_grid_vec_partial_int
    integer, dimension(:,:,:), allocatable :: d_grid_int

    ! Add routine to call stack
    call init_routine( routine_name)

    d_grid_vec_partial( grid%n1: grid%n2, 1:field%third_dimension_val%n) = &
      field%d_nih( grid%n1: grid%n2, 1:field%third_dimension_val%n)

    where (d_grid_vec_partial)
      d_grid_vec_partial_int = 1
    elsewhere
      d_grid_vec_partial_int = 0
    end where

    if (par%primary) then
      allocate( d_grid_int( grid%nx, grid%ny, field%third_dimension_val%n))
    else
      allocate( d_grid_int( 0,0,0))
    end if
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial_int, d_grid_int)

    call write_var_primary( filename, ncid, id_var, d_grid_int)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_grid_logical_3D

  subroutine write_to_netcdf_grid_int_3D( field, grid, filename, ncid, id_var)

    ! In/output variables:
    class(type_field_int_3D), intent(in) :: field
    type(type_grid),          intent(in) :: grid
    character(len=*),         intent(in) :: filename
    integer,                  intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_netcdf_grid_int_3D'
    integer, dimension(grid%n1:grid%n2,1:field%third_dimension_val%n)  :: d_grid_vec_partial
    integer, dimension(:,:,:), allocatable :: d_grid

    ! Add routine to call stack
    call init_routine( routine_name)

    d_grid_vec_partial( grid%n1: grid%n2, 1:field%third_dimension_val%n) = &
      field%d_nih( grid%n1: grid%n2, 1:field%third_dimension_val%n)

    if (par%primary) then
      allocate( d_grid( grid%nx, grid%ny, field%third_dimension_val%n))
    else
      allocate( d_grid( 0,0,0))
    end if
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    call write_var_primary( filename, ncid, id_var, d_grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_grid_int_3D

  subroutine write_to_netcdf_grid_dp_3D( field, grid, filename, ncid, id_var)

    ! In/output variables:
    class(type_field_dp_3D), intent(in) :: field
    type(type_grid),         intent(in) :: grid
    character(len=*),        intent(in) :: filename
    integer,                 intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_netcdf_grid_dp_3D'
    real(dp), dimension(grid%n1:grid%n2,1:field%third_dimension_val%n)  :: d_grid_vec_partial
    real(dp), dimension(:,:,:), allocatable :: d_grid

    ! Add routine to call stack
    call init_routine( routine_name)

    d_grid_vec_partial( grid%n1: grid%n2, 1:field%third_dimension_val%n) = &
      field%d_nih( grid%n1: grid%n2, 1:field%third_dimension_val%n)

    if (par%primary) then
      allocate( d_grid( grid%nx, grid%ny, field%third_dimension_val%n))
    else
      allocate( d_grid( 0,0,0))
    end if
    call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)

    call write_var_primary( filename, ncid, id_var, d_grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_grid_dp_3D

  ! ===== Mesh-based fields
  ! =======================

  subroutine write_to_netcdf_mesh( field, mesh, filename, ncid, id_var)

    ! In/output variables:
    class(atype_field), intent(in) :: field
    type(type_mesh),    intent(in) :: mesh
    character(len=*),   intent(in) :: filename
    integer,            intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_netcdf_mesh'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f => field)
    class default
      call crash('invalid field type')
    class is (type_field_logical_2D)
      call write_to_netcdf_mesh_logical_2D( f, mesh, filename, ncid, id_var)
    class is (type_field_int_2D)
      call write_to_netcdf_mesh_int_2D    ( f, mesh, filename, ncid, id_var)
    class is (type_field_dp_2D)
      call write_to_netcdf_mesh_dp_2D     ( f, mesh, filename, ncid, id_var)
    class is (type_field_logical_3D)
      call write_to_netcdf_mesh_logical_3D( f, mesh, filename, ncid, id_var)
    class is (type_field_int_3D)
      call write_to_netcdf_mesh_int_3D    ( f, mesh, filename, ncid, id_var)
    class is (type_field_dp_3D)
      call write_to_netcdf_mesh_dp_3D     ( f, mesh, filename, ncid, id_var)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_mesh

  subroutine write_to_netcdf_mesh_logical_2D( field, mesh, filename, ncid, id_var)

    ! In/output variables:
    class(type_field_logical_2D), intent(in) :: field
    type(type_mesh),              intent(in) :: mesh
    character(len=*),             intent(in) :: filename
    integer,                      intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'write_to_netcdf_mesh_logical_2D'
    logical, dimension(:), allocatable :: d_tot
    integer, dimension(:), allocatable :: d_tot_int

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) then
      allocate( d_tot    ( field%pai_val%n))
      allocate( d_tot_int( field%pai_val%n))
      call gather_dist_shared_to_primary( field%pai_val, field%d_nih, d_tot=d_tot)
      where (d_tot)
        d_tot_int = 1
      elsewhere
        d_tot_int = 0
      end where
    else
      allocate( d_tot( 0))
      call gather_dist_shared_to_primary( field%pai_val, field%d_nih)
    end if

    call write_var_primary( filename, ncid, id_var, d_tot_int)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_mesh_logical_2D

  subroutine write_to_netcdf_mesh_int_2D( field, mesh, filename, ncid, id_var)

    ! In/output variables:
    class(type_field_int_2D), intent(in) :: field
    type(type_mesh),          intent(in) :: mesh
    character(len=*),         intent(in) :: filename
    integer,                  intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'write_to_netcdf_mesh_int_2D'
    integer, dimension(:), allocatable :: d_tot

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) then
      allocate( d_tot( field%pai_val%n))
      call gather_dist_shared_to_primary( field%pai_val, field%d_nih, d_tot=d_tot)
    else
      allocate( d_tot( 0))
      call gather_dist_shared_to_primary( field%pai_val, field%d_nih)
    end if

    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_mesh_int_2D

  subroutine write_to_netcdf_mesh_dp_2D( field, mesh, filename, ncid, id_var)

    ! In/output variables:
    class(type_field_dp_2D), intent(in) :: field
    type(type_mesh),         intent(in) :: mesh
    character(len=*),        intent(in) :: filename
    integer,                 intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'write_to_netcdf_mesh_dp_2D'
    real(dp), dimension(:), allocatable :: d_tot

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) then
      allocate( d_tot( field%pai_val%n))
      call gather_dist_shared_to_primary( field%pai_val, field%d_nih, d_tot=d_tot)
    else
      allocate( d_tot( 0))
      call gather_dist_shared_to_primary( field%pai_val, field%d_nih)
    end if

    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_mesh_dp_2D

  subroutine write_to_netcdf_mesh_logical_3D( field, mesh, filename, ncid, id_var)

    ! In/output variables:
    class(type_field_logical_3D), intent(in) :: field
    type(type_mesh),              intent(in) :: mesh
    character(len=*),             intent(in) :: filename
    integer,                      intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'write_to_netcdf_mesh_logical_3D'
    logical, dimension(:,:), allocatable :: d_tot
    integer, dimension(:,:), allocatable :: d_tot_int

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) then
      allocate( d_tot    ( field%pai_val%n, field%third_dimension_val%n))
      allocate( d_tot_int( field%pai_val%n, field%third_dimension_val%n))
      call gather_dist_shared_to_primary( field%pai_val, field%third_dimension_val%n, &
        field%d_nih, d_tot=d_tot)
      where (d_tot)
        d_tot_int = 1
      elsewhere
        d_tot_int = 0
      end where
    else
      allocate( d_tot( 0,0))
      call gather_dist_shared_to_primary( field%pai_val, field%third_dimension_val%n, &
        field%d_nih)
    end if

    call write_var_primary( filename, ncid, id_var, d_tot_int)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_mesh_logical_3D

  subroutine write_to_netcdf_mesh_int_3D( field, mesh, filename, ncid, id_var)

    ! In/output variables:
    class(type_field_int_3D), intent(in) :: field
    type(type_mesh),          intent(in) :: mesh
    character(len=*),         intent(in) :: filename
    integer,                  intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'write_to_netcdf_mesh_int_3D'
    integer, dimension(:,:), allocatable :: d_tot

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) then
      allocate( d_tot( field%pai_val%n, field%third_dimension_val%n))
      call gather_dist_shared_to_primary( field%pai_val, field%third_dimension_val%n, &
        field%d_nih, d_tot=d_tot)
    else
      allocate( d_tot( 0,0))
      call gather_dist_shared_to_primary( field%pai_val, field%third_dimension_val%n, &
        field%d_nih)
    end if

    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_mesh_int_3D

  subroutine write_to_netcdf_mesh_dp_3D( field, mesh, filename, ncid, id_var)

    ! In/output variables:
    class(type_field_dp_3D), intent(in) :: field
    type(type_mesh),         intent(in) :: mesh
    character(len=*),        intent(in) :: filename
    integer,                 intent(in) :: ncid, id_var

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'write_to_netcdf_mesh_dp_3D'
    real(dp), dimension(:,:), allocatable :: d_tot

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) then
      allocate( d_tot( field%pai_val%n, field%third_dimension_val%n))
      call gather_dist_shared_to_primary( field%pai_val, field%third_dimension_val%n, &
        field%d_nih, d_tot=d_tot)
    else
      allocate( d_tot( 0,0))
      call gather_dist_shared_to_primary( field%pai_val, field%third_dimension_val%n, &
        field%d_nih)
    end if

    call write_var_primary( filename, ncid, id_var, d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf_mesh_dp_3D

  ! ===== Utilities
  ! ===============

  function get_var_type( field) result( var_type)

    ! In/output variables:
    class(atype_field), intent(in) :: field
    integer                        :: var_type

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'get_var_type'

    ! Add routine to call stack
    call init_routine( routine_name)

    var_type = -42

    select type (f => field)
    class default
      call crash('invalid field type')
    class is (type_field_logical_2D)
      var_type = NF90_INT
    class is (type_field_logical_3D)
      var_type = NF90_INT
    class is (type_field_int_2D)
      var_type = NF90_INT
    class is (type_field_int_3D)
      var_type = NF90_INT
    class is (type_field_dp_2D)
      var_type = NF90_DOUBLE
    class is (type_field_dp_3D)
      var_type = NF90_DOUBLE
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function get_var_type

  function get_dim_ids( field, filename, ncid) result( dim_ids)

    ! In/output variables:
    class(atype_field),     intent(in) :: field
    character(len=*),       intent(in) :: filename
    integer,                intent(in) :: ncid
    integer, dimension(:), allocatable :: dim_ids

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'get_dim_ids'
    integer, dimension(:), allocatable :: dim_ids_grid

    ! Add routine to call stack
    call init_routine( routine_name)

    dim_ids_grid = get_dim_ids_grid( field, filename, ncid)

    select type (f => field)
    class default
      call crash('invalid field type')
    class is (atype_field_2D)
      dim_ids = dim_ids_grid
    class is (atype_field_3D)
      dim_ids = [dim_ids_grid, get_dim_id_third_dimension( f, filename, ncid)]
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function get_dim_ids

  function get_dim_ids_grid( field, filename, ncid) result( dim_ids)

    ! In/output variables:
    class(atype_field),     intent(in) :: field
    character(len=*),       intent(in) :: filename
    integer,                intent(in) :: ncid
    integer, dimension(:), allocatable :: dim_ids

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'get_dim_ids_grid'
    integer                        :: id_dim_vi, id_dim_ti, id_dim_ei, id_dim_x, id_dim_y

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (g => field%grid())
    class default
      call crash('invalid field%grid type')

    class is (type_grid)

      call inquire_dim_multopt( filename, ncid, field_name_options_x, id_dim_x)
      call inquire_dim_multopt( filename, ncid, field_name_options_y, id_dim_y)
      if (id_dim_x == -1) call crash('couldnt find x dimension in file ' // trim( filename))
      if (id_dim_y == -1) call crash('couldnt find y dimension in file ' // trim( filename))

      dim_ids = [id_dim_x, id_dim_y]

    class is (type_mesh)

      if (field%Arakawa_grid() == Arakawa_grid%a()) then
        call inquire_dim_multopt( filename, ncid, 'vi', id_dim_vi)
        if (id_dim_vi == -1) call crash('couldnt find vi dimension in file ' // trim( filename))
        dim_ids = [id_dim_vi]
      elseif (field%Arakawa_grid() == Arakawa_grid%b()) then
        call inquire_dim_multopt( filename, ncid, 'ti', id_dim_ti)
        if (id_dim_ti == -1) call crash('couldnt find ti dimension in file ' // trim( filename))
        dim_ids = [id_dim_ti]
      elseif (field%Arakawa_grid() == Arakawa_grid%c()) then
        call inquire_dim_multopt( filename, ncid, 'ei', id_dim_ei)
        if (id_dim_ei == -1) call crash('couldnt find ei dimension in file ' // trim( filename))
        dim_ids = [id_dim_ei]
      else
        call crash('invalid Arakawa grid')
      end if

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function get_dim_ids_grid

  function get_dim_id_third_dimension( field, filename, ncid) result( dim_id)

    ! In/output variables:
    class(atype_field_3D),   intent(in) :: field
    character(len=*),       intent(in) :: filename
    integer,                intent(in) :: ncid
    integer                            :: dim_id

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'get_dim_id_third_dimension'
    integer                        :: dim_length

    ! Add routine to call stack
    call init_routine( routine_name)

    call inquire_dim( filename, ncid, trim( field%third_dimension_val%name), dim_length, dim_id)
    if (dim_id == -1) call crash('couldnt find dimension ' // trim( field%third_dimension_val%name) // &
      ' in file ' // trim( filename))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function get_dim_id_third_dimension

end submodule fields_basic_submod_write_to_netcdf