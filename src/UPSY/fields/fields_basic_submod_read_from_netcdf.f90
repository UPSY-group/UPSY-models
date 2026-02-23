submodule (fields_basic) fields_basic_submod_read_from_netcdf

  use netcdf_io_main
  use mpi_distributed_memory_grid, only: distribute_gridded_data_from_primary
  use mpi_distributed_memory, only: distribute_from_primary

contains

  subroutine read_from_netcdf( self, filename, ncid)
    ! NOTE: assumes the NetCDF file is defined on
    !       the same grid/mesh as the field

    ! In/output variables:
    class(atype_field), intent(inout) :: self
    character(len=*),   intent(in   ) :: filename
    integer,            intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'read_from_netcdf'
    integer                        :: id_var

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Check if a variable for this field exists in the netcdf file
    call inquire_var( filename, ncid, self%name(), id_var)
    if (id_var == -1) call crash('variable "' // trim( self%name()) // &
      '" not found in file "' // trim( filename) // '"')

    select type (f => self)
    class default
      call crash('invalid field type')
    class is (atype_field_2D)
    class is (atype_field_3D)
      if (.not. f%third_dimension_val%is_dim_and_var_in_netcdf( filename, ncid)) then
        call crash('third dimension of field "' // trim( self%name()) // &
          '" not found in file "' // trim( filename) // '"')
      end if
    end select

    select type (g => self%grid())
    class default
      call crash('invalid field%grid type')
    class is (type_grid)
      call read_from_netcdf_grid( self, g, filename, ncid)
    class is (type_mesh)
      call read_from_netcdf_mesh( self, filename, ncid)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf

  ! ===== Grid-based fields
  ! =======================

  subroutine read_from_netcdf_grid( field, grid, filename, ncid)

    ! In/output variables:
    class(atype_field), intent(inout) :: field
    type(type_grid),    intent(in   ) :: grid
    character(len=*),   intent(in   ) :: filename
    integer,            intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'read_from_netcdf_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f => field)
    class default
      call crash('invalid field type')
    class is (type_field_logical_2D)
      call read_from_netcdf_grid_logical_2D( f, grid, filename, ncid)
    class is (type_field_logical_3D)
      call read_from_netcdf_grid_logical_3D( f, grid, filename, ncid)
    class is (type_field_int_2D)
      call read_from_netcdf_grid_int_2D( f, grid, filename, ncid)
    class is (type_field_int_3D)
      call read_from_netcdf_grid_int_3D( f, grid, filename, ncid)
    class is (type_field_dp_2D)
      call read_from_netcdf_grid_dp_2D( f, grid, filename, ncid)
    class is (type_field_dp_3D)
      call read_from_netcdf_grid_dp_3D( f, grid, filename, ncid)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_grid

  subroutine read_from_netcdf_grid_logical_2D( field, grid, filename, ncid)

    ! In/output variables:
    type(type_field_logical_2D), intent(inout) :: field
    type(type_grid),             intent(in   ) :: grid
    character(len=*),            intent(in   ) :: filename
    integer,                     intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'read_from_netcdf_grid_logical_2D'
    type(type_par_arr_info)              :: pai
    integer                              :: id_var
    logical, dimension(:,:), allocatable :: d_grid_tot
    logical, dimension(:), pointer       :: d_grid_vec_partial

    ! Add routine to call stack
    call init_routine( routine_name)

    pai = field%pai()
    d_grid_vec_partial => field%d_nih( pai%i1: pai%i2)

    call inquire_var( filename, ncid, field%name(), id_var)

    if (par%primary) then
      allocate( d_grid_tot( grid%nx, grid%ny))
    else
      allocate( d_grid_tot( 0,0))
    end if

    call read_var_primary( filename, ncid, id_var, d_grid_tot)
    call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_grid_logical_2D

  subroutine read_from_netcdf_grid_logical_3D( field, grid, filename, ncid)

    ! In/output variables:
    type(type_field_logical_3D), intent(inout) :: field
    type(type_grid),             intent(in   ) :: grid
    character(len=*),            intent(in   ) :: filename
    integer,                     intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'read_from_netcdf_grid_logical_3D'
    type(type_par_arr_info)                :: pai
    integer                                :: id_var
    type(type_third_dimension)             :: field_third_dimension
    logical, dimension(:,:,:), allocatable :: d_grid_tot
    logical, dimension(:,:), pointer       :: d_grid_vec_partial

    ! Add routine to call stack
    call init_routine( routine_name)

    pai = field%pai()
    field_third_dimension = field%third_dimension()
    d_grid_vec_partial => field%d_nih( pai%i1: pai%i2, 1: field_third_dimension%n)

    call inquire_var( filename, ncid, field%name(), id_var)

    if (par%primary) then
      allocate( d_grid_tot( grid%nx, grid%ny, field_third_dimension%n))
    else
      allocate( d_grid_tot( 0,0,0))
    end if

    call read_var_primary( filename, ncid, id_var, d_grid_tot)
    call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_grid_logical_3D

  subroutine read_from_netcdf_grid_int_2D( field, grid, filename, ncid)

    ! In/output variables:
    type(type_field_int_2D), intent(inout) :: field
    type(type_grid),         intent(in   ) :: grid
    character(len=*),        intent(in   ) :: filename
    integer,                 intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'read_from_netcdf_grid_int_2D'
    type(type_par_arr_info)              :: pai
    integer                              :: id_var
    integer, dimension(:,:), allocatable :: d_grid_tot
    integer, dimension(:), pointer       :: d_grid_vec_partial

    ! Add routine to call stack
    call init_routine( routine_name)

    pai = field%pai()
    d_grid_vec_partial => field%d_nih( pai%i1: pai%i2)

    call inquire_var( filename, ncid, field%name(), id_var)

    if (par%primary) then
      allocate( d_grid_tot( grid%nx, grid%ny))
    else
      allocate( d_grid_tot( 0,0))
    end if

    call read_var_primary( filename, ncid, id_var, d_grid_tot)
    call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_grid_int_2D

  subroutine read_from_netcdf_grid_int_3D( field, grid, filename, ncid)

    ! In/output variables:
    type(type_field_int_3D), intent(inout) :: field
    type(type_grid),         intent(in   ) :: grid
    character(len=*),        intent(in   ) :: filename
    integer,                 intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'read_from_netcdf_grid_int_3D'
    type(type_par_arr_info)                :: pai
    integer                                :: id_var
    type(type_third_dimension)             :: field_third_dimension
    integer, dimension(:,:,:), allocatable :: d_grid_tot
    integer, dimension(:,:), pointer       :: d_grid_vec_partial

    ! Add routine to call stack
    call init_routine( routine_name)

    pai = field%pai()
    field_third_dimension = field%third_dimension()
    d_grid_vec_partial => field%d_nih( pai%i1: pai%i2, 1: field_third_dimension%n)

    call inquire_var( filename, ncid, field%name(), id_var)

    if (par%primary) then
      allocate( d_grid_tot( grid%nx, grid%ny, field_third_dimension%n))
    else
      allocate( d_grid_tot( 0,0,0))
    end if

    call read_var_primary( filename, ncid, id_var, d_grid_tot)
    call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_grid_int_3D

  subroutine read_from_netcdf_grid_dp_2D( field, grid, filename, ncid)

    ! In/output variables:
    type(type_field_dp_2D), intent(inout) :: field
    type(type_grid),        intent(in   ) :: grid
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_from_netcdf_grid_dp_2D'
    type(type_par_arr_info)               :: pai
    integer                               :: id_var
    real(dp), dimension(:,:), allocatable :: d_grid_tot
    real(dp), dimension(:), pointer       :: d_grid_vec_partial

    ! Add routine to call stack
    call init_routine( routine_name)

    pai = field%pai()
    d_grid_vec_partial => field%d_nih( pai%i1: pai%i2)

    call inquire_var( filename, ncid, field%name(), id_var)

    if (par%primary) then
      allocate( d_grid_tot( grid%nx, grid%ny))
    else
      allocate( d_grid_tot( 0,0))
    end if

    call read_var_primary( filename, ncid, id_var, d_grid_tot)
    call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_grid_dp_2D

  subroutine read_from_netcdf_grid_dp_3D( field, grid, filename, ncid)

    ! In/output variables:
    type(type_field_dp_3D), intent(inout) :: field
    type(type_grid),        intent(in   ) :: grid
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'read_from_netcdf_grid_dp_3D'
    type(type_par_arr_info)                 :: pai
    integer                                 :: id_var
    type(type_third_dimension)              :: field_third_dimension
    real(dp), dimension(:,:,:), allocatable :: d_grid_tot
    real(dp), dimension(:,:), pointer       :: d_grid_vec_partial

    ! Add routine to call stack
    call init_routine( routine_name)

    pai = field%pai()
    field_third_dimension = field%third_dimension()
    d_grid_vec_partial => field%d_nih( pai%i1: pai%i2, 1: field_third_dimension%n)

    call inquire_var( filename, ncid, field%name(), id_var)

    if (par%primary) then
      allocate( d_grid_tot( grid%nx, grid%ny, field_third_dimension%n))
    else
      allocate( d_grid_tot( 0,0,0))
    end if

    call read_var_primary( filename, ncid, id_var, d_grid_tot)
    call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_grid_dp_3D

  ! ===== Mesh-based fields
  ! =======================

  subroutine read_from_netcdf_mesh( field, filename, ncid)

    ! In/output variables:
    class(atype_field), intent(inout) :: field
    character(len=*),   intent(in   ) :: filename
    integer,            intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'read_from_netcdf_mesh'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f => field)
    class default
      call crash('invalid field type')
    class is (type_field_logical_2D)
      call read_from_netcdf_mesh_logical_2D( f, filename, ncid)
    class is (type_field_logical_3D)
      call read_from_netcdf_mesh_logical_3D( f, filename, ncid)
    class is (type_field_int_2D)
      call read_from_netcdf_mesh_int_2D( f, filename, ncid)
    class is (type_field_int_3D)
      call read_from_netcdf_mesh_int_3D( f, filename, ncid)
    class is (type_field_dp_2D)
      call read_from_netcdf_mesh_dp_2D( f, filename, ncid)
    class is (type_field_dp_3D)
      call read_from_netcdf_mesh_dp_3D( f, filename, ncid)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_mesh

  subroutine read_from_netcdf_mesh_logical_2D( field, filename, ncid)

    ! In/output variables:
    type(type_field_logical_2D), intent(inout) :: field
    character(len=*),            intent(in   ) :: filename
    integer,                     intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'read_from_netcdf_mesh_logical_2D'
    type(type_par_arr_info)            :: pai
    integer                            :: id_var
    logical, dimension(:), allocatable :: d_tot
    logical, dimension(:), pointer     :: d_partial

    ! Add routine to call stack
    call init_routine( routine_name)

    pai = field%pai()
    d_partial => field%d_nih( pai%i1: pai%i2)

    call inquire_var( filename, ncid, field%name(), id_var)

    if (par%primary) then
      allocate( d_tot( pai%n))
    else
      allocate( d_tot( 0))
    end if

    call read_var_primary( filename, ncid, id_var, d_tot)
    call distribute_from_primary( d_partial, d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_mesh_logical_2D

  subroutine read_from_netcdf_mesh_logical_3D( field, filename, ncid)

    ! In/output variables:
    type(type_field_logical_3D), intent(inout) :: field
    character(len=*),            intent(in   ) :: filename
    integer,                     intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'read_from_netcdf_mesh_logical_3D'
    type(type_par_arr_info)              :: pai
    integer                              :: id_var
    type(type_third_dimension)           :: field_third_dimension
    logical, dimension(:,:), allocatable :: d_tot
    logical, dimension(:,:), pointer     :: d_partial

    ! Add routine to call stack
    call init_routine( routine_name)

    pai = field%pai()
    field_third_dimension = field%third_dimension()
    d_partial => field%d_nih( pai%i1: pai%i2, 1: field_third_dimension%n)

    call inquire_var( filename, ncid, field%name(), id_var)

    if (par%primary) then
      allocate( d_tot( pai%n, field_third_dimension%n))
    else
      allocate( d_tot( 0,0))
    end if

    call read_var_primary( filename, ncid, id_var, d_tot)
    call distribute_from_primary( d_partial, d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_mesh_logical_3D

  subroutine read_from_netcdf_mesh_int_2D( field, filename, ncid)

    ! In/output variables:
    type(type_field_int_2D), intent(inout) :: field
    character(len=*),        intent(in   ) :: filename
    integer,                 intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'read_from_netcdf_mesh_int_2D'
    type(type_par_arr_info)            :: pai
    integer                            :: id_var
    integer, dimension(:), allocatable :: d_tot
    integer, dimension(:), pointer     :: d_partial

    ! Add routine to call stack
    call init_routine( routine_name)

    pai = field%pai()
    d_partial => field%d_nih( pai%i1: pai%i2)

    call inquire_var( filename, ncid, field%name(), id_var)

    if (par%primary) then
      allocate( d_tot( pai%n))
    else
      allocate( d_tot( 0))
    end if

    call read_var_primary( filename, ncid, id_var, d_tot)
    call distribute_from_primary( d_partial, d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_mesh_int_2D

  subroutine read_from_netcdf_mesh_int_3D( field, filename, ncid)

    ! In/output variables:
    type(type_field_int_3D), intent(inout) :: field
    character(len=*),        intent(in   ) :: filename
    integer,                 intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'read_from_netcdf_mesh_int_3D'
    type(type_par_arr_info)              :: pai
    integer                              :: id_var
    type(type_third_dimension)           :: field_third_dimension
    integer, dimension(:,:), allocatable :: d_tot
    integer, dimension(:,:), pointer     :: d_partial

    ! Add routine to call stack
    call init_routine( routine_name)

    pai = field%pai()
    field_third_dimension = field%third_dimension()
    d_partial => field%d_nih( pai%i1: pai%i2, 1: field_third_dimension%n)

    call inquire_var( filename, ncid, field%name(), id_var)

    if (par%primary) then
      allocate( d_tot( pai%n, field_third_dimension%n))
    else
      allocate( d_tot( 0,0))
    end if

    call read_var_primary( filename, ncid, id_var, d_tot)
    call distribute_from_primary( d_partial, d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_mesh_int_3D

  subroutine read_from_netcdf_mesh_dp_2D( field, filename, ncid)

    ! In/output variables:
    type(type_field_dp_2D), intent(inout) :: field
    character(len=*),        intent(in   ) :: filename
    integer,                 intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'read_from_netcdf_mesh_dp_2D'
    type(type_par_arr_info)             :: pai
    integer                             :: id_var
    real(dp), dimension(:), allocatable :: d_tot
    real(dp), dimension(:), pointer     :: d_partial

    ! Add routine to call stack
    call init_routine( routine_name)

    pai = field%pai()
    d_partial => field%d_nih( pai%i1: pai%i2)

    call inquire_var( filename, ncid, field%name(), id_var)

    if (par%primary) then
      allocate( d_tot( pai%n))
    else
      allocate( d_tot( 0))
    end if

    call read_var_primary( filename, ncid, id_var, d_tot)
    call distribute_from_primary( d_partial, d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_mesh_dp_2D

  subroutine read_from_netcdf_mesh_dp_3D( field, filename, ncid)

    ! In/output variables:
    type(type_field_dp_3D), intent(inout) :: field
    character(len=*),        intent(in   ) :: filename
    integer,                 intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'read_from_netcdf_mesh_dp_3D'
    type(type_par_arr_info)               :: pai
    integer                               :: id_var
    type(type_third_dimension)            :: field_third_dimension
    real(dp), dimension(:,:), allocatable :: d_tot
    real(dp), dimension(:,:), pointer     :: d_partial

    ! Add routine to call stack
    call init_routine( routine_name)

    pai = field%pai()
    field_third_dimension = field%third_dimension()
    d_partial => field%d_nih( pai%i1: pai%i2, 1: field_third_dimension%n)

    call inquire_var( filename, ncid, field%name(), id_var)

    if (par%primary) then
      allocate( d_tot( pai%n, field_third_dimension%n))
    else
      allocate( d_tot( 0,0))
    end if

    call read_var_primary( filename, ncid, id_var, d_tot)
    call distribute_from_primary( d_partial, d_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_mesh_dp_3D

end submodule fields_basic_submod_read_from_netcdf