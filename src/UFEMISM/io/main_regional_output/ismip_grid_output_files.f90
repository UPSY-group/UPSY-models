module ismip_grid_output_files

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
  use UPSY_main, only: UPSY
  use precisions, only: dp 
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use model_configuration, only: C
  use region_types, only: type_model_region
  use grid_types, only: type_grid
  use netcdf_io_main
  use ice_mass_and_fluxes, only: calc_ice_margin_fluxes
  use remapping_main, only: map_from_mesh_vertices_to_xy_grid_2D, &
    map_from_mesh_vertices_to_xy_grid_3D, map_from_mesh_vertices_to_xy_grid_2D_minval, &
    map_from_mesh_triangles_to_xy_grid_2D, map_from_mesh_triangles_to_xy_grid_3D
  use mpi_distributed_memory, only: gather_to_all

  implicit none

  private

  public :: create_ISMIP_regional_output_file_grid, write_to_ISMIP_regional_output_file_grid

contains

  subroutine write_to_ISMIP_regional_output_files_grid( region)
    !< Write to ISMIP regional output NetCDF files - grid version

    ! In/output variables:
    type(type_model_region), intent(in   ) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_ISMIP_regional_output_files_grid'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no ISMIP output should be created, do nothing 
    if (.not. C%do_create_ismip_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to ISMIP grid output files' // '...'

    ! Write to all files
    call write_to_single_ISMIP_regional_output_file_grid( region, 'lithk')
    ! TODO add the rest

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_ISMIP_regional_output_files_grid

  subroutine write_to_single_ISMIP_regional_output_file_grid( region, fieldname)
    !< Write to a single ISMIP regional output NetCDF file - grid version

    ! In/output variables:
    type(type_model_region), intent(in   ) :: region
    character(len=*),        intent(in   ) :: fieldname

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_single_ISMIP_regional_output_file_grid'
    integer                        :: ncid
    character(len=*)               :: filename

    ! Add routine to path
    call init_routine( routine_name)

    ! Get filename
    filename = get_ISMIP_filename( region, fieldname)    

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( region%output_filename_grid, ncid)

    ! TODO get fieldtype and write proper cftime with/without bounds
    ! write the time to the file
    call write_time_to_file( region%output_filename_grid, ncid, region%time)

    call write_to_ISMIP_regional_output_file_grid_field( region, grid, filename, fieldname)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_single_ISMIP_regional_output_file_grid

  subroutine write_to_ISMIP_regional_output_file_grid_field( region, grid, filename, ncid, choice_output_field)
    !< Write a single field to the main regional output NetCDF file - grid version

    ! In/output variables:
    type(type_model_region), intent(in   ) :: region
    type(type_grid),         intent(in   ) :: grid
    character(len=*),        intent(in   ) :: filename
    integer,                 intent(in   ) :: ncid
    character(len=*),        intent(in   ) :: choice_output_field

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_ISMIP_regional_output_file_grid_field'
    real(dp), dimension(:),   allocatable :: d_mesh_vec_partial_2D
    real(dp), dimension(:),   allocatable :: d_grid_vec_partial_2D
    !real(dp), dimension(:),   allocatable :: mask_int
    !integer                               :: vi, ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    !allocate( d_mesh_vec_partial_2D(    region%mesh%vi1:region%mesh%vi2))
    allocate( d_grid_vec_partial_2D(         grid%n_loc                ))
    !allocate( mask_int( region%mesh%vi1:region%mesh%vi2), source = 0._dp)

    ! Add the specified data field to the file
    select case (choice_output_field)
      case default
        call crash('unknown choice_output_field "' // trim( choice_output_field) // '"')
      case ('none')
        ! Do nothing
      case ('lithk')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, grid, C%output_dir, region%ice%Hi, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( grid, filename, ncid, 'lithk', d_grid_vec_partial_2D)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_ISMIP_regional_output_file_grid_field

  subroutine create_ISMIP_regional_output_file_grids( region)
    !< Create all ISMIP regional output NetCDF files - grid version

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_ISMIP_regional_output_file_grids'

    ! Add routine to path
    call init_routine( routine_name)

    ! if no ISMIP output should be created, do nothing
    if (.not. C%do_create_ismip_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Create all grid files
    call create_single_ISMIP_regional_output_file_grid( region, 'lithk')
    ! TODO add the rest

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_ISMIP_regional_output_file_grids

  subroutine create_single_ISMIP_regional_output_file_grid( region, fieldname)
    !< Create a single ISMIP regional output NetCDF file - grid version

    ! In/output variables:
    type(type_model_region), intent(inout) :: region
    character(len=*),        intent(in   ) :: fieldname

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_single_ISMIP_regional_output_file_grid'
    integer                        :: ncid
    character(len=*),              :: filename
    character(len=*),              :: long_name
    character(len=*),              :: units
    character(len=2),              :: fieldtype        ! 'ST' (state) or 'FL' (flux)
    character(len=*), parameter    :: iprecision = 'single'
    logical, parameter             :: do_compress = .true.

    ! Add routine to path
    call init_routine( routine_name)

    ! Get filename of this field
    filename = get_ISMIP_filename( region, fieldname)    

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating grid output file "' // &
      UPSY%stru%colour_string( trim( filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( filename, ncid)

    ! Set up the grid in the file
    call setup_xy_grid_in_netcdf_file( filename, ncid, region%output_grid)

    ! Add time dimension to the file
    call add_time_dimension_to_file( filename, ncid)

    ! Get fieldtype
    fieldtype = get_ISMIP_fieldtype( fieldname)

    ! TODO if fieldtype == 'FL', add_bnds_dimension_to_file

    ! Gather metadata
    call get_ISMIP_metadata( fieldname, long_name, units)

    ! Add the field
    call add_field_grid_dp_2D( filename, ncid, fieldname, precision = iprecision, do_compress = do_compress, long_name = long_name, units = units)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_single_ISMIP_regional_output_file_grid

  subroutine get_ISMIP_metadata( fieldname, long_name, units)
    !< Get the metadata of an ISMIP field

    ! In/output variables:
    type(type_model_region), intent(inout) :: region
    character(len=*),        intent(in   ) :: fieldname

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'get_ISMIP_metadata'
    character(len=*),              :: long_name
    character(len=*),              :: units

    ! Add routine to path
    call init_routine( routine_name)

    select case (fieldname)
      case default
        call crash('invalid ISMIP fieldname "' // trim( fieldname) // '"')
      case ('lithk')
        long_name = 'land_ice_thickness'
        units = 'm'
      ! TODO add the rest
    end select

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_single_ISMIP_regional_output_file_grid

  function get_ISMIP_filename( region, fieldname) result(filename)
    ! Define the filename for this specific field

    type(type_model_region)     :: region
    character(len=*)            :: fieldname
    character(len=*)            :: filename

    ! Local variables
    character(len=3)            :: is_name

    ! Convert region name to ISMIP ice sheet name
    select case (region%name)
      case default
        call crash('invalid region name for ISMIP output "' // trim( region%name) // '"')
      case ('ANT')
        is_name = 'AIS'
      case ('GRL')
        is_name = 'GIS'
    end select

    filename = trim( C%output_dir) // trim(C%ismip_folder) // '/' // trim(fieldname) // '_' // trim(is_name) // &
      '_' // trim(C%ismip_group_name) // '_' // trim(C%ismip_model_name) // '_' // trim(C%ismip_experiment_name) // '.nc'

  end function get_ISMIP_filename

  function get_ISMIP_fieldtype( fieldname) result (fieldtype)
    ! Determine whether this is a ST (state) or FL (flux) variable

    character(len=*)           :: fieldname
    character(len=*)           :: fieldtype

    select case (fieldname)
      case default
        call crash('unknown fieldname for ISMIP output "' // trim(fieldname) // '"')
      case ('lithk')
        fieldtype = 'ST'
      ! TODO add the rest
    end select

  end function get_ISMIP_fieldtype

end module ismip_grid_output_files
