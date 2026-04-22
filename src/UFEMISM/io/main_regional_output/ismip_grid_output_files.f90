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
  use ismip_output_types, only: type_ismip_grid_output, type_ismip_gridded_field

  implicit none

  private

  public :: create_ISMIP_regional_output_files_grid, write_to_ISMIP_regional_output_files_grid

contains

  subroutine write_to_ISMIP_regional_output_files_grid( region)
    !< Write to ISMIP regional output NetCDF files - grid version

    ! In/output variables:
    type(type_model_region), intent(in   ) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_ISMIP_regional_output_files_grid'

    ! Add routine to path
    call init_routine( routine_name)

    ! if no ISMIP output should be created, do nothing 
    if (.not. C%do_create_ismip_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to ISMIP grid output files' // '...'

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_ISMIP_regional_output_files_grid

  subroutine create_ISMIP_regional_output_files_grid( region)
    ! MAIN creation
    !< Create all ISMIP regional output NetCDF files - grid version

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_ISMIP_regional_output_files_grid'

    ! Add routine to path
    call init_routine( routine_name)

    ! if no ISMIP output should be created, do nothing
    if (.not. C%do_create_ismip_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Initialise ISMIP_grid_output
    call initialise_ismip_grid_output( region)

    ! Create folder
    if (par%primary) call system('mkdir ' // trim(C%output_dir) // trim(C%ismip_folder))

    ! Create all grid files
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%lithk)
    ! TODO add the rest

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_ISMIP_regional_output_files_grid

  subroutine create_single_ISMIP_regional_output_file_grid( ismip_grid_output, field)
    !< Create a single ISMIP regional output NetCDF file - grid version

    ! In/output variables:
    type(type_ismip_grid_output),      intent(inout) :: ismip_grid_output
    type(type_ismip_gridded_field),    intent(in   ) :: field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_single_ISMIP_regional_output_file_grid'
    integer                        :: ncid
    character(len=*), parameter    :: iprecision = 'single'
    logical, parameter             :: do_compress = .true.

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating ISMIP grid output file "' // &
      UPSY%stru%colour_string( trim( field%filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( field%filename, ncid)

    ! Set up the grid in the file
    call setup_xy_grid_in_netcdf_file( field%filename, ncid, ismip_grid_output%grid)

    ! Add time dimension to the file
    ! TODO change to cftime
    call add_time_dimension_to_file( field%filename, ncid)

    ! TODO if fieldtype == 'FL', add_bnds_dimension_to_file

    ! Add the field
    call add_field_grid_dp_2D( field%filename, ncid, field%name, precision = iprecision, & 
      do_compress = do_compress, long_name = field%long_name, units = field%units)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_single_ISMIP_regional_output_file_grid

  subroutine initialise_ISMIP_grid_output( region)
    ! Initialise the ISMIP grid type

    ! NOTE: should only be called by the primary (?)

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_ISMIP_grid_output'

    ! Add routine to path
    call init_routine( routine_name)

    ! Convert region name to ISMIP ice sheet name
    select case (region%name)
      case default
        call crash('invalid region name for ISMIP output "' // trim( region%name) // '"')
      case ('ANT')
        region%ismip_grid_output%IS_name = 'AIS'
      case ('GRL')
        region%ismip_grid_output%IS_name = 'GIS'
    end select

    ! Copy the output grid
    region%ismip_grid_output%grid = region%output_grid

    call initialise_ISMIP_field( region%ismip_grid_output, region%ismip_grid_output%lithk, 'lithk', 'land_ice_thickness', 'm', 'ST')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ISMIP_grid_output

  subroutine initialise_ISMIP_field( ismip_grid_output, field, name, long_name, units, fieldtype)
    ! Initialise a single field

    type(type_ismip_grid_output),      intent(inout) :: ismip_grid_output
    type(type_ismip_gridded_field),    intent(inout) :: field
    character(len=*),               intent(in   ) :: name
    character(len=*),               intent(in   ) :: long_name
    character(len=*),               intent(in   ) :: units
    character(len=2),                  intent(in   ) :: fieldtype

    ! Local variables
    character(len=1024), parameter :: routine_name = 'initialise_ISMIP_field'

    ! Add routine to path
    call init_routine( routine_name)

    ! Inherit metadata
    field%name      = name
    field%long_name = long_name
    field%units     = units
    field%fieldtype = fieldtype

    ! Define the filename for this field
    field%filename = trim( C%output_dir) // trim(C%ismip_folder) // '/' // trim(field%name) // '_' // &
      trim(ismip_grid_output%IS_name) // '_' // trim(C%ismip_group_name) // '_' // trim(C%ismip_model_name) // & 
      '_' // trim(C%ismip_exp_name) // '.nc'

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Defined name "' // &
      UPSY%stru%colour_string( trim( name), 'light blue') // '"...'

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ISMIP_field

end module ismip_grid_output_files
