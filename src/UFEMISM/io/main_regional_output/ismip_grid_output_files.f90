module ismip_grid_output_files

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
  use UPSY_main, only: UPSY
  use parameters
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

    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%lithk)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%orog)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%topg)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%base)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%xvelsurf)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%yvelsurf)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%zvelsurf)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%xvelbase)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%yvelbase)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%zvelbase)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%xvelmean)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%yvelmean)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%strbasemag)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_ISMIP_regional_output_files_grid

  subroutine write_to_single_ISMIP_regional_output_file_grid( region, field)
    !< Write to single ISMIP regional output NetCDF file - grid version

    ! In/output variables:
    type(type_model_region),           intent(in   ) :: region
    type(type_ismip_gridded_field),    intent(in   ) :: field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_single_ISMIP_regional_output_file_grid'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( field%filename, ncid)

    ! write the time to the file
    call write_time_to_file( field%filename, ncid, region%time)

    ! write the data to the file
    call write_to_ISMIP_regional_output_file_grid_field( region, field, ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_single_ISMIP_regional_output_file_grid

  subroutine write_to_ISMIP_regional_output_file_grid_field( region, field, ncid)
    !< Write a single field to the ISMIP regional output NetCDF file - grid version

    ! In/output variables:
    type(type_model_region),           intent(in   ) :: region
    type(type_ismip_gridded_field),    intent(in   ) :: field
    integer,                           intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_ISMIP_regional_output_file_grid_field'
    real(dp), dimension(:),   allocatable :: d_grid_vec_partial_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate memory
    allocate( d_grid_vec_partial_2D( region%output_grid%n_loc ))

    ! Add the specified data field to the file
    select case (field%name)
      case default
        call crash('unknown choice_output_field "' // trim( field%name) // '"')
      case ('lithk')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%Hi, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
      case ('orog')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%Hs, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
      case ('topg')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%Hb, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
      case ('base')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%Hib, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
      case ('xvelsurf')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%u_surf, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('yvelsurf')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%v_surf, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('zvelsurf')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%w_surf, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('xvelbase')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%u_base, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('yvelbase')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%v_base, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('zvelbase')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%w_base, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('xvelmean')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%u_vav, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('yvelmean')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%v_vav, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('strbasemag')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%basal_shear_stress, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_ISMIP_regional_output_file_grid_field

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
    if (par%primary) call system('mkdir ' // trim(C%output_dir) // trim(C%ismip_exp_name))

    ! Create all grid files
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%lithk)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%orog)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%topg)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%base)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%xvelsurf)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%yvelsurf)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%zvelsurf)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%xvelbase)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%yvelbase)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%zvelbase)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%xvelmean)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%yvelmean)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%strbasemag)
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
    call initialise_ISMIP_field( region%ismip_grid_output, region%ismip_grid_output%orog,  'orog' , 'surface_altitude', 'm', 'ST')
    call initialise_ISMIP_field( region%ismip_grid_output, region%ismip_grid_output%topg,  'topg' , 'bedrock_altitude', 'm', 'ST')
    call initialise_ISMIP_field( region%ismip_grid_output, region%ismip_grid_output%base,  'base' , 'base_altitude', 'm', 'ST')
    call initialise_ISMIP_field( region%ismip_grid_output, region%ismip_grid_output%xvelsurf, 'xvelsurf' , 'land_ice_surface_x_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region%ismip_grid_output, region%ismip_grid_output%yvelsurf, 'yvelsurf' , 'land_ice_surface_y_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region%ismip_grid_output, region%ismip_grid_output%zvelsurf, 'zvelsurf' , 'land_ice_surface_upward_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region%ismip_grid_output, region%ismip_grid_output%xvelbase, 'xvelbase' , 'land_ice_basal_x_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region%ismip_grid_output, region%ismip_grid_output%yvelbase, 'yvelbase' , 'land_ice_basal_y_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region%ismip_grid_output, region%ismip_grid_output%zvelbase, 'zvelbase' , 'land_ice_surface_upward_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region%ismip_grid_output, region%ismip_grid_output%xvelmean, 'xvelmean' , 'land_ice_vertical_mean_x_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region%ismip_grid_output, region%ismip_grid_output%yvelmean, 'yvelmean' , 'land_ice_vertical_mean_y_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region%ismip_grid_output, region%ismip_grid_output%strbasemag, 'strbasemag' , 'land_ice_basal_drag', 'Pa', 'ST')

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
    field%filename = trim( C%output_dir) // trim(C%ismip_exp_name) // '/' // trim(field%name) // '_' // &
      trim(ismip_grid_output%IS_name) // '_' // trim(C%ismip_group_name) // '_' // trim(C%ismip_model_name) // & 
      '_' // trim(C%ismip_exp_name) // '.nc'

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Defined name "' // &
      UPSY%stru%colour_string( trim( name), 'light blue') // '"...'

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ISMIP_field

end module ismip_grid_output_files
