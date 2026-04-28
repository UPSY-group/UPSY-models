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
    map_from_mesh_triangles_to_xy_grid_2D
  use mpi_distributed_memory, only: gather_to_all
  use ismip_output_types, only: type_ismip_grid_output, type_ismip_gridded_field

  implicit none

  private

  public :: create_ISMIP_regional_output_files_grid, write_to_ISMIP_regional_output_files_grid, &
    accumulate_ISMIP_flux_fields

contains

  subroutine accumulate_ISMIP_flux_fields( region)
    !< Write to ISMIP regional output NetCDF files - grid version

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'accumulate_ISMIP_flux_fields'
    real( dp)                      :: deltat     ! 

    ! Add routine to path
    call init_routine( routine_name)

    ! if no ISMIP output should be created, do nothing 
    if (.not. C%do_create_ismip_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Get delta t since last current time
    deltat = region%time - region%ismip_grid_output%t_curr

    ! Accumulate all FL fields (except non-varying such as geothermal heat flux)
    call accumulate_single_ISMIP_flux_field( region, region%ismip_grid_output%acabf, deltat)

    ! Update current time
    region%ismip_grid_output%t_curr = region%time

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine accumulate_ISMIP_flux_fields

  subroutine accumulate_single_ISMIP_flux_field( region, field, deltat)
    !< Write to ISMIP regional output NetCDF files - grid version

    ! In/output variables:
    type(type_model_region),           intent(inout) :: region
    type(type_ismip_gridded_field),    intent(inout) :: field
    real(dp),                          intent(in   ) :: deltat

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'accumulate_single_ISMIP_flux_field'
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: SMB_loc
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    select case (field%name)
      case default
        call crash('unknown ISMIP field name "' // trim( field%name) // '"')
      case ('acabf')
        ! For hybrid memory reasons, make a local copy (is this necessary?)
        SMB_loc( region%mesh%vi1: region%mesh%vi2) = region%SMB%SMB( region%mesh%vi1: region%mesh%vi2)
        do vi = region%mesh%vi1, region%mesh%vi2
          if (region%ice%Hi( vi) > 0) then
            ! Accumulate values only where ice is present
            field%accum( vi) = field%accum( vi) + SMB_loc( vi) * ice_density / sec_per_year * deltat
          end if
        end do
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine accumulate_single_ISMIP_flux_field

  subroutine write_to_ISMIP_regional_output_files_grid( region)
    !< Write to ISMIP regional output NetCDF files - grid version

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

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

    ! Basic topography
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%lithk)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%orog)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%topg)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%base)

    ! Geothermal heat flux
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%hfgeoubed)

    ! Surface and basal mass balance
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%acabf)

    ! Thickness tendency


    ! Velocities
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%xvelsurf)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%yvelsurf)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%zvelsurf)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%xvelbase)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%yvelbase)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%zvelbase)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%xvelmean)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%yvelmean)

    ! Temperatures
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%litemptop)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%litempbotgr)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%litempbotfl)

    ! Basal drag
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%strbasemag)

    ! Lateral mass balance

    ! Area fractions
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%sftgif)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%sftgrf)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%sftflf)

    ! Set previous time step to current
    region%ismip_grid_output%t_prev = region%ismip_grid_output%t_curr

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_ISMIP_regional_output_files_grid

  subroutine write_to_single_ISMIP_regional_output_file_grid( region, field)
    !< Write to single ISMIP regional output NetCDF file - grid version

    ! In/output variables:
    type(type_model_region),           intent(inout) :: region
    type(type_ismip_gridded_field),    intent(inout) :: field

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

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_single_ISMIP_regional_output_file_grid

  subroutine write_to_ISMIP_regional_output_file_grid_field( region, field, ncid)
    !< Write a single field to the ISMIP regional output NetCDF file - grid version

    ! In/output variables:
    type(type_model_region),           intent(inout) :: region
    type(type_ismip_gridded_field),    intent(inout) :: field
    integer,                           intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_ISMIP_regional_output_file_grid_field'
    real(dp), dimension(:),   allocatable :: d_grid_vec_partial_2D, d_mesh_vec_partial_2D
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: SMB_loc
    integer                               :: vi
    real(dp)                              :: deltat

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine deltat since previous writing (should be equal to 0 (first time step) or dt_ismip_output)
    deltat = region%ismip_grid_output%t_curr - region%ismip_grid_output%t_prev

    ! Allocate memory
    allocate( d_grid_vec_partial_2D( region%output_grid%n_loc ))

    ! Add the specified data field to the file
    select case (field%name)
      case default
        call crash('unknown choice_output_field "' // trim( field%name) // '"')

      ! Basic topography (ST)
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

      ! Geothermal heat flux (FL)
      case ('hfgeoubed')
        ! This is always a snapshot
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        ! First timeframe, no accumulation yet. Following protocol, using snapshot field instead
        d_mesh_vec_partial_2D = region%ice%geothermal_heat_flux / sec_per_year
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
        deallocate( d_mesh_vec_partial_2D)

      ! Surface and basal mass balance (FL)
      case ('acabf')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        if (field%is_initial) then
          ! First timeframe, no accumulation yet. Following protocol, using snapshot field instead
          SMB_loc( region%mesh%vi1: region%mesh%vi2) = region%SMB%SMB( region%mesh%vi1: region%mesh%vi2)
          d_mesh_vec_partial_2D = SMB_loc * ice_density / sec_per_year
          field%is_initial = .false.
        else
          d_mesh_vec_partial_2D = field%accum / deltat
          field%accum( region%mesh%vi1: region%mesh%vi2) = 0._dp
        end if
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
        deallocate( d_mesh_vec_partial_2D)

      ! Thickness tendency (FL)


      ! Velocities (ST)
      case ('xvelsurf')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%u_surf_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('yvelsurf')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%v_surf_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('zvelsurf')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%w_surf, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('xvelbase')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%u_base_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('yvelbase')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%v_base_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('zvelbase')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%w_base, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('xvelmean')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%u_vav_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)
      case ('yvelmean')
        call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%v_vav_b, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D / sec_per_year)

      ! Temperatures
      case ('litemptop')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%Ti( :,1), d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
      case ('litempbotgr')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        do vi = region%mesh%vi1, region%mesh%vi2
          if (region%ice%fraction_gr( vi) > 0._dp) then
            d_mesh_vec_partial_2D( vi) = region%ice%Ti( vi, region%mesh%nz)
          else
            d_mesh_vec_partial_2D( vi) = NaN
          end if
        end do
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
        deallocate( d_mesh_vec_partial_2D)
      case ('litempbotfl')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        do vi = region%mesh%vi1, region%mesh%vi2
          if (region%ice%mask_floating_ice( vi)) then
            d_mesh_vec_partial_2D( vi) = region%ice%Ti( vi, region%mesh%nz)
          else
            d_mesh_vec_partial_2D( vi) = NaN
          end if
        end do
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
        deallocate( d_mesh_vec_partial_2D)

      ! Basal drag
      case ('strbasemag')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%basal_shear_stress, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)

      ! Lateral mass balance

      ! Area fractions
      case ('sftgif')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%fraction_margin, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
      case ('sftgrf')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%fraction_gr, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
      case ('sftflf')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        ! TODO check whether this is appropriate
        d_mesh_vec_partial_2D = region%ice%fraction_margin( region%mesh%vi1: region%mesh%vi2) - region%ice%fraction_gr( region%mesh%vi1: region%mesh%vi2)
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
        deallocate( d_mesh_vec_partial_2D)
    end select

    ! Clean up memory
    deallocate( d_grid_vec_partial_2D)

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
    if (par%primary) call system('mkdir ' // trim(region%ismip_grid_output%folder))

    ! Create all grid files

    ! Basic topography
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%lithk)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%orog)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%topg)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%base)

    ! Geothermal heat flux
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%hfgeoubed)

    ! Surface and basal mass balance
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%acabf)


    ! Thickness tendency


    ! Velocities
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%xvelsurf)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%yvelsurf)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%zvelsurf)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%xvelbase)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%yvelbase)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%zvelbase)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%xvelmean)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%yvelmean)

    ! Temperatures
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%litemptop)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%litempbotgr)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%litempbotfl)

    ! Basal drag
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%strbasemag)

    ! Lateral mass balance

    ! Area fractions
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%sftgif)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%sftgrf)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%sftflf)

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

    ! Initialise the time trackers
    region%ismip_grid_output%t_prev = C%start_time_of_run
    region%ismip_grid_output%t_curr = C%start_time_of_run

    ! Basic topography
    call initialise_ISMIP_field( region, region%ismip_grid_output%lithk, 'lithk', 'land_ice_thickness', 'm', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%orog,  'orog' , 'surface_altitude', 'm', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%topg,  'topg' , 'bedrock_altitude', 'm', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%base,  'base' , 'base_altitude', 'm', 'ST')

    ! Geothermal heat flux
    call initialise_ISMIP_field( region, region%ismip_grid_output%hfgeoubed, 'hfgeoubed' , 'upward_geothermal_heat_flux_in_land_ice', 'W m-2', 'FL')

    ! Surface and basal mass balances
    call initialise_ISMIP_field( region, region%ismip_grid_output%acabf, 'acabf' , 'land_ice_surface_specific_mass_balance_flux', 'kg m-2 s-1', 'FL')

    ! Thickness tendency

    ! Velocities
    call initialise_ISMIP_field( region, region%ismip_grid_output%xvelsurf, 'xvelsurf' , 'land_ice_surface_x_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%yvelsurf, 'yvelsurf' , 'land_ice_surface_y_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%zvelsurf, 'zvelsurf' , 'land_ice_surface_upward_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%xvelbase, 'xvelbase' , 'land_ice_basal_x_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%yvelbase, 'yvelbase' , 'land_ice_basal_y_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%zvelbase, 'zvelbase' , 'land_ice_basal_upward_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%xvelmean, 'xvelmean' , 'land_ice_vertical_mean_x_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%yvelmean, 'yvelmean' , 'land_ice_vertical_mean_y_velocity', 'm s-1', 'ST')

    ! Temperatures
    call initialise_ISMIP_field( region, region%ismip_grid_output%litemptop, 'litemptop' , 'temperature_at_top_of_ice_sheet_model', 'K', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%litempbotgr, &
      'litempbotgr' , 'temperature_at_base_of_ice_sheet_model', 'K', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%litempbotfl, &
      'litempbotfl' , 'temperature_at_base_of_ice_sheet_model', 'K', 'ST')

    ! Basal drag
    call initialise_ISMIP_field( region, region%ismip_grid_output%strbasemag,  'strbasemag' , 'land_ice_basal_drag', 'Pa', 'ST')

    ! Lateral mass balance

    ! Area fractions
    call initialise_ISMIP_field( region, region%ismip_grid_output%sftgif, 'sftgif' , 'land_ice_area_fraction', '', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%sftgrf, 'sftgrf' , 'grounded_ice_sheet_area_fraction', '', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%sftflf, 'sftflf' , 'floating_ice_shelf_area_fraction', '', 'ST')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ISMIP_grid_output

  subroutine initialise_ISMIP_field( region, field, name, long_name, units, fieldtype)
    ! Initialise a single field

    type(type_model_region),           intent(inout) :: region
    type(type_ismip_gridded_field),    intent(inout) :: field
    character(len=*),                  intent(in   ) :: name
    character(len=*),                  intent(in   ) :: long_name
    character(len=*),                  intent(in   ) :: units
    character(len=2),                  intent(in   ) :: fieldtype

    ! Local variables
    character(len=1024), parameter :: routine_name = 'initialise_ISMIP_field'
    integer                        :: res_int
    character(len=2)               :: res_str

    ! Add routine to path
    call init_routine( routine_name)

    ! Inherit metadata
    field%name      = name
    field%long_name = long_name
    field%units     = units
    field%fieldtype = fieldtype

    ! Convert grid resolution to string
    res_int = int(C%dx_output_grid_ANT / 1000._dp)
    write(res_str, '(I2.2)') res_int

    ! Define the name of the subfolder
    region%ismip_grid_output%folder = trim( C%output_dir) // trim( C%ismip_exp_name) // '_' // trim(res_str) // '/'

    ! Define the filename for this field
    field%filename = trim( region%ismip_grid_output%folder) // trim(field%name) // '_' // &
      trim(region%ismip_grid_output%IS_name) // '_' // trim(C%ismip_group_name) // '_' // &
      trim(C%ismip_model_name) // '_' // trim(C%ismip_exp_name) // '.nc'

    ! Allocate fields for accumulation during each timestep for flux fields
    if (field%fieldtype == 'FL') then
      allocate(field%accum( region%mesh%vi1: region%mesh%vi2), source=0._dp)
      field%is_initial = .true.
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ISMIP_field

  ! TODO remap_ISMIP_fields -> reallocate_bounds

end module ismip_grid_output_files
