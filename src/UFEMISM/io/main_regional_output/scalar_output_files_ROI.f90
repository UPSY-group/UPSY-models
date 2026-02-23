module scalar_output_files_ROI

  use mpi_basic, only: par
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning
  use model_configuration, only: C
  use region_types, only: type_model_region
  use netcdf_io_main
  use reallocate_mod
  use netcdf_basic
  use netcdf, only: NF90_DOUBLE, NF90_UNLIMITED

  implicit none

  private

  public :: create_scalar_regional_output_file_ROI, buffer_scalar_output_ROI, write_to_scalar_regional_output_file_ROI, &
           create_ISMIP_scalar_regional_output_file_ROI, buffer_ISMIP_scalar_output_ROI, write_to_ISMIP_scalar_regional_output_file_ROI, &
            write_buffer_to_scalar_file_single_variable_ROI

  interface write_buffer_to_scalar_file_single_variable_ROI
    procedure :: write_buffer_to_scalar_file_single_variable_int_ROI
    procedure :: write_buffer_to_scalar_file_single_variable_dp_ROI
  end interface write_buffer_to_scalar_file_single_variable_ROI

contains

  subroutine write_to_scalar_regional_output_file_ROI( region)
    !< Write to the scalar regional output NetCDF file

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_scalar_regional_output_file_ROI'
    character(len=1024)            :: filename
    integer                        :: ncid, n, id_dim_time, ti, i_ROI

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    do i_ROI=1, region%nROI

      ! Print to terminal
      if (par%primary) write(0,'(A)') '   Writing to ROI scalar output file "' // &
        UPSY%stru%colour_string( trim( region%output_filenames_scalar_ROI(i_ROI)), 'light blue') // '"...'

      ! Shorthand for variable names
      filename = region%output_filenames_scalar_ROI(i_ROI)
      n        = region%scalars_ROI(i_ROI)%buffer%n

      ! Open the NetCDF file
      call open_existing_netcdf_file_for_writing( filename, ncid)

      ! Inquire number of timeframes already present in the file
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

      ! Write the time to the file
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'time',              region%scalars_ROI(i_ROI)%buffer%time,              n, ti+1)

      ! Integrated ice geometry
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'ice_area',          region%scalars_ROI(i_ROI)%buffer%ice_area,          n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'ice_volume',        region%scalars_ROI(i_ROI)%buffer%ice_volume,        n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'ice_volume_af',     region%scalars_ROI(i_ROI)%buffer%ice_volume_af,     n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'ice_area_PD',       region%scalars_ROI(i_ROI)%buffer%ice_area_PD,       n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'ice_volume_PD',     region%scalars_ROI(i_ROI)%buffer%ice_volume_PD,     n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'ice_volume_af_PD',  region%scalars_ROI(i_ROI)%buffer%ice_volume_af_PD,  n, ti+1)

      ! Integrated ice shelf geometry
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'ice_shelf_area',    region%scalars_ROI(i_ROI)%buffer%ice_shelf_area,    n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'ice_shelf_volume',  region%scalars_ROI(i_ROI)%buffer%ice_shelf_volume,  n, ti+1)


      ! Integrated mass fluxes
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'SMB_total',         region%scalars_ROI(i_ROI)%buffer%SMB_total,         n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'SMB_gr',            region%scalars_ROI(i_ROI)%buffer%SMB_gr,            n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'SMB_fl',            region%scalars_ROI(i_ROI)%buffer%SMB_fl,            n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'SMB_land',          region%scalars_ROI(i_ROI)%buffer%SMB_land,          n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'SMB_ocean',         region%scalars_ROI(i_ROI)%buffer%SMB_ocean,         n, ti+1)

      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'BMB_total',         region%scalars_ROI(i_ROI)%buffer%BMB_total,         n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'BMB_gr',            region%scalars_ROI(i_ROI)%buffer%BMB_gr,            n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'BMB_fl',            region%scalars_ROI(i_ROI)%buffer%BMB_fl,            n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'BMB_land',          region%scalars_ROI(i_ROI)%buffer%BMB_land,          n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'BMB_ocean',         region%scalars_ROI(i_ROI)%buffer%BMB_ocean,         n, ti+1)

      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'LMB_total',         region%scalars_ROI(i_ROI)%buffer%LMB_total,         n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'LMB_gr',            region%scalars_ROI(i_ROI)%buffer%LMB_gr,            n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'LMB_fl',            region%scalars_ROI(i_ROI)%buffer%LMB_fl,            n, ti+1)

      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'AMB_total',         region%scalars_ROI(i_ROI)%buffer%AMB_total,         n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'AMB_gr',            region%scalars_ROI(i_ROI)%buffer%AMB_gr,            n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'AMB_fl',            region%scalars_ROI(i_ROI)%buffer%AMB_fl,            n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'AMB_land',          region%scalars_ROI(i_ROI)%buffer%AMB_land,          n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'AMB_ocean',         region%scalars_ROI(i_ROI)%buffer%AMB_ocean,         n, ti+1)

      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'gl_flux',           region%scalars_ROI(i_ROI)%buffer%gl_flux,           n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'cf_gr_flux',        region%scalars_ROI(i_ROI)%buffer%cf_gr_flux,        n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'cf_fl_flux',        region%scalars_ROI(i_ROI)%buffer%cf_fl_flux,        n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'margin_land_flux',  region%scalars_ROI(i_ROI)%buffer%margin_land_flux,  n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'margin_ocean_flux', region%scalars_ROI(i_ROI)%buffer%margin_ocean_flux, n, ti+1)

      ! ! Numerical stability info
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'dt_ice',            region%scalars_ROI(i_ROI)%buffer%dt_ice,            n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'n_visc_its',        region%scalars_ROI(i_ROI)%buffer%n_visc_its,        n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'n_Axb_its',         region%scalars_ROI(i_ROI)%buffer%n_Axb_its,         n, ti+1)

      ! Reset buffer
      region%scalars_ROI(i_ROI)%buffer%n = 0

      ! Close the file
      call close_netcdf_file( ncid)

    end do ! i_ROI=1, region%nROI

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_scalar_regional_output_file_ROI

  subroutine create_scalar_regional_output_file_ROI( region)
    !< Create the scalar regional output NetCDF file

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_scalar_regional_output_file_ROI'
    character(len=1024)            :: filename_base, filename
    integer                        :: ncid, i_ROI

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    do i_ROI = 1, region%nROI

      ! Get the filename
      filename = region%output_filenames_scalar_ROI(i_ROI)

      ! Print to terminal
      if (par%primary) write(0,'(A)') '   Creating scalar output file "' // &
        UPSY%stru%colour_string( trim( filename), 'light blue') // '"...'

      ! Create the NetCDF file
      call create_new_netcdf_file_for_writing( filename, ncid)

      ! Add time, zeta, and month dimensions+variables to the file
      call add_time_dimension_to_file( filename, ncid)

      ! Integrated ice geometry
      call add_field_dp_0D( filename, ncid, 'ice_area',          long_name = 'Total ice area', units = 'm^2')
      call add_field_dp_0D( filename, ncid, 'ice_volume',        long_name = 'Total ice volume', units = 'm s.l.e.')
      call add_field_dp_0D( filename, ncid, 'ice_volume_af',     long_name = 'Total ice volume above floatation', units = 'm s.l.e.')

      call add_field_dp_0D( filename, ncid, 'ice_area_PD',       long_name = 'Total ice area for present-day', units = 'm^2')
      call add_field_dp_0D( filename, ncid, 'ice_volume_PD',     long_name = 'Total ice volume for present-day', units = 'm s.l.e.')
      call add_field_dp_0D( filename, ncid, 'ice_volume_af_PD',  long_name = 'Total ice volume above floatation for present-day', units = 'm s.l.e.')

      ! Integrated ice shelf geometry
      call add_field_dp_0D( filename, ncid, 'ice_shelf_area',    long_name = 'Total ice shelf area', units = 'm^2')
      call add_field_dp_0D( filename, ncid, 'ice_shelf_volume',  long_name = 'Total ice shelf volume', units = 'm^3')

      ! Integrated mass fluxes
      call add_field_dp_0D( filename, ncid, 'SMB_total',         long_name = 'Area-integrated total SMB', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'SMB_gr',            long_name = 'Area-integrated ice sheet SMB', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'SMB_fl',            long_name = 'Area-integrated ice shelf SMB', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'SMB_land',          long_name = 'Area-integrated ice-free land SMB', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'SMB_ocean',         long_name = 'Area-integrated ice-free ocean SMB', units = 'Gt yr^-1')

      call add_field_dp_0D( filename, ncid, 'BMB_total',         long_name = 'Area-integrated total BMB', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'BMB_gr',            long_name = 'Area-integrated ice sheet BMB', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'BMB_fl',            long_name = 'Area-integrated ice shelf BMB', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'BMB_land',          long_name = 'Area-integrated ice-free land BMB', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'BMB_ocean',         long_name = 'Area-integrated ice-free ocean BMB', units = 'Gt yr^-1')

      call add_field_dp_0D( filename, ncid, 'LMB_total',         long_name = 'Area-integrated total LMB', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'LMB_gr',            long_name = 'Area-integrated ice sheet LMB', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'LMB_fl',            long_name = 'Area-integrated ice shelf LMB', units = 'Gt yr^-1')

      call add_field_dp_0D( filename, ncid, 'AMB_total',         long_name = 'Area-integrated total additional MB from other sources', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'AMB_gr',            long_name = 'Area-integrated ice sheet additional MB from other sources', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'AMB_fl',            long_name = 'Area-integrated ice shelf additional MB from other sources', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'AMB_land',          long_name = 'Area-integrated ice-free land additional MB from other sources', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'AMB_ocean',         long_name = 'Area-integrated ice-free ocean additional MB from other sources', units = 'Gt yr^-1')

      call add_field_dp_0D( filename, ncid, 'gl_flux',           long_name = 'Total lateral grounding line flux', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'cf_gr_flux',        long_name = 'Total lateral grounded calving front flux', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'cf_fl_flux',        long_name = 'Total lateral floating calving front flux', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'margin_land_flux',  long_name = 'Total lateral flux exiting the ice margin into ground', units = 'Gt yr^-1')
      call add_field_dp_0D( filename, ncid, 'margin_ocean_flux', long_name = 'Total lateral flux exiting the ice margin into water', units = 'Gt yr^-1')

      ! Numerical stability info
      call add_field_dp_0D(  filename, ncid, 'dt_ice',           long_name = 'Ice-dynamical time step', units = 'yr')
      call add_field_int_0D( filename, ncid, 'n_visc_its',       long_name = 'Number of non-linear viscosity iterations')
      call add_field_int_0D( filename, ncid, 'n_Axb_its',        long_name = 'Number of iterations in iterative solver for linearised momentum balance')

      ! Allocate memory to buffer scalar output data between output writing intervals
      call allocate_scalar_output_buffer_ROI( region, i_ROI)

      ! Close the file
      call close_netcdf_file( ncid)
    end do ! i_ROI = 1, region%nROI

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_scalar_regional_output_file_ROI

  subroutine allocate_scalar_output_buffer_ROI( region, i_ROI)
    !< Allocate memory to buffer the scalar output data between output writing intervals

    ! In/output variables:
    type(type_model_region), intent(inout) :: region
    integer,                 intent(in)    :: i_ROI

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_scalar_output_buffer_ROI'
    integer                        :: n_mem

    ! Add routine to path
    call init_routine( routine_name)

    region%scalars_ROI(i_ROI)%buffer%n_mem = 0
    region%scalars_ROI(i_ROI)%buffer%n     = 0

    ! Only allocate memory for this on the primary
    if (par%primary) then

      n_mem = 1000
      region%scalars_ROI(i_ROI)%buffer%n_mem = n_mem
      region%scalars_ROI(i_ROI)%buffer%n     = 0

      allocate( region%scalars_ROI(i_ROI)%buffer%time             ( n_mem), source = 0._dp)

      allocate( region%scalars_ROI(i_ROI)%buffer%ice_area         ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ice_volume       ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ice_volume_af    ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ice_area_PD      ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ice_volume_PD    ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ice_volume_af_PD ( n_mem), source = 0._dp)

      allocate( region%scalars_ROI(i_ROI)%buffer%ice_shelf_area   ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ice_shelf_volume ( n_mem), source = 0._dp)

      allocate( region%scalars_ROI(i_ROI)%buffer%SMB_total        ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%SMB_gr           ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%SMB_fl           ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%SMB_land         ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%SMB_ocean        ( n_mem), source = 0._dp)

      allocate( region%scalars_ROI(i_ROI)%buffer%BMB_total        ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%BMB_gr           ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%BMB_fl           ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%BMB_land         ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%BMB_ocean        ( n_mem), source = 0._dp)

      allocate( region%scalars_ROI(i_ROI)%buffer%LMB_total        ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%LMB_gr           ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%LMB_fl           ( n_mem), source = 0._dp)

      allocate( region%scalars_ROI(i_ROI)%buffer%AMB_total        ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%AMB_gr           ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%AMB_fl           ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%AMB_land         ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%AMB_ocean        ( n_mem), source = 0._dp)

      allocate( region%scalars_ROI(i_ROI)%buffer%gl_flux          ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%cf_gr_flux       ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%cf_fl_flux       ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%margin_land_flux ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%margin_ocean_flux( n_mem), source = 0._dp)

      allocate( region%scalars_ROI(i_ROI)%buffer%dt_ice           ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%n_visc_its       ( n_mem), source = 0)
      allocate( region%scalars_ROI(i_ROI)%buffer%n_Axb_its        ( n_mem), source = 0)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_scalar_output_buffer_ROI

  subroutine buffer_scalar_output_ROI( region)
    !< Buffer the scalar output data between output writing intervals

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'buffer_scalar_output'
    integer                        :: n, i_ROI

    ! Add routine to path
    call init_routine( routine_name)

    ! Only the primary does this
    if (par%primary) then

      do i_ROI=1, region%nROI

        ! Increase timeframe count
        region%scalars_ROI(i_ROI)%buffer%n = region%scalars_ROI(i_ROI)%buffer%n + 1
        n = region%scalars_ROI(i_ROI)%buffer%n

        ! Extend buffer memory if necessary
        if (n > region%scalars_ROI(i_ROI)%buffer%n_mem - 10) call extend_scalar_output_buffer_ROI( region)

        ! Store new timeframe in buffer
        region%scalars_ROI(i_ROI)%buffer%time             ( n) = region%time

        region%scalars_ROI(i_ROI)%buffer%ice_area         ( n) = region%scalars_ROI(i_ROI)%ice_area
        region%scalars_ROI(i_ROI)%buffer%ice_volume       ( n) = region%scalars_ROI(i_ROI)%ice_volume
        region%scalars_ROI(i_ROI)%buffer%ice_volume_af    ( n) = region%scalars_ROI(i_ROI)%ice_volume_af
        region%scalars_ROI(i_ROI)%buffer%ice_area_PD      ( n) = region%scalars_ROI(i_ROI)%ice_area_PD
        region%scalars_ROI(i_ROI)%buffer%ice_volume_PD    ( n) = region%scalars_ROI(i_ROI)%ice_volume_PD
        region%scalars_ROI(i_ROI)%buffer%ice_volume_af_PD ( n) = region%scalars_ROI(i_ROI)%ice_volume_af_PD

        region%scalars_ROI(i_ROI)%buffer%ice_shelf_area   ( n) = region%scalars_ROI(i_ROI)%ice_shelf_area
        region%scalars_ROI(i_ROI)%buffer%ice_shelf_volume ( n) = region%scalars_ROI(i_ROI)%ice_shelf_volume

        region%scalars_ROI(i_ROI)%buffer%SMB_total        ( n) = region%scalars_ROI(i_ROI)%SMB_total
        region%scalars_ROI(i_ROI)%buffer%SMB_gr           ( n) = region%scalars_ROI(i_ROI)%SMB_gr
        region%scalars_ROI(i_ROI)%buffer%SMB_fl           ( n) = region%scalars_ROI(i_ROI)%SMB_fl
        region%scalars_ROI(i_ROI)%buffer%SMB_land         ( n) = region%scalars_ROI(i_ROI)%SMB_land
        region%scalars_ROI(i_ROI)%buffer%SMB_ocean        ( n) = region%scalars_ROI(i_ROI)%SMB_ocean

        region%scalars_ROI(i_ROI)%buffer%BMB_total        ( n) = region%scalars_ROI(i_ROI)%BMB_total
        region%scalars_ROI(i_ROI)%buffer%BMB_gr           ( n) = region%scalars_ROI(i_ROI)%BMB_gr
        region%scalars_ROI(i_ROI)%buffer%BMB_fl           ( n) = region%scalars_ROI(i_ROI)%BMB_fl
        region%scalars_ROI(i_ROI)%buffer%BMB_land         ( n) = region%scalars_ROI(i_ROI)%BMB_land
        region%scalars_ROI(i_ROI)%buffer%BMB_ocean        ( n) = region%scalars_ROI(i_ROI)%BMB_ocean

        region%scalars_ROI(i_ROI)%buffer%LMB_total        ( n) = region%scalars_ROI(i_ROI)%LMB_total
        region%scalars_ROI(i_ROI)%buffer%LMB_gr           ( n) = region%scalars_ROI(i_ROI)%LMB_gr
        region%scalars_ROI(i_ROI)%buffer%LMB_fl           ( n) = region%scalars_ROI(i_ROI)%LMB_fl

        region%scalars_ROI(i_ROI)%buffer%AMB_total        ( n) = region%scalars_ROI(i_ROI)%AMB_total
        region%scalars_ROI(i_ROI)%buffer%AMB_gr           ( n) = region%scalars_ROI(i_ROI)%AMB_gr
        region%scalars_ROI(i_ROI)%buffer%AMB_fl           ( n) = region%scalars_ROI(i_ROI)%AMB_fl
        region%scalars_ROI(i_ROI)%buffer%AMB_land         ( n) = region%scalars_ROI(i_ROI)%AMB_land
        region%scalars_ROI(i_ROI)%buffer%AMB_ocean        ( n) = region%scalars_ROI(i_ROI)%AMB_ocean

        region%scalars_ROI(i_ROI)%buffer%gl_flux          ( n) = region%scalars_ROI(i_ROI)%gl_flux
        region%scalars_ROI(i_ROI)%buffer%cf_gr_flux       ( n) = region%scalars_ROI(i_ROI)%cf_gr_flux
        region%scalars_ROI(i_ROI)%buffer%cf_fl_flux       ( n) = region%scalars_ROI(i_ROI)%cf_fl_flux
        region%scalars_ROI(i_ROI)%buffer%margin_land_flux ( n) = region%scalars_ROI(i_ROI)%margin_land_flux
        region%scalars_ROI(i_ROI)%buffer%margin_ocean_flux( n) = region%scalars_ROI(i_ROI)%margin_ocean_flux

        region%scalars_ROI(i_ROI)%buffer%dt_ice           ( n) = region%ice%dt_ice
        region%scalars_ROI(i_ROI)%buffer%n_visc_its       ( n) = region%ice%n_visc_its
        region%scalars_ROI(i_ROI)%buffer%n_Axb_its        ( n) = region%ice%n_Axb_its
      end do ! i_ROI=1, region%nROI
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine buffer_scalar_output_ROI

  subroutine extend_scalar_output_buffer_ROI( region)
    !< Extend memory to buffer the scalar output data between output writing intervals
    !
    ! NOTE: should only be called by the primary!

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'extend_scalar_output_buffer_ROI'
    integer                        :: n_mem, i_ROI

    ! Add routine to path
    call init_routine( routine_name)

    do i_ROI=1, region%nROI

      n_mem = region%scalars%buffer%n_mem * 2
      region%scalars%buffer%n_mem = n_mem

      call reallocate( region%scalars%buffer%time             , n_mem, source = 0._dp)

      call reallocate( region%scalars_ROI(i_ROI)%buffer%ice_area         , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%ice_volume       , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%ice_volume_af    , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%ice_area_PD      , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%ice_volume_PD    , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%ice_volume_af_PD , n_mem, source = 0._dp)

      call reallocate( region%scalars_ROI(i_ROI)%buffer%ice_shelf_area   , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%ice_shelf_volume , n_mem, source = 0._dp)

      call reallocate( region%scalars_ROI(i_ROI)%buffer%SMB_total        , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%SMB_gr           , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%SMB_fl           , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%SMB_land         , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%SMB_ocean        , n_mem, source = 0._dp)

      call reallocate( region%scalars_ROI(i_ROI)%buffer%BMB_total        , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%BMB_gr           , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%BMB_fl           , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%BMB_land         , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%BMB_ocean        , n_mem, source = 0._dp)

      call reallocate( region%scalars_ROI(i_ROI)%buffer%LMB_total        , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%LMB_gr           , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%LMB_fl           , n_mem, source = 0._dp)

      call reallocate( region%scalars_ROI(i_ROI)%buffer%AMB_total        , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%AMB_gr           , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%AMB_fl           , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%AMB_land         , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%AMB_ocean        , n_mem, source = 0._dp)

      call reallocate( region%scalars_ROI(i_ROI)%buffer%gl_flux          , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%cf_gr_flux       , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%cf_fl_flux       , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%margin_land_flux , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%margin_ocean_flux, n_mem, source = 0._dp)

      call reallocate( region%scalars_ROI(i_ROI)%buffer%dt_ice           , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%n_visc_its       , n_mem, source = 0)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%n_Axb_its        , n_mem, source = 0)

    end do ! i_ROI=1, region%nROI

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extend_scalar_output_buffer_ROI

  subroutine write_buffer_to_scalar_file_single_variable_int_ROI( filename, ncid, var_name, d, n, ti)
    !< Write buffered scalar data of a single variable to the scalar output file

    ! In/output variables:
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: var_name
    integer,  dimension(:), intent(in   ) :: d
    integer,                intent(in   ) :: n
    integer,                intent(in   ) :: ti

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_buffer_to_scalar_file_single_variable_int_ROI'
    integer                        :: id_var
    integer, dimension(1)          :: start, count
    integer,  dimension(n)         :: d_to_write

    ! Add routine to path
    call init_routine( routine_name)

    call inquire_var( filename, ncid, var_name, id_var)

    start = ti
    count = n
    d_to_write = d(1:n)

    call write_var_primary(  filename, ncid, id_var, d_to_write, start = start, count = count)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_buffer_to_scalar_file_single_variable_int_ROI

  subroutine write_buffer_to_scalar_file_single_variable_dp_ROI( filename, ncid, var_name, d, n, ti)
    !< Write buffered scalar data of a single variable to the scalar output file

    ! In/output variables:
    character(len=*),       intent(in   ) :: filename
    integer,                intent(in   ) :: ncid
    character(len=*),       intent(in   ) :: var_name
    real(dp), dimension(:), intent(in   ) :: d
    integer,                intent(in   ) :: n
    integer,                intent(in   ) :: ti

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_buffer_to_scalar_file_single_variable_dp_ROI'
    integer                        :: id_var
    integer, dimension(1)          :: start, count
    real(dp), dimension(n)         :: d_to_write

    ! Add routine to path
    call init_routine( routine_name)

    call inquire_var( filename, ncid, var_name, id_var)

    start = ti
    count = n
    d_to_write = d(1:n)

    call write_var_primary(  filename, ncid, id_var, d_to_write, start = start, count = count)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_buffer_to_scalar_file_single_variable_dp_ROI



  subroutine write_to_ISMIP_scalar_regional_output_file_ROI( region)
    !< Write to the scalar regional output NetCDF file

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_ISMIP_scalar_regional_output_file_ROI'
    character(len=1024)            :: filename
    integer                        :: ncid, n, id_dim_time, ti, i_ROI

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_ISMIP_output) then
      call finalise_routine( routine_name)
      return
    end if

    do i_ROI=1, region%nROI

      ! Print to terminal
      if (par%primary) write(0,'(A)') '   Writing to ROI ISMIP scalar output file "' // &
        UPSY%stru%colour_string( trim( region%output_filenames_ISMIP_scalar_ROI(i_ROI)), 'light blue') // '"...'

      ! Shorthand for variable names
      filename = region%output_filenames_ISMIP_scalar_ROI(i_ROI)
      n        = region%scalars_ROI(i_ROI)%buffer%ismip%n

      ! Open the NetCDF file
      call open_existing_netcdf_file_for_writing( filename, ncid)

      ! Inquire number of timeframes already present in the file
      call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = ti)

      ! Write the time to the file
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'time',               region%scalars_ROI(i_ROI)%buffer%ismip%time,               n, ti+1)

      ! Integrated ice geometry
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'lim',                region%scalars_ROI(i_ROI)%buffer%ismip%lim,                n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'limnsw',             region%scalars_ROI(i_ROI)%buffer%ismip%limnsw,             n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'iareagr',            region%scalars_ROI(i_ROI)%buffer%ismip%iareagr,            n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'iareafl',            region%scalars_ROI(i_ROI)%buffer%ismip%iareafl,            n, ti+1)

      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'tendacabf',          region%scalars_ROI(i_ROI)%buffer%ismip%tendacabf,          n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'tendlibmassbf',      region%scalars_ROI(i_ROI)%buffer%ismip%tendlibmassbf,      n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'tendlibmassbffl',    region%scalars_ROI(i_ROI)%buffer%ismip%tendlibmassbffl,    n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'tendlicalvf',        region%scalars_ROI(i_ROI)%buffer%ismip%tendlicalvf,        n, ti+1)
      call write_buffer_to_scalar_file_single_variable_ROI( filename, ncid, 'tendlifmassbf',      region%scalars_ROI(i_ROI)%buffer%ismip%tendlifmassbf,      n, ti+1)

      ! Reset buffer
      region%scalars_ROI(i_ROI)%buffer%ismip%n = 0

      ! Close the file
      call close_netcdf_file( ncid)

    end do ! i_ROI=1, region%nROI

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_ISMIP_scalar_regional_output_file_ROI

  subroutine create_ISMIP_scalar_regional_output_file_ROI( region)
    !< Create the scalar regional output NetCDF file

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_ISMIP_scalar_regional_output_file_ROI'
    character(len=1024)            :: filename_base, filename
    integer                        :: ncid, i_ROI
    integer                        :: id_dim_time
    integer                        :: id_var_time

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_ISMIP_output) then
      call finalise_routine( routine_name)
      return
    end if

    do i_ROI = 1, region%nROI

      ! Get the filename
      filename = region%output_filenames_ISMIP_scalar_ROI(i_ROI)

      ! Print to terminal
      if (par%primary) write(0,'(A)') '   Creating scalar output file "' // &
        UPSY%stru%colour_string( trim( filename), 'light blue') // '"...'

      ! Create the NetCDF file
      call create_new_netcdf_file_for_writing( filename, ncid)

      ! Add time dimension to the file manually because of different units
      ! call add_time_dimension_to_file( filename, ncid)
      call create_dimension(   filename, ncid, get_first_option_from_list( field_name_options_time), NF90_UNLIMITED, id_dim_time)
      call create_variable(    filename, ncid, get_first_option_from_list( field_name_options_time), NF90_DOUBLE, (/ id_dim_time /), id_var_time)
      call add_attribute_char( filename, ncid, id_var_time, 'long_name', 'Time')
      call add_attribute_char( filename, ncid, id_var_time, 'units', 'days')

      ! Integrated ice geometry
      call add_field_dp_0D( filename, ncid, 'lim',        long_name = 'land_ice_mass',                          units = 'kg')
      call add_field_dp_0D( filename, ncid, 'limnsw',     long_name = 'land_ice_mass_not_displacing_sea_water', units = 'kg')
      call add_field_dp_0D( filename, ncid, 'iareagr',    long_name = 'grounded_ice_sheet_area',                units = 'm2')
      call add_field_dp_0D( filename, ncid, 'iareafl',    long_name = 'floating_ice_sheet_area',                units = 'm2')

      call add_field_dp_0D( filename, ncid, 'tendacabf',       long_name = 'tendency_of_land_ice_mass_due_to_surface_mass_balance',          units = 'kg s-1')
      call add_field_dp_0D( filename, ncid, 'tendlibmassbf',   long_name = 'tendency_of_land_ice_mass_due_to_basal_mass_balance',            units = 'kg s-1')
      call add_field_dp_0D( filename, ncid, 'tendlibmassbffl', long_name = 'tendency_of_land_ice_mass_due_to_basal_mass_balance',            units = 'kg s-1')
      call add_field_dp_0D( filename, ncid, 'tendlicalvf',     long_name = 'tendency_of_land_ice_mass_due_to_calving',                       units = 'kg s-1')
      call add_field_dp_0D( filename, ncid, 'tendlifmassbf',   long_name = 'tendency_of_land_ice_mass_due_to_calving_and_ice_front_melting', units = 'kg s-1')

      ! Allocate memory to buffer scalar output data between output writing intervals
      call allocate_ISMIP_scalar_output_buffer_ROI( region, i_ROI)

      ! Close the file
      call close_netcdf_file( ncid)
    end do ! i_ROI = 1, region%nROI

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_ISMIP_scalar_regional_output_file_ROI

  subroutine allocate_ISMIP_scalar_output_buffer_ROI( region, i_ROI)
    !< Allocate memory to buffer the scalar output data between output writing intervals

    ! In/output variables:
    type(type_model_region), intent(inout) :: region
    integer,                 intent(in)    :: i_ROI

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_ISMIP_scalar_output_buffer_ROI'
    integer                        :: n_mem

    ! Add routine to path
    call init_routine( routine_name)

    region%scalars_ROI(i_ROI)%buffer%ismip%n_mem = 0
    region%scalars_ROI(i_ROI)%buffer%ismip%n     = 0

    ! Only allocate memory for this on the primary
    if (par%primary) then

      n_mem = 1000
      region%scalars_ROI(i_ROI)%buffer%ismip%n_mem = n_mem
      region%scalars_ROI(i_ROI)%buffer%ismip%n     = 0

      allocate( region%scalars_ROI(i_ROI)%buffer%ismip%time            ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ismip%lim             ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ismip%limnsw          ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ismip%iareagr         ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ismip%iareafl         ( n_mem), source = 0._dp)

      allocate( region%scalars_ROI(i_ROI)%buffer%ismip%tendacabf       ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ismip%tendlibmassbf   ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ismip%tendlibmassbffl ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ismip%tendlicalvf     ( n_mem), source = 0._dp)
      allocate( region%scalars_ROI(i_ROI)%buffer%ismip%tendlifmassbf   ( n_mem), source = 0._dp)

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_ISMIP_scalar_output_buffer_ROI

  subroutine buffer_ISMIP_scalar_output_ROI( region)
    !< Buffer the scalar output data between output writing intervals

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'buffer_buffer_ISMIP_scalar_output_ROIscalar_output'
    integer                        :: n, i_ROI

    ! Add routine to path
    call init_routine( routine_name)

    ! Only the primary does this
    if (par%primary) then

      do i_ROI=1, region%nROI

        ! Increase timeframe count
        region%scalars_ROI(i_ROI)%buffer%ismip%n = region%scalars_ROI(i_ROI)%buffer%ismip%n + 1
        n = region%scalars_ROI(i_ROI)%buffer%ismip%n

        ! Extend buffer memory if necessary
        if (n > region%scalars_ROI(i_ROI)%buffer%ismip%n_mem - 10) call extend_scalar_output_buffer_ROI( region)

        ! Store new timeframe in buffer
        region%scalars_ROI(i_ROI)%buffer%ismip%time    ( n) = region%time * 360._dp ! need to convert to number of days for ISMIP

        region%scalars_ROI(i_ROI)%buffer%ismip%lim     ( n) = region%scalars_ROI(i_ROI)%ismip%lim
        region%scalars_ROI(i_ROI)%buffer%ismip%limnsw  ( n) = region%scalars_ROI(i_ROI)%ismip%limnsw
        region%scalars_ROI(i_ROI)%buffer%ismip%iareagr ( n) = region%scalars_ROI(i_ROI)%ismip%iareagr
        region%scalars_ROI(i_ROI)%buffer%ismip%iareafl ( n) = region%scalars_ROI(i_ROI)%ismip%iareafl

        region%scalars_ROI(i_ROI)%buffer%ismip%tendacabf       ( n) = region%scalars_ROI(i_ROI)%ismip%tendacabf
        region%scalars_ROI(i_ROI)%buffer%ismip%tendlibmassbf   ( n) = region%scalars_ROI(i_ROI)%ismip%tendlibmassbf
        region%scalars_ROI(i_ROI)%buffer%ismip%tendlibmassbffl ( n) = region%scalars_ROI(i_ROI)%ismip%tendlibmassbffl
        region%scalars_ROI(i_ROI)%buffer%ismip%tendlicalvf     ( n) = region%scalars_ROI(i_ROI)%ismip%tendlicalvf
        region%scalars_ROI(i_ROI)%buffer%ismip%tendlifmassbf   ( n) = region%scalars_ROI(i_ROI)%ismip%tendlifmassbf

      end do ! i_ROI=1, region%nROI
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine buffer_ISMIP_scalar_output_ROI

  subroutine extend_ISMIP_scalar_output_buffer_ROI( region)
    !< Extend memory to buffer the scalar output data between output writing intervals
    !
    ! NOTE: should only be called by the primary!

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'extend_ISMIP_scalar_output_buffer_ROI'
    integer                        :: n_mem, i_ROI

    ! Add routine to path
    call init_routine( routine_name)

    do i_ROI=1, region%nROI

      n_mem = region%scalars%buffer%ismip%n_mem * 2
      region%scalars%buffer%ismip%n_mem = n_mem

      call reallocate( region%scalars_ROI(i_ROI)%buffer%ismip%time            , n_mem, source = 0._dp)

      call reallocate( region%scalars_ROI(i_ROI)%buffer%ismip%lim             , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%ismip%limnsw          , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%ismip%iareagr         , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%ismip%iareafl         , n_mem, source = 0._dp)

      call reallocate( region%scalars_ROI(i_ROI)%buffer%ismip%tendacabf       , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%ismip%tendlibmassbf   , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%ismip%tendlibmassbffl , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%ismip%tendlicalvf     , n_mem, source = 0._dp)
      call reallocate( region%scalars_ROI(i_ROI)%buffer%ismip%tendlifmassbf   , n_mem, source = 0._dp)

    end do ! i_ROI=1, region%nROI

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine extend_ISMIP_scalar_output_buffer_ROI


end module scalar_output_files_ROI
