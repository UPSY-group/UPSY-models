module ismip_output_files

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
  use UPSY_main, only: UPSY
  use parameters
  use precisions, only: dp
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use model_configuration, only: C
  use region_types, only: type_model_region
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use netcdf_io_main
  use ice_mass_and_fluxes, only: calc_ISMIP_fluxes
  use remapping_main, only: map_from_mesh_vertices_to_xy_grid_2D, &
    map_from_mesh_triangles_to_xy_grid_2D
  use mpi_distributed_memory, only: gather_to_all
  use ismip_output_types, only: type_ismip_output, type_ismip_gridded_field, type_ismip_scalar
  use mesh_zeta, only: vertical_average
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_local
  use netcdf, only: NF90_GLOBAL
  use reallocate_mod, only: reallocate_bounds

  implicit none

  private

  public :: create_ISMIP_regional_output_files, write_to_ISMIP_regional_output_files, &
    accumulate_ISMIP_flux_fields, remap_ISMIP_output

  interface create_single_ISMIP_regional_output_file
    module procedure create_single_ISMIP_regional_output_file_grid
    module procedure create_single_ISMIP_regional_output_file_scalar
  end interface

  interface write_to_single_ISMIP_regional_output_file
    module procedure write_to_single_ISMIP_regional_output_file_grid
  end interface

  interface write_to_file_scalar
    module procedure write_to_file_scalar_ST
    module procedure write_to_file_scalar_FL
  end interface

  interface initialise_ISMIP_field
    module procedure initialise_ISMIP_field_grid
    module procedure initialise_ISMIP_field_scalar
  end interface

contains

  subroutine accumulate_ISMIP_flux_fields( region)
    !< Write to ISMIP regional output NetCDF files - grid version

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter                       :: routine_name = 'accumulate_ISMIP_flux_fields'
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: calving_flux
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: gl_flux
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: SMB_loc
    logical, dimension(region%mesh%vi1:region%mesh%vi2)  :: mask_ice
    real( dp)                                            :: deltat
    integer                                              :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! if no ISMIP output should be created, do nothing
    if (.not. C%do_create_ismip_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Compute the calving flux
    call calc_ISMIP_fluxes( region%mesh, region%ice, calving_flux, gl_flux)

    ! Copy SMB for hybrid memory reasons
    SMB_loc( region%mesh%vi1: region%mesh%vi2) = region%SMB%SMB( region%mesh%vi1: region%mesh%vi2)

    ! Determine the mask where ice is present
    do vi = region%mesh%vi1, region%mesh%vi2
      if (region%ice%Hi( vi) > 0._dp) then
        mask_ice( vi) = .true.
      else
        mask_ice( vi) = .false.
      end if
    end do

    ! Get delta t since last current time
    deltat = region%time - region%ismip_output%t_curr

    ! Accumulate FL fields. Exceptions:
    ! - dlithkdt (accumulation is used to store the previous written thickness
    ! - hfgeoubed (doesn't vary, so just writing out the initial snapshot)
    ! - lifmassbf (not computed yet, so just spitting out zeros)
    call accumulate_single_ISMIP_flux_field( region, SMB_loc, mask_ice, deltat, &
      region%ismip_output%tendacabf, field = region%ismip_output%acabf)
    call accumulate_single_ISMIP_flux_field( region, region%BMB%BMB, region%ice%mask_grounded_ice, deltat, &
      region%ismip_output%tendlibmassbfgr, field = region%ismip_output%libmassbfgr)
    call accumulate_single_ISMIP_flux_field( region, region%BMB%BMB, region%ice%mask_floating_ice, deltat, &
      region%ismip_output%tendlibmassbffl, field = region%ismip_output%libmassbffl)
    call accumulate_single_ISMIP_flux_field( region, calving_flux, mask_ice, deltat, &
      region%ismip_output%tendlicalvf, field = region%ismip_output%licalvf)
    call accumulate_single_ISMIP_flux_field( region, gl_flux, mask_ice, deltat, &
      region%ismip_output%tendligroundf)

    ! Update current time
    region%ismip_output%t_curr = region%time

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine accumulate_ISMIP_flux_fields

  subroutine accumulate_single_ISMIP_flux_field( region, d_partial, mask, deltat, scalar, field)
    !< Write to ISMIP regional output NetCDF files - grid version

    ! In/output variables:
    type(type_model_region),                              intent(inout) :: region
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2), intent(in   ) :: d_partial
    logical,  dimension(region%mesh%vi1:region%mesh%vi2), intent(in   ) :: mask
    real(dp),                                             intent(in   ) :: deltat
    type(type_ismip_scalar),                              intent(inout) :: scalar
    type(type_ismip_gridded_field), optional,             intent(inout) :: field

    ! Local variables:
    character(len=1024), parameter                       :: routine_name = 'accumulate_single_ISMIP_flux_field'
    integer                                              :: vi, ierr
    real(dp)                                             :: scalar_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Accumulate field if requested
    if (present( field)) then
      do vi = region%mesh%vi1, region%mesh%vi2
        if (mask( vi)) then
          field%accum( vi) = field%accum( vi) + d_partial( vi) * ice_density / sec_per_year * deltat
        end if
      end do
    end if

    ! Initialise local scalar
    scalar_loc = 0._dp

    ! Compute integrated scalar
    do vi = region%mesh%vi1, region%mesh%vi2
      if (mask( vi)) then
        scalar_loc = scalar_loc + d_partial( vi) * ice_density * region%mesh%A( vi) / sec_per_year * deltat
      end if
    end do

    ! Add together values from each process
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalar_loc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Accumulate
    scalar%accum = scalar%accum + scalar_loc

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine accumulate_single_ISMIP_flux_field

  subroutine write_to_ISMIP_regional_output_files( region)
    !< Write to ISMIP regional output NetCDF files

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter                       :: routine_name = 'write_to_ISMIP_regional_output_files'
    logical,  dimension(region%mesh%vi1:region%mesh%vi2) :: mask_ice_a
    integer                                              :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! if no ISMIP output should be created, do nothing
    if (.not. C%do_create_ismip_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Determine masks
    do vi = region%mesh%vi1, region%mesh%vi2
      if (region%ice%Hi( vi) > 0._dp) then
        mask_ice_a( vi) = .true.
      else
        mask_ice_a( vi) = .false.
      end if
    end do

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to ISMIP output files' // '...'

    ! Basic topography
    call write_to_file_grid_ST_a( region, region%ismip_output%lithk, region%ice%Hi, vmin=0._dp)
    call write_to_file_grid_ST_a( region, region%ismip_output%orog, region%ice%Hs, vmin=0._dp)
    call write_to_file_grid_ST_a( region, region%ismip_output%topg, region%ice%Hb)
    call write_to_file_grid_ST_a( region, region%ismip_output%base, region%ice%Hib)

    ! Geothermal heat flux
    call write_to_single_ISMIP_regional_output_file( region, region%ismip_output%hfgeoubed)

    ! Surface and basal mass balance
    call write_to_file_grid_FL( region, region%ismip_output%acabf)
    call write_to_file_grid_FL( region, region%ismip_output%libmassbfgr)
    call write_to_file_grid_FL( region, region%ismip_output%libmassbffl)

    ! Thickness tendency
    call write_to_single_ISMIP_regional_output_file( region, region%ismip_output%dlithkdt)

    ! Velocities
    call write_to_file_grid_ST_b( region, region%ismip_output%xvelsurf, region%ice%u_surf_b / sec_per_year)
    call write_to_file_grid_ST_b( region, region%ismip_output%yvelsurf, region%ice%v_surf_b / sec_per_year)
    call write_to_file_grid_ST_a( region, region%ismip_output%zvelsurf, region%ice%w_surf   / sec_per_year)
    call write_to_file_grid_ST_b( region, region%ismip_output%xvelbase, region%ice%u_base_b / sec_per_year)
    call write_to_file_grid_ST_b( region, region%ismip_output%yvelbase, region%ice%v_base_b / sec_per_year)
    call write_to_file_grid_ST_a( region, region%ismip_output%zvelbase, region%ice%w_base   / sec_per_year)
    call write_to_file_grid_ST_b( region, region%ismip_output%xvelmean, region%ice%u_vav_b  / sec_per_year)
    call write_to_file_grid_ST_b( region, region%ismip_output%yvelmean, region%ice%v_vav_b  / sec_per_year)

    ! Temperatures
    call write_to_file_grid_ST_a( region, region%ismip_output%litemptop, region%ice%Ti( :, 1))
    call write_to_single_ISMIP_regional_output_file( region, region%ismip_output%litempavg)
    call write_to_single_ISMIP_regional_output_file( region, region%ismip_output%litempgradgr)
    call write_to_single_ISMIP_regional_output_file( region, region%ismip_output%litempgradfl)
    call write_to_single_ISMIP_regional_output_file( region, region%ismip_output%litempbotgr)
    call write_to_single_ISMIP_regional_output_file( region, region%ismip_output%litempbotfl)

    ! Basal drag
    call write_to_file_grid_ST_a( region, region%ismip_output%strbasemag, region%ice%basal_shear_stress)

    ! Lateral mass balance
    call write_to_single_ISMIP_regional_output_file( region, region%ismip_output%licalvf)
    call write_to_single_ISMIP_regional_output_file( region, region%ismip_output%lifmassbf)

    ! Area fractions
    call write_to_file_grid_ST_a( region, region%ismip_output%sftgif, region%ice%fraction_margin, vmin=0._dp, vmax=1._dp)
    call write_to_file_grid_ST_a( region, region%ismip_output%sftgrf, region%ice%fraction_gr, vmin=0._dp, vmax=1._dp)
    call write_to_file_grid_ST_a( region, region%ismip_output%sftflf, region%ice%fraction_margin - region%ice%fraction_gr, vmin=0._dp, vmax=1._dp)

    ! === Scalars ===

    ! State with provided inputfields and optional masks
    call write_to_file_scalar( region, region%ismip_output%lim, region%ice%Hi * ice_density)
    call write_to_file_scalar( region, region%ismip_output%limnsw, region%ice%TAF * ice_density, mask=mask_ice_a)
    call write_to_file_scalar( region, region%ismip_output%iareagr, region%ice%fraction_gr, mask=mask_ice_a)
    call write_to_file_scalar( region, region%ismip_output%iareafl, 1._dp-region%ice%fraction_gr, mask=mask_ice_a)

    ! Fluxes with provided initial scalar values in Gt/yr
    call write_to_file_scalar( region, region%ismip_output%tendacabf, region%scalars%SMB_gr + region%scalars%SMB_fl)
    call write_to_file_scalar( region, region%ismip_output%tendlibmassbfgr, region%scalars%BMB_gr)
    call write_to_file_scalar( region, region%ismip_output%tendlibmassbffl, region%scalars%BMB_fl)
    call write_to_file_scalar( region, region%ismip_output%tendlicalvf, region%scalars%margin_ocean_flux)
    call write_to_file_scalar( region, region%ismip_output%tendlifmassbf, 0._dp)
    call write_to_file_scalar( region, region%ismip_output%tendligroundf, region%scalars%gl_flux)

    ! Set previous time step to current
    region%ismip_output%t_prev = region%ismip_output%t_curr

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_ISMIP_regional_output_files

  subroutine write_to_single_ISMIP_regional_output_file_grid( region, field)
    !< Write to single ISMIP regional output NetCDF file

    ! In/output variables:
    type(type_model_region),           intent(inout) :: region
    type(type_ismip_gridded_field),    intent(inout) :: field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_single_ISMIP_regional_output_file_grid'
    integer                        :: ncid
    character(len=16)              :: nt_str

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( field%filename, ncid)

    ! write the time to the file
    select case (field%fieldtype)
      case default
        call crash('invalid fieldtype for CFtime writing "' // trim(field%fieldtype) //  '"')
      case ('FL')
        call write_cftime_to_file( field%filename, ncid, region%time, with_bounds = .true.)
      case ('ST')
        call write_cftime_to_file( field%filename, ncid, region%time, with_bounds = .false.)
    end select

    ! Update the time counter attribute
    field%nt = field%nt + 1
    ! Convert counter to string
    write(nt_str, '(I16)') field%nt
    nt_str = adjustl(nt_str)
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'nt', trim(nt_str))

    ! write the data to the file
    call write_to_ISMIP_regional_output_file_field_grid( region, field, ncid)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_single_ISMIP_regional_output_file_grid

  subroutine write_to_file_grid_FL( region, field)
    !< Write to single ISMIP regional output NetCDF file

    ! In/output variables:
    type(type_model_region),                              intent(inout) :: region
    type(type_ismip_gridded_field),                       intent(inout) :: field

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_file_grid_FL'
    integer                               :: ncid
    character(len=16)                     :: nt_str
    real(dp)                              :: deltat
    real(dp), dimension(:),   allocatable :: d_grid_vec_partial_2D, d_mesh_vec_partial_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( field%filename, ncid)

    ! write the time to the file
    call write_cftime_to_file( field%filename, ncid, region%time, with_bounds = .true.)

    ! Update the time counter attribute
    field%nt = field%nt + 1
    ! Convert counter to string
    write(nt_str, '(I16)') field%nt
    nt_str = adjustl(nt_str)
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'nt', trim(nt_str))

    ! Determine deltat since previous writing (should be equal to 0 (first time step) or dt_ismip_output)
    deltat = region%ismip_output%t_curr - region%ismip_output%t_prev

    ! Allocate memory
    allocate( d_grid_vec_partial_2D( region%output_grid%n_loc ))
    allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))

    ! Determine values
    if (field%is_initial) then
      ! First timeframe, no accumulation yet. Following protocol, using snapshot field instead
      d_mesh_vec_partial_2D = field%accum
      field%is_initial = .false.
    else
      d_mesh_vec_partial_2D = field%accum / deltat
    end if

    ! Restore accumulation to 0
    field%accum( region%mesh%vi1: region%mesh%vi2) = 0._dp

    ! Map from mesh to grid
    call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)

    ! Write gridded field to file
    call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)

    ! Clean up memory
    deallocate( d_grid_vec_partial_2D)
    deallocate( d_mesh_vec_partial_2D)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_file_grid_FL

  subroutine write_to_file_grid_ST_a( region, field, inputfield, vmin, vmax)
    !< Write STATE gridded mesh field to single ISMIP regional output NetCDF file

    ! In/output variables:
    type(type_model_region),                              intent(inout) :: region
    type(type_ismip_gridded_field),                       intent(inout) :: field
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2), intent(in   ) :: inputfield
    real(dp), optional,                                   intent(in   ) :: vmin
    real(dp), optional,                                   intent(in   ) :: vmax

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_file_grid_ST_a'
    integer                               :: ncid
    character(len=16)                     :: nt_str
    real(dp)                              :: deltat
    real(dp), dimension(:),   allocatable :: d_grid_vec_partial_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( field%filename, ncid)

    ! write the time to the file
    call write_cftime_to_file( field%filename, ncid, region%time, with_bounds = .false.)

    ! Update the time counter attribute
    field%nt = field%nt + 1
    ! Convert counter to string
    write(nt_str, '(I16)') field%nt
    nt_str = adjustl(nt_str)
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'nt', trim(nt_str))

    ! Determine deltat since previous writing (should be equal to 0 (first time step) or dt_ismip_output)
    deltat = region%ismip_output%t_curr - region%ismip_output%t_prev

    ! Allocate memory
    allocate( d_grid_vec_partial_2D( region%output_grid%n_loc ))

    ! Map from mesh to grid
    call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, inputfield, d_grid_vec_partial_2D)

    ! Enforce bounds
    if (present( vmin)) then
      d_grid_vec_partial_2D = max(vmin, d_grid_vec_partial_2D)
    end if
    if (present( vmax)) then
      d_grid_vec_partial_2D = min(vmax, d_grid_vec_partial_2D)
    end if

    ! Write gridded field to file
    call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)

    ! Clean up memory
    deallocate( d_grid_vec_partial_2D)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_file_grid_ST_a

  subroutine write_to_file_grid_ST_b( region, field, inputfield)
    !< Write STATE gridded mesh field to single ISMIP regional output NetCDF file

    ! In/output variables:
    type(type_model_region),                              intent(inout) :: region
    type(type_ismip_gridded_field),                       intent(inout) :: field
    real(dp), dimension(region%mesh%ti1:region%mesh%ti2), intent(in   ) :: inputfield

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_file_grid_ST_b'
    integer                               :: ncid
    character(len=16)                     :: nt_str
    real(dp)                              :: deltat
    real(dp), dimension(:),   allocatable :: d_grid_vec_partial_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( field%filename, ncid)

    ! write the time to the file
    call write_cftime_to_file( field%filename, ncid, region%time, with_bounds = .false.)

    ! Update the time counter attribute
    field%nt = field%nt + 1
    ! Convert counter to string
    write(nt_str, '(I16)') field%nt
    nt_str = adjustl(nt_str)
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'nt', trim(nt_str))

    ! Determine deltat since previous writing (should be equal to 0 (first time step) or dt_ismip_output)
    deltat = region%ismip_output%t_curr - region%ismip_output%t_prev

    ! Allocate memory
    allocate( d_grid_vec_partial_2D( region%output_grid%n_loc ))

    ! Map from mesh to grid
    call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, inputfield, d_grid_vec_partial_2D)

    ! Write gridded field to file
    call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)

    ! Clean up memory
    deallocate( d_grid_vec_partial_2D)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_file_grid_ST_b

  subroutine write_to_file_scalar_ST( region, scalar, inputfield, mask)
    !< Write to STATE scalar to single ISMIP regional output NetCDF file

    ! In/output variables:
    type(type_model_region),                                        intent(inout) :: region
    type(type_ismip_scalar),                                        intent(inout) :: scalar
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2),           intent(in   ) :: inputfield
    logical,  dimension(region%mesh%vi1:region%mesh%vi2), optional, intent(in   ) :: mask

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_file_scalar_ST'
    integer                        :: ncid
    character(len=16)              :: nt_str
    real(dp)                       :: deltat
    real(dp)                       :: scalar_loc
    integer                        :: vi, ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Check wheter scalar is state
    if (scalar%fieldtype == 'FL') call crash('Flux variable should get initval')

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( scalar%filename, ncid)

    ! write the time to the file
    call write_cftime_to_file( scalar%filename, ncid, region%time, with_bounds = .false.)

    ! Update the time counter attribute
    scalar%nt = scalar%nt + 1
    ! Convert counter to string
    write(nt_str, '(I16)') scalar%nt
    nt_str = adjustl(nt_str)
    call add_attribute_char( scalar%filename, ncid, NF90_GLOBAL, 'nt', trim(nt_str))

    ! Determine deltat since previous writing (should be equal to 0 (first time step) or dt_ismip_output)
    deltat = region%ismip_output%t_curr - region%ismip_output%t_prev

    ! Accumulate scalar values per process
    scalar_loc = 0._dp
    do vi = region%mesh%vi1, region%mesh%vi2
      if (.not. present(mask) .or. mask( vi)) then
        ! Add value if no mask is provided, or if the mask is true
        scalar_loc = scalar_loc + inputfield( vi) * region%mesh%A( vi)
      end if
    end do

    ! Accumulate sum over all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, scalar_loc, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Write scalar to file
    call write_to_field_multopt_dp_0D( scalar%filename, ncid, scalar%name, scalar_loc)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_file_scalar_ST

  subroutine write_to_file_scalar_FL( region, scalar, initval)
    !< Write to FLUX scalar to single ISMIP regional output NetCDF file

    ! In/output variables:
    type(type_model_region),  intent(inout) :: region
    type(type_ismip_scalar),  intent(inout) :: scalar
    real(dp),                 intent(in   ) :: initval ! Gt/y

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_file_scalar_FL'
    integer                        :: ncid
    character(len=16)              :: nt_str
    real(dp)                       :: deltat
    real(dp)                       :: scalar_loc
    integer                        :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Check wheter scalar is flux
    if (scalar%fieldtype == 'ST') call crash('State variable should get inputfield')

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( scalar%filename, ncid)

    call write_cftime_to_file( scalar%filename, ncid, region%time, with_bounds = .true.)

    ! Update the time counter attribute
    scalar%nt = scalar%nt + 1
    ! Convert counter to string
    write(nt_str, '(I16)') scalar%nt
    nt_str = adjustl(nt_str)
    call add_attribute_char( scalar%filename, ncid, NF90_GLOBAL, 'nt', trim(nt_str))

    ! Determine deltat since previous writing (should be equal to 0 (first time step) or dt_ismip_output)
    deltat = region%ismip_output%t_curr - region%ismip_output%t_prev

    ! Determine scalar value
    if (scalar%is_initial) then
      ! First trimeframe, no accumulation yet, provide initial value instead
      scalar_loc = initval * 1.e12_dp / sec_per_year
      scalar%is_initial = .false.
    else
      scalar_loc = scalar%accum / deltat
      scalar%accum = 0._dp
    end if

    ! Write scalar to file
    call write_to_field_multopt_dp_0D( scalar%filename, ncid, scalar%name, scalar_loc)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_file_scalar_FL

  subroutine write_to_ISMIP_regional_output_file_field_grid( region, field, ncid)
    !< Write a single field to the ISMIP regional output NetCDF file

    ! In/output variables:
    type(type_model_region),           intent(inout) :: region
    type(type_ismip_gridded_field),    intent(inout) :: field
    integer,                           intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_ISMIP_regional_output_file_field_grid'
    real(dp), dimension(:),   allocatable :: d_grid_vec_partial_2D, d_mesh_vec_partial_2D
    real(dp), dimension(:),   allocatable :: dTdzeta
    real(dp), dimension(C%nz)             :: d_zeta_temp
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: SMB_loc, calving_flux, gl_flux
    integer                               :: vi
    real(dp)                              :: deltat

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine deltat since previous writing (should be equal to 0 (first time step) or dt_ismip_output)
    deltat = region%ismip_output%t_curr - region%ismip_output%t_prev

    ! Allocate memory
    allocate( d_grid_vec_partial_2D( region%output_grid%n_loc ))

    ! Add the specified data field to the file
    select case (field%name)
      case default
        call crash('unknown choice_output_field "' // trim( field%name) // '"')

      ! Geothermal heat flux (FL)
      case ('hfgeoubed')
        ! This is always a snapshot
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        ! First timeframe, no accumulation yet. Following protocol, using snapshot field instead
        d_mesh_vec_partial_2D = region%ice%geothermal_heat_flux / sec_per_year
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
        deallocate( d_mesh_vec_partial_2D)

      ! Thickness tendency (FL)
      case ('dlithkdt')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        if (field%is_initial) then
          ! First timeframe, spit out 0
          d_mesh_vec_partial_2D( :) = 0._dp
          field%is_initial = .false.
        else
          d_mesh_vec_partial_2D = (region%ice%Hi - field%accum) / (region%time - region%ismip_output%t_prev) / sec_per_year
          field%accum = region%ice%Hi
        end if
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
        deallocate( d_mesh_vec_partial_2D)

      ! Temperatures
      case ('litempavg')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        do vi = region%mesh%vi1, region%mesh%vi2
          if (region%ice%Hi( vi) > 0._dp) then
            d_zeta_temp = region%ice%Ti( vi, :)
            d_mesh_vec_partial_2D( vi) = vertical_average( region%mesh%zeta, d_zeta_temp)
          else
            d_mesh_vec_partial_2D( vi) = NaN
          end if
        end do
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
        deallocate( d_mesh_vec_partial_2D)

      case ('litempgradgr')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        allocate( dTdzeta( 1:region%mesh%nz))
        do vi = region%mesh%vi1, region%mesh%vi2
          if (region%ice%mask_grounded_ice( vi)) then
            d_zeta_temp = region%ice%Ti( vi, :)
            call multiply_CSR_matrix_with_vector_local( region%mesh%M_ddzeta_k_k_1D, d_zeta_temp, dTdzeta)
            d_mesh_vec_partial_2D( vi) = -1._dp / region%ice%Hi( vi) * dTdzeta( region%mesh%nz)
          else
            d_mesh_vec_partial_2D( vi) = NaN
          end if
        end do
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
        deallocate( d_mesh_vec_partial_2D)
        deallocate( dTdzeta)

      case ('litempgradfl')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        allocate( dTdzeta( 1:region%mesh%nz))
        do vi = region%mesh%vi1, region%mesh%vi2
          if (region%ice%mask_floating_ice( vi)) then
            call multiply_CSR_matrix_with_vector_local( region%mesh%M_ddzeta_k_k_1D, region%ice%Ti( vi, :), dTdzeta)
            d_mesh_vec_partial_2D( vi) = -1._dp / region%ice%Hi( vi) * dTdzeta( region%mesh%nz)
          else
            d_mesh_vec_partial_2D( vi) = NaN
          end if
        end do
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
        deallocate( d_mesh_vec_partial_2D)
        deallocate( dTdzeta)

      case ('litempbotgr')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        do vi = region%mesh%vi1, region%mesh%vi2
          if (region%ice%mask_grounded_ice( vi)) then
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

      ! Lateral mass balance
      case ('licalvf')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        if (field%is_initial) then
          ! First timeframe, no accumulation yet. Following protocol, using snapshot field instead
          call calc_ISMIP_fluxes( region%mesh, region%ice, calving_flux, gl_flux)
          d_mesh_vec_partial_2D = calving_flux * ice_density / sec_per_year
          field%is_initial = .false.
        else
          d_mesh_vec_partial_2D = field%accum / deltat
          field%accum( region%mesh%vi1: region%mesh%vi2) = 0._dp
        end if
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
        deallocate( d_mesh_vec_partial_2D)

      case ('lifmassbf')
        ! Undefined, so just write out zeros
        d_grid_vec_partial_2D( :) = 0._dp
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)

    end select

    ! Clean up memory
    deallocate( d_grid_vec_partial_2D)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_ISMIP_regional_output_file_field_grid

  subroutine create_ISMIP_regional_output_files( region)
    ! MAIN creation
    !< Create all ISMIP regional output NetCDF files

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_ISMIP_regional_output_files'

    ! Add routine to path
    call init_routine( routine_name)

    ! if no ISMIP output should be created, do nothing
    if (.not. C%do_create_ismip_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Initialise ISMIP_output
    call initialise_ismip_output( region)

    ! Create folder
    if (par%primary) call system('mkdir ' // trim(region%ismip_output%folder))

    ! Create all grid files

    ! Basic topography
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%lithk)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%orog)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%topg)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%base)

    ! Geothermal heat flux
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%hfgeoubed)

    ! Surface and basal mass balance
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%acabf)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%libmassbfgr)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%libmassbffl)

    ! Thickness tendency
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%dlithkdt)

    ! Velocities
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%xvelsurf)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%yvelsurf)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%zvelsurf)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%xvelbase)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%yvelbase)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%zvelbase)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%xvelmean)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%yvelmean)

    ! Temperatures
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%litemptop)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%litempavg)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%litempgradgr)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%litempgradfl)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%litempbotgr)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%litempbotfl)

    ! Basal drag
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%strbasemag)

    ! Lateral mass balance
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%licalvf)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%lifmassbf)

    ! Area fractions
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%sftgif)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%sftgrf)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%sftflf)

    ! === Scalars ===

    ! State
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%lim)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%limnsw)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%iareagr)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%iareafl)

    ! Fluxes
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%tendacabf)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%tendlibmassbfgr)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%tendlibmassbffl)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%tendlicalvf)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%tendlifmassbf)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%tendligroundf)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_ISMIP_regional_output_files

  subroutine create_single_ISMIP_regional_output_file_grid( ismip_output, field)
    !< Create a single ISMIP regional output NetCDF file - grid version

    ! In/output variables:
    type(type_ismip_output),      intent(inout) :: ismip_output
    type(type_ismip_gridded_field),    intent(in   ) :: field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_single_ISMIP_regional_output_file_grid'
    integer                        :: ncid
    character(len=*), parameter    :: iprecision = 'single'
    logical, parameter             :: do_compress = .true.
    integer                        :: res_int
    character(len=16)              :: res_str
    character(len=1024)            :: title

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating ISMIP grid output file "' // &
      UPSY%stru%colour_string( trim( field%filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( field%filename, ncid)

    ! Set up the grid in the file
    call setup_xy_grid_in_netcdf_file( field%filename, ncid, ismip_output%grid, do_include_lonlat = .false.)

    ! Add time dimension to the file
    select case (field%fieldtype)
      case default
        call crash('invalid fieldtype for CFtime writing "' // trim(field%fieldtype) //  '"')
      case ('FL')
        call add_cftime_dimension_to_file( field%filename, ncid, with_bounds = .true.)
      case ('ST')
        call add_cftime_dimension_to_file( field%filename, ncid, with_bounds = .false.)
    end select

    ! Add the field
    call add_field_grid_dp_2D( field%filename, ncid, field%name, precision = iprecision, &
      do_compress = do_compress, long_name = field%long_name, units = field%units, &
      standard_name = field%standard_name)

    ! Convert grid resolution and start/end times to string
    res_int = int(C%dx_output_grid_ANT)
    write(res_str, '(I16)') res_int
    res_str = adjustl(res_str)

    ! Make a title
    title = 'ISMIP7 output by UFEMISM - ' // trim(field%name)

    ! Add attributes
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'title', trim(title))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'Conventions', trim(C%ismip_conventions))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'grid_type', trim(ismip_output%IS_name))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'grid_resolution', trim(res_str) // 'm')
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'group', trim(C%ismip_group_name))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'model', trim(C%ismip_model_name))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'scenario', trim(C%ismip_scenario_name))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'contact_name', trim(C%ismip_contact_name))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'contact_email', trim(C%ismip_contact_email))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'crs', trim(ismip_output%crs))

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_single_ISMIP_regional_output_file_grid

  subroutine create_single_ISMIP_regional_output_file_scalar( ismip_output, scalar)
    !< Create a single ISMIP regional output NetCDF file - grid version

    ! In/output variables:
    type(type_ismip_output),      intent(inout) :: ismip_output
    type(type_ismip_scalar),      intent(in   ) :: scalar

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_single_ISMIP_regional_output_file_scalar'
    integer                        :: ncid
    character(len=*), parameter    :: iprecision = 'single'
    logical, parameter             :: do_compress = .true.
    integer                        :: res_int
    character(len=16)              :: res_str
    character(len=1024)            :: title

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating ISMIP grid output file "' // &
      UPSY%stru%colour_string( trim( scalar%filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( scalar%filename, ncid)

    ! Add time dimension to the file
    select case (scalar%fieldtype)
      case default
        call crash('invalid fieldtype for CFtime writing "' // trim(scalar%fieldtype) //  '"')
      case ('FL')
        call add_cftime_dimension_to_file( scalar%filename, ncid, with_bounds = .true.)
      case ('ST')
        call add_cftime_dimension_to_file( scalar%filename, ncid, with_bounds = .false.)
    end select

    ! Add the field
    call add_field_dp_0D( scalar%filename, ncid, scalar%name, precision = iprecision, &
      do_compress = do_compress, long_name = scalar%long_name, units = scalar%units, &
      standard_name = scalar%standard_name)

    ! Convert grid resolution and start/end times to string
    res_int = int(C%dx_output_grid_ANT)
    write(res_str, '(I16)') res_int
    res_str = adjustl(res_str)

    ! Make a title
    title = 'ISMIP7 output by UFEMISM - ' // trim(scalar%name)

    ! Add attributes
    call add_attribute_char( scalar%filename, ncid, NF90_GLOBAL, 'title', trim(title))
    call add_attribute_char( scalar%filename, ncid, NF90_GLOBAL, 'Conventions', trim(C%ismip_conventions))
    call add_attribute_char( scalar%filename, ncid, NF90_GLOBAL, 'grid_type', trim(ismip_output%IS_name))
    call add_attribute_char( scalar%filename, ncid, NF90_GLOBAL, 'grid_resolution', trim(res_str) // 'm')
    call add_attribute_char( scalar%filename, ncid, NF90_GLOBAL, 'group', trim(C%ismip_group_name))
    call add_attribute_char( scalar%filename, ncid, NF90_GLOBAL, 'model', trim(C%ismip_model_name))
    call add_attribute_char( scalar%filename, ncid, NF90_GLOBAL, 'scenario', trim(C%ismip_scenario_name))
    call add_attribute_char( scalar%filename, ncid, NF90_GLOBAL, 'contact_name', trim(C%ismip_contact_name))
    call add_attribute_char( scalar%filename, ncid, NF90_GLOBAL, 'contact_email', trim(C%ismip_contact_email))
    call add_attribute_char( scalar%filename, ncid, NF90_GLOBAL, 'crs', trim(ismip_output%crs))

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_single_ISMIP_regional_output_file_scalar

  subroutine initialise_ISMIP_output( region)
    ! Initialise the ISMIP grid type

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter                       :: routine_name = 'initialise_ISMIP_output'
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: SMB_loc, BMB_gr_loc, BMB_fl_loc
    integer                                              :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Convert region name to ISMIP ice sheet name
    select case (region%name)
      case default
        call crash('invalid region name for ISMIP output "' // trim( region%name) // '"')
      case ('ANT')
        region%ismip_output%IS_name = 'AIS'
        region%ismip_output%crs = 'EPSG:3031'
      case ('GRL')
        region%ismip_output%IS_name = 'GIS'
        region%ismip_output%crs = 'EPSG:3413'
    end select

    ! Copy the output grid
    region%ismip_output%grid = region%output_grid

    ! Initialise the time trackers
    region%ismip_output%t_prev = C%start_time_of_run
    region%ismip_output%t_curr = C%start_time_of_run

    ! Basic topography
    call initialise_ISMIP_field( region, region%ismip_output%lithk, 'lithk', 'Ice thickness', 'land_ice_thickness', 'm', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%orog,  'orog' , 'Surface elevation', 'surface_altitude', 'm', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%topg,  'topg' , 'Bedrock elevation', 'bedrock_altitude', 'm', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%base,  'base' , 'Ice base elevation', 'base_altitude', 'm', 'ST')

    ! Geothermal heat flux
    call initialise_ISMIP_field( region, region%ismip_output%hfgeoubed, 'hfgeoubed' , &
      'Geothermal heat flux', 'upward_geothermal_heat_flux_in_land_ice', 'W m-2', 'FL')

    ! Surface and basal mass balances
    SMB_loc( region%mesh%vi1: region%mesh%vi2) = region%SMB%SMB( region%mesh%vi1: region%mesh%vi2)
    call initialise_ISMIP_field( region, region%ismip_output%acabf, 'acabf' , &
      'Surface mass balance flux', 'land_ice_surface_specific_mass_balance_flux', 'kg m-2 s-1', 'FL', &
      initfield = SMB_loc * ice_density / sec_per_year)

    BMB_gr_loc( :) = 0._dp
    BMB_fl_loc( :) = 0._dp
    do vi = region%mesh%vi1, region%mesh%vi2
      if (region%ice%mask_grounded_ice( vi)) BMB_gr_loc( vi) = region%BMB%BMB( vi)
      if (region%ice%mask_floating_ice( vi)) BMB_fl_loc( vi) = region%BMB%BMB( vi)
    end do
    call initialise_ISMIP_field( region, region%ismip_output%libmassbfgr, 'libmassbfgr' , &
      'Basal mass balance flux beneath grounded ice', 'land_ice_basal_specific_mass_balance_flux', 'kg m-2 s-1', 'FL', &
      initfield = BMB_gr_loc * ice_density / sec_per_year)
    call initialise_ISMIP_field( region, region%ismip_output%libmassbffl, 'libmassbffl' , &
      'Basal mass balance flux beneath floating ice', 'land_ice_basal_specific_mass_balance_flux', 'kg m-2 s-1', 'FL', &
      initfield = BMB_fl_loc * ice_density / sec_per_year)

    ! Thickness tendency
    call initialise_ISMIP_field( region, region%ismip_output%dlithkdt, 'dlithkdt' , &
      'Ice thickness imbalance', 'tendency_of_land_ice_thickness', 'm s-1', 'FL')
    ! Use accum of thickness tendency to store previous ice thickness
    region%ismip_output%dlithkdt%accum = region%ice%Hi

    ! Velocities
    call initialise_ISMIP_field( region, region%ismip_output%xvelsurf, 'xvelsurf' , &
      'Surface velocity in x', 'land_ice_surface_x_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%yvelsurf, 'yvelsurf' , &
      'Surface velocity in y', 'land_ice_surface_y_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%zvelsurf, 'zvelsurf' , &
      'Surface velocity in z', 'land_ice_surface_upward_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%xvelbase, 'xvelbase' , &
      'Basal velocity in x', 'land_ice_basal_x_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%yvelbase, 'yvelbase' , &
      'Basal velocity in y', 'land_ice_basal_y_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%zvelbase, 'zvelbase' , &
      'Basal velocity in z', 'land_ice_basal_upward_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%xvelmean, 'xvelmean' , &
      'Mean velocity in x', 'land_ice_vertical_mean_x_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%yvelmean, 'yvelmean' , &
      'Mean velocity in y', 'land_ice_vertical_mean_y_velocity', 'm s-1', 'ST')

    ! Temperatures
    call initialise_ISMIP_field( region, region%ismip_output%litemptop, 'litemptop' , &
      'Surface temperature', 'temperature_at_top_of_ice_sheet_model', 'K', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%litempavg, 'litempavg' , &
      'Depth average temperature', 'land_ice_temperature', 'K', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%litempgradgr, 'litempgradgr' , &
      'Vertical Basal temperature gradient beneath grounded ice sheet', &
      'temperature_gradient_at_base_of_ice_sheet_model', 'K m-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%litempgradfl, 'litempgradfl' , &
      'Vertical Basal temperature gradient beneath floating ice sheet', &
      'temperature_gradient_at_base_of_ice_sheet_model', 'K m-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%litempbotgr, 'litempbotgr' , &
      'Basal temperature beneath grounded ice', 'temperature_at_base_of_ice_sheet_model', 'K', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%litempbotfl, 'litempbotfl' , &
      'Basal temperature beneath floating ice', 'temperature_at_base_of_ice_sheet_model', 'K', 'ST')

    ! Basal drag
    call initialise_ISMIP_field( region, region%ismip_output%strbasemag, 'strbasemag' , &
      'Basal drag', 'land_ice_basal_drag', 'Pa', 'ST')

    ! Lateral mass balance
    call initialise_ISMIP_field( region, region%ismip_output%licalvf, 'licalvf' , &
      'Calving flux', 'land_ice_specific_mass_flux_due_to_calving', 'kg m-2 s-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_output%lifmassbf,   'lifmassbf' , &
      'Ice front melt flux', 'land_ice_specific_mass_flux_due_to_ice_front_melting', 'kg m-2 s-1', 'FL')

    ! Area fractions
    call initialise_ISMIP_field( region, region%ismip_output%sftgif, 'sftgif' , &
      'Land ice area fraction', 'land_ice_area_fraction', '1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%sftgrf, 'sftgrf' , &
      'Grounded ice sheet area fraction', 'grounded_ice_sheet_area_fraction', '1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%sftflf, 'sftflf' , &
      'Floating ice sheet area fraction', 'floating_ice_shelf_area_fraction', '1', 'ST')

    ! === Scalars ===

    ! State
    call initialise_ISMIP_field( region, region%ismip_output%lim, 'lim', &
      'Total ice mass', 'land_ice_mass', 'kg', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%limnsw, 'limnsw', &
      'Mass above floatation', 'land_ice_mass_not_displacing_sea_water', 'kg', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%iareagr, 'iareagr', &
      'Grounded ice area', 'grounded_ice_sheet_area', 'm^2', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%iareafl, 'iareafl', &
      'Floating ice area', 'floating_ice_shelf_area', 'm^2', 'ST')

    ! Fluxes
    call initialise_ISMIP_field( region, region%ismip_output%tendacabf, 'tendacabf', &
      'Total SMB flux', 'tendency_of_land_ice_mass_due_to_surface_mass_balance', 'kg s^-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_output%tendlibmassbfgr, 'tendlibmassbfgr', &
      'Total BMB flux beneath grounded ice', 'tendency_of_land_ice_mass_due_to_basal_mass_balance', 'kg s^-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_output%tendlibmassbffl, 'tendlibmassbffl', &
      'Total BMB flux beneath floating ice', 'tendency_of_land_ice_mass_due_to_basal_mass_balance', 'kg s^-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_output%tendlicalvf, 'tendlicalvf', &
      'Total calving flux', 'tendency_of_land_ice_mass_due_to_calving', 'kg s^-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_output%tendlifmassbf, 'tendlifmassbf', &
      'Total ice front melting flux', 'tendency_of_land_ice_mass_due_to_ice_front_melting', 'kg s^-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_output%tendligroundf, 'tendligroundf', &
      'Total grounding line flux', 'tbd', 'kg s^-1', 'FL')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ISMIP_output

  subroutine initialise_ISMIP_field_grid( region, field, name, long_name, standard_name, units, fieldtype, initfield)
    ! Initialise a single field

    type(type_model_region),                                        intent(inout) :: region
    type(type_ismip_gridded_field),                                 intent(inout) :: field
    character(len=*),                                               intent(in   ) :: name
    character(len=*),                                               intent(in   ) :: long_name
    character(len=*),                                               intent(in   ) :: standard_name
    character(len=*),                                               intent(in   ) :: units
    character(len=2),                                               intent(in   ) :: fieldtype
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2), optional, intent(in   ) :: initfield

    ! Local variables
    character(len=1024), parameter :: routine_name = 'initialise_ISMIP_field_grid'
    character(len=4)               :: start_year, end_year

    ! Add routine to path
    call init_routine( routine_name)

    ! Inherit metadata
    field%name          = name
    field%long_name     = long_name
    field%standard_name = standard_name
    field%units         = units
    field%fieldtype     = fieldtype

    ! Convert grid resolution and start/end times to string
    write(start_year, '(I4)') int(C%start_time_of_run)
    write(end_year, '(I4)') int(C%end_time_of_run)

    ! Define the name of the subfolder
    region%ismip_output%folder = trim( C%output_dir) // 'CORE' // '/'

    ! Define the filename for this field
    field%filename = trim( region%ismip_output%folder) // trim(field%name) // '_' // &
      trim(region%ismip_output%IS_name) // '_' // trim(C%ismip_group_name) // '_' // &
      trim(C%ismip_model_name) // '_' // trim(C%ismip_member_id) // '_' // &
      trim(C%ismip_esm_name) // '_' // trim(C%ismip_forcing_member_id) // '_' // &
      trim(C%ismip_scenario_name) // '_' // trim(C%ismip_counter) // '_' // &
      trim(start_year) // '-' // trim(end_year) // '.nc'

    ! Allocate fields for accumulation during each timestep for flux fields
    if (field%fieldtype == 'FL') then
      allocate(field%accum( region%mesh%vi1: region%mesh%vi2), source = 0._dp)
      field%is_initial = .true.
      if (present(initfield)) then
        field%accum( region%mesh%vi1:region%mesh%vi2) = initfield( region%mesh%vi1:region%mesh%vi2)
      end if
    end if

    ! Initialise the counter for number of time values
    field%nt = 0

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ISMIP_field_grid

  subroutine initialise_ISMIP_field_scalar( region, scalar, name, long_name, standard_name, units, fieldtype)
    ! Initialise a single field

    type(type_model_region),           intent(inout) :: region
    type(type_ismip_scalar),           intent(inout) :: scalar
    character(len=*),                  intent(in   ) :: name
    character(len=*),                  intent(in   ) :: long_name
    character(len=*),                  intent(in   ) :: standard_name
    character(len=*),                  intent(in   ) :: units
    character(len=2),                  intent(in   ) :: fieldtype

    ! Local variables
    character(len=1024), parameter :: routine_name = 'initialise_ISMIP_field_scalar'
    character(len=4)               :: start_year, end_year

    ! Add routine to path
    call init_routine( routine_name)

    ! Inherit metadata
    scalar%name          = name
    scalar%long_name     = long_name
    scalar%standard_name = standard_name
    scalar%units         = units
    scalar%fieldtype     = fieldtype

    ! Convert start/end times to string
    write(start_year, '(I4)') int(C%start_time_of_run)
    write(end_year, '(I4)') int(C%end_time_of_run)

    ! Define the name of the subfolder
    region%ismip_output%folder = trim( C%output_dir) // 'CORE' // '/'

    ! Define the filename for this field
    scalar%filename = trim( region%ismip_output%folder) // trim(scalar%name) // '_' // &
      trim(region%ismip_output%IS_name) // '_' // trim(C%ismip_group_name) // '_' // &
      trim(C%ismip_model_name) // '_' // trim(C%ismip_member_id) // '_' // &
      trim(C%ismip_esm_name) // '_' // trim(C%ismip_forcing_member_id) // '_' // &
      trim(C%ismip_scenario_name) // '_' // trim(C%ismip_counter) // '_' // &
      trim(start_year) // '-' // trim(end_year) // '.nc'

    ! Initialise
    if ((scalar%fieldtype == 'FL')) then
      scalar%accum = 0._dp
      scalar%is_initial = .true.
    end if

    ! Initialise the counter for number of time values
    scalar%nt = 0

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ISMIP_field_scalar

  subroutine remap_ISMIP_output( mesh_old, mesh_new, ice, ismip_output)
    ! Reallocate the accumulated fields and redefine Hi_prev

    type(type_mesh),                   intent(in   ) :: mesh_old
    type(type_mesh),                   intent(in   ) :: mesh_new
    type(type_ice_model),              intent(in   ) :: ice
    type(type_ismip_output),      intent(inout) :: ismip_output

    ! Local variables
    character(len=1024), parameter :: routine_name = 'remap_ISMIP_output'

    ! Add routine to path
    call init_routine( routine_name)

    call reallocate_bounds( ismip_output%acabf%accum, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( ismip_output%libmassbfgr%accum, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( ismip_output%libmassbffl%accum, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( ismip_output%licalvf%accum, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( ismip_output%lifmassbf%accum, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( ismip_output%dlithkdt%accum, mesh_new%vi1, mesh_new%vi2)

    ! Use accum to store current Hi for dHidt
    ismip_output%dlithkdt%accum = ice%Hi

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_ISMIP_output

end module ismip_output_files
