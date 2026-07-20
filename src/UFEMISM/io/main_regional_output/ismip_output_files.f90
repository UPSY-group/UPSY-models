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
  use netcdf, only: NF90_GLOBAL, NF90_FILL_DOUBLE
  use reallocate_mod, only: reallocate_bounds

  implicit none

  private

  public :: create_ISMIP_regional_output_files, write_to_ISMIP_regional_output_files, &
    accumulate_ISMIP_flux_fields, remap_ISMIP_output

  interface create_single_ISMIP_regional_output_file
    module procedure create_single_ISMIP_regional_output_file_grid
    module procedure create_single_ISMIP_regional_output_file_scalar
  end interface

  interface write_to_file
    module procedure write_to_file_grid_ST
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
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: SMB_loc, BMB_gr, BMB_fl
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

    ! Extract SMB only over ice-covered fraction
    do vi = region%mesh%vi1, region%mesh%vi2
      SMB_loc( vi) = SMB_loc( vi) * min(1._dp,max(0._dp,(region%ice%fraction_margin( vi))))
    end do

    ! Determine BMB_gr and BMB_fl, dependent on subgrid scheme
    do vi = region%mesh%vi1, region%mesh%vi2
      select case (C%choice_BMB_subgrid)
        case default
          call crash('unknown choice_BMB_subgrid "' // C%choice_BMB_subgrid // '"')
        case ('FCMP')
          if (region%ice%mask_floating_ice( vi) .or. region%ice%mask_gl_fl( vi)) then
            BMB_fl( vi) = region%BMB%BMB_shelf( vi)
            BMB_gr( vi) = 0._dp
          elseif (region%ice%mask_grounded_ice( vi) .or. region%ice%mask_gl_gr( vi)) then
            BMB_fl( vi) = 0._dp
            BMB_gr( vi) = region%BMB%BMB_sheet( vi)
          else
            BMB_fl( vi) = 0._dp
            BMB_gr( vi) = 0._dp
          end if
        case ('NMP')
          if (region%ice%mask_floating_ice( vi) .and. region%ice%fraction_gr( vi) == 0._dp) then
            BMB_fl( vi) = region%BMB%BMB_shelf( vi)
            BMB_gr( vi) = 0._dp
          elseif (region%ice%fraction_gr( vi) > 0._dp) then
            BMB_fl( vi) = 0._dp
            BMB_gr( vi) = region%BMB%BMB_sheet( vi)
          else
            BMB_fl( vi) = 0._dp
            BMB_gr( vi) = 0._dp
          end if
        case ('PMP')
          if (region%ice%mask_floating_ice( vi) .or. region%ice%mask_grounded_ice( vi)) then
            BMB_fl( vi) = (1._dp - region%ice%fraction_gr( vi)) * region%BMB%BMB_shelf( vi)
            BMB_gr( vi) = region%ice%fraction_gr( vi) * region%BMB%BMB_sheet( vi)
          else
            BMB_fl( vi) = 0._dp
            BMB_gr( vi) = 0._dp
          end if
      end select
    end do

    ! Get delta t since last current time
    deltat = region%time - region%ismip_output%t_curr

    ! Accumulate regular FL fields.
    call accumulate_single_ISMIP_flux_field( region, SMB_loc, deltat, &
      region%ismip_output%tendacabf, region%ismip_output%acabf)
    call accumulate_single_ISMIP_flux_field( region, BMB_gr, deltat, &
      region%ismip_output%tendlibmassbfgr, region%ismip_output%libmassbfgr)
    call accumulate_single_ISMIP_flux_field( region, BMB_fl, deltat, &
      region%ismip_output%tendlibmassbffl, region%ismip_output%libmassbffl)
    call accumulate_single_ISMIP_flux_field( region, calving_flux, deltat, &
      region%ismip_output%tendlicalvf, region%ismip_output%licalvf)
    call accumulate_single_ISMIP_flux_field( region, -gl_flux, deltat, &
      region%ismip_output%tendligroundf, region%ismip_output%ligroundf)

    ! === Exceptions ===
    ! Geothermal heat flux: Store snapshot in accum
    region%ismip_output%hfgeoubed%accum( region%mesh%vi1:region%mesh%vi2) = &
      region%ice%geothermal_heat_flux( region%mesh%vi1:region%mesh%vi2) / sec_per_year

    ! Frontal melting: just spit out zeros throughout

    ! Ice thickness tendency: accumulation used to store previously written thickness

    ! Update current time
    region%ismip_output%t_curr = region%time

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine accumulate_ISMIP_flux_fields

  subroutine accumulate_single_ISMIP_flux_field( region, d_partial, deltat, scalar, field)
    !< Write to ISMIP regional output NetCDF files - grid version

    ! In/output variables:
    type(type_model_region),                              intent(inout) :: region
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2), intent(in   ) :: d_partial
    real(dp),                                             intent(in   ) :: deltat
    type(type_ismip_scalar),                              intent(inout) :: scalar
    type(type_ismip_gridded_field),                       intent(inout) :: field

    ! Local variables:
    character(len=1024), parameter                       :: routine_name = 'accumulate_single_ISMIP_flux_field'
    integer                                              :: vi, ierr
    real(dp)                                             :: scalar_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise local scalar
    scalar_loc = 0._dp

    do vi = region%mesh%vi1, region%mesh%vi2
      ! Accumulate field
      field%accum( vi) = field%accum( vi) + d_partial( vi) * ice_density / sec_per_year * deltat
      ! Add to integrated scalar
      scalar_loc = scalar_loc + d_partial( vi) * ice_density * region%mesh%A( vi) / sec_per_year * deltat
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
    logical,  dimension(region%mesh%ti1:region%mesh%ti2) :: mask_gr_b
    real(dp),  dimension(region%mesh%vi1:region%mesh%vi2) :: T_vav, TF
    real(dp), dimension(C%nz)                            :: d_zeta_temp
    integer                                              :: vi, ti

    ! Add routine to path
    call init_routine( routine_name)

    ! if no ISMIP output should be created, do nothing
    if (.not. C%do_create_ismip_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Determine masks
    do vi = region%mesh%vi1, region%mesh%vi2
      if (region%ice%geom%Hi( vi) > 0._dp) then
        mask_ice_a( vi) = .true.
      else
        mask_ice_a( vi) = .false.
      end if
    end do

    do ti = region%mesh%ti1, region%mesh%ti2
      if (region%ice%fraction_gr_b( ti) > 0._dp) then
        mask_gr_b( ti) = .true.
      else
        mask_gr_b( ti) = .false.
      end if
    end do

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to ISMIP output files' // '...'

    ! Basic topography
    call write_to_file( region, region%ismip_output%lithk, inputfield_a=region%ice%geom%Hi, vmin=0._dp)
    call write_to_file( region, region%ismip_output%orog,  inputfield_a=region%ice%Hs, vmin=0._dp)
    call write_to_file( region, region%ismip_output%topg,  inputfield_a=region%ice%geom%Hb)
    call write_to_file( region, region%ismip_output%base,  inputfield_a=region%ice%Hib)

    ! Geothermal heat flux
    call write_to_file_grid_FL( region, region%ismip_output%hfgeoubed, vmin=0._dp)

    ! Surface and basal mass balance
    call write_to_file_grid_FL( region, region%ismip_output%acabf)
    call write_to_file_grid_FL( region, region%ismip_output%libmassbfgr)
    call write_to_file_grid_FL( region, region%ismip_output%libmassbffl)

    ! Thickness tendency
    do vi = region%mesh%vi1, region%mesh%vi2
      region%ismip_output%dlithkdt%accum( vi) = (region%ice%geom%Hi( vi) - region%ismip_output%dlithkdt%accum( vi)) / sec_per_year
    end do
    call write_to_file_grid_FL( region, region%ismip_output%dlithkdt)
    region%ismip_output%dlithkdt%accum( region%mesh%vi1:region%mesh%vi2) = region%ice%geom%Hi( region%mesh%vi1:region%mesh%vi2)

    ! Velocities
    call write_to_file( region, region%ismip_output%xvelsurf, inputfield_b=region%ice%u_surf_b / sec_per_year)
    call write_to_file( region, region%ismip_output%yvelsurf, inputfield_b=region%ice%v_surf_b / sec_per_year)
    call write_to_file( region, region%ismip_output%zvelsurf, inputfield_a=region%ice%w_3D( :, 1) / sec_per_year)
    call write_to_file( region, region%ismip_output%xvelbase, inputfield_b=region%ice%u_base_b / sec_per_year)
    call write_to_file( region, region%ismip_output%yvelbase, inputfield_b=region%ice%v_base_b / sec_per_year)
    call write_to_file( region, region%ismip_output%zvelbase, inputfield_a=region%ice%w_3D( :, C%nz) / sec_per_year)
    call write_to_file( region, region%ismip_output%xvelmean, inputfield_b=region%ice%u_vav_b  / sec_per_year)
    call write_to_file( region, region%ismip_output%yvelmean, inputfield_b=region%ice%v_vav_b  / sec_per_year)

    ! Temperatures
    call write_to_file( region, region%ismip_output%litemptop, inputfield_a=region%ice%Ti( :, 1), mask_a=mask_ice_a)

    do vi = region%mesh%vi1, region%mesh%vi2
      if (mask_ice_a( vi)) then
        ! Get vertical T profile
        d_zeta_temp = region%ice%Ti( vi, :)
        ! Compute vertical average
        T_vav( vi) = vertical_average( region%mesh%zeta, d_zeta_temp)
      else
        ! No ice here, return NaN
        T_vav( vi) = NaN
      end if
    end do
    call write_to_file( region, region%ismip_output%litempavg,    inputfield_a=T_vav)
    call write_to_file( region, region%ismip_output%litempbotgr,  inputfield_a=region%ice%Ti( :, C%nz), mask_a=region%ice%mask_grounded_ice)
    call write_to_file( region, region%ismip_output%litempbotfl,  inputfield_a=region%ice%Ti( :, C%nz), mask_a=region%ice%mask_floating_ice)

    ! Basal drag
    call write_to_file( region, region%ismip_output%strbasemag, inputfield_b=region%ice%basal_shear_stress, mask_b=mask_gr_b, vmin=0._dp)

    ! Lateral mass balance
    call write_to_file_grid_FL( region, region%ismip_output%licalvf, vmax=0._dp)
    call write_to_file_grid_FL( region, region%ismip_output%ligroundf)
    call write_to_file_grid_FL( region, region%ismip_output%lifmassbf)

    ! Area fractions
    call write_to_file( region, region%ismip_output%sftgif, inputfield_a=region%ice%fraction_margin, vmin=0._dp, vmax=1._dp)
    call write_to_file( region, region%ismip_output%sftgrf, inputfield_a=region%ice%fraction_gr, vmin=0._dp, vmax=1._dp)
    call write_to_file( region, region%ismip_output%sftflf, inputfield_a=region%ice%fraction_margin - region%ice%fraction_gr, &
      vmin=0._dp, vmax=1._dp)

    ! Other stuff
    do vi = region%mesh%vi1, region%mesh%vi2
      ! Determine thermal forcing
      if (C%choice_BMB_model_ANT == 'laddie') then
        TF( vi) = region%BMB%laddie%now%T( vi) - region%BMB%laddie%T_freeze( vi)
      else
        TF( vi) = region%ocean%T_draft( vi) - region%ocean%T_freezing_point( vi)
      end if
    end do
    call write_to_file( region, region%ismip_output%tfbase, inputfield_a=TF, mask_a=region%ice%mask_floating_ice)

    ! === Scalars ===

    ! State with provided inputfields and optional masks
    call write_to_file( region, region%ismip_output%lim, region%ice%geom%Hi * ice_density)
    call write_to_file( region, region%ismip_output%limnsw, max(0._dp,region%ice%TAF) * ice_density, mask=mask_ice_a)
    call write_to_file( region, region%ismip_output%iareagr, region%ice%fraction_gr, mask=mask_ice_a)
    call write_to_file( region, region%ismip_output%iareafl, 1._dp-region%ice%fraction_gr, mask=mask_ice_a)

    ! Fluxes
    call write_to_file( region, region%ismip_output%tendacabf)
    call write_to_file( region, region%ismip_output%tendlibmassbfgr)
    call write_to_file( region, region%ismip_output%tendlibmassbffl)
    call write_to_file( region, region%ismip_output%tendlicalvf)
    call write_to_file( region, region%ismip_output%tendlifmassbf)
    call write_to_file( region, region%ismip_output%tendligroundf)

    ! Set previous time step to current
    region%ismip_output%t_prev = region%ismip_output%t_curr

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_ISMIP_regional_output_files

  subroutine write_to_file_grid_FL( region, field, vmin, vmax)
    !< Write to single ISMIP regional output NetCDF file

    ! In/output variables:
    type(type_model_region),                              intent(inout) :: region
    type(type_ismip_gridded_field),                       intent(inout) :: field
    real(dp),                                   optional, intent(in   ) :: vmin
    real(dp),                                   optional, intent(in   ) :: vmax

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
    d_mesh_vec_partial_2D = field%accum / deltat

    ! Restore accumulation to 0 if required
    field%accum( region%mesh%vi1: region%mesh%vi2) = 0._dp

    ! Map from mesh to grid
    call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)

    ! Enforce bounds
    if (present( vmin)) then
      d_grid_vec_partial_2D = max(vmin, d_grid_vec_partial_2D)
    end if
    if (present( vmax)) then
      d_grid_vec_partial_2D = min(vmax, d_grid_vec_partial_2D)
    end if

    ! Convert nans to fill values
    where (isnan( d_grid_vec_partial_2D))
      d_grid_vec_partial_2D = NF90_FILL_DOUBLE
    end where

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

  subroutine write_to_file_grid_ST( region, field, inputfield_a, inputfield_b, vmin, vmax, mask_a, mask_b)
    !< Write STATE gridded mesh field to single ISMIP regional output NetCDF file

    ! In/output variables:
    type(type_model_region),                                        intent(inout) :: region
    type(type_ismip_gridded_field),                                 intent(inout) :: field
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2), optional, intent(in   ) :: inputfield_a
    real(dp), dimension(region%mesh%ti1:region%mesh%ti2), optional, intent(in   ) :: inputfield_b
    real(dp),                                             optional, intent(in   ) :: vmin
    real(dp),                                             optional, intent(in   ) :: vmax
    logical, dimension(region%mesh%vi1:region%mesh%vi2),  optional, intent(in   ) :: mask_a
    logical, dimension(region%mesh%ti1:region%mesh%ti2),  optional, intent(in   ) :: mask_b

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_file_grid_ST'
    integer                               :: ncid, vi, ti
    character(len=16)                     :: nt_str
    real(dp)                              :: deltat
    real(dp), dimension(:),  allocatable :: d_grid_vec_partial_2D, d_mesh_vec_partial_2D, mask_grid

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

    if (present(inputfield_b)) then
      allocate( d_mesh_vec_partial_2D( region%mesh%ti1:region%mesh%ti2))

      ! Apply mask if requested
      do ti = region%mesh%ti1, region%mesh%ti2
        if (.not. present(mask_b)) then
          ! Add value if no mask is provided, or if the mask is true
          d_mesh_vec_partial_2D( ti) = inputfield_b( ti)
        elseif (mask_b( ti)) then
          ! Add value if no mask is provided, or if the mask is true
          d_mesh_vec_partial_2D( ti) = inputfield_b( ti)
        else
          ! Only used for strbasemag. Set to 0 over floating regions, rather than masking
          d_mesh_vec_partial_2D( ti) = 0._dp
        end if
      end do

      ! Map from mesh triangles to grid
      call map_from_mesh_triangles_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
    elseif (present(inputfield_a)) then
      allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))

      ! Apply mask if requested
      do vi = region%mesh%vi1, region%mesh%vi2
        if (.not. present(mask_a)) then
          ! Add value if no mask is provided, or if the mask is true
          d_mesh_vec_partial_2D( vi) = inputfield_a( vi)
        elseif (mask_a( vi)) then
          ! Add value if no mask is provided, or if the mask is true
          d_mesh_vec_partial_2D( vi) = inputfield_a( vi)
        else
          d_mesh_vec_partial_2D( vi) = NaN
        end if
      end do

      ! Map from mesh vertices to grid
      call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
    else
      call crash('write_to_file_grid_ST requires either inputfield_a or inputfield_b')
    end if

    ! Enforce bounds
    if (present( vmin)) then
      d_grid_vec_partial_2D = max(vmin, d_grid_vec_partial_2D)
    end if
    if (present( vmax)) then
      d_grid_vec_partial_2D = min(vmax, d_grid_vec_partial_2D)
    end if

    ! For velocity fields, mask regions with ice fraction < 0.5
    if (present(inputfield_b)) then
      allocate( mask_grid( region%output_grid%n_loc ))
      call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%fraction_margin, mask_grid)
      where (mask_grid < 0.5_dp)
        d_grid_vec_partial_2D = NF90_FILL_DOUBLE
      end where
      deallocate( mask_grid)
    end if

    ! Convert nans to fill values
    where (isnan( d_grid_vec_partial_2D))
      d_grid_vec_partial_2D = NF90_FILL_DOUBLE
    end where

    ! Write gridded field to file
    call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)

    ! Clean up memory
    deallocate( d_grid_vec_partial_2D)
    if (allocated( d_mesh_vec_partial_2D)) deallocate( d_mesh_vec_partial_2D)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_file_grid_ST

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
    ! Apply mask if requested
    do vi = region%mesh%vi1, region%mesh%vi2
      if (.not. present(mask)) then
        ! Add value if no mask is provided, or if the mask is true
        scalar_loc = scalar_loc + inputfield( vi) * region%mesh%A( vi)
      elseif (mask( vi)) then
        ! Add value if no mask is provided, or if the mask is true
        scalar_loc = scalar_loc + inputfield( vi) * region%mesh%A( vi)
      else
        ! Do nothing
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

  subroutine write_to_file_scalar_FL( region, scalar)
    !< Write to FLUX scalar to single ISMIP regional output NetCDF file

    ! In/output variables:
    type(type_model_region),  intent(inout) :: region
    type(type_ismip_scalar),  intent(inout) :: scalar

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
    scalar_loc = scalar%accum / deltat
    scalar%accum = 0._dp

    ! Write scalar to file
    call write_to_field_multopt_dp_0D( scalar%filename, ncid, scalar%name, scalar_loc)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_file_scalar_FL

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
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%litempbotgr)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%litempbotfl)

    ! Basal drag
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%strbasemag)

    ! Lateral mass balance
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%licalvf)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%ligroundf)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%lifmassbf)

    ! Area fractions
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%sftgif)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%sftgrf)
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%sftflf)

    ! Other stuff
    call create_single_ISMIP_regional_output_file( region%ismip_output, region%ismip_output%tfbase)

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

    ! Define the name of the subfolder
    region%ismip_output%folder = trim( C%output_dir) // trim(C%ismip_counter) // '/'

    ! Basic topography
    call initialise_ISMIP_field( region, region%ismip_output%lithk, 'lithk', 'Ice thickness', 'land_ice_thickness', 'm', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%orog,  'orog' , 'Surface elevation', 'surface_altitude', 'm', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%topg,  'topg' , 'Bedrock elevation', 'bedrock_altitude', 'm', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%base,  'base' , 'Ice base elevation', 'base_altitude', 'm', 'ST')

    ! Geothermal heat flux
    call initialise_ISMIP_field( region, region%ismip_output%hfgeoubed, 'hfgeoubed' , &
      'Geothermal heat flux', 'upward_geothermal_heat_flux_in_land_ice', 'W m-2', 'FL')

    ! Surface and basal mass balances
    call initialise_ISMIP_field( region, region%ismip_output%acabf, 'acabf' , &
      'Surface mass balance flux', 'land_ice_surface_specific_mass_balance_flux', 'kg m-2 s-1', 'FL')

    call initialise_ISMIP_field( region, region%ismip_output%libmassbfgr, 'libmassbfgr' , &
      'Basal mass balance flux beneath grounded ice', 'land_ice_basal_specific_mass_balance_flux', 'kg m-2 s-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_output%libmassbffl, 'libmassbffl' , &
      'Basal mass balance flux beneath floating ice', 'land_ice_basal_specific_mass_balance_flux', 'kg m-2 s-1', 'FL')

    ! Thickness tendency
    call initialise_ISMIP_field( region, region%ismip_output%dlithkdt, 'dlithkdt' , &
      'Ice thickness imbalance', 'tendency_of_land_ice_thickness', 'm s-1', 'FL')
    ! Store snapshot ice thickness in accum
    region%ismip_output%dlithkdt%accum( region%mesh%vi1:region%mesh%vi2) = region%ice%geom%Hi( region%mesh%vi1:region%mesh%vi2)

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
    call initialise_ISMIP_field( region, region%ismip_output%ligroundf, 'ligroundf' , &
      'Calving flux', 'land_ice_specific_grounding_line_flux', 'kg m-2 s-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_output%lifmassbf,   'lifmassbf' , &
      'Ice front melt flux', 'land_ice_specific_mass_flux_due_to_ice_front_melting', 'kg m-2 s-1', 'FL')

    ! Area fractions
    call initialise_ISMIP_field( region, region%ismip_output%sftgif, 'sftgif' , &
      'Land ice area fraction', 'land_ice_area_fraction', '1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%sftgrf, 'sftgrf' , &
      'Grounded ice sheet area fraction', 'grounded_ice_sheet_area_fraction', '1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_output%sftflf, 'sftflf' , &
      'Floating ice sheet area fraction', 'floating_ice_shelf_area_fraction', '1', 'ST')

    ! Other stuff
    call initialise_ISMIP_field( region, region%ismip_output%tfbase, 'tfbase' , &
      'Thermal forcing under floating ice shelves', '', 'K', 'ST')

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
      'Total SMB flux', 'tendency_of_land_ice_mass_due_to_surface_mass_balance', 'kg s-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_output%tendlibmassbfgr, 'tendlibmassbfgr', &
      'Total BMB flux beneath grounded ice', 'tendency_of_land_ice_mass_due_to_basal_mass_balance', 'kg s-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_output%tendlibmassbffl, 'tendlibmassbffl', &
      'Total BMB flux beneath floating ice', 'tendency_of_land_ice_mass_due_to_basal_mass_balance', 'kg s-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_output%tendlicalvf, 'tendlicalvf', &
      'Total calving flux', 'tendency_of_land_ice_mass_due_to_calving', 'kg s-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_output%tendlifmassbf, 'tendlifmassbf', &
      'Total ice front melting flux', 'tendency_of_land_ice_mass_due_to_ice_front_melting', 'kg s-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_output%tendligroundf, 'tendligroundf', &
      'Total grounding line flux', 'tbd', 'kg s-1', 'FL')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ISMIP_output

  subroutine initialise_ISMIP_field_grid( region, field, name, long_name, standard_name, units, fieldtype)
    ! Initialise a single field

    type(type_model_region),        intent(inout) :: region
    type(type_ismip_gridded_field), intent(inout) :: field
    character(len=*),               intent(in   ) :: name
    character(len=*),               intent(in   ) :: long_name
    character(len=*),               intent(in   ) :: standard_name
    character(len=*),               intent(in   ) :: units
    character(len=2),               intent(in   ) :: fieldtype

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

    ! Convert grid resolution and start/end times to string.
    ! Offset end_year of 1 year to indicate last completed year
    write(start_year, '(I4)') int(C%start_time_of_run)
    write(end_year, '(I4)') int(C%end_time_of_run - 1)

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
    ! Offset end_year of 1 year to indicate last completed year
    write(start_year, '(I4)') int(C%start_time_of_run)
    write(end_year, '(I4)') int(C%end_time_of_run - 1)

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
    call reallocate_bounds( ismip_output%ligroundf%accum, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( ismip_output%lifmassbf%accum, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( ismip_output%dlithkdt%accum, mesh_new%vi1, mesh_new%vi2)

    ! Use accum to store current Hi for dHidt
    ismip_output%dlithkdt%accum = ice%geom%Hi

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_ISMIP_output

end module ismip_output_files
