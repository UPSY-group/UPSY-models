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
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use netcdf_io_main
  use ice_mass_and_fluxes, only: calc_ice_margin_fluxes
  use remapping_main, only: map_from_mesh_vertices_to_xy_grid_2D, &
    map_from_mesh_triangles_to_xy_grid_2D
  use mpi_distributed_memory, only: gather_to_all
  use ismip_output_types, only: type_ismip_grid_output, type_ismip_gridded_field
  use mesh_zeta, only: vertical_average
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_local
  use netcdf, only: NF90_GLOBAL
  use reallocate_mod, only: reallocate_bounds

  implicit none

  private

  public :: create_ISMIP_regional_output_files_grid, write_to_ISMIP_regional_output_files_grid, &
    accumulate_ISMIP_flux_fields, remap_ISMIP_grid_output

contains

  subroutine accumulate_ISMIP_flux_fields( region)
    !< Write to ISMIP regional output NetCDF files - grid version

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=1024), parameter                       :: routine_name = 'accumulate_ISMIP_flux_fields'
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: calving_flux
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
    call calc_ice_margin_fluxes( region%mesh, region%ice, calving_flux)

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
    deltat = region%time - region%ismip_grid_output%t_curr

    ! Accumulate FL fields. Exceptions:
    ! - dlithkdt (accumulation is used to store the previous written thickness
    ! - hfgeoubed (doesn't vary, so just writing out the initial snapshot)
    ! - lifmassbf (not computed yet, so just spitting out zeros)
    call accumulate_single_ISMIP_flux_field( region, &
      region%ismip_grid_output%acabf, SMB_loc, mask_ice, deltat)
    call accumulate_single_ISMIP_flux_field( region, &
      region%ismip_grid_output%libmassbfgr, region%BMB%BMB, region%ice%mask_grounded_ice, deltat)
    call accumulate_single_ISMIP_flux_field( region, &
      region%ismip_grid_output%libmassbffl, region%BMB%BMB, region%ice%mask_floating_ice,  deltat)
    call accumulate_single_ISMIP_flux_field( region, &
      region%ismip_grid_output%licalvf, calving_flux, mask_ice, deltat)

    ! Update current time
    region%ismip_grid_output%t_curr = region%time

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine accumulate_ISMIP_flux_fields

  subroutine accumulate_single_ISMIP_flux_field( region, field, d_partial, mask, deltat)
    !< Write to ISMIP regional output NetCDF files - grid version

    ! In/output variables:
    type(type_model_region),                              intent(inout) :: region
    type(type_ismip_gridded_field),                       intent(inout) :: field
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2), intent(in   ) :: d_partial
    logical,  dimension(region%mesh%vi1:region%mesh%vi2), intent(in   ) :: mask
    real(dp),                                             intent(in   ) :: deltat

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'accumulate_single_ISMIP_flux_field'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = region%mesh%vi1, region%mesh%vi2
      if (mask( vi)) then
        field%accum( vi) = field%accum( vi) + d_partial( vi) * ice_density / sec_per_year * deltat
      end if
    end do

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
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%libmassbfgr)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%libmassbffl)

    ! Thickness tendency
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%dlithkdt)

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
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%litempavg)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%litempgradgr)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%litempgradfl)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%litempbotgr)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%litempbotfl)

    ! Basal drag
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%strbasemag)

    ! Lateral mass balance
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%licalvf)
    call write_to_single_ISMIP_regional_output_file_grid( region, region%ismip_grid_output%lifmassbf)

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
    real(dp), dimension(:),   allocatable :: dTdzeta
    real(dp), dimension(region%mesh%vi1:region%mesh%vi2) :: SMB_loc, calving_flux
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
        ! Prevent negative values due to rounding errors
        d_grid_vec_partial_2D = max(0._dp, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
      case ('orog')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%Hs, d_grid_vec_partial_2D)
        ! Prevent negative values due to rounding errors
        d_grid_vec_partial_2D = max(0._dp, d_grid_vec_partial_2D)
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

      case ('libmassbfgr')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        if (field%is_initial) then
          ! First timeframe, no accumulation yet. Following protocol, using snapshot field instead
          do vi = region%mesh%vi1, region%mesh%vi2
            if (region%ice%mask_grounded_ice( vi)) then
              d_mesh_vec_partial_2D( vi) = region%BMB%BMB( vi) * ice_density / sec_per_year
            else
              d_mesh_vec_partial_2D( vi) = NaN
            end if
          end do 
          field%is_initial = .false.
        else
          d_mesh_vec_partial_2D = field%accum / deltat
          field%accum( region%mesh%vi1: region%mesh%vi2) = 0._dp
        end if
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
        deallocate( d_mesh_vec_partial_2D)

      case ('libmassbffl')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        if (field%is_initial) then
          ! First timeframe, no accumulation yet. Following protocol, using snapshot field instead
          do vi = region%mesh%vi1, region%mesh%vi2
            if (region%ice%mask_floating_ice( vi)) then
              d_mesh_vec_partial_2D( vi) = region%BMB%BMB( vi) * ice_density / sec_per_year
            else
              d_mesh_vec_partial_2D( vi) = NaN
            end if
          end do 
          field%is_initial = .false.
        else
          d_mesh_vec_partial_2D = field%accum / deltat
          field%accum( region%mesh%vi1: region%mesh%vi2) = 0._dp
        end if
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
          d_mesh_vec_partial_2D = (region%ice%Hi - field%accum) / (region%time - region%ismip_grid_output%t_prev) / sec_per_year
          field%accum = region%ice%Hi
        end if
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
        deallocate( d_mesh_vec_partial_2D)

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

      case ('litempavg')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        do vi = region%mesh%vi1, region%mesh%vi2
          if (region%ice%Hi( vi) > 0._dp) then
            d_mesh_vec_partial_2D( vi) = vertical_average( region%mesh%zeta, region%ice%Ti( vi, :))
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
            ! TODO replace matrix with first-order once available
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

      case ('litempgradfl')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        allocate( dTdzeta( 1:region%mesh%nz))
        do vi = region%mesh%vi1, region%mesh%vi2
          if (region%ice%mask_floating_ice( vi)) then
            ! TODO replace matrix with first-order once available
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

      ! Basal drag
      case ('strbasemag')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%basal_shear_stress, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)

      ! Lateral mass balance
      case ('licalvf')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        if (field%is_initial) then
          ! First timeframe, no accumulation yet. Following protocol, using snapshot field instead
          call calc_ice_margin_fluxes( region%mesh, region%ice, calving_flux)
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

      ! Area fractions
      case ('sftgif')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%fraction_margin, d_grid_vec_partial_2D)
        ! Prevent values outside [0,1] due to rounding errors
        d_grid_vec_partial_2D = min(1._dp,max(0._dp, d_grid_vec_partial_2D))
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
      case ('sftgrf')
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, region%ice%fraction_gr, d_grid_vec_partial_2D)
        ! Prevent values outside [0,1] due to rounding errors
        d_grid_vec_partial_2D = min(1._dp,max(0._dp, d_grid_vec_partial_2D))
        call write_to_field_multopt_grid_dp_2D( region%output_grid, field%filename, ncid, field%name, d_grid_vec_partial_2D)
      case ('sftflf')
        allocate( d_mesh_vec_partial_2D( region%mesh%vi1:region%mesh%vi2))
        ! TODO check whether this is appropriate
        d_mesh_vec_partial_2D = region%ice%fraction_margin( region%mesh%vi1: region%mesh%vi2) - region%ice%fraction_gr( region%mesh%vi1: region%mesh%vi2)
        call map_from_mesh_vertices_to_xy_grid_2D( region%mesh, region%output_grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        ! Prevent values outside [0,1] due to rounding errors
        d_grid_vec_partial_2D = min(1._dp,max(0._dp, d_grid_vec_partial_2D))
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
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%libmassbfgr)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%libmassbffl)

    ! Thickness tendency
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%dlithkdt)

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
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%litempavg)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%litempgradgr)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%litempgradfl)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%litempbotgr)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%litempbotfl)

    ! Basal drag
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%strbasemag)

    ! Lateral mass balance
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%licalvf)
    call create_single_ISMIP_regional_output_file_grid( region%ismip_grid_output, region%ismip_grid_output%lifmassbf)

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
    call setup_xy_grid_in_netcdf_file( field%filename, ncid, ismip_grid_output%grid)

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
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'grid_type', trim(ismip_grid_output%IS_name))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'grid_resolution', trim(res_str) // 'm')
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'group', trim(C%ismip_group_name))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'model', trim(C%ismip_model_name))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'scenario', trim(C%ismip_scenario_name))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'contact_name', trim(C%ismip_contact_name))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'contact_email', trim(C%ismip_contact_email))
    call add_attribute_char( field%filename, ncid, NF90_GLOBAL, 'crs', trim(ismip_grid_output%crs))

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
        region%ismip_grid_output%crs = 'EPSG:3031'
      case ('GRL')
        region%ismip_grid_output%IS_name = 'GIS'
        region%ismip_grid_output%crs = 'EPSG:3413'
    end select

    ! Copy the output grid
    region%ismip_grid_output%grid = region%output_grid

    ! Initialise the time trackers
    region%ismip_grid_output%t_prev = C%start_time_of_run
    region%ismip_grid_output%t_curr = C%start_time_of_run

    ! Basic topography
    call initialise_ISMIP_field( region, region%ismip_grid_output%lithk, 'lithk', 'Ice thickness', 'land_ice_thickness', 'm', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%orog,  'orog' , 'Surface elevation', 'surface_altitude', 'm', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%topg,  'topg' , 'Bedrock elevation', 'bedrock_altitude', 'm', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%base,  'base' , 'Ice base elevation', 'base_altitude', 'm', 'ST')

    ! Geothermal heat flux
    call initialise_ISMIP_field( region, region%ismip_grid_output%hfgeoubed, 'hfgeoubed' , &
      'Geothermal heat flux', 'upward_geothermal_heat_flux_in_land_ice', 'W m-2', 'FL')

    ! Surface and basal mass balances
    call initialise_ISMIP_field( region, region%ismip_grid_output%acabf, 'acabf' , &
      'Surface mass balance flux', 'land_ice_surface_specific_mass_balance_flux', 'kg m-2 s-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_grid_output%libmassbfgr, 'libmassbfgr' , &
      'Basal mass balance flux beneath grounded ice', 'land_ice_basal_specific_mass_balance_flux', 'kg m-2 s-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_grid_output%libmassbffl, 'libmassbffl' , &
      'Basal mass balance flux beneath floating ice', 'land_ice_basal_specific_mass_balance_flux', 'kg m-2 s-1', 'FL')

    ! Thickness tendency
    call initialise_ISMIP_field( region, region%ismip_grid_output%dlithkdt, 'dlithkdt' , &
      'Ice thickness imbalance', 'tendency_of_land_ice_thickness', 'm s-1', 'FL')
    ! Use accum of thickness tendency to store previous ice thickness
    region%ismip_grid_output%dlithkdt%accum = region%ice%Hi

    ! Velocities
    call initialise_ISMIP_field( region, region%ismip_grid_output%xvelsurf, 'xvelsurf' , &
      'Surface velocity in x', 'land_ice_surface_x_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%yvelsurf, 'yvelsurf' , &
      'Surface velocity in y', 'land_ice_surface_y_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%zvelsurf, 'zvelsurf' , &
      'Surface velocity in z', 'land_ice_surface_upward_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%xvelbase, 'xvelbase' , &
      'Basal velocity in x', 'land_ice_basal_x_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%yvelbase, 'yvelbase' , &
      'Basal velocity in y', 'land_ice_basal_y_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%zvelbase, 'zvelbase' , &
      'Basal velocity in z', 'land_ice_basal_upward_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%xvelmean, 'xvelmean' , &
      'Mean velocity in x', 'land_ice_vertical_mean_x_velocity', 'm s-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%yvelmean, 'yvelmean' , &
      'Mean velocity in y', 'land_ice_vertical_mean_y_velocity', 'm s-1', 'ST')

    ! Temperatures
    call initialise_ISMIP_field( region, region%ismip_grid_output%litemptop, 'litemptop' , &
      'Surface temperature', 'temperature_at_top_of_ice_sheet_model', 'K', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%litempavg, 'litempavg' , &
      'Depth average temperature', 'land_ice_temperature', 'K', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%litempgradgr, 'litempgradgr' , &
      'Vertical Basal temperature gradient beneath grounded ice sheet', & 
      'temperature_gradient_at_base_of_ice_sheet_model', 'K m-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%litempgradfl, 'litempgradfl' , &
      'Vertical Basal temperature gradient beneath floating ice sheet', & 
      'temperature_gradient_at_base_of_ice_sheet_model', 'K m-1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%litempbotgr, 'litempbotgr' , &
      'Basal temperature beneath grounded ice', 'temperature_at_base_of_ice_sheet_model', 'K', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%litempbotfl, 'litempbotfl' , &
      'Basal temperature beneath floating ice', 'temperature_at_base_of_ice_sheet_model', 'K', 'ST')

    ! Basal drag
    call initialise_ISMIP_field( region, region%ismip_grid_output%strbasemag, 'strbasemag' , &
      'Basal drag', 'land_ice_basal_drag', 'Pa', 'ST')

    ! Lateral mass balance
    call initialise_ISMIP_field( region, region%ismip_grid_output%licalvf, 'licalvf' , &
      'Calving flux', 'land_ice_specific_mass_flux_due_to_calving', 'kg m-2 s-1', 'FL')
    call initialise_ISMIP_field( region, region%ismip_grid_output%lifmassbf,   'lifmassbf' , &
      'Ice front melt flux', 'land_ice_specific_mass_flux_due_to_ice_front_melting', 'kg m-2 s-1', 'FL')

    ! Area fractions
    call initialise_ISMIP_field( region, region%ismip_grid_output%sftgif, 'sftgif' , &
      'Land ice area fraction', 'land_ice_area_fraction', '1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%sftgrf, 'sftgrf' , &
      'Grounded ice sheet area fraction', 'grounded_ice_sheet_area_fraction', '1', 'ST')
    call initialise_ISMIP_field( region, region%ismip_grid_output%sftflf, 'sftflf' , &
      'Floating ice sheet area fraction', 'floating_ice_shelf_area_fraction', '1', 'ST')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ISMIP_grid_output

  subroutine initialise_ISMIP_field( region, field, name, long_name, standard_name, units, fieldtype)
    ! Initialise a single field

    type(type_model_region),           intent(inout) :: region
    type(type_ismip_gridded_field),    intent(inout) :: field
    character(len=*),                  intent(in   ) :: name
    character(len=*),                  intent(in   ) :: long_name
    character(len=*),                  intent(in   ) :: standard_name
    character(len=*),                  intent(in   ) :: units
    character(len=2),                  intent(in   ) :: fieldtype

    ! Local variables
    character(len=1024), parameter :: routine_name = 'initialise_ISMIP_field'
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
    region%ismip_grid_output%folder = trim( C%output_dir) // 'CORE' // '/'

    ! Define the filename for this field
    field%filename = trim( region%ismip_grid_output%folder) // trim(field%name) // '_' // &
      trim(region%ismip_grid_output%IS_name) // '_' // trim(C%ismip_group_name) // '_' // &
      trim(C%ismip_model_name) // '_' // trim(C%ismip_member_id) // '_' // &
      trim(C%ismip_esm_name) // '_' // trim(C%ismip_forcing_member_id) // '_' // &
      trim(C%ismip_scenario_name) // '_' // trim(C%ismip_counter) // '_' // &
      trim(start_year) // '-' // trim(end_year) // '.nc'

    ! Allocate fields for accumulation during each timestep for flux fields
    if ((field%fieldtype == 'FL') .and. (field%name /= 'hfgeoubed')) then
      allocate(field%accum( region%mesh%vi1: region%mesh%vi2), source=0._dp)
      field%is_initial = .true.
    end if

    ! Initialise the counter for number of time values
    field%nt = 0

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ISMIP_field

  subroutine remap_ISMIP_grid_output( mesh_old, mesh_new, ice, ismip_grid_output)
    ! Reallocate the accumulated fields and redefine Hi_prev

    type(type_mesh),                   intent(in   ) :: mesh_old
    type(type_mesh),                   intent(in   ) :: mesh_new
    type(type_ice_model),              intent(in   ) :: ice
    type(type_ismip_grid_output),      intent(inout) :: ismip_grid_output

    ! Local variables
    character(len=1024), parameter :: routine_name = 'remap_ISMIP_grid_output'

    ! Add routine to path
    call init_routine( routine_name)

    call reallocate_bounds( ismip_grid_output%acabf%accum, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( ismip_grid_output%libmassbfgr%accum, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( ismip_grid_output%libmassbffl%accum, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( ismip_grid_output%licalvf%accum, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( ismip_grid_output%lifmassbf%accum, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( ismip_grid_output%dlithkdt%accum, mesh_new%vi1, mesh_new%vi2)

    ! Use accum to store current Hi for dHidt
    ismip_grid_output%dlithkdt%accum = ice%Hi

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_ISMIP_grid_output

end module ismip_grid_output_files
