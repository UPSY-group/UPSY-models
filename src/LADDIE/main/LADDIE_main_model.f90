module LADDIE_main_model

  ! The main stand-alone model of LADDIE

! ===== Preamble =====
! ====================

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, MPI_WTIME
  use precisions, only: dp
  use mpi_basic, only: par, sync 
  use control_resources_and_error_messaging, only: happy, warning, crash, init_routine, finalise_routine, &
    colour_string, str2int, int2str, insert_val_into_string_dp
  use model_configuration                                    , ONLY: C
  use parameters
  use reference_geometry_types, only: type_reference_geometry
  use laddie_model_types, only: type_laddie_model, type_laddie_timestep
  use laddie_forcing_types, only: type_laddie_forcing
  use laddie_main_utils, only: update_laddie_masks, extrapolate_laddie_variables, repartition_laddie
  use laddie_utilities, only: compute_ambient_TS, allocate_laddie_model, allocate_laddie_timestep, &
    map_H_a_b, map_H_a_c, allocate_laddie_forcing
  use laddie_operators, only: update_laddie_operators
  use laddie_mesh_output, only: create_laddie_mesh_output_file, write_to_laddie_mesh_output_file
  use laddie_grid_output, only: create_laddie_output_file_grid, write_to_laddie_output_file_grid
  use laddie_scalar_output, only: create_laddie_scalar_output_file, write_to_laddie_scalar_output_file, &
    buffer_laddie_scalars
  use laddie_integration, only: integrate_euler, integrate_fbrk3, integrate_lfra, move_laddie_timestep
  use laddie_physics, only: compute_subglacial_discharge, compute_SGD_at_transects
  use mesh_types, only: type_mesh
  use ocean_main, only: initialise_ocean_model, run_ocean_model
  use netcdf_io_main
  use mesh_creation_main, only: create_mesh_from_gridded_geometry, create_mesh_from_meshed_geometry, write_mesh_success
  use grid_basic, only: setup_square_grid
  use mesh_output_files, only: create_main_regional_output_file_mesh, write_to_main_regional_output_file_mesh
  use grid_output_files, only: create_main_regional_output_file_grid, write_to_main_regional_output_file_grid, &
    create_main_regional_output_file_grid_ROI, write_to_main_regional_output_file_grid_ROI
  use scalar_output_files, only: create_scalar_regional_output_file, buffer_scalar_output, write_to_scalar_regional_output_file
  use mesh_ROI_polygons
  use apply_maps, only: clear_all_maps_involving_this_mesh
  use mesh_memory, only: deallocate_mesh
  use mesh_halo_exchange, only: exchange_halos
  use mesh_repartitioning, only: repartition_mesh
  use checksum_mod, only: checksum

  implicit none

contains

! ===== Main routine =====
! ========================

  subroutine run_laddie_model( mesh, laddie, forcing, time, is_initial, is_standalone)
    ! Integrate the model until t_end

    ! In/output variables
    type(type_mesh),           intent(in   ) :: mesh 
    type(type_laddie_model),   intent(inout) :: laddie
    type(type_laddie_forcing), intent(inout) :: forcing
    real(dp),                  intent(in   ) :: time 
    logical,                   intent(in   ) :: is_initial
    logical,                   intent(in   ) :: is_standalone

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_laddie_model'
    type(type_mesh)                :: mesh_repartitioned

    ! Add routine to path
    call init_routine( routine_name)

    if (C%do_repartition_laddie) then
      ! Repartition the mesh so each process has (approximately)
      ! the same number of ice shelf vertices/triangles
      call repartition_mesh( mesh, mesh_repartitioned, laddie%mask_a, laddie%mask_b)

      ! Repartition laddie
      call repartition_laddie( mesh, mesh_repartitioned, laddie, forcing)

      ! Run laddie on the repartitioned mesh
      call run_laddie_model_leg( mesh_repartitioned, laddie, forcing, time, is_initial, is_standalone)

      ! Un-repartition laddie
      call repartition_laddie( mesh_repartitioned, mesh, laddie, forcing)
    else
      ! Run laddie on the original mesh
      call run_laddie_model_leg( mesh, laddie, forcing, time, is_initial, is_standalone)
    end if

    ! Clean up after yourself
    call deallocate_mesh( mesh_repartitioned)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_laddie_model

  subroutine run_laddie_model_leg( mesh, laddie, forcing, time, is_initial, is_standalone)
    ! Run one leg of the laddie model

    ! In/output variables
    type(type_mesh),           intent(in   ) :: mesh
    type(type_laddie_model),   intent(inout) :: laddie
    type(type_laddie_forcing), intent(in   ) :: forcing
    real(dp),                  intent(in   ) :: time
    logical,                   intent(in   ) :: is_initial
    logical,                   intent(in   ) :: is_standalone

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_laddie_model_leg'
    integer                        :: vi, ti, iyear
    real(dp)                       :: tl               ! [s] Laddie time
    real(dp)                       :: dt               ! [s] Laddie time step
    real(dp)                       :: duration         ! [days] Duration of run
    real(dp)                       :: ref_time         ! [s] Reference time for writing
    real(dp), parameter            :: time_relax_laddie = 0.02_dp ! [days]
    real(dp), parameter            :: fac_dt_relax = 3.0_dp ! Reduction factor of time step
    real(dp)                       :: time_to_write    ! [days]
    real(dp)                       :: last_write_time  ! [days]
    real(dp)                       :: time_to_write_mesh    ! [days]
    real(dp)                       :: last_write_time_mesh  ! [days]
    real(dp)                       :: time_to_write_grid    ! [days]
    real(dp)                       :: last_write_time_grid  ! [days]
    real(dp)                       :: fac_time_varying_SGD

    ! Add routine to path
    call init_routine( routine_name)

    ! == Preparation ==
    ! =================

    ! Extrapolate data into new cells
    call extrapolate_laddie_variables( mesh, laddie, forcing)

    ! == Update masks ==
    call update_laddie_masks( mesh, laddie, forcing)

    ! == Compute SGD ==
    select case (C%choice_laddie_SGD)
      case default
        call crash('unknown choice_laddie_SGD "' // trim( C%choice_laddie_SGD) // '"!')
      case ('none')
        ! Do nothing
      case ('idealised','read_from_file')
        if (time >= C%start_time_of_applying_SGD) then 
          ! Compute SGD
          call compute_subglacial_discharge( mesh, laddie, forcing)
          ! Check if apply boost
          if (C%apply_SGD_boost_every_10_years) then
              ! Shift by +1 year so the window moves to 9–10, 19–20, 29–30, ...
              iyear = floor(time + 1.0_dp)
              if (mod(iyear, 10) == 0) then
                  laddie%SGD = laddie%SGD * C%SGD_boost_factor
              end if
          end if
        else
          ! Set SGD to zero
          laddie%SGD = 0._dp
        end if
      case ('read_transects')
        ! Compute SGD from transects
        call compute_SGD_at_transects(mesh, laddie, forcing)
        ! Check if apply boost
        if (C%apply_SGD_boost_every_10_years) then
            ! Shift by +1 year so the window moves to 9–10, 19–20, 29–30, ...
            iyear = floor(time + 1.0_dp)
            if (mod(iyear, 10) == 0) then
                laddie%SGD = laddie%SGD * C%SGD_boost_factor
            end if
        end if
    end select

    ! Set values to zero if outside laddie mask
    do vi = mesh%vi1, mesh%vi2
      if (.not. laddie%mask_a( vi)) then
        laddie%now%H( vi)     = 0.0_dp
        laddie%now%T( vi)     = 0.0_dp
        laddie%now%S( vi)     = 0.0_dp
        laddie%melt( vi)  = 0.0_dp
        laddie%entr( vi)  = 0.0_dp
      end if
    end do

    call checksum( laddie%now%H, 'laddie%now%H', mesh%pai_V)
    call checksum( laddie%now%T, 'laddie%now%T', mesh%pai_V)
    call checksum( laddie%now%S, 'laddie%now%S', mesh%pai_V)
    call checksum( laddie%melt , 'laddie%melt' , mesh%pai_V)
    call checksum( laddie%entr , 'laddie%entr' , mesh%pai_V)

    do ti = mesh%ti1, mesh%ti2
      if (.not. laddie%mask_b( ti)) then
        laddie%now%U( ti)     = 0.0_dp
        laddie%now%V( ti)     = 0.0_dp
        laddie%now%H_b( ti)   = 0.0_dp
      end if
    end do

    call checksum( laddie%now%U  , 'laddie%now%U'  , mesh%pai_Tri)
    call checksum( laddie%now%V  , 'laddie%now%V'  , mesh%pai_Tri)
    call checksum( laddie%now%H_b, 'laddie%now%H_b', mesh%pai_Tri)

    ! Simply set H_c zero everywhere, will be recomputed through mapping later
    laddie%now%H_c( mesh%ei1:mesh%ei2) = 0.0_dp

    ! == Update operators ==
    call update_laddie_operators( mesh, laddie)

    ! == Main time loop ==
    ! ====================

    ! Determine run duration and apply offset for initial run
    if (is_initial) then
      duration = C%time_duration_laddie_init
      ref_time = time*sec_per_year - duration*sec_per_day
    else
      duration = C%time_duration_laddie
      ref_time = time*sec_per_year
    end if

    tl = 0.0_dp
    last_write_time = 0.0_dp
    last_write_time_mesh = 0.0_dp
    last_write_time_grid = 0.0_dp
    time_to_write = C%time_interval_scalar_output
    time_to_write_mesh = C%dt_output
    time_to_write_grid = C%dt_output_grid

    ! Perform first integration with half the time step for LFRA scheme
    dt = C%dt_laddie / fac_dt_relax
    if (C%choice_laddie_integration_scheme == 'lfra') then
      call integrate_lfra( mesh, laddie, forcing, tl, time, dt)
    end if

    do while (tl < duration * sec_per_day)

      ! Set time step
      if (tl < time_relax_laddie * sec_per_day) then
        ! Relaxation, take short time step
        dt = C%dt_laddie / fac_dt_relax
      else
        ! Regular timestep
        dt = C%dt_laddie
      end if

      select case(C%choice_laddie_integration_scheme)
        case default
          call crash('unknown choice_laddie_integration_scheme "' // trim( C%choice_laddie_integration_scheme) // '"')
        case ('euler')
          call integrate_euler( mesh, laddie, forcing, tl, time, dt)
        case ('fbrk3')
          call integrate_fbrk3( mesh, laddie, forcing, tl, time, dt)
        case ('lfra')
          call integrate_lfra( mesh, laddie, forcing, tl, time, 2*dt)
      end select

      ! Write to output
      if (is_standalone) then
        ! Always include scalar output
        call buffer_laddie_scalars( mesh, laddie, ref_time + tl)
  
        ! Write scalars if required
        if (tl > time_to_write * sec_per_day) then
          call write_to_laddie_scalar_output_file( laddie)
          last_write_time = time_to_write
          time_to_write = time_to_write + C%time_interval_scalar_output
        end if

        ! Write mesh if required
        if (tl > time_to_write_mesh * sec_per_day) then
          call write_to_laddie_mesh_output_file( mesh, laddie, forcing, ref_time + tl, is_standalone)
          last_write_time_mesh = time_to_write_mesh
          time_to_write_mesh = time_to_write_mesh + C%dt_output
        end if

        ! Write grid if required
        if (tl > time_to_write_grid * sec_per_day) then
          call write_to_laddie_output_file_grid( mesh, laddie, forcing, ref_time + tl)
          last_write_time_grid = time_to_write_grid
          time_to_write_grid = time_to_write_grid + C%dt_output_grid
        end if

      else
        if (C%do_write_laddie_output_fields) then
          ! Write if required
          if (tl > time_to_write_mesh * sec_per_day) then
            call write_to_laddie_mesh_output_file( mesh, laddie, forcing, ref_time + tl, is_standalone)
            last_write_time_mesh = time_to_write_mesh
            time_to_write_mesh = time_to_write_mesh + C%dt_output
          end if
        end if
  
        if (C%do_write_laddie_output_scalar) then
          call buffer_laddie_scalars( mesh, laddie, ref_time + tl)
  
          ! Write if required
          if (tl > time_to_write * sec_per_day) then
            call write_to_laddie_scalar_output_file( laddie)
            last_write_time = time_to_write
            time_to_write = time_to_write + C%time_interval_scalar_output
          end if
        end if
      end if

    end do !do while (tl < C%time_duration_laddie)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_laddie_model_leg

! ===== Model initialisation =====
! ================================

  subroutine initialise_laddie_model( mesh, laddie, forcing, is_standalone)
    ! Initialise the model

    implicit none

    ! In/output variables:
    type(type_mesh)                                    , intent(in)    :: mesh
    type(type_laddie_model)                            , intent(inout) :: laddie
    type(type_laddie_forcing)                          , intent(in)    :: forcing
    logical                                            , intent(in)    :: is_standalone

    ! Local variables:
    character(len=256), parameter                                      :: routine_name = 'initialise_laddie_model'
    character(len=256), parameter                                      :: grid_name = 'square_grid'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) write (0,'(A)') '   Initialising laddie model... '

    ! Allocate variables
    call allocate_laddie_model( mesh, laddie)

    ! == Update masks ==
    call update_laddie_masks( mesh, laddie, forcing)

    ! == Update operators ==
    call update_laddie_operators( mesh, laddie)

    ! Initialise requested timesteps
    call initialise_laddie_model_timestep( mesh, laddie, forcing, laddie%now)

    select case(C%choice_laddie_integration_scheme)
      case default
        call crash('unknown choice_laddie_integration_scheme "' // trim( C%choice_laddie_integration_scheme) // '"')
      case ('euler')
        call initialise_laddie_model_timestep( mesh, laddie, forcing, laddie%np1)
      case ('fbrk3')
        call initialise_laddie_model_timestep( mesh, laddie, forcing, laddie%np13)
        call initialise_laddie_model_timestep( mesh, laddie, forcing, laddie%np12)
        call initialise_laddie_model_timestep( mesh, laddie, forcing, laddie%np1)
      case ('lfra')
        call crash('LeapFrog RobertAsselin scheme does not work yet, use euler or fbrk3')
        call initialise_laddie_model_timestep( mesh, laddie, forcing, laddie%nm1)
        call initialise_laddie_model_timestep( mesh, laddie, forcing, laddie%np1)
    end select

    ! Create output file
    if (is_standalone) then
      ! Always include scalar output
      call create_laddie_scalar_output_file( laddie)
      ! Always include mesh output
      call create_laddie_mesh_output_file( mesh, laddie, is_standalone)
      ! Always include grid output
      call setup_square_grid( grid_name, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, C%dx_output_grid_ANT, laddie%output_grid, &
        lambda_M = mesh%lambda_M, phi_M = mesh%phi_M, beta_stereo = mesh%beta_stereo)
      call create_laddie_output_file_grid( mesh, laddie, forcing)
    else
      if (C%do_write_laddie_output_fields) call create_laddie_mesh_output_file( mesh, laddie, is_standalone)
      if (C%do_write_laddie_output_scalar) call create_laddie_scalar_output_file( laddie)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_laddie_model

  subroutine initialise_laddie_model_timestep( mesh, laddie, forcing, npx)
    ! Initialise the laddie model for given timestep

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(inout) :: laddie
    type(type_laddie_forcing),              intent(in   ) :: forcing
    type(type_laddie_timestep),             intent(inout) :: npx

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'initialise_laddie_model_timestep'
    integer                                               :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate timestep
    call allocate_laddie_timestep( mesh, npx)

    ! Layer thickness
    do vi = mesh%vi1, mesh%vi2
       if (laddie%mask_a( vi)) then
         npx%H( vi)      = C%laddie_initial_thickness
       end if
    end do
    call exchange_halos( mesh, npx%H)
    call checksum( npx%H, 'npx%H', mesh%pai_V)

    ! Layer thickness on b and c grid
    call map_H_a_b( mesh, laddie, npx%H, npx%H_b)
    call map_H_a_c( mesh, laddie, npx%H, npx%H_c)
    call checksum( npx%H_b, 'npx%H_b', mesh%pai_Tri)
    call checksum( npx%H_c, 'npx%H_c', mesh%pai_E)

    ! Initialise ambient T and S
    call compute_ambient_TS( mesh, laddie, forcing, npx%H)

    ! Initialise main T and S
    do vi = mesh%vi1, mesh%vi2
       if (laddie%mask_a( vi)) then
         npx%T( vi)      = laddie%T_amb( vi) + C%laddie_initial_T_offset
         npx%S( vi)      = laddie%S_amb( vi) + C%laddie_initial_S_offset
       end if
    end do
    call checksum( npx%T, 'npx%T', mesh%pai_V)
    call checksum( npx%S, 'npx%S', mesh%pai_V)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_laddie_model_timestep

end module LADDIE_main_model
