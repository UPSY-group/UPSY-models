MODULE ocean_main

  ! The main ocean model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  use UPSY_main, only: UPSY
  USE mpi_basic                                              , ONLY: par, sync
  USE call_stack_and_comp_time_tracking                  , ONLY: crash, init_routine, finalise_routine
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  use grid_types, only: type_grid
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ocean_model_types                                      , ONLY: type_ocean_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE ocean_utilities                                        , ONLY: initialise_ocean_vertical_grid, calc_ocean_temperature_at_shelf_base, calc_ocean_freezing_point_at_shelf_base, interpolate_ocean_depth
  USE ocean_realistic                                        , ONLY: initialise_ocean_model_realistic, run_ocean_model_realistic, remap_ocean_model_realistic
  USE ocean_idealised                                        , ONLY: initialise_ocean_model_idealised, run_ocean_model_idealised
  use netcdf_io_main
  use ocean_snapshot_nudge2D, only: initialise_ocean_model_snapshot_nudge2D, run_ocean_model_snapshot_nudge2D
  use reference_geometry_types, only: type_reference_geometry
  use ocean_snapshot_plus_anomalies, only: initialise_ocean_model_snapshot_plus_anomalies, &
    run_ocean_model_snapshot_plus_anomalies

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_ocean_model( mesh, grid_smooth, ice, ocean, region_name, time)
    ! Calculate the ocean

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    type(type_grid),                        intent(in   ) :: grid_smooth
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_ocean_model'
    CHARACTER(LEN=256)                                    :: choice_ocean_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if we need to calculate a new ocean
    IF (C%do_asynchronous_ocean) THEN
      ! Asynchronous coupling: do not calculate a new ocean in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next ocean time step
      IF (time == ocean%t_next) THEN
        ! Go on to calculate a new ocean
        ocean%t_next = time + C%dt_ocean
      ELSEIF (time > ocean%t_next) THEN
        ! This should not be possible
        CALL crash('overshot the ocean time step')
      ELSE
        ! It is not yet time to calculate a new ocean
        CALL finalise_routine( routine_name)
        RETURN
      END IF

    ELSE ! IF (C%do_asynchronous_ocean) THEN
      ! Synchronous coupling: calculate a new ocean in every model loop
      ocean%t_next = time + C%dt_ocean
    END IF

    ! Determine which ocean model to run for this region
    select case (region_name)
    case default
      call crash('unknown region_name "' // region_name // '"')
    case ('NAM')
      choice_ocean_model = C%choice_ocean_model_NAM
    case ('EAS')
      choice_ocean_model = C%choice_ocean_model_EAS
    case ('GRL')
      choice_ocean_model = C%choice_ocean_model_GRL
    case ('ANT')
      choice_ocean_model = C%choice_ocean_model_ANT
    end select

    ! Run the chosen ocean model
    select case( choice_ocean_model)
    case default
      call crash('unknown choice_ocean_model "' // trim( choice_ocean_model) // '"')
    case( 'none')
      ! No need to do anything
    case( 'idealised')
      call run_ocean_model_idealised( mesh, ice, ocean)
    case( 'realistic')
      call limit_ocean_supercooling( mesh, ice, ocean)
      call run_ocean_model_realistic( mesh, ice, ocean, time)
    case( 'snapshot+nudge2D')
      call run_ocean_model_snapshot_nudge2D( mesh, grid_smooth, ice, ocean, time)
    case( 'snapshot_plus_anomalies')
      call run_ocean_model_snapshot_plus_anomalies( mesh, ocean, time)
    end select

    ! Compute secondary variables
    CALL calc_ocean_temperature_at_shelf_base(    mesh, ice, ocean)
    CALL calc_ocean_freezing_point_at_shelf_base( mesh, ice, ocean)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model

  SUBROUTINE initialise_ocean_model( mesh, ice, ocean, region_name, start_time_of_run, refgeo_PD, refgeo_init)
    ! Initialise the ocean model

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(OUT)   :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: start_time_of_run
    type(type_reference_geometry),          intent(in   ) :: refgeo_PD, refgeo_init

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_ocean_model'
    CHARACTER(LEN=256)                                    :: choice_ocean_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising ocean model...'

    ! Determine which ocean model to run for this region
    select case (region_name)
    case default
      call crash('unknown region_name "' // region_name // '"')
    case ('NAM')
      choice_ocean_model = C%choice_ocean_model_NAM
    case ('EAS')
      choice_ocean_model = C%choice_ocean_model_EAS
    case ('GRL')
      choice_ocean_model = C%choice_ocean_model_GRL
    case ('ANT')
      choice_ocean_model = C%choice_ocean_model_ANT
    end select

    ! Initialise vertical grid: C%z_ocean and C%nz_ocean
    CALL initialise_ocean_vertical_grid

    ! Allocate memory for main variables
    ALLOCATE( ocean%T( mesh%vi1:mesh%vi2,C%nz_ocean))
    ALLOCATE( ocean%S( mesh%vi1:mesh%vi2,C%nz_ocean))
    ocean%T = 0._dp
    ocean%S = 0._dp

    ! Allocate memory for secondary variables
    ALLOCATE( ocean%T_draft(          mesh%vi1:mesh%vi2))
    ALLOCATE( ocean%T_freezing_point( mesh%vi1:mesh%vi2))
    ocean%T_draft          = 0._dp
    ocean%T_freezing_point = 0._dp

    ! Set time of next calculation to start time
    ocean%t_next = C%start_time_of_run

    ! Determine which ocean model to initialise
    select case( choice_ocean_model)
    case default
      call crash('unknown choice_ocean_model "' // trim( choice_ocean_model) // '"')
    case( 'none')
    case( 'idealised')
      call initialise_ocean_model_idealised( mesh, ocean)
    case( 'realistic')
      call initialise_ocean_model_realistic( mesh, ice, ocean, region_name, start_time_of_run)
    case( 'snapshot+nudge2D')
      call initialise_ocean_model_snapshot_nudge2D( mesh, ocean%snapshot_nudge2D, region_name, refgeo_PD, refgeo_init)
    case( 'snapshot_plus_anomalies')
      call initialise_ocean_model_snapshot_plus_anomalies( mesh, ocean%snapshot_plus_anomalies)
    end select

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ocean_model

  SUBROUTINE write_to_restart_file_ocean_model( mesh, ocean, region_name, time)
    ! Write to the restart file for the ocean model

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(IN)    :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'write_to_restart_file_ocean_model'
    CHARACTER(LEN=256)                                    :: choice_ocean_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which ocean model to run for this region
    select case (region_name)
    case default
      call crash('unknown region_name "' // region_name // '"')
    case ('NAM')
      choice_ocean_model = C%choice_ocean_model_NAM
    case ('EAS')
      choice_ocean_model = C%choice_ocean_model_EAS
    case ('GRL')
      choice_ocean_model = C%choice_ocean_model_GRL
    case ('ANT')
      choice_ocean_model = C%choice_ocean_model_ANT
    end select

    ! Write to the restart file of the chosen ocean model
    select case(choice_ocean_model)
    case default
      call crash('unknown choice_ocean_model "' // trim( choice_ocean_model) // '"')
    case( 'none', 'idealised', 'snapshot+nudge2D', 'snapshot_plus_anomalies')
      ! No need to do anything
    case( 'realistic')
      call write_to_restart_file_ocean_model_region( mesh, ocean, region_name, time)
    end select

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_ocean_model

  SUBROUTINE write_to_restart_file_ocean_model_region( mesh, ocean, region_name, time)
    ! Write to the restart NetCDF file for the ocean model

    ! In/output variables:
    TYPE(type_mesh),        INTENT(IN)    :: mesh
    TYPE(type_ocean_model), INTENT(IN)    :: ocean
    CHARACTER(LEN=3),       INTENT(IN)    :: region_name
    REAL(dp),               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER         :: routine_name = 'write_to_restart_file_ocean_model_region'
    INTEGER                               :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%primary) WRITE(0,'(A)') '   Writing to ocean restart file "' // &
      UPSY%stru%colour_string( TRIM( ocean%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( ocean%restart_filename, ncid)

    ! Write the time to the file
    CALL write_time_to_file( ocean%restart_filename, ncid, time)

    ! ! Write the velocity fields to the file
    CALL write_to_field_multopt_mesh_dp_3D_ocean( mesh, ocean%restart_filename, ncid, 'T', ocean%T)
    CALL write_to_field_multopt_mesh_dp_3D_ocean( mesh, ocean%restart_filename, ncid, 'S', ocean%S)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_ocean_model_region

  SUBROUTINE create_restart_file_ocean_model( mesh, ocean, region_name)
    ! Create the restart file for the ocean model

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_file_ocean_model'
    CHARACTER(LEN=256)                                    :: choice_ocean_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which ocean model to run for this region
    select case (region_name)
    case default
      call crash('unknown region_name "' // region_name // '"')
    case ('NAM')
      choice_ocean_model = C%choice_ocean_model_NAM
    case ('EAS')
      choice_ocean_model = C%choice_ocean_model_EAS
    case ('GRL')
      choice_ocean_model = C%choice_ocean_model_GRL
    case ('ANT')
      choice_ocean_model = C%choice_ocean_model_ANT
    end select

    ! Create the restart file of the chosen ocean model
    select case(choice_ocean_model)
    case default
      call crash('unknown choice_ocean_model "' // trim( choice_ocean_model) // '"')
    case( 'none', 'idealised', 'snapshot+nudge2D', 'snapshot_plus_anomalies')
      ! No need to do anything
    case( 'realistic')
      call create_restart_file_ocean_model_region( mesh, ocean, region_name)
    end select

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_ocean_model

  SUBROUTINE create_restart_file_ocean_model_region( mesh, ocean, region_name)
    ! Create a restart NetCDF file for the ocean submodel
    ! Includes generation of the procedural filename (e.g. "restart_ocean_00001.nc")

    ! In/output variables:
    TYPE(type_mesh),        INTENT(IN)    :: mesh
    TYPE(type_ocean_model), INTENT(INOUT) :: ocean
    CHARACTER(LEN=3),       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER         :: routine_name = 'create_restart_file_ocean_model_region'
    CHARACTER(LEN=256)                    :: filename_base
    INTEGER                               :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'restart_ocean_' // region_name
    CALL generate_filename_XXXXXdotnc( filename_base, ocean%restart_filename)

    ! Print to terminal
    IF (par%primary) WRITE(0,'(A)') '   Creating ocean model restart file "' // &
      UPSY%stru%colour_string( TRIM( ocean%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( ocean%restart_filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( ocean%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    CALL add_time_dimension_to_file( ocean%restart_filename, ncid)

    ! Add a depth dimension to the file
    CALL add_depth_dimension_to_file( ocean%restart_filename, ncid, C%z_ocean)

    ! Add the data fields to the file
    CALL add_field_mesh_dp_3D_ocean( ocean%restart_filename, ncid, 'T' , long_name = 'Ocean temperatures', units = 'degrees C')
    CALL add_field_mesh_dp_3D_ocean( ocean%restart_filename, ncid, 'S' , long_name = 'Ocean salinity', units = 'PSU')

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_ocean_model_region

  SUBROUTINE remap_ocean_model( mesh_old, mesh_new, ice, ocean, region_name, time)
    ! Remap the ocean model

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_ocean_model'
    CHARACTER(LEN=256)                                    :: choice_ocean_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '    Remapping ocean model data to the new mesh...'

    ! Determine which ocean model to run for this region
    select case (region_name)
    case default
      call crash('unknown region_name "' // region_name // '"')
    case ('NAM')
      choice_ocean_model = C%choice_ocean_model_NAM
    case ('EAS')
      choice_ocean_model = C%choice_ocean_model_EAS
    case ('GRL')
      choice_ocean_model = C%choice_ocean_model_GRL
    case ('ANT')
      choice_ocean_model = C%choice_ocean_model_ANT
    end select

    ! Reallocate memory for main variables
    CALL reallocate_bounds( ocean%T, mesh_new%vi1, mesh_new%vi2, C%nz_ocean)
    CALL reallocate_bounds( ocean%S, mesh_new%vi1, mesh_new%vi2, C%nz_ocean)

    ! Reallocate memory for secondary variables
    CALL reallocate_bounds( ocean%T_draft,          mesh_new%vi1, mesh_new%vi2)
    CALL reallocate_bounds( ocean%T_freezing_point, mesh_new%vi1, mesh_new%vi2)

    ! Determine which ocean model to remap
    select case (choice_ocean_model)
    case default
      call crash('unknown choice_ocean_model "' // TRIM( choice_ocean_model) // '"')
    case ('none')
      ! No need to do anything
    case ('idealised')
      call initialise_ocean_model_idealised( mesh_new, ocean)
    case ('realistic')
      call remap_ocean_model_realistic( mesh_old, mesh_new, ice, ocean, region_name, time)
    end select

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_ocean_model

  SUBROUTINE limit_ocean_supercooling( mesh, ice, ocean)
  ! Calculate the ocean freezing point at 6 km depth, taking it as the coldest temperature we will allow the ocean to be
    ! so we can prevent unrealistically high refreezing due to uniformly applied deltaT values that might not be representative of a certain region

    implicit none

    ! In/output variables
    type(type_mesh),                    intent(in)    :: mesh
    type(type_ice_model),               intent(in)    :: ice
    type(type_ocean_model),             intent(inout) :: ocean

    ! Local variables:
    character(len=256), parameter                     :: routine_name = 'limit_ocean_supercooling'
    integer                                           :: vi, z
    real(dp), parameter                               :: depth_limit = 6E3_dp ! depth of freezing point that we want to limit the temperature to be
    real(dp)                                          :: S0                   ! Practical salinity [PSU]
    real(dp), parameter                               :: lambda1 = -0.0575_dp ! Liquidus slope                [degC PSU^-1] (Favier et al. (2019), Table 2)
    real(dp), parameter                               :: lambda2 = 0.0832_dp  ! Liquidus intercept            [degC]        (Favier et al. (2019), Table 2)
    real(dp), parameter                               :: lambda3 = 7.59E-4_dp ! Liquidus pressure coefficient [degC m^-1]   (Favier et al. (2019), Table 2)
    real(dp)                                          :: T_supercooled_limit

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
    do z = 1, C%nz_ocean
      ! Find salinity at 6 km depth
      call interpolate_ocean_depth( C%nz_ocean, C%z_ocean, ocean%S( vi,:), depth_limit, S0)
      ! Calculate ocean freezing temperature (Favier et al. (2019), Eq. 3) in degrees Celsius
      T_supercooled_limit = lambda1 * S0 + lambda2 - lambda3 * depth_limit

      if (ocean%T(vi, z) < T_supercooled_limit) then
        ocean%T(vi, z) = T_supercooled_limit
      end if
    end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  END SUBROUTINE limit_ocean_supercooling

END MODULE ocean_main
