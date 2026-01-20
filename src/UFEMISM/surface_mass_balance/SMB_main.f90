module SMB_main

  ! The main SMB model module.

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use mpi_f08, only: MPI_WIN
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use grid_basic, only: type_grid
  use ice_model_types, only: type_ice_model
  use climate_model_types, only: type_climate_model
  use SMB_idealised, only: type_SMB_model_idealised
  ! use SMB_prescribed, only: type_SMB_model_prescribed
  ! use SMB_reconstructed, only: type_SMB_model_reconstructed
  ! use SMB_IMAU_ITM, only: type_SMB_model_IMAU_ITM
  ! use SMB_snapshot_plus_anomalies, only: type_SMB_model_snapshot_plus_anomalies
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use reallocate_dist_shared_mod, only: reallocate_dist_shared
  use mesh_ROI_polygons, only: calc_polygon_Patagonia
  use plane_geometry, only: is_in_polygon
  use mesh_data_smoothing, only: smooth_Gaussian
  use netcdf_io_main

  implicit none

  type type_SMB_model
    ! The surface mass balance model

    ! Main data fields
    real(dp), dimension(:), contiguous, pointer :: SMB                       ! Yearly  SMB (m)
    type(MPI_WIN) :: wSMB

    ! Sub-models
    type(type_SMB_model_idealised)               :: idealised
    ! type(type_SMB_model_prescribed)              :: prescribed
    ! type(type_SMB_model_reconstructed)           :: reconstructed
    ! type(type_SMB_model_IMAU_ITM)                :: IMAUITM
    ! type(type_SMB_model_snapshot_plus_anomalies) :: snapshot_plus_anomalies

    ! Timestepping
    real(dp)                                     :: t_next

    ! Metadata
    character(:), allocatable                    :: restart_filename          ! Name for generated restart file

  end type type_SMB_model

contains

! ===== Main routines =====
! =========================

  SUBROUTINE run_SMB_model( mesh, grid_smooth, ice, climate, SMB, region_name, time)
    ! Calculate the surface mass balance

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_grid),                        INTENT(IN)    :: grid_smooth
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(IN)    :: climate
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model
    integer                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if we need to calculate a new SMB
    IF (C%do_asynchronous_SMB) THEN
      ! Asynchronous coupling: do not calculate a new SMB in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next SMB time step
      IF (time == SMB%t_next) THEN
        ! Go on to calculate a new SMB
        SMB%t_next = time + C%dt_SMB
      ELSEIF (time > SMB%t_next) THEN
        ! This should not be possible
        CALL crash('overshot the SMB time step')
      ELSE
        ! It is not yet time to calculate a new SMB
        CALL finalise_routine( routine_name)
        RETURN
      END IF

    ELSE ! IF (C%do_asynchronous_SMB) THEN
      ! Synchronous coupling: calculate a new SMB in every model loop
      SMB%t_next = time + C%dt_SMB
    END IF

    ! Determine which SMB model to run for this region
    SELECT CASE (region_name)
    CASE ('NAM')
      choice_SMB_model = C%choice_SMB_model_NAM
    CASE ('EAS')
      choice_SMB_model = C%choice_SMB_model_EAS
    CASE ('GRL')
      choice_SMB_model = C%choice_SMB_model_GRL
    CASE ('ANT')
      choice_SMB_model = C%choice_SMB_model_ANT
    CASE DEFAULT
      CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Run the chosen SMB model
    SELECT CASE (choice_SMB_model)
    CASE DEFAULT
      CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')

    CASE ('uniform')
      SMB%SMB( mesh%vi1: mesh%vi2) = C%uniform_SMB

    CASE ('idealised')
      call SMB%idealised%run( SMB%idealised%ct_run( &
        ice, climate, grid_smooth, time, region_name))
      do vi = mesh%vi1, mesh%vi2
        SMB%SMB( vi) = SMB%idealised%SMB( vi)
      end do

    ! CASE ('prescribed')
    !   call SMB%prescribed%run( mesh, region_name, time)
    !   do vi = mesh%vi1, mesh%vi2
    !     SMB%SMB( vi) = SMB%prescribed%SMB( vi)
    !   end do

    ! CASE ('reconstructed')
    !   call SMB%reconstructed%run( mesh, grid_smooth, ice, region_name, time)
    !   do vi = mesh%vi1, mesh%vi2
    !     SMB%SMB( vi) = SMB%reconstructed%SMB( vi)
    !   end do

    ! CASE ('IMAU-ITM')
    !   call SMB%IMAUITM%run( mesh, ice, climate)
    !   do vi = mesh%vi1, mesh%vi2
    !     SMB%SMB( vi) = sum( SMB%IMAUITM%SMB_monthly( vi,:))
    !   end do

    ! case ('snapshot_plus_anomalies')
    !   call SMB%snapshot_plus_anomalies%run( mesh, time)
    !   do vi = mesh%vi1, mesh%vi2
    !     SMB%SMB( vi) = SMB%snapshot_plus_anomalies%SMB( vi)
    !   end do

    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model

  SUBROUTINE initialise_SMB_model( mesh, ice, SMB, region_name)
    ! Initialise the SMB model

    ! In- and output variables
    TYPE(type_mesh),                        intent(in   ) :: mesh
    TYPE(type_ice_model),                   intent(in   ) :: ice
    TYPE(type_SMB_model),                   intent(inout) :: SMB
    CHARACTER(LEN=3),                       intent(in   ) :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising surface mass balance model...'

    ! Determine which SMB model to initialise for this region
    SELECT CASE (region_name)
    CASE ('NAM')
      choice_SMB_model = C%choice_SMB_model_NAM
    CASE ('EAS')
      choice_SMB_model = C%choice_SMB_model_EAS
    CASE ('GRL')
      choice_SMB_model = C%choice_SMB_model_GRL
    CASE ('ANT')
      choice_SMB_model = C%choice_SMB_model_ANT
    CASE DEFAULT
      CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    IF (par%primary)  WRITE(*,"(A)") '   Initialising SMB model ' // &
      colour_string( trim( choice_SMB_model),'light blue') // '...'

    ! Allocate memory for main variables
    call allocate_dist_shared( SMB%SMB, SMB%wSMB, mesh%pai_V%n_nih)
    SMB%SMB( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => SMB%SMB

    ! Set time of next calculation to start time
    SMB%t_next = C%start_time_of_run

    ! Determine which SMB model to initialise
    SELECT CASE (choice_SMB_model)
    CASE ('uniform')
      SMB%SMB( mesh%vi1: mesh%vi2) = C%uniform_SMB
    CASE ('idealised')
      call SMB%idealised%allocate( SMB%idealised%ct_allocate( mesh))
      call SMB%idealised%initialise( SMB%idealised%ct_initialise( region_name))
    ! CASE ('prescribed')
    !   call SMB%prescribed%init( mesh, region_name)
    ! CASE ('reconstructed')
    !   call SMB%reconstructed%init( mesh)
    ! CASE ('IMAU-ITM')
    !   call SMB%IMAUITM%init( mesh, ice, region_name)
    ! case ('snapshot_plus_anomalies')
    !   call SMB%snapshot_plus_anomalies%init( mesh)
    CASE DEFAULT
      CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SMB_model

  ! SUBROUTINE write_to_restart_file_SMB_model( mesh, SMB, region_name, time)
  !   ! Write to the restart file for the SMB model

  !   IMPLICIT NONE

  !   ! In/output variables:
  !   TYPE(type_mesh),                        INTENT(IN)    :: mesh
  !   TYPE(type_SMB_model),                   INTENT(IN)    :: SMB
  !   CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
  !   REAL(dp),                               INTENT(IN)    :: time

  !   ! Local variables:
  !   CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'write_to_restart_file_SMB_model'
  !   CHARACTER(LEN=256)                                    :: choice_SMB_model

  !   ! Add routine to path
  !   CALL init_routine( routine_name)

  !   ! Determine which SMB model to use for this region
  !   SELECT CASE (region_name)
  !   CASE ('NAM')
  !     choice_SMB_model = C%choice_SMB_model_NAM
  !   CASE ('EAS')
  !     choice_SMB_model = C%choice_SMB_model_EAS
  !   CASE ('GRL')
  !     choice_SMB_model = C%choice_SMB_model_GRL
  !   CASE ('ANT')
  !     choice_SMB_model = C%choice_SMB_model_ANT
  !   CASE DEFAULT
  !     CALL crash('unknown region_name "' // region_name // '"')
  !   END SELECT

  !   ! Write to the restart file of the chosen SMB model
  !   SELECT CASE (choice_SMB_model)
  !   CASE DEFAULT
  !     CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
  !   CASE ('uniform', &
  !         'idealised', &
  !         'prescribed', &
  !         'reconstructed', &
  !         'snapshot_plus_anomalies')
  !     ! No need to do anything
  !   CASE ('IMAU-ITM')
  !     call write_to_restart_file_SMB_model_region(mesh, SMB, region_name, time)
  !   END SELECT

  !   ! Finalise routine path
  !   CALL finalise_routine( routine_name)

  ! END SUBROUTINE write_to_restart_file_SMB_model

  ! SUBROUTINE write_to_restart_file_SMB_model_region(mesh, SMB, region_name, time)
  !   ! Write to the restart NetCDF file for the SMB model

  !   ! In/output variables:
  !   TYPE(type_mesh),          INTENT(IN) :: mesh
  !   TYPE(type_SMB_model),     INTENT(IN) :: SMB
  !   CHARACTER(LEN=3),         INTENT(IN) :: region_name
  !   REAL(dp),                 INTENT(IN) :: time

  !   ! Local variables:
  !   CHARACTER(LEN=256), PARAMETER        :: routine_name = 'write_to_restart_file_SMB_model_region'
  !   INTEGER                              :: ncid

  !   ! Add routine to path
  !   CALL init_routine( routine_name)

  !   ! If no NetCDF output should be created, do nothing
  !   IF (.NOT. C%do_create_netcdf_output) THEN
  !     CALL finalise_routine( routine_name)
  !     RETURN
  !   END IF

  !   ! Print to terminal
  !   IF (par%primary) WRITE(0,'(A)') '   Writing to SMB restart file "' // &
  !     colour_string( TRIM( SMB%restart_filename), 'light blue') // '"...'

  !   ! Open the NetCDF file
  !   CALL open_existing_netcdf_file_for_writing( SMB%restart_filename, ncid)

  !   ! Write the time to the file
  !   CALL write_time_to_file( SMB%restart_filename, ncid, time)
  !   ! month dimension is already written when adding it to file

  !   ! ! Write the SMB fields to the file
  !   ! TODO: do we need to check if IMAUITM is being used before writing these files?
  !   CALL write_to_field_multopt_mesh_dp_2D_monthly( mesh, SMB%restart_filename, ncid, 'SMB_monthly', SMB%IMAUITM%SMB_monthly)
  !   CALL write_to_field_multopt_mesh_dp_2D_monthly( mesh, SMB%restart_filename, ncid, 'FirnDepth', SMB%IMAUITM%FirnDepth)
  !   CALL write_to_field_multopt_mesh_dp_2D(         mesh, SMB%restart_filename, ncid, 'MeltPreviousYear', SMB%IMAUITM%MeltPreviousYear)
  !   CALL write_to_field_multopt_mesh_dp_2D(         mesh, SMB%restart_filename, ncid, 'SMB', SMB%SMB)

  !   ! Close the file
  !   CALL close_netcdf_file( ncid)

  !   ! Finalise routine path
  !   CALL finalise_routine( routine_name)

  ! END SUBROUTINE write_to_restart_file_SMB_model_region

  ! SUBROUTINE create_restart_file_SMB_model( mesh, SMB, region_name)
  !   ! Create the restart file for the SMB model

  !   IMPLICIT NONE

  !   ! In/output variables:
  !   TYPE(type_mesh),                        INTENT(IN)    :: mesh
  !   TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
  !   CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

  !   ! Local variables:
  !   CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_file_SMB_model'
  !   CHARACTER(LEN=256)                                    :: choice_SMB_model

  !   ! Add routine to path
  !   CALL init_routine( routine_name)

  !   ! Determine which SMB model to use for this region
  !   SELECT CASE (region_name)
  !   CASE ('NAM')
  !     choice_SMB_model = C%choice_SMB_model_NAM
  !   CASE ('EAS')
  !     choice_SMB_model = C%choice_SMB_model_EAS
  !   CASE ('GRL')
  !     choice_SMB_model = C%choice_SMB_model_GRL
  !   CASE ('ANT')
  !     choice_SMB_model = C%choice_SMB_model_ANT
  !   CASE DEFAULT
  !     CALL crash('unknown region_name "' // region_name // '"')
  !   END SELECT

  !   ! Create the restart file of the chosen SMB model
  !   SELECT CASE (choice_SMB_model)
  !   CASE ('uniform', &
  !         'idealised', &
  !         'prescribed', &
  !         'reconstructed', &
  !         'snapshot_plus_anomalies')
  !     ! No need to do anything
  !   CASE ('IMAU-ITM')
  !     call create_restart_file_SMB_model_region(mesh, SMB, region_name)
  !   CASE DEFAULT
  !     CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
  !   END SELECT

  !   ! Finalise routine path
  !   CALL finalise_routine( routine_name)

  ! END SUBROUTINE create_restart_file_SMB_model

  ! SUBROUTINE create_restart_file_SMB_model_region(mesh, SMB, region_name)
  !   ! Create the restart file for the SMB model

  !   IMPLICIT NONE

  !   ! In/output variables:
  !   TYPE(type_mesh),                        INTENT(IN)    :: mesh
  !   TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
  !   CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

  !   ! Local variables:
  !   CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_file_SMB_model'
  !   CHARACTER(LEN=256)                                    :: choice_SMB_model
  !   CHARACTER(LEN=256)                                    :: filename_base
  !   INTEGER                                               :: ncid

  !   ! Add routine to path
  !   CALL init_routine( routine_name)

  !   ! If no NetCDF output should be created, do nothing
  !   IF (.NOT. C%do_create_netcdf_output) THEN
  !     CALL finalise_routine( routine_name)
  !     RETURN
  !   END IF

  !   ! Set the filename
  !   filename_base = TRIM( C%output_dir) // 'restart_SMB_' // region_name
  !   CALL generate_filename_XXXXXdotnc( filename_base, SMB%restart_filename)

  !   ! Print to terminal
  !   IF (par%primary) WRITE(0,'(A)') '   Creating SMB model restart file "' // &
  !     colour_string( TRIM( SMB%restart_filename), 'light blue') // '"...'

  !   ! Create the NetCDF file
  !   CALL create_new_netcdf_file_for_writing( SMB%restart_filename, ncid)

  !   ! Set up the mesh in the file
  !   CALL setup_mesh_in_netcdf_file( SMB%restart_filename, ncid, mesh)

  !   ! Add a time dimension to the file
  !   CALL add_month_dimension_to_file( SMB%restart_filename, ncid)
  !   CALL add_time_dimension_to_file( SMB%restart_filename, ncid)


  !   ! Add the data fields to the file
  !   CALL add_field_mesh_dp_2D( SMB%restart_filename, ncid, 'SMB',                 long_name = 'Surface mass balance',            units = 'm/yr')
  !   CALL add_field_mesh_dp_2D_monthly( SMB%restart_filename, ncid, 'SMB_monthly', long_name = 'monthly Surface mass balance',    units = 'm/yr')
  !   CALL add_field_mesh_dp_2D_monthly( SMB%restart_filename, ncid, 'FirnDepth',   long_name = 'Firn Depth',                      units = 'm')
  !   CALL add_field_mesh_dp_2D( SMB%restart_filename, ncid, 'MeltPreviousYear',    long_name = 'Total melt in the previous year', units = 'm w.e.')

  !   ! Close the file
  !   CALL close_netcdf_file( ncid)

  !   ! Finalise routine path
  !   CALL finalise_routine( routine_name)

  ! END SUBROUTINE create_restart_file_SMB_model_region

  SUBROUTINE remap_SMB_model( mesh_old, mesh_new, SMB, region_name)
    ! Remap the SMB model

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_SMB_model),                   INTENT(inout) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '    Remapping surface mass balance model data to the new mesh...'

    ! Determine which SMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Reallocate memory for main variables
    call reallocate_dist_shared( SMB%SMB, SMB%wSMB, mesh_new%pai_V%n_nih)

    ! Determine which SMB model to initialise
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        call SMB%idealised%remap( SMB%idealised%ct_remap( mesh_new))
      ! CASE ('prescribed')
      !   call SMB%prescribed%remap( mesh_new, region_name)
      ! CASE ('reconstructed')
      !   call SMB%reconstructed%remap( mesh_new)
      ! CASE ('IMAU-ITM')
        ! call SMB%IMAUITM%remap( mesh_old, mesh_new)
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_SMB_model

end module SMB_main
