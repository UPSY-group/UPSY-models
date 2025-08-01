MODULE BMB_laddie

  ! LADDIE1.0 model (python implementation)
  ! Called when BMB model config is 'laddie_py'
  ! Lambert et al. (2023) doi: 10.5194/tc-17-3203-2023
  ! To use this option, get the code here: https://github.com/erwinlambert/laddie

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE BMB_model_types                                        , ONLY: type_BMB_model
  use netcdf_io_main

  IMPLICIT NONE


CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_BMB_model_laddie( mesh, ice, BMB, time, do_hybrid)
    ! Calculate the basal mass balance
    !
    ! Call the external LADDIE model

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB
    REAL(dp),                               INTENT(IN)    :: time
    LOGICAL,                                INTENT(IN)    :: do_hybrid  ! If .TRUE. run laddie within ROI, and a different BMB model outside ROI. If .FALSE. run laddie for the full domain.

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_BMB_model_laddie'
    CHARACTER(LEN=256)                                    :: filename_BMB_laddie_output
    CHARACTER(LEN=256)                                    :: filename_laddieready
    LOGICAL                                               :: found_laddie_file
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2)                :: temporary_BMB
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (time > C%start_time_of_run) THEN
      ! Compute and read laddie values for current geometry

      ! Define filename of BMB output from LADDIE
      filename_BMB_laddie_output      = TRIM(C%fixed_output_dir) // '/laddie_output/BMB_latest_UFE.nc'

      ! Run LADDIE
      IF (par%primary) THEN
        ! Different commands are needed to run laddie on different systems. local_mac or slurm_HPC (the latter is used for snellius)
        IF (C%choice_BMB_laddie_system == 'local_mac') THEN
          CALL system('cd ' // TRIM(C%dir_BMB_laddie_model) // ';' // TRIM(C%conda_activate_prompt) // '; python3 runladdie.py ' // TRIM(C%filename_BMB_laddie_configname) // '; conda deactivate')
        ELSEIF (C%choice_BMB_laddie_system == 'slurm_HPC') THEN
          CALL system('cd ' // TRIM(C%dir_BMB_laddie_model) // '; srun --ntasks=1 --exact --overlap --cpu-bind=cores python3 runladdie.py ' // TRIM(C%filename_BMB_laddie_configname))
        ELSE
          CALL crash('C%BMB_laddie_system not recognized, should be "local_mac" or "slurm_HPC".')
        END IF
      END IF

      ! Other processes wait for primary process to finish
      CALL sync

      ! Let UFEMISM sleep until LADDIE is finished
      CALL wait_for_laddie_to_finish( filename_laddieready, found_laddie_file)

      IF (do_hybrid) THEN
        ! Only apply values to ROI
        CALL read_field_from_file_2D( filename_BMB_laddie_output, 'BMBext', mesh, C%output_dir, temporary_BMB)

        DO vi = mesh%vi1, mesh%vi2
          IF (ice%mask_ROI(vi)) THEN
            ! Initialise all values in ROI at zero
            BMB%BMB( vi) = 0._dp

            ! Copy values from temporary BMB to all floating / grounding cells within ROI
            IF (ice%mask_floating_ice( vi) .OR. ice%mask_gl_fl( vi) .OR. ice%mask_gl_gr( vi)) THEN
              BMB%BMB_shelf( vi) = temporary_BMB(vi)
              BMB%BMB( vi)       = BMB%BMB_shelf( vi)
            END IF

          END IF
        END DO

      ELSE
        ! Apply values to all ice shelf grid cells
        CALL read_field_from_file_2D( filename_BMB_laddie_output, 'BMBext', mesh, C%output_dir, BMB%BMB_shelf)

      END IF

    ELSE ! (time = C%start_time_of_run)
      ! Read values from laddie initial output

      IF (do_hybrid) THEN
        ! Only apply values to ROI
        CALL read_field_from_file_2D( C%filename_BMB_laddie_initial_output, 'BMBext', mesh, C%output_dir, temporary_BMB)

        DO vi = mesh%vi1, mesh%vi2
          IF (ice%mask_ROI(vi)) THEN
            ! Initialise all values in ROI at zero
            BMB%BMB( vi) = 0._dp

            ! Copy values from temporary BMB to all floating / grounding cells within ROI
            IF (ice%mask_floating_ice( vi) .OR. ice%mask_gl_fl( vi) .OR. ice%mask_gl_gr( vi)) THEN
              BMB%BMB_shelf( vi) = temporary_BMB(vi)
              BMB%BMB( vi)       = BMB%BMB_shelf( vi)
            END IF

          END IF
        END DO

      ELSE
        ! Apply values to all ice shelf grid cells
        CALL read_field_from_file_2D( C%filename_BMB_laddie_initial_output, 'BMBext', mesh, C%output_dir, BMB%BMB_shelf)

      END IF

    END IF ! (time > C%start_time_of_run)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_BMB_model_laddie

  SUBROUTINE initialise_BMB_model_laddie( mesh, BMB)
    ! Initialise the BMB model
    !
    ! Call the external LADDIE model

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_BMB_model_laddie'
    LOGICAL                                               :: file_exists

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (par%primary) THEN

      ! Check if LADDIE restartfile exists to copy into laddie_output directory
      INQUIRE( EXIST = file_exists, FILE = TRIM( C%filename_BMB_laddie_initial_restart))
      IF (.NOT. file_exists) THEN
        CALL crash('file "' // TRIM( C%filename_BMB_laddie_initial_restart) // '" not found!')
      END IF

      ! Create folder for laddie output in fixed_output_dir,
      CALL system('mkdir ' // TRIM(C%fixed_output_dir) // '/laddie_output')

      ! Copy initial restart file to restart_latest.nc
      CALL system('cp ' // TRIM(C%filename_BMB_laddie_initial_restart) // ' ' // TRIM(C%fixed_output_dir) // '/laddie_output/restart_latest.nc')

    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_BMB_model_laddie

  SUBROUTINE remap_BMB_model_laddie( mesh, BMB)
    ! Remap the BMB model
    !
    ! Simply reread the latest BMB output

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_BMB_model),                   INTENT(INOUT) :: BMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_BMB_model_laddie'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Read in BMB data from LADDIE
    CALL read_field_from_file_2D( C%filename_BMB_laddie_initial_output, 'BMBext', mesh, C%output_dir, BMB%BMB_shelf)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_BMB_model_laddie

  SUBROUTINE wait_for_laddie_to_finish( filename_laddieready, found_laddie_file)
    ! Let model sleep until a filename laddieready is detected
    !
    ! This filename is created by LADDIE when it's run is finished

    ! In/output variables:
    CHARACTER(LEN=256),                       INTENT(IN)  :: filename_laddieready
    LOGICAL,                                  INTENT(OUT) :: found_laddie_file

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'wait_for_laddie_to_finish'
    CHARACTER(LEN=256)                                    :: laddieready

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Define laddieready filename
    laddieready = TRIM(C%fixed_output_dir) // '/laddie_output/laddieready'   ! Dummy file to indicate LADDIE run is finished

    ! Check whether new LADDIE file is available, if not sleep
    found_laddie_file = .FALSE.

    ! Start sleep loop
    DO WHILE (.NOT. found_laddie_file)
      IF (par%primary) THEN
        print*, 'Looking for LADDIE file...'
      END IF
      CALL sync

      INQUIRE( EXIST = found_laddie_file, FILE = laddieready)

      CALL SLEEP(1)

    END DO ! End sleep loop if laddieready is found

    ! Remove laddieready file
    IF (par%primary) THEN
      CALL system('rm ' // TRIM(laddieready))
    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE wait_for_laddie_to_finish

END MODULE BMB_laddie
