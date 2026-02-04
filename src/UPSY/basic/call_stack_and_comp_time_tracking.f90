MODULE call_stack_and_comp_time_tracking

  ! Keep track of which subroutine has control, so that error messages can actually tell
  ! you where they originated.
  !
  ! Also keep track of how much computation time each subroutine uses.

! ===== Preamble =====
! ====================

  use basic_program_info, only: program_name, routine_path, n_MPI_windows_used
  use mpi_f08, only: MPI_WTIME
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use string_module, only: colour_string, insert_val_into_string_dp, insert_val_into_string_int
  use crash_mod, only: crash, warning, happy

  IMPLICIT NONE

! ===== Global variables =====
! ============================

  LOGICAL             :: do_display_path   = .FALSE.

  TYPE subroutine_resource_tracker
    character(len=:), allocatable :: routine_path
    real(dp)                      :: tstart, tcomp
    integer                       :: n_MPI_windows_used_initial
    integer                       :: n_MPI_windows_used_final
  END TYPE subroutine_resource_tracker

  TYPE( subroutine_resource_tracker), DIMENSION(:), ALLOCATABLE :: resource_tracker

CONTAINS

! ===== Control & resource tracking ======
! ========================================

  SUBROUTINE initialise_control_and_resource_tracker
    ! Initialise the control and resource tracker

#if (DO_RESOURCE_TRACKING)

    ! Local variables:
    INTEGER :: i

    ALLOCATE( resource_tracker( 4000))

    ! Initialise values
    DO i = 1, size( resource_tracker,1)
      resource_tracker( i)%routine_path = 'subroutine_placeholder'
      resource_tracker( i)%tstart       = 0._dp
      resource_tracker( i)%tcomp        = 0._dp
    END DO

#endif

    ! Initialise the routine path
    routine_path = program_name

    ! Initialise the shared memory leak tracker
    n_MPI_windows_used = 0

  END SUBROUTINE initialise_control_and_resource_tracker

  SUBROUTINE init_routine( routine_name, do_track_resource_use)
    ! Initialise a subroutine; update the routine path

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    CHARACTER(LEN=*)                                 , INTENT(IN)    :: routine_name
    LOGICAL                                  , OPTIONAL, INTENT(IN)    :: do_track_resource_use

#if (DO_RESOURCE_TRACKING)

    ! Local variables:
    INTEGER :: len_path, len_name
    INTEGER :: i
    LOGICAL :: do_track_resource_use_loc

    len_path = LEN(      routine_path)
    len_name = LEN_TRIM( routine_name)

    ! Append this routine to the routine path
    routine_path = routine_path // '/' // trim( routine_name)

    ! If so specified, print the current routine path to the terminal (useful for debugging)
    IF (do_display_path) THEN
      IF (par%primary) WRITE(0,'(A)') '   Initialising ' // routine_path
    END IF

    ! Check if resource use for this subroutine should be tracked
    ! (especially for the NetCDF routines we don't want to do this, as there are
    ! a great many of them and the resource tracker output file will become annoyingly big)

    IF (PRESENT( do_track_resource_use)) THEN
      do_track_resource_use_loc = do_track_resource_use
    ELSE
      do_track_resource_use_loc = .TRUE.
    END IF

    ! Don't use the resource tracker for the test programs
    ! (as those use wayyy more combinations than the actual model,
    !  quickly overflowing the tracker)
    if (index( routine_path, 'UPSY_unit_tests') /= 0 .or. &
        index( routine_path, 'UPSY_multinode_unit_tests') /= 0 .or. &
        index( routine_path, 'UPSY_component_tests') /= 0 .or. &
        index( routine_path, 'UFEMISM_program/run_all_unit_tests') /= 0.or. &
        index( routine_path, 'UFEMISM_program/run_all_component_tests') /= 0.or. &
        index( routine_path, 'LADDIE/run_laddie_unit_tests') /= 0) then
      do_track_resource_use_loc = .false.
    end if

    IF (do_track_resource_use_loc) THEN

      ! Initialise the computation time tracker
      CALL find_subroutine_in_resource_tracker( i)
      resource_tracker( i)%tstart = MPI_WTIME()

      ! Check how many MPI windows are in use when this subroutine is initialised
      resource_tracker( i)%n_MPI_windows_used_initial = n_MPI_windows_used

    ELSE

      routine_path = routine_path // '_NOTRACK'

    END IF

#endif

  END SUBROUTINE init_routine

  SUBROUTINE finalise_routine( routine_name, n_extra_MPI_windows_expected)
    ! Finalise; remove the current routine name from the routine path

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    CHARACTER(LEN=*)                                 , INTENT(IN)    :: routine_name
    integer, optional, intent(in) :: n_extra_MPI_windows_expected

#if (DO_RESOURCE_TRACKING)

    ! Local variables:
    LOGICAL                                                            :: do_track_resource_use
    INTEGER                                                            :: len_path_tot, i
    REAL(dp)                                                           :: dt
    integer :: n_extra_MPI_windows_expected_

    ! If so specified, print the current routine path to the terminal (useful for debugging)
    IF (do_display_path) THEN
      IF (par%primary) WRITE(0,'(A)') '   Finalising   ' // routine_path
    END IF

    ! Check if resource use should be tracked for this subroutine
    i = INDEX( routine_path, '_NOTRACK')
    IF ( i == 0) THEN
      do_track_resource_use = .TRUE.
    ELSE
      do_track_resource_use = .FALSE.
    END IF

    IF (do_track_resource_use) THEN
      ! Resource use for this subroutine should be tracked

      ! Add computation time to the resource tracker
      CALL find_subroutine_in_resource_tracker( i)
      dt = MPI_WTIME() - resource_tracker( i)%tstart
      resource_tracker( i)%tcomp = resource_tracker( i)%tcomp + dt

      ! Check how many MPI windows are in use when this subroutine is finalised
      resource_tracker( i)%n_MPI_windows_used_final = n_MPI_windows_used

      if (present( n_extra_MPI_windows_expected)) then
        n_extra_MPI_windows_expected_ = n_extra_MPI_windows_expected
      else
        n_extra_MPI_windows_expected_ = 0
      end if

      ! If it is larger than at the start, mention this
      if (resource_tracker( i)%n_MPI_windows_used_final > &
        resource_tracker( i)%n_MPI_windows_used_initial + n_extra_MPI_windows_expected_) then

        ! Some exceptions where extra used MPI windows should not be interpreted as a memory leak
        if (index( routine_path, 'UFEMISM/initialise') == 0 .and. &
            index( routine_path, 'LADDIE_program/initialise') == 0 .and. &
            index( routine_path, 'allocate_dist_shared') == 0 .and. &
            index( routine_path, 'UPSY_unit_tests') == 0 .and. &
            index( routine_path, 'UPSY_multinode_unit_tests') == 0 .and. &
            index( routine_path, 'UPSY_component_tests') == 0 .and. &
            index( routine_path, 'run_all_unit_tests') == 0 .and. &
            index( routine_path, 'run_laddie_unit_tests') == 0 .and. &
            index( routine_path, 'run_all_multinode_unit_tests') == 0 .and. &
            index( routine_path, 'run_all_component_tests') == 0) then

          call warning('shared memory was allocated but not freed, possibly memory leak ' // &
            '(n_init = {int_01}, n_final = {int_02})', &
            int_01 = resource_tracker( i)%n_MPI_windows_used_initial, &
            int_02 = resource_tracker( i)%n_MPI_windows_used_final)
        end if
      end if

      ! Find where in the string exactly the current routine name is located
      i = INDEX( routine_path, trim( routine_name), back = .true.)

      IF (i == 0) THEN
        CALL crash('BEEP - finalise_routine - routine_name "' // TRIM( routine_name) // &
          '" not found in routine_path = "' // routine_path // '"!')
      END IF

      ! Remove the current routine name from the routine path
      routine_path = routine_path( 1:i-2)

    ELSE ! IF (do_track_resource_use) THEN
      ! Resource use for this subroutine should not be tracked

      ! Find where in the string exactly the current routine name is located
      i = INDEX( routine_path, TRIM( routine_name) // '_NOTRACK', back = .true.)

      IF (i == 0) THEN
        CALL crash('finalise_routine - routine_name "' // TRIM( routine_name) // &
          '" not found in routine_path = "' // routine_path // '"!')
      END IF

      ! Remove the current routine name from the routine path
      routine_path = routine_path( 1:i-2)

    END IF ! IF (do_track_resource_use) THEN

#endif

  END SUBROUTINE finalise_routine

  SUBROUTINE find_subroutine_in_resource_tracker( i)
    ! Find the current subroutine in the resource tracker. If it's not there yet, add it.

    IMPLICIT NONE

    ! In/output variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER                                            , INTENT(OUT)   :: i

#if (DO_RESOURCE_TRACKING)

    ! Local variables:
    INTEGER                                                            :: n

    n = SIZE( resource_tracker)

    DO i = 1, n
      IF     (resource_tracker( i)%routine_path == routine_path) THEN
        ! The current subroutine is listed at this position in the resource tracker
        RETURN
      ELSEIF (resource_tracker( i)%routine_path == 'subroutine_placeholder') THEN
        ! We've checked all listed subroutines and haven't found the current one; add it
        resource_tracker( i)%routine_path = routine_path
        RETURN
      END IF
    END DO

    ! If we've reached this point, then the resource tracker is overflowing
    CALL crash('Resource tracker overflows! Allocate more memory for it in call_stack_and_comp_time_tracking/initialise_control_and_resource_tracker!')

#endif

  END SUBROUTINE find_subroutine_in_resource_tracker

  SUBROUTINE reset_resource_tracker
    ! Reset the computation times and maximum memory use for all subroutines in the resource tracker

    IMPLICIT NONE

#if (DO_RESOURCE_TRACKING)

    ! Local variables:
    ! REAL(dp), DIMENSION(:,:,:,:), ALLOCATABLE, OPTIONAL, INTENT(INOUT) :: i
    INTEGER                                                            :: i,n


    n = SIZE( resource_tracker)

    DO i = 1, n
      resource_tracker( i)%tstart      = 0._dp
      resource_tracker( i)%tcomp       = 0._dp
    END DO

#endif

  END SUBROUTINE reset_resource_tracker

END MODULE call_stack_and_comp_time_tracking
