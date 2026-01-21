module SMB_snapshot_plus_anomalies

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mesh_types, only: type_mesh
  use climate_model_types, only: type_climate_model
  use parameters, only: NaN
  use model_configuration, only: C
  use netcdf_io_main
  use mpi_f08, only: MPI_WIN, MPI_BCAST, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use SMB_basic, only: atype_SMB_model, type_SMB_model_context_allocate, &
    type_SMB_model_context_initialise, type_SMB_model_context_run, type_SMB_model_context_remap
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension

  implicit none

  private

  public :: type_SMB_model_snapshot_plus_anomalies

  type, extends(atype_SMB_model) :: type_SMB_model_snapshot_plus_anomalies
    !< The snapshot+anomalies SMB model (monthly T2m, annual SMB)

      ! Baseline climate
      real(dp), dimension(:,:), contiguous, pointer :: T2m_baseline => null()
      real(dp), dimension(:  ), contiguous, pointer :: SMB_baseline => null()
      type(MPI_WIN) :: wT2m_baseline, wSMB_baseline

      ! Two anomaly timeframes enveloping the current model time
      real(dp)                                      :: anomaly_t0
      real(dp), dimension(:  ), contiguous, pointer :: T2m_anomaly_0 => null()
      real(dp), dimension(:  ), contiguous, pointer :: SMB_anomaly_0 => null()
      type(MPI_WIN) :: wT2m_anomaly_0, wSMB_anomaly_0

      real(dp)                                      :: anomaly_t1
      real(dp), dimension(:  ), contiguous, pointer :: T2m_anomaly_1 => null()
      real(dp), dimension(:  ), contiguous, pointer :: SMB_anomaly_1 => null()
      type(MPI_WIN) :: wT2m_anomaly_1, wSMB_anomaly_1

      ! Time-weighted anomaly
      real(dp), dimension(:  ), contiguous, pointer :: T2m_anomaly => null()
      real(dp), dimension(:  ), contiguous, pointer :: SMB_anomaly => null()
      type(MPI_WIN) :: wT2m_anomaly, wSMB_anomaly

      ! Applied climate
      real(dp), dimension(:,:), contiguous, pointer :: T2m    ! = baseline + anomaly
      type(MPI_WIN) :: wT2m

    contains

      procedure, public :: allocate_SMB_model   => allocate_SMB_model_snapshot_plus_anomalies
      procedure, public :: initialise_SMB_model => initialise_SMB_model_snapshot_plus_anomalies
      procedure, public :: run_SMB_model        => run_SMB_model_snapshot_plus_anomalies_abs
      procedure, public :: remap_SMB_model      => remap_SMB_model_snapshot_plus_anomalies

      procedure, private :: run_SMB_model_snapshot_plus_anomalies
      procedure, private :: run_SMB_model_snapshot_plus_anomalies_climate
      procedure, private :: update_timeframes

  end type type_SMB_model_snapshot_plus_anomalies

contains

  subroutine allocate_SMB_model_snapshot_plus_anomalies( self, context)

    ! In/output variables:
    class(type_SMB_model_snapshot_plus_anomalies), intent(inout) :: self
    type(type_SMB_model_context_allocate),         intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SMB_model_snapshot_plus_anomalies'

    ! Add routine to path
    call init_routine( routine_name)

    call self%set_name('SMB_model_snapshot_plus_anomalies')

    ! Create all model fields
    ! =======================

    ! Baseline climate
    call self%create_field( self%T2m_baseline, self%wT2m_baseline, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m_baseline', &
      long_name = 'baseline monthly 2-m air temperature', &
      units     = 'K')
    call self%create_field( self%SMB_baseline, self%wSMB_baseline, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'SMB_baseline', &
      long_name = 'baseline surface mass balance', &
      units     = 'm yr^-1')

    ! Two anomaly snapshots enveloping the current model time
    call self%create_field( self%T2m_anomaly_0, self%wT2m_anomaly_0, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'T2m_anomaly_0', &
      long_name = 'previous annual 2-m air temperature anomaly', &
      units     = 'K')
    call self%create_field( self%SMB_anomaly_0, self%wSMB_anomaly_0, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'SMB_anomaly_0', &
      long_name = 'previous surface mass balance anomaly', &
      units     = 'm yr^-1')
    call self%create_field( self%T2m_anomaly_1, self%wT2m_anomaly_1, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'T2m_anomaly_1', &
      long_name = 'next annual 2-m air temperature anomaly', &
      units     = 'K')
    call self%create_field( self%SMB_anomaly_1, self%wSMB_anomaly_1, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'SMB_anomaly_1', &
      long_name = 'next surface mass balance anomaly', &
      units     = 'm yr^-1')

    ! Time-weighted anomaly
    call self%create_field( self%T2m_anomaly, self%wT2m_anomaly, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'T2m_anomaly', &
      long_name = 'annual 2-m air temperature anomaly', &
      units     = 'K')
    call self%create_field( self%SMB_anomaly, self%wSMB_anomaly, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'SMB_anomaly', &
      long_name = 'surface mass balance anomaly', &
      units     = 'm yr^-1')

    ! Applied climate
    call self%create_field( self%T2m, self%wT2m, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m', &
      long_name = 'monthly 2-m air temperature', &
      units     = 'K')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_snapshot_plus_anomalies

  subroutine initialise_SMB_model_snapshot_plus_anomalies( self, context)

    ! In/output variables:
    class(type_SMB_model_snapshot_plus_anomalies), intent(inout) :: self
    type(type_SMB_model_context_initialise),       intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_snapshot_plus_anomalies'

    ! Add routine to path
    call init_routine( routine_name)

    ! Read baseline snapshot
    call read_field_from_file_2D_monthly( C%SMB_snp_p_anml_filename_snapshot_T2m, 'T2m', &
      self%mesh, C%output_dir, self%T2m_baseline)
    call read_field_from_file_2D( C%SMB_snp_p_anml_filename_snapshot_SMB, 'SMB', &
      self%mesh, C%output_dir, self%SMB_baseline)

    ! Initialise anomaly timeframes
    self%anomaly_t0 = C%start_time_of_run - 200._dp
    self%anomaly_t1 = C%start_time_of_run - 100._dp
    call self%update_timeframes( self%mesh, C%start_time_of_run)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_snapshot_plus_anomalies

  subroutine run_SMB_model_snapshot_plus_anomalies_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_snapshot_plus_anomalies), intent(inout) :: self
    type(type_SMB_model_context_run),              intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_snapshot_plus_anomalies_abs'

    ! Add routine to path
    call init_routine( routine_name)

    call self%run_SMB_model_snapshot_plus_anomalies( self%mesh, context%time)
    call self%run_SMB_model_snapshot_plus_anomalies_climate( self%mesh, context%climate, context%time)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_snapshot_plus_anomalies_abs

  subroutine remap_SMB_model_snapshot_plus_anomalies( self, context)

    ! In/output variables:
    class(type_SMB_model_snapshot_plus_anomalies), intent(inout) :: self
    type(type_SMB_model_context_remap),            intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_SMB_model_snapshot_plus_anomalies'

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap all model fields
    ! ======================

    ! Baseline climate
    call self%remap_field( context%mesh_new, 'T2m_baseline', self%T2m_baseline)
    call self%remap_field( context%mesh_new, 'SMB_baseline', self%SMB_baseline)

    ! Two anomaly snapshots enveloping the current model time
    call self%remap_field( context%mesh_new, 'T2m_anomaly_0', self%T2m_anomaly_0)
    call self%remap_field( context%mesh_new, 'SMB_anomaly_0', self%SMB_anomaly_0)
    call self%remap_field( context%mesh_new, 'T2m_anomaly_1', self%T2m_anomaly_1)
    call self%remap_field( context%mesh_new, 'SMB_anomaly_1', self%SMB_anomaly_1)

    ! ! Time-weighted anomaly
    call self%remap_field( context%mesh_new, 'T2m_anomaly', self%T2m_anomaly)
    call self%remap_field( context%mesh_new, 'SMB_anomaly', self%SMB_anomaly)

    ! Applied climate
    call self%remap_field( context%mesh_new, 'T2m', self%T2m)

    ! Set the timestamps of the timeframes so that they
    ! will be updated the next time the model is run
    self%anomaly_t0 = C%start_time_of_run - 200._dp
    self%anomaly_t1 = C%start_time_of_run - 100._dp

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_snapshot_plus_anomalies

  subroutine run_SMB_model_snapshot_plus_anomalies( self, mesh, time)

    ! In/output variables:
    class(type_SMB_model_snapshot_plus_anomalies), intent(inout) :: self
    type(type_mesh),                               intent(in   ) :: mesh
    real(dp),                                      intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_snapshot_plus_anomalies'
    real(dp)                       :: w0, w1

    ! Add routine to path
    call init_routine( routine_name)

    ! If the current model time falls outside the enveloping window
    ! of the two timeframes that have been read, update them
    if (time < self%anomaly_t0 .or. &
        time > self%anomaly_t1) then
      call self%update_timeframes( mesh, time)
    end if

    ! Interpolate between the two timeframes to find the applied anomaly
    w0 = (self%anomaly_t1 - time) / &
         (self%anomaly_t1 - self%anomaly_t0)
    w1 = 1._dp - w0

    self%SMB_anomaly( mesh%vi1: mesh%vi2) = &
      w0 * self%SMB_anomaly_0( mesh%vi1: mesh%vi2) + &
      w1 * self%SMB_anomaly_1( mesh%vi1: mesh%vi2)

    ! Add anomaly to snapshot to find the applied SMB
    self%SMB( mesh%vi1: mesh%vi2) = &
      self%SMB_baseline( mesh%vi1: mesh%vi2) + &
      self%SMB_anomaly( mesh%vi1: mesh%vi2)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_snapshot_plus_anomalies

  subroutine run_SMB_model_snapshot_plus_anomalies_climate( self, mesh, climate, time)

    ! In/output variables:
    class(type_SMB_model_snapshot_plus_anomalies), intent(inout) :: self
    type(type_mesh),                               intent(in   ) :: mesh
    type(type_climate_model),                      intent(in   ) :: climate
    real(dp),                                      intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_snapshot_plus_anomalies_climate'
    real(dp)                       :: w0, w1
    integer                        :: m

    ! Add routine to path
    call init_routine( routine_name)

    ! If the current model time falls outside the enveloping window
    ! of the two timeframes that have been read, update them
    if (time < self%anomaly_t0 .or. &
        time > self%anomaly_t1) then
      call self%update_timeframes( mesh, time)
    end if

    ! Interpolate between the two timeframes to find the applied anomaly
    w0 = (self%anomaly_t1 - time) / &
         (self%anomaly_t1 - self%anomaly_t0)
    w1 = 1._dp - w0

    self%T2m_anomaly = &
      w0 * self%T2m_anomaly_0 + &
      w1 * self%T2m_anomaly_1

    ! Add anomaly to snapshot to find the applied temperature
    do m = 1, 12
      self%T2m( mesh%vi1:mesh%vi2,m) = &
        self%T2m_baseline( mesh%vi1:mesh%vi2,m) + &
        self%T2m_anomaly ( mesh%vi1:mesh%vi2  )
    end do

    ! ! Copy to climate model
    ! climate%T2m( mesh%vi1:mesh%vi2,:) = self%T2m( mesh%vi1:mesh%vi2,:)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_snapshot_plus_anomalies_climate

  subroutine update_timeframes( self, mesh, time)

    ! In/output variables:
    class(type_SMB_model_snapshot_plus_anomalies), intent(inout) :: self
    type(type_mesh),                               intent(in   ) :: mesh
    real(dp),                                      intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'update_timeframes'
    character(len=1024)                 :: filename
    integer                             :: ncid, id_dim_time, nt, id_var_time, ierr
    real(dp), dimension(:), allocatable :: time_from_file
    integer                             :: ti0, ti1

    ! Add routine to path
    call init_routine( routine_name)

    filename = trim( C%SMB_snp_p_anml_filename_anomalies)

    ! Read time variable from the file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call check_time( filename, ncid)
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)
    call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)
    allocate( time_from_file( nt))
    call read_var_primary( filename, ncid, id_var_time, time_from_file)
    call MPI_BCAST( time_from_file(:), nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call close_netcdf_file( ncid)

    ! Find the two timeframes
    if (time < time_from_file( 1)) then
      if (par%primary) call warning('model time before start of anomaly file; using first timeframe')
      ti0 = 1
      ti1 = 2
    elseif (time > time_from_file( size( time_from_file,1))) then
      if (par%primary) call warning('model time beyond end of anomaly file; using last timeframe')
      ti0 = size( time_from_file,1) - 1
      ti1 = size( time_from_file,1)
    else
      ti0 = 1
      ti1 = 2
      do while (time_from_file( ti1) < time .and. ti1 < size( time_from_file,1))
        ti0 = ti1
        ti1 = ti1 + 1
      end do
    end if

    self%anomaly_t0 = time_from_file( ti0)
    self%anomaly_t1 = time_from_file( ti1)

    ! Read the two timeframes
    call read_field_from_file_2D( filename, 'T2m_anomaly', &
      mesh, C%output_dir, self%T2m_anomaly_0, &
      time_to_read = self%anomaly_t0)
    call read_field_from_file_2D( filename, 'T2m_anomaly', &
      mesh, C%output_dir, self%T2m_anomaly_1, &
      time_to_read = self%anomaly_t1)
    call read_field_from_file_2D( filename, 'SMB_anomaly', &
      mesh, C%output_dir, self%SMB_anomaly_0, &
      time_to_read = self%anomaly_t0)
    call read_field_from_file_2D( filename, 'SMB_anomaly', &
      mesh, C%output_dir, self%SMB_anomaly_1, &
      time_to_read = self%anomaly_t1)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_timeframes

end module SMB_snapshot_plus_anomalies
