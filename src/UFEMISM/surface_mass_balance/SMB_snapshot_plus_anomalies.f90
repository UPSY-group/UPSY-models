module SMB_snapshot_plus_anomalies

  use mpi_basic, only: par
  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use mesh_types, only: type_mesh
  use SMB_model_basic, only: atype_SMB_model, type_SMB_model_context_allocate, &
    type_SMB_model_context_initialise, type_SMB_model_context_run, &
    type_SMB_model_context_remap
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use netcdf_io_main, only: open_existing_netcdf_file_for_reading, check_time, &
    inquire_dim_multopt, inquire_var_multopt, read_var_primary, close_netcdf_file, &
    field_name_options_time, read_field_from_file_2D, read_field_from_file_2D_monthly
  use mpi_f08, only: MPI_WIN, MPI_BCAST, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
  use climate_model_types, only: type_climate_model

  implicit none

  private

  public :: type_SMB_model_snp_p_anml

  type, extends(atype_SMB_model) :: type_SMB_model_snp_p_anml

    ! Baseline climate
    real(dp), dimension(:,:), contiguous, pointer :: T2m_baseline
    real(dp), dimension(:  ), contiguous, pointer :: SMB_baseline
    type(MPI_WIN) :: wT2m_baseline, wSMB_baseline

    ! Two anomaly timeframes enveloping the current model time
    real(dp)                                      :: anomaly_t0
    real(dp), dimension(:  ), contiguous, pointer :: T2m_anomaly_0
    real(dp), dimension(:  ), contiguous, pointer :: SMB_anomaly_0
    type(MPI_WIN) :: wT2m_anomaly_0, wSMB_anomaly_0

    real(dp)                                      :: anomaly_t1
    real(dp), dimension(:  ), contiguous, pointer :: T2m_anomaly_1
    real(dp), dimension(:  ), contiguous, pointer :: SMB_anomaly_1
    type(MPI_WIN) :: wT2m_anomaly_1, wSMB_anomaly_1

    ! Time-weighted anomaly
    real(dp), dimension(:  ), contiguous, pointer :: T2m_anomaly
    real(dp), dimension(:  ), contiguous, pointer :: SMB_anomaly
    type(MPI_WIN) :: wT2m_anomaly, wSMB_anomaly

    ! Applied climate
    real(dp), dimension(:,:), contiguous, pointer :: T2m    ! = baseline + anomaly
    type(MPI_WIN) :: wT2m

    contains

      procedure, public :: allocate_SMB_model   => allocate_SMB_model_snp_p_anml_abs
      procedure, public :: deallocate_SMB_model => deallocate_SMB_model_snp_p_anml_abs
      procedure, public :: initialise_SMB_model => initialise_SMB_model_snp_p_anml_abs
      procedure, public :: run_SMB_model        => run_SMB_model_snp_p_anml_abs
      procedure, public :: remap_SMB_model      => remap_SMB_model_snp_p_anml_abs

      procedure, private :: allocate_SMB_model_snp_p_anml
      procedure, private :: initialise_SMB_model_snp_p_anml
      procedure, private :: run_SMB_model_snp_p_anml
      procedure, private :: update_timeframes

  end type type_SMB_model_snp_p_anml

contains

  subroutine allocate_SMB_model_snp_p_anml_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_snp_p_anml),               intent(inout) :: self
    type(type_SMB_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SMB_model_snp_p_anml_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call self%allocate_SMB_model_snp_p_anml

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_snp_p_anml_abs

  subroutine deallocate_SMB_model_snp_p_anml_abs( self)

    ! In/output variables:
    class(type_SMB_model_snp_p_anml), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'deallocate_SMB_model_snp_p_anml'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_SMB_model_snp_p_anml_abs

  subroutine initialise_SMB_model_snp_p_anml_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_snp_p_anml),                 intent(inout) :: self
    type(type_SMB_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_snp_p_anml_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call self%initialise_SMB_model_snp_p_anml( self%mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_snp_p_anml_abs

  subroutine run_SMB_model_snp_p_anml_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_snp_p_anml),          intent(inout) :: self
    type(type_SMB_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_snp_p_anml_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call self%run_SMB_model_snp_p_anml( self%mesh, context%climate, context%time)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_snp_p_anml_abs

  subroutine remap_SMB_model_snp_p_anml_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_snp_p_anml),           intent(inout) :: self
    type(type_SMB_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_SMB_model_snp_p_anml_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap data fields
    call self%remap_field( context%mesh_new, 'T2m_baseline' , self%T2m_baseline)
    call self%remap_field( context%mesh_new, 'SMB_baseline' , self%SMB_baseline)
    call self%remap_field( context%mesh_new, 'T2m_anomaly_0', self%T2m_anomaly_0)
    call self%remap_field( context%mesh_new, 'SMB_anomaly_0', self%SMB_anomaly_0)
    call self%remap_field( context%mesh_new, 'T2m_anomaly_1', self%T2m_anomaly_1)
    call self%remap_field( context%mesh_new, 'SMB_anomaly_1', self%SMB_anomaly_1)
    call self%remap_field( context%mesh_new, 'T2m_anomaly'  , self%T2m_anomaly)
    call self%remap_field( context%mesh_new, 'SMB_anomaly'  , self%SMB_anomaly)
    call self%remap_field( context%mesh_new, 'T2m'          , self%T2m)

    ! Re-initialise and update timeframes
    call self%initialise_SMB_model_snp_p_anml( self%mesh)
    call self%update_timeframes( self%mesh, context%time)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_snp_p_anml_abs



  subroutine allocate_SMB_model_snp_p_anml( self)

    ! In/output variables:
    class(type_SMB_model_snp_p_anml), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SMB_model_snp_p_anml'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Baseline climate
    call self%create_field( self%T2m_baseline, self%wT2m_baseline, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m_baseline', &
      long_name = 'baseline monthly 2-m air temperature', &
      units     = 'K', &
      remap_method = 'reallocate')

    call self%create_field( self%SMB_baseline, self%wSMB_baseline, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'SMB_baseline', &
      long_name = 'baseline surface mass balance', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    ! Two anomaly timeframes enveloping the current model time
    call self%create_field( self%T2m_anomaly_0, self%wT2m_anomaly_0, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'T2m_anomaly_0', &
      long_name = 'first 2-m air temperature anomaly timeframe', &
      units     = 'K', &
      remap_method = 'reallocate')

    call self%create_field( self%SMB_anomaly_0, self%wSMB_anomaly_0, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'SMB_anomaly_0', &
      long_name = 'first surface mass balance anomaly timeframe', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%T2m_anomaly_1, self%wT2m_anomaly_1, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'T2m_anomaly_1', &
      long_name = 'second 2-m air temperature anomaly timeframe', &
      units     = 'K', &
      remap_method = 'reallocate')

    call self%create_field( self%SMB_anomaly_1, self%wSMB_anomaly_1, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'SMB_anomaly_1', &
      long_name = 'second surface mass balance anomaly timeframe', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    ! Time-weighted anomaly
    call self%create_field( self%T2m_anomaly, self%wT2m_anomaly, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'T2m_anomaly', &
      long_name = '2-m air temperature anomaly', &
      units     = 'K', &
      remap_method = 'reallocate')

    call self%create_field( self%SMB_anomaly, self%wSMB_anomaly, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'SMB_anomaly', &
      long_name = 'surface mass balance anomaly', &
      units     = 'm yr^-1', &
      remap_method = 'reallocate')

    ! Applied climate
    call self%create_field( self%T2m, self%wT2m, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m', &
      long_name = 'applied monthly 2-m air temperature', &
      units     = 'K', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_snp_p_anml

  subroutine initialise_SMB_model_snp_p_anml( self, mesh)

    ! In/output variables:
    class(type_SMB_model_snp_p_anml), intent(inout) :: self
    type(type_mesh),                  intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_snp_p_anml'

    ! Add routine to path
    call init_routine( routine_name)

    ! Read baseline snapshot
    call read_field_from_file_2D_monthly( C%SMB_snp_p_anml_filename_snapshot_T2m, 'T2m', &
      mesh, C%output_dir, self%T2m_baseline)
    call read_field_from_file_2D( C%SMB_snp_p_anml_filename_snapshot_SMB, 'SMB', &
      mesh, C%output_dir, self%SMB_baseline)

    ! Initialise anomaly timeframes
    self%anomaly_t0 = C%start_time_of_run - 200._dp
    self%anomaly_t1 = C%start_time_of_run - 100._dp
    call self%update_timeframes( mesh, C%start_time_of_run)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_snp_p_anml

  subroutine run_SMB_model_snp_p_anml( self, mesh, climate, time)

    ! In/output variables:
    class(type_SMB_model_snp_p_anml), intent(inout) :: self
    type(type_mesh),                  intent(in   ) :: mesh
    type(type_climate_model),         intent(inout) :: climate
    real(dp),                         intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_snp_p_anml'
    real(dp)                       :: w0, w1
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! If the current model time falls outside the enveloping window
    ! of the two timeframes that have been read, update them
    if (time < self%anomaly_t0 .or. time > self%anomaly_t1) then
      call self%update_timeframes( mesh, time)
    end if

    ! Interpolate between the two timeframes to find the applied anomaly
    w0 = (self%anomaly_t1 - time) / (self%anomaly_t1 - self%anomaly_t0)
    w1 = 1._dp - w0

    do vi = mesh%vi1, mesh%vi2

      ! Note that the baseline and the applied temperature are monthly, but the anomaly is annual
      self%T2m_anomaly( vi) = w0 * self%T2m_anomaly_0( vi) + w1 * self%T2m_anomaly_1( vi)
      self%SMB_anomaly( vi) = w0 * self%SMB_anomaly_0( vi) + w1 * self%SMB_anomaly_1( vi)

      ! Add anomaly to snapshot to find the applied temperature and SMB
      self%T2m( vi,:) = self%T2m_baseline( vi,:) + self%T2m_anomaly( vi)
      self%SMB( vi  ) = self%SMB_baseline( vi  ) + self%SMB_anomaly( vi)

      ! Copy T2m to climate model
      climate%T2m( vi,:) = self%T2m( vi,:)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_snp_p_anml

  subroutine update_timeframes( self, mesh, time)

    ! In/output variables:
    class(type_SMB_model_snp_p_anml), intent(inout) :: self
    type(type_mesh),                  intent(in   ) :: mesh
    real(dp),                         intent(in   ) :: time

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