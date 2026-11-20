module climate_snapshot_plus_anomalies
  ! This climate model takes a snapshot as in the other "snapshot-type" climate models, 
  ! and modified it via time-and-space-varying temperature and precipitation "anomalies" over it. 
  ! Anomalies are read from a file, and are a difference for temperatures and ratio for precipitation 
  ! (therefore avoiding negative values).
  ! In short, the climate field is computed as follows:
  ! Temperature   = snapshot(x,y,m) + dT(x,y,m,time)
  ! Precipitation = snapshot(x,y,m) * rP(x,y,m,time)
  !
  ! As in the other snapshot models, values are downscaled to account for the difference in elevation between the
  ! GCM ice surface from which climate variables come, and the modelled ice surface

  ! ===== Preamble =====
  ! ====================

  use precisions,                        only: dp
  use upsy_main,                         only: upsy
  use mpi_basic,                         only: par, sync
   use mpi_f08,                          only: MPI_DOUBLE_PRECISION, MPI_COMM_WORLD, MPI_BCAST
  use call_stack_and_comp_time_tracking, only: crash, init_routine, finalise_routine, warning
  use model_configuration,               only: C
  use parameters
  use mpi_f08,                           only: MPI_WIN
  use mesh_basic,                        only: type_mesh
  use ice_geometry_model_data,           only: atype_ice_geometry_model_data
  use reference_geometry_types,          only: type_reference_geometry
  use Arakawa_grid_mod,                  only: Arakawa_grid
  use fields_dimensions,                 only: third_dimension
  use global_forcings_main
  use netcdf_io_main
  use netcdf_basic
  use reallocate_mod,                    only: reallocate_bounds
  use climate_model_basic,               only: atype_climate_model
  use series_utilities,                  only: update_timeframes_from_record, interpolate_value_from_forcing_record

  implicit none

  private
  public :: type_climate_model_snapshot_plus_anomalies

    type, extends(atype_climate_model) :: type_climate_model_snapshot_plus_anomalies

      character(len=1024)                              :: filename_climate_anomalies

      real(dp), dimension(:,:), contiguous, pointer :: T2m_baseline    => null()   !< [K]                      Baseline monthly mean 2-m air temperature
      real(dp), dimension(:,:), contiguous, pointer :: Precip_baseline => null()   !< [m.w.e. month^-1]        Baseline monthly total precipitation
      real(dp), dimension(:  ), contiguous, pointer :: Hs_baseline     => null()   !< [m w.r.t. PD sea level]  Baseline surface elevation
      type(MPI_WIN) :: wT2m_baseline, wPrecip_baseline, wHs_baseline

      ! Insolation
      real(dp), dimension(:,:), contiguous, pointer :: Q_TOA           => null()
      type(MPI_WIN) :: wQ_TOA

      real(dp), dimension(:,:), contiguous, pointer :: ins_Q_TOA0  => null()
      real(dp), dimension(:,:), contiguous, pointer :: ins_Q_TOA1  => null()
      type(MPI_WIN) :: wins_Q_TOA0, wins_Q_TOA1

      ! Wind
      real(dp), dimension(:,:), contiguous, pointer :: Wind_LR         => null()   ! [m s^-1]     10m wind velocity in x-direction
      real(dp), dimension(:,:), contiguous, pointer :: Wind_DU         => null()   ! [m s^-1]     10m wind velocity in y-direction
      type(MPI_WIN) :: wWind_LR, wWind_DU

      ! Region-specific info
      character(len=1024)                           :: filename_climate_snapshot
      logical                                       :: do_lapse_rate_corrections
      real(dp)                                      :: lapse_rate_temp
      real(dp)                                      :: lapse_rate_precip
      logical                                       :: has_insolation

      integer,                    allocatable       :: ins_nyears
      integer,                    allocatable       :: ins_nlat,ins_nlon
      real(dp), dimension(:    ), allocatable       :: ins_time
      real(dp), dimension(:    ), allocatable       :: ins_lat
      real(dp),                   allocatable       :: ins_t0, ins_t1
      integer,                    allocatable       :: ins_ti0,ins_ti1
      
      ! anomalies
      real(dp), dimension(:,:), contiguous, pointer :: T2m_anomaly       => null()  ! Monthly mean 2m air temperature anomaly (K)
      real(dp), dimension(:,:), contiguous, pointer :: Precip_anomaly    => null()  ! Monthly mean precipitation anomaly (m)
      type(MPI_WIN) :: wT2m_anomaly, wPrecip_anomaly

      real(dp), dimension(:,:), contiguous, pointer :: T2m_anomaly_0     => null()  ! Prev. timeframe of monthly mean 2m air temperature anomaly (K)
      real(dp), dimension(:,:), contiguous, pointer :: Precip_anomaly_0  => null()  ! Prev. timeframe of monthly mean precipitation anomaly (m)
      type(MPI_WIN) :: wT2m_anomaly_0, wPrecip_anomaly_0

      real(dp), dimension(:,:), contiguous, pointer :: T2m_anomaly_1     => null()   ! Next timeframe of monthly mean 2m air temperature anomaly (K)
      real(dp), dimension(:,:), contiguous, pointer :: Precip_anomaly_1  => null()   ! Next timeframe of monthly mean precipitation anomaly (m)
      type(MPI_WIN) :: wT2m_anomaly_1, wPrecip_anomaly_1
      REAL(dp)                                      :: anomaly_t0, anomaly_t1        ! times of each timeframe

      ! Reference absorbed insolation (for GCM snapshots), or insolation at model time for the applied climate
      real(dp), dimension(:,:),   allocatable       :: Albedo                        ! Monthly mean surface albedo (calculated using our own SMB scheme for consistency)
      real(dp), dimension(:  ),   allocatable       :: I_abs                         ! Total yearly absorbed insolation, used in the climate matrix for interpolation

      contains
        procedure, public :: allocate_climate_model   => climate_model_snapshot_plus_anomalies_allocate
        procedure, public :: deallocate_climate_model => climate_model_snapshot_plus_anomalies_deallocate
        procedure, public :: initialise_climate_model => climate_model_snapshot_plus_anomalies_initialise
        procedure, public :: run_climate_model        => climate_model_snapshot_plus_anomalies_run
        procedure, public :: remap_climate_model      => climate_model_snapshot_plus_anomalies_remap
        procedure, public :: get_climate_model_name

        procedure, private :: apply_geometry_downscaling_corrections

        procedure, private :: update_climate_timeframes
        procedure, private :: interpolate_between_climate_anomaly_timeframes

        procedure, private :: initialise_insolation_forcing
        procedure, private :: get_insolation_at_time
        procedure, private :: update_insolation_timeframes_from_file

    end type type_climate_model_snapshot_plus_anomalies

  


contains

! ===== Main routines =====
! =========================
  subroutine climate_model_snapshot_plus_anomalies_allocate( self)

    ! In/output variables:
    class(type_climate_model_snapshot_plus_anomalies), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_snapshot_plus_anomalies_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%create_field( self%T2m_baseline, self%wT2m_baseline, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m_baseline', &
      long_name = 'Baseline monthly mean 2-m air temperature', &
      units     = 'K', &
      remap_method = 'reallocate')

    call self%create_field( self%Precip_baseline, self%wPrecip_baseline, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Precip_baseline', &
      long_name = 'Baseline monthly total precipitation', &
      units     = 'm.w.e. month^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%Hs_baseline, self%wHs_baseline, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Hs_baseline', &
      long_name = 'Baseline Surface elevation', &
      units     = 'm', &
      remap_method = 'reallocate')

    call self%create_field( self%Q_TOA, self%wQ_TOA, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Q_TOA', &
      long_name = 'Monthly insolation', &
      units     = 'W m^-2', &
      remap_method = 'reallocate')

    call self%create_field( self%ins_Q_TOA0, self%wins_Q_TOA0, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'ins_Q_TOA0', &
      long_name = 'Monthly insolation (prev timeframe)', &
      units     = 'W m^-2', &
      remap_method = 'reallocate')

    call self%create_field( self%ins_Q_TOA1, self%wins_Q_TOA1, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'ins_Q_TOA1', &
      long_name = 'Monthly insolation (next timeframe)', &
      units     = 'W m^-2', &
      remap_method = 'reallocate')

    call self%create_field( self%Wind_LR, self%wWind_LR, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Wind_LR', &
      long_name = 'Monthly 10m wind velocity in x-direction', &
      units     = 'm s^-1', &
      remap_method = 'reallocate')

    call self%create_field( self%Wind_DU, self%wWind_DU, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Wind_DU', &
      long_name = 'Monthly 10m wind velocity in y-direction', &
      units     = 'm s^-1', &
      remap_method = 'reallocate')

    ! Anomalies
    call self%create_field( self%T2m_anomaly, self%wT2m_anomaly, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m_anomaly', &
      long_name = 'Monthly mean 2-m air temperature anomaly', &
      units     = 'K', &
      remap_method = 'reallocate')

    call self%create_field( self%Precip_anomaly, self%wPrecip_anomaly, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Precip_anomaly', &
      long_name = 'Monthly total precipitation anomaly (ratio)', &
      units     = '-', &
      remap_method = 'reallocate')

    call self%create_field( self%T2m_anomaly_0, self%wT2m_anomaly_0, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m_anomaly_0', &
      long_name = 'Monthly mean 2-m air temperature anomaly (prev frame)', &
      units     = 'K', &
      remap_method = 'reallocate')

    call self%create_field( self%Precip_anomaly_0, self%wPrecip_anomaly_0, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Precip_anomaly_0', &
      long_name = 'Monthly total precipitation anomaly (ratio, prev frame)', &
      units     = '-', &
      remap_method = 'reallocate')

    call self%create_field( self%T2m_anomaly_1, self%wT2m_anomaly_1, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m_anomaly_1', &
      long_name = 'Monthly mean 2-m air temperature anomaly (next frame)', &
      units     = 'K', &
      remap_method = 'reallocate')

    call self%create_field( self%Precip_anomaly_1, self%wPrecip_anomaly_1, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Precip_anomaly_1', &
      long_name = 'Monthly total precipitation anomaly (ratio, next frame)', &
      units     = '-', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_plus_anomalies_allocate

  subroutine climate_model_snapshot_plus_anomalies_deallocate( self)

    ! In/output variables:
    class(type_climate_model_snapshot_plus_anomalies), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_snapshot_plus_anomalies_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to climate model snapshot_plus_transient_deltaT
    nullify( self%T2m_baseline)
    nullify( self%Precip_baseline)
    nullify( self%Hs_baseline)

    ! Insolation
    nullify( self%Q_TOA)
    nullify( self%ins_Q_TOA0)
    nullify( self%ins_Q_TOA1)

    ! Winds
    nullify( self%Wind_LR)
    nullify( self%Wind_DU)

    ! Anomalies
    nullify( self%T2m_anomaly)
    nullify( self%Precip_anomaly)
    nullify( self%T2m_anomaly_0)
    nullify( self%Precip_anomaly_0)
    nullify( self%T2m_anomaly_1)
    nullify( self%Precip_anomaly_1)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_plus_anomalies_deallocate

  subroutine climate_model_snapshot_plus_anomalies_initialise( self, geom, refgeo_PD, refgeo_init)
    ! Initialise the climate model
    !
    ! In- and output variables
    class(type_climate_model_snapshot_plus_anomalies),        intent(inout) :: self
    class(atype_ice_geometry_model_data),                     intent(in   ) :: geom
    type(type_reference_geometry),                            intent(in   ) :: refgeo_PD
    type(type_reference_geometry),                            intent(in   ) :: refgeo_init

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'climate_model_snapshot_plus_anomalies_initialise'
    integer                                               :: vi, m
    character(len=256)                                    :: filename_climate_snapshot
    real(dp)                                              :: timeframe_init_insolation, w0, w1

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen realistic climate model
    self%has_insolation = .FALSE.

    ! Read single-time data from external file
    ! Determine which climate model to initialise for this region
    select case(self%region_name())
    case ('NAM')
        filename_climate_snapshot       = C%climate_snp_p_anml_filename_snapshot_NAM
        self%filename_climate_anomalies = C%climate_snp_p_anml_filename_anomalies_NAM
        self%do_lapse_rate_corrections  = C%do_lapse_rate_corrections_NAM
        self%lapse_rate_precip          = C%precip_CC_correction_NAM
        self%lapse_rate_temp            = C%lapse_rate_temp_NAM
        if (C%choice_SMB_model_NAM == 'IMAU-ITM' .or. C%choice_SMB_model_NAM == 'ITM_v2') then
           self%has_insolation = .TRUE.
        end if
    case ('EAS')
        filename_climate_snapshot       = C%climate_snp_p_anml_filename_snapshot_EAS
        self%filename_climate_anomalies = C%climate_snp_p_anml_filename_anomalies_EAS
        self%do_lapse_rate_corrections  = C%do_lapse_rate_corrections_EAS
        self%lapse_rate_precip          = C%precip_CC_correction_EAS
        self%lapse_rate_temp            = C%lapse_rate_temp_EAS
        if (C%choice_SMB_model_EAS == 'IMAU-ITM' .or. C%choice_SMB_model_EAS == 'ITM_v2') then
           self%has_insolation = .TRUE.
        end if
    case ('GRL')
        filename_climate_snapshot       = C%climate_snp_p_anml_filename_snapshot_GRL
        self%filename_climate_anomalies = C%climate_snp_p_anml_filename_anomalies_GRL
        self%do_lapse_rate_corrections  = C%do_lapse_rate_corrections_GRL
        self%lapse_rate_precip          = C%precip_CC_correction_GRL
        self%lapse_rate_temp            = C%lapse_rate_temp_GRL
        if (C%choice_SMB_model_GRL == 'IMAU-ITM' .or. C%choice_SMB_model_GRL == 'ITM_v2') then
           self%has_insolation = .TRUE.
        end if
    case ('ANT')
        filename_climate_snapshot       = C%climate_snp_p_anml_filename_snapshot_ANT
        self%filename_climate_anomalies = C%climate_snp_p_anml_filename_anomalies_ANT
        self%do_lapse_rate_corrections  = C%do_lapse_rate_corrections_ANT
        self%lapse_rate_precip          = C%precip_CC_correction_ANT
        self%lapse_rate_temp            = C%lapse_rate_temp_ANT
        if (C%choice_SMB_model_ANT == 'IMAU-ITM' .or. C%choice_SMB_model_ANT == 'ITM_v2') then
           self%has_insolation = .TRUE.
        end if
    case default
      CALL crash('unknown region_name "' // self%region_name() // '"')
    end select

    ! Read in baseline variables
    CALL read_field_from_file_2D(         filename_climate_snapshot, 'Hs'    , self%mesh, C%output_dir, self%Hs_baseline)
    CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'T2m'   , self%mesh, C%output_dir, self%T2m_baseline)
    CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'Precip', self%mesh, C%output_dir, self%Precip_baseline)
    call read_field_from_file_2D_monthly( filename_climate_snapshot, 'uas'   , self%mesh, C%output_dir, self%Wind_LR)
    call read_field_from_file_2D_monthly( filename_climate_snapshot, 'vas'   , self%mesh, C%output_dir, self%Wind_DU)

    ! Read in anomalies

    ! Initialise anomaly timeframes
    self%anomaly_t0 = C%start_time_of_run - 200._dp
    self%anomaly_t1 = C%start_time_of_run - 100._dp
    call self%update_climate_timeframes(C%start_time_of_run)

    ! Interpolate between the two timeframes to find the applied anomaly
    call self%interpolate_between_climate_anomaly_timeframes(C%start_time_of_run)

    ! correct for mismatch between GCM and ISM ice surfaces
    call self%apply_geometry_downscaling_corrections(geom)

    ! Initialises the insolation (if needed)
    if (self%has_insolation) then
        call self%initialise_insolation_forcing()
        call self%get_insolation_at_time(C%start_time_of_run)
    end if

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  end subroutine climate_model_snapshot_plus_anomalies_initialise

  subroutine climate_model_snapshot_plus_anomalies_run( self, geom, time)
    ! Calculate the climate
    !
    ! Use a snapshot plus anomalies (temp difference, precip ratio)

    ! In/output variables:
    class(type_climate_model_snapshot_plus_anomalies), intent(inout) :: self
    class(atype_ice_geometry_model_data),              intent(in   ) :: geom
    real(dp),                                          intent(in   ) :: time

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'climate_model_snapshot_plus_anomalies_run'

    ! Add routine to path
    call init_routine( routine_name)

    ! Run the chosen snapshot+deltaT(time) climate model
    ! If the current model time falls outside the enveloping window
    ! of the two timeframes that have been read, update them
    if (time < self%anomaly_t0 .or. time > self%anomaly_t1) then
      call self%update_climate_timeframes(time)
    end if

    ! Interpolate between the two timeframes to find the applied anomaly
    call self%interpolate_between_climate_anomaly_timeframes(time)

    ! Update temperature and precipitation fields based on the mismatch between
    ! the ice sheet surface elevation in the forcing climate and the model's ice sheet surface elevation
    call self%apply_geometry_downscaling_corrections(geom)

    ! if needed for IMAU-ITM or climate matrix, we need to update insolation
    if (self%has_insolation) then
      call self%get_insolation_at_time(time)
    end if

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  end subroutine climate_model_snapshot_plus_anomalies_run

  subroutine climate_model_snapshot_plus_anomalies_remap(self, mesh_new)
    ! In/out variables
    class(type_climate_model_snapshot_plus_anomalies), intent(inout) :: self
    type(type_mesh), target,                           intent(in   ) :: mesh_new
    

    ! Local variables
    character(LEN=256), parameter                         :: routine_name = 'climate_model_snapshot_plus_anomalies_remap'
    
    ! Add routine to path
    call init_routine( routine_name)

    call self%remap_field( mesh_new, 'T2m_baseline',    self%T2m_baseline)
    call self%remap_field( mesh_new, 'Precip_baseline', self%Precip_baseline)
    call self%remap_field( mesh_new, 'Hs_baseline',     self%Hs_baseline)

    call self%remap_field( mesh_new, 'Q_TOA' ,          self%Q_TOA)
    call self%remap_field( mesh_new, 'ins_Q_TOA0',      self%ins_Q_TOA0)
    call self%remap_field( mesh_new, 'ins_Q_TOA1',      self%ins_Q_TOA1)

    call self%remap_field( mesh_new, 'Wind_LR',      self%Wind_LR)
    call self%remap_field( mesh_new, 'Wind_DU',      self%Wind_DU)

    call self%remap_field( mesh_new, 'T2m_anomaly',    self%T2m_anomaly)
    call self%remap_field( mesh_new, 'Precip_anomaly', self%Precip_anomaly)

    call self%remap_field( mesh_new, 'T2m_anomaly_0',    self%T2m_anomaly_0)
    call self%remap_field( mesh_new, 'Precip_anomaly_0', self%Precip_anomaly_0)
    call self%remap_field( mesh_new, 'T2m_anomaly_1',    self%T2m_anomaly_1)
    call self%remap_field( mesh_new, 'Precip_anomaly_1', self%Precip_anomaly_1)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_plus_anomalies_remap

  subroutine update_climate_timeframes(self, time)

    ! In/output variables:
    class(type_climate_model_snapshot_plus_anomalies), intent(inout) :: self
    real(dp),                                          intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'update_climate_timeframes'
    integer                             :: ncid, id_dim_time, nt, id_var_time, ierr
    real(dp), dimension(:), allocatable :: time_from_file
    integer                             :: ti0, ti1

    ! Add routine to path
    call init_routine( routine_name)

    ! Read time variable from the file
    call open_existing_netcdf_file_for_reading( self%filename_climate_anomalies, ncid)
    call check_time( self%filename_climate_anomalies, ncid)
    call inquire_dim_multopt( self%filename_climate_anomalies, ncid, field_name_options_time, id_dim_time, dim_length = nt)
    call inquire_var_multopt( self%filename_climate_anomalies, ncid, field_name_options_time, id_var_time)
    allocate( time_from_file( nt))
    call read_var_primary( self%filename_climate_anomalies, ncid, id_var_time, time_from_file)
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
    call read_field_from_file_2D_monthly( self%filename_climate_anomalies, 'T2m_anomaly', &
      self%mesh, C%output_dir, self%T2m_anomaly_0, &
      time_to_read = self%anomaly_t0)
    call read_field_from_file_2D_monthly( self%filename_climate_anomalies, 'T2m_anomaly', &
      self%mesh, C%output_dir, self%T2m_anomaly_1, &
      time_to_read = self%anomaly_t1)

    call read_field_from_file_2D_monthly( self%filename_climate_anomalies, 'Precip_anomaly', &
      self%mesh, C%output_dir, self%Precip_anomaly_0, &
      time_to_read = self%anomaly_t0)
    call read_field_from_file_2D_monthly( self%filename_climate_anomalies, 'Precip_anomaly', &
      self%mesh, C%output_dir, self%Precip_anomaly_1, &
      time_to_read = self%anomaly_t1)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_climate_timeframes

  subroutine interpolate_between_climate_anomaly_timeframes(self, time)
    ! Calculates anomalies by interpolating between the two
    ! timeframes to find the applied anomaly, then applies it to the climate variables


    ! In/output variables:
    class(type_climate_model_snapshot_plus_anomalies), intent(inout) :: self
    real(dp),                                          intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'interpolate_between_climate_anomaly_timeframes'
    real(dp)                            :: w0, w1
    integer                             :: vi, m

    ! Add routine to path
    call init_routine( routine_name)

    ! We need to make sure that the weight interpolation does not produce absurd anomaly values
    ! capping when the simulation time happens before/after all time frames
    if (time < self%anomaly_t0 .and. time < self%anomaly_t1) then
      call warning('Simulation time is before any available anomaly timeframe. Using first timeframe...')
      w0 = 1._dp
    elseif (time > self%anomaly_t0 .and. time > self%anomaly_t1) then
      call warning('Simulation time is after any available anomaly timeframe. Using last timeframe...')
      w0 = 0._dp
    else
      w0 = (self%anomaly_t1 - time) / (self%anomaly_t1 - self%anomaly_t0)
    end if

    w1 = 1._dp - w0

    do m = 1,12
      do vi = self%mesh%vi1, self%mesh%vi2

        ! Note that the baseline and the applied temperature and precip are monthly, but the anomaly is annual
        self%T2m_anomaly( vi,m)    = w0 * self%T2m_anomaly_0( vi,m)    + w1 * self%T2m_anomaly_1( vi,m)
        self%Precip_anomaly( vi,m) = w0 * self%Precip_anomaly_0( vi,m) + w1 * self%Precip_anomaly_1( vi,m)

        ! Add/multiply anomaly to T2m/Precip snapshots 
        ! to find the applied temperature and precipitation
        self%T2m( vi,m)     = self%T2m_baseline( vi,m)     + self%T2m_anomaly( vi,m)
        self%Precip( vi,m)  = self%Precip_baseline( vi,m)  * self%Precip_anomaly( vi,m)
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine interpolate_between_climate_anomaly_timeframes

  subroutine apply_geometry_downscaling_corrections(self, geom)
    ! Applies corrections for temperature and precipitation
    ! to correct for the mismatch between T and P at the forcing's ice surface elevation and the model's ice surface elevation
    
    ! In this version, since T2m and Precip have been computed with their anomalies before
    ! the downscaling function is called (and are always updated prior to its call),
    ! we use T2m and Precip directly instead of their baseline

    class(type_climate_model_snapshot_plus_anomalies), intent(inout) :: self
    class(atype_ice_geometry_model_data),              intent(in   ) :: geom

    ! local variables
    character(len=256), parameter                        :: routine_name = 'apply_geometry_downscaling_corrections'
    integer                                              :: vi, m
    real(dp)                                             :: deltaH, deltaT, deltaP
    real(dp), dimension(:,:), allocatable                :: T_inv, T_inv_ref

    ! Add routine to path
    CALL init_routine( routine_name)

    
    allocate( T_inv     (self%mesh%vi1:self%mesh%vi2, 12))
    allocate( T_inv_ref (self%mesh%vi1:self%mesh%vi2, 12))

    do vi = self%mesh%vi1, self%mesh%vi2

    ! we only apply corrections where it is not open ocean
    if (geom%mask_icefree_ocean( vi) .eqv. .FALSE.) then
        deltaT  = (geom%Hs( vi) - self%Hs_baseline( vi)) * (-1._dp * abs(self%lapse_rate_temp))
        do m = 1, 12
        ! Do corrections - based on Eq. 11 of Albrecht et al. (2020; TC) for PISM
        self%T2m( vi, m) = self%T2m( vi, m) + deltaT

        ! Calculate inversion-layer temperatures
        T_inv_ref( vi, m) = 88.9_dp + 0.67_dp *  self%T2m( vi, m)
        T_inv(     vi, m) = 88.9_dp + 0.67_dp * (self%T2m( vi, m) - self%lapse_rate_temp * (geom%Hs( vi) - self%Hs_baseline( vi)))
        
        ! Correct precipitation based on a simple Clausius-Clapeyron method (Jouzel & Merlivat, 1984; Huybrechts, 2002)
        ! Same as implemented in IMAU-ICE
        self%Precip( vi, m) = self%Precip( vi, m) * (T_inv_ref( vi, m) / T_inv( vi, m))**2 * EXP(22.47_dp * (T0 / T_inv_ref( vi, m) - T0 / T_inv( vi, m)))
        end do ! m
    end if
    end do ! vi

    deallocate(T_inv)
    deallocate(T_inv_ref)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_geometry_downscaling_corrections


  ! ==================================================
  ! Insolation
  ! ==================================================
  subroutine initialise_insolation_forcing(self)
    ! initialise the insolation series in the forcing structure

    ! In/output variables:
    class(type_climate_model_snapshot_plus_anomalies),     intent(inout) :: self

    ! Local variables:
    character(len=256), parameter                    :: routine_name = 'initialise_insolation_forcing'
    character(len=256)                               :: str
    integer                                          :: ncid
    real(dp)                                         :: closest_t0

    ! Add routine to path
    call init_routine( routine_name)

    if     (C%choice_insolation_forcing == 'none') then
      ! No insolation included, likely because we're running an idealised-geometry experiment
    elseif (C%choice_insolation_forcing == 'static' .or. &
            C%choice_insolation_forcing == 'realistic') then

      ! Initialise insolation
      ! The times at which we have insolation fields from Laskar, between which we'll interpolate
      ! to find the insolation at model time (assuming that t0 <= model_time <= t1)
      if (par%primary)   write(0,*) ' initialising insolation data from ', trim(C%filename_insolation), '...'

      ! Memory allocation
      allocate(self%ins_t0)
      allocate(self%ins_t1)
      allocate(self%ins_ti0)
      allocate(self%ins_ti1)
      allocate(self%ins_nlat)
      allocate(self%ins_nlon)
      allocate(self%Albedo    (self%mesh%vi1:self%mesh%vi2,12)) ! TODO: do we need to allocate this here?
      allocate(self%I_abs     (self%mesh%vi1:self%mesh%vi2))    ! TODO: do we need to allocate this here?
      allocate(self%ins_lat   (self%ins_nlat))
      self%ins_t0     = C%start_time_of_run
      self%ins_t1     = C%start_time_of_run
      self%ins_nlat   = 181
      self%ins_nlon   = 360
      self%ins_lat    = 0._dp
      self%ins_Q_TOA0 = 0._dp
      self%ins_Q_TOA1 = 0._dp
      self%Q_TOA      = 0._dp

      ! find the closest timeframe to the start of the run
      call read_field_from_file_0D( C%filename_insolation, field_name_options_time, closest_t0, time_to_read = C%start_time_of_run)

      if (C%start_time_of_run >= closest_t0) then
        if (par%primary) WRITE(0,*) '     start time is after ins_t0, reading one step further for ins_t1...'
        self%ins_t0 = closest_t0
        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, self%ins_t1, time_to_read = C%start_time_of_run+1000._dp)
      else
        ! otherwise we read one record before for t1
        if (par%primary) WRITE(0,*) '     start time is before closest ins_t0, reading one step earlier for t0, and using that one for t1...'
        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, self%ins_t0, time_to_read = C%start_time_of_run-1000._dp)
        self%ins_t1 = closest_t0
      end if

      ! Read the fields at ins_t0
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, self%mesh, C%output_dir, self%ins_Q_TOA0, time_to_read = self%ins_t0)
      ! Read the fields at ins_t1
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, self%mesh, C%output_dir, self%ins_Q_TOA1, time_to_read = self%ins_t1)

    else
      call crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_insolation_forcing

  subroutine get_insolation_at_time(self, time)
    ! Get monthly insolation at time t on the regional grid

    ! In/output variables:
    class(type_climate_model_snapshot_plus_anomalies), intent(inout) :: self
    real(dp),                                               intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                    :: routine_name = 'get_insolation_at_time'
    real(dp)                                         :: time_applied
    integer                                          :: vi,m
    real(dp)                                         :: wt0, wt1

    ! Add routine to path
    call init_routine( routine_name)

    time_applied = 0._dp

    ! Safety
    if     (C%choice_insolation_forcing == 'none') then
      call crash('insolation should not be used when choice_insolation_forcing = "none"!')
    elseif (C%choice_insolation_forcing == 'static') then
      time_applied = C%static_insolation_time
    elseif (C%choice_insolation_forcing == 'realistic') then
      time_applied = time
    else
      call crash('unknown choice_insolation_forcing "' // trim( C%choice_insolation_forcing) // '"!')
    end if

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    if (time_applied < self%ins_t0 .or. time_applied > self%ins_t1) then
      if (par%primary)  write(0,*) '   Model time is out of the current insolation timeframes. Updating timeframes...'
      call self%update_insolation_timeframes_from_file(time_applied)
    end if

    ! Calculate timeframe interpolation weights (plus safety checks for when the extend beyond the record)
    if (self%ins_t1 == self%ins_t0) then
      wt0 = 0._dp
      wt1 = 1._dp
    else
      if (time_applied > self%ins_t1) then
        wt0 = 0._dp
      elseif (time_applied < self%ins_t0) then
        wt0 = 1._dp
      else
        wt0 = (self%ins_t1 - time_applied) / (self%ins_t1 - self%ins_t0)
      end if
      wt1 = 1._dp - wt0
    end if

    ! Interpolate the two timeframes
    do vi = self%mesh%vi1, self%mesh%vi2
      do m = 1, 12
        self%Q_TOA(vi, m) = wt0 * self%ins_Q_TOA0(vi, m) + wt1 * self%ins_Q_TOA1(vi, m)
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine get_insolation_at_time

  subroutine update_insolation_timeframes_from_file(self,time)
    ! Read the NetCDF file containing the insolation forcing data. Only read the time frames enveloping the current
    ! coupling timestep to save on memory usage. Only done by master.

    ! NOTE: assumes time in forcing file is in kyr
    class(type_climate_model_snapshot_plus_anomalies), intent(inout) :: self
    real(dp),                                               intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                      :: routine_name = 'update_insolation_timeframes_from_file'
    integer                                            :: ti0, ti1, ncid
    character(len=256)                                 :: str

    ! Add routine to path
    call init_routine( routine_name)

    if     (C%choice_insolation_forcing == 'none') then
      call crash('insolation should not be used when choice_insolation_forcing = "none"!')
    elseif (C%choice_insolation_forcing == 'static' .or. &
            C%choice_insolation_forcing == 'realistic') then

      ! Update insolation
      call read_field_from_file_0D( C%filename_insolation, field_name_options_time, self%ins_t0, time_to_read = time)

      if (time >= self%ins_t0) then
        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, self%ins_t1, time_to_read = time+1000._dp)
      else
        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, self%ins_t1, time_to_read = time)
        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, self%ins_t0, time_to_read = time-1000._dp)
      end if

      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, self%mesh, C%output_dir, self%ins_Q_TOA0, time_to_read = self%ins_t0)
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, self%mesh, C%output_dir, self%ins_Q_TOA1, time_to_read = self%ins_t1)

      !call warning('insolation timeframes at t = {dp_01} are ins_t0={dp_02} and ins_t1={dp_03}', dp_01 =  time, dp_02 = self%ins_t0, dp_03 = self%%ins_t1)

    else
      call crash('unknown choice_insolation_forcing "' // trim( C%choice_insolation_forcing) // '"!')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_insolation_timeframes_from_file

  function get_climate_model_name( self) result( climate_model_name)
    class(type_climate_model_snapshot_plus_anomalies), intent(in) :: self
    character(len=:), allocatable :: climate_model_name
    climate_model_name = 'snapshot_plus_anomalies'
  end function get_climate_model_name

end module climate_snapshot_plus_anomalies
