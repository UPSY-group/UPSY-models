module climate_snapshot_plus_uniform_deltaT

  ! Snapshot + dT climate model

! ===== Preamble =====
! ====================

  use precisions,                        only: dp
  use upsy_main,                         only: upsy
  use mpi_basic,                         only: par, sync
  use call_stack_and_comp_time_tracking, only: crash, init_routine, finalise_routine, warning
  use model_configuration,               only: C
  use parameters
  use mpi_f08,                           only: MPI_WIN
  use mesh_basic,                        only: type_mesh
  use ice_geometry_model_data,           only: atype_ice_geometry_model_data
  use reference_geometry_types,          only: type_reference_geometry
  !   use climate_model_types,               only: type_climate_model, type_climate_model_snapshot
  use Arakawa_grid_mod,                  only: Arakawa_grid
  use fields_dimensions,                 only: third_dimension
  use global_forcings_main
  use netcdf_io_main
  use netcdf_basic
  !   use climate_model_utilities,           only: get_insolation_at_time
  use reallocate_mod,                    only: reallocate_bounds
  use climate_model_basic,               only: atype_climate_model

  implicit none

  private
  public :: type_climate_model_snapshot_plus_uniform_deltaT

    type, extends(atype_climate_model) :: type_climate_model_snapshot_plus_uniform_deltaT

        real(dp), dimension(:,:), contiguous, pointer :: T2m_baseline    => null()   !< [K]                      Baseline monthly mean 2-m air temperature
        real(dp), dimension(:,:), contiguous, pointer :: Precip_baseline => null()   !< [m.w.e. month^-1]        Baseline monthly total precipitation
        real(dp), dimension(:  ), contiguous, pointer :: Hs_baseline     => null()   !< [m w.r.t. PD sea level]  Baseline surface elevation
        type(MPI_WIN) :: wT2m_baseline, wPrecip_baseline, wHs_baseline

        real(dp), dimension(:,:), contiguous, pointer :: ins_Q_TOA0  => null()
        real(dp), dimension(:,:), contiguous, pointer :: ins_Q_TOA1  => null()
        type(MPI_WIN) :: wins_Q_TOA0, wins_Q_TOA1
        
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
        real(dp)                                      :: deltaT

        ! Reference absorbed insolation (for GCM snapshots), or insolation at model time for the applied climate
        real(dp), dimension(:,:),   allocatable       :: Albedo                        ! Monthly mean surface albedo (calculated using our own SMB scheme for consistency)
        real(dp), dimension(:  ),   allocatable       :: I_abs                         ! Total yearly absorbed insolation, used in the climate matrix for interpolation
        
        contains

          procedure, public :: allocate_climate_model   => climate_model_snapshot_plus_uniform_deltaT_allocate
          procedure, public :: deallocate_climate_model => climate_model_snapshot_plus_uniform_deltaT_deallocate
          procedure, public :: initialise_climate_model => climate_model_snapshot_plus_uniform_deltaT_initialise
          procedure, public :: run_climate_model        => climate_model_snapshot_plus_uniform_deltaT_run
          procedure, public :: remap_climate_model      => climate_model_snapshot_plus_uniform_deltaT_remap
          procedure, public :: get_climate_model_name

          procedure, private :: apply_precipitation_CC_correction
          procedure, private :: apply_geometry_downscaling_corrections

          procedure, private :: initialise_insolation_forcing
          procedure, private :: get_insolation_at_time
          procedure, private :: update_insolation_timeframes_from_file

    end type type_climate_model_snapshot_plus_uniform_deltaT

contains

! ===== Main routines =====
! =========================
  subroutine climate_model_snapshot_plus_uniform_deltaT_allocate( self)

    ! In/output variables:
    class(type_climate_model_snapshot_plus_uniform_deltaT), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_snapshot_plus_uniform_deltaT_allocate'

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


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_plus_uniform_deltaT_allocate

  subroutine climate_model_snapshot_plus_uniform_deltaT_deallocate( self)

    ! In/output variables:
    class(type_climate_model_snapshot_plus_uniform_deltaT), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_snapshot_plus_uniform_deltaT_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to climate model snapshot_plus_uniform_deltaT
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

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_plus_uniform_deltaT_deallocate

  subroutine climate_model_snapshot_plus_uniform_deltaT_initialise(self, geom, refgeo_PD, refgeo_init)
    ! Initialise the climate model
    !
    ! Use a realistic climate scheme
    ! In- and output variables
    class(type_climate_model_snapshot_plus_uniform_deltaT), intent(inout) :: self
    class(atype_ice_geometry_model_data),                   intent(in   ) :: geom
    type(type_reference_geometry),                          intent(in   ) :: refgeo_PD
    type(type_reference_geometry),                          intent(in   ) :: refgeo_init

    ! local variables:
    character(len=256), parameter                         :: routine_name = 'climate_model_snapshot_plus_uniform_deltaT_initialise'
    integer                                               :: vi, m
    character(len=256)                                    :: filename_climate_snapshot
    real(dp)                                              :: timeframe_init_insolation

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%primary)  write(*,"(a)") '     initialising snapshot+uniform_dT climate model...'

    ! Run the chosen realistic climate model
    self%has_insolation = .FALSE.

    ! Determine which region to initialise
    select case (self%region_name())
      case('NAM')
          filename_climate_snapshot      = C%filename_climate_snapshot_unif_dT_NAM
          self%do_lapse_rate_corrections = C%do_lapse_rate_corrections_NAM
          self%lapse_rate_precip         = C%lapse_rate_precip_NAM
          self%lapse_rate_temp           = C%lapse_rate_temp_NAM
          self%deltaT                    = C%uniform_deltaT_NAM
          if (C%choice_SMB_model_NAM == 'IMAU-ITM' .or. C%choice_SMB_model_NAM == 'ITM_v2') then
           self%has_insolation = .TRUE.
        end if
      case('EAS')
          filename_climate_snapshot      = C%filename_climate_snapshot_unif_dT_EAS
          self%do_lapse_rate_corrections = C%do_lapse_rate_corrections_EAS
          self%lapse_rate_precip         = C%lapse_rate_precip_EAS
          self%lapse_rate_temp           = C%lapse_rate_temp_EAS
          self%deltaT                    = C%uniform_deltaT_EAS
          if (C%choice_SMB_model_EAS == 'IMAU-ITM' .or. C%choice_SMB_model_EAS == 'ITM_v2') then
           self%has_insolation = .TRUE.
        end if
      case('GRL')
          filename_climate_snapshot      = C%filename_climate_snapshot_unif_dT_GRL
          self%do_lapse_rate_corrections = C%do_lapse_rate_corrections_GRL
          self%lapse_rate_precip         = C%lapse_rate_precip_GRL
          self%lapse_rate_temp           = C%lapse_rate_temp_GRL
          self%deltaT                    = C%uniform_deltaT_GRL
          if (C%choice_SMB_model_GRL == 'IMAU-ITM' .or. C%choice_SMB_model_GRL == 'ITM_v2') then
           self%has_insolation = .TRUE.
        end if
      case('ANT')
          filename_climate_snapshot      = C%filename_climate_snapshot_unif_dT_ANT
          self%do_lapse_rate_corrections = C%do_lapse_rate_corrections_ANT
          self%lapse_rate_precip         = C%lapse_rate_precip_ANT
          self%lapse_rate_temp           = C%lapse_rate_temp_ANT
          self%deltaT                    = C%uniform_deltaT_ANT
          if (C%choice_SMB_model_ANT == 'IMAU-ITM' .or. C%choice_SMB_model_ANT == 'ITM_v2') then
           self%has_insolation = .TRUE.
        end if
      case default
        call crash('unknown region_name "' // self%region_name() // '"')
    end select

    call read_field_from_file_2D(         filename_climate_snapshot, 'Hs'    , self%mesh, C%output_dir, self%Hs_baseline)
    call read_field_from_file_2D_monthly( filename_climate_snapshot, 'T2m'   , self%mesh, C%output_dir, self%T2m_baseline)
    call read_field_from_file_2d_monthly( filename_climate_snapshot, 'Precip', self%mesh, C%output_dir, self%Precip_baseline)

    ! Adding deltaT to the baseline temperature field (uniform in space and time)
    do vi = self%mesh%vi1, self%mesh%vi2
    do m = 1, 12
        self%T2m_baseline( vi, m)   = self%T2m_baseline( vi, m) + self%deltaT
    end do
    end do

    ! Scale Baseline Precipitation based on the CC lapse-rate for the applied deltaT
    call self%apply_precipitation_CC_correction()

    ! apply corrections due to the mismatch between GCM ice surface and ISM ice surface
    call self%apply_geometry_downscaling_corrections(geom)

    ! Initialises the insolation (if needed)
    if (self%has_insolation) then
        call self%initialise_insolation_forcing()
        if (C%start_time_of_run < 0._dp) then
            timeframe_init_insolation = C%start_time_of_run
        else
            timeframe_init_insolation = 0._dp
        end if
        call self%get_insolation_at_time(timeframe_init_insolation)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_plus_uniform_deltaT_initialise

  subroutine climate_model_snapshot_plus_uniform_deltaT_run(self, geom, time)
    ! Calculate the climate
    !
    ! Use a snapshot plus a prescribed uniform deltaT

    ! In/output variables:
    class(type_climate_model_snapshot_plus_uniform_deltaT), intent(inout) :: self
    class(atype_ice_geometry_model_data),                   intent(in   ) :: geom
    real(dp),                                               intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'climate_model_snapshot_plus_uniform_deltaT_run'

    ! Add routine to path
    call init_routine( routine_name)

    
    ! Update temperature and precipitation fields based on the mismatch between
    ! the ice sheet surface elevation in the forcing climate and the model's ice sheet surface elevation
    call self%apply_geometry_downscaling_corrections(geom)

    ! if needed for an ITM or climate matrix, we need to update insolation
    if (self%has_insolation) then
      call self%get_insolation_at_time(time)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_plus_uniform_deltaT_run

  
  subroutine climate_model_snapshot_plus_uniform_deltaT_remap(self, mesh_new, refgeo_init, refgeo_PD, geom)
    ! In/out variables
    class(type_climate_model_snapshot_plus_uniform_deltaT), intent(inout) :: self
    type(type_mesh), target,                                intent(in   ) :: mesh_new
    type(type_reference_geometry), target,                  intent(in   ) :: refgeo_init, refgeo_PD
    class(atype_ice_geometry_model_data),                   intent(in   ) :: geom

    ! Local variables
    character(LEN=256), parameter                         :: routine_name = 'climate_model_snapshot_plus_uniform_deltaT_remap'

    ! Add routine to path
    call init_routine( routine_name)

    ! TODO: what is the proper way to remap these variables?
    call self%remap_field( mesh_new, 'T2m_baseline',    self%T2m_baseline)
    call self%remap_field( mesh_new, 'Precip_baseline', self%Precip_baseline)
    call self%remap_field( mesh_new, 'Hs_baseline',     self%Hs_baseline)

    call self%remap_field( mesh_new, 'Q_TOA' ,          self%Q_TOA)
    call self%remap_field( mesh_new, 'ins_Q_TOA0',      self%ins_Q_TOA0)
    call self%remap_field( mesh_new, 'ins_Q_TOA1',      self%ins_Q_TOA1)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_plus_uniform_deltaT_remap


  ! ==================================================
  ! Temperature and precipitation corrections
  ! ==================================================
  subroutine apply_precipitation_CC_correction(self)
  
    ! Applies a simple Clausius-Clapeyron correction to temperature based on the prescribed deltaT
    class(type_climate_model_snapshot_plus_uniform_deltaT), intent(inout) :: self

    ! Local Variables
    character(len=256), parameter                        :: routine_name = 'apply_precipitation_CC_correction'
    integer                                              :: vi, m

    ! Add routine to path
    call init_routine( routine_name)

    do vi = self%mesh%vi1, self%mesh%vi2
        do m = 1, 12
            ! Precip(\Delta T) = Precip(PD) \times 1.068^{\Delta T}
            self%Precip_baseline( vi, m) = self%Precip_baseline( vi, m) * self%lapse_rate_precip**self%deltaT
        end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_precipitation_CC_correction

  subroutine apply_geometry_downscaling_corrections(self, geom)
    ! Applies corrections for temperature and precipitation
    ! to correct for the mismatch between T and P at the forcing's ice surface elevation and the model's ice surface elevation

    class(type_climate_model_snapshot_plus_uniform_deltaT), intent(inout) :: self
    class(atype_ice_geometry_model_data),                   intent(in   )  :: geom

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
        self%T2m( vi, m) = self%T2m_baseline( vi, m) + deltaT

        ! Calculate inversion-layer temperatures
        T_inv_ref( vi, m) = 88.9_dp + 0.67_dp *  self%T2m( vi, m)
        T_inv(     vi, m) = 88.9_dp + 0.67_dp * (self%T2m( vi, m) - self%lapse_rate_temp * (geom%Hs( vi) - self%Hs_baseline( vi)))
        
        ! Correct precipitation based on a simple Clausius-Clapeyron method (Jouzel & Merlivat, 1984; Huybrechts, 2002)
        ! Same as implemented in IMAU-ICE
        self%Precip( vi, m) = self%Precip_baseline( vi, m) * (T_inv_ref( vi, m) / T_inv( vi, m))**2 * EXP(22.47_dp * (T0 / T_inv_ref( vi, m) - T0 / T_inv( vi, m)))
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
    class(type_climate_model_snapshot_plus_uniform_deltaT),     intent(inout) :: self

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
    class(type_climate_model_snapshot_plus_uniform_deltaT), intent(inout) :: self
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
    class(type_climate_model_snapshot_plus_uniform_deltaT), intent(inout) :: self
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
    class(type_climate_model_snapshot_plus_uniform_deltaT), intent(in) :: self
    character(len=:), allocatable :: climate_model_name
    climate_model_name = 'snapshot_plus_uniform_deltaT'
  end function get_climate_model_name

end module climate_snapshot_plus_uniform_deltaT
