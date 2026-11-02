module climate_snapshot

  ! "Snapshot of a climate" model

! ===== Preamble =====
! ====================

  use precisions,                        only: dp
  use upsy_main,                         only: upsy
  use mpi_basic,                         only: par, sync
  use call_stack_and_comp_time_tracking, only: crash, init_routine, finalise_routine, warning
  use model_configuration,               only: C
  use parameters
  use mesh_basic,                        only: type_mesh
  use ice_geometry_model_data,           only: atype_ice_geometry_model_data
  use climate_model_types,               only: type_climate_model, type_climate_model_snapshot
  use global_forcing_types,              only: type_global_forcing
  use global_forcings_main
  use netcdf_io_main
  use netcdf_basic
  use climate_model_utilities,           only: get_insolation_at_time
  use reallocate_mod,                    only: reallocate_bounds

  implicit none

  private
  public :: type_climate_model_snapshot

    type, extends(atype_climate_model) :: type_climate_model_snapshot

        real(dp), dimension(:  ), contiguous, pointer :: Hs     => null()   !< [m w.r.t. PD sea level]  Baseline surface elevation
        type(MPI_WIN) :: wHs

        real(dp), dimension(:,:), contiguous, pointer :: Q_TOA  => null()
        type(MPI_WIN) :: wQ_TOA
        
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
        real(dp), dimension(:,:  ), allocatable       :: ins_Q_TOA0, ins_Q_TOA1

        ! Reference absorbed insolation (for GCM snapshots), or insolation at model time for the applied climate
        real(dp), dimension(:,:),   allocatable       :: Q_TOA                         ! Monthly mean insolation at the top of the atmosphere (W/m2) (taken from the prescribed insolation solution at orbit_time)
        real(dp), dimension(:,:),   allocatable       :: Albedo                        ! Monthly mean surface albedo (calculated using our own SMB scheme for consistency)
        real(dp), dimension(:  ),   allocatable       :: I_abs                         ! Total yearly absorbed insolation, used in the climate matrix for interpolation
        
        contains

            procedure, public :: allocate_climate_model   => climate_model_snapshot_allocate
            procedure, public :: deallocate_climate_model => climate_model_snapshot_deallocate
            procedure, public :: initialise_climate_model => climate_model_snapshot_initialise
            procedure, public :: run_climate_model        => climate_model_snapshot_run
            procedure, public :: remap_climate_model      => climate_model_snapshot_remap

            procedure, public :: get_climate_model_name

            procedure, private :: apply_geometry_downscaling_corrections
            procedure, private :: initialise_insolation_forcing

    end type type_climate_model_snapshot
  

contains

! ===== Main routines =====
! =========================

subroutine climate_model_snapshot_allocate( self)

    ! In/output variables:
    class(type_climate_model_snapshot), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_snapshot_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is specific to the prescribed SMB model
    call self%create_field( self%Hs, self%wHs, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Hs', &
      long_name = 'Surface elevation', &
      units     = 'm', &
      remap_method = 'reallocate')

    call self%create_field( self%Q_TOA, self%wQ_TOA, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Q_TOA', &
      long_name = 'Monthly insolation', &
      units     = 'W m^-2', &
      remap_method = 'reallocate')


    ! Remove routine from call stack
    call finalise_routine( routine_name)

end subroutine climate_model_snapshot_allocate

  subroutine climate_model_snapshot_deallocate( self)

    ! In/output variables:
    class(type_climate_model_snapshot), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_snapshot_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to SMB model prescribed
    nullify( self%Hs)
    nullify( self%Q_TOA)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_deallocate


  subroutine climate_model_snapshot_initialise(self, geom)
    ! Initialise the climate model
    !
    ! Use a realistic climate scheme

    ! In- and output variables
    class(type_climate_model_snapshot),     intent(inout) :: self
    class(atype_ice_geometry_model_data),   intent(in   ) :: geom

    ! local variables:
    character(len=256), parameter                         :: routine_name = 'initialise_climate_model_realistic'
    character(len=256)                                    :: filename_climate_snapshot
    logical                                               :: do_lapse_rates
    real(dp)                                              :: timeframe_init_insolation

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    if (par%primary)  write(*,"(a)") '     initialising realistic climate model "' // &
      UPSY%stru%colour_string( trim( C%choice_climate_model_realistic),'light blue') // '"...'

    ! Run the chosen realistic climate model
    self%has_insolation = .FALSE.
    
    ! Determine which climate model to initialise for this region
    select case (self%region_name())
      case('NAM')
        filename_climate_snapshot      = C%filename_climate_snapshot_NAM
        self%do_lapse_rate_corrections = C%do_lapse_rate_corrections_NAM
        self%lapse_rate_precip         = C%lapse_rate_precip_NAM
        self%lapse_rate_temp           = C%lapse_rate_temp_NAM
        if (C%choice_SMB_model_NAM == 'IMAU-ITM' .or. C%choice_SMB_model_NAM == 'ITM_v2') then
           self%%has_insolation = .TRUE.
        end if
      case ('EAS')
        filename_climate_snapshot      = C%filename_climate_snapshot_EAS
        self%do_lapse_rate_corrections = C%do_lapse_rate_corrections_EAS
        self%lapse_rate_precip         = C%lapse_rate_precip_EAS
        self%lapse_rate_temp           = C%lapse_rate_temp_EAS
        if (C%choice_SMB_model_NAM == 'IMAU-ITM' .or. C%choice_SMB_model_NAM == 'ITM_v2') then
           self%has_insolation = .TRUE.
        end if
      case('GRL')
        filename_climate_snapshot      = C%filename_climate_snapshot_GRL
        self%do_lapse_rate_corrections = C%do_lapse_rate_corrections_GRL
        self%lapse_rate_precip         = C%lapse_rate_precip_GRL
        self%lapse_rate_temp           = C%lapse_rate_temp_GRL
        if (C%choice_SMB_model_NAM == 'IMAU-ITM' .or. C%choice_SMB_model_NAM == 'ITM_v2') then
           self%has_insolation = .TRUE.
        end if
      case('ANT') then
        filename_climate_snapshot      = C%filename_climate_snapshot_ANT
        self%do_lapse_rate_corrections = C%do_lapse_rate_corrections_ANT
        self%lapse_rate_precip         = C%lapse_rate_precip_ANT
        self%lapse_rate_temp           = C%lapse_rate_temp_ANT
        if (C%choice_SMB_model_NAM == 'IMAU-ITM' .or. C%choice_SMB_model_NAM == 'ITM_v2') then
           self%has_insolation = .TRUE.
        end if
      case default
        call crash('unknown region_name "' // self%region_name() // '"')
    end select

    call read_field_from_file_2D(         filename_climate_snapshot, 'Hs'    , self%mesh, C%output_dir, self%Hs)
    call read_field_from_file_2D_monthly( filename_climate_snapshot, 'T2m'   , self%mesh, C%output_dir, self%T2m)
    call read_field_from_file_2D_monthly( filename_climate_snapshot, 'Precip', self%mesh, C%output_dir, self%Precip)

    call self%apply_geometry_downscaling_corrections(geom)

    ! Initialises the insolation (if needed)
    if (self%has_insolation) then
        if (C%choice_insolation_forcing == 'none') then
            call crash('Chosen climate or SMB model cannot be used with choice_insolation_forcing = "none"!')
        else
            call self%initialise_insolation_forcing()
            if (C%start_time_of_run < 0._dp) then
                timeframe_init_insolation = C%start_time_of_run
            else
                timeframe_init_insolation = 0._dp
            end if
            call self%get_insolation_at_time(timeframe_init_insolation)
        end if
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_initialise

  subroutine climate_model_snapshot_run(self, geom, forcing, time)
    ! Calculate the climate
    !
    ! Use an realistic climate scheme

    ! In/output variables:
    class(type_climate_model_snapshot),     intent(inout) :: self
    class(atype_ice_geometry_model_data),   intent(in   ) :: geom
    type(type_global_forcing),              intent(in)    :: forcing
    real(dp),                               intent(in)    :: time

    ! local variables:
    character(len=256), parameter                         :: routine_name = 'climate_model_realistic_run'

    ! add routine to path
    call init_routine( routine_name)

    ! run the chosen realistic climate model

    ! update temperature and precipitation fields based on the mismatch between
    ! the ice sheet surface elevation in the forcing climate and the model's ice sheet surface elevation
    call self%apply_geometry_downscaling_corrections(geom)

    ! if needed for imau-itm or climate matrix, we need to update insolation
    if (self%has_insolation) then
      call self%get_insolation_at_time(self, time)
    end if

    ! finalise routine path
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_run

  subroutine apply_geometry_downscaling_corrections(self, geom)
    ! Applies the lapse rate corrections for temperature and precipitation
    ! to correct for the mismatch between T and P at the forcing's ice surface elevation and the model's ice surface elevation

    class(type_climate_model_snapshot),     intent(inout) :: self
    class(atype_ice_geometry_model_data),  intent(in   )  :: geom

    ! local variables
    character(len=256), parameter                        :: routine_name = 'apply_geometry_downscaling_corrections'
    integer                                              :: vi, m
    real(dp)                                             :: deltaH, deltaT, deltaP
    real(dp), dimension(:,:), allocatable                :: T_inv, T_inv_ref

    ! Add routine to path
    CALL init_routine( routine_name)

    
    allocate( T_inv     (self%mesh%vi1:self%mesh%vi2, 12))
    allocate( T_inv_ref (self%mesh%vi1:self%mesh%vi2, 12))

    do vi = mesh%vi1, mesh%vi2

    ! we only apply corrections where it is not open ocean
    if (geom%mask_icefree_ocean( vi) .eqv. .FALSE.) then
        deltaT  = (geom%Hs( vi) - self%Hs( vi)) * (-1._dp * abs(self%lapse_rate_temp))
        do m = 1, 12
        ! Do corrections - based on Eq. 11 of Albrecht et al. (2020; TC) for PISM
        self%T2m( vi, m) = self%T2m( vi, m) + deltaT

        ! Calculate inversion-layer temperatures
        T_inv_ref( vi, m) = 88.9_dp + 0.67_dp *  self%T2m( vi, m)
        T_inv(     vi, m) = 88.9_dp + 0.67_dp * (self%T2m( vi, m) - self%lapse_rate_temp * (geom%Hs( vi) - self%Hs( vi)))
        
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

  ! == Insolation

  subroutine climate_model_snapshot_remap(self, geom, time)
    ! In/out variables
    class(type_climate_model_snapshot),     intent(inout) :: self
    class(atype_ice_geometry_model_data),   intent(in   ) :: geom
    real(dp),                               intent(in)    :: time

    ! Local variables
    character(LEN=256), parameter                         :: routine_name = 'remap_climate_realistic'

    ! Add routine to path
    call init_routine( routine_name)

    ! TODO: what is the proper way to remap these variables?

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine climate_model_snapshot_remap

  subroutine initialise_insolation_forcing(self)
    ! initialise the insolation series in the forcing structure

    ! In/output variables:
    class(type_climate_model_snapshot),     intent(inout) :: self

    ! Local variables:
    character(len=256), parameter                    :: routine_name = 'initialise_insolation_forcing'
    character(len=256)                               :: str
    integer                                          :: ncid
    real(dp)                                         :: closest_t0

    ! Add routine to path
    CALL init_routine( routine_name)

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
      allocate(self%ins_Q_TOA0         (mesh%vi1:mesh%vi2,12))
      allocate(self%ins_Q_TOA1         (mesh%vi1:mesh%vi2,12))
      allocate(self%Albedo             (mesh%vi1:mesh%vi2,12))
      allocate(self%I_abs              (mesh%vi1:mesh%vi2))
      self%ins_t0     = C%start_time_of_run
      self%ins_t1     = C%start_time_of_run
      self%ins_nlat   = 181
      self%ins_nlon   = 360
      allocate(self%ins_lat            (   self%ins_nlat))
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
        snapshot%ins_t1 = closest_t0
      end if

      ! Read the fields at ins_t0
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, selfmesh, C%output_dir, self%ins_Q_TOA0, time_to_read = self%ins_t0)
      ! Read the fields at ins_t1
      call read_field_from_file_2D_monthly( C%filename_insolation, field_name_options_insolation, selfmesh, C%output_dir, self%ins_Q_TOA1, time_to_read = self%ins_t1)

    else
      call crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_insolation_forcing

  subroutine remap_insolation(self, geom)
    ! In/out variables
    class(type_climate_model_snapshot),     intent(inout) :: self
    class(atype_ice_geometry_model_data),   intent(in   ) :: geom

    ! Local variables
    character(LEN=256), parameter                         :: routine_name = 'remap_insolation'

    ! Add routine to path
    call init_routine( routine_name)

    ! TODO: is this the proper way to do it?
    call self%initialise_insolation_forcing()
    call self%get_insolation_at_time(time)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_insolation

end module climate_snapshot
