module ISMIP7_climate

  ! The ISMIP7 protocol provides separate NetCDF files for each year, with a single file contaning 12 monthly fields
  ! for that year. There are separate files for each variable, namely:
  !
  ! - tas        : the absolute 2m air temp (in K)
  ! - tas-anomaly: the 2m air temp anomaly relative to the 1960-1989 baseline
  ! - dtsdz      : the vertical surface temp gradient, to correct for differences in elevation between
  !                the observed and modelled ice-sheet geometry
  ! - pr         : the absolute precipitation (in kg m^-2 s^-1)
  ! - pr-anomaly : the precipitation anomaly relative to the 1960-1989 baseline
  !
  ! The directory structure for the forcing data is:
  !
  !  ../
  !    path/
  !      to/
  !        base_folder/
  !          tas/
  !            version/
  !              tas_somethingsomethingsomething_2015.nc
  !              tas_somethingsomethingsomething_2016.nc
  !              tas_somethingsomethingsomething_2017.nc
  !              ...
  !          tas-anomaly/
  !            version/
  !              tas-anomaly_somethingsomethingsomething_2015.nc
  !              tas-anomaly_somethingsomethingsomething_2016.nc
  !              tas-anomaly_somethingsomethingsomething_2017.nc
  !              ...
  !          dtsdz/
  !            version/
  !              dtsdz_somethingsomethingsomething_2015.nc
  !              dtsdz_somethingsomethingsomething_2016.nc
  !              dtsdz_somethingsomethingsomething_2017.nc
  !              ...
  !          pr/
  !            version/
  !              pr_somethingsomethingsomething_2015.nc
  !              pr_somethingsomethingsomething_2016.nc
  !              pr_somethingsomethingsomething_2017.nc
  !              ...
  !          pr-anomaly/
  !            version/
  !              pr-anomaly_somethingsomethingsomething_2015.nc
  !              pr-anomaly_somethingsomethingsomething_2016.nc
  !              pr-anomaly_somethingsomethingsomething_2017.nc
  !              ...
  !
  ! In the config, you only need to provide the path/to/base_folder and the version; UFEMISM will take
  ! care of the rest, assuming the directory structure is as expected.

  use mpi_basic, only: par
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: warning, crash
  use mpi_f08, only: MPI_WIN
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use netcdf_io_main, only: read_field_from_file_2D_monthly, read_field_from_file_2D
  use climate_model_basic, only: atype_climate_model
  use ISMIP7_forcing_field_types, only: type_ISMIP7_forcing_field_monthly, type_ISMIP7_forcing_field_yearly

  implicit none

  private

  public :: type_climate_model_ISMIP7

  type, extends(atype_climate_model) :: type_climate_model_ISMIP7
    !< Variables and functions that are specific to the ISMIP7 climate model

      ! Baseline climate and surface elevation
      real(dp), dimension(:,:), contiguous, pointer :: T2m_baseline    => null()   !< [K]                      Baseline monthly mean 2-m air temperature
      real(dp), dimension(:,:), contiguous, pointer :: Precip_baseline => null()   !< [m.w.e. month^-1]        Baseline monthly total precipitation
      real(dp), dimension(:  ), contiguous, pointer :: Hs_baseline     => null()   !< [m w.r.t. PD sea level]  Baseline surface elevation
      type(MPI_WIN) :: wT2m_baseline, wPrecip_baseline, wHs_baseline

      ! ISMIP7-style input forcing fields
      type(type_ISMIP7_forcing_field_monthly) :: tas                 !< [K]                GCM-derived monthly mean surface air temperature
      type(type_ISMIP7_forcing_field_monthly) :: tas_anomaly         !< [K]                GCM-derived monthly mean surface air temperature anomaly
      type(type_ISMIP7_forcing_field_monthly) :: pr                  !< [m.w.e. month^-1]  GCM-derived monthly total precipitation
      type(type_ISMIP7_forcing_field_monthly) :: pr_anomaly          !< [m.w.e. month^-1]  GCM-derived monthly total precipitation anomaly
      type(type_ISMIP7_forcing_field_yearly)  :: dtsdz               !< [K m^-1]           GCM-derived anual mean vertical temperature gradient

      ! Elevation-based temperature correction
      real(dp), dimension(:  ), contiguous, pointer :: delta_z    => null()
      real(dp), dimension(:  ), contiguous, pointer :: delta_ts   => null()
      type(MPI_WIN) :: wdelta_z, wdelta_ts

    contains

      procedure, public :: allocate   => climate_model_ISMIP7_allocate
      procedure, public :: deallocate => climate_model_ISMIP7_deallocate
      procedure, public :: initialise => climate_model_ISMIP7_initialise
      procedure, public :: run        => climate_model_ISMIP7_run
      procedure, public :: remap      => climate_model_ISMIP7_remap

      procedure, private :: initialise_climate_baseline_fixed

  end type type_climate_model_ISMIP7

contains

  subroutine climate_model_ISMIP7_allocate( self, region_name, mesh)

    ! In/output variables:
    class(type_climate_model_ISMIP7), intent(inout) :: self
    character(len=*),                 intent(in   ) :: region_name
    type(type_mesh), target,          intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_ISMIP7_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all climate models
    call self%allocate_climate_model( 'climate_ISMIP7', region_name, mesh)

    ! Allocate all the stuff that is specific to the ISMIP7 climate model

    ! Allocate fields and baseline climate
    select case (C%climate_ISMIP7_choice_baseline)
    case default
      call crash('invalid climate_ISMIP7_choice_baseline "' // trim( C%climate_ISMIP7_choice_baseline) // '"')
    case ('yearly')

      ! Allocate monthly climate (as ISMIP7 forcing fields)
      call self%tas%allocate( self, 'tas', 'Monthly mean 2-m air temperature', 'K')
      call self%pr%allocate ( self, 'pr' , 'Monthly total precipitation', 'm.w.e. month^-1')

    case ('fixed')

      ! Allocate baseline climate
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

      ! Allocate anomalies (as ISMIP7 forcing fields)
      call self%tas_anomaly%allocate( self, 'tas-anomaly', 'Monthly mean 2-m air temperature anomaly', 'K')
      call self%pr_anomaly%allocate ( self, 'pr-anomaly' , 'Monthly total precipitation anomaly', 'm.w.e. month^-1')

    end select

    ! Initialise vertical gradient (as ISMIP7 forcing field)
    call self%dtsdz%allocate( self, 'dtsdz', 'Vertical temperature gradient', 'K/m')

    ! Elevation-based temperature correction
    call self%create_field( self%Hs_baseline, self%wHs_baseline, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Hs_baseline', &
      long_name = 'Baseline surface elevation', &
      units     = 'm', &
      remap_method = 'reallocate')

    call self%create_field( self%delta_z, self%wdelta_z, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'delta_z', &
      long_name = 'Elevation difference w.r.t. baseline', &
      units     = 'm', &
      remap_method = 'reallocate')

    call self%create_field( self%delta_ts, self%wdelta_ts, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'delta_ts', &
      long_name = 'Temperature correction difference w.r.t. baseline', &
      units     = 'K', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_ISMIP7_allocate

  subroutine climate_model_ISMIP7_deallocate( self)

    ! In/output variables:
    class(type_climate_model_ISMIP7), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_ISMIP7_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is common to all climate models
    call self%deallocate_climate_model()

    ! Deallocate all the stuff that is specific to climate model ISMIP7

    ! Baseline climate and surface elevation
    nullify( self%T2m_baseline)
    nullify( self%Precip_baseline)
    nullify( self%Hs_baseline)

    ! ISMIP7-style input forcing fields
    ! call self%tas%deallocate()
    ! call self%tas_anomaly%deallocate()
    ! call self%pr%deallocate()
    ! call self%pr_anomaly%deallocate()
    ! call self%dtsdz%deallocate()

    ! Elevation-based temperature correction
    nullify( self%delta_z )
    nullify( self%delta_ts)


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_ISMIP7_deallocate

  subroutine climate_model_ISMIP7_initialise( self, refgeo_PD, refgeo_init)

    ! In/output variables:
    class(type_climate_model_ISMIP7), intent(inout) :: self
    type(type_reference_geometry),    intent(in   ) :: refgeo_PD
    type(type_reference_geometry),    intent(in   ) :: refgeo_init

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'climate_model_ISMIP7_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Initialise all the stuff that is common to all climate models
    call self%initialise_climate_model()

    ! Initialise all the stuff that is specific to climate model ISMIP7

    ! Initialise fields and baseline climate
    select case (C%climate_ISMIP7_choice_baseline)
    case default
      call crash('invalid climate_ISMIP7_choice_baseline "' // trim( C%climate_ISMIP7_choice_baseline) // '"')
    case ('yearly')

      call self%tas%initialise( C%climate_ISMIP7_forcing_foldername, C%climate_ISMIP7_forcing_version, self%mesh)
      call self%pr%initialise ( C%climate_ISMIP7_forcing_foldername, C%climate_ISMIP7_forcing_version, self%mesh)

    case ('fixed')

      call self%initialise_climate_baseline_fixed()

      call self%tas_anomaly%initialise( C%climate_ISMIP7_forcing_foldername, C%climate_ISMIP7_forcing_version, self%mesh)
      call self%pr_anomaly%initialise ( C%climate_ISMIP7_forcing_foldername, C%climate_ISMIP7_forcing_version, self%mesh)

    end select

    ! Initialise vertical gradient
    call self%dtsdz%initialise( C%climate_ISMIP7_forcing_foldername, C%climate_ISMIP7_forcing_version, self%mesh)

    ! Initialise the baseline surface elevation
    select case (C%climate_ISMIP7_choice_refgeo)
    case default
      call crash('invalid climate_ISMIP7_choice_refgeo "' // trim( C%climate_ISMIP7_choice_refgeo) // '"')
    case ('init')
      self%Hs_baseline( self%mesh%vi1:self%mesh%vi2) = refgeo_init%Hs( self%mesh%vi1: self%mesh%vi2)
    case ('PD')
      self%Hs_baseline( self%mesh%vi1:self%mesh%vi2) = refgeo_PD%Hs( self%mesh%vi1: self%mesh%vi2)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_ISMIP7_initialise

  subroutine initialise_climate_baseline_fixed( self)

    ! In/output variables:
    class(type_climate_model_ISMIP7), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'initialise_climate_baseline_fixed'
    character(len=:), allocatable :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) then
      write(0,*) '   Reading ISMIP7 climate baseline from file: ', &
        UPSY%stru%colour_string( trim( C%climate_ISMIP7_filename_baseline), 'light blue')
    end if

    ! Read the fixed baseline climate
    call read_field_from_file_2D_monthly( C%climate_ISMIP7_filename_baseline, 'T2m'   , self%mesh, C%output_dir, self%T2m_baseline)
    call read_field_from_file_2D_monthly( C%climate_ISMIP7_filename_baseline, 'Precip', self%mesh, C%output_dir, self%Precip_baseline)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_climate_baseline_fixed

  subroutine climate_model_ISMIP7_run( self, ice, time)

    ! In/output variables:
    class(type_climate_model_ISMIP7), intent(inout) :: self
    type(type_ice_model),             intent(in   ) :: ice
    real(dp),                         intent(in   ) :: time

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_ISMIP7_run'
    logical                     :: do_run_climate_model
    integer                     :: vi, mi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is common to all climate models
    call self%run_climate_model( time, do_run_climate_model)
    if (.not. do_run_climate_model) then
      call finalise_routine( routine_name)
      return
    end if

    ! Run all the stuff that is specific to climate model ISMIP7

    ! Calculate elevation-based T2m correction
    call self%dtsdz%update_and_interpolate( self%mesh, time)

    do vi = self%mesh%vi1, self%mesh%vi2
      self%delta_z ( vi) = ice%geom%Hs( vi) - self%Hs_baseline ( vi)
      self%delta_ts( vi) = self%delta_z( vi) * self%dtsdz%val_interp( vi)
    end do

    ! Calculate monthly climate
    select case (C%climate_ISMIP7_choice_baseline)
    case default
      call crash('invalid climate_ISMIP7_choice_baseline "' // trim( C%climate_ISMIP7_choice_baseline) // '"')
    case ('yearly')

      ! Update and interpolate timeframes
      call self%tas%update_and_interpolate( self%mesh, time)
      call self%pr%update_and_interpolate ( self%mesh, time)

      ! Calculate monthly climate
      do vi = self%mesh%vi1, self%mesh%vi2
        do mi = 1, 12
          self%T2m   ( vi, mi) = self%tas%val_interp( vi, mi) + self%delta_ts( vi)
          self%Precip( vi, mi) = self%pr%val_interp ( vi, mi)
        end do
      end do

    case ('fixed')

      ! Update and interpolate timeframes
      call self%tas_anomaly%update_and_interpolate( self%mesh, time)
      call self%pr_anomaly%update_and_interpolate ( self%mesh, time)

      ! Calculate monthly climate
      do vi = self%mesh%vi1, self%mesh%vi2
        do mi = 1, 12
          self%T2m   ( vi, mi) =             self%T2m_baseline   ( vi, mi) + self%tas_anomaly%val_interp( vi, mi) + self%delta_ts( vi)
          self%Precip( vi, mi) = max( 0._dp, self%Precip_baseline( vi, mi) + self%pr_anomaly%val_interp ( vi, mi))   ! Must be positive
        end do
      end do

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_ISMIP7_run

  subroutine climate_model_ISMIP7_remap( self, mesh_new)

    ! In/output variables:
    class(type_climate_model_ISMIP7), intent(inout) :: self
    type(type_mesh), target,          intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_ISMIP7_remap'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remap all the stuff that is common to all climate models
    call self%remap_climate_model( mesh_new)

    ! Remap all the stuff that is specific to climate model ISMIP7

    call crash('remapping not yet supported for ISMIP7 climate forcing')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine climate_model_ISMIP7_remap

end module ISMIP7_climate
