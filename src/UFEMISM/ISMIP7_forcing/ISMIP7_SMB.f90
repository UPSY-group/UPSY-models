module ISMIP7_SMB

  ! The ISMIP7 protocol provides separate NetCDF files for each year, with a single file contaning 12 monthly fields
  ! for that year. There are separate files for each variable, namely:
  !
  ! - acabf        : the absolute SMB (in kg m^-2 s^-1)
  ! - acabf-anomaly: the SMB anomaly relative to the 1960-1989 baseline
  ! - dacabfdz     : the vertical SMB gradient, to correct for differences in elevation between
  !                  the observed and modelled ice-sheet geometry
  !
  ! The directory structure for the forcing data is:
  !
  !  ../
  !    path/
  !      to/
  !        base_folder/
  !          acabf/
  !            version/
  !              acabf_somethingsomethingsomething_2015.nc
  !              acabf_somethingsomethingsomething_2016.nc
  !              acabf_somethingsomethingsomething_2017.nc
  !              ...
  !          acabf-anomaly/
  !            version/
  !              acabf-anomaly_somethingsomethingsomething_2015.nc
  !              acabf-anomaly_somethingsomethingsomething_2016.nc
  !              acabf-anomaly_somethingsomethingsomething_2017.nc
  !              ...
  !          dacabfdz/
  !            version/
  !              dacabfdz_somethingsomethingsomething_2015.nc
  !              dacabfdz_somethingsomethingsomething_2016.nc
  !              dacabfdz_somethingsomethingsomething_2017.nc
  !              ...
  !
  ! In the config, you only need to provide the path/to/base_folder and the version; UFEMISM will take
  ! care of the rest, assuming the directory structure is as expected.

  use mpi_basic, only: par
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use SMB_model_basic, only: atype_SMB_model, &
    type_SMB_model_context_run, &
    type_SMB_model_context_remap
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use mpi_f08, only: MPI_WIN
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use parameters, only: sec_per_year, ice_density, NaN, freshwater_density
  use ISMIP7_forcing_field_types, only: type_ISMIP7_forcing_field_monthly, type_ISMIP7_forcing_field_yearly
  use netcdf_io_main, only: read_field_from_file_2d

  implicit none

  private

  public :: type_SMB_model_ISMIP7

  type, extends(atype_SMB_model) :: type_SMB_model_ISMIP7
    !< Variables and functions that are specific to the ISMIP7 SMB model

      ! Baseline SMB and surface elevation
      real(dp), dimension(:), contiguous, pointer   :: SMB_baseline  => null()   !< [m.i.e. yr^-1]           Baseline yearly total SMB
      real(dp), dimension(:), contiguous, pointer   :: SMB_offset    => null()   !< [m.i.e. yr^-1]           Offset in SMB to shift baseline period to 1960-1989
      real(dp), dimension(:), contiguous, pointer   :: Hs_baseline   => null()   !< [m w.r.t. PD sea level]  Baseline surface elevation
      type(MPI_WIN) :: wSMB_baseline, wHs_baseline, wSMB_offset

      ! ISMIP7-style input forcing fields
      type(type_ISMIP7_forcing_field_monthly)       :: acabf                     !< [m.i.e. month^_1]        Monthly total SMB          (scaled from SI units to ice model units upon reading)
      type(type_ISMIP7_forcing_field_monthly)       :: acabf_anomaly             !< [m.i.e. month^_1]        Monthly total SMB anomaly  (scaled from SI units to ice model units upon reading)
      type(type_ISMIP7_forcing_field_yearly)        :: dacabfdz                  !< [m.i.e. yr^-1 m^-1]      Vertical SMB gradient      (scaled from SI units to ice model units upon reading)

      ! Elevation-induced SMB change
      real(dp), dimension(:), contiguous, pointer   :: delta_z       => null()   !< [m]                      Surface elevation change w.r.t. baseline
      real(dp), dimension(:), contiguous, pointer   :: delta_SMB     => null()   !< [m.i.e. yr^-1]           Elevation-induced SMB change
      type(MPI_WIN) :: wdelta_z, wdelta_SMB

      ! Monthly SMB
      real(dp), dimension(:,:), contiguous, pointer :: SMB_monthly   => null()   !< [m.i.e. month^-1]        Monthly total SMB
      type(MPI_WIN) :: wSMB_monthly

    contains

      procedure, public :: allocate   => SMB_model_ISMIP7_allocate
      procedure, public :: deallocate => SMB_model_ISMIP7_deallocate
      procedure, public :: initialise => SMB_model_ISMIP7_initialise
      procedure, public :: run_SMB_model        => run_SMB_model_ISMIP7_abs
      procedure, public :: remap_SMB_model      => remap_SMB_model_ISMIP7_abs

      procedure, private :: run_SMB_model_ISMIP7
      ! procedure, private :: remap_SMB_model_ISMIP7

      procedure, private :: initialise_SMB_baseline_fixed

  end type type_SMB_model_ISMIP7

contains

  subroutine run_SMB_model_ISMIP7_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_ISMIP7),             intent(inout) :: self
    type(type_SMB_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_SMB_model_ISMIP7_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call self%run_SMB_model_ISMIP7( self%mesh, context%ice, context%time)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_ISMIP7_abs

  subroutine remap_SMB_model_ISMIP7_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_ISMIP7),               intent(inout) :: self
    type(type_SMB_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_SMB_model_ISMIP7_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    call crash('remapping not yet supported for ISMIP7 SMB forcing')
    ! call self%remap_SMB_model_ISMIP7( context%mesh_new)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_ISMIP7_abs



  subroutine SMB_model_ISMIP7_allocate( self, region_name, mesh)

    ! In/output variables:
    class(type_SMB_model_ISMIP7), intent(inout) :: self
    character(len=*),             intent(in   ) :: region_name
    type(type_mesh), target,      intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_ISMIP7_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all SMB models
    call self%allocate_SMB_model( 'SMB_ISMIP7', region_name, mesh)

    ! Allocate all the stuff that is specific to the ISMIP7 SMB model

    ! Allocate fields and baseline SMB
    select case (C%SMB_ISMIP7_choice_SMB_baseline)
    case default
      call crash('invalid SMB_ISMIP7_choice_SMB_baseline "' // trim( C%SMB_ISMIP7_choice_SMB_baseline) // '"')
    case ('yearly')

      ! Allocate monthly climate (as ISMIP7 forcing fields)
      call self%acabf%allocate( self, 'acabf', 'Monthly total SMB', 'm.i.e. month^-1')

    case ('fixed')

      ! Allocate baseline climate
      call self%create_field( self%SMB_baseline, self%wSMB_baseline, &
        self%mesh, Arakawa_grid%a(), &
        name      = 'SMB_baseline', &
        long_name = 'Baseline yearly total SMB', &
        units     = 'm.i.e. yr^-1', &
        remap_method = 'reallocate')

      if (C%SMB_ISMIP7_choice_SMB_offset) then
        ! SMB offset
        call self%create_field( self%SMB_offset, self%wSMB_offset, &
          self%mesh, Arakawa_grid%a(), &
          name      = 'SMB_offset', &
          long_name = 'Offset to shift SMB baseline period', &
          units     = 'm.i.e. yr^-1', &
          remap_method = 'reallocate')
      end if

      ! Allocate anomalies (as ISMIP7 forcing fields)
      call self%acabf_anomaly%allocate( self, 'acabf-anomaly', 'Monthly total SMB anomaly', 'm.i.e. month^-1')

    end select

    ! Initialise vertical gradient (as ISMIP7 forcing field)
    call self%dacabfdz%allocate( self, 'dacabfdz', 'Vertical SMB gradient', 'm.i.e. yr^-1 m^-1')

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
      long_name = 'Surface elevation change w.r.t. baseline', &
      units     = 'm', &
      remap_method = 'reallocate')

    call self%create_field( self%delta_SMB, self%wdelta_SMB, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'delta_SMB', &
      long_name = 'Elevation-induced SMB change', &
      units     = 'm.i.e. yr^-1', &
      remap_method = 'reallocate')

    ! Monthly SMB
    ! ===========

    call self%create_field( self%SMB_monthly, self%wSMB_monthly, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'SMB_monthly', &
      long_name = 'Monthly total SMB', &
      units     = 'm.i.e. month^-1')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine SMB_model_ISMIP7_allocate

  subroutine SMB_model_ISMIP7_deallocate( self)

    ! In/output variables:
    class(type_SMB_model_ISMIP7), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_ISMIP7_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is common to all SMB models
    call self%deallocate_SMB_model()

    ! Deallocate all the stuff that is specific to SMB model ISMIP7

    ! Baseline SMB and surface elevation
    nullify( self%SMB_baseline)
    nullify( self%SMB_offset)
    nullify( self%Hs_baseline)

    ! ! ISMIP7-style input forcing fields
    ! call self%acabf%deallocate()
    ! call self%acabf_anomaly%deallocate()
    ! call self%dacabfdz%deallocate()

    ! Elevation-induced SMB change
    nullify( self%delta_z)
    nullify( self%delta_SMB)

    ! Monthly SMB
    nullify( self%SMB_monthly)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine SMB_model_ISMIP7_deallocate

  subroutine SMB_model_ISMIP7_initialise( self, ice, refgeo_init, refgeo_PD)

    ! In/output variables
    class(type_SMB_model_ISMIP7),  intent(inout) :: self
    type(type_ice_model),          intent(in   ) :: ice
    type(type_reference_geometry), intent(in   ) :: refgeo_init
    type(type_reference_geometry), intent(in   ) :: refgeo_PD

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_ISMIP7_initialise'

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise all the stuff that is common to all SMB models
    call self%initialise_SMB_model()

    ! Deallocate all the stuff that is specific to SMB model ISMIP7

    ! Initialise fields and baseline climate
    select case (C%SMB_ISMIP7_choice_SMB_baseline)
    case default
      call crash('invalid SMB_ISMIP7_choice_SMB_baseline "' // trim( C%SMB_ISMIP7_choice_SMB_baseline) // '"')
    case ('yearly')

      call self%acabf%initialise( C%SMB_ISMIP7_forcing_foldername, C%SMB_ISMIP7_forcing_version, self%mesh)

    case ('fixed')

      call self%initialise_SMB_baseline_fixed()

      call self%acabf_anomaly%initialise( C%SMB_ISMIP7_forcing_foldername, C%SMB_ISMIP7_forcing_version, self%mesh)

    end select

    ! Initialise vertical gradient
    call self%dacabfdz%initialise( C%SMB_ISMIP7_forcing_foldername, C%SMB_ISMIP7_forcing_version, self%mesh)

    ! Initialise the baseline surface elevation
    select case (C%SMB_ISMIP7_choice_refgeo)
    case default
      call crash('invalid SMB_ISMIP7_choice_refgeo "' // trim( C%SMB_ISMIP7_choice_refgeo) // '"')
    case ('init')
      self%Hs_baseline( self%mesh%vi1:self%mesh%vi2) = refgeo_init%Hs( self%mesh%vi1: self%mesh%vi2)
    case ('PD')
      self%Hs_baseline( self%mesh%vi1:self%mesh%vi2) = refgeo_PD%Hs( self%mesh%vi1: self%mesh%vi2)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine SMB_model_ISMIP7_initialise

  subroutine initialise_SMB_baseline_fixed( self)

    ! In/output variables
    class(type_SMB_model_ISMIP7), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_SMB_baseline_fixed'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then
      write(0,*) '   Reading ISMIP7 SMB baseline from file: ', &
        UPSY%stru%colour_string( trim( C%SMB_ISMIP7_filename_SMB_baseline_fixed), 'light blue')
    end if

    call read_field_from_file_2D( trim( C%SMB_ISMIP7_filename_SMB_baseline_fixed), 'SMB||surface_mass_balance', &
      self%mesh, C%output_dir, self%SMB_baseline)

    if (C%SMB_ISMIP7_choice_SMB_offset) then
      if (par%primary) then
        write(0,*) '   Reading SMB offset from file: ', &
          UPSY%stru%colour_string( trim( C%SMB_ISMIP7_filename_SMB_offset), 'light blue')
      end if

      call read_field_from_file_2D( trim( C%SMB_ISMIP7_filename_SMB_offset), &
        'SMB||surface_mass_balance', self%mesh, C%output_dir, self%SMB_offset)

      ! Subtract offset from SMB_baseline
      self%SMB_baseline( self%mesh%vi1:self%mesh%vi2) = self%SMB_baseline( self%mesh%vi1:self%mesh%vi2) &
        - self%SMB_offset( self%mesh%vi1:self%mesh%vi2)

    end if

    ! Convert from [m.w.e. yr^-1] to [m.i.e. yr^-1]
    self%SMB_baseline( self%mesh%vi1:self%mesh%vi2) = self%SMB_baseline( self%mesh%vi1:self%mesh%vi2) &
      * freshwater_density / ice_density

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_baseline_fixed

  subroutine run_SMB_model_ISMIP7( self, mesh, ice, time)

    ! In/output variables:
    class(type_SMB_model_ISMIP7), intent(inout) :: self
    type(type_mesh),              intent(in   ) :: mesh
    type(type_ice_model),         intent(in   ) :: ice
    real(dp),                     intent(in   ) :: time

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_SMB_model_ISMIP7'
    real(dp)                    :: delta_z
    integer                     :: vi, mi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Calculate elevation-based T2m correction
    call self%dacabfdz%update_and_interpolate( mesh, time)

    do vi = mesh%vi1, mesh%vi2
      self%delta_z  ( vi) = ice%Hs( vi) - self%Hs_baseline ( vi)
      self%delta_SMB( vi) = self%delta_z( vi) * self%dacabfdz%val_interp( vi)
    end do

    ! Calculate monthly climate
    select case (C%SMB_ISMIP7_choice_SMB_baseline)
    case default
      call crash('invalid SMB_ISMIP7_choice_SMB_baseline "' // trim( C%SMB_ISMIP7_choice_SMB_baseline) // '"')
    case ('yearly')

      ! Update and interpolate timeframes
      call self%acabf%update_and_interpolate( mesh, time)

      ! Calculate monthly SMB
      do vi = mesh%vi1, mesh%vi2
        do mi = 1, 12
                                      ! Divide delta by 12 to convert from [m.i.e. yr^-1] to [m.i.e. month^-1]
          self%SMB_monthly( vi, mi) = self%acabf%val_interp( vi, mi) + self%delta_SMB( vi) / 12._dp
        end do
      end do

    case ('fixed')

      ! Update and interpolate timeframes
      call self%acabf_anomaly%update_and_interpolate( mesh, time)

      ! Calculate monthly climate
      do vi = mesh%vi1, mesh%vi2
        do mi = 1, 12
                                      ! Divide baseline and delta by 12 to convert from [m.i.e. yr^-1] to [m.i.e. month^-1]
          self%SMB_monthly( vi, mi) = (self%SMB_baseline( vi) + self%delta_SMB( vi)) / 12._dp + self%acabf_anomaly%val_interp( vi, mi)
        end do
      end do

    end select

    ! Calculate yearly SMB by summing all monthly values
    do vi = mesh%vi1, mesh%vi2
      self%SMB( vi) = sum( self%SMB_monthly( vi,:))
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_ISMIP7

end module ISMIP7_SMB
