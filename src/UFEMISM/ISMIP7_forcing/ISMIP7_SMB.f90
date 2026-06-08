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
  use SMB_model_basic, only: atype_SMB_model, type_SMB_model_context_allocate, &
    type_SMB_model_context_initialise, type_SMB_model_context_run, &
    type_SMB_model_context_remap
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use mpi_f08, only: MPI_WIN
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use netcdf_io_main, only: read_field_from_file_2D, read_time_from_file, &
    open_existing_netcdf_file_for_reading, setup_xy_grid_from_file, close_netcdf_file, &
    inquire_var, read_var_primary, inquire_fill_value
  use basic_model_utilities, only: list_files_in_folder
  use parameters, only: sec_per_year, ice_density, NaN
  use grid_types, only: type_grid
  use remapping_types, only: type_map
  use remapping_grid_to_mesh_vertices, only: create_map_from_xy_grid_to_mesh_vertices
  use mpi_distributed_memory_grid, only: distribute_gridded_data_from_primary
  use apply_maps, only: apply_map_xy_grid_to_mesh_2D, apply_map_xy_grid_to_mesh_3D
  use smooth_gridded_data, only: extrapolate_fillvalue_Gaussian_grid
  use dist_to_hybrid_mod, only: dist_to_hybrid

  implicit none

  private

  public :: type_SMB_model_ISMIP7

  type type_SMB_ISMIP7_timeframe
    real(dp)                                      :: time
    real(dp), dimension(:,:), contiguous, pointer :: acabf         => null()   !< [m yr^-1]      monthly surface mass balance flux          (scaled from SI units to ice model units upon reading)
    real(dp), dimension(:,:), contiguous, pointer :: acabf_anomaly => null()   !< [m yr^-1]      monthly surface mass balance flux anomaly  (scaled from SI units to ice model units upon reading)
    real(dp), dimension(:  ), contiguous, pointer :: dacabfdz      => null()   !< [m yr^-1 m^-1]         surface mass balance flux gradient (scaled from SI units to ice model units upon reading)
    type(MPI_WIN) :: wacabf, wacabf_anomaly, wdacabfdz
  end type type_SMB_ISMIP7_timeframe

  type, extends(atype_SMB_model) :: type_SMB_model_ISMIP7

      ! Mapping object to remap data from the ISMIP7 grid to the UFEMISM mesh
      type(type_grid) :: grid_raw
      type(type_map)  :: map

      ! Baseline SMB and surface elevation
      real(dp), dimension(:), contiguous, pointer :: SMB_baseline => null()   !< [m.i.e./yr]             baseline annual mean SMB
      real(dp), dimension(:), contiguous, pointer :: Hs_baseline  => null()   !< [m w.r.t. PD sea level] baseline surface elevation
      type(MPI_WIN) :: wSMB_baseline, wHs_baseline

      ! List of forcing files and their timestamps
      character(len=1024), dimension(:), allocatable :: filenames_acabf
      character(len=1024), dimension(:), allocatable :: filenames_acabf_anomaly
      character(len=1024), dimension(:), allocatable :: filenames_dacabfdz
      real(dp),            dimension(:), allocatable :: timestamps

      ! Timeframes of forcng fields enveloping the current model time
      type(type_SMB_ISMIP7_timeframe) :: timeframe_before
      type(type_SMB_ISMIP7_timeframe) :: timeframe_after
      type(type_SMB_ISMIP7_timeframe) :: timeframe_interp

      ! Elevation-induced SMB change
      real(dp), dimension(:), contiguous, pointer :: delta_z    => null()   !< [m]       surface elevation change w.r.t. baseline
      real(dp), dimension(:), contiguous, pointer :: delta_SMB  => null()   !< [m yr^-1] elevation-induced change in SMB w.r.t. baseline
      type(MPI_WIN) :: wdelta_z, wdelta_SMB

      ! Monthly SMB
      real(dp), dimension(:,:), contiguous, pointer :: SMB_monthly => null()   !< [m yr^-1] monthly surface mass balance
      type(MPI_WIN) :: wSMB_monthly

    contains

      procedure, public :: allocate_SMB_model   => allocate_SMB_model_ISMIP7_abs
      procedure, public :: deallocate_SMB_model => deallocate_SMB_model_ISMIP7_abs
      procedure, public :: initialise_SMB_model => initialise_SMB_model_ISMIP7_abs
      procedure, public :: run_SMB_model        => run_SMB_model_ISMIP7_abs
      procedure, public :: remap_SMB_model      => remap_SMB_model_ISMIP7_abs

      procedure, private :: allocate_SMB_model_ISMIP7
      procedure, private :: initialise_SMB_model_ISMIP7
      procedure, private :: run_SMB_model_ISMIP7
      procedure, private :: remap_SMB_model_ISMIP7

      procedure, private :: allocate_SMB_model_ISMIP7_timeframe
      procedure, private :: initialise_SMB_baseline_fixed
      procedure, private :: initialise_Hs_baseline
      procedure, private :: initialise_lists_of_files_and_timestamps
      procedure, private :: initialise_list_of_files_and_timestamps
      procedure, private :: initialise_remapping_object
      procedure, private :: update_timeframes
      procedure, private :: update_timeframe
      procedure, private :: read_monthly_data_from_ISMIP7_forcing_file
      procedure, private :: read_annual_mean_data_from_ISMIP7_forcing_file
      procedure, private :: reallocate_timeframe

  end type type_SMB_model_ISMIP7

contains

  subroutine allocate_SMB_model_ISMIP7_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_ISMIP7),                  intent(inout) :: self
    type(type_SMB_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_SMB_model_ISMIP7_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%allocate_SMB_model_ISMIP7

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_ISMIP7_abs

  subroutine deallocate_SMB_model_ISMIP7_abs( self)

    ! In/output variables:
    class(type_SMB_model_ISMIP7), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_SMB_model_ISMIP7_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_SMB_model_ISMIP7_abs

  subroutine initialise_SMB_model_ISMIP7_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_ISMIP7),                    intent(inout) :: self
    type(type_SMB_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_SMB_model_ISMIP7_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call self%initialise_SMB_model_ISMIP7( self%mesh, context%refgeo_init, context%refgeo_PD, self%region_name())

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_ISMIP7_abs

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

    call self%remap_SMB_model_ISMIP7( context%mesh_new, context%time, &
      context%region_name, context%refgeo_init, context%refgeo_PD, context%ice)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_ISMIP7_abs



  subroutine allocate_SMB_model_ISMIP7( self)

    ! In/output variables
    class(type_SMB_model_ISMIP7), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_SMB_model_ISMIP7'

    ! Add routine to path
    call init_routine( routine_name)

    ! Baseline
    ! ========

    call self%create_field( self%SMB_baseline, self%wSMB_baseline, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'SMB_baseline', &
      long_name = 'Baseline SMB', &
      units     = 'm yr^-1')

    call self%create_field( self%Hs_baseline, self%wHs_baseline, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Hs_baseline', &
      long_name = 'Baseline surface elevation', &
      units     = 'm')

    ! Timeframes
    ! ==========

    call self%allocate_SMB_model_ISMIP7_timeframe( self%timeframe_before, 'before')
    call self%allocate_SMB_model_ISMIP7_timeframe( self%timeframe_after , 'after')
    call self%allocate_SMB_model_ISMIP7_timeframe( self%timeframe_interp, 'interp')

    ! Elevation-induced SMB change
    ! ============================

    call self%create_field( self%delta_z, self%wdelta_z, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'delta_z', &
      long_name = 'Surface elevation change w.r.t. baseline', &
      units     = 'm')

    call self%create_field( self%delta_SMB, self%wdelta_SMB, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'delta_SMB', &
      long_name = 'elevation-induced change in SMB w.r.t. baseline', &
      units     = 'm yr^-1')

    ! Monthly SMB
    ! ===========

    call self%create_field( self%SMB_monthly, self%wSMB_monthly, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'SMB_monthly', &
      long_name = 'monthly surface mass balance', &
      units     = 'm yr^-1')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_ISMIP7

  subroutine allocate_SMB_model_ISMIP7_timeframe( self, timeframe, name_postfix)

    ! In/output variables
    class(type_SMB_model_ISMIP7),    intent(inout) :: self
    type(type_SMB_ISMIP7_timeframe), intent(inout) :: timeframe
    character(len=*),                intent(in   ) :: name_postfix

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_SMB_model_ISMIP7_timeframe'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%SMB_ISMIP7_choice_SMB_baseline)
    case default
      call crash('invalid SMB_ISMIP7_choice_SMB_baseline "' // trim( C%SMB_ISMIP7_choice_SMB_baseline) // '"')
    case ('yearly')
      call self%create_field( timeframe%acabf, timeframe%wacabf, &
        self%mesh, Arakawa_grid%a(), third_dimension%month(), &
        name      = 'acabf_' // trim( name_postfix), &
        long_name = 'monthly surface mass balance flux', &
        units     = 'm yr^-1')
    case ('fixed')
      call self%create_field( timeframe%acabf_anomaly, timeframe%wacabf_anomaly, &
        self%mesh, Arakawa_grid%a(), third_dimension%month(), &
        name      = 'acabf_anomaly_' // trim( name_postfix), &
        long_name = 'monthly surface mass balance flux anomaly', &
        units     = 'm yr^-1')
    end select

    call self%create_field( timeframe%dacabfdz, timeframe%wdacabfdz, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'dacabfdz_' // trim( name_postfix), &
      long_name = 'monthly surface mass balance flux gradient', &
      units     = 'm yr^-1 m^-1')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_ISMIP7_timeframe

  subroutine initialise_SMB_model_ISMIP7( self, mesh, refgeo_init, refgeo_PD, region_name)

    ! In/output variables
    class(type_SMB_model_ISMIP7),  intent(inout) :: self
    type(type_mesh),               intent(in   ) :: mesh
    type(type_reference_geometry), intent(in   ) :: refgeo_init, refgeo_PD
    character(len=3),              intent(in   ) :: region_name

    ! Local variables:
    character(len=*), parameter            :: routine_name = 'initialise_SMB_model_ISMIP7'
    type(type_reference_geometry), pointer :: refgeo
    real(dp), dimension(mesh%vi1:mesh%vi2) :: dummy

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise the baseline SMB
    select case (C%SMB_ISMIP7_choice_SMB_baseline)
    case default
      call crash('invalid SMB_ISMIP7_choice_SMB_baseline "' // trim( C%SMB_ISMIP7_choice_SMB_baseline) // '"')
    case ('yearly')
      ! No need to do anything
    case ('fixed')
      call self%initialise_SMB_baseline_fixed( mesh, region_name)
    end select

    ! Initialise the baseline surface elevation
    select case (C%SMB_ISMIP7_choice_refgeo)
    case default
      call crash('invalid SMB_ISMIP7_choice_refgeo "' // trim( C%SMB_ISMIP7_choice_refgeo) // '"')
    case ('init')
      call self%initialise_Hs_baseline( mesh, refgeo_init)
    case ('PD')
      call self%initialise_Hs_baseline( mesh, refgeo_PD)
    end select

    ! Initialise lists of forcing files and their timestamps
    call self%initialise_lists_of_files_and_timestamps
    call self%initialise_remapping_object( mesh)

    ! Set the timestamps of the two timeframes so that the first time the model is run,
    ! it will automatically read and update them
    self%timeframe_before%time = C%start_time_of_run - 2._dp
    self%timeframe_after%time  = C%start_time_of_run - 1._dp

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_ISMIP7

  subroutine initialise_remapping_object( self, mesh)

    ! In/output variables
    class(type_SMB_model_ISMIP7), intent(inout) :: self
    type(type_mesh),              intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'initialise_remapping_object'
    character(len=:), allocatable :: filename
    integer                       :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! Read the grid from the first listed input file
    filename = trim( self%filenames_dacabfdz(1))
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_xy_grid_from_file( filename, ncid, self%grid_raw)
    call close_netcdf_file( ncid)

    ! Calculate the remapping operator
    self%grid_raw%name = 'grid_ISMIP7_SMB_files'
    call create_map_from_xy_grid_to_mesh_vertices( self%grid_raw, mesh, C%output_dir, self%map, '2nd_order_conservative')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_remapping_object

  subroutine initialise_SMB_baseline_fixed( self, mesh, region_name)

    ! In/output variables
    class(type_SMB_model_ISMIP7), intent(inout) :: self
    type(type_mesh),              intent(in   ) :: mesh
    character(len=3),             intent(in   ) :: region_name

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_SMB_baseline_fixed'

    ! Add routine to path
    call init_routine( routine_name)

    call read_field_from_file_2D( trim( C%SMB_ISMIP7_filename_SMB_baseline_fixed), 'SMB||surface_mass_balance', &
      mesh, C%output_dir, self%SMB_baseline)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_baseline_fixed

  subroutine initialise_Hs_baseline( self, mesh, refgeo)

    ! In/output variables
    class(type_SMB_model_ISMIP7),  intent(inout) :: self
    type(type_mesh),               intent(in   ) :: mesh
    type(type_reference_geometry), intent(in   ) :: refgeo

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_Hs_baseline'

    ! Add routine to path
    call init_routine( routine_name)

    self%Hs_baseline( mesh%vi1:mesh%vi2) = refgeo%Hs( mesh%vi1: mesh%vi2)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_Hs_baseline

  subroutine initialise_lists_of_files_and_timestamps( self)

    ! In/output variables
    class(type_SMB_model_ISMIP7), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_lists_of_files_and_timestamps'

    ! Add routine to path
    call init_routine( routine_name)

    call self%initialise_list_of_files_and_timestamps( self%filenames_acabf        , self%timestamps, 'acabf')
    call self%initialise_list_of_files_and_timestamps( self%filenames_acabf_anomaly, self%timestamps, 'acabf-anomaly')
    call self%initialise_list_of_files_and_timestamps( self%filenames_dacabfdz     , self%timestamps, 'dacabfdz')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_lists_of_files_and_timestamps

  subroutine initialise_list_of_files_and_timestamps( self, filenames, timestamps, var_name)

    ! In/output variables
    class(type_SMB_model_ISMIP7),                   intent(inout) :: self
    character(len=1024), dimension(:), allocatable, intent(inout) :: filenames
    real(dp),            dimension(:), allocatable, intent(inout) :: timestamps
    character(len=*),                               intent(in   ) :: var_name

    ! Local variables:
    character(len=*), parameter                    :: routine_name = 'initialise_list_of_files_and_timestamps'
    character(len=:), allocatable                  :: foldername
    integer                                        :: i
    real(dp)                                       :: year

    ! Add routine to path
    call init_routine( routine_name)

    if (allocated( filenames )) deallocate( filenames)
    if (allocated( timestamps)) deallocate( timestamps)

    ! Construct foldernames
    foldername = trim( C%SMB_ISMIP7_forcing_foldername) // '/' // var_name // '/' // trim( C%SMB_ISMIP7_forcing_version)

    call list_files_in_folder( foldername, filenames, var_name)
    if (size( filenames,1) == 0) call crash('could not find any valid NetCDF files in directory "' // trim( foldername) // '"')

    ! Read timestamps
    allocate( timestamps( size( filenames,1)))
    do i = 1, size( filenames,1)
      call read_year_from_netcdf_filename( filenames( i), year)
      ! Add half a year, so that we have the timestamp at the middle of the year rather than the start
      timestamps( i) = year + 0.5_dp
    end do

    ! Append foldername to filenames
    do i = 1, size( filenames,1)
      filenames( i) = trim( foldername) // '/' // trim( filenames( i))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_list_of_files_and_timestamps

  subroutine read_year_from_netcdf_filename( filename, year)
    ! Rather than trying to read from the file's time dimension, which would mean
    ! dealing with whatever calendar it uses, just read it from the filename

    ! In/output variables
    character(len=*), intent(in   ) :: filename
    real(dp),         intent(  out) :: year

    ! Local variables
    character(len=*), parameter :: routine_name = 'read_year_from_netcdf_filename'
    character(len=4) :: year_str
    integer          :: year_int, stat

    ! Add routine to path
    call init_routine( routine_name)

    ! Since the files are called e.g. 'acabf_AIS_CESM2-WACCM_ssp585_SDBN1-8000m_v2_2015.nc',
    ! just assumed that the last few characters spell out the year

    year_str = filename( len_trim( filename) - 6 : len_trim( filename) - 3)
    year_int = UPSY%stru%str2int( year_str, stat)
    if (stat /= 0) call crash('could not read year from filename "' // trim( filename) // '"')

    year = real( year_int,dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_year_from_netcdf_filename

  subroutine run_SMB_model_ISMIP7( self, mesh, ice, time)

    ! In/output variables:
    class(type_SMB_model_ISMIP7), intent(inout) :: self
    type(type_mesh),              intent(in   ) :: mesh
    type(type_ice_model),         intent(in   ) :: ice
    real(dp),                     intent(in   ) :: time

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_SMB_model_ISMIP7'
    real(dp)                    :: w_before, w_after, delta_z
    integer                     :: vi, mi

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Update timeframes if necessary
    if (time < self%timeframe_before%time .or. time > self%timeframe_after%time) then
      call self%update_timeframes( mesh, time)
    end if

    ! Interpolate timeframes
    w_before = (self%timeframe_after%time - time) / (self%timeframe_after%time - self%timeframe_before%time)
    w_after  = 1._dp - w_before

    select case (C%SMB_ISMIP7_choice_SMB_baseline)
    case default
      call crash('invalid SMB_ISMIP7_choice_SMB_baseline "' // trim( C%SMB_ISMIP7_choice_SMB_baseline) // '"')
    case ('yearly')
      do vi = mesh%vi1, mesh%vi2
        self%timeframe_interp%acabf( vi,:) = &
          w_before * self%timeframe_before%acabf( vi,:) + &
          w_after  * self%timeframe_after%acabf ( vi,:)
        self%timeframe_interp%dacabfdz( vi) = &
          w_before * self%timeframe_before%dacabfdz( vi) + &
          w_after  * self%timeframe_after%dacabfdz ( vi)
      end do
    case ('fixed')
      do vi = mesh%vi1, mesh%vi2
        self%timeframe_interp%acabf_anomaly( vi,:) = &
          w_before * self%timeframe_before%acabf_anomaly( vi,:) + &
          w_after  * self%timeframe_after%acabf_anomaly ( vi,:)
        self%timeframe_interp%dacabfdz( vi) = &
          w_before * self%timeframe_before%dacabfdz( vi) + &
          w_after  * self%timeframe_after%dacabfdz ( vi)
      end do
    end select

    ! Calculate elevation-based SMB correction
    do vi = mesh%vi1, mesh%vi2
      self%delta_z( vi) = ice%Hs( vi) - self%Hs_baseline( vi)
      self%delta_SMB( vi) = self%delta_z( vi) * self%timeframe_interp%dacabfdz( vi)
    end do

    ! Calculate monthly SMB
    select case (C%SMB_ISMIP7_choice_SMB_baseline)
    case default
      call crash('invalid SMB_ISMIP7_choice_SMB_baseline "' // trim( C%SMB_ISMIP7_choice_SMB_baseline) // '"')
    case ('yearly')
      do vi = mesh%vi1, mesh%vi2
        do mi = 1, 12
          self%SMB_monthly( vi,mi) = self%timeframe_interp%acabf( vi,mi) + self%delta_SMB( vi)
        end do
      end do
    case ('fixed')
      do vi = mesh%vi1, mesh%vi2
        do mi = 1, 12
          self%SMB_monthly( vi,mi) = self%SMB_baseline( vi) + self%timeframe_interp%acabf_anomaly( vi,mi) + self%delta_SMB( vi)
        end do
      end do
    end select

    ! Calculate applied SMB
    do vi = mesh%vi1, mesh%vi2
      self%SMB( vi) = sum( self%SMB_monthly( vi,:)) / 12._dp
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_ISMIP7

  subroutine update_timeframes( self, mesh, time)

    ! In/output variables:
    class(type_SMB_model_ISMIP7), intent(inout) :: self
    type(type_mesh),              intent(in   ) :: mesh
    real(dp),                     intent(in   ) :: time

    ! Local variables:
    character(len=*), parameter :: routine_name = 'update_timeframes'
    real(dp)                    :: forcing_time_min, forcing_time_max
    integer                     :: ti0, ti1

    ! Add routine to call stack
    call init_routine( routine_name)

    forcing_time_min = minval( self%timestamps)
    forcing_time_max = maxval( self%timestamps)

    if (time < forcing_time_min) then
      ! Set forcing equal to the first timeframe

      call self%update_timeframe( mesh, self%timeframe_before, 1)
      call self%update_timeframe( mesh, self%timeframe_after , 1)
      self%timeframe_before%time = C%start_time_of_run
      self%timeframe_after%time  = forcing_time_min

    elseif (time > forcing_time_max) then
      ! Set forcing equal to the last timeframe

      call self%update_timeframe( mesh, self%timeframe_before, size( self%timestamps))
      call self%update_timeframe( mesh, self%timeframe_after , size( self%timestamps))
      self%timeframe_before%time = forcing_time_max
      self%timeframe_after%time  = C%end_time_of_run

    else
      ! Read both timeframes

      ti1 = 2
      do while (self%timestamps( ti1) < time)
        ti1 = ti1 + 1
      end do
      ti0 = ti1 - 1

      call self%update_timeframe( mesh, self%timeframe_before, ti0)
      call self%update_timeframe( mesh, self%timeframe_after , ti1)
      self%timeframe_before%time = self%timestamps( ti0)
      self%timeframe_after%time  = self%timestamps( ti1)

    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine update_timeframes

  subroutine update_timeframe( self, mesh, timeframe, ti)
    ! Read data from the ti'th NetCDF files to this timeframe

    ! In/output variables:
    class(type_SMB_model_ISMIP7),    intent(in   ) :: self
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_SMB_ISMIP7_timeframe), intent(inout) :: timeframe
    integer,                         intent(in   ) :: ti

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'update_timeframe'
    character(len=:), allocatable :: filename_acabf
    character(len=:), allocatable :: filename_acabf_anomaly
    character(len=:), allocatable :: filename_dacabfdz

    ! Add routine to call stack
    call init_routine( routine_name)

    select case (C%SMB_ISMIP7_choice_SMB_baseline)
    case default
      call crash('invalid SMB_ISMIP7_choice_SMB_baseline "' // trim( C%SMB_ISMIP7_choice_SMB_baseline) // '"')
    case ('yearly')

      filename_acabf         = self%filenames_acabf        ( ti)
      filename_dacabfdz      = self%filenames_dacabfdz     ( ti)

      if (par%primary) then
        write(0,*) '   Reading ISMIP7 SMB forcing data from files:'
        write(0,*) '    ', UPSY%stru%colour_string( trim( filename_acabf)        , 'light blue')
        write(0,*) '    ', UPSY%stru%colour_string( trim( filename_dacabfdz)     , 'light blue')
      end if

      call self%read_monthly_data_from_ISMIP7_forcing_file    ( mesh, filename_acabf   , 'acabf'        , timeframe%acabf)
      call self%read_annual_mean_data_from_ISMIP7_forcing_file( mesh, filename_dacabfdz, 'dacabfdz'     , timeframe%dacabfdz)

      ! Convert from SI units (kg m^-2 s^-1) to ice model units (m yr^-1)
      timeframe%acabf        ( mesh%vi1: mesh%vi2,:) = timeframe%acabf        ( mesh%vi1: mesh%vi2,:) * sec_per_year / ice_density
      timeframe%dacabfdz     ( mesh%vi1: mesh%vi2  ) = timeframe%dacabfdz     ( mesh%vi1: mesh%vi2  ) * sec_per_year / ice_density

    case ('fixed')

      filename_acabf_anomaly = self%filenames_acabf_anomaly( ti)
      filename_dacabfdz      = self%filenames_dacabfdz     ( ti)

      if (par%primary) then
        write(0,*) '   Reading ISMIP7 SMB forcing data from files:'
        write(0,*) '    ', UPSY%stru%colour_string( trim( filename_acabf_anomaly), 'light blue')
        write(0,*) '    ', UPSY%stru%colour_string( trim( filename_dacabfdz)     , 'light blue')
      end if

      call self%read_monthly_data_from_ISMIP7_forcing_file    ( mesh, filename_acabf_anomaly, 'acabf-anomaly', timeframe%acabf_anomaly)
      call self%read_annual_mean_data_from_ISMIP7_forcing_file( mesh, filename_dacabfdz     , 'dacabfdz'     , timeframe%dacabfdz)

      ! Convert from SI units (kg m^-2 s^-1) to ice model units (m yr^-1)
      timeframe%acabf_anomaly( mesh%vi1: mesh%vi2,:) = timeframe%acabf_anomaly( mesh%vi1: mesh%vi2,:) * sec_per_year / ice_density
      timeframe%dacabfdz     ( mesh%vi1: mesh%vi2  ) = timeframe%dacabfdz     ( mesh%vi1: mesh%vi2  ) * sec_per_year / ice_density

    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine update_timeframe

  subroutine read_monthly_data_from_ISMIP7_forcing_file( self, mesh, filename, var_name, d)
    ! Since these files don't have an UPSY-style month dimension,
    ! we have to work around this a little bit...

    ! In/output variables:
    class(type_SMB_model_ISMIP7), intent(in   ) :: self
    type(type_mesh),              intent(in   ) :: mesh
    character(len=*),             intent(in   ) :: filename, var_name
    real(dp), dimension(:,:),     intent(  out) :: d

    ! Local variables:
    character(len=*), parameter               :: routine_name = 'read_monthly_data_from_ISMIP7_forcing_file'
    real(dp), dimension(:,:,:), allocatable   :: d_grid_tot
    integer                                   :: ncid, id_var
    real(dp)                                  :: fill_value
    real(dp), dimension(:,:  ), allocatable   :: d_grid_vec_partial
    real(dp)                                  :: sigma
    real(dp), dimension(mesh%vi1:mesh%vi2,12) :: d_dist

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Read raw gridded data to the primary
    if (par%primary) then
      allocate( d_grid_tot( self%grid_raw%nx, self%grid_raw%ny, 12), source = NaN)
    else
      allocate( d_grid_tot(0,0,0))
    end if

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call inquire_fill_value( filename, ncid, var_name, fill_value)
    call inquire_var( filename, ncid, var_name, id_var)
    call read_var_primary( filename, ncid, id_var, d_grid_tot)
    call close_netcdf_file( ncid)

    ! Distribute gridded data to the processes
    allocate( d_grid_vec_partial( self%grid_raw%pai%i1: self%grid_raw%pai%i2, 12))
    call distribute_gridded_data_from_primary( self%grid_raw, d_grid_vec_partial, d_grid_tot)
    deallocate( d_grid_tot)

    ! Extrapolate data into fill_value cells
    sigma = self%grid_raw%dx * 2._dp
    call extrapolate_fillvalue_Gaussian_grid( self%grid_raw, d_grid_vec_partial, sigma, fill_value)

    ! Remap data to mesh
    call apply_map_xy_grid_to_mesh_3D( self%grid_raw, mesh, self%map, d_grid_vec_partial, d_dist)
    call dist_to_hybrid( mesh%pai_V, 12, d_dist, d)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_monthly_data_from_ISMIP7_forcing_file

  subroutine read_annual_mean_data_from_ISMIP7_forcing_file( self, mesh, filename, var_name, d)
    ! Even though it's a file with only one frame, it still has a time dimension...

    ! In/output variables:
    class(type_SMB_model_ISMIP7), intent(in   ) :: self
    type(type_mesh),              intent(in   ) :: mesh
    character(len=*),             intent(in   ) :: filename, var_name
    real(dp), dimension(:),       intent(  out) :: d

    ! Local variables:
    character(len=*), parameter             :: routine_name = 'read_annual_mean_data_from_ISMIP7_forcing_file'
    real(dp), dimension(:,:,:), allocatable :: d_grid_tot_with_time
    real(dp), dimension(:,:  ), allocatable :: d_grid_tot
    integer                                 :: ncid, id_var
    real(dp)                                :: fill_value
    real(dp), dimension(:    ), allocatable :: d_grid_vec_partial
    real(dp)                                :: sigma
    real(dp), dimension(mesh%vi1:mesh%vi2)  :: d_dist

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Read raw gridded data to the primary
    if (par%primary) then
      allocate( d_grid_tot_with_time( self%grid_raw%nx, self%grid_raw%ny, 1), source = NaN)
      allocate( d_grid_tot          ( self%grid_raw%nx, self%grid_raw%ny   ), source = NaN)
    else
      allocate( d_grid_tot_with_time(0,0,0))
      allocate( d_grid_tot          (0,0  ))
    end if

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call inquire_fill_value( filename, ncid, var_name, fill_value)
    call inquire_var( filename, ncid, var_name, id_var)
    call read_var_primary( filename, ncid, id_var, d_grid_tot_with_time)
    if (par%primary) d_grid_tot = d_grid_tot_with_time( :,:,1)
    deallocate( d_grid_tot_with_time)
    call close_netcdf_file( ncid)

    ! Distribute gridded data to the processes
    allocate( d_grid_vec_partial( self%grid_raw%pai%i1: self%grid_raw%pai%i2))
    call distribute_gridded_data_from_primary( self%grid_raw, d_grid_vec_partial, d_grid_tot)
    deallocate( d_grid_tot)

    ! Extrapolate data into fill_value cells
    sigma = self%grid_raw%dx * 2._dp
    call extrapolate_fillvalue_Gaussian_grid( self%grid_raw, d_grid_vec_partial, sigma, fill_value)

    ! Remap data to mesh
    call apply_map_xy_grid_to_mesh_2D( self%grid_raw, mesh, self%map, d_grid_vec_partial, d_dist)
    call dist_to_hybrid( mesh%pai_V, d_dist, d)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_annual_mean_data_from_ISMIP7_forcing_file

  subroutine remap_SMB_model_ISMIP7( self, mesh_new, time, region_name, refgeo_init, refgeo_PD, ice)

    ! In/output variables
    class(type_SMB_model_ISMIP7),  intent(inout) :: self
    type(type_mesh),               intent(in   ) :: mesh_new
    real(dp),                      intent(in   ) :: time
    character(len=3),              intent(in   ) :: region_name
    type(type_reference_geometry), intent(in   ) :: refgeo_init, refgeo_PD
    type(type_ice_model),          intent(in   ) :: ice

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_SMB_model_ISMIP7'

    ! Add routine to path
    call init_routine( routine_name)

    ! Reallocate all the fields
    ! =========================

    ! Baseline SMB and surface elevation
    call self%reallocate_field( mesh_new, 'SMB_baseline', self%SMB_baseline)
    call self%reallocate_field( mesh_new, 'Hs_baseline' , self%Hs_baseline )

    ! Timeframes of forcng fields enveloping the current model time
    call self%reallocate_timeframe( self%timeframe_before, 'before', mesh_new)
    call self%reallocate_timeframe( self%timeframe_after , 'after' , mesh_new)
    call self%reallocate_timeframe( self%timeframe_interp, 'interp', mesh_new)

      ! ! Elevation-induced SMB change
    call self%reallocate_field( mesh_new, 'delta_z'  , self%delta_z)
    call self%reallocate_field( mesh_new, 'delta_SMB', self%delta_SMB)

      ! ! Monthly SMB
    call self%reallocate_field( mesh_new, 'SMB_monthly'  , self%SMB_monthly)

    ! Re-initialise and run the SMB model to fill in all the fields
    ! =============================================================

    call self%initialise_SMB_model_ISMIP7( mesh_new, refgeo_init, refgeo_PD, region_name)
    call self%run_SMB_model_ISMIP7( mesh_new, ice, time)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_ISMIP7

  subroutine reallocate_timeframe( self, timeframe, name_postfix, mesh_new)

    ! In/output variables
    class(type_SMB_model_ISMIP7),    intent(inout) :: self
    type(type_SMB_ISMIP7_timeframe), intent(inout) :: timeframe
    character(len=*),                intent(in   ) :: name_postfix
    type(type_mesh),                 intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'reallocate_timeframe'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%SMB_ISMIP7_choice_SMB_baseline)
    case default
      call crash('invalid SMB_ISMIP7_choice_SMB_baseline "' // trim( C%SMB_ISMIP7_choice_SMB_baseline) // '"')
    case ('yearly')

      call self%reallocate_field( mesh_new, 'acabf_'         // trim( name_postfix), timeframe%acabf)
      call self%reallocate_field( mesh_new, 'dacabfdz_'      // trim( name_postfix), timeframe%dacabfdz)

    case ('fixed')

      call self%reallocate_field( mesh_new, 'acabf_anomaly_' // trim( name_postfix), timeframe%acabf_anomaly)
      call self%reallocate_field( mesh_new, 'dacabfdz_'      // trim( name_postfix), timeframe%dacabfdz)

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine reallocate_timeframe

end module ISMIP7_SMB
