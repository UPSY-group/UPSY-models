module SMB_ISMIP7

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
  use netcdf_io_main, only: read_field_from_file_2D, read_field_from_file_2D_monthly
  use basic_model_utilities, only: list_files_in_folder

  implicit none

  private

  public :: type_SMB_model_ISMIP7

  type, extends(atype_SMB_model) :: type_SMB_model_ISMIP7

      ! Main data fields
      real(dp), dimension(:), contiguous, pointer :: SMB_baseline => null()   !< Baseline annual mean SMB   [m.i.e./yr]
      real(dp), dimension(:), contiguous, pointer :: Hs_baseline  => null()   !< Baseline surface elevation [m w.r.t. PD sea level]
      type(MPI_WIN) :: wSMB_baseline, wHs_baseline

      ! List of forcing files and their timestamps
      character(len=1024), dimension(:), allocatable :: filenames_acabf
      character(len=1024), dimension(:), allocatable :: filenames_acabf_anomaly
      character(len=1024), dimension(:), allocatable :: filenames_dacabfdz
      real(dp),            dimension(:), allocatable :: timestamps_acabf
      real(dp),            dimension(:), allocatable :: timestamps_acabf_anomaly
      real(dp),            dimension(:), allocatable :: timestamps_dacabfdz

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

      procedure, private :: initialise_SMB_baseline_fixed
      procedure, private :: initialise_Hs_baseline
      procedure, private :: initialise_lists_of_files_and_timestamps
      procedure, private :: initialise_list_of_files_and_timestamps

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
    call self%initialise_SMB_model_ISMIP7( self%mesh, context%ice, context%refgeo_init, context%refgeo_PD, self%region_name())

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
    call self%run_SMB_model_ISMIP7( self%mesh, context%ice, self%region_name())

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

    call self%remap_SMB_model_ISMIP7( context%mesh_new)

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

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_ISMIP7

  subroutine initialise_SMB_model_ISMIP7( self, mesh, ice, refgeo_init, refgeo_PD, region_name)

    ! In/output variables
    class(type_SMB_model_ISMIP7),  intent(inout) :: self
    type(type_mesh),               intent(in   ) :: mesh
    type(type_ice_model),          intent(in   ) :: ice
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
      call self%initialise_SMB_baseline_fixed( mesh, ice, region_name)
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

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_ISMIP7

  subroutine initialise_SMB_baseline_fixed( self, mesh, ice, region_name)

    ! In/output variables
    class(type_SMB_model_ISMIP7), intent(inout) :: self
    type(type_mesh),              intent(in   ) :: mesh
    type(type_ice_model),         intent(in   ) :: ice
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

    call self%initialise_list_of_files_and_timestamps( self%filenames_acabf        , self%timestamps_acabf        , 'acabf')
    call self%initialise_list_of_files_and_timestamps( self%filenames_acabf_anomaly, self%timestamps_acabf_anomaly, 'acabf-anomaly')
    call self%initialise_list_of_files_and_timestamps( self%filenames_dacabfdz     , self%timestamps_dacabfdz     , 'dacabfdz')

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
    character(len=1024), dimension(:), allocatable :: list_of_filenames
    integer                                        :: i
    real(dp)                                       :: year

    ! Add routine to path
    call init_routine( routine_name)

    if (allocated( filenames )) deallocate( filenames)
    if (allocated( timestamps)) deallocate( timestamps)

    ! Construct foldernames
    foldername = trim( C%SMB_ISMIP7_forcing_foldername) // '/' // var_name // '/' // trim( C%SMB_ISMIP7_forcing_version)

    call list_files_in_folder( foldername, list_of_filenames)
    call take_only_valid_netcdf_files_from_list( list_of_filenames, filenames, var_name)
    if (size( filenames,1) == 0) call crash('could not find any valid NetCDF files in directory "' // trim( foldername) // '"')

    allocate( timestamps( size( filenames,1)))
    do i = 1, size( filenames,1)
      call read_year_from_netcdf_filename( filenames( i), year)
      ! Add half a year, so that we have the timestamp at the middle of the year rather than the start
      timestamps( i) = year + 0.5_dp
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_list_of_files_and_timestamps

  subroutine take_only_valid_netcdf_files_from_list( filenames_all, filenames, var_name)

    ! In/output variables
    character(len=1024), dimension(:),              intent(in   ) :: filenames_all
    character(len=1024), dimension(:), allocatable, intent(inout) :: filenames
    character(len=*),                               intent(in   ) :: var_name

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'take_only_valid_netcdf_files_from_list'
    character(len=:), allocatable :: filename
    integer                       :: i, n_valid_netcdf_files, j

    ! Add routine to path
    call init_routine( routine_name)

    n_valid_netcdf_files = 0
    do i = 1, size( filenames_all)
      filename = filenames_all( i)
      if (UPSY%stru%startswith( trim( filename), trim( var_name), case_sensitive = .false.)) then
        n_valid_netcdf_files = n_valid_netcdf_files + 1
      end if
    end do

    allocate( filenames( n_valid_netcdf_files))

    j = 0
    do i = 1, size( filenames_all)
      filename = filenames_all( i)
      if (UPSY%stru%startswith( trim( filename), trim( var_name), case_sensitive = .false.)) then
        j = j + 1
        filenames( j) = filename
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine take_only_valid_netcdf_files_from_list

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

  subroutine run_SMB_model_ISMIP7( self, mesh, ice, region_name)

    ! In/output variables:
    class(type_SMB_model_ISMIP7), intent(inout) :: self
    type(type_mesh),              intent(in   ) :: mesh
    type(type_ice_model),         intent(in   ) :: ice
    character(len=3),             intent(in   ) :: region_name

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_SMB_model_ISMIP7'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! DENK DROM
    call self%write_to_restart_file( C%output_dir)
    call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_ISMIP7

  subroutine remap_SMB_model_ISMIP7( self, mesh_new)

    ! In/output variables
    class(type_SMB_model_ISMIP7), intent(inout) :: self
    type(type_mesh),              intent(in   ) :: mesh_new

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_SMB_model_ISMIP7'

    ! Add routine to path
    call init_routine( routine_name)

    call self%remap_field( mesh_new, 'SMB_baseline', self%SMB_baseline)
    call self%remap_field( mesh_new, 'Hs_baseline' , self%Hs_baseline )

    call crash('whoopsiedaisy')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_ISMIP7

end module SMB_ISMIP7