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

  implicit none

  private

  public :: type_SMB_model_ISMIP7

  type, extends(atype_SMB_model) :: type_SMB_model_ISMIP7

      ! Main data fields
      real(dp), dimension(:), contiguous, pointer :: SMB_baseline => null()   !< Baseline annual mean SMB   [m.i.e./yr]
      real(dp), dimension(:), contiguous, pointer :: Hs_baseline  => null()   !< Baseline surface elevation [m w.r.t. PD sea level]
      type(MPI_WIN) :: wSMB_baseline, wHs_baseline

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
      call initialise_SMB_baseline_fixed( self, mesh, ice, region_name)
    end select

    ! Initialise the baseline surface elevation
    select case (C%SMB_ISMIP7_choice_refgeo)
    case default
      call crash('invalid SMB_ISMIP7_choice_refgeo "' // trim( C%SMB_ISMIP7_choice_refgeo) // '"')
    case ('init')
      call initialise_Hs_baseline( self, mesh, refgeo_init)
    case ('PD')
      call initialise_Hs_baseline( self, mesh, refgeo_PD)
    end select

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