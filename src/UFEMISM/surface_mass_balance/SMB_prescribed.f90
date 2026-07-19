module SMB_prescribed

  use parameters, only: ice_density, freshwater_density
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use SMB_model_basic, only: atype_SMB_model, type_SMB_model_context_remap
  use mpi_basic, only: par
  use netcdf_io_main, only: read_field_from_file_2D
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use climate_model_types, only: type_climate_model
  use grid_types, only: type_grid

  implicit none

  private

  public :: type_SMB_model_prescribed

  type, extends(atype_SMB_model) :: type_SMB_model_prescribed
    !< Variables and functions that are specific to the prescribed SMB model

    contains

      procedure, public :: allocate   => SMB_model_prescribed_allocate
      procedure, public :: deallocate => SMB_model_prescribed_deallocate
      procedure, public :: initialise => SMB_model_prescribed_initialise
      procedure, public :: run        => SMB_model_prescribed_run
      procedure, public :: remap_SMB_model      => remap_SMB_model_prescribed_abs

      procedure, private :: initialise_SMB_model_prescribed_notime

  end type type_SMB_model_prescribed

contains

  subroutine SMB_model_prescribed_allocate( self, region_name, mesh)

    ! In/output variables:
    class(type_SMB_model_prescribed), intent(inout) :: self
    character(len=*),                 intent(in   ) :: region_name
    type(type_mesh), target,          intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_prescribed_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all SMB models
    call self%allocate_SMB_model( 'SMB_prescribed', region_name, mesh)

    ! Allocate all the stuff that is specific to the prescribed SMB model

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine SMB_model_prescribed_allocate

  subroutine SMB_model_prescribed_deallocate( self)

    ! In/output variables:
    class(type_SMB_model_prescribed), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_prescribed_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is common to all SMB models
    call self%deallocate_SMB_model()

    ! Deallocate all the stuff that is specific to SMB model prescribed


    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine SMB_model_prescribed_deallocate

  subroutine SMB_model_prescribed_initialise( self, ice, refgeo_init, refgeo_PD)

    ! In/output variables
    class(type_SMB_model_prescribed), intent(inout) :: self
    type(type_ice_model),             intent(in   ) :: ice
    type(type_reference_geometry),    intent(in   ) :: refgeo_init
    type(type_reference_geometry),    intent(in   ) :: refgeo_PD

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_prescribed_initialise'
    character(:), allocatable   :: choice_SMB_prescribed

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise all the stuff that is common to all SMB models
    call self%initialise_SMB_model()

    ! Initialise all the stuff that is specific to SMB model prescribed

    ! Determine the type of prescribed SMB forcing for this region
    select case (self%region_name())
    case default
      call crash('unknown region_name "' // trim( self%region_name()) // '"!')
    case ('NAM')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_NAM)
    case ('EAS')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_EAS)
    case ('GRL')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_GRL)
    case ('ANT')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_ANT)
    end select

    ! Initialised the chosen type of prescribed SMB forcing
    select case (choice_SMB_prescribed)
    case default
      call crash('unknown choice_SMB_prescribed "' // trim( choice_SMB_prescribed) // '"!')
    case ('SMB_no_time')
      call self%initialise_SMB_model_prescribed_notime()
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine SMB_model_prescribed_initialise

  subroutine initialise_SMB_model_prescribed_notime( self)
    ! Prescribe SMB from a file without a time dimension

    ! In/output variables
    class(type_SMB_model_prescribed), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_SMB_model_prescribed_notime'
    character(:), allocatable   :: filename_SMB_prescribed
    real(dp)                    :: timeframe_SMB_prescribed

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine filename for this model region
    select case (self%region_name())
    case default
      call crash('unknown region_name "' // trim( self%region_name()) // '"!')
    case ('NAM')
      filename_SMB_prescribed  = trim( C%filename_SMB_prescribed_NAM)
      timeframe_SMB_prescribed = C%timeframe_SMB_prescribed_NAM
    case ('EAS')
      filename_SMB_prescribed  = trim( C%filename_SMB_prescribed_EAS)
      timeframe_SMB_prescribed = C%timeframe_SMB_prescribed_EAS
    case ('GRL')
      filename_SMB_prescribed  = trim( C%filename_SMB_prescribed_GRL)
      timeframe_SMB_prescribed = C%timeframe_SMB_prescribed_GRL
    case ('ANT')
      filename_SMB_prescribed  = trim( C%filename_SMB_prescribed_ANT)
      timeframe_SMB_prescribed = C%timeframe_SMB_prescribed_ANT
    end select

    ! Print to terminal
    if (par%primary)  write(*,"(A)") '   Initialising SMB from file "' // &
      UPSY%stru%colour_string( trim( filename_SMB_prescribed),'light blue') // '"...'

    ! Read SMB from file
    if (timeframe_SMB_prescribed == 1E9_dp) then
      ! Assume the file has no time dimension
      call read_field_from_file_2D( filename_SMB_prescribed, &
        'SMB||surface_mass_balance||', self%mesh, C%output_dir, self%SMB)
    else
      ! Assume the file has a time dimension, and read the specified timeframe
      call read_field_from_file_2D( filename_SMB_prescribed, &
        'SMB||surface_mass_balance||', self%mesh, C%output_dir, self%SMB, &
        time_to_read = timeframe_SMB_prescribed)
    end if

    ! Convert from [m.w.e. yr^-1] to [m.i.e. yr^-1]
    self%SMB( self%mesh%vi1:self%mesh%vi2) = self%SMB( self%mesh%vi1:self%mesh%vi2) &
      * freshwater_density / ice_density

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_prescribed_notime

  subroutine SMB_model_prescribed_run( self, time, ice, climate, grid_smooth)

    ! In/output variables:
    class(type_SMB_model_prescribed), intent(inout) :: self
    real(dp),                         intent(in   ) :: time
    type(type_ice_model),             intent(in   ) :: ice
    type(type_climate_model),         intent(inout) :: climate
    type(type_grid),                  intent(in   ) :: grid_smooth

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_prescribed_run'
    logical                     :: do_run_SMB_model
    character(:), allocatable   :: choice_SMB_prescribed

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is common to all SMB models
    call self%run_SMB_model( time, do_run_SMB_model)
    if (.not. do_run_SMB_model) then
      call finalise_routine( routine_name)
      return
    end if

    ! Run all the stuff that is specific to SMB model idealised

    ! Determine the type of prescribed SMB forcing for this region
    select case (self%region_name())
    case default
      call crash('unknown region_name "' // trim( self%region_name()) // '"!')
    case ('NAM')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_NAM)
    case ('EAS')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_EAS)
    case ('GRL')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_GRL)
    case ('ANT')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_ANT)
    end select

    ! Initialised the chosen type of prescribed SMB forcing
    select case (choice_SMB_prescribed)
    case default
      call crash('unknown choice_SMB_prescribed "' // trim( choice_SMB_prescribed) // '"!')
    case ('SMB_no_time')
      ! No need to do anything
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine SMB_model_prescribed_run

  subroutine remap_SMB_model_prescribed_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_prescribed),           intent(inout) :: self
    type(type_SMB_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_SMB_model_prescribed_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Re-initialise to read and remap the SMB from the input file again
    call self%initialise( context%ice, context%refgeo_init, context%refgeo_PD)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_prescribed_abs

end module SMB_prescribed
