module SMB_prescribed

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: crash, warning, happy, init_routine, &
    finalise_routine, colour_string
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use climate_model_types, only: type_climate_model
  use netcdf_io_main
  use SMB_basic, only: atype_SMB_model, type_SMB_model_context_allocate, type_SMB_model_context_initialise, &
    type_SMB_model_context_run, type_SMB_model_context_remap

  implicit none

  private

  public :: type_SMB_model_prescribed

  type, extends(atype_SMB_model) :: type_SMB_model_prescribed

    contains

      procedure, public :: allocate_SMB_model   => allocate_SMB_model_prescribed
      procedure, public :: initialise_SMB_model => initialise_SMB_model_prescribed
      procedure, public :: run_SMB_model        => run_SMB_model_prescribed
      procedure, public :: remap_SMB_model      => remap_SMB_model_prescribed

      procedure, private :: run_SMB_model_prescribed_notime
      procedure, private :: fill_initialised_fields

  end type type_SMB_model_prescribed

contains

  subroutine allocate_SMB_model_prescribed( self, context)

    ! In/output variables:
    class(type_SMB_model_prescribed),      intent(inout) :: self
    type(type_SMB_model_context_allocate), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SMB_model_prescribed'

    ! Add routine to path
    call init_routine( routine_name)

    call self%set_name('SMB_model_prescribed')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_prescribed

  subroutine initialise_SMB_model_prescribed( self, context)

    ! In/output variables:
    class(type_SMB_model_prescribed),        intent(inout) :: self
    type(type_SMB_model_context_initialise), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_prescribed'

    ! Add routine to path
    call init_routine( routine_name)

    call self%fill_initialised_fields( self%mesh, context%region_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_prescribed

  subroutine run_SMB_model_prescribed( self, context)

    ! In/output variables:
    class(type_SMB_model_prescribed), intent(inout) :: self
    type(type_SMB_model_context_run), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_prescribed'
    character(:), allocatable      :: choice_SMB_prescribed

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine the type of prescribed SMB forcing for this region
    select case (context%region_name)
    case default
      call crash('unknown region_name "' // trim( context%region_name) // '"!')
    case ('NAM')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_NAM)
    case ('EAS')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_EAS)
    case ('GRL')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_GRL)
    case ('ANT')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_ANT)
    end select

    ! Run the chosen type of prescribed SMB forcing
    select case (choice_SMB_prescribed)
    case default
      call crash('unknown choice_SMB_prescribed "' // trim( choice_SMB_prescribed) // '"!')
    case ('SMB_no_time')
      ! SMB only, no time
      call self%run_SMB_model_prescribed_notime( self%mesh)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_prescribed

  subroutine remap_SMB_model_prescribed( self, context)

    ! In/output variables:
    class(type_SMB_model_prescribed),   intent(inout) :: self
    type(type_SMB_model_context_remap), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_SMB_model_prescribed'

    ! Add routine to path
    call init_routine( routine_name)

    call self%fill_initialised_fields( self%mesh, context%region_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_prescribed

  ! == SMB only, no time
  ! ====================

  subroutine run_SMB_model_prescribed_notime( self, mesh)

    ! In/output variables:
    class(type_SMB_model_prescribed), intent(inout) :: self
    type(type_mesh),                  intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_prescribed_notime'

    ! Add routine to path
    call init_routine( routine_name)

    ! No need to do anything, as the SMB was already read during initialisation

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_prescribed_notime

  subroutine fill_initialised_fields( self, mesh, region_name)

    ! In- and output variables
    class(type_SMB_model_prescribed), intent(inout) :: self
    type(type_mesh),                  intent(in   ) :: mesh
    character(len=3),                 intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'fill_initialised_fields'
    character(:), allocatable      :: choice_SMB_prescribed
    character(len=1024)            :: filename_SMB_prescribed
    real(dp)                       :: timeframe_SMB_prescribed

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine the type of prescribed SMB forcing for this region
    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name) // '"!')
    case ('NAM')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_NAM)
    case ('EAS')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_EAS)
    case ('GRL')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_GRL)
    case ('ANT')
      choice_SMB_prescribed  = trim( C%choice_SMB_prescribed_ANT)
    end select

    ! Safety
    select case (choice_SMB_prescribed)
    case default
      call crash('unknown choice_SMB_prescribed "' // trim( choice_SMB_prescribed) // '"!')
    case ('SMB_no_time')
    end select

    ! Determine filename for this model region
    select case (region_name)
    case default
      call crash('unknown region_name "' // TRIM( region_name) // '"!')
    case ('NAM')
      filename_SMB_prescribed  = C%filename_SMB_prescribed_NAM
      timeframe_SMB_prescribed = C%timeframe_SMB_prescribed_NAM
    case ('EAS')
      filename_SMB_prescribed  = C%filename_SMB_prescribed_EAS
      timeframe_SMB_prescribed = C%timeframe_SMB_prescribed_EAS
    case ('GRL')
      filename_SMB_prescribed  = C%filename_SMB_prescribed_GRL
      timeframe_SMB_prescribed = C%timeframe_SMB_prescribed_GRL
    case ('ANT')
      filename_SMB_prescribed  = C%filename_SMB_prescribed_ANT
      timeframe_SMB_prescribed = C%timeframe_SMB_prescribed_ANT
    end select

    ! Print to terminal
    if (par%primary)  write(*,"(A)") '   Initialising SMB from file "' // &
      colour_string( trim( filename_SMB_prescribed),'light blue') // '"...'

    ! Read SMB from file
    if (timeframe_SMB_prescribed == 1E9_dp) then
      ! Assume the file has no time dimension
      call read_field_from_file_2D( filename_SMB_prescribed, &
        'SMB||surface_mass_balance||', mesh, C%output_dir, self%SMB)
    else
      ! Assume the file has a time dimension, and read the specified timeframe
      call read_field_from_file_2D( filename_SMB_prescribed, &
        'SMB||surface_mass_balance||', mesh, C%output_dir, self%SMB, time_to_read = timeframe_SMB_prescribed)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine fill_initialised_fields

end module SMB_prescribed
