module SMB_prescribed

  use parameters, only: pi
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use SMB_model_basic, only: atype_SMB_model, type_SMB_model_context_allocate, &
    type_SMB_model_context_initialise, type_SMB_model_context_run, &
    type_SMB_model_context_remap
  use mpi_basic, only: par
  use netcdf_io_main, only: read_field_from_file_2D

  implicit none

  private

  public :: type_SMB_model_prescribed

  type, extends(atype_SMB_model) :: type_SMB_model_prescribed

    contains

      procedure, public :: allocate_SMB_model   => allocate_SMB_model_prescribed_abs
      procedure, public :: deallocate_SMB_model => deallocate_SMB_model_prescribed_abs
      procedure, public :: initialise_SMB_model => initialise_SMB_model_prescribed_abs
      procedure, public :: run_SMB_model        => run_SMB_model_prescribed_abs
      procedure, public :: remap_SMB_model      => remap_SMB_model_prescribed_abs

      procedure, private :: run_SMB_model_prescribed
      procedure, private :: initialise_SMB_model_prescribed
      procedure, private :: initialise_SMB_model_prescribed_notime

  end type type_SMB_model_prescribed

contains

  subroutine allocate_SMB_model_prescribed_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_prescribed),              intent(inout) :: self
    type(type_SMB_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SMB_model_prescribed_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_prescribed_abs

  subroutine deallocate_SMB_model_prescribed_abs( self)

    ! In/output variables:
    class(type_SMB_model_prescribed), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'deallocate_SMB_model_prescribed_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_SMB_model_prescribed_abs

  subroutine initialise_SMB_model_prescribed_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_prescribed),                intent(inout) :: self
    type(type_SMB_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_prescribed_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call self%initialise_SMB_model_prescribed( self%mesh, self%region_name())

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_prescribed_abs

  subroutine run_SMB_model_prescribed_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_prescribed),         intent(inout) :: self
    type(type_SMB_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_prescribed_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call self%run_SMB_model_prescribed( self%region_name())

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_prescribed_abs

  subroutine remap_SMB_model_prescribed_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_prescribed),           intent(inout) :: self
    type(type_SMB_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_SMB_model_prescribed_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Re-initialise to read and remap the SMB from the input file again
    call self%initialise_SMB_model_prescribed( self%mesh, self%region_name())

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_prescribed_abs



  subroutine initialise_SMB_model_prescribed( self, mesh, region_name)

    ! In/output variables
    class(type_SMB_model_prescribed), intent(inout) :: self
    type(type_mesh),                  intent(in   ) :: mesh
    character(len=3),                 intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_prescribed'
    character(:), allocatable      :: choice_SMB_prescribed

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

    ! Initialised the chosen type of prescribed SMB forcing
    select case (choice_SMB_prescribed)
    case default
      call crash('unknown choice_SMB_prescribed "' // trim( choice_SMB_prescribed) // '"!')
    case ('SMB_no_time')
      call self%initialise_SMB_model_prescribed_notime( mesh, region_name)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_prescribed

  subroutine initialise_SMB_model_prescribed_notime( self, mesh, region_name)
    ! Prescribe SMB from a file without a time dimension

    ! In/output variables
    class(type_SMB_model_prescribed), intent(inout) :: self
    type(type_mesh),                  intent(in   ) :: mesh
    character(len=3),                 intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_prescribed_notime'
    character(:), allocatable      :: filename_SMB_prescribed
    real(dp)                       :: timeframe_SMB_prescribed

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine filename for this model region
    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name) // '"!')
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
        'SMB||surface_mass_balance||', mesh, C%output_dir, self%SMB)
    else
      ! Assume the file has a time dimension, and read the specified timeframe
      call read_field_from_file_2D( filename_SMB_prescribed, &
        'SMB||surface_mass_balance||', mesh, C%output_dir, self%SMB, &
        time_to_read = timeframe_SMB_prescribed)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_prescribed_notime

  subroutine run_SMB_model_prescribed( self, region_name)

    ! In/output variables:
    class(type_SMB_model_prescribed), intent(inout) :: self
    character(len=3),                 intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_prescribed'
    character(:), allocatable      :: choice_SMB_prescribed

    ! Add routine to call stack
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

    ! Initialised the chosen type of prescribed SMB forcing
    select case (choice_SMB_prescribed)
    case default
      call crash('unknown choice_SMB_prescribed "' // trim( choice_SMB_prescribed) // '"!')
    case ('SMB_no_time')
      ! No need to do anything
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_prescribed

end module SMB_prescribed