MODULE laddie_hydrology

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: crash, init_routine, finalise_routine
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use mesh_utilities, only: find_containing_vertex
  use laddie_forcing_types, only: type_laddie_forcing, type_transect_SGD
  use transects_main
  use UPSY_main, only: UPSY

  implicit none

  private

  public :: initialise_transects_SGD

contains

  subroutine initialise_transects_SGD( mesh, forcing)

    ! In/output variables
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_laddie_forcing),       intent(inout) :: forcing

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'initialise_transects_SGD'
    character(len=1024)                            :: transects_str
    character(len=1024), dimension(:), allocatable :: transect_strs
    integer                                        :: it

    ! Add routine to path
    call init_routine( routine_name)

    transects_str = C%transects_SGD

    ! If empty, finalise routine
    if (transects_str == '') then
      allocate( forcing%transects( 0))
      call finalise_routine( routine_name)
      return
    end if

    ! Separate by double vertical bars
    call UPSY%stru%separate_strings_by_double_vertical_bars( transects_str, transect_strs)

    ! Allocate the right size
    allocate( forcing%transects( size(transect_strs)))

    ! Initialise each transact
    do it = 1, size(forcing%transects)
      call initialise_transect_SGD( mesh, forcing%transects(it), transect_strs(it))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_transects_SGD


  subroutine initialise_transect_SGD( mesh, transect, transect_str)

    ! In/output variables
    type(type_mesh),         intent(in   ) :: mesh
    type(type_transect_SGD), intent(  out) :: transect
    character(len=*),        intent(in   ) :: transect_str

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'initialise_transect_SGD'
    character(len=1024)                   :: source, name, filename
    real(dp), dimension(:,:), allocatable :: waypoints
    real(dp), dimension(2) :: p
    real(dp)                              :: flux_strength
    integer                               :: index_point, vi

    ! Add routine to path
    call init_routine( routine_name)

    call parse_transect_str_SGD( transect_str, source, name, filename, flux_strength)

    transect%name            = name
    transect%flux_strength   = flux_strength

    if (par%primary) write(0,*) '  Initialising SGD transect ', &
      UPSY%stru%colour_string( trim( transect%name),'light blue'), '...'

    select case (source)
    case default
      call crash('invalid transect source "' // trim( source) // '"')
    case ('hardcoded')
      ! Not installed currently
    case ('read_from_file')
      call initialise_transect_waypoints_from_file( filename, waypoints)
    end select

    ! Create transect vertices on 500m resolution to ensure the transect is fully connected
    call calc_transect_vertices_from_waypoints( transect, waypoints, 100._dp)

    index_point = 1
    ! Find vertex in mesh that contains vertex in transect
    do vi = 1, transect%nV
       p = transect%V(vi,:)
       call find_containing_vertex( mesh, p, index_point)

       ! Save index_points of mesh to transect type
       transect%index_point( vi) = index_point

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_transect_SGD

  subroutine parse_transect_str_SGD( transect_str, source, name, filename, flux_strength)

    ! The transect should be specified by a name (either one of the hard-coded options,
    ! or a direction to an external file), and a flux strength, e.g.:
    !
    !   file:transect_channel1.txt,SF=10
    !
    ! NOTE: This currently only works with file input (did not test it for hardcoded transects).

    ! In/output variables:
    character(len=*), intent(in   ) :: transect_str
    character(len=*), intent(  out) :: source, name, filename
    real(dp),         intent(  out) :: flux_strength

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'parse_transect_str_SGD'
    integer                        :: i

    ! Add routine to path
    call init_routine( routine_name)

    i = index( transect_str,',SF=')
    if (i==0) call crash('invalid transect string "' // trim(transect_str) // '" - could not find resolution specification!')

    name = transect_str( 1:i-1)
    read( transect_str( i+4:len_trim(transect_str)),*) flux_strength

    ! Separate source and name (and filename)
    i = index(name,'file:')
    if (i==0) then
      source = 'hardcoded'
    elseif (i==1) then
      source = 'read_from_file'
      filename = name(6:len_trim(name))
      i = index( filename,'/',back=.true.)
      name = filename(i+1:len_trim(filename)-4)
    else
      call crash('invalid transect string "' // trim(transect_str))
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine parse_transect_str_SGD

END MODULE laddie_hydrology

