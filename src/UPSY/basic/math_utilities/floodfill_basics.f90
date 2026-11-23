module floodfill_basics

  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash

  implicit none

  private

  public :: type_floodfill

  type :: type_floodfill

      integer                              :: n               !              Total number of nodes
      integer                              :: nC_max          !              Maximum number of connections per node
      integer, dimension(:  ), allocatable :: nC              ! [n-by-1]     Number of connections per node
      integer, dimension(:,:), allocatable :: C               ! [n-by-nCmax] Indices of connected nodes
      integer, dimension(:  ), allocatable :: map             ! [n-by-1]     Flood-fill map
      integer, dimension(:  ), allocatable :: stack_front     ! [n-by-1]     Stack of elements at the front
      integer                              :: stack_front_n   !              Number of elements in that stack
      integer, dimension(:  ), allocatable :: stack_checked   ! [n-by-1]     Stack of elements that have been checked
      integer                              :: stack_checked_n !              Number of elements in that stack

    contains

      procedure, public :: create   => create_floodfill
      procedure, public :: initialise_floodfill_at_single_node
      procedure, public :: take_last_element_from_front_stack
      procedure, public :: add_unchecked_neighbours_to_front_stack
      procedure, public :: cleanup_map_and_stacks

      final :: destroy_floodfill

  end type type_floodfill

  integer, parameter :: map_unchecked = 0
  integer, parameter :: map_front     = 1
  integer, parameter :: map_checked   = 2

contains

  subroutine create_floodfill( self, n, nC, C)

    ! In/output variables:
    class(type_floodfill),   intent(inout) :: self
    integer,                 intent(in   ) :: n
    integer, dimension(:  ), intent(in   ) :: nC
    integer, dimension(:,:), intent(in   ) :: C

    ! Local variables:
    character(len=*), parameter :: routine_name = 'create_floodfill'

    ! Add routine to path
    call init_routine( routine_name)

    self%n      = n
    self%nC_max = maxval( nC)

    ! Safety
    if (size( nC,1) /= n) call crash('size( nC,1) /= n')
    if (size(  C,1) /= n) call crash('size( C,1) /= n')
    if (size(  C,2) < self%nC_max) call crash('size(  C,2) < nC_max')

    allocate( self%nC( n), source = nC)
    allocate( self%C( n, self%nC_max), source = C( :,1:self%nC_max))

    allocate( self%map( n), source = map_unchecked)
    allocate( self%stack_front( n))
    self%stack_front_n = 0
    allocate( self%stack_checked( n))
    self%stack_checked_n = 0

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_floodfill

  subroutine initialise_floodfill_at_single_node( self, i)

    ! In/output variables:
    class(type_floodfill),   intent(inout) :: self
    integer,                 intent(in   ) :: i

    ! Local variables:
    character(len=*), parameter :: routine_name = 'initialise_floodfill_at_single_node'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (i < 1 .or. i > self%n) call crash('invalid node index')

    call self%cleanup_map_and_stacks()

    self%stack_front_n = 1
    self%stack_front( self%stack_front_n) = i
    self%map( i) = map_front

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_floodfill_at_single_node

  subroutine take_last_element_from_front_stack( self, i)

    ! In/output variables:
    class(type_floodfill),   intent(inout) :: self
    integer,                 intent(  out) :: i

    i = self%stack_front( self%stack_front_n)
    self%stack_front_n = self%stack_front_n - 1

    ! Mark it as checked
    self%stack_checked_n = self%stack_checked_n + 1
    self%stack_checked( self%stack_checked_n) = i
    self%map( i) = map_checked

  end subroutine take_last_element_from_front_stack

  subroutine add_unchecked_neighbours_to_front_stack( self, i)

    ! In/output variables:
    class(type_floodfill),   intent(inout) :: self
    integer,                 intent(in   ) :: i

    ! Local variables:
    integer :: ci, j

    do ci = 1, self%nC( i)
      j = self%C( i, ci)
      if (self%map( j) == 0) then
        self%map( j) = map_front
        self%stack_front_n = self%stack_front_n + 1
        self%stack_front( self%stack_front_n) = j
      end if
    end do

  end subroutine add_unchecked_neighbours_to_front_stack

  subroutine cleanup_map_and_stacks( self)

    ! In/output variables:
    class(type_floodfill),   intent(inout) :: self

    ! Local variables:
    integer :: ii, i

    ! Front stack
    do ii = 1, self%stack_front_n
      i = self%stack_front( ii)
      self%map( i) = map_unchecked
    end do
    self%stack_front_n = 0

    ! Checked stack
    do ii = 1, self%stack_checked_n
      i = self%stack_checked( ii)
      self%map( i) = map_unchecked
    end do
    self%stack_checked_n = 0

  end subroutine cleanup_map_and_stacks

  subroutine destroy_floodfill( self)

    ! In/output variables:
    type(type_floodfill), intent(inout) :: self

    self%n = 0
    self%nC_max = 0
    self%stack_front_n = 0
    self%stack_checked_n = 0

    if (allocated( self%nC)) deallocate( self%nC)
    if (allocated( self%C)) deallocate( self%C)
    if (allocated( self%map)) deallocate( self%map)
    if (allocated( self%stack_front)) deallocate( self%stack_front)
    if (allocated( self%stack_checked)) deallocate( self%stack_checked)

  end subroutine destroy_floodfill

end module floodfill_basics
