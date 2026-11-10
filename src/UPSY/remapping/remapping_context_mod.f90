module remapping_context_mod

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use basic_grid_mod, only: atype_grid
  use Arakawa_grid_mod, only: type_Arakawa_grid, Arakawa_grid
  use remapping_methods, only: type_remapping_method, remapping_method
  use parallel_array_info_type, only: type_par_arr_info
  use xy_grid_basic, only: type_xy_grid
  use lonlat_grid_basic, only: type_lonlat_grid
  use mesh_basic, only: type_mesh

  implicit none

  private

  public :: type_remapping_context

  type :: type_remapping_context

      class(atype_grid), pointer  :: grid_src          ! Source      grid (x/y-grid, lon/lat-grid, or mesh)
      type(type_Arakawa_grid)     :: Arakawa_grid_src  ! Source      Arakawa grid
      type(type_par_arr_info)     :: pai_src           ! Source      grid data parallelisation info
      class(atype_grid), pointer  :: grid_dst          ! Destination grid
      type(type_Arakawa_grid)     :: Arakawa_grid_dst  ! Destination Arakawa grid
      type(type_par_arr_info)     :: pai_dst           ! Destination grid data parallelisation info
      type(type_remapping_method) :: method            ! Remapping method (nearest-neighbour, linear, 1st order cons., 2nd order cons.)

    contains

      procedure, public :: create => create_remapping_context
      procedure, public :: eq     => test_remapping_context_equality
      generic           :: operator(==) => eq

  end type type_remapping_context

contains

  subroutine create_remapping_context( self, grid_src, Arakawa_grid_src, grid_dst, Arakawa_grid_dst, method)

    ! In/output variables:
    class(type_remapping_context), intent(inout) :: self
    class(atype_grid), target,     intent(in   ) :: grid_src
    type(type_Arakawa_grid),       intent(in   ) :: Arakawa_grid_src
    class(atype_grid), target,     intent(in   ) :: grid_dst
    type(type_Arakawa_grid),       intent(in   ) :: Arakawa_grid_dst
    type(type_remapping_method),   intent(in   ) :: method

    ! Local variables:
    character(len=*), parameter  :: routine_name = 'create_remapping_context'

    ! Add routine to path
    call init_routine( routine_name)

    self%grid_src         => grid_src
    self%Arakawa_grid_src =  Arakawa_grid_src
    self%grid_dst         => grid_dst
    self%Arakawa_grid_dst =  Arakawa_grid_dst
    self%method           =  method

    call set_remapping_context_pai( grid_src, Arakawa_grid_src, self%pai_src)
    call set_remapping_context_pai( grid_dst, Arakawa_grid_dst, self%pai_dst)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_remapping_context

  subroutine set_remapping_context_pai( grid, grid_Arakawa_grid, pai)

    ! In/output variables:
    class(atype_grid),       intent(in   ) :: grid
    type(type_Arakawa_grid), intent(in   ) :: grid_Arakawa_grid
    type(type_par_arr_info), intent(  out) :: pai

    ! Local variables:
    character(len=*), parameter  :: routine_name = 'set_remapping_context_pai'

    ! Add routine to path
    call init_routine( routine_name)

    select type (g => grid)
    class default
      call crash('unsupported grid type')
    class is (type_xy_grid)

      if (grid_Arakawa_grid == Arakawa_grid%a()) then
        pai = g%pai
      else
        call crash('unsupported Arakawa grid')
      end if

    class is (type_lonlat_grid)

      call crash('fixme')

    class is (type_mesh)

      if     (grid_Arakawa_grid == Arakawa_grid%a()) then
        pai = g%pai_V
      elseif (grid_Arakawa_grid == Arakawa_grid%b()) then
        pai = g%pai_Tri
      elseif (grid_Arakawa_grid == Arakawa_grid%c()) then
        pai = g%pai_E
      else
        call crash('unsupported Arakawa grid')
      end if

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine set_remapping_context_pai

  pure function test_remapping_context_equality( self, other) result( res)
    class(type_remapping_context), intent(in) :: self
    class(type_remapping_context), intent(in) :: other
    logical                                   :: res
    res = &
      self%grid_src%name    == other%grid_src%name .and. &
      self%Arakawa_grid_src == other%Arakawa_grid_src .and. &
      self%grid_dst%name    == other%grid_dst%name .and. &
      self%Arakawa_grid_dst == other%Arakawa_grid_dst .and. &
      self%method           == other%method
  end function test_remapping_context_equality

end module remapping_context_mod