module remapping_basic

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use basic_grid_mod, only: atype_grid
  use Arakawa_grid_mod, only: type_Arakawa_grid
  use remapping_context_mod, only: type_remapping_context
  use remapping_calc_CSR_operators, only: calc_CSR_remapping_operator
  use remapping_methods, only: type_remapping_method, remapping_method
  use ugrid_basic, only: type_ugrid
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_1D_wrapper, &
    multiply_CSR_matrix_with_vector_2D_wrapper

  implicit none

  private

  public :: remap

  ! Interface bundling the 2-D and 3-D remapping procedures
  interface remap
    procedure :: remap_2D
    procedure :: remap_3D
  end interface remap

  ! Type of a remapping object
  type :: type_map

      type(type_remapping_context) :: ct
      type(type_CSR_matrix_dp)     :: M

    contains

      procedure, public :: create => create_map
      generic,   public :: apply => &
        apply_map_2D, &
        apply_map_3D
      procedure, public :: apply_map_2D
      procedure, public :: apply_map_3D

  end type type_map

  ! Collection of all the maps that have already been calculated
  type(type_map), dimension(:), allocatable, target :: Atlas

contains

  subroutine remap_2D( grid_src, Arakawa_grid_src, d_src, grid_dst, Arakawa_grid_dst, d_dst, method)

    ! In/output variables:
    class(atype_grid), target,   intent(in   ) :: grid_src
    type(type_Arakawa_grid),     intent(in   ) :: Arakawa_grid_src
    real(dp), dimension(:),      intent(in   ) :: d_src
    class(atype_grid), target,   intent(in   ) :: grid_dst
    type(type_Arakawa_grid),     intent(in   ) :: Arakawa_grid_dst
    real(dp), dimension(:),      intent(  out) :: d_dst
    type(type_remapping_method), intent(in   ) :: method

    ! Local variables:
    character(len=*), parameter  :: routine_name = 'remap_2D'
    type(type_remapping_context) :: ct
    type(type_map), pointer      :: map

    ! Add routine to path
    call init_routine( routine_name)

    call ct%create( grid_src, Arakawa_grid_src, grid_dst, Arakawa_grid_dst, method)
    call find_in_or_add_map_to_Atlas( ct, map)
    call map%apply( d_src, d_dst)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_2D

  subroutine remap_3D( grid_src, Arakawa_grid_src, d_src, grid_dst, Arakawa_grid_dst, d_dst, method)

    ! In/output variables:
    class(atype_grid), target,   intent(in   ) :: grid_src
    type(type_Arakawa_grid),     intent(in   ) :: Arakawa_grid_src
    real(dp), dimension(:,:),    intent(in   ) :: d_src
    class(atype_grid), target,   intent(in   ) :: grid_dst
    type(type_Arakawa_grid),     intent(in   ) :: Arakawa_grid_dst
    real(dp), dimension(:,:),    intent(  out) :: d_dst
    type(type_remapping_method), intent(in   ) :: method

    ! Local variables:
    character(len=*), parameter  :: routine_name = 'remap_3D'
    type(type_remapping_context) :: ct
    type(type_map), pointer      :: map

    ! Add routine to path
    call init_routine( routine_name)

    call ct%create( grid_src, Arakawa_grid_src, grid_dst, Arakawa_grid_dst, method)
    call find_in_or_add_map_to_Atlas( ct, map)
    call map%apply( d_src, d_dst)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_3D

  subroutine find_in_or_add_map_to_Atlas( ct, map)
    ! Return a pointer to an existing or newly created mapping object in the Atlas

    ! In/output variables:
    type(type_remapping_context), intent(in   ) :: ct
    type(type_map), pointer,      intent(  out) :: map

    ! Local variables:
    character(len=*), parameter               :: routine_name = 'find_in_or_add_map_to_Atlas'
    integer                                   :: mi
    logical                                   :: found_it
    type(type_map), dimension(:), allocatable :: Atlas_old

    ! Add routine to path
    call init_routine( routine_name)

    if (.not. allocated( Atlas)) then
      ! Allocate room for a single map in the Atlas, then create that map
      allocate( Atlas( 1))
      map => Atlas( 1)
      call map%create( ct)

    else

      ! See if the correct map already exists in the Atlas
      found_it = .false.
      do mi = 1, size( Atlas)
        if (Atlas(mi)%ct == ct) then
          found_it = .true.
          map => Atlas( mi)
        end if
      end do

      ! If not, extend the Atlas and add a new mapping object to it
      if (.not. found_it) then
        allocate( Atlas_old( size( Atlas)), source = Atlas)
        deallocate( Atlas)
        allocate( Atlas( size( Atlas_old)+1))
        Atlas( 1: size( Atlas_old)) = Atlas_old
        deallocate( Atlas_old)
        map => Atlas( size( Atlas))
        call map%create( ct)
      end if

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine find_in_or_add_map_to_Atlas

  subroutine create_map( self, ct)
    ! Create a new mapping object matching the specified grids and remapping method

    ! In/output variables:
    class(type_map),              intent(inout) :: self
    type(type_remapping_context), intent(in   ) :: ct

    ! Local variables:
    character(len=*), parameter :: routine_name = 'create_map'
    type(type_ugrid)            :: ugrid_src
    type(type_ugrid)            :: ugrid_dst

    ! Add routine to path
    call init_routine( routine_name)

    self%ct%grid_src         => ct%grid_src
    self%ct%Arakawa_grid_src =  ct%Arakawa_grid_src
    self%ct%grid_dst         => ct%grid_dst
    self%ct%Arakawa_grid_dst =  ct%Arakawa_grid_dst
    self%ct%method           =  ct%method

    ! Create ugrids from source and destination grids
    call ugrid_src%build_from_grid( self%ct%grid_src, self%ct%Arakawa_grid_src)
    call ugrid_dst%build_from_grid( self%ct%grid_dst, self%ct%Arakawa_grid_dst)

    ! Calculate remapping operator between the two ugrids
    call calc_CSR_remapping_operator( ugrid_src, ugrid_dst, self%ct%method, self%M)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_map

  subroutine apply_map_2D( self, d_src, d_dst)
    ! Apply a remapping object as a sparse matrix multiplication

    ! In/output variables:
    class(type_map),        intent(in   ) :: self
    real(dp), dimension(:), intent(in   ) :: d_src
    real(dp), dimension(:), intent(  out) :: d_dst

    ! Local variables:
    character(len=*), parameter :: routine_name = 'apply_map_2D'

    ! Add routine to path
    call init_routine( routine_name)

    call multiply_CSR_matrix_with_vector_1D_wrapper( self%M, &
      self%ct%pai_src, d_src, self%ct%pai_dst, d_dst)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_2D

  subroutine apply_map_3D( self, d_src, d_dst)
    ! Apply a remapping object as a sparse matrix multiplication

    ! In/output variables:
    class(type_map),          intent(in   ) :: self
    real(dp), dimension(:,:), intent(in   ) :: d_src
    real(dp), dimension(:,:), intent(  out) :: d_dst

    ! Local variables:
    character(len=*), parameter :: routine_name = 'apply_map_3D'

    ! Add routine to path
    call init_routine( routine_name)

    call multiply_CSR_matrix_with_vector_2D_wrapper( self%M, &
      self%ct%pai_src, d_src, self%ct%pai_dst, d_dst)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_3D

end module remapping_basic