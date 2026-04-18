module reallocate_dist_shared_mod

  use precisions, only: dp
  use mpi_basic, only: par, sync_node
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use deallocate_dist_shared_mod, only: deallocate_dist_shared
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: reallocate_dist_shared

  interface reallocate_dist_shared
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object
    procedure :: reallocate_dist_shared_logical_1D
    procedure :: reallocate_dist_shared_logical_2D
    procedure :: reallocate_dist_shared_logical_3D
    procedure :: reallocate_dist_shared_int_1D
    procedure :: reallocate_dist_shared_int_2D
    procedure :: reallocate_dist_shared_int_3D
    procedure :: reallocate_dist_shared_dp_1D
    procedure :: reallocate_dist_shared_dp_2D
    procedure :: reallocate_dist_shared_dp_3D
    procedure :: reallocate_dist_shared_complex_1D
    procedure :: reallocate_dist_shared_complex_2D
    procedure :: reallocate_dist_shared_complex_3D
  end interface reallocate_dist_shared

contains

  subroutine reallocate_dist_shared_logical_1D( p, win, i1_new, i2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    logical, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                  intent(inout) :: win        !< Corresponding MPI window
    integer,                        intent(in   ) :: i1_new, i2_new     !< Bounds of the memory to be allocated

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_dist_shared_logical_1D'
    integer                        :: i1_old, i2_old, n1_old, n1_new
    logical, dimension(:), pointer :: p_temp => null()
    type(MPI_WIN)                  :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    i1_old = lbound( p,1)
    i2_old = ubound( p,1)
    n1_old = size( p,1)
    n1_new = i2_new + 1 - i1_new

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, i1_old, i2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, i1_new, i2_new)

    ! Copy data there
    if (par%node_primary) p     ( i1_new: i1_new - 1 + min( n1_old, n1_new)) = &
                          p_temp( i1_old: i1_old - 1 + min( n1_old, n1_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_logical_1D

  subroutine reallocate_dist_shared_logical_2D( p, win, i1_new, i2_new, j1_new, j2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    logical, dimension(:,:), pointer, intent(inout) :: p                  !< Pointer to memory
    type(MPI_WIN),                    intent(inout) :: win                !< Corresponding MPI window
    integer,                          intent(in   ) :: i1_new, i2_new, j1_new, j2_new     !< Bounds of the memory to be allocated

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'reallocate_dist_shared_logical_2D'
    integer                          :: i1_old, i2_old, j1_old, j2_old, n1_old, n2_old, n1_new, n2_new
    logical, dimension(:,:), pointer :: p_temp => null()
    type(MPI_WIN)                    :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    i1_old = lbound( p,1)
    i2_old = ubound( p,1)
    j1_old = lbound( p,2)
    j2_old = ubound( p,2)
    n1_old = size( p,1)
    n2_old = size( p,2)
    n1_new = i2_new + 1 - i1_new
    n2_new = j2_new + 1 - j1_new

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, i1_old, i2_old, j1_old, j2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, i1_new, i2_new, j1_new, j2_new)

    ! Copy data there
    if (par%node_primary) p( i1_new: i1_new - 1 + min( n1_old, n1_new), j1_new: j1_new - 1 + min( n2_old, n2_new)) = &
      p_temp( i1_old: i1_old - 1 + min( n1_old, n1_new), j1_old: j1_old - 1 + min( n2_old, n2_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_logical_2D

  subroutine reallocate_dist_shared_logical_3D( p, win, i1_new, i2_new, j1_new, j2_new, k1_new, k2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    logical, dimension(:,:,:), pointer, intent(inout) :: p                          !< Pointer to memory
    type(MPI_WIN),                      intent(inout) :: win                        !< Corresponding MPI window
    integer,                            intent(in   ) :: i1_new, i2_new, j1_new, j2_new, k1_new, k2_new     !< Bounds of the memory to be allocated

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'reallocate_dist_shared_logical_3D'
    integer                            :: i1_old, i2_old, j1_old, j2_old, k1_old, k2_old, n1_old, n2_old, n3_old, n1_new, n2_new, n3_new
    logical, dimension(:,:,:), pointer :: p_temp => null()
    type(MPI_WIN)                      :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    i1_old = lbound( p,1)
    i2_old = ubound( p,1)
    j1_old = lbound( p,2)
    j2_old = ubound( p,2)
    k1_old = lbound( p,3)
    k2_old = ubound( p,3)
    n1_old = size( p,1)
    n2_old = size( p,2)
    n3_old = size( p,3)
    n1_new = i2_new + 1 - i1_new
    n2_new = j2_new + 1 - j1_new
    n3_new = k2_new + 1 - k1_new

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, i1_old, i2_old, j1_old, j2_old, k1_old, k2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, i1_new, i2_new, j1_new, j2_new, k1_new, k2_new)

    ! Copy data there
    if (par%node_primary) p( i1_new: i1_new - 1 + min( n1_old, n1_new), j1_new: j1_new - 1 + min( n2_old, n2_new), k1_new: k1_new - 1 + min( n3_old, n3_new)) = &
      p_temp( i1_old: i1_old - 1 + min( n1_old, n1_new), j1_old: j1_old - 1 + min( n2_old, n2_new), k1_old: k1_old - 1 + min( n3_old, n3_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_logical_3D

  subroutine reallocate_dist_shared_int_1D( p, win, i1_new, i2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    integer, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                  intent(inout) :: win        !< Corresponding MPI window
    integer,                        intent(in   ) :: i1_new, i2_new     !< Bounds of the memory to be allocated

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_dist_shared_int_1D'
    integer                        :: i1_old, i2_old, n1_old, n1_new
    integer, dimension(:), pointer :: p_temp => null()
    type(MPI_WIN)                  :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    i1_old = lbound( p,1)
    i2_old = ubound( p,1)
    n1_old = size( p,1)
    n1_new = i2_new + 1 - i1_new

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, i1_old, i2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, i1_new, i2_new)

    ! Copy data there
    if (par%node_primary) p( i1_new: i1_new - 1 + min( n1_old, n1_new)) = p_temp( i1_old: i1_old - 1 + min( n1_old, n1_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_int_1D

  subroutine reallocate_dist_shared_int_2D( p, win, i1_new, i2_new, j1_new, j2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    integer, dimension(:,:), pointer, intent(inout) :: p                  !< Pointer to memory
    type(MPI_WIN),                    intent(inout) :: win                !< Corresponding MPI window
    integer,                          intent(in   ) :: i1_new, i2_new, j1_new, j2_new     !< Bounds of the memory to be allocated

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'reallocate_dist_shared_int_2D'
    integer                          :: i1_old, i2_old, j1_old, j2_old, n1_old, n2_old, n1_new, n2_new
    integer, dimension(:,:), pointer :: p_temp => null()
    type(MPI_WIN)                    :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    i1_old = lbound( p,1)
    i2_old = ubound( p,1)
    j1_old = lbound( p,2)
    j2_old = ubound( p,2)
    n1_old = size( p,1)
    n2_old = size( p,2)
    n1_new = i2_new + 1 - i1_new
    n2_new = j2_new + 1 - j1_new

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, i1_old, i2_old, j1_old, j2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, i1_new, i2_new, j1_new, j2_new)

    ! Copy data there
    if (par%node_primary) p( i1_new: i1_new - 1 + min( n1_old, n1_new), j1_new: j1_new - 1 + min( n2_old, n2_new)) = &
      p_temp( i1_old: i1_old - 1 + min( n1_old, n1_new), j1_old: j1_old - 1 + min( n2_old, n2_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_int_2D

  subroutine reallocate_dist_shared_int_3D( p, win, i1_new, i2_new, j1_new, j2_new, k1_new, k2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    integer, dimension(:,:,:), pointer, intent(inout) :: p                          !< Pointer to memory
    type(MPI_WIN),                      intent(inout) :: win                        !< Corresponding MPI window
    integer,                            intent(in   ) :: i1_new, i2_new, j1_new, j2_new, k1_new, k2_new     !< Bounds of the memory to be allocated

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'reallocate_dist_shared_int_3D'
    integer                            :: i1_old, i2_old, j1_old, j2_old, k1_old, k2_old, n1_old, n2_old, n3_old, n1_new, n2_new, n3_new
    integer, dimension(:,:,:), pointer :: p_temp => null()
    type(MPI_WIN)                      :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    i1_old = lbound( p,1)
    i2_old = ubound( p,1)
    j1_old = lbound( p,2)
    j2_old = ubound( p,2)
    k1_old = lbound( p,3)
    k2_old = ubound( p,3)
    n1_old = size( p,1)
    n2_old = size( p,2)
    n3_old = size( p,3)
    n1_new = i2_new + 1 - i1_new
    n2_new = j2_new + 1 - j1_new
    n3_new = k2_new + 1 - k1_new

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, i1_old, i2_old, j1_old, j2_old, k1_old, k2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, i1_new, i2_new, j1_new, j2_new, k1_new, k2_new)

    ! Copy data there
    if (par%node_primary) p( i1_new: i1_new - 1 + min( n1_old, n1_new), j1_new: j1_new - 1 + min( n2_old, n2_new), k1_new: k1_new - 1 + min( n3_old, n3_new)) = &
      p_temp( i1_old: i1_old - 1 + min( n1_old, n1_new), j1_old: j1_old - 1 + min( n2_old, n2_new), k1_old: k1_old - 1 + min( n3_old, n3_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_int_3D

  subroutine reallocate_dist_shared_dp_1D( p, win, i1_new, i2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    real(dp), dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                   intent(inout) :: win        !< Corresponding MPI window
    integer,                         intent(in   ) :: i1_new, i2_new     !< Bounds of the memory to be allocated

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'reallocate_dist_shared_dp_1D'
    integer                         :: i1_old, i2_old, n1_old, n1_new
    real(dp), dimension(:), pointer :: p_temp => null()
    type(MPI_WIN)                   :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    i1_old = lbound( p,1)
    i2_old = ubound( p,1)
    n1_old = size( p,1)
    n1_new = i2_new + 1 - i1_new

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, i1_old, i2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, i1_new, i2_new)

    ! Copy data there
    if (par%node_primary) p( i1_new: i1_new - 1 + min( n1_old, n1_new)) = p_temp( i1_old: i1_old - 1 + min( n1_old, n1_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_dp_1D

  subroutine reallocate_dist_shared_dp_2D( p, win, i1_new, i2_new, j1_new, j2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    real(dp), dimension(:,:), pointer, intent(inout) :: p                  !< Pointer to memory
    type(MPI_WIN),                     intent(inout) :: win                !< Corresponding MPI window
    integer,                           intent(in   ) :: i1_new, i2_new, j1_new, j2_new     !< Bounds of the memory to be allocated

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'reallocate_dist_shared_dp_2D'
    integer                           :: i1_old, i2_old, j1_old, j2_old, n1_old, n2_old, n1_new, n2_new
    real(dp), dimension(:,:), pointer :: p_temp => null()
    type(MPI_WIN)                     :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    i1_old = lbound( p,1)
    i2_old = ubound( p,1)
    j1_old = lbound( p,2)
    j2_old = ubound( p,2)
    n1_old = size( p,1)
    n2_old = size( p,2)
    n1_new = i2_new + 1 - i1_new
    n2_new = j2_new + 1 - j1_new

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, i1_old, i2_old, j1_old, j2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, i1_new, i2_new, j1_new, j2_new)

    ! Copy data there
    if (par%node_primary) p( i1_new: i1_new - 1 + min( n1_old, n1_new), j1_new: j1_new - 1 + min( n2_old, n2_new)) = &
      p_temp( i1_old: i1_old - 1 + min( n1_old, n1_new), j1_old: j1_old - 1 + min( n2_old, n2_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_dp_2D

  subroutine reallocate_dist_shared_dp_3D( p, win, i1_new, i2_new, j1_new, j2_new, k1_new, k2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    real(dp), dimension(:,:,:), pointer, intent(inout) :: p                          !< Pointer to memory
    type(MPI_WIN),                       intent(inout) :: win                        !< Corresponding MPI window
    integer,                             intent(in   ) :: i1_new, i2_new, j1_new, j2_new, k1_new, k2_new     !< Bounds of the memory to be allocated

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'reallocate_dist_shared_dp_3D'
    integer                             :: i1_old, i2_old, j1_old, j2_old, k1_old, k2_old, n1_old, n2_old, n3_old, n1_new, n2_new, n3_new
    real(dp), dimension(:,:,:), pointer :: p_temp => null()
    type(MPI_WIN)                       :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    i1_old = lbound( p,1)
    i2_old = ubound( p,1)
    j1_old = lbound( p,2)
    j2_old = ubound( p,2)
    k1_old = lbound( p,3)
    k2_old = ubound( p,3)
    n1_old = size( p,1)
    n2_old = size( p,2)
    n3_old = size( p,3)
    n1_new = i2_new + 1 - i1_new
    n2_new = j2_new + 1 - j1_new
    n3_new = k2_new + 1 - k1_new

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, i1_old, i2_old, j1_old, j2_old, k1_old, k2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, i1_new, i2_new, j1_new, j2_new, k1_new, k2_new)

    ! Copy data there
    if (par%node_primary) p( i1_new: i1_new - 1 + min( n1_old, n1_new), j1_new: j1_new - 1 + min( n2_old, n2_new), k1_new: k1_new - 1 + min( n3_old, n3_new)) = &
      p_temp( i1_old: i1_old - 1 + min( n1_old, n1_new), j1_old: j1_old - 1 + min( n2_old, n2_new), k1_old: k1_old - 1 + min( n3_old, n3_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_dp_3D

  subroutine reallocate_dist_shared_complex_1D( p, win, i1_new, i2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    complex(dp), dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                      intent(inout) :: win        !< Corresponding MPI window
    integer,                            intent(in   ) :: i1_new, i2_new     !< Bounds of the memory to be allocated

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'reallocate_dist_shared_complex_1D'
    integer                            :: i1_old, i2_old, n1_old, n1_new
    complex(dp), dimension(:), pointer :: p_temp => null()
    type(MPI_WIN)                      :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    i1_old = lbound( p,1)
    i2_old = ubound( p,1)
    n1_old = size( p,1)
    n1_new = i2_new + 1 - i1_new

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, i1_old, i2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, i1_new, i2_new)

    ! Copy data there
    if (par%node_primary) p( i1_new: i1_new - 1 + min( n1_old, n1_new)) = p_temp( i1_old: i1_old - 1 + min( n1_old, n1_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_complex_1D

  subroutine reallocate_dist_shared_complex_2D( p, win, i1_new, i2_new, j1_new, j2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    complex(dp), dimension(:,:), pointer, intent(inout) :: p                  !< Pointer to memory
    type(MPI_WIN),                        intent(inout) :: win                !< Corresponding MPI window
    integer,                              intent(in   ) :: i1_new, i2_new, j1_new, j2_new     !< Bounds of the memory to be allocated

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'reallocate_dist_shared_complex_2D'
    integer                              :: i1_old, i2_old, j1_old, j2_old, n1_old, n2_old, n1_new, n2_new
    complex(dp), dimension(:,:), pointer :: p_temp => null()
    type(MPI_WIN)                        :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    i1_old = lbound( p,1)
    i2_old = ubound( p,1)
    j1_old = lbound( p,2)
    j2_old = ubound( p,2)
    n1_old = size( p,1)
    n2_old = size( p,2)
    n1_new = i2_new + 1 - i1_new
    n2_new = j2_new + 1 - j1_new

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, i1_old, i2_old, j1_old, j2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, i1_new, i2_new, j1_new, j2_new)

    ! Copy data there
    if (par%node_primary) p( i1_new: i1_new - 1 + min( n1_old, n1_new), j1_new: j1_new - 1 + min( n2_old, n2_new)) = &
      p_temp( i1_old: i1_old - 1 + min( n1_old, n1_new), j1_old: j1_old - 1 + min( n2_old, n2_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_complex_2D

  subroutine reallocate_dist_shared_complex_3D( p, win, i1_new, i2_new, j1_new, j2_new, k1_new, k2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    complex(dp), dimension(:,:,:), pointer, intent(inout) :: p                          !< Pointer to memory
    type(MPI_WIN),                          intent(inout) :: win                        !< Corresponding MPI window
    integer,                                intent(in   ) :: i1_new, i2_new, j1_new, j2_new, k1_new, k2_new     !< Bounds of the memory to be allocated

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'reallocate_dist_shared_complex_3D'
    integer                                :: i1_old, i2_old, j1_old, j2_old, k1_old, k2_old, n1_old, n2_old, n3_old, n1_new, n2_new, n3_new
    complex(dp), dimension(:,:,:), pointer :: p_temp => null()
    type(MPI_WIN)                          :: w_temp

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    i1_old = lbound( p,1)
    i2_old = ubound( p,1)
    j1_old = lbound( p,2)
    j2_old = ubound( p,2)
    k1_old = lbound( p,3)
    k2_old = ubound( p,3)
    n1_old = size( p,1)
    n2_old = size( p,2)
    n3_old = size( p,3)
    n1_new = i2_new + 1 - i1_new
    n2_new = j2_new + 1 - j1_new
    n3_new = k2_new + 1 - k1_new

    ! Allocate temporary memory, copy data there
    call allocate_dist_shared( p_temp, w_temp, i1_old, i2_old, j1_old, j2_old, k1_old, k2_old)
    if (par%node_primary) p_temp = p
    call sync_node

    ! Deallocate p
    call deallocate_dist_shared( p, win)

    ! Allocate new memory
    call allocate_dist_shared( p, win, i1_new, i2_new, j1_new, j2_new, k1_new, k2_new)

    ! Copy data there
    if (par%node_primary) p( i1_new: i1_new - 1 + min( n1_old, n1_new), j1_new: j1_new - 1 + min( n2_old, n2_new), k1_new: k1_new - 1 + min( n3_old, n3_new)) = &
      p_temp( i1_old: i1_old - 1 + min( n1_old, n1_new), j1_old: j1_old - 1 + min( n2_old, n2_new), k1_old: k1_old - 1 + min( n3_old, n3_new))
    call sync_node

    ! Deallocate temporary memory
    call deallocate_dist_shared( p_temp, w_temp)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_complex_3D

end module reallocate_dist_shared_mod
