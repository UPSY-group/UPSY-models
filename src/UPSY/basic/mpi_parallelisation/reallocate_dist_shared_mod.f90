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

  subroutine reallocate_dist_shared_logical_1D( p, win, n1_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    logical, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                  intent(inout) :: win        !< Corresponding MPI window
    integer,                        intent(in   ) :: n1_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=*), parameter :: routine_name = 'reallocate_dist_shared_logical_1D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    call deallocate_dist_shared( p, win)
    call allocate_dist_shared( p, win, n1_new)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_logical_1D

  subroutine reallocate_dist_shared_logical_2D( p, win, n1_new, n2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    logical, dimension(:,:), pointer, intent(inout) :: p                  !< Pointer to memory
    type(MPI_WIN),                    intent(inout) :: win                !< Corresponding MPI window
    integer,                          intent(in   ) :: n1_new, n2_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'reallocate_dist_shared_logical_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    call deallocate_dist_shared( p, win)
    call allocate_dist_shared( p, win, n1_new, n2_new)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_logical_2D

  subroutine reallocate_dist_shared_logical_3D( p, win, n1_new, n2_new, n3_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    logical, dimension(:,:,:), pointer, intent(inout) :: p                          !< Pointer to memory
    type(MPI_WIN),                      intent(inout) :: win                        !< Corresponding MPI window
    integer,                            intent(in   ) :: n1_new, n2_new, n3_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=*), parameter     :: routine_name = 'reallocate_dist_shared_logical_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    call deallocate_dist_shared( p, win)
    call allocate_dist_shared( p, win, n1_new, n2_new, n3_new)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_logical_3D

  subroutine reallocate_dist_shared_int_1D( p, win, n1_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    integer, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                  intent(inout) :: win        !< Corresponding MPI window
    integer,                        intent(in   ) :: n1_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=*), parameter :: routine_name = 'reallocate_dist_shared_int_1D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    call deallocate_dist_shared( p, win)
    call allocate_dist_shared( p, win, n1_new)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_int_1D

  subroutine reallocate_dist_shared_int_2D( p, win, n1_new, n2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    integer, dimension(:,:), pointer, intent(inout) :: p                  !< Pointer to memory
    type(MPI_WIN),                    intent(inout) :: win                !< Corresponding MPI window
    integer,                          intent(in   ) :: n1_new, n2_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'reallocate_dist_shared_int_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    call deallocate_dist_shared( p, win)
    call allocate_dist_shared( p, win, n1_new, n2_new)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_int_2D

  subroutine reallocate_dist_shared_int_3D( p, win, n1_new, n2_new, n3_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    integer, dimension(:,:,:), pointer, intent(inout) :: p                          !< Pointer to memory
    type(MPI_WIN),                      intent(inout) :: win                        !< Corresponding MPI window
    integer,                            intent(in   ) :: n1_new, n2_new, n3_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=*), parameter     :: routine_name = 'reallocate_dist_shared_int_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    call deallocate_dist_shared( p, win)
    call allocate_dist_shared( p, win, n1_new, n2_new, n3_new)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_int_3D

  subroutine reallocate_dist_shared_dp_1D( p, win, n1_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    real(dp), dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                   intent(inout) :: win        !< Corresponding MPI window
    integer,                         intent(in   ) :: n1_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=*), parameter  :: routine_name = 'reallocate_dist_shared_dp_1D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    call deallocate_dist_shared( p, win)
    call allocate_dist_shared( p, win, n1_new)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_dp_1D

  subroutine reallocate_dist_shared_dp_2D( p, win, n1_new, n2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    real(dp), dimension(:,:), pointer, intent(inout) :: p                  !< Pointer to memory
    type(MPI_WIN),                     intent(inout) :: win                !< Corresponding MPI window
    integer,                           intent(in   ) :: n1_new, n2_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=*), parameter    :: routine_name = 'reallocate_dist_shared_dp_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    call deallocate_dist_shared( p, win)
    call allocate_dist_shared( p, win, n1_new, n2_new)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_dp_2D

  subroutine reallocate_dist_shared_dp_3D( p, win, n1_new, n2_new, n3_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    real(dp), dimension(:,:,:), pointer, intent(inout) :: p                          !< Pointer to memory
    type(MPI_WIN),                       intent(inout) :: win                        !< Corresponding MPI window
    integer,                             intent(in   ) :: n1_new, n2_new, n3_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=*), parameter      :: routine_name = 'reallocate_dist_shared_dp_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    call deallocate_dist_shared( p, win)
    call allocate_dist_shared( p, win, n1_new, n2_new, n3_new)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_dp_3D

  subroutine reallocate_dist_shared_complex_1D( p, win, n1_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    complex(dp), dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                      intent(inout) :: win        !< Corresponding MPI window
    integer,                            intent(in   ) :: n1_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=*), parameter     :: routine_name = 'reallocate_dist_shared_complex_1D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    call deallocate_dist_shared( p, win)
    call allocate_dist_shared( p, win, n1_new)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_complex_1D

  subroutine reallocate_dist_shared_complex_2D( p, win, n1_new, n2_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    complex(dp), dimension(:,:), pointer, intent(inout) :: p                  !< Pointer to memory
    type(MPI_WIN),                        intent(inout) :: win                !< Corresponding MPI window
    integer,                              intent(in   ) :: n1_new, n2_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=*), parameter       :: routine_name = 'reallocate_dist_shared_complex_2D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    call deallocate_dist_shared( p, win)
    call allocate_dist_shared( p, win, n1_new, n2_new)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_complex_2D

  subroutine reallocate_dist_shared_complex_3D( p, win, n1_new, n2_new, n3_new)
    !< Reallocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    complex(dp), dimension(:,:,:), pointer, intent(inout) :: p                          !< Pointer to memory
    type(MPI_WIN),                          intent(inout) :: win                        !< Corresponding MPI window
    integer,                                intent(in   ) :: n1_new, n2_new, n3_new     !< Dimension(s) of memory to be allocated

    ! Local variables:
    character(len=*), parameter :: routine_name = 'reallocate_dist_shared_complex_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (.not. associated(p)) call crash('pointer not associated with any memory')

    call deallocate_dist_shared( p, win)
    call allocate_dist_shared( p, win, n1_new, n2_new, n3_new)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine reallocate_dist_shared_complex_3D

end module reallocate_dist_shared_mod
