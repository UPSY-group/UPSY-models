module allocate_dist_shared_mod

  use precisions, only: dp
  use mpi_basic, only: par, sync_node
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, n_MPI_windows_used
  use, intrinsic :: iso_c_binding, only: c_ptr, c_f_pointer
  use mpi_f08, only: MPI_WIN, MPI_ADDRESS_KIND, MPI_WIN_ALLOCATE_SHARED, MPI_INFO_NULL, &
    MPI_WIN_SHARED_QUERY

  implicit none

  private

  public :: allocate_dist_shared

  interface allocate_dist_shared
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object
    procedure :: allocate_dist_shared_logical_1D
    procedure :: allocate_dist_shared_logical_2D
    procedure :: allocate_dist_shared_logical_3D
    procedure :: allocate_dist_shared_int_1D
    procedure :: allocate_dist_shared_int_2D
    procedure :: allocate_dist_shared_int_3D
    procedure :: allocate_dist_shared_dp_1D
    procedure :: allocate_dist_shared_dp_2D
    procedure :: allocate_dist_shared_dp_3D
    procedure :: allocate_dist_shared_complex_1D
    procedure :: allocate_dist_shared_complex_2D
    procedure :: allocate_dist_shared_complex_3D
  end interface allocate_dist_shared

contains

  subroutine allocate_dist_shared_logical_1D( p, win, bounds_dim1)
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    logical, dimension(:), pointer, intent(inout) :: p           !< Pointer to memory
    type(MPI_WIN),                  intent(inout) :: win         !< Corresponding MPI window
    integer, dimension(2),          intent(in   ) :: bounds_dim1 !< [lower bound, upper bound] of memory to be allocated

    ! Local variables:
    character(len=*), parameter    :: routine_name = 'allocate_dist_shared_logical_1D'
    integer                        :: ierr
    integer                        :: n1
    integer(kind=MPI_ADDRESS_KIND) :: windowsize
    integer                        :: disp_unit
    type(c_ptr)                    :: baseptr
    logical, dimension(:), pointer :: p_raw

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (associated(p)) call crash('pointer is already/still associated with memory')

    n1 = bounds_dim1(2) + 1 - bounds_dim1(1)

    if (par%node_primary) then
      windowsize = n1*4_MPI_ADDRESS_KIND
      disp_unit  = 4
    else
      windowsize = 0_MPI_ADDRESS_KIND
      disp_unit  = 1
    end if

    call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
      MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

    if (.not. par%primary) then
      ! Get the baseptr, size and disp_unit values of the primary's memory space.
      call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    end if

    ! Associate temporary pointer with this memory space.
    call c_f_pointer( baseptr, p_raw, [n1])

    ! Associate actual pointer (with correct bounds)
    p( bounds_dim1(1):bounds_dim1(2)) => p_raw

    ! Initialise
    if (par%node_primary) p = .false.
    call sync_node

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used + 1

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine allocate_dist_shared_logical_1D

  subroutine allocate_dist_shared_logical_2D( p, win, bounds_dim1, bounds_dim2)
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    logical, dimension(:,:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                    intent(inout) :: win        !< Corresponding MPI window
    integer, dimension(2),            intent(in   ) :: bounds_dim1, bounds_dim2 !< [lower bound, upper bound] of memory to be allocated

    ! Local variables:
    character(len=*), parameter      :: routine_name = 'allocate_dist_shared_logical_2D'
    integer                          :: ierr
    integer                          :: n1, n2
    integer(kind=MPI_ADDRESS_KIND)   :: windowsize
    integer                          :: disp_unit
    type(c_ptr)                      :: baseptr
    logical, dimension(:,:), pointer :: p_raw

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (associated(p)) call crash('pointer is already/still associated with memory')

    n1 = bounds_dim1(2) + 1 - bounds_dim1(1)
    n2 = bounds_dim2(2) + 1 - bounds_dim2(1)

    if (par%node_primary) then
      windowsize = n1*n2*4_MPI_ADDRESS_KIND
      disp_unit  = 4
    else
      windowsize = 0_MPI_ADDRESS_KIND
      disp_unit  = 1
    end if

    call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
      MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

    if (.not. par%primary) then
      ! Get the baseptr, size and disp_unit values of the primary's memory space.
      call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    end if

    ! Associate temporary pointer with this memory space.
    call c_f_pointer( baseptr, p_raw, [n1, n2])

    ! Associate actual pointer (with correct bounds)
    p( bounds_dim1(1):bounds_dim1(2), bounds_dim2(1):bounds_dim2(2)) => p_raw

    ! Initialise
    if (par%node_primary) p = .false.
    call sync_node

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used + 1

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine allocate_dist_shared_logical_2D

  subroutine allocate_dist_shared_logical_3D( p, win, bounds_dim1, bounds_dim2, bounds_dim3)
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    logical, dimension(:,:,:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                      intent(inout) :: win        !< Corresponding MPI window
    integer, dimension(2),              intent(in   ) :: bounds_dim1, bounds_dim2, bounds_dim3 !< [lower bound, upper bound] of memory to be allocated

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'allocate_dist_shared_logical_3D'
    integer                            :: ierr
    integer                            :: n1, n2, n3
    integer(kind=MPI_ADDRESS_KIND)     :: windowsize
    integer                            :: disp_unit
    type(c_ptr)                        :: baseptr
    logical, dimension(:,:,:), pointer :: p_raw

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (associated(p)) call crash('pointer is already/still associated with memory')

    n1 = bounds_dim1(2) + 1 - bounds_dim1(1)
    n2 = bounds_dim2(2) + 1 - bounds_dim2(1)
    n3 = bounds_dim3(2) + 1 - bounds_dim3(1)

    if (par%node_primary) then
      windowsize = n1*n2*n3*4_MPI_ADDRESS_KIND
      disp_unit  = 4
    else
      windowsize = 0_MPI_ADDRESS_KIND
      disp_unit  = 1
    end if

    call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
      MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

    if (.not. par%primary) then
      ! Get the baseptr, size and disp_unit values of the primary's memory space.
      call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    end if

    ! Associate temporary pointer with this memory space.
    call c_f_pointer( baseptr, p_raw, [n1, n2, n3])

    ! Associate actual pointer (with correct bounds)
    p( bounds_dim1(1):bounds_dim1(2), bounds_dim2(1):bounds_dim2(2), bounds_dim3(1):bounds_dim3(2)) => p_raw

    ! Initialise
    if (par%node_primary) p = .false.
    call sync_node

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used + 1

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine allocate_dist_shared_logical_3D

  subroutine allocate_dist_shared_int_1D( p, win, bounds_dim1)
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    integer, dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                  intent(inout) :: win        !< Corresponding MPI window
    integer, dimension(2),          intent(in   ) :: bounds_dim1 !< [lower bound, upper bound] of memory to be allocated

    ! Local variables:
    character(len=*), parameter    :: routine_name = 'allocate_dist_shared_int_1D'
    integer                        :: ierr
    integer                        :: n1
    integer(kind=MPI_ADDRESS_KIND) :: windowsize
    integer                        :: disp_unit
    type(c_ptr)                    :: baseptr
    integer, dimension(:), pointer :: p_raw

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (associated(p)) call crash('pointer is already/still associated with memory')

    n1 = bounds_dim1(2) + 1 - bounds_dim1(1)

    if (par%node_primary) then
      windowsize = n1*4_MPI_ADDRESS_KIND
      disp_unit  = 4
    else
      windowsize = 0_MPI_ADDRESS_KIND
      disp_unit  = 1
    end if

    call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
      MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

    if (.not. par%primary) then
      ! Get the baseptr, size and disp_unit values of the primary's memory space.
      call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    end if

    ! Associate temporary pointer with this memory space.
    call c_f_pointer( baseptr, p_raw, [n1])

    ! Associate actual pointer (with correct bounds)
    p( bounds_dim1(1):bounds_dim1(2)) => p_raw

    ! Initialise
    if (par%node_primary) p = 0
    call sync_node

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used + 1

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine allocate_dist_shared_int_1D

  subroutine allocate_dist_shared_int_2D( p, win, bounds_dim1, bounds_dim2)
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    integer, dimension(:,:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                    intent(inout) :: win        !< Corresponding MPI window
    integer, dimension(2),            intent(in   ) :: bounds_dim1, bounds_dim2 !< [lower bound, upper bound] of memory to be allocated

    ! Local variables:
    character(len=*), parameter      :: routine_name = 'allocate_dist_shared_int_2D'
    integer                          :: ierr
    integer                          :: n1, n2
    integer(kind=MPI_ADDRESS_KIND)   :: windowsize
    integer                          :: disp_unit
    type(c_ptr)                      :: baseptr
    integer, dimension(:,:), pointer :: p_raw

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (associated(p)) call crash('pointer is already/still associated with memory')

    n1 = bounds_dim1(2) + 1 - bounds_dim1(1)
    n2 = bounds_dim2(2) + 1 - bounds_dim2(1)

    if (par%node_primary) then
      windowsize = n1*n2*4_MPI_ADDRESS_KIND
      disp_unit  = 4
    else
      windowsize = 0_MPI_ADDRESS_KIND
      disp_unit  = 1
    end if

    call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
      MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

    if (.not. par%primary) then
      ! Get the baseptr, size and disp_unit values of the primary's memory space.
      call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    end if

    ! Associate temporary pointer with this memory space.
    call c_f_pointer( baseptr, p_raw, [n1, n2])

    ! Associate actual pointer (with correct bounds)
    p( bounds_dim1(1):bounds_dim1(2), bounds_dim2(1):bounds_dim2(2)) => p_raw

    ! Initialise
    if (par%node_primary) p = 0
    call sync_node

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used + 1

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine allocate_dist_shared_int_2D

  subroutine allocate_dist_shared_int_3D( p, win, bounds_dim1, bounds_dim2, bounds_dim3)
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    integer, dimension(:,:,:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                      intent(inout) :: win        !< Corresponding MPI window
    integer, dimension(2),              intent(in   ) :: bounds_dim1, bounds_dim2, bounds_dim3 !< [lower bound, upper bound] of memory to be allocated

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'allocate_dist_shared_int_3D'
    integer                            :: ierr
    integer                            :: n1, n2, n3
    integer(kind=MPI_ADDRESS_KIND)     :: windowsize
    integer                            :: disp_unit
    type(c_ptr)                        :: baseptr
    integer, dimension(:,:,:), pointer :: p_raw

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (associated(p)) call crash('pointer is already/still associated with memory')

    n1 = bounds_dim1(2) + 1 - bounds_dim1(1)
    n2 = bounds_dim2(2) + 1 - bounds_dim2(1)
    n3 = bounds_dim3(2) + 1 - bounds_dim3(1)

    if (par%node_primary) then
      windowsize = n1*n2*n3*4_MPI_ADDRESS_KIND
      disp_unit  = 4
    else
      windowsize = 0_MPI_ADDRESS_KIND
      disp_unit  = 1
    end if

    call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
      MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

    if (.not. par%primary) then
      ! Get the baseptr, size and disp_unit values of the primary's memory space.
      call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    end if

    ! Associate temporary pointer with this memory space.
    call c_f_pointer( baseptr, p_raw, [n1, n2, n3])

    ! Associate actual pointer (with correct bounds)
    p( bounds_dim1(1):bounds_dim1(2), bounds_dim2(1):bounds_dim2(2), bounds_dim3(1):bounds_dim3(2)) => p_raw

    ! Initialise
    if (par%node_primary) p = 0
    call sync_node

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used + 1

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine allocate_dist_shared_int_3D

  subroutine allocate_dist_shared_dp_1D( p, win, bounds_dim1)
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    real(dp), dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                   intent(inout) :: win        !< Corresponding MPI window
    integer, dimension(2),           intent(in   ) :: bounds_dim1 !< [lower bound, upper bound] of memory to be allocated

    ! Local variables:
    character(len=*), parameter     :: routine_name = 'allocate_dist_shared_dp_1D'
    integer                         :: ierr
    integer                         :: n1
    integer(kind=MPI_ADDRESS_KIND)  :: windowsize
    integer                         :: disp_unit
    type(c_ptr)                     :: baseptr
    real(dp), dimension(:), pointer :: p_raw

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (associated(p)) call crash('pointer is already/still associated with memory')

    n1 = bounds_dim1(2) + 1 - bounds_dim1(1)

    if (par%node_primary) then
      windowsize = n1*8_MPI_ADDRESS_KIND
      disp_unit  = 4
    else
      windowsize = 0_MPI_ADDRESS_KIND
      disp_unit  = 1
    end if

    call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
      MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

    if (.not. par%primary) then
      ! Get the baseptr, size and disp_unit values of the primary's memory space.
      call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    end if

    ! Associate temporary pointer with this memory space.
    call c_f_pointer( baseptr, p_raw, [n1])

    ! Associate actual pointer (with correct bounds)
    p( bounds_dim1(1):bounds_dim1(2)) => p_raw

    ! Initialise
    if (par%node_primary) p = 0._dp
    call sync_node

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used + 1

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine allocate_dist_shared_dp_1D

  subroutine allocate_dist_shared_dp_2D( p, win, bounds_dim1, bounds_dim2)
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    real(dp), dimension(:,:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                     intent(inout) :: win        !< Corresponding MPI window
    integer, dimension(2),             intent(in   ) :: bounds_dim1, bounds_dim2 !< [lower bound, upper bound] of memory to be allocated

    ! Local variables:
    character(len=*), parameter       :: routine_name = 'allocate_dist_shared_dp_2D'
    integer                           :: ierr
    integer                           :: n1, n2
    integer(kind=MPI_ADDRESS_KIND)    :: windowsize
    integer                           :: disp_unit
    type(c_ptr)                       :: baseptr
    real(dp), dimension(:,:), pointer :: p_raw

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (associated(p)) call crash('pointer is already/still associated with memory')

    n1 = bounds_dim1(2) + 1 - bounds_dim1(1)
    n2 = bounds_dim2(2) + 1 - bounds_dim2(1)

    if (par%node_primary) then
      windowsize = n1*n2*8_MPI_ADDRESS_KIND
      disp_unit  = 4
    else
      windowsize = 0_MPI_ADDRESS_KIND
      disp_unit  = 1
    end if

    call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
      MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

    if (.not. par%primary) then
      ! Get the baseptr, size and disp_unit values of the primary's memory space.
      call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    end if

    ! Associate temporary pointer with this memory space.
    call c_f_pointer( baseptr, p_raw, [n1, n2])

    ! Associate actual pointer (with correct bounds)
    p( bounds_dim1(1):bounds_dim1(2), bounds_dim2(1):bounds_dim2(2)) => p_raw

    ! Initialise
    if (par%node_primary) p = 0._dp
    call sync_node

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used + 1

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine allocate_dist_shared_dp_2D

  subroutine allocate_dist_shared_dp_3D( p, win, bounds_dim1, bounds_dim2, bounds_dim3)
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    real(dp), dimension(:,:,:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                       intent(inout) :: win        !< Corresponding MPI window
    integer, dimension(2),               intent(in   ) :: bounds_dim1, bounds_dim2, bounds_dim3 !< [lower bound, upper bound] of memory to be allocated

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'allocate_dist_shared_dp_3D'
    integer                             :: ierr
    integer                             :: n1, n2, n3
    integer(kind=MPI_ADDRESS_KIND)      :: windowsize
    integer                             :: disp_unit
    type(c_ptr)                         :: baseptr
    real(dp), dimension(:,:,:), pointer :: p_raw

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (associated(p)) call crash('pointer is already/still associated with memory')

    n1 = bounds_dim1(2) + 1 - bounds_dim1(1)
    n2 = bounds_dim2(2) + 1 - bounds_dim2(1)
    n3 = bounds_dim3(2) + 1 - bounds_dim3(1)

    if (par%node_primary) then
      windowsize = n1*n2*n3*8_MPI_ADDRESS_KIND
      disp_unit  = 4
    else
      windowsize = 0_MPI_ADDRESS_KIND
      disp_unit  = 1
    end if

    call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
      MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

    if (.not. par%primary) then
      ! Get the baseptr, size and disp_unit values of the primary's memory space.
      call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    end if

    ! Associate temporary pointer with this memory space.
    call c_f_pointer( baseptr, p_raw, [n1, n2, n3])

    ! Associate actual pointer (with correct bounds)
    p( bounds_dim1(1):bounds_dim1(2), bounds_dim2(1):bounds_dim2(2), bounds_dim3(1):bounds_dim3(2)) => p_raw

    ! Initialise
    if (par%node_primary) p = 0._dp
    call sync_node

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used + 1

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine allocate_dist_shared_dp_3D

  subroutine allocate_dist_shared_complex_1D( p, win, bounds_dim1)
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    complex(dp), dimension(:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                      intent(inout) :: win        !< Corresponding MPI window
    integer, dimension(2),              intent(in   ) :: bounds_dim1 !< [lower bound, upper bound] of memory to be allocated

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'allocate_dist_shared_complex_1D'
    integer                            :: ierr
    integer                            :: n1
    integer(kind=MPI_ADDRESS_KIND)     :: windowsize
    integer                            :: disp_unit
    type(c_ptr)                        :: baseptr
    complex(dp), dimension(:), pointer :: p_raw

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (associated(p)) call crash('pointer is already/still associated with memory')

    n1 = bounds_dim1(2) + 1 - bounds_dim1(1)

    if (par%node_primary) then
      windowsize = n1*16_MPI_ADDRESS_KIND
      disp_unit  = 4
    else
      windowsize = 0_MPI_ADDRESS_KIND
      disp_unit  = 1
    end if

    call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
      MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

    if (.not. par%primary) then
      ! Get the baseptr, size and disp_unit values of the primary's memory space.
      call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    end if

    ! Associate temporary pointer with this memory space.
    call c_f_pointer( baseptr, p_raw, [n1])

    ! Associate actual pointer (with correct bounds)
    p( bounds_dim1(1):bounds_dim1(2)) => p_raw

    ! Initialise
    if (par%node_primary) p = 0._dp
    call sync_node

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used + 1

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine allocate_dist_shared_complex_1D

  subroutine allocate_dist_shared_complex_2D( p, win, bounds_dim1, bounds_dim2)
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    complex(dp), dimension(:,:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                        intent(inout) :: win        !< Corresponding MPI window
    integer, dimension(2),                intent(in   ) :: bounds_dim1, bounds_dim2 !< [lower bound, upper bound] of memory to be allocated

    ! Local variables:
    character(len=*), parameter          :: routine_name = 'allocate_dist_shared_complex_2D'
    integer                              :: ierr
    integer                              :: n1, n2
    integer(kind=MPI_ADDRESS_KIND)       :: windowsize
    integer                              :: disp_unit
    type(c_ptr)                          :: baseptr
    complex(dp), dimension(:,:), pointer :: p_raw

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (associated(p)) call crash('pointer is already/still associated with memory')

    n1 = bounds_dim1(2) + 1 - bounds_dim1(1)
    n2 = bounds_dim2(2) + 1 - bounds_dim2(1)

    if (par%node_primary) then
      windowsize = n1*n2*16_MPI_ADDRESS_KIND
      disp_unit  = 4
    else
      windowsize = 0_MPI_ADDRESS_KIND
      disp_unit  = 1
    end if

    call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
      MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

    if (.not. par%primary) then
      ! Get the baseptr, size and disp_unit values of the primary's memory space.
      call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    end if

    ! Associate temporary pointer with this memory space.
    call c_f_pointer( baseptr, p_raw, [n1, n2])

    ! Associate actual pointer (with correct bounds)
    p( bounds_dim1(1):bounds_dim1(2), bounds_dim2(1):bounds_dim2(2)) => p_raw

    ! Initialise
    if (par%node_primary) p = 0._dp
    call sync_node

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used + 1

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine allocate_dist_shared_complex_2D

  subroutine allocate_dist_shared_complex_3D( p, win, bounds_dim1, bounds_dim2, bounds_dim3)
    !< Allocate hybrid distributed/shared memory, with an associated MPI window object

    ! In/output variables:
    complex(dp), dimension(:,:,:), pointer, intent(inout) :: p          !< Pointer to memory
    type(MPI_WIN),                          intent(inout) :: win        !< Corresponding MPI window
    integer, dimension(2),                  intent(in   ) :: bounds_dim1, bounds_dim2, bounds_dim3 !< [lower bound, upper bound] of memory to be allocated

    ! Local variables:
    character(len=*), parameter            :: routine_name = 'allocate_dist_shared_complex_3D'
    integer                                :: ierr
    integer                                :: n1, n2, n3
    integer(kind=MPI_ADDRESS_KIND)         :: windowsize
    integer                                :: disp_unit
    type(c_ptr)                            :: baseptr
    complex(dp), dimension(:,:,:), pointer :: p_raw

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (associated(p)) call crash('pointer is already/still associated with memory')

    n1 = bounds_dim1(2) + 1 - bounds_dim1(1)
    n2 = bounds_dim2(2) + 1 - bounds_dim2(1)
    n3 = bounds_dim3(2) + 1 - bounds_dim3(1)

    if (par%node_primary) then
      windowsize = n1*n2*n3*16_MPI_ADDRESS_KIND
      disp_unit  = 4
    else
      windowsize = 0_MPI_ADDRESS_KIND
      disp_unit  = 1
    end if

    call MPI_WIN_ALLOCATE_SHARED( windowsize, disp_unit, &
      MPI_INFO_NULL, par%mpi_comm_node, baseptr, win, ierr)

    if (.not. par%primary) then
      ! Get the baseptr, size and disp_unit values of the primary's memory space.
      call MPI_WIN_SHARED_QUERY( win, 0, windowsize, disp_unit, baseptr, ierr)
    end if

    ! Associate temporary pointer with this memory space.
    call c_f_pointer( baseptr, p_raw, [n1, n2, n3])

    ! Associate actual pointer (with correct bounds)
    p( bounds_dim1(1):bounds_dim1(2), bounds_dim2(1):bounds_dim2(2), bounds_dim3(1):bounds_dim3(2)) => p_raw

    ! Initialise
    if (par%node_primary) p = 0._dp
    call sync_node

    ! Update the shared memory leak tracker
    n_MPI_windows_used = n_MPI_windows_used + 1

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine allocate_dist_shared_complex_3D

end module allocate_dist_shared_mod
