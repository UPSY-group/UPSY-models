module CSR_matrix_mod

  ! The Compressed Sparse Row matrix type

  use precisions, only: dp
  use parallel_array_info_type, only: type_par_arr_info
  use mpi_f08, only: MPI_ALLGATHER, MPI_INTEGER, MPI_ANY_TAG, MPI_COMM_WORLD, MPI_SEND, MPI_RECV, &
    MPI_STATUS, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_ALLREDUCE, MPI_MIN, MPI_MAX, MPI_IN_PLACE, &
    MPI_LOGICAL, MPI_LOR
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use reallocate_mod, only: reallocate
  use netcdf, only: NF90_INT, NF90_DOUBLE, NF90_CREATE, NF90_NOCLOBBER, NF90_NETCDF4, NF90_DEF_DIM, &
    NF90_DEF_VAR, NF90_INT, NF90_DOUBLE, NF90_PUT_VAR, NF90_OPEN, NF90_CLOSE, NF90_NOWRITE, &
    NF90_INQ_DIMID, NF90_INQUIRE_DIMENSION, NF90_INQ_VARID, NF90_GET_VAR
  use netcdf_basic_wrappers, only: create_dimension, create_variable, handle_netcdf_error, &
    create_and_write_to_scalar_variable_dist_int, read_scalar_variable_dist_int, &
    create_new_netcdf_file_for_writing, close_netcdf_file
  use netcdf_write_var_primary, only: write_var_primary
  use string_module, only: endswith, insert_val_into_string_int
  use crash_mod, only: warning, crash

  implicit none

  private

  public :: type_CSR_matrix_dp

  ! The basic CSR matrix type
  type type_CSR_matrix_dp
    ! Compressed Sparse Row (CSR) format matrix

      integer                             :: m,n         ! A = [m-by-n]
      integer                             :: nnz_max     ! Maximum number of non-zero entries in A
      integer                             :: nnz         ! Actual  number of non-zero entries in A
      integer,  dimension(:), allocatable :: ptr         ! Row start indices
      integer,  dimension(:), allocatable :: ind         ! Column indices
      real(dp), dimension(:), allocatable :: val         ! Values

      ! Parallelisation
      integer :: m_loc,  i1,      i2      ! Rows    owned by each process
      integer :: m_node, i1_node, i2_node ! Rows    owned by each node
      integer :: n_loc,  j1,      j2      ! Columns owned by each process
      integer :: n_node, j1_node, j2_node ! Columns owned by each node

      type(type_par_arr_info) :: pai_x
      type(type_par_arr_info) :: pai_y

      logical :: is_finalised = .false.
      integer :: j_min_node, j_max_node   ! Range of rows of x needed by this node to compute y = A*x
      integer :: needs_x_tot = -1         ! Whether y = A*x can be evaluated with just x_nih, or if it needs x_tot (0: can use x_nih; 1: needs x_tot; -1: not checked yet)

    contains

      private

      procedure, public  :: allocate               => allocate_matrix_CSR_dist
      procedure, public  :: allocate_loc           => allocate_matrix_CSR_loc
      procedure, public  :: deallocate             => deallocate_matrix_CSR_dist
      procedure, public  :: duplicate              => duplicate_matrix_CSR_dist
      procedure, public  :: add_entry              => add_entry_CSR_dist
      procedure, public  :: add_empty_row          => add_empty_row_CSR_dist
      procedure, public  :: gather_to_primary      => gather_CSR_dist_to_primary
      procedure, public  :: read_single_row        => read_single_row_CSR_dist
      procedure, public  :: finalise               => finalise_matrix_CSR_dist
      procedure, public  :: set_diagonal_to_one_and_rest_of_row_to_zero
      procedure, public  :: save_as_NetCDF         => save_CSR_matrix_as_NetCDF
      procedure, public  :: write_to_NetCDF        => write_CSR_matrix_to_NetCDF
      procedure, public  :: write_to_dist_NetCDFs  => write_CSR_matrix_to_dist_NetCDFs
      procedure, public  :: read_from_dist_NetCDFs => read_CSR_matrix_from_dist_NetCDFs

      procedure, private :: extend                 => extend_matrix_CSR_dist
      procedure, private :: crop                   => crop_matrix_CSR_dist
      procedure, private :: calc_j_node_range
      procedure, private :: set_row_to_value
      procedure, private :: set_row_diag_to_val

      generic,   public  :: operator(==) => eq
      procedure, private :: eq => test_CSR_matrix_equality

      final              :: deallocate_matrix_CSR_dist_final

  end type type_CSR_matrix_dp

contains

  ! ===== CSR matrices in distributed memory =====
  ! ==============================================

  subroutine allocate_matrix_CSR_dist( A, m_glob, n_glob, m_loc, n_loc, nnz_max_proc, pai_x, pai_y)
    ! Allocate memory for a CSR-format sparse m-by-n matrix A

    ! In- and output variables:
    class(type_CSR_matrix_dp),         intent(inout) :: A
    integer,                           intent(in   ) :: m_glob, n_glob, m_loc, n_loc, nnz_max_proc
    type(type_par_arr_info), optional, intent(in   ) :: pai_x, pai_y

    ! Local variables:
    character(len=*), parameter :: routine_name = 'allocate_matrix_CSR_dist'
    integer                     :: ierr
    integer, dimension(par%n)   :: m_loc_all, n_loc_all

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Matrix dimensions
    A%m       = m_glob
    A%n       = n_glob
    A%m_loc   = m_loc
    A%n_loc   = n_loc
    A%nnz_max = nnz_max_proc
    A%nnz     = 0

    ! Partition rows and columns over the processes
    call MPI_ALLGATHER( m_loc, 1, MPI_integer, m_loc_all, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( n_loc, 1, MPI_integer, n_loc_all, 1, MPI_integer, MPI_COMM_WORLD, ierr)

    ! Safety
    if (sum( m_loc_all) /= m_glob) call crash('sum of numbers of local rows doesnt match number of global rows!')
    if (sum( n_loc_all) /= n_glob) call crash('sum of numbers of local columns doesnt match number of global columns!')

    A%i1 = 1 + sum( m_loc_all( 1:par%i  ))
    A%i2 = 1 + sum( m_loc_all( 1:par%i+1))-1
    A%j1 = 1 + sum( n_loc_all( 1:par%i  ))
    A%j2 = 1 + sum( n_loc_all( 1:par%i+1))-1

    ! Range owned by this node
    call MPI_ALLREDUCE( A%i1, A%i1_node, 1, MPI_INTEGER, MPI_MIN, par%mpi_comm_node, ierr)
    call MPI_ALLREDUCE( A%i2, A%i2_node, 1, MPI_INTEGER, MPI_MAX, par%mpi_comm_node, ierr)
    A%m_node = A%i2_node + 1 - A%i1_node

    call MPI_ALLREDUCE( A%j1, A%j1_node, 1, MPI_INTEGER, MPI_MIN, par%mpi_comm_node, ierr)
    call MPI_ALLREDUCE( A%j2, A%j2_node, 1, MPI_INTEGER, MPI_MAX, par%mpi_comm_node, ierr)
    A%n_node = A%j2_node + 1 - A%j1_node

    ! Allocate memory
    allocate( A%ptr( A%i1: A%i2+1), source = 1)
    allocate( A%ind( A%nnz_max), source = 0    )
    allocate( A%val( A%nnz_max), source = 0._dp)

    ! Parallel array info
    if (present( pai_x)) then
      A%pai_x = pai_x
    else
      A%pai_x%n       = A%n

      A%pai_x%n_loc   = A%n_loc
      A%pai_x%i1      = A%j1
      A%pai_x%i2      = A%j2

      A%pai_x%n_node  = A%n_node
      A%pai_x%i1_node = A%j1_node
      A%pai_x%i2_node = A%j2_node

      A%pai_x%n_nih   = A%pai_x%n_node
      A%pai_x%i1_nih  = A%pai_x%i1_node
      A%pai_x%i2_nih  = A%pai_x%i2_node

      A%pai_x%n_hle   = 0
      A%pai_x%i1_hle  = 0
      A%pai_x%i2_hle  = -1

      A%pai_x%n_hli   = 0
      A%pai_x%i1_hli  = 0
      A%pai_x%i2_hli  = -1

      A%pai_x%n_hre   = 0
      A%pai_x%i1_hre  = 0
      A%pai_x%i2_hre  = -1

      A%pai_x%n_hri   = 0
      A%pai_x%i1_hri  = 0
      A%pai_x%i2_hri  = -1
    end if

    if (present( pai_y)) then
      A%pai_y = pai_y
    else
      A%pai_y%n       = A%m

      A%pai_y%n_loc   = A%m_loc
      A%pai_y%i1      = A%i1
      A%pai_y%i2      = A%i2

      A%pai_y%n_node  = A%m_node
      A%pai_y%i1_node = A%i1_node
      A%pai_y%i2_node = A%i2_node

      A%pai_y%n_nih   = A%pai_y%n_node
      A%pai_y%i1_nih  = A%pai_y%i1_node
      A%pai_y%i2_nih  = A%pai_y%i2_node

      A%pai_y%n_hle   = 0
      A%pai_y%i1_hle  = 0
      A%pai_y%i2_hle  = -1

      A%pai_y%n_hli   = 0
      A%pai_y%i1_hli  = 0
      A%pai_y%i2_hli  = -1

      A%pai_y%n_hre   = 0
      A%pai_y%i1_hre  = 0
      A%pai_y%i2_hre  = -1

      A%pai_y%n_hri   = 0
      A%pai_y%i1_hri  = 0
      A%pai_y%i2_hri  = -1
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_matrix_CSR_dist

  subroutine deallocate_matrix_CSR_dist( A)
    !< Deallocate memory for a CSR-format sparse matrix A

    ! In- and output variables:
    class(type_CSR_matrix_dp), intent(inout) :: A

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_matrix_CSR_dist'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Matrix dimensions
    A%m       = 0
    A%n       = 0
    A%nnz_max = 0
    A%nnz     = 0

    ! Parallelisation ranges
    A%m_loc   = 0
    A%i1      = 0
    A%i2      = 0

    A%n_loc   = 0
    A%j1      = 0
    A%j2      = 0

    A%m_node  = 0
    A%i1_node = 0
    A%i2_node = 0

    A%n_node  = 0
    A%j1_node = 0
    A%j2_node = 0

    if (allocated( A%ptr)) deallocate( A%ptr)
    if (allocated( A%ind)) deallocate( A%ind)
    if (allocated( A%val)) deallocate( A%val)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_matrix_CSR_dist

  subroutine deallocate_matrix_CSR_dist_final( A)
    type(type_CSR_matrix_dp), intent(inout) :: A
    call A%deallocate
  end subroutine deallocate_matrix_CSR_dist_final

  subroutine duplicate_matrix_CSR_dist( A, B)
    ! Duplicate the CSR-format sparse matrix A

    ! In- and output variables:
    class(type_CSR_matrix_dp), intent(in   ) :: A
    type( type_CSR_matrix_dp), intent(inout) :: B

    ! Local variables:
    character(len=*), parameter :: routine_name = 'duplicate_matrix_CSR_dist'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate any memory that was allcoated for B
    call B%deallocate

    ! Matrix dimensions
    B%m       = A%m
    B%n       = A%n
    B%nnz_max = A%nnz_max
    B%nnz     = A%nnz

    ! Parallelisation ranges
    B%m_loc   = A%m_loc
    B%i1      = A%i1
    B%i2      = A%i2

    B%n_loc   = A%n_loc
    B%j1      = A%j1
    B%j2      = A%j2

    B%m_node  = A%m_node
    B%i1_node = A%i1_node
    B%i2_node = A%i2_node

    B%n_node  = A%n_node
    B%j1_node = A%j1_node
    B%j2_node = A%j2_node

    B%pai_x = A%pai_x
    B%pai_y = A%pai_y

    ! Allocate memory
    allocate( B%ptr( B%i1: B%i2+1    ), source = 1    )
    allocate( B%ind( B%nnz_max), source = 0    )
    allocate( B%val( B%nnz_max), source = 0._dp)

    ! Copy data
    B%ptr = A%ptr
    B%ind = A%ind
    B%val = A%val

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine duplicate_matrix_CSR_dist

  subroutine add_entry_CSR_dist( A, i, j, v)
    ! Add value v to row i, column j of CSR-formatted matrix A
    !
    ! NOTE: assumes all rows before i are finished and nothing exists yet for rows after i!

    ! In- and output variables:
    class(type_CSR_matrix_dp), intent(inout) :: A
    integer,                   intent(in   ) :: i,j
    real(dp),                  intent(in   ) :: v

    ! Safety
    if (i < A%i1 .or. i > A%i2) call crash('out of ownership range!')

    ! Increase number of non-zeros
    A%nnz = A%nnz + 1

    ! List entry
    A%ind( A%nnz) = j
    A%val( A%nnz) = v

    ! Update pointer list
    A%ptr( i+1) = A%nnz+1

    ! Extend memory if necessary
    if (A%nnz > A%nnz_max - 10) call A%extend( 1000)

  end subroutine add_entry_CSR_dist

  subroutine add_empty_row_CSR_dist( A, i)
    ! Add an empty row i to CSR-formatted matrix A
    !
    ! NOTE: assumes all rows before i are finished and nothing exists yet for rows after i!

    ! In- and output variables:
    class(type_CSR_matrix_dp), intent(inout) :: A
    integer,                   intent(in)    :: i

    ! Safety
    if (i < A%i1 .or. i > A%i2) call crash('out of ownership range!')

    ! Update pointer list
    A%ptr( i+1) = A%nnz+1

  end subroutine add_empty_row_CSR_dist

  subroutine extend_matrix_CSR_dist( A, nnz_extra)
    ! Extend memory for a CSR-format sparse m-by-n matrix A
    class(type_CSR_matrix_dp), intent(inout) :: A
    integer,                   intent(in   ) :: nnz_extra
    A%nnz_max = A%nnz + nnz_extra
    call reallocate( A%ind, A%nnz_max)
    call reallocate( A%val, A%nnz_max)
  end subroutine extend_matrix_CSR_dist

  subroutine crop_matrix_CSR_dist( A)
    ! Crop memory for a CSR-format sparse m-by-n matrix A
    class(type_CSR_matrix_dp), intent(inout) :: A
    A%nnz_max = A%nnz
    call reallocate( A%ind, A%nnz_max)
    call reallocate( A%val, A%nnz_max)
  end subroutine crop_matrix_CSR_dist

  subroutine gather_CSR_dist_to_primary( A, A_tot)
    ! Gather a CSR-format sparse m-by-n matrix A that is distributed over the processes, to the primary

    ! In- and output variables:
    class(type_CSR_matrix_dp), intent(in   ) :: A
    type(type_CSR_matrix_dp),  intent(inout) :: A_tot

    ! Local variables:
    character(len=*), parameter :: routine_name = 'gather_CSR_dist_to_primary'
    integer                     :: ierr
    type(MPI_STATUS)            :: recv_status
    integer,  dimension(par%n)  :: m_glob_all, n_glob_all, m_loc_all, n_loc_all
    integer                     :: nnz_tot
    integer                     :: p
    integer                     :: row, k1, k2, k, col
    real(dp)                    :: val
    type(type_CSR_matrix_dp)    :: A_proc

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Gather dimensions
    call MPI_ALLGATHER( A%m    , 1, MPI_integer, m_glob_all, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( A%n    , 1, MPI_integer, n_glob_all, 1, MPI_integer, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( A%m_loc, 1, MPI_integer, m_loc_all , 1, MPI_integer, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( A%n_loc, 1, MPI_integer, n_loc_all , 1, MPI_integer, MPI_COMM_WORLD, ierr)

    call MPI_ALLREDUCE( A%nnz, nnz_tot, 1, MPI_integer, MPI_sum, MPI_COMM_WORLD, ierr)

    ! Safety - check if dimensions match
    if (any( m_glob_all /= A%m)) call crash('global numbers of rows do not match across the processes!')
    if (any( n_glob_all /= A%n)) call crash('global numbers of columns do not match across the processes!')
    if (sum( m_loc_all) /= A%m ) call crash('local numbers of rows do not add up across the processes!')
    if (sum( n_loc_all) /= A%n ) call crash('local numbers of columns do not add up across the processes!')

    ! Allocate memory
    if (par%primary) then
      A_tot%m       = A%m
      A_tot%m_loc   = A%m
      A_tot%i1      = 1
      A_tot%i2      = A%m
      A_tot%n       = A%n
      A_tot%n_loc   = A%n
      A_tot%j1      = 1
      A_tot%j2      = A%n
      A_tot%nnz     = 0
      A_tot%nnz_max = nnz_tot
      allocate( A_tot%ptr( A%m+1)  , source = 1    )
      allocate( A_tot%ind( nnz_tot), source = 0    )
      allocate( A_tot%val( nnz_tot), source = 0._dp)
    else
      A_tot%m       = A%m
      A_tot%m_loc   = 0
      A_tot%i1      = 1
      A_tot%i2      = 0
      A_tot%n       = A%n
      A_tot%n_loc   = 0
      A_tot%j1      = 1
      A_tot%j2      = 0
      A_tot%nnz     = 0
      A_tot%nnz_max = 0
    end if

    ! Start with the primary's own data
    if (par%primary) then
      do row = A%i1, A%i2
        k1 = A%ptr( row)
        k2 = A%ptr( row+1) - 1
        do k = k1, k2
          col = A%ind( k)
          val = A%val( k)
          call add_entry_CSR_dist( A_tot, row, col, val)
        end do
      end do
    end if

    ! Collect data from the other processes
    do p = 1, par%n-1

      if     (par%i == p) then

        ! Send matrix metadata to primary
        call MPI_Send( A%m      , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%m_loc  , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%i1     , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%i2     , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%n      , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%n_loc  , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%j1     , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%j2     , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%nnz    , 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%nnz_max, 1, MPI_integer, 0, 0, MPI_COMM_WORLD, ierr)

        ! Send matrix data to primary
        call MPI_Send( A%ptr, A%m_loc+1, MPI_integer         , 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%ind, A%nnz_max, MPI_integer         , 0, 0, MPI_COMM_WORLD, ierr)
        call MPI_Send( A%val, A%nnz_max, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierr)

      elseif (par%primary) then

        ! Receive matrix metadata from process
        call MPI_RECV( A_proc%m      , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%m_loc  , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%i1     , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%i2     , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%n      , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%n_loc  , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%j1     , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%j2     , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%nnz    , 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%nnz_max, 1, MPI_integer, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)

        ! Allocate memory
        allocate( A_proc%ptr( A_proc%i1: A_proc%i2+1), source = 0    )
        allocate( A_proc%ind( A_proc%nnz_max), source = 0    )
        allocate( A_proc%val( A_proc%nnz_max), source = 0._dp)

        ! Receive matrix data from process
        call MPI_RECV( A_proc%ptr, A_proc%m_loc+1, MPI_integer         , p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%ind, A_proc%nnz_max, MPI_integer         , p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)
        call MPI_RECV( A_proc%val, A_proc%nnz_max, MPI_DOUBLE_PRECISION, p, MPI_any_TAG, MPI_COMM_WORLD, recv_status, ierr)

        ! Write to total matrix
        do row = A_proc%i1, A_proc%i2
          k1 = A_proc%ptr( row)
          k2 = A_proc%ptr( row+1) - 1
          do k = k1, k2
            col = A_proc%ind( k)
            val = A_proc%val( k)
            call add_entry_CSR_dist( A_tot, row, col, val)
          end do
        end do

        ! Clean up after yourself
        deallocate( A_proc%ptr)
        deallocate( A_proc%ind)
        deallocate( A_proc%val)

      end if ! if     (par%i == p) then
      call sync

    end do ! do p = 1, par%n-1

    if (par%primary) call A_tot%crop

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_CSR_dist_to_primary

  subroutine read_single_row_CSR_dist( A, i, ind, val, nnz)
    ! Read the coefficients of a single row of A

    ! In- and output variables:
    class(type_CSR_matrix_dp),  intent(in   ) :: A
    integer,                    intent(in   ) :: i
    integer,  dimension(:    ), intent(  out) :: ind
    real(dp), dimension(:    ), intent(  out) :: val
    integer,                    intent(  out) :: nnz

    ! Local variables:
    integer :: k1,k2

    ! Safety
    if (i < A%i1 .or. i > A%i2) call crash('row {int_01} is not owned by process {int_02}!', int_01 = i, int_02 = par%i)

    k1 = A%ptr( i)
    k2 = A%ptr( i+1) - 1

    nnz = k2 + 1 - k1

    ind( 1:nnz) = A%ind( k1:k2)
    val( 1:nnz) = A%val( k1:k2)

  end subroutine read_single_row_CSR_dist

  subroutine finalise_matrix_CSR_dist( A)

    ! In- and output variables:
    class(type_CSR_matrix_dp), intent(inout) :: A

    ! Local variables:
    character(len=*), parameter :: routine_name = 'finalise_matrix_CSR_dist'

    ! Add routine to call stack
    call init_routine( routine_name)

    call A%crop
    call A%calc_j_node_range

    A%is_finalised = .true.

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine finalise_matrix_CSR_dist

  subroutine calc_j_node_range( A)

    ! In- and output variables:
    class(type_CSR_matrix_dp), intent(inout) :: A

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_j_node_range'
    integer                        :: i, k1, k2, k, j, ierr
    logical                        :: needs_x_tot

    ! Add routine to call stack
    call init_routine( routine_name)

    A%j_min_node =  huge( A%j_min_node)
    A%j_max_node = -huge( A%j_max_node)

    do i = A%i1, A%i2

      k1 = A%ptr( i)
      k2 = A%ptr( i+1)-1

      do k = k1, k2
        j = A%ind( k)
        A%j_min_node = min( A%j_min_node, j)
        A%j_max_node = max( A%j_max_node, j)
      end do

    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, A%j_min_node, 1, MPI_INTEGER, MPI_MIN, par%mpi_comm_node, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, A%j_max_node, 1, MPI_INTEGER, MPI_MAX, par%mpi_comm_node, ierr)

    needs_x_tot = A%j_min_node < A%pai_x%i1_nih .or. A%j_max_node > A%pai_x%i2_nih
    call MPI_ALLREDUCE( MPI_IN_PLACE, needs_x_tot, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ierr)
    if (needs_x_tot) then
      A%needs_x_tot = 1
    else
      A%needs_x_tot = 0
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_j_node_range

  subroutine set_diagonal_to_one_and_rest_of_row_to_zero( A, i)
    class(type_CSR_matrix_dp), intent(inout) :: A
    integer,                   intent(in   ) :: i
    call set_row_to_value(    A, i, 0._dp)
    call set_row_diag_to_val( A, i, 1._dp)
  end subroutine set_diagonal_to_one_and_rest_of_row_to_zero

  subroutine set_row_to_value( A, i, val)
    !< Set A(i,:) to val
    class(type_CSR_matrix_dp), intent(inout) :: A
    integer,                   intent(in   ) :: i
    real(dp),                  intent(in   ) :: val
    integer :: k
    do k = A%ptr( i), A%ptr( i+1) - 1
      A%val( k) = val
    end do
  end subroutine set_row_to_value

  subroutine set_row_diag_to_val( A, i, val)
    !< Set A(i,i) to val
    class(type_CSR_matrix_dp), intent(inout) :: A
    integer,                   intent(in   ) :: i
    real(dp),                  intent(in   ) :: val
    integer :: k
    do k = A%ptr( i), A%ptr( i+1) - 1
      if (A%ind( k) == i) then
        A%val( k) = val
      end if
    end do
  end subroutine set_row_diag_to_val

  subroutine save_CSR_matrix_as_NetCDF( A, output_dir, filename)
    !< Save a CSR matrix as a NetCDF file

    ! In- and output variables:
    class(type_CSR_matrix_dp), intent(in) :: A                !< The CSR matrix
    character(len=*),          intent(in) :: output_dir       !< The name of the output directory
    character(len=*),          intent(in) :: filename         !< The name of the file

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'save_CSR_matrix_as_NetCDF'
    character(len=:), allocatable :: filename_sans_nc, filename_with_nc, filename_including_output_dir
    integer                       :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    if (endswith( filename, '.nc', case_sensitive = .false.)) then
      filename_sans_nc = filename( 1: len_trim( filename)-3)
    else
      filename_sans_nc = filename( 1: len_trim( filename))
    end if
    filename_with_nc = filename_sans_nc // '.nc'

    filename_including_output_dir = trim( output_dir) // '/' // trim( filename_with_nc)

    call create_new_netcdf_file_for_writing( filename_including_output_dir, ncid)
    call write_CSR_matrix_to_NetCDF( A, filename_including_output_dir, ncid)
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine save_CSR_matrix_as_NetCDF

  subroutine write_CSR_matrix_to_NetCDF( A, filename, ncid)
    !< Write a CSR matrix to a NetCDF file (or a group therein)

    ! In- and output variables:
    class(type_CSR_matrix_dp), intent(in) :: A           !< The CSR matrix
    character(len=*),          intent(in) :: filename    !< The name of the file (only used for error messaging)
    integer,                   intent(in) :: ncid        !< ID of a file, or a group within a file

    ! Local variables:
    character(len=*), parameter :: routine_name = 'write_CSR_matrix_to_NetCDF'
    integer                     :: id_dim_m, id_dim_mp1, id_dim_n, id_dim_nnz, id_var_ptr, id_var_ind, id_var_val
    type(type_CSR_matrix_dp)    :: A_tot

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety
    if (.not. A%is_finalised) then
      call crash('A is not finalised')
    end if

    call A%gather_to_primary( A_tot)

    ! Create dimensions
    call create_dimension( filename, ncid, 'm'     , A_tot%m  , id_dim_m  )
    call create_dimension( filename, ncid, 'mplus1', A_tot%m+1, id_dim_mp1)
    call create_dimension( filename, ncid, 'n'     , A_tot%n  , id_dim_n  )
    call create_dimension( filename, ncid, 'nnz'   , A_tot%nnz, id_dim_nnz)

    ! Create variables
    call create_variable( filename, ncid, 'ptr', NF90_INT   , [id_dim_mp1], id_var_ptr)
    call create_variable( filename, ncid, 'ind', NF90_INT   , [id_dim_nnz], id_var_ind)
    call create_variable( filename, ncid, 'val', NF90_DOUBLE, [id_dim_nnz], id_var_val)

    ! Write variables
    call write_var_primary( filename, ncid, id_var_ptr, A_tot%ptr)
    call write_var_primary( filename, ncid, id_var_ind, A_tot%ind)
    call write_var_primary( filename, ncid, id_var_val, A_tot%val)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_CSR_matrix_to_NetCDF

  subroutine write_CSR_matrix_to_dist_NetCDFs( A, output_dir, filename_base)
    !< Write a CSR matrix to multiple NetCDF files (one for each process), including all the parallellisation info

    ! In- and output variables:
    class(type_CSR_matrix_dp), intent(in) :: A                !< The CSR matrix
    character(len=*),          intent(in) :: output_dir       !< The name of the output directory
    character(len=*),          intent(in) :: filename_base    !< The name of the file

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'write_CSR_matrix_to_dist_NetCDFs'
    character(len=:), allocatable :: filename_base_sans_nc, filename_template, filename
    integer                       :: ncid
    integer                       :: id_dim_m_locp1, id_dim_nnz
    integer                       :: id_var_ptr, id_var_ind, id_var_val

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Safety
    if (.not. A%is_finalised) then
      call crash('A is not finalised')
    end if

    ! Create a process-specific filename
    if (endswith( filename_base, '.nc', case_sensitive = .false.)) then
      filename_base_sans_nc = filename_base( 1: len_trim( filename_base)-3)
    else
      filename_base_sans_nc = filename_base( 1: len_trim( filename_base))
    end if

    filename_template = trim( filename_base_sans_nc) // '_proc_{proc}.nc'
    filename = trim( output_dir) // '/' // trim( insert_val_into_string_int( filename_template, '{proc}', par%i))

    ! Create a NetCDF file for each process
    call handle_netcdf_error( NF90_CREATE( filename, ior( NF90_NOCLOBBER, NF90_NETCDF4), ncid), filename = filename)

    ! Create dimensions and variables for the CSR matrix
    call handle_netcdf_error( NF90_DEF_DIM( ncid, 'm_locplus1', A%m_loc+1, id_dim_m_locp1), filename = filename, dimvarname = 'm_locplus1')
    call handle_netcdf_error( NF90_DEF_DIM( ncid, 'nnz'       , A%nnz    , id_dim_nnz    ), filename = filename, dimvarname = 'nnz')

    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'm'          , A%m          )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'n'          , A%n          )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'nnz'        , A%nnz        )

    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'm_loc'      , A%m_loc      )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'i1'         , A%i1         )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'i2'         , A%i2         )

    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'm_node'     , A%m_node     )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'i1_node'    , A%i1_node    )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'i2_node'    , A%i2_node    )

    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'n_loc'      , A%n_loc      )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'j1'         , A%j1         )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'j2'         , A%j2         )

    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'n_node'     , A%n_node     )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'j1_node'    , A%j1_node    )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'j2_node'    , A%j2_node    )

    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'j_min_node' , A%j_min_node )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'j_max_node' , A%j_max_node )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, 'needs_x_tot', A%needs_x_tot)

    call A%pai_x%setup_in_netcdf_file( filename, ncid, 'pai_x')
    call A%pai_y%setup_in_netcdf_file( filename, ncid, 'pai_y')

    call handle_netcdf_error( NF90_DEF_VAR( ncid, name = 'ptr', xtype = NF90_INT, &
      dimids = [id_dim_m_locp1], varid = id_var_ptr), filename = filename, dimvarname = 'ptr')
    call handle_netcdf_error( NF90_DEF_VAR( ncid, name = 'ind', xtype = NF90_INT, &
      dimids = [id_dim_nnz    ], varid = id_var_ind), filename = filename, dimvarname = 'ind')
    call handle_netcdf_error( NF90_DEF_VAR( ncid, name = 'val', xtype = NF90_DOUBLE, &
      dimids = [id_dim_nnz    ], varid = id_var_val), filename = filename, dimvarname = 'val')

    ! Write data to the NetCDF file
    call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var_ptr, A%ptr), filename = filename, dimvarname = 'ptr')
    call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var_ind, A%ind), filename = filename, dimvarname = 'ind')
    call handle_netcdf_error( NF90_PUT_VAR( ncid, id_var_val, A%val), filename = filename, dimvarname = 'val')

    call handle_netcdf_error( NF90_CLOSE( ncid))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_CSR_matrix_to_dist_NetCDFs

  subroutine read_CSR_matrix_from_dist_NetCDFs( A, foldername, filename_base)
    !< Read a CSR matrix from multiple NetCDF files (one for each process)

    ! In- and output variables:
    class(type_CSR_matrix_dp), intent(inout) :: A                !< The CSR matrix
    character(len=*),          intent(in   ) :: foldername       !< The name of the directory containing the NetCDF files
    character(len=*),          intent(in   ) :: filename_base    !< The base name of the file

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'read_CSR_matrix_from_dist_NetCDFs'
    character(len=:), allocatable :: filename_base_sans_nc, filename_template, filename
    logical                       :: file_exists
    integer                       :: ncid
    integer                       :: id_var_ptr, id_var_ind, id_var_val

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate any memory that was already allocated for A
    call A%deallocate

    ! Create a process-specific filename
    if (endswith( filename_base, '.nc', case_sensitive = .false.)) then
      filename_base_sans_nc = filename_base( 1: len_trim( filename_base)-3)
    else
      filename_base_sans_nc = filename_base( 1: len_trim( filename_base))
    end if

    filename_template = trim( filename_base_sans_nc) // '_proc_{proc}.nc'
    filename = trim( foldername) // '/' // trim( insert_val_into_string_int( filename_template, '{proc}', par%i))

    ! Check if this file actually exists
    inquire( exist = file_exists, file = trim( filename))
    if (.not. file_exists) call crash('file "' // trim( filename) // '" not found!')

    ! Open the NetCDF file with read-only access
    call handle_netcdf_error( NF90_OPEN( trim( filename), NF90_NOWRITE, ncid), filename = filename)

    ! Read matrix data
    call read_scalar_variable_dist_int( filename, ncid, 'm'          , A%m          )
    call read_scalar_variable_dist_int( filename, ncid, 'n'          , A%n          )
    call read_scalar_variable_dist_int( filename, ncid, 'nnz'        , A%nnz        )
    A%nnz_max = A%nnz

    call read_scalar_variable_dist_int( filename, ncid, 'm_loc'      , A%m_loc      )
    call read_scalar_variable_dist_int( filename, ncid, 'i1'         , A%i1         )
    call read_scalar_variable_dist_int( filename, ncid, 'i2'         , A%i2         )

    call read_scalar_variable_dist_int( filename, ncid, 'm_node'     , A%m_node     )
    call read_scalar_variable_dist_int( filename, ncid, 'i1_node'    , A%i1_node    )
    call read_scalar_variable_dist_int( filename, ncid, 'i2_node'    , A%i2_node    )

    call read_scalar_variable_dist_int( filename, ncid, 'n_loc'      , A%n_loc      )
    call read_scalar_variable_dist_int( filename, ncid, 'j1'         , A%j1         )
    call read_scalar_variable_dist_int( filename, ncid, 'j2'         , A%j2         )

    call read_scalar_variable_dist_int( filename, ncid, 'n_node'     , A%n_node     )
    call read_scalar_variable_dist_int( filename, ncid, 'j1_node'    , A%j1_node    )
    call read_scalar_variable_dist_int( filename, ncid, 'j2_node'    , A%j2_node    )

    call read_scalar_variable_dist_int( filename, ncid, 'j_min_node' , A%j_min_node )
    call read_scalar_variable_dist_int( filename, ncid, 'j_max_node' , A%j_max_node )
    call read_scalar_variable_dist_int( filename, ncid, 'needs_x_tot', A%needs_x_tot)

    call A%pai_x%read_from_netcdf_file( filename, ncid, 'pai_x')
    call A%pai_y%read_from_netcdf_file( filename, ncid, 'pai_y')

    ! Allocate memory for A
    allocate( A%ptr( A%i1: A%i2+1), source = 1)
    allocate( A%ind( A%nnz_max), source = 0    )
    allocate( A%val( A%nnz_max), source = 0._dp)

    ! Read matrix data from the NetCDF file
    call handle_netcdf_error( NF90_INQ_VARID( ncid, 'ptr', id_var_ptr))
    call handle_netcdf_error( NF90_GET_VAR( ncid, id_var_ptr, A%ptr), filename = filename, dimvarname = 'ptr')
    call handle_netcdf_error( NF90_INQ_VARID( ncid, 'ind', id_var_ind))
    call handle_netcdf_error( NF90_GET_VAR( ncid, id_var_ind, A%ind), filename = filename, dimvarname = 'ind')
    call handle_netcdf_error( NF90_INQ_VARID( ncid, 'val', id_var_val))
    call handle_netcdf_error( NF90_GET_VAR( ncid, id_var_val, A%val), filename = filename, dimvarname = 'val')

    call handle_netcdf_error( NF90_CLOSE( ncid))

    A%is_finalised = .true.

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_CSR_matrix_from_dist_NetCDFs

  function test_CSR_matrix_equality( A, B) result( res)

    class(type_CSR_matrix_dp), intent(in) :: A
    class(type_CSR_matrix_dp), intent(in) :: B
    logical                               :: res

    integer :: ierr

    res = &
      A%is_finalised .and. &
      B%is_finalised .and. &
      A%m       == B%m       .and. &
      A%n       == B%n       .and. &
      A%nnz_max == B%nnz_max .and. &
      A%nnz     == B%nnz     .and. &
      A%m_loc   == B%m_loc   .and. &
      A%i1      == B%i1      .and. &
      A%i2      == B%i2      .and. &
      A%m_node  == B%m_node  .and. &
      A%i1_node == B%i1_node .and. &
      A%i2_node == B%i2_node .and. &
      A%n_loc   == B%n_loc   .and. &
      A%j1      == B%j1      .and. &
      A%j2      == B%j2      .and. &
      A%n_node  == B%n_node  .and. &
      A%j1_node == B%j1_node .and. &
      A%j2_node == B%j2_node .and. &
      A%j_min_node == B%j_min_node .and. &
      A%j_max_node == B%j_max_node .and. &
      A%needs_x_tot == B%needs_x_tot .and. &
      A%pai_x == B%pai_x .and. &
      A%pai_y == B%pai_y .and. &
      all( A%ptr == B%ptr) .and. &
      all( A%ind == B%ind) .and. &
      all( A%val == B%val)

  end function test_CSR_matrix_equality

  ! ===== CSR matrices in local memory =====
  ! ========================================

  subroutine allocate_matrix_CSR_loc( A, m, n, nnz_max)
    ! Allocate memory for a CSR-format sparse m-by-n matrix A

    ! In- and output variables:
    class(type_CSR_matrix_dp), intent(inout) :: A
    integer,                   intent(in   ) :: m, n, nnz_max

    ! Local variables:
    character(*), parameter:: routine_name = 'allocate_matrix_CSR_loc'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Matrix dimensions
    A%m       = m
    A%n       = n
    A%m_loc   = m
    A%n_loc   = n
    A%nnz_max = nnz_max
    A%nnz     = 0

    A%i1 = 1
    A%i2 = m
    A%j1 = 1
    A%j2 = n

    ! Allocate memory
    allocate( A%ptr( A%m+1    ), source = 1    )
    allocate( A%ind( A%nnz_max), source = 0    )
    allocate( A%val( A%nnz_max), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_matrix_CSR_loc

end module CSR_matrix_mod
