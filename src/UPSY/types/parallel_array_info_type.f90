module parallel_array_info_type

  use netcdf_basic_wrappers, only: handle_netcdf_error
  use netcdf, only: NF90_DEF_VAR, NF90_INT
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use netcdf_basic_wrappers, only: create_and_write_to_scalar_variable_dist_int, &
    read_scalar_variable_dist_int
  use crash_mod, only: crash
  use precisions, only: dp

  implicit none

  private

  public :: type_par_arr_info

  type type_par_arr_info

    integer :: n                             ! Global number of elements
    integer :: i1,      i2,       n_loc      ! Range owned by this process
    integer :: i1_node, i2_node,  n_node     ! Range owned by this shared-memory node
    integer :: i1_nih,  i2_nih,   n_nih      ! Size of shared-memory array on this node, including exterior halos
    integer :: i1_hle,  i2_hle,   n_hle      ! Range of left  exterior halo
    integer :: i1_hli,  i2_hli,   n_hli      ! Range of left  interior halo
    integer :: i1_hre,  i2_hre,   n_hre      ! Range of right exterior halo
    integer :: i1_hri,  i2_hri,   n_hri      ! Range of right interior halo

  contains

    generic,   public  :: operator(==) => eq
    procedure, private :: eq => test_pai_equality

    procedure, public  :: setup_in_netcdf_file
    procedure, public  :: read_from_netcdf_file

    generic,   public  :: is_hybrid => &
      is_hybrid_logical_2D, &
      is_hybrid_logical_3D, &
      is_hybrid_int_2D, &
      is_hybrid_int_3D, &
      is_hybrid_dp_2D, &
      is_hybrid_dp_3D
    procedure, private :: is_hybrid_logical_2D
    procedure, private :: is_hybrid_logical_3D
    procedure, private :: is_hybrid_int_2D
    procedure, private :: is_hybrid_int_3D
    procedure, private :: is_hybrid_dp_2D
    procedure, private :: is_hybrid_dp_3D

  end type type_par_arr_info

contains

  pure function test_pai_equality( pai1, pai2) result( res)
    class(type_par_arr_info), intent(in) :: pai1
    class(type_par_arr_info), intent(in) :: pai2
    logical                                 :: res
    res = &
      pai1%n       == pai2%n       .and. &
      pai1%i1      == pai2%i1      .and. &
      pai1%i2      == pai2%i2      .and. &
      pai1%n_loc   == pai2%n_loc   .and. &
      pai1%i1_node == pai2%i1_node .and. &
      pai1%i2_node == pai2%i2_node .and. &
      pai1%n_node  == pai2%n_node  .and. &
      pai1%i1_nih  == pai2%i1_nih  .and. &
      pai1%i2_nih  == pai2%i2_nih  .and. &
      pai1%n_nih   == pai2%n_nih   .and. &
      pai1%i1_hle  == pai2%i1_hle  .and. &
      pai1%i2_hle  == pai2%i2_hle  .and. &
      pai1%n_hle   == pai2%n_hle   .and. &
      pai1%i1_hli  == pai2%i1_hli  .and. &
      pai1%i2_hli  == pai2%i2_hli  .and. &
      pai1%n_hli   == pai2%n_hli   .and. &
      pai1%i1_hre  == pai2%i1_hre  .and. &
      pai1%i2_hre  == pai2%i2_hre  .and. &
      pai1%n_hre   == pai2%n_hre   .and. &
      pai1%i1_hri  == pai2%i1_hri  .and. &
      pai1%i2_hri  == pai2%i2_hri  .and. &
      pai1%n_hri   == pai2%n_hri
  end function test_pai_equality

  subroutine setup_in_netcdf_file( self, filename, ncid, pai_name)
    !< Write parallel array info to an existing, open NetCDF file

    ! In/output variables:
    class(type_par_arr_info), intent(in) :: self
    character(len=*),         intent(in) :: filename
    integer,                  intent(in) :: ncid
    character(len=*),         intent(in) :: pai_name

    ! Local variables:

    ! Local variables:
    character(len=*), parameter :: routine_name = 'setup_parallel_array_info_in_netcdf_file'

    ! Add routine to call stack
    call init_routine( routine_name)

    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n'      , self%n      )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1'     , self%i1     )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2'     , self%i2     )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_loc'  , self%n_loc  )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1_node', self%i1_node)
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2_node', self%i2_node)
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_node' , self%n_node )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1_nih' , self%i1_nih )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2_nih' , self%i2_nih )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_nih'  , self%n_nih  )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1_hle' , self%i1_hle )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2_hle' , self%i2_hle )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_hle'  , self%n_hle  )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1_hli' , self%i1_hli )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2_hli' , self%i2_hli )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_hli'  , self%n_hli  )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1_hre' , self%i1_hre )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2_hre' , self%i2_hre )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_hre'  , self%n_hre  )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1_hri' , self%i1_hri )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2_hri' , self%i2_hri )
    call create_and_write_to_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_hri'  , self%n_hri  )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_in_netcdf_file

  subroutine read_from_netcdf_file( self, filename, ncid, pai_name)
    !< Read parallel array info from a NetCDF file

    ! In/output variables:
    class(type_par_arr_info), intent(  out) :: self
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    character(len=*),         intent(in   ) :: pai_name

    ! Local variables:

    ! Local variables:
    character(len=*), parameter :: routine_name = 'read_parallel_array_info_from_netcdf_file'

    ! Add routine to call stack
    call init_routine( routine_name)

    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n'      , self%n      )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1'     , self%i1     )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2'     , self%i2     )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_loc'  , self%n_loc  )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1_node', self%i1_node)
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2_node', self%i2_node)
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_node' , self%n_node )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1_nih' , self%i1_nih )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2_nih' , self%i2_nih )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_nih'  , self%n_nih  )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1_hle' , self%i1_hle )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2_hle' , self%i2_hle )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_hle'  , self%n_hle  )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1_hli' , self%i1_hli )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2_hli' , self%i2_hli )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_hli'  , self%n_hli  )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1_hre' , self%i1_hre )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2_hre' , self%i2_hre )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_hre'  , self%n_hre  )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i1_hri' , self%i1_hri )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_i2_hri' , self%i2_hri )
    call read_scalar_variable_dist_int( filename, ncid, trim( pai_name) // '_n_hri'  , self%n_hri  )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_file

  function is_hybrid_logical_2D( self, d) result( isso)
    class(type_par_arr_info), intent(in) :: self
    logical, dimension(:),    intent(in) :: d
    logical :: isso
    isso = .false.
    if (size( d,1) == self%n_loc) then
      isso = .false.
    elseif (size( d,1) == self%n_nih) then
      isso = .true.
    else
      call crash('array doesnt meet expected size for either distributed or hybrid distributed/shared memory')
    end if
  end function is_hybrid_logical_2D

  function is_hybrid_logical_3D( self, d) result( isso)
    class(type_par_arr_info), intent(in) :: self
    logical, dimension(:,:),  intent(in) :: d
    logical :: isso
    isso = .false.
    if (size( d,1) == self%n_loc) then
      isso = .false.
    elseif (size( d,1) == self%n_nih) then
      isso = .true.
    else
      call crash('array doesnt meet expected size for either distributed or hybrid distributed/shared memory')
    end if
  end function is_hybrid_logical_3D

  function is_hybrid_int_2D( self, d) result( isso)
    class(type_par_arr_info), intent(in) :: self
    integer, dimension(:),    intent(in) :: d
    logical :: isso
    isso = .false.
    if (size( d,1) == self%n_loc) then
      isso = .false.
    elseif (size( d,1) == self%n_nih) then
      isso = .true.
    else
      call crash('array doesnt meet expected size for either distributed or hybrid distributed/shared memory')
    end if
  end function is_hybrid_int_2D

  function is_hybrid_int_3D( self, d) result( isso)
    class(type_par_arr_info), intent(in) :: self
    integer, dimension(:,:),  intent(in) :: d
    logical :: isso
    isso = .false.
    if (size( d,1) == self%n_loc) then
      isso = .false.
    elseif (size( d,1) == self%n_nih) then
      isso = .true.
    else
      call crash('array doesnt meet expected size for either distributed or hybrid distributed/shared memory')
    end if
  end function is_hybrid_int_3D

  function is_hybrid_dp_2D( self, d) result( isso)
    class(type_par_arr_info), intent(in) :: self
    real(dp), dimension(:),   intent(in) :: d
    logical :: isso
    isso = .false.
    if (size( d,1) == self%n_loc) then
      isso = .false.
    elseif (size( d,1) == self%n_nih) then
      isso = .true.
    else
      call crash('array doesnt meet expected size for either distributed or hybrid distributed/shared memory')
    end if
  end function is_hybrid_dp_2D

  function is_hybrid_dp_3D( self, d) result( isso)
    class(type_par_arr_info), intent(in) :: self
    real(dp), dimension(:,:), intent(in) :: d
    logical :: isso
    isso = .false.
    if (size( d,1) == self%n_loc) then
      isso = .false.
    elseif (size( d,1) == self%n_nih) then
      isso = .true.
    else
      call crash('array doesnt meet expected size for either distributed or hybrid distributed/shared memory')
    end if
  end function is_hybrid_dp_3D

end module parallel_array_info_type
