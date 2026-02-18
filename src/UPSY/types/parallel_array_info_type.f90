module parallel_array_info_type

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

end module parallel_array_info_type
