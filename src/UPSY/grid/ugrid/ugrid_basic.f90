module ugrid_basic

  use precisions, only: dp
  use ugrid_data, only: atype_ugrid_data

  implicit none

  private

  public :: atype_ugrid

  type, extends(atype_ugrid_data) :: atype_ugrid
  end type atype_ugrid

end module ugrid_basic
