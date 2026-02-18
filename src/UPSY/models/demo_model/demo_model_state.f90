module demo_model_state

  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_demo_model_state

  type :: type_demo_model_state

      ! Some ice-model-esque data fields
      real(dp), dimension(:  ), contiguous, pointer :: H        => null()
      real(dp), dimension(:,:), contiguous, pointer :: u_3D     => null()
      real(dp), dimension(:,:), contiguous, pointer :: v_3D     => null()
      logical,  dimension(:  ), contiguous, pointer :: mask_ice => null()
      real(dp), dimension(:,:), contiguous, pointer :: T2m      => null()
      type(MPI_WIN) :: wH, wu_3D, wv_3D, wmask_ice, wT2m

  end type type_demo_model_state

end module demo_model_state