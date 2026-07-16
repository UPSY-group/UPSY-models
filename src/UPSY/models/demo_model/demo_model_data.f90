module demo_model_data

  use models_basic, only: atype_model
  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_demo_model_data

  type, abstract, extends(atype_model) :: atype_demo_model_data
    !< Variables that are common to all demo models,
    !< which we want other models to be able to access

    ! Some ice-model-esque data fields
    real(dp), dimension(:  ), contiguous, pointer :: H        => null()
    real(dp), dimension(:,:), contiguous, pointer :: u_3D     => null()
    real(dp), dimension(:,:), contiguous, pointer :: v_3D     => null()
    logical,  dimension(:  ), contiguous, pointer :: mask_ice => null()
    real(dp), dimension(:,:), contiguous, pointer :: T2m      => null()
    type(MPI_WIN) :: wH, wu_3D, wv_3D, wmask_ice, wT2m

  end type atype_demo_model_data

end module demo_model_data