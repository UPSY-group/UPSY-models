module SMB_basic
  !< The basic surface mass balance model

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mpi_f08, only: MPI_WIN
  use UPSY_main, only: atype_model
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid

  implicit none

  private

  public :: atype_SMB_model

  type, abstract, extends(atype_model) :: atype_SMB_model
    !< The basic surface mass balance model

      ! Main data fields
      real(dp), dimension(:), contiguous, pointer :: SMB   !< [m] Net annual  SMB
      type(MPI_WIN) :: wSMB

      ! Timestepping
      real(dp) :: t_next

    contains

      private

      procedure, public :: init_common

  end type atype_SMB_model

contains

  subroutine init_common( self, mesh)

    ! In/output variables:
    class(atype_SMB_model), intent(inout) :: self
    type(type_mesh),        intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_SMB_model_init_common'

    ! Add routine to path
    call init_routine( routine_name)

    call self%set_grid( mesh)

    call self%create_field( self%SMB, self%wSMB, &
      mesh, Arakawa_grid%a(), &
      name      = 'SMB', &
      long_name = 'surface mass balance', &
      units     = 'm yr^-1')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine init_common

end module SMB_basic