module ocean_ismip

  use precisions, only: dp
  use UPSY_main, only: UPSY
  use parameters, only: NaN
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ocean_model_types, only: type_ocean_model, type_ocean_model_ismip
  use netcdf_io_main
  use mpi_f08, only: MPI_BCAST, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD

  implicit none

  private

  public :: initialise_ocean_model_ismip, run_ocean_model_ismip

contains

  subroutine run_ocean_model_ismip( mesh, ocean, time)

    ! In/output variables:
    type(type_mesh),        intent(in   ) :: mesh
    type(type_ocean_model), intent(inout) :: ocean
    real(dp),               intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_ocean_model_ismip'
    real(dp)                       :: w0, w1

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_ocean_model_ismip

  subroutine initialise_ocean_model_ismip( mesh, ismip)

    ! In/output variables:
    type(type_mesh),              intent(in   ) :: mesh
    type(type_ocean_model_ismip), intent(inout) :: ismip

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_ocean_model_ismip'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate memory

    allocate (ismip%T_t0( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate (ismip%S_t0( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)

    allocate (ismip%T_t1( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate (ismip%S_t1( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model_ismip


end module ocean_ismip
