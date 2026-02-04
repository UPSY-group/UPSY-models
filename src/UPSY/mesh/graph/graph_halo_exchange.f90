module graph_halo_exchange

  ! See graph_parallelisation for an explanation of the halos.

  use precisions, only: dp
  use graph_types, only: type_graph
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mpi_basic, only: par, sync
  use halo_exchange_mod, only: basic_halo_exchange

  implicit none

  private

  public :: exchange_halos

  interface exchange_halos
    procedure :: exchange_halos_logical
    procedure :: exchange_halos_int
    procedure :: exchange_halos_int_3D
    procedure :: exchange_halos_dp
    procedure :: exchange_halos_dp_3D
  end interface exchange_halos

contains

  subroutine exchange_halos_logical( graph, d_nih)

    ! In/output variables:
    type(type_graph),      intent(in   ) :: graph
    logical, dimension(:), intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_logical'

    ! Add routine to path
    call init_routine( routine_name)

    ! If running on one node, do nothing
    if (par%n_nodes == 1) then
      call sync
      call finalise_routine( routine_name)
      return
    end if

    call basic_halo_exchange( graph%pai, d_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_logical

  subroutine exchange_halos_int( graph, d_nih)

    ! In/output variables:
    type(type_graph),      intent(in   ) :: graph
    integer, dimension(:), intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_int'

    ! Add routine to path
    call init_routine( routine_name)

    ! If running on one node, do nothing
    if (par%n_nodes == 1) then
      call sync
      call finalise_routine( routine_name)
      return
    end if

    call basic_halo_exchange( graph%pai, d_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_int

  subroutine exchange_halos_int_3D( graph, d_nih)

    ! In/output variables:
    type(type_graph),        intent(in   ) :: graph
    integer, dimension(:,:), intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_int_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! If running on one node, do nothing
    if (par%n_nodes == 1) then
      call sync
      call finalise_routine( routine_name)
      return
    end if

    call basic_halo_exchange( graph%pai, size( d_nih,2), d_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_int_3D

  subroutine exchange_halos_dp( graph, d_nih)

    ! In/output variables:
    type(type_graph),       intent(in   ) :: graph
    real(dp), dimension(:), intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_dp'

    ! Add routine to path
    call init_routine( routine_name)

    ! If running on one node, do nothing
    if (par%n_nodes == 1) then
      call sync
      call finalise_routine( routine_name)
      return
    end if

    call basic_halo_exchange( graph%pai, d_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_dp

  subroutine exchange_halos_dp_3D( graph, d_nih)

    ! In/output variables:
    type(type_graph),         intent(in   ) :: graph
    real(dp), dimension(:,:), intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'exchange_halos_dp_3D'

    ! Add routine to path
    call init_routine( routine_name)

    ! If running on one node, do nothing
    if (par%n_nodes == 1) then
      call sync
      call finalise_routine( routine_name)
      return
    end if

    call basic_halo_exchange( graph%pai, size( d_nih,2), d_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine exchange_halos_dp_3D

end module graph_halo_exchange
