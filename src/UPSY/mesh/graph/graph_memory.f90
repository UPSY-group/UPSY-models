module graph_memory

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use graph_types, only: type_graph
  use reallocate_mod, only: reallocate
  use mpi_distributed_shared_memory, only: deallocate_dist_shared

  implicit none

  private

  public :: allocate_graph_primary, crop_graph_primary, deallocate_graph

contains

  subroutine allocate_graph_primary( graph, n, nV, nTri, nE, nC_mem)

    ! In/output variables:
    type(type_graph), intent(inout) :: graph
    integer,          intent(in   ) :: n, nV, nTri, nE, nC_mem

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'allocate_graph_primary'

    ! Add routine to path
    call init_routine( routine_name)

    graph%n      = 0
    graph%nC_mem = nC_mem

    ! Mapping between parent mesh and graph
    allocate( graph%ni2vi     ( n        ), source = 0)
    allocate( graph%ni2ti     ( n        ), source = 0)
    allocate( graph%ni2ei     ( n        ), source = 0)
    allocate( graph%vi2ni     ( nV       ), source = 0)
    allocate( graph%ti2ni     ( nTri     ), source = 0)
    allocate( graph%ei2ni     ( nE       ), source = 0)

    ! Node coordinates and connectivity
    allocate( graph%V         ( n, 2     ), source = 0._dp)
    allocate( graph%nC        ( n        ), source = 0)
    allocate( graph%C         ( n, nC_mem), source = 0)

    ! border nodes
    allocate( graph%is_border  ( n        ), source = .false.)
    allocate( graph%border_nhat( n, 2     ), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_graph_primary

  subroutine crop_graph_primary( graph)

    ! In/output variables:
    type(type_graph), intent(inout) :: graph

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'crop_graph_primary'

    ! Add routine to path
    call init_routine( routine_name)

    ! Mapping between parent mesh and graph
    call reallocate( graph%ni2vi, graph%n)
    call reallocate( graph%ni2ti, graph%n)
    call reallocate( graph%ni2ei, graph%n)

    ! Node coordinates and connectivity
    call reallocate( graph%V         , graph%n, 2           )
    call reallocate( graph%nC        , graph%n              )
    call reallocate( graph%C         , graph%n, graph%nC_mem)

    ! border nodes
    call reallocate( graph%is_border  , graph%n              )
    call reallocate( graph%border_nhat, graph%n, 2           )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine crop_graph_primary

  subroutine deallocate_graph( graph)

    ! In/output variables:
    type(type_graph), intent(inout) :: graph

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'deallocate_graph'

    ! Add routine to path
    call init_routine( routine_name)

    ! Mapping between parent mesh and graph
    deallocate( graph%ni2vi)
    deallocate( graph%ni2ti)
    deallocate( graph%ni2ei)
    deallocate( graph%vi2ni)
    deallocate( graph%ti2ni)
    deallocate( graph%ei2ni)

    ! Node coordinates and connectivity
    deallocate( graph%V )
    deallocate( graph%nC)
    deallocate( graph%C )

    ! border nodes
    deallocate( graph%is_border  )
    deallocate( graph%border_nhat)

    ! Parallelisation ranges
    deallocate( graph%owning_process)
    deallocate( graph%owning_node   )

    call deallocate_dist_shared( graph%buffer1_g_nih, graph%wbuffer1_g_nih)
    call deallocate_dist_shared( graph%buffer2_g_nih, graph%wbuffer2_g_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_graph

end module graph_memory