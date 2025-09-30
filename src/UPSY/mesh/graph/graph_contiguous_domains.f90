module graph_contiguous_domains

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use graph_types, only: type_graph
  use sorting, only: quick_n_dirty_sort

  implicit none

  private

  public :: enforce_contiguous_process_domains_graph

contains

  !> Shuffle nodes so that the nodes owned by each process form
  !> a contiguous domain with short borders, to minimise the data volumes for halo exchanges.
  subroutine enforce_contiguous_process_domains_graph( graph)

    ! In/output variables:
    type(type_graph), intent(inout) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'enforce_contiguous_process_domains_graph'
    real(dp), dimension(graph%n)               :: xx
    integer,  dimension(graph%n)               :: ni_new2ni_old, ni_old2ni_new
    integer                                    :: ni_old, ni_new, ci, nj_old, nj_new, ni, vi, ti, ei
    real(dp), dimension(graph%n,2)             :: V_old
    integer,  dimension(graph%n)               :: nC_old
    integer,  dimension(graph%n, graph%nC_mem) :: C_old
    integer,  dimension(graph%n)               :: ni_old2vi, ni_old2ti, ni_old2ei
    integer,  dimension(size(graph%vi2ni))     :: vi2ni_old
    integer,  dimension(size(graph%ti2ni))     :: ti2ni_old
    integer,  dimension(size(graph%ei2ni))     :: ei2ni_old
    logical,  dimension(graph%n)               :: is_border_old
    real(dp), dimension(graph%n,2)             :: border_nhat_old

    ! Add routine to path
    call init_routine( routine_name)

    ! Sort nodes by x-coordinate
    xx = graph%V( :,1)
    call quick_n_dirty_sort( xx, ni_new2ni_old)

    ! Calculate translation table in the opposite direction
    do ni_new = 1, graph%n
      ni_old = ni_new2ni_old( ni_new)
      ni_old2ni_new( ni_old) = ni_new
    end do

    ! Save old translation tables
    ni_old2vi       = graph%ni2vi
    ni_old2ti       = graph%ni2ti
    ni_old2ei       = graph%ni2ei
    vi2ni_old       = graph%vi2ni
    ti2ni_old       = graph%ti2ni
    ei2ni_old       = graph%ei2ni
    V_old           = graph%V
    nC_old          = graph%nC
    C_old           = graph%C
    is_border_old   = graph%is_border
    border_nhat_old = graph%border_nhat

    ! Reset node data
    graph%ni2vi       = 0
    graph%ni2ti       = 0
    graph%ni2ei       = 0
    graph%vi2ni       = 0
    graph%ti2ni       = 0
    graph%ei2ni       = 0
    graph%V           = 0._dp
    graph%nC          = 0
    graph%C           = 0
    graph%is_border   = .false.
    graph%border_nhat = 0._dp

    ! Fill in new node data
    do ni_new = 1, graph%n

      ! This new node corresponds to this old node
      ni_old = ni_new2ni_old( ni_new)

      ! ni2vi, vi2ni
      vi = ni_old2vi( ni_old)
      graph%ni2vi( ni_new) = vi
      if (vi > 0) graph%vi2ni( vi) = ni_new

      ! ni2ti, ti2ni
      ti = ni_old2ti( ni_old)
      graph%ni2ti( ni_new) = ti
      if (ti > 0) graph%ti2ni( ti) = ni_new

      ! ni2ei, ei2ni
      ei = ni_old2ei( ni_old)
      graph%ni2ei( ni_new) = ei
      if (ei > 0) graph%ei2ni( ei) = ni_new

      ! V
      graph%V( ni_new,:) = V_old( ni_old,:)

      ! nC
      graph%nC( ni_new) = nC_old( ni_old)

      ! C
      do ci = 1, nC_old( ni_old)
        nj_old = C_old( ni_old, ci)
        nj_new = ni_old2ni_new( nj_old)
        graph%C( ni_new,ci) = nj_new
      end do

      ! is_border
      graph%is_border( ni_new) = is_border_old( ni_old)

      ! border_nhat
      graph%border_nhat( ni_new,:) = border_nhat_old( ni_old,:)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine enforce_contiguous_process_domains_graph

end module graph_contiguous_domains