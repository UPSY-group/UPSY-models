module graph_parallelisation

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use graph_types, only: type_graph
  use mesh_parallelisation, only: determine_ownership_ranges_equal
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared
  use mpi_basic, only: par, sync
  use mpi_f08, only: MPI_ALLGATHER, MPI_INTEGER, MPI_COMM_WORLD

  implicit none

  private

  public :: setup_graph_parallelisation

contains

  subroutine setup_graph_parallelisation( graph, nz)
    !< Setup the parallelisation of memory and halos on them graph

    ! In/output variables:
    type(type_graph), intent(inout) :: graph
    integer,          intent(in   ) :: nz

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_graph_parallelisation'

    ! Add routine to path
    call init_routine( routine_name)

    if (allocated( graph%owning_process)) deallocate( graph%owning_process)
    if (allocated( graph%owning_node   )) deallocate( graph%owning_node   )

    allocate( graph%owning_process( graph%n))
    allocate( graph%owning_node   ( graph%n))

    call determine_ownership_ranges_equal( graph%n, graph%ni1, graph%ni2, graph%n_loc, &
      graph%pai%i1_node, graph%pai%i2_node, graph%pai%n_node, &
      graph%owning_process, graph%owning_node)

    graph%pai%i1    = graph%ni1
    graph%pai%i2    = graph%ni2
    graph%pai%n_loc = graph%n_loc

    ! Determine all halos
    call determine_halos( graph)

    ! Allocate buffer shared memory for e.g. matrix multiplications
    if (associated( graph%buffer1_g_nih )) call deallocate_dist_shared( graph%buffer1_g_nih , graph%wbuffer1_g_nih )
    if (associated( graph%buffer2_g_nih )) call deallocate_dist_shared( graph%buffer2_g_nih , graph%wbuffer2_g_nih )
    if (associated( graph%buffer1_gk_nih)) call deallocate_dist_shared( graph%buffer1_gk_nih, graph%wbuffer1_gk_nih)
    if (associated( graph%buffer2_gk_nih)) call deallocate_dist_shared( graph%buffer2_gk_nih, graph%wbuffer2_gk_nih)
    call allocate_dist_shared( graph%buffer1_g_nih , graph%wbuffer1_g_nih , graph%pai%n_nih    )
    call allocate_dist_shared( graph%buffer2_g_nih , graph%wbuffer2_g_nih , graph%pai%n_nih    )
    call allocate_dist_shared( graph%buffer1_gk_nih, graph%wbuffer1_gk_nih, graph%pai%n_nih, nz)
    call allocate_dist_shared( graph%buffer2_gk_nih, graph%wbuffer2_gk_nih, graph%pai%n_nih, nz)
    graph%buffer1_g_nih(  graph%pai%i1_nih:graph%pai%i2_nih      ) => graph%buffer1_g_nih
    graph%buffer2_g_nih(  graph%pai%i1_nih:graph%pai%i2_nih      ) => graph%buffer2_g_nih
    graph%buffer1_gk_nih( graph%pai%i1_nih:graph%pai%i2_nih, 1:nz) => graph%buffer1_gk_nih
    graph%buffer2_gk_nih( graph%pai%i1_nih:graph%pai%i2_nih, 1:nz) => graph%buffer2_gk_nih

    ! call print_parallelisation_info( graph)

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 4)

  end subroutine setup_graph_parallelisation

  subroutine determine_halos( graph)
    !< Determine all the halos

    ! In/output variables:
    type(type_graph), intent(inout) :: graph

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'determine_halos'
    integer                        :: node_ID_left, node_ID_right

    ! Add routine to path
    call init_routine( routine_name)

    ! Left ("west")
    if (par%node_ID == 0) then
      ! There is no node to the left

      graph%pai%i1_nih  = graph%pai%i1_node
      graph%pai%n_hle   = 0
      graph%pai%i1_hle  =  0
      graph%pai%i2_hle  = -1
      graph%pai%n_hli   = 0
      graph%pai%i1_hli  =  0
      graph%pai%i2_hli  = -1

    else
      node_ID_left = par%node_ID - 1

      call determine_halo_range_a( graph, node_ID_left, par%node_ID , graph%pai%i1_hle, graph%pai%i2_hle)
      call determine_halo_range_a( graph, par%node_ID , node_ID_left, graph%pai%i1_hli, graph%pai%i2_hli)
      graph%pai%n_hle = graph%pai%i2_hle + 1 - graph%pai%i1_hle
      graph%pai%n_hli = graph%pai%i2_hli + 1 - graph%pai%i1_hli
      graph%pai%i1_nih = graph%pai%i1_hle

    end if

    ! Right ("east")
    if (par%node_ID == par%n_nodes-1) then
      ! There is no node to the right

      graph%pai%i2_nih  = graph%pai%i2_node
      graph%pai%n_hre   = 0
      graph%pai%i1_hre  =  0
      graph%pai%i2_hre  = -1
      graph%pai%n_hri   = 0
      graph%pai%i1_hri  =  0
      graph%pai%i2_hri  = -1

    else
      node_ID_right = par%node_ID + 1

      call determine_halo_range_a( graph, node_ID_right, par%node_ID  , graph%pai%i1_hre, graph%pai%i2_hre)
      call determine_halo_range_a( graph, par%node_ID  , node_ID_right, graph%pai%i1_hri, graph%pai%i2_hri)
      graph%pai%n_hre  = graph%pai%i2_hre + 1 - graph%pai%i1_hre
      graph%pai%n_hri  = graph%pai%i2_hri + 1 - graph%pai%i1_hri
      graph%pai%i2_nih = graph%pai%i2_hre

    end if

    graph%pai%n     = graph%n
    graph%pai%n_nih = graph%pai%i2_nih + 1 - graph%pai%i1_nih

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_halos

  subroutine determine_halo_range_a( graph, node_in, node_next_to, ni1, ni2)
    !< Find all graph nodes that are owned by HPC node_in but adjacent to HPC node_next_to,
    !< and return their range ni1-ni2

    ! In/output variables:
    type(type_graph), intent(in   ) :: graph
    integer,          intent(in   ) :: node_in, node_next_to
    integer,          intent(  out) :: ni1, ni2

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'determine_halo_range_a'
    integer                        :: ni, ci, nj
    logical                        :: has_neighbours_in_node_next_to

    ! Add routine to path
    call init_routine( routine_name)

    ni1 = graph%n
    ni2 = 1
    do ni = 1, graph%n
      if (graph%owning_node( ni) == node_in) then
        ! This graph node is owned by HPC node_in


        ! Check if this graph node has any neighbouring graph nodes owned by HPC node_next_to
        has_neighbours_in_node_next_to  = .false.
        do ci = 1, graph%nC( ni)
          nj = graph%C( ni,ci)
          if (graph%owning_node( nj) == node_next_to) then
            has_neighbours_in_node_next_to = .true.
            exit
          end if
        end do

        if (has_neighbours_in_node_next_to) then
          ni1 = min( ni1, ni)
          ni2 = max( ni2, ni)
        end if

      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_halo_range_a

  subroutine print_parallelisation_info( graph)

    ! In/output variables:
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'print_parallelisation_info'
    integer, dimension(0:par%n-1) :: node, process, &
      ni1, ni2, ni1_node, ni2_node, &
      ni1_nih, ni2_nih, ni1_hle, ni2_hle, n_hle, ni1_hli, ni2_hli, n_hli, ni1_hre, ni2_hre, n_hre, ni1_hri, ni2_hri, n_hri
    integer :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    call MPI_ALLGATHER( par%node_ID  , 1, MPI_INTEGER, node    , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( par%i        , 1, MPI_INTEGER, process , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER( graph%ni1     , 1, MPI_INTEGER, ni1     , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%ni2     , 1, MPI_INTEGER, ni2     , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    call MPI_ALLGATHER( graph%pai%i1_node, 1, MPI_INTEGER, ni1_node, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%i2_node, 1, MPI_INTEGER, ni2_node, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%i1_nih , 1, MPI_INTEGER, ni1_nih , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%i2_nih , 1, MPI_INTEGER, ni2_nih , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%n_hle  , 1, MPI_INTEGER, n_hle   , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%i1_hle , 1, MPI_INTEGER, ni1_hle , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%i2_hle , 1, MPI_INTEGER, ni2_hle , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%n_hli  , 1, MPI_INTEGER, n_hli   , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%i1_hli , 1, MPI_INTEGER, ni1_hli , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%i2_hli , 1, MPI_INTEGER, ni2_hli , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%n_hre  , 1, MPI_INTEGER, n_hre   , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%i1_hre , 1, MPI_INTEGER, ni1_hre , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%i2_hre , 1, MPI_INTEGER, ni2_hre , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%n_hri  , 1, MPI_INTEGER, n_hri   , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%i1_hri , 1, MPI_INTEGER, ni1_hri , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    call MPI_ALLGATHER( graph%pai%i2_hri , 1, MPI_INTEGER, ni2_hri , 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    if (par%primary) then
      write(0,'(A,7I8)') 'HPC node     :', node
      write(0,'(A,7I8)') 'HPC process  :', process
      write(0,*) ''
      write(0,*) '-= graph nodes =-'
      write(0,*) ''
      write(0,'(A,I8)')  'n = ', graph%n
      write(0,*) ''
      write(0,'(A,7I8)') 'ni1      :', ni1
      write(0,'(A,7I8)') 'ni2      :', ni2
      write(0,'(A,7I8)') 'ni1_node :', ni1_node
      write(0,'(A,7I8)') 'ni2_node :', ni2_node
      write(0,'(A,7I8)') 'ni1_nih  :', ni1_nih
      write(0,'(A,7I8)') 'ni2_nih  :', ni2_nih
      write(0,*) ''
      write(0,'(A,7I8)') 'n_hle    :', n_hle
      write(0,'(A,7I8)') 'ni1_hle  :', ni1_hle
      write(0,'(A,7I8)') 'ni2_hle  :', ni2_hle
      write(0,*) ''
      write(0,'(A,7I8)') 'n_hli    :', n_hli
      write(0,'(A,7I8)') 'ni1_hli  :', ni1_hli
      write(0,'(A,7I8)') 'ni2_hli  :', ni2_hli
      write(0,*) ''
      write(0,'(A,7I8)') 'n_hre    :', n_hre
      write(0,'(A,7I8)') 'ni1_hre  :', ni1_hre
      write(0,'(A,7I8)') 'ni2_hre  :', ni2_hre
      write(0,*) ''
      write(0,'(A,7I8)') 'n_hri    :', n_hri
      write(0,'(A,7I8)') 'ni1_hri  :', ni1_hri
      write(0,'(A,7I8)') 'ni2_hri  :', ni2_hri
    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine print_parallelisation_info

end module graph_parallelisation