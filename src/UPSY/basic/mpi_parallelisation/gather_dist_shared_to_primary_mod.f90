module gather_dist_shared_to_primary_mod

  ! Gather a hybrid distributed/shared array including halos to the primary

  use assertions_basic
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mpi_f08, only: MPI_INTEGER, MPI_ALLGATHER, MPI_GATHERV, MPI_COMM_WORLD, &
    MPI_DOUBLE_PRECISION, MPI_LOGICAL, MPI_DOUBLE_COMPLEX
  use parallel_array_info_type, only: type_par_arr_info

  implicit none

  private

  public :: gather_dist_shared_to_primary

  interface gather_dist_shared_to_primary
    procedure :: gather_dist_shared_to_primary_logical_1D
    procedure :: gather_dist_shared_to_primary_logical_2D
    procedure :: gather_dist_shared_to_primary_logical_3D
    procedure :: gather_dist_shared_to_primary_int_1D
    procedure :: gather_dist_shared_to_primary_int_2D
    procedure :: gather_dist_shared_to_primary_int_3D
    procedure :: gather_dist_shared_to_primary_dp_1D
    procedure :: gather_dist_shared_to_primary_dp_2D
    procedure :: gather_dist_shared_to_primary_dp_3D
    procedure :: gather_dist_shared_to_primary_complex_1D
    procedure :: gather_dist_shared_to_primary_complex_2D
    procedure :: gather_dist_shared_to_primary_complex_3D
  end interface gather_dist_shared_to_primary

contains

  subroutine gather_dist_shared_to_primary_logical_1D( pai, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                           intent(in   ) :: pai
    logical, dimension(pai%i1_nih:pai%i2_nih), target, intent(in   ) :: d_nih
    logical, dimension(1:pai%n), optional, target,     intent(  out) :: d_tot

    ! Local variables:
    character(len=*), parameter    :: routine_name = 'gather_dist_shared_to_primary_logical_1D'
    logical, dimension(:), pointer :: d_loc
    integer                        :: ierr, i
    integer, dimension(1:par%n)    :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    call assert( ((par%primary .and. present( d_tot)) .or. &
      (.not. par%primary .and. .not. present( d_tot))), 'd_tot should only be present on primary')

    ! Make sure all other processes have finished writing
    call sync

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_nih
      call finalise_routine( routine_name)
      return
    end if

    ! Determine ranges owned by each process
    call MPI_ALLGATHER( pai%n_loc, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    if( sum( counts) /= pai%n) call crash('inconsistency between pai%n_loc and pai%n')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    do i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    end do

    ! Gather data to the primary
    d_loc( pai%i1:pai%i2) => d_nih( pai%i1:pai%i2)
    call MPI_GATHERV( d_loc, pai%n_loc, MPI_LOGICAL, &
      d_tot, counts, displs, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_logical_1D

  subroutine gather_dist_shared_to_primary_logical_2D( pai, nz, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                intent(in   ) :: pai
    integer,                                                intent(in   ) :: nz
    logical, dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(in   ) :: d_nih
    logical, dimension(1:pai%n,1:nz), optional, target,     intent(  out) :: d_tot

    ! Local variables:
    character(len=*), parameter    :: routine_name = 'gather_dist_shared_to_primary_logical_2D'
    logical, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                        :: k

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        d_tot_1D => d_tot(:,k)
        call gather_dist_shared_to_primary_logical_1D( pai, d_nih_1D, d_tot_1D)
      end do

    else

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        call gather_dist_shared_to_primary_logical_1D( pai, d_nih_1D)
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_logical_2D

  subroutine gather_dist_shared_to_primary_logical_3D( pai, nz, nl, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                     intent(in   ) :: pai
    integer,                                                     intent(in   ) :: nz, nl
    logical, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(in   ) :: d_nih
    logical, dimension(1:pai%n,1:nz,1:nl), optional, target,     intent(  out) :: d_tot

    ! Local variables:
    character(len=*), parameter    :: routine_name = 'gather_dist_shared_to_primary_logical_3D'
    logical, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                        :: k,l

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          d_tot_1D => d_tot(:,k,l)
          call gather_dist_shared_to_primary_logical_1D( pai, d_nih_1D, d_tot_1D)
        end do
      end do

    else

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          call gather_dist_shared_to_primary_logical_1D( pai, d_nih_1D)
        end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_logical_3D

  subroutine gather_dist_shared_to_primary_int_1D( pai, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                           intent(in   ) :: pai
    integer, dimension(pai%i1_nih:pai%i2_nih), target, intent(in   ) :: d_nih
    integer, dimension(1:pai%n), optional, target,     intent(  out) :: d_tot

    ! Local variables:
    character(len=*), parameter    :: routine_name = 'gather_dist_shared_to_primary_int_1D'
    integer, dimension(:), pointer :: d_loc
    integer                        :: ierr, i
    integer, dimension(1:par%n)    :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    call assert( ((par%primary .and. present( d_tot)) .or. &
      (.not. par%primary .and. .not. present( d_tot))), 'd_tot should only be present on primary')

    ! Make sure all other processes have finished writing
    call sync

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_nih
      call finalise_routine( routine_name)
      return
    end if

    ! Determine ranges owned by each process
    call MPI_ALLGATHER( pai%n_loc, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    if( sum( counts) /= pai%n) call crash('inconsistency between pai%n_loc and pai%n')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    do i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    end do

    ! Gather data to the primary
    d_loc( pai%i1:pai%i2) => d_nih( pai%i1:pai%i2)
    call MPI_GATHERV( d_loc, pai%n_loc, MPI_INTEGER, &
      d_tot, counts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_int_1D

  subroutine gather_dist_shared_to_primary_int_2D( pai, nz, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                intent(in   ) :: pai
    integer,                                                intent(in   ) :: nz
    integer, dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(in   ) :: d_nih
    integer, dimension(1:pai%n,1:nz), optional, target,     intent(  out) :: d_tot

    ! Local variables:
    character(len=*), parameter    :: routine_name = 'gather_dist_shared_to_primary_int_2D'
    integer, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                        :: k

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        d_tot_1D => d_tot(:,k)
        call gather_dist_shared_to_primary_int_1D( pai, d_nih_1D, d_tot_1D)
      end do

    else

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        call gather_dist_shared_to_primary_int_1D( pai, d_nih_1D)
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_int_2D

  subroutine gather_dist_shared_to_primary_int_3D( pai, nz, nl, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                     intent(in   ) :: pai
    integer,                                                     intent(in   ) :: nz, nl
    integer, dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(in   ) :: d_nih
    integer, dimension(1:pai%n,1:nz,1:nl), optional, target,     intent(  out) :: d_tot

    ! Local variables:
    character(len=*), parameter    :: routine_name = 'gather_dist_shared_to_primary_int_3D'
    integer, dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                        :: k,l

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          d_tot_1D => d_tot(:,k,l)
          call gather_dist_shared_to_primary_int_1D( pai, d_nih_1D, d_tot_1D)
        end do
      end do

    else

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          call gather_dist_shared_to_primary_int_1D( pai, d_nih_1D)
        end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_int_3D

  subroutine gather_dist_shared_to_primary_dp_1D( pai, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                            intent(in   ) :: pai
    real(dp), dimension(pai%i1_nih:pai%i2_nih), target, intent(in   ) :: d_nih
    real(dp), dimension(1:pai%n), optional, target,     intent(  out) :: d_tot

    ! Local variables:
    character(len=*), parameter     :: routine_name = 'gather_dist_shared_to_primary_dp_1D'
    real(dp), dimension(:), pointer :: d_loc
    integer                         :: ierr, i
    integer, dimension(1:par%n)     :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    call assert( ((par%primary .and. present( d_tot)) .or. &
      (.not. par%primary .and. .not. present( d_tot))), 'd_tot should only be present on primary')

    ! Make sure all other processes have finished writing
    call sync

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_nih
      call finalise_routine( routine_name)
      return
    end if

    ! Determine ranges owned by each process
    call MPI_ALLGATHER( pai%n_loc, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    if( sum( counts) /= pai%n) call crash('inconsistency between pai%n_loc and pai%n')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    do i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    end do

    ! Gather data to the primary
    d_loc( pai%i1:pai%i2) => d_nih( pai%i1:pai%i2)
    call MPI_GATHERV( d_loc, pai%n_loc, MPI_DOUBLE_PRECISION, &
      d_tot, counts, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_dp_1D

  subroutine gather_dist_shared_to_primary_dp_2D( pai, nz, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                 intent(in   ) :: pai
    integer,                                                 intent(in   ) :: nz
    real(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(in   ) :: d_nih
    real(dp), dimension(1:pai%n,1:nz), optional, target,     intent(  out) :: d_tot

    ! Local variables:
    character(len=*), parameter     :: routine_name = 'gather_dist_shared_to_primary_dp_2D'
    real(dp), dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                         :: k

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        d_tot_1D => d_tot(:,k)
        call gather_dist_shared_to_primary_dp_1D( pai, d_nih_1D, d_tot_1D)
      end do

    else

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        call gather_dist_shared_to_primary_dp_1D( pai, d_nih_1D)
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_dp_2D

  subroutine gather_dist_shared_to_primary_dp_3D( pai, nz, nl, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                      intent(in   ) :: pai
    integer,                                                      intent(in   ) :: nz, nl
    real(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(in   ) :: d_nih
    real(dp), dimension(1:pai%n,1:nz,1:nl), optional, target,     intent(  out) :: d_tot

    ! Local variables:
    character(len=*), parameter     :: routine_name = 'gather_dist_shared_to_primary_dp_3D'
    real(dp), dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                         :: k,l

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          d_tot_1D => d_tot(:,k,l)
          call gather_dist_shared_to_primary_dp_1D( pai, d_nih_1D, d_tot_1D)
        end do
      end do

    else

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          call gather_dist_shared_to_primary_dp_1D( pai, d_nih_1D)
        end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_dp_3D

  subroutine gather_dist_shared_to_primary_complex_1D( pai, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                               intent(in   ) :: pai
    complex(dp), dimension(pai%i1_nih:pai%i2_nih), target, intent(in   ) :: d_nih
    complex(dp), dimension(1:pai%n), optional, target,     intent(  out) :: d_tot

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'gather_dist_shared_to_primary_complex_1D'
    complex(dp), dimension(:), pointer :: d_loc
    integer                            :: ierr, i
    integer, dimension(1:par%n)        :: counts, displs

    ! Add routine to path
    call init_routine( routine_name)

    call assert( ((par%primary .and. present( d_tot)) .or. &
      (.not. par%primary .and. .not. present( d_tot))), 'd_tot should only be present on primary')

    ! Make sure all other processes have finished writing
    call sync

    ! Exception when we're running on a single node
    if (par%n_nodes == 1) then
      if (par%primary) d_tot = d_nih
      call finalise_routine( routine_name)
      return
    end if

    ! Determine ranges owned by each process
    call MPI_ALLGATHER( pai%n_loc, 1, MPI_INTEGER, counts, 1, MPI_INTEGER, MPI_COMM_WORLD, ierr)
    if( sum( counts) /= pai%n) call crash('inconsistency between pai%n_loc and pai%n')

    ! Calculate displacements for MPI_GATHERV
    displs( 1) = 0
    do i = 2, par%n
      displs( i) = displs( i-1) + counts( i-1)
    end do

    ! Gather data to the primary
    d_loc( pai%i1:pai%i2) => d_nih( pai%i1:pai%i2)
    call MPI_GATHERV( d_loc, pai%n_loc, MPI_DOUBLE_COMPLEX, &
      d_tot, counts, displs, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD, ierr)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_complex_1D

  subroutine gather_dist_shared_to_primary_complex_2D( pai, nz, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                    intent(in   ) :: pai
    integer,                                                    intent(in   ) :: nz
    complex(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz), target, intent(in   ) :: d_nih
    complex(dp), dimension(1:pai%n,1:nz), optional, target,     intent(  out) :: d_tot

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'gather_dist_shared_to_primary_complex_2D'
    complex(dp), dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                            :: k

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        d_tot_1D => d_tot(:,k)
        call gather_dist_shared_to_primary_complex_1D( pai, d_nih_1D, d_tot_1D)
      end do

    else

      do k = 1, nz
        d_nih_1D => d_nih(:,k)
        call gather_dist_shared_to_primary_complex_1D( pai, d_nih_1D)
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_complex_2D

  subroutine gather_dist_shared_to_primary_complex_3D( pai, nz, nl, d_nih, d_tot)

    ! In/output variables:
    type(type_par_arr_info),                                         intent(in   ) :: pai
    integer,                                                         intent(in   ) :: nz, nl
    complex(dp), dimension(pai%i1_nih:pai%i2_nih,1:nz,1:nl), target, intent(in   ) :: d_nih
    complex(dp), dimension(1:pai%n,1:nz,1:nl), optional, target,     intent(  out) :: d_tot

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'gather_dist_shared_to_primary_complex_3D'
    complex(dp), dimension(:), pointer :: d_nih_1D, d_tot_1D
    integer                            :: k,l

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          d_tot_1D => d_tot(:,k,l)
          call gather_dist_shared_to_primary_complex_1D( pai, d_nih_1D, d_tot_1D)
        end do
      end do

    else

      do k = 1, nz
        do l = 1, nl
          d_nih_1D => d_nih(:,k,l)
          call gather_dist_shared_to_primary_complex_1D( pai, d_nih_1D)
        end do
      end do

    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine gather_dist_shared_to_primary_complex_3D

end module gather_dist_shared_to_primary_mod
