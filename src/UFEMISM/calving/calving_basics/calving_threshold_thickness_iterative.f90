module calving_threshold_thickness_iterative

  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use ice_geometry_basics, only: is_floating
  use mpi_basic, only: par
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD
  use mpi_distributed_memory, only: gather_to_all, gather_to_primary, distribute_from_primary

  implicit none

  private

  public :: apply_calving_threshold_thickness_iterative

contains

  subroutine apply_calving_threshold_thickness_iterative( mesh, Hb, SL, Hi)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hb
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: SL
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: Hi

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'apply_calving_threshold_thickness_iterative'
    logical, dimension(:), allocatable :: mask_ice, mask_open_ocean, mask_front
    integer                            :: nit, nV_calved
    integer                            :: vi
    real(dp)                           :: threshold_thickness
    integer                            :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    allocate( mask_ice       ( mesh%vi1:mesh%vi2))
    allocate( mask_open_ocean( mesh%vi1:mesh%vi2))
    allocate( mask_front     ( mesh%vi1:mesh%vi2))

    nV_calved = 1
    nit = 0
    do while (nV_calved > 0)

      ! Safety
      nit = nit + 1
      if (nit > mesh%nV) call crash('programming error - iterative open-ocean-ice-front-only threshold-thickness calving broke down')

      ! Calculate ice mask and open ocean mask
      call calc_mask_ice                 ( mesh, Hi, mask_ice)
      call calc_mask_open_ocean_floodfill( mesh, mask_ice, mask_open_ocean)
      call calc_mask_ice_front           ( mesh, mask_ice, mask_open_ocean, mask_front)

      nV_calved = 0
      do vi = mesh%vi1, mesh%vi2

        if ((mask_front( vi) .or. mask_open_ocean( vi))) then
          ! Calving is only allowed at the ice front adjacent to the open ocean

          if (is_floating( Hi( vi), Hb( vi), SL( vi))) then
            threshold_thickness = C%calving_threshold_thickness_shelf
          else
            threshold_thickness = C%calving_threshold_thickness_sheet
          end if

          if (Hi( vi) < threshold_thickness) then
            Hi( vi) = 0._dp
            nV_calved = nV_calved + 1
          end if

        end if

      end do

      ! Check if any processes applied any calving
      call MPI_ALLREDUCE( MPI_IN_PLACE, nV_calved, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine apply_calving_threshold_thickness_iterative

  subroutine calc_mask_ice( mesh, Hi, mask_ice)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hi
    logical,  dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_ice

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_mask_ice'
    integer                     :: vi

    ! Add routine to call stack
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      mask_ice( vi) = Hi( vi) > 0._dp
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calc_mask_ice

  subroutine calc_mask_open_ocean_floodfill( mesh, mask_ice, mask_open_ocean)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_ice
    logical, dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_open_ocean

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'calc_mask_open_ocean_floodfill'
    logical, dimension(:), allocatable :: mask_ice_tot, mask_open_ocean_tot
    integer, dimension(:), allocatable :: map, stack
    integer                            :: stackN
    integer                            :: vi, ci, vj

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) then
      allocate( mask_ice_tot       ( mesh%nV), source = .false.)
      allocate( mask_open_ocean_tot( mesh%nV), source = .false.)
      call gather_to_primary( mask_ice, mask_ice_tot)
    else
      allocate( mask_ice_tot       ( 0))
      allocate( mask_open_ocean_tot( 0))
      call gather_to_primary( mask_ice)
    end if

    ! Let the primary do the work
    if (par%primary) then

      allocate( map  ( mesh%nV), source = 0)
      allocate( stack( mesh%nV), source = 0)
      stackN = 0

      ! Initialise the stack with all ice-free border vertices
      do vi = 1, mesh%nV
        if (.not. mask_ice_tot( vi) .and. mesh%VBI( vi) > 0) then
          map( vi) = 1
          stackN = stackN + 1
          stack( stackN) = vi
        end if
      end do

      ! Expand the open ocean inward to the ice front
      do while (stackN > 0)

        ! Take the last vertex from the stack
        vi = stack( stackN)
        stackN = stackN - 1

        ! Mark it on the map
        map( vi) = 2
        mask_open_ocean_tot( vi) = .true.

        ! Add its non-marked neighbours to the stack
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          if (.not. mask_ice_tot( vj) .and. map( vj) == 0) then
            map( vj) = 1
            stackN = stackN + 1
            stack( stackN) = vj
          end if
        end do

      end do

    end if

    ! Distribute the result from the primary
    if (par%primary) then
      call distribute_from_primary( mask_open_ocean, mask_open_ocean_tot)
    else
      call distribute_from_primary( mask_open_ocean)
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calc_mask_open_ocean_floodfill

  subroutine calc_mask_ice_front( mesh, mask_ice, mask_open_ocean, mask_front)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_ice
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_open_ocean
    logical, dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_front

    ! Local variables:
    character(len=*), parameter        :: routine_name = 'calc_mask_ice_front'
    logical, dimension(:), allocatable :: mask_open_ocean_tot
    integer                            :: vi, ci, vj

    ! Add routine to call stack
    call init_routine( routine_name)

    allocate( mask_open_ocean_tot( mesh%nV), source = .false.)
    call gather_to_all( mask_open_ocean, mask_open_ocean_tot)

    do vi = mesh%vi1, mesh%vi2
      mask_front( vi) = .false.
      if (mask_ice( vi)) then
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          if (mask_open_ocean_tot( vj)) then
            mask_front( vi) = .true.
            exit
          end if
        end do
      end if
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine calc_mask_ice_front

end module calving_threshold_thickness_iterative
