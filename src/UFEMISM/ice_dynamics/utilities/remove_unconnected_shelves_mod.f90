module remove_unconnected_shelves_mod

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use mpi_distributed_memory, only: gather_to_primary
  use mpi_basic, only: par, sync
  use mpi_distributed_memory, only: distribute_from_primary
  use ice_geometry_basics, only: is_floating
  use netcdf_io_main
  use model_configuration, only: C

  implicit none

  private

  public :: remove_unconnected_shelves

contains

  subroutine remove_unconnected_shelves( mesh, Hb, SL, Hi)
    !< Detect shelves that are not connected to grounded ice,
    !< and remove them (i.e., set ice thickness there to zero)

    ! In- and output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), intent(in   ) :: Hb
    real(dp), dimension(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), intent(in   ) :: SL
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: Hi

    ! Local variables:
    character(len=*), parameter            :: routine_name = 'remove_unconnected_shelves'
    logical,  dimension(:), allocatable    :: mask_grounded_ice_tot
    logical,  dimension(:), allocatable    :: mask_floating_ice_tot
    integer                                :: vi, ci, vj
    integer,  dimension(:), allocatable    :: map, stack
    integer                                :: stackN
    logical                                :: is_next_to_grounded_ice
    logical,  dimension(:), allocatable    :: mask_ice_connected_to_grounded_ice_tot
    logical,  dimension(mesh%vi1:mesh%vi2) :: mask_ice_connected_to_grounded_ice
    real(dp) :: x0, y0, r

    ! Add routine to path
    call init_routine(routine_name)

    call calc_ice_masks( mesh, Hi, Hb, SL, mask_grounded_ice_tot, mask_floating_ice_tot)

    ! Flood-fill stuff cannot be parallellised; let the primary do the work
    if (par%primary) then

      allocate( map  ( mesh%nV), source = 0)
      allocate( stack( mesh%nV), source = 0)
      stackN = 0

      ! Start by marking all grounded ice on the map
      do vi = 1, mesh%nV
        if (mask_grounded_ice_tot( vi)) map( vi) = 2
      end do

      ! Initialise the stack with all shelf-next-to-sheet vertices
      do vi = 1, mesh%nV
        if (mask_floating_ice_tot( vi)) then
          is_next_to_grounded_ice = .false.
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            if (mask_grounded_ice_tot( vj)) then
              is_next_to_grounded_ice = .true.
              exit
            end if
          end do
          if (is_next_to_grounded_ice) then
            map( vi) = 1
            stackN = stackN + 1
            stack( stackN) = vi
          end if
        end if
      end do

      ! Expand the stack outward, marking only floating-ice vertices
      do while (stackN > 0)

        ! Take the last element from the stack
        vi = stack( stackN)
        stackN = stackN - 1

        ! If it is floating, mark it on the map, and add its non-mapped, non-stacked neighbours to the stack
        if (mask_floating_ice_tot( vi)) then
          map( vi) = 2
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            if (map( vj) == 0) then
              map( vj) = 1
              stackN = stackN + 1
              stack( stackN) = vj
            end if
          end do
        end if

      end do

      ! Fill in the ice-connected-to-grounded-ice-mask
      allocate( mask_ice_connected_to_grounded_ice_tot( mesh%nV), source = .false.)
      do vi = 1, mesh%nV
        if (map( vi) == 2) mask_ice_connected_to_grounded_ice_tot( vi) = .true.
      end do

    end if
    call sync

    call distribute_from_primary( mask_ice_connected_to_grounded_ice, d_tot=mask_ice_connected_to_grounded_ice_tot)

    ! Remove ice that is not connected to grounded ice
    do vi = mesh%vi1, mesh%vi2
      if (.not. mask_ice_connected_to_grounded_ice( vi)) Hi( vi) = 0._dp
    end do

    ! Finalise routine path
    call finalise_routine(routine_name)

  end subroutine remove_unconnected_shelves

  subroutine calc_ice_masks( mesh, Hi, Hb, SL, mask_grounded_ice_tot, mask_floating_ice_tot)

    ! In- and output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), intent(in   ) :: Hb
    real(dp), dimension(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), intent(in   ) :: SL
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hi
    logical,  dimension(:), allocatable,    intent(inout) :: mask_grounded_ice_tot
    logical,  dimension(:), allocatable,    intent(inout) :: mask_floating_ice_tot

    ! Local variables:
    character(len=*), parameter            :: routine_name = 'calc_ice_masks'
    logical,  dimension(mesh%vi1:mesh%vi2) :: mask_grounded_ice
    logical,  dimension(mesh%vi1:mesh%vi2) :: mask_floating_ice
    integer                                :: vi

    ! Add routine to path
    call init_routine(routine_name)

    mask_grounded_ice = .false.
    mask_floating_ice = .false.

    do vi = mesh%vi1, mesh%vi2
      if (Hi( vi) > 0._dp) then
        if (is_floating( Hi( vi), Hb( vi), SL( vi))) then
          mask_floating_ice( vi) = .true.
        else
          mask_grounded_ice( vi) = .true.
        end if
      end if
    end do

    if (par%primary) then
      allocate( mask_grounded_ice_tot( mesh%nV))
      allocate( mask_floating_ice_tot( mesh%nV))
      call gather_to_primary( mask_grounded_ice, d_tot = mask_grounded_ice_tot)
      call gather_to_primary( mask_floating_ice, d_tot = mask_floating_ice_tot)
    else
      call gather_to_primary( mask_grounded_ice)
      call gather_to_primary( mask_floating_ice)
    end if

    ! Finalise routine path
    call finalise_routine(routine_name)

  end subroutine calc_ice_masks

end module remove_unconnected_shelves_mod
