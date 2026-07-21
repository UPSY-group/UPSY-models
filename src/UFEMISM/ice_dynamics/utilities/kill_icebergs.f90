! module remove_icebergs
!   use precisions, only: dp
!   use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
!   use mesh_types, only: type_mesh
!   use ice_geometry_calculations, only: type_ice_geometry
!   use mpi_distributed_memory, only: gather_to_all
!   use mpi_basic, only: par
!   use mpi_distributed_memory, only: distribute_from_primary

!   implicit none

!   private

!   public :: kill_icebergs

! contains

!   subroutine kill_icebergs(mesh, icegeom, Hi, dHi_dt)
!     !< Icebergs slow down the model so they are removed using the flood fill method.

!     ! In- and output variables:
!     type(type_ice_geometry),                intent(in   ) :: icegeom
!     type(type_mesh),                        intent(in   ) :: mesh
!     real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: Hi
!     real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: dHi_dt

!     ! Local variables:
!     character(len=1024), parameter         :: routine_name = 'kill_icebergs'
!     real(dp), dimension(mesh%vi1:mesh%vi2) :: iceberg_mask, iceberg_mask_fin
!     real(dp), dimension(mesh%nV)           :: iceberg_mask_tot
!     integer, dimension(mesh%nV)            :: map
!     integer, dimension(mesh%nV)            :: stack
!     integer                                :: stackN
!     integer                                :: vi, ci, vj
!     integer                                :: ierr


!     ! Add routine to path
!     call init_routine(routine_name)

!     map    = 0
!     stack  = 0
!     stackN = 0

!     ! Initialise mask
!     iceberg_mask = 0.0_dp

!     ! Associate each element (open ocean, floating ice and grounded ice) to a value for the mask
!     do vi = mesh%vi1, mesh%vi2
!         if (icegeom%mask_grounded_ice(vi) .eqv. .TRUE.) then
!         iceberg_mask(vi) = 1.0 ! grounded
!         else if (icegeom%mask_floating_ice(vi) .eqv. .TRUE.) then
!         iceberg_mask(vi) = 2.0 ! floating
!         else
!         iceberg_mask(vi) = 0.0 ! ocean
!         end if
!     end do

!     ! Gather all so the flood fill method has access to all neighbours
!     call gather_to_all( iceberg_mask, iceberg_mask_tot)

!     ! Run on primary core
!     if (par%primary) then

!         ! Fill the stack with all shelf-next-to-sheet points
!         do vi = 1, mesh%nV
!           if (iceberg_mask_tot(vi) == 1.0_dp) then                           ! if grounded
!             do ci = 1, mesh%nC(vi)                                           ! check neighbours
!                vj = mesh%C(vi, ci)
!                 if (iceberg_mask_tot(vj) == 2.0_dp .and. map(vj) == 0) then  ! if neighbour is floating and not visited yet
!                   stackN = stackN + 1                                        ! add to list of floating points that are connected to the grounded part
!                   stack(stackN) = vj
!                   map(vj) = 1_dp
!                 end if
!             end do
!           end if
!         end do

!         ! Start finding the neighbours that are floating using the flood fill method
!         do while (stackN > 0)
!             vi = stack(stackN)
!             stackN = stackN - 1
!             ! Check all neighbors
!             do ci = 1, mesh%nC(vi)
!                 vj = mesh%C(vi, ci)
!                 if ( map(vj) .ne. 1_dp) then                                  ! hasnt been visited yet
!                     if (iceberg_mask_tot(vj) == 2.0_dp) then                  ! if it's floating add to stack
!                         map(vj) = 1_dp
!                         stackN = stackN + 1
!                         stack(stackN) = vj
!                     end if
!                 end if
!             end do
!         end do

!         ! Mark the other floating points as icebergs
!         do vi = 1, mesh%nV
!             if (iceberg_mask_tot(vi) == 2.0_dp .and. map(vi) == 0) then        ! is floating but not marked as connected to grounded ice
!                 iceberg_mask_tot(vi) = 3.0_dp
!             end if
!         end do
!     end if

!     call distribute_from_primary(iceberg_mask_fin, iceberg_mask_tot)

!     ! Remove ice corresponding to icebergs
!     do vi = mesh%vi1, mesh%vi2
!         if (iceberg_mask_fin(vi) == 3.0_dp) then
!             Hi(vi) = 0.0_dp
!             dHi_dt(vi) = 0.0_dp
!         end if
!     end do

!     call finalise_routine(routine_name)

!   end subroutine kill_icebergs

! end module remove_icebergs
