module calving_threshold_thickness_iterative

  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use ice_geometry_model_basic, only: type_ice_geometry_model
  use mpi_basic, only: par
  use mpi_distributed_memory, only: gather_to_primary, distribute_from_primary

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
    character(len=*), parameter                :: routine_name = 'apply_calving_threshold_thickness_iterative'
    type(type_ice_geometry_model), allocatable :: geom
    logical, dimension(mesh%vi1:mesh%vi2)      :: mask_Hi_below_threshold
    logical, dimension(:), allocatable         :: mask_Hi_below_threshold_tot
    logical, dimension(:), allocatable         :: mask_calve_tot
    integer                                    :: vi
    integer, dimension(:), allocatable         :: map, stack
    integer                                    :: stackN
    integer                                    :: ci, vj
    logical, dimension(mesh%vi1:mesh%vi2)      :: mask_calve

    ! Add routine to call stack
    call init_routine( routine_name)

    allocate( geom)
    call geom%allocate( 'temp_name', mesh)
    geom%Hi( mesh%vi1:mesh%vi2) = Hi
    geom%Hb( mesh%vi1:mesh%vi2) = Hb
    geom%SL( mesh%vi1:mesh%vi2) = SL
    call geom%determine_masks()
    call geom%calc_effective_thickness()

    ! Determine which vertices have an ice thickness below their relevant threshold
    mask_Hi_below_threshold = .false.
    do vi = mesh%vi1, mesh%vi2
      if ((geom%mask_grounded_ice( vi) .and. geom%Hi_eff( vi) < C%calving_threshold_thickness_sheet) .or. &
          (geom%mask_floating_ice( vi) .and. geom%Hi_eff( vi) < C%calving_threshold_thickness_shelf) .or. &
          geom%mask_icefree_land( vi) .or. geom%mask_icefree_ocean( vi)) then
        mask_Hi_below_threshold( vi) = .true.
      end if
    end do
    if (par%primary) then
      allocate( mask_Hi_below_threshold_tot( 1:mesh%nV))
      call gather_to_primary( mask_Hi_below_threshold, d_tot = mask_Hi_below_threshold_tot)
    else
      allocate( mask_Hi_below_threshold_tot( 0))
      call gather_to_primary( mask_Hi_below_threshold)
    end if

    ! Calve inward from the open ocean
    if (par%primary) then

      allocate( mask_calve_tot( mesh%nV), source = .false.)

      ! Initialise flood-fill map and stack
      allocate( map  ( mesh%nV), source = 0)
      allocate( stack( mesh%nV), source = 0)
      stackN = 0

      ! Start at the domain border
      do vi = 1, mesh%nV
        if (mesh%VBI( vi) > 0 .and. mask_Hi_below_threshold_tot( vi)) then
          map( vi) = 1
          stackN = stackN + 1
          stack( stackN) = vi
        end if
      end do

      ! Flood-fill
      do while (stackN > 0)

        ! Take the last element from the stack
        vi = stack( stackN)
        stackN = stackN - 1

        ! Mark it on the map
        map( vi) = 2
        mask_calve_tot( vi) = .true.

        ! Add its non-marked neighbours to the stack
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          if (map( vj) == 0 .and. mask_Hi_below_threshold_tot( vj)) then
            map( vj) = 1
            stackN = stackN + 1
            stack( stackN) = vj
          end if
        end do

      end do

    end if
    if (par%primary) then
      call distribute_from_primary( mask_calve, d_tot = mask_calve_tot)
    else
      call distribute_from_primary( mask_calve)
    end if

    ! Apply the calving mask
    do vi = mesh%vi1, mesh%vi2
      if (mask_calve( vi)) Hi( vi) = 0._dp
    end do

    ! Clean up after yourself
    call geom%deallocate()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine apply_calving_threshold_thickness_iterative

end module calving_threshold_thickness_iterative
