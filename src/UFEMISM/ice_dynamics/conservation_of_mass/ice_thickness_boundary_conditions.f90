module ice_thickness_boundary_conditions

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_geometry_basics, only: ice_surface_elevation, Hi_from_Hb_Hs_and_SL
  use mpi_distributed_memory, only: gather_to_all
  use conservation_of_mass_utilities, only: calc_n_interior_neighbours

  implicit none

  private

  public :: apply_ice_thickness_BC_explicit

contains

  subroutine apply_ice_thickness_BC_explicit( mesh, mask_noice, Hb, SL, Hi)
    !< Apply boundary conditions to the ice thickness on the domain border directly

    ! In/output variables:
    type(type_mesh),                                          intent(in   ) :: mesh
    logical,  dimension(mesh%vi1:mesh%vi2),                   intent(in   ) :: mask_noice
    real(dp), dimension(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), intent(in   ) :: Hb
    real(dp), dimension(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), intent(in   ) :: SL
    real(dp), dimension(:), target,                           intent(inout) :: Hi

    ! Local variables:
    character(len=*), parameter            :: routine_name = 'apply_ice_thickness_BC_explicit'
    integer,  dimension(mesh%vi1:mesh%vi2) :: n_interior_neighbours
    real(dp), dimension(mesh%vi1:mesh%vi2) :: Hs_tplusdt
    integer                                :: vi
    real(dp), dimension(mesh%nV)           :: Hs_tplusdt_tot
    logical,  dimension(mesh%nV)           :: mask_noice_tot
    character(len=:), allocatable          :: BC_H
    integer                                :: ci,vj
    real(dp)                               :: Hs_sum, Hs_av
    real(dp), dimension(:), pointer        :: Hi_nih, Hi_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Support both distributed and hybrid distributed/shared versions of Hi
    if (size( Hi,1) == mesh%pai_V%n_loc) then
      Hi_loc( mesh%vi1:mesh%vi2) => Hi
    elseif (size( Hi,1) == mesh%pai_V%n_nih) then
      Hi_nih( mesh%pai_V%i1_nih:mesh%pai_V%i2_nih) => Hi
      Hi_loc( mesh%vi1:mesh%vi2) => Hi_nih( mesh%vi1:mesh%vi2)
    else
      call crash('Hi has invalid size')
    end if

    ! Calculate Hs( t+dt)
    do vi = mesh%vi1, mesh%vi2
      Hs_tplusdt( vi) = ice_surface_elevation( Hi_loc( vi), Hb( vi), SL( vi))
    end do

    ! Gather global data fields
    call gather_to_all( Hs_tplusdt, Hs_tplusdt_tot)
    call gather_to_all( mask_noice, mask_noice_tot)

    call calc_n_interior_neighbours( mesh, mask_noice, n_interior_neighbours)

    ! == First pass: set values of border vertices to mean of interior neighbours
    !    ...for those border vertices that actually have interior neighbours.
    ! ===========================================================================

    do vi = mesh%vi1, mesh%vi2

      if     (mesh%VBI( vi) == 1 .or. mesh%VBI( vi) == 2) then
        ! Northern domain border
        BC_H = C%BC_H_north
      elseif (mesh%VBI( vi) == 3 .or. mesh%VBI( vi) == 4) then
        ! Eastern domain border
        BC_H = C%BC_H_east
      elseif (mesh%VBI( vi) == 5 .or. mesh%VBI( vi) == 6) then
        ! Southern domain border
        BC_H = C%BC_H_south
      elseif (mesh%VBI( vi) == 7 .or. mesh%VBI( vi) == 8) then
        ! Western domain border
        BC_H = C%BC_H_west
      else
        ! Free vertex
        cycle
      end if

      select case (BC_H)
      case default
        call crash('unknown BC_H "' // trim( BC_H) // '"')
      case ('zero')
        ! Set ice thickness to zero here

        Hi_loc( vi) = 0._dp

      case ('infinite')
        ! Set H on this vertex equal to the average value on its neighbours

        if (n_interior_neighbours( vi) > 0) then

          Hs_sum = 0._dp
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            if (mesh%VBI( vj) == 0 .and. .not. mask_noice_tot( vj)) then
              Hs_sum = Hs_sum + Hs_tplusdt_tot( vj)
            end if
          end do
          Hs_av = Hs_sum / real( n_interior_neighbours( vi),dp)

          Hs_tplusdt( vi) = max( Hb( vi), Hs_av)
          Hi_loc( vi) = Hi_from_Hb_Hs_and_SL( Hb( vi), Hs_tplusdt( vi), SL( vi))

        end if

      end select

    end do

    ! == Second pass: set values of border vertices to mean of all neighbours
    !    ...for those border vertices that have no interior neighbours.
    ! =======================================================================

    ! Gather global data fields again
    call gather_to_all( Hs_tplusdt, Hs_tplusdt_tot)

    do vi = mesh%vi1, mesh%vi2

      if     (mesh%VBI( vi) == 1 .or. mesh%VBI( vi) == 2) then
        ! Northern domain border
        BC_H = C%BC_H_north
      elseif (mesh%VBI( vi) == 3 .or. mesh%VBI( vi) == 4) then
        ! Eastern domain border
        BC_H = C%BC_H_east
      elseif (mesh%VBI( vi) == 5 .or. mesh%VBI( vi) == 6) then
        ! Southern domain border
        BC_H = C%BC_H_south
      elseif (mesh%VBI( vi) == 7 .or. mesh%VBI( vi) == 8) then
        ! Western domain border
        BC_H = C%BC_H_west
      else
        ! Free vertex
        cycle
      end if

      select case (BC_H)
      case default
        call crash('unknown BC_H "' // trim( BC_H) // '"')
      case ('zero')
        ! Set ice thickness to zero here

        Hi_loc( vi) = 0._dp

      case ('infinite')
        ! Set H on this vertex equal to the average value on its neighbours

        if (n_interior_neighbours( vi) == 0) then

          Hs_sum = 0._dp
          do ci = 1, mesh%nC( vi)
            vj = mesh%C( vi,ci)
            Hs_sum = Hs_sum + Hs_tplusdt_tot( vj)
          end do
          Hs_av = Hs_sum / real( mesh%nC( vi),dp)

          Hs_tplusdt( vi) = max( Hb( vi), Hs_av)
          Hi_loc( vi) = Hi_from_Hb_Hs_and_SL( Hb( vi), Hs_tplusdt( vi), SL( vi))

        end if

      end select

    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_ice_thickness_BC_explicit

end module ice_thickness_boundary_conditions
