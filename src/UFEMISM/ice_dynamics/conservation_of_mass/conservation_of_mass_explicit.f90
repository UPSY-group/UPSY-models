module conservation_of_mass_explicit

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use ice_geometry_basics, only: ice_surface_elevation, Hi_from_Hb_Hs_and_SL
  use mpi_distributed_memory, only: gather_to_all
  use conservation_of_mass_utilities, only: calc_ice_flux_divergence_matrix_upwind, &
    calc_flux_limited_timestep, apply_mask_noice_direct, calc_n_interior_neighbours
  use CSR_matrix_vector_multiplication, only: multiply_csr_matrix_with_vector_1d_wrapper
  use checksum_mod, only: checksum
  use ice_thickness_boundary_conditions, only: apply_ice_thickness_BC_explicit

  implicit none

  private

  public :: calc_dHi_dt_explicit

contains

  subroutine calc_dHi_dt_explicit( mesh, geom, u_vav_perp, SMB, BMB, LMB, AMB, &
    mask_noice, dt, dHi_dt, Hi_tplusdt, divQ, dHi_dt_target, BC_prescr_mask, BC_prescr_Hi)
    !< Calculate ice thickness rates of change (dH/dt)
    !< Use a time-explicit discretisation scheme for the ice fluxes

    ! The ice continuity equation (alternatively known as the ice thickness equation,
    ! or just conservation of mass) reads:
    !
    !     [1] dH/dt = -div( Q) + m
    !
    ! Here, Q is the horizontal ice flux vector, and m is the net mass balance.
    !
    ! We define a matrix operator M_divQ that can be multiplied with the ice thickness
    ! to produce the flux divergence:
    !
    !     [2] div( Q) = M_divQ * H
    !
    ! Substituting [2] into [1] yields:
    !
    !     [3] dH/dt = -M_divQ H + m
    !
    ! Using a time-explicit discretisation scheme so that H on the right-hand side
    ! is defined at time t yields:
    !
    !     [4] (H( t+dt) - H( t)) / dt = -M_divQ H( t) + m

    ! In/output variables:
    type(type_mesh),                                     intent(in   )           :: mesh                  ! [-]       The model mesh
    class(atype_ice_geometry_model_data),                intent(in   )           :: geom                  !           The ice-sheet geometry
    real(dp), dimension(mesh%vi1:mesh%vi2, mesh%nC_mem), intent(in   )           :: u_vav_perp                ! [m yr^-1] Vertically-averaged ice velocity components perpendicular to Voronoi cell boundaries
    real(dp), dimension(mesh%vi1:mesh%vi2),              intent(in   )           :: SMB                   ! [m yr^-1] Surface mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2),              intent(in   )           :: BMB                   ! [m yr^-1] Basal   mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2),              intent(in   )           :: LMB                   ! [m yr^-1] Lateral mass balance
    real(dp), dimension(mesh%vi1:mesh%vi2),              intent(inout)           :: AMB                   ! [m yr^-1] Artificial mass balance
    logical,  dimension(mesh%vi1:mesh%vi2),              intent(in   )           :: mask_noice            ! [-]       Mask of vertices where no ice is allowed
    real(dp),                                            intent(inout)           :: dt                    ! [dt]      Time step
    real(dp), dimension(mesh%vi1:mesh%vi2),              intent(  out)           :: dHi_dt                ! [m yr^-1] Ice thickness rate of change
    real(dp), dimension(mesh%vi1:mesh%vi2),              intent(  out)           :: Hi_tplusdt            ! [m]       Ice thickness at time t + dt
    real(dp), dimension(mesh%vi1:mesh%vi2),              intent(  out)           :: divQ                  ! [m yr^-1] Horizontal ice flux divergence
    real(dp), dimension(mesh%vi1:mesh%vi2),              intent(in   )           :: dHi_dt_target         ! [m yr^-1] Target ice thickness rate of change
    integer,  dimension(mesh%vi1:mesh%vi2),              intent(in   ), optional :: BC_prescr_mask        ! [-]       Mask of vertices where thickness is prescribed
    real(dp), dimension(mesh%vi1:mesh%vi2),              intent(in   ), optional :: BC_prescr_Hi          ! [m]       Prescribed thicknesses

    ! Local variables:
    character(len=*), parameter  :: routine_name = 'calc_dHi_dt_explicit'
    type(type_CSR_matrix_dp)     :: M_divQ
    real(dp)                     :: dt_max
    integer                      :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the ice flux divergence matrix M_divQ using an upwind scheme
    call calc_ice_flux_divergence_matrix_upwind( mesh, u_vav_perp, geom%fraction_margin, M_divQ)

    ! Calculate the ice flux divergence div(Q)
    call multiply_CSR_matrix_with_vector_1D_wrapper( M_divQ, &
      mesh%pai_V, geom%Hi, mesh%pai_V, divQ, &
      buffer_xx_nih = mesh%buffer1_d_a_nih, buffer_yy_nih = mesh%buffer2_d_a_nih)

    ! Calculate rate of ice thickness change dHi/dt
    do vi = mesh%vi1, mesh%vi2
      dHi_dt( vi) = -divQ( vi) + geom%fraction_margin( vi) * (SMB( vi) + BMB( vi) - dHi_dt_target( vi)) + LMB( vi)
    end do

    ! Store this value in the artificial mass balance field
    AMB( mesh%vi1: mesh%vi2) = dHi_dt( mesh%vi1: mesh%vi2)

    ! Calculate largest time step possible based on dHi_dt
    call calc_flux_limited_timestep( mesh, geom, dHi_dt, dt_max)

    ! Constrain dt based on new limit
    dt = min( dt, dt_max)

    ! Calculate ice thickness at t+dt
    Hi_tplusdt( mesh%vi1: mesh%vi2) = max( 0._dp, geom%Hi( mesh%vi1: mesh%vi2) + dHi_dt( mesh%vi1: mesh%vi2) * dt)

    ! Apply boundary conditions at the domain border
    call apply_ice_thickness_BC_explicit( mesh, mask_noice, geom%Hb, geom%SL, Hi_tplusdt)

    ! Set predicted ice thickness to prescribed values where told to do so
    if (present( BC_prescr_mask) .or. present( BC_prescr_Hi)) then
      ! Safety
      if (.not. (present( BC_prescr_mask) .and. present( BC_prescr_Hi))) then
        call crash('need to provide prescribed both Hi and mask!')
      end if
      do vi = mesh%vi1, mesh%vi2
        if (BC_prescr_mask( vi) == 1) then
          Hi_tplusdt( vi) = max( 0._dp, BC_prescr_Hi( vi))
        end if
      end do
    end if

    ! Enforce Hi = 0 where told to do so
    call apply_mask_noice_direct( mesh, mask_noice, Hi_tplusdt)

    ! Recalculate dH/dt, accounting for limit of no negative ice thickness
    dHi_dt( mesh%vi1: mesh%vi2) = (Hi_tplusdt( mesh%vi1: mesh%vi2) - geom%Hi( mesh%vi1: mesh%vi2)) / dt

    ! Remove the final dH/dt field, which now includes some
    ! artificial ice modifications, from the original field
    ! stored in the AMB field. Any residuals will represent
    ! the component of the original dH/dt that was removed.
    ! The negative of this we call artificial mass balance.
    AMB( mesh%vi1: mesh%vi2) = dHi_dt( mesh%vi1: mesh%vi2) - AMB( mesh%vi1: mesh%vi2)

    call checksum( mesh%pai_V, AMB       , 'AMB')
    call checksum( mesh%pai_V, dHi_dt    , 'dHi_dt')
    call checksum( mesh%pai_V, Hi_tplusdt, 'Hi_tplusdt')
    call checksum( mesh%pai_V, divQ      , 'divQ')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_dHi_dt_explicit

end module conservation_of_mass_explicit
