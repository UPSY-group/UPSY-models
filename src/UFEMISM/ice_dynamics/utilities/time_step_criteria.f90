module time_step_criteria

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_MIN
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use ice_velocity_model_data, only: atype_ice_velocity_model_data
  use mpi_distributed_memory, only: gather_to_all
  use mpi_distributed_shared_memory, only: gather_dist_shared_to_all
  use map_velocities_to_c_grid, only: map_velocities_from_b_to_c_2D

  implicit none

  private

  public :: calc_critical_timestep_adv

contains

  subroutine calc_critical_timestep_adv( mesh, geom, vel, dt_crit_adv)
    !< Calculate the critical time step for advective ice flow (CFL criterion)

    ! In- and output variables:
    type(type_mesh),                      intent(in   ) :: mesh
    class(atype_ice_geometry_model_data), intent(in   ) :: geom
    class(atype_ice_velocity_model_data), intent(in   ) :: vel
    real(dp),                             intent(  out) :: dt_crit_adv

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_critical_timestep_adv'
    real(dp), dimension(mesh%nV)           :: Hi_tot
    logical,  dimension(mesh%nV)           :: mask_floating_ice_tot
    real(dp), dimension(mesh%ei1:mesh%ei2) :: u_vav_c, v_vav_c
    real(dp), dimension(mesh%nE)           :: u_vav_c_tot, v_vav_c_tot
    integer                                :: ei, vi, vj
    real(dp)                               :: dist, dt
    real(dp), parameter                    :: dt_correction_factor = 0.9_dp ! Make actual applied time step a little bit smaller, just to be sure.
    integer                                :: ierr

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather global ice thickness
    call gather_dist_shared_to_all( mesh%pai_V, geom%Hi, Hi_tot)
    call gather_dist_shared_to_all( mesh%pai_V, geom%mask_floating_ice, mask_floating_ice_tot)

    ! Calculate vertically averaged ice velocities on the edges
    call map_velocities_from_b_to_c_2D( mesh, vel%u_vav_b, vel%v_vav_b, u_vav_c, v_vav_c)
    call gather_to_all( u_vav_c, u_vav_c_tot)
    call gather_to_all( v_vav_c, v_vav_c_tot)

    ! Initialise time step with maximum allowed value
    dt_crit_adv = C%dt_ice_max

    do ei = mesh%ei1, mesh%ei2

      ! Only check at ice-covered vertices
      vi = mesh%EV( ei,1)
      vj = mesh%EV( ei,2)
      if (Hi_tot( vi) == 0._dp .OR. Hi_tot( vj) == 0._dp) CYCLE

      if (C%do_grounded_only_adv_dt) then
        ! Only check grounded vertices
        if (mask_floating_ice_tot( vi) .OR. mask_floating_ice_tot( vj)) CYCLE
      end if

      dist = norm2( mesh%V( vi,:) - mesh%V( vj,:))
      dt = dist / max( 0.1_dp, abs( u_vav_c_tot( ei)) + abs( v_vav_c_tot( ei))) * dt_correction_factor
      dt_crit_adv = min( dt_crit_adv, dt)

    end do

    call MPI_ALLREDUCE( MPI_IN_PLACE, dt_crit_adv, 1, MPI_doUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
    dt_crit_adv = min( C%dt_ice_max, dt_crit_adv)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_critical_timestep_adv

end module time_step_criteria
