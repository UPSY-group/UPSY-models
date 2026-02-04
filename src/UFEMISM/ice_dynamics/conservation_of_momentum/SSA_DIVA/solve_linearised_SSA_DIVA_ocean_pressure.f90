module solve_linearised_SSA_DIVA_ocean_pressure

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use ice_model_types, only: type_ice_velocity_solver_DIVA_graphs
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist, &
    finalise_matrix_CSR_dist
  use mpi_distributed_shared_memory, only: gather_dist_shared_to_all
  use petsc_basic, only: solve_matrix_equation_CSR_PETSc
  use graph_types, only: type_graph_pair

  use netcdf_io_main

  implicit none

  private

  public :: solve_SSA_DIVA_linearised_ocean_pressure

contains

  subroutine solve_SSA_DIVA_linearised_ocean_pressure( DIVA, &
    PETSc_rtol, PETSc_abstol, n_Axb_its)
    !< Solve the linearised SSA

    ! In/output variables:
    type(type_ice_velocity_solver_DIVA_graphs), intent(inout) :: DIVA
    real(dp),                                   intent(in   ) :: PETSc_rtol, PETSc_abstol
    integer,                                    intent(  out) :: n_Axb_its             ! Number of iterations used in the iterative solver

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'solve_SSA_DIVA_linearised_ocean_pressure'
    integer                             :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    type(type_sparse_matrix_CSR_dp)     :: A_CSR
    real(dp), dimension(:), allocatable :: bb_buv
    real(dp), dimension(:), allocatable :: uv_buv
    integer                             :: row_niuv, ni, uv

    real(dp), dimension(:), pointer :: d_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Store the previous solution
    call gather_dist_shared_to_all( DIVA%graphs%graph_b%pai, DIVA%u_vav_b, DIVA%u_b_prev)
    call gather_dist_shared_to_all( DIVA%graphs%graph_b%pai, DIVA%v_vav_b, DIVA%v_b_prev)

    ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
    ! =============================================================================

    ! Matrix size
    ncols           = DIVA%graphs%graph_b%n     * 2      ! from
    ncols_loc       = DIVA%graphs%graph_b%n_loc * 2
    nrows           = DIVA%graphs%graph_b%n     * 2      ! to
    nrows_loc       = DIVA%graphs%graph_b%n_loc * 2
    nnz_est_proc    = DIVA%graphs%M2_ddx_b_b%nnz * 4

    call allocate_matrix_CSR_dist( A_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Allocate memory for the load vector and the solution
    allocate( bb_buv( DIVA%graphs%graph_b%ni1*2-1: DIVA%graphs%graph_b%ni2*2))
    allocate( uv_buv( DIVA%graphs%graph_b%ni1*2-1: DIVA%graphs%graph_b%ni2*2))

    ! Fill in the current velocity solution
    do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2

      ! u
      row_niuv = DIVA%graphs%biuv2n( ni,1)
      uv_buv( row_niuv) = DIVA%u_vav_b( ni)

      ! v
      row_niuv = DIVA%graphs%biuv2n( ni,2)
      uv_buv( row_niuv) = DIVA%v_vav_b( ni)

    end do

    ! == Construct the stiffness matrix for the linearised SSA
    ! ========================================================

    do row_niuv = A_CSR%i1, A_CSR%i2

      ni = DIVA%graphs%n2biuv( row_niuv,1)
      uv = DIVA%graphs%n2biuv( row_niuv,2)

      ! if (BC_prescr_mask_b( ti) == 1) then
      !   ! Dirichlet boundary condition; velocities are prescribed for this triangle

      !   ! Stiffness matrix: diagonal element set to 1
      !   call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

      !   ! Load vector: prescribed velocity
      !   if     (uv == 1) then
      !     bb( row_tiuv) = BC_prescr_u_b( ti)
      !   elseif (uv == 2) then
      !     bb( row_tiuv) = BC_prescr_v_b( ti)
      !   else
      !     call crash('uv can only be 1 or 2!')
      !   end if

      if (DIVA%graphs%graph_b%is_border( ni)) then
        ! Ice margin: apply ocean pressure boundary condition

        call calc_SSA_DIVA_stiffness_matrix_row_BC_ice_front( DIVA%graphs, &
          DIVA%N_b, DIVA%tau_ox_b, DIVA%tau_oy_b, A_CSR, bb_buv, row_niuv)

      else
        ! No boundary conditions apply; solve the SSA

        if (C%do_include_SSADIVA_crossterms) then
          ! Calculate matrix coefficients for the full SSA/DIVA
          call calc_SSA_DIVA_stiffness_matrix_row_free( DIVA%graphs, &
            DIVA%N_b, DIVA%dN_dx_b, DIVA%dN_dy_b, &
            DIVA%beta_eff_b, &
            DIVA%tau_dx_b, DIVA%tau_dy_b, &
            A_CSR, bb_buv, row_niuv)
        else
          ! Calculate matrix coefficients for the SSA sans the gradients of the effective viscosity (the "cross-terms")
          call calc_SSA_DIVA_sans_stiffness_matrix_row_free( DIVA%graphs, &
            DIVA%N_b, &
            DIVA%beta_eff_b, &
            DIVA%tau_dx_b, DIVA%tau_dy_b, &
            A_CSR, bb_buv, row_niuv)
        end if

      end if

    end do

    call finalise_matrix_CSR_dist( A_CSR)

    ! ! == Solve the matrix equation
    ! ! ============================

    ! Use PETSc to solve the matrix equation
    call solve_matrix_equation_CSR_PETSc( A_CSR, bb_buv, uv_buv, PETSc_rtol, PETSc_abstol, &
      n_Axb_its, PETSc_KSPtype = 'bicg', PETSc_PCtype = 'bjacobi')

    ! Disentangle the u and v components of the velocity solution
    do ni = DIVA%graphs%graph_b%ni1, DIVA%graphs%graph_b%ni2

      ! u
      row_niuv = DIVA%graphs%biuv2n( ni,1)
      DIVA%u_vav_b( ni) = uv_buv( row_niuv)

      ! v
      row_niuv = DIVA%graphs%biuv2n( ni,2)
      DIVA%v_vav_b( ni) = uv_buv( row_niuv)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSA_DIVA_linearised_ocean_pressure

  subroutine calc_SSA_DIVA_stiffness_matrix_row_free( graphs, N_b, dN_dx_b, dN_dy_b, &
    basal_friction_coefficient_b, tau_dx_b, tau_dy_b, A_CSR, bb_b, row_niuv)
    !< Add coefficients to this matrix row to represent the linearised SSA

    ! The SSA reads;
    !
    !   d/dx [ 2 N ( 2 du/dx + dv/dy )] + d/dy [ N ( du/dy + dv/dx)] - beta_b u = -tau_dx
    !
    !   d/dy [ 2 N ( 2 dv/dy + du/dx )] + d/dx [ N ( dv/dx + du/dy)] - beta_b v = -tau_dy
    !
    ! Using the chain rule, this expands to read:
    !
    !   4 N d2u/dx2 + 4 dN/dx du/dx + 2 N d2v/dxdy + 2 dN/dx dv/dy + ...
    !     N d2u/dy2 +   dN/dy du/dy +   N d2v/dxdy +   dN/dy dv/dx - beta_b u = -tau_dx
    !
    !   4 N d2v/dy2 + 4 dN/dy dv/dy + 2 N d2u/dxdy + 2 dN/dy du/dx + ...
    !     N d2v/dx2 +   dN/dx dv/dx +   N d2u/dxdy +   dN/dx du/dy - beta_b v = -tau_dy
    !
    ! Rearranging to gather the terms involving u and v gives:
    !
    !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_b u + ...
    !   3 N d2v/dxdy + 2 dN/dx dv/dy +             dN/dy dv/dx = -tau_dx
    !
    !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_b v + ...
    !   3 N d2u/dxdy + 2 dN/dy du/dx +             dN/dx du/dy = -tau_dy

    ! In/output variables:
    type(type_graph_pair),                                                    intent(in   ) :: graphs
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: N_b, dN_dx_b, dN_dy_b
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: basal_friction_coefficient_b
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: tau_dx_b, tau_dy_b
    type(type_sparse_matrix_CSR_dp),                                          intent(inout) :: A_CSR
    real(dp), dimension(graphs%graph_b%ni1*2-1: graphs%graph_b%ni2*2),        intent(inout) :: bb_b
    integer,                                                                  intent(in   ) :: row_niuv

    ! Local variables:
    integer                             :: ni, uv
    real(dp)                            :: N, dN_dx, dN_dy, basal_friction_coefficient, tau_dx, tau_dy
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_ddx_val
    real(dp), dimension(:), allocatable :: single_row_ddy_val
    real(dp), dimension(:), allocatable :: single_row_d2dx2_val
    real(dp), dimension(:), allocatable :: single_row_d2dxdy_val
    real(dp), dimension(:), allocatable :: single_row_d2dy2_val
    integer                             :: single_row_nnz
    real(dp)                            :: Au, Av
    integer                             :: k, nj, col_nju, col_njv

    ! Relevant indices for this graph node
    ni = graphs%n2biuv( row_niuv,1)
    uv = graphs%n2biuv( row_niuv,2)

    ! N, dN/dx, dN/dy, basal_friction_coefficient_b, tau_dx, and tau_dy on this graph node
    N                          = N_b                         ( ni)
    dN_dx                      = dN_dx_b                     ( ni)
    dN_dy                      = dN_dy_b                     ( ni)
    basal_friction_coefficient = basal_friction_coefficient_b( ni)
    tau_dx                     = tau_dx_b                    ( ni)
    tau_dy                     = tau_dy_b                    ( ni)

    ! allocate memory for single matrix rows
    allocate( single_row_ind(        graphs%graph_a%nC_mem*2))
    allocate( single_row_ddx_val(    graphs%graph_a%nC_mem*2))
    allocate( single_row_ddy_val(    graphs%graph_a%nC_mem*2))
    allocate( single_row_d2dx2_val(  graphs%graph_a%nC_mem*2))
    allocate( single_row_d2dxdy_val( graphs%graph_a%nC_mem*2))
    allocate( single_row_d2dy2_val(  graphs%graph_a%nC_mem*2))

    ! Read coefficients of the operator matrices
    call read_single_row_CSR_dist( graphs%M2_ddx_b_b   , ni, single_row_ind, single_row_ddx_val   , single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_ddy_b_b   , ni, single_row_ind, single_row_ddy_val   , single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_d2dx2_b_b , ni, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_d2dxdy_b_b, ni, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_d2dy2_b_b , ni, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    if (uv == 1) then
      ! x-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        nj      = single_row_ind( k)
        col_nju = graphs%biuv2n( nj,1)
        col_njv = graphs%biuv2n( nj,2)

        !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_b u + ...
        !   3 N d2v/dxdy + 2 dN/dx dv/dy +             dN/dy dv/dx = -tau_dx

        ! Combine the mesh operators
        Au = 4._dp * N     * single_row_d2dx2_val(  k) + &  ! 4  N    d2u/dx2
             4._dp * dN_dx * single_row_ddx_val(    k) + &  ! 4 dN/dx du/dx
                     N     * single_row_d2dy2_val(  k) + &  !    N    d2u/dy2
                     dN_dy * single_row_ddy_val(    k)      !   dN/dy du/dy
        if (nj == ni) Au = Au - basal_friction_coefficient  ! - beta_b u

        Av = 3._dp * N     * single_row_d2dxdy_val( k) + &  ! 3  N    d2v/dxdy
             2._dp * dN_dx * single_row_ddy_val(    k) + &  ! 2 dN/dx dv/dy
                     dN_dy * single_row_ddx_val(    k)      !   dN/dy dv/dx

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_niuv, col_nju, Au)
        call add_entry_CSR_dist( A_CSR, row_niuv, col_njv, Av)

      end do

      ! Load vector
      bb_b( row_niuv) = -tau_dx

    elseif (uv == 2) then
      ! y-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        nj      = single_row_ind( k)
        col_nju = graphs%biuv2n( nj,1)
        col_njv = graphs%biuv2n( nj,2)

        !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_b v + ...
        !   3 N d2u/dxdy + 2 dN/dy du/dx +             dN/dx du/dy = -tau_dy

        ! Combine the mesh operators
        Av = 4._dp * N     * single_row_d2dy2_val(  k) + &  ! 4  N    d2v/dy2
             4._dp * dN_dy * single_row_ddy_val(    k) + &  ! 4 dN/dy dv/dy
                     N     * single_row_d2dx2_val(  k) + &  !    N    d2v/dx2
                     dN_dx * single_row_ddx_val(    k)      !   dN/dx dv/dx
        if (nj == ni) Av = Av - basal_friction_coefficient  ! - beta_b v

        Au = 3._dp * N     * single_row_d2dxdy_val( k) + &  ! 3  N    d2u/dxdy
             2._dp * dN_dy * single_row_ddx_val(    k) + &  ! 2 dN/dy du/dx
                     dN_dx * single_row_ddy_val(    k)      !   dN/dx du/dy

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_niuv, col_nju, Au)
        call add_entry_CSR_dist( A_CSR, row_niuv, col_njv, Av)

      end do

      ! Load vector
      bb_b( row_niuv) = -tau_dy

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_SSA_DIVA_stiffness_matrix_row_free

  subroutine calc_SSA_DIVA_sans_stiffness_matrix_row_free( graphs, N_b, basal_friction_coefficient_b, &
    tau_dx_b, tau_dy_b, A_CSR, bb_b, row_niuv)
    !< Add coefficients to this matrix row to represent the linearised SSA
    !< sans the gradients of the effective viscosity (the "cross-terms")

    ! The SSA reads;
    !
    !   d/dx [ 2 N ( 2 du/dx + dv/dy )] + d/dy [ N ( du/dy + dv/dx)] - beta_b u = -tau_dx
    !
    !   d/dy [ 2 N ( 2 dv/dy + du/dx )] + d/dx [ N ( dv/dx + du/dy)] - beta_b v = -tau_dy
    !
    ! Using the chain rule, this expands to read:
    !
    !   4 N d2u/dx2 + 4 dN/dx du/dx + 2 N d2v/dxdy + 2 dN/dx dv/dy + ...
    !     N d2u/dy2 +   dN/dy du/dy +   N d2v/dxdy +   dN/dy dv/dx - beta_b u = -tau_dx
    !
    !   4 N d2v/dy2 + 4 dN/dy dv/dy + 2 N d2u/dxdy + 2 dN/dy du/dx + ...
    !     N d2v/dx2 +   dN/dx dv/dx +   N d2u/dxdy +   dN/dx du/dy - beta_b v = -tau_dy
    !
    ! The "sans" approximation neglects the gradients dN/dx, dN/dy of N:
    !
    !   4 N d2u/dx2 + N d2u/dy2 + 3 N d2v/dxdy - beta_b u = -tau_dx
    !   4 N d2v/dy2 + N d2v/dx2 + 3 N d2u/dxdy - beta_b v = -tau_dy
    !
    ! Dividing both sides by N yields:
    !
    !   4 d2u/dx2 + d2u/dy2 + 3 d2v/dxdy - beta_b u / N = -tau_dx / N
    !   4 d2v/dy2 + d2v/dx2 + 3 d2u/dxdy - beta_b v / N = -tau_dy / N
    !
    ! Note that there is no clear mathematical or physical reason why this should be allowed.
    ! However, while I (Tijn Berends, 2023) have found a few cases where there are noticeable
    ! differences in the solutions (e.g. ISMIP-HOM experiments with high strain rates),
    ! most of the time the difference with respect to the full SSA/DIVA is very small.
    ! The "sans" option makes the solver quite a lot more stable and therefore faster.
    ! Someone really ought to perform some proper experiments to determine whether or not
    ! this should be the default.

    ! In/output variables:
    type(type_graph_pair),                                                    intent(in   ) :: graphs
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: N_b
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: basal_friction_coefficient_b
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: tau_dx_b, tau_dy_b
    type(type_sparse_matrix_CSR_dp),                                          intent(inout) :: A_CSR
    real(dp), dimension(graphs%graph_b%ni1*2-1: graphs%graph_b%ni2*2),        intent(inout) :: bb_b
    integer,                                                                  intent(in   ) :: row_niuv

    ! Local variables:
    integer                             :: ni, uv
    real(dp)                            :: N, basal_friction_coefficient, tau_dx, tau_dy
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_d2dx2_val
    real(dp), dimension(:), allocatable :: single_row_d2dxdy_val
    real(dp), dimension(:), allocatable :: single_row_d2dy2_val
    integer                             :: single_row_nnz
    real(dp)                            :: Au, Av
    integer                             :: k, nj, col_nju, col_njv

    ! Relevant indices for this triangle
    ni = graphs%n2biuv( row_niuv,1)
    uv = graphs%n2biuv( row_niuv,2)

    ! N, beta_b, tau_dx, and tau_dy on this graph node
    N                          = N_b                         ( ni)
    basal_friction_coefficient = basal_friction_coefficient_b( ni)
    tau_dx                     = tau_dx_b                    ( ni)
    tau_dy                     = tau_dy_b                    ( ni)

    ! allocate memory for single matrix rows
    allocate( single_row_ind(        graphs%graph_a%nC_mem*2))
    allocate( single_row_d2dx2_val(  graphs%graph_a%nC_mem*2))
    allocate( single_row_d2dxdy_val( graphs%graph_a%nC_mem*2))
    allocate( single_row_d2dy2_val(  graphs%graph_a%nC_mem*2))

    ! Read coefficients of the operator matrices
    call read_single_row_CSR_dist( graphs%M2_d2dx2_b_b , ni, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_d2dxdy_b_b, ni, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_d2dy2_b_b , ni, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    if (uv == 1) then
      ! x-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        nj      = single_row_ind( k)
        col_nju = graphs%biuv2n( nj,1)
        col_njv = graphs%biuv2n( nj,2)

        !   4 d2u/dx2 + d2u/dy2 + 3 d2v/dxdy - beta_b u / N = -tau_dx / N

        ! Combine the mesh operators
        Au = 4._dp * single_row_d2dx2_val(  k) + &             ! 4 d2u/dx2
                     single_row_d2dy2_val(  k)                 !   d2u/dy2
        if (nj == ni) Au = Au - basal_friction_coefficient / N ! - beta_b u / N

        Av = 3._dp * single_row_d2dxdy_val( k)                 ! 3 d2v/dxdy

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_niuv, col_nju, Au)
        call add_entry_CSR_dist( A_CSR, row_niuv, col_njv, Av)

      end do

      ! Load vector
      bb_b( row_niuv) = -tau_dx / N

    elseif (uv == 2) then
      ! y-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        nj      = single_row_ind( k)
        col_nju = graphs%biuv2n( nj,1)
        col_njv = graphs%biuv2n( nj,2)

        !   4 d2v/dy2 + d2v/dx2 + 3 d2u/dxdy - beta_b v / N = -tau_dy / N

        ! Combine the mesh operators
        Av = 4._dp * single_row_d2dy2_val(  k) + &             ! 4 d2v/dy2
                     single_row_d2dx2_val(  k)                 !   d2v/dx2
        if (nj == ni) Av = Av - basal_friction_coefficient / N ! - beta_b v / N

        Au = 3._dp * single_row_d2dxdy_val( k)                 ! 3 d2u/dxdy

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_niuv, col_nju, Au)
        call add_entry_CSR_dist( A_CSR, row_niuv, col_njv, Av)

      end do

      ! Load vector
      bb_b( row_niuv) = -tau_dy / N

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_SSA_DIVA_sans_stiffness_matrix_row_free

  subroutine calc_SSA_DIVA_stiffness_matrix_row_BC_ice_front( graphs, N_b, tau_ox_b, tau_oy_b, &
    A_CSR, bb_b, row_niuv)
    ! The ocean-pressure boundary condition to the SSA at the ice front reads;
    !
    !   2 N ( 2 du/dx + dv/dy ) n_x + N ( du/dy + dv/dx) n_y = tau_ox
    !
    !   2 N ( 2 dv/dy + du/dx ) n_y + N ( dv/dx + du/dy) n_x = tau_oy
    !
    ! (See also Robinson et al., 2020, Eq. 19)
    !
    ! Rearranging to gather the terms involving u and v gives:
    !
    !   4 N n_x du/dx + N n_y du/dy + 2 N n_x dv/dy + N n_y dv/dx = tau_ox
    !
    !   4 N n_y dv/dy + N n_x dv/dx + 2 N n_y du/dx + N n_x du/dy = tau_oy
    !
    ! The ocean back pressure at the ice front is defined as:
    !
    !   tau_ox = (1/2 rho_i g H^2 - 1/2 rho_sw g Ho^2) n_x
    !   tau_oy = (1/2 rho_i g H^2 - 1/2 rho_sw g Ho^2) n_y
    !
    ! The height of the water column in contact with the ice front is defined as:
    !
    !   Ho = min( max( SL - Hb, 0), rho_i / rho_sw H)

    ! In/output variables:
    type(type_graph_pair),                                                    intent(in   ) :: graphs
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: N_b, tau_ox_b, tau_oy_b
    type(type_sparse_matrix_CSR_dp),                                          intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2),                                   intent(inout) :: bb_b
    integer,                                                                  intent(in   ) :: row_niuv

    ! Local variables:
    integer                             :: ni, uv
    real(dp)                            :: N, tau_ox, tau_oy, n_x, n_y
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_ddx_val
    real(dp), dimension(:), allocatable :: single_row_ddy_val
    integer                             :: single_row_nnz
    real(dp)                            :: Au, Av
    integer                             :: k, nj, col_nju, col_njv

    ! Relevant indices for this graph node
    ni = graphs%n2biuv( row_niuv,1)
    uv = graphs%n2biuv( row_niuv,2)

    ! N, Hi, Ho, and outward unit normal vector on this graph node
    N       = N_b      ( ni)
    tau_ox = tau_ox_b( ni)
    tau_oy = tau_oy_b( ni)
    n_x     = graphs%graph_b%border_nhat( ni,1)
    n_y     = graphs%graph_b%border_nhat( ni,2)

    ! Allocate memory for single matrix rows
    allocate( single_row_ind(     graphs%graph_a%nC_mem*2))
    allocate( single_row_ddx_val( graphs%graph_a%nC_mem*2))
    allocate( single_row_ddy_val( graphs%graph_a%nC_mem*2))

    ! Read coefficients of the operator matrices
    call read_single_row_CSR_dist( graphs%M_ddx_b_b, ni, single_row_ind, single_row_ddx_val, single_row_nnz)
    call read_single_row_CSR_dist( graphs%M_ddy_b_b, ni, single_row_ind, single_row_ddy_val, single_row_nnz)

    if (uv == 1) then
      ! x-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        nj      = single_row_ind( k)
        col_nju = graphs%biuv2n( nj,1)
        col_njv = graphs%biuv2n( nj,2)

        !   4 N n_x du/dx + N n_y du/dy + 2 N n_x dv/dy + N n_y dv/dx = tau_ox

        ! Combine the mesh operators
        Au = 4._dp * N * n_x * single_row_ddx_val( k) + &  ! 4 N n_x du/dx
                     N * n_y * single_row_ddy_val( k)      !   N n_y du/dy

        Av = 2._dp * N * n_x * single_row_ddy_val( k) + &  ! 2 N n_x dv/dy
                     N * n_y * single_row_ddx_val( k)      !   N n_y dv/dx

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_niuv, col_nju, Au)
        call add_entry_CSR_dist( A_CSR, row_niuv, col_njv, Av)

      end do

      ! Load vector
      bb_b( row_niuv) = tau_ox

    elseif (uv == 2) then
      ! y-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        nj      = single_row_ind( k)
        col_nju = graphs%biuv2n( nj,1)
        col_njv = graphs%biuv2n( nj,2)

        !   4 N n_y dv/dy + N n_x dv/dx + 2 N n_y du/dx + N n_x du/dy = tau_oy

        ! Combine the mesh operators
        Av = 4._dp * N * n_y * single_row_ddy_val( k) + &  ! 4 N n_y dv/dy
                     N * n_x * single_row_ddx_val( k)      !   N n_x dv/dx

        Au = 2._dp * N * n_y * single_row_ddx_val( k) + &  ! 2 N n_y du/dx
                     N * n_x * single_row_ddy_val( k)      !   N n_x du/dy

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_niuv, col_nju, Au)
        call add_entry_CSR_dist( A_CSR, row_niuv, col_njv, Av)

      end do

      ! Load vector
      bb_b( row_niuv) = tau_oy

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_SSA_DIVA_stiffness_matrix_row_BC_ice_front

end module solve_linearised_SSA_DIVA_ocean_pressure
