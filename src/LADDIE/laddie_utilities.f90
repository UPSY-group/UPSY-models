MODULE laddie_utilities

  ! Utilities for the laddie model

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE call_stack_and_comp_time_tracking                  , ONLY: crash, init_routine, finalise_routine, warning
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE laddie_model_types                                     , ONLY: type_laddie_model, type_laddie_timestep
  use laddie_forcing_types, only: type_laddie_forcing
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE ocean_utilities                                        , ONLY: interpolate_ocean_depth
  USE mpi_distributed_memory                                 , ONLY: gather_to_all
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_1D
  use mesh_integrate_over_domain, only: average_over_domain
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared
  use checksum_mod, only: checksum

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE compute_ambient_TS( mesh, laddie, forcing, Hstar)
    ! Compute T and S of ambient ocean water at the depth of LADDIE's layer bottom
    ! through vertical interpolation

    ! In- and output variables

    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
    type(type_laddie_forcing),              intent(in)    :: forcing
    REAL(dp), DIMENSION(mesh%pai_V%i1_nih:mesh%pai_V%i2_nih), INTENT(IN)    :: Hstar

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_ambient_TS'
    INTEGER                                               :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Get T and S at layer base
    DO vi = mesh%vi1, mesh%vi2
       IF (laddie%mask_a( vi)) THEN
         CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, forcing%T_ocean( vi,:), Hstar( vi) - forcing%Hib( vi), laddie%T_amb( vi))
         CALL interpolate_ocean_depth( C%nz_ocean, C%z_ocean, forcing%S_ocean( vi,:), Hstar( vi) - forcing%Hib( vi), laddie%S_amb( vi))
       END IF
    END DO
    call checksum( mesh%pai_V, laddie%T_amb, 'laddie%T_amb')
    call checksum( mesh%pai_V, laddie%S_amb, 'laddie%S_amb')

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE compute_ambient_TS

  subroutine map_H_a_b( mesh, laddie, H_a, H_b)
    ! Map layer thickness from a to b grid, accounting for BCs

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(in)    :: laddie
    real(dp), dimension(:),                 intent(in)    :: H_a
    real(dp), dimension(:),                 intent(inout) :: H_b

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'map_H_a_b'

    ! Add routine to path
    call init_routine( routine_name)

    call sync

    call multiply_CSR_matrix_with_vector_1D( laddie%M_map_H_a_b, &
      mesh%pai_V, H_a, mesh%pai_Tri, H_b)
    call checksum( mesh%pai_Tri, H_b, 'H_b')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_H_a_b

  subroutine map_H_a_c( mesh, laddie, H_a, H_c)
    ! Map layer thickness from a to c grid, accounting for BCs

    ! In- and output variables

    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_model),                intent(in)    :: laddie
    real(dp), dimension(:),                 intent(in)    :: H_a
    real(dp), dimension(:),                 intent(inout) :: H_c

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'map_H_a_c'

    ! Add routine to path
    call init_routine( routine_name)

    call sync

    call multiply_CSR_matrix_with_vector_1D( laddie%M_map_H_a_c, &
      mesh%pai_V, H_a, mesh%pai_E, H_c)
    call checksum( mesh%pai_E, H_c, 'H_c')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_H_a_c

  SUBROUTINE allocate_laddie_model( mesh, laddie)
    ! Allocate variables of the laddie model

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_model),                INTENT(INOUT) :: laddie

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'allocate_laddie_model'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Thickness
    if (associated( laddie%dH_dt )) call deallocate_dist_shared( laddie%dH_dt , laddie%wdH_dt )
    call allocate_dist_shared( laddie%dH_dt         , laddie%wdH_dt         , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [m]             change

    ! Temperatures
    if (associated( laddie%T_amb    )) call deallocate_dist_shared( laddie%T_amb    , laddie%wT_amb    )
    if (associated( laddie%T_base   )) call deallocate_dist_shared( laddie%T_base   , laddie%wT_base   )
    if (associated( laddie%T_freeze )) call deallocate_dist_shared( laddie%T_freeze , laddie%wT_freeze )
    call allocate_dist_shared( laddie%T_amb         , laddie%wT_amb         , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [degC]          Temperature layer bottom
    call allocate_dist_shared( laddie%T_base        , laddie%wT_base        , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [degC]          Temperature ice shelf base
    call allocate_dist_shared( laddie%T_freeze      , laddie%wT_freeze      , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [degC]          Temperature freezing

    ! Salinities
    if (associated( laddie%S_amb  )) call deallocate_dist_shared( laddie%S_amb  , laddie%wS_amb  )
    if (associated( laddie%S_base )) call deallocate_dist_shared( laddie%S_base , laddie%wS_base )
    call allocate_dist_shared( laddie%S_amb         , laddie%wS_amb         , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [PSU]           Salinity layer bottom
    call allocate_dist_shared( laddie%S_base        , laddie%wS_base        , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [PSU]           Salinity ice shelf base

    ! Densities and buoyancies
    if (associated( laddie%rho         )) call deallocate_dist_shared( laddie%rho         , laddie%wrho         )
    if (associated( laddie%rho_amb     )) call deallocate_dist_shared( laddie%rho_amb     , laddie%wrho_amb     )
    if (associated( laddie%drho_amb    )) call deallocate_dist_shared( laddie%drho_amb    , laddie%wdrho_amb    )
    if (associated( laddie%Hdrho_amb   )) call deallocate_dist_shared( laddie%Hdrho_amb   , laddie%wHdrho_amb   )
    if (associated( laddie%Hdrho_amb_b )) call deallocate_dist_shared( laddie%Hdrho_amb_b , laddie%wHdrho_amb_b )
    if (associated( laddie%drho_base   )) call deallocate_dist_shared( laddie%drho_base   , laddie%wdrho_base   )
    call allocate_dist_shared( laddie%rho           , laddie%wrho           , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [kg m^-3]       Layer density
    call allocate_dist_shared( laddie%rho_amb       , laddie%wrho_amb       , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [kg m^-3]       Ambient water density
    call allocate_dist_shared( laddie%drho_amb      , laddie%wdrho_amb      , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! []              Buoyancy at layer bottom
    call allocate_dist_shared( laddie%Hdrho_amb     , laddie%wHdrho_amb     , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! []              Depth-integrated buoyancy at layer bottom
    call allocate_dist_shared( laddie%Hdrho_amb_b   , laddie%wHdrho_amb_b   , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    ! []              Depth-integrated buoyancy at layer bottom
    call allocate_dist_shared( laddie%drho_base     , laddie%wdrho_base     , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! []              Buoyancy at ice base

    ! Friction velocity
    if (associated( laddie%u_star )) call deallocate_dist_shared( laddie%u_star , laddie%wu_star )
    call allocate_dist_shared( laddie%u_star        , laddie%wu_star        , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [m s^-1]        Friction velocity

    ! Physical parameter fields
    if (associated( laddie%gamma_T )) call deallocate_dist_shared( laddie%gamma_T , laddie%wgamma_T )
    if (associated( laddie%gamma_S )) call deallocate_dist_shared( laddie%gamma_S , laddie%wgamma_S )
    if (associated( laddie%A_h     )) call deallocate_dist_shared( laddie%A_h     , laddie%wA_h     )
    if (associated( laddie%K_h     )) call deallocate_dist_shared( laddie%K_h     , laddie%wK_h     )
    call allocate_dist_shared( laddie%gamma_T       , laddie%wgamma_T       , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! []              Turbulent heat exchange coefficient
    call allocate_dist_shared( laddie%gamma_S       , laddie%wgamma_S       , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! []              Turbulent salt exchange coefficient
    call allocate_dist_shared( laddie%A_h           , laddie%wA_h           , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    ! [m^2 s^-1]      Horizontal laplacian viscosity
    call allocate_dist_shared( laddie%K_h           , laddie%wK_h           , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [m^2 s^-1]      Horizontal diffusivity

    ! Vertical rates
    if (associated( laddie%melt      )) call deallocate_dist_shared( laddie%melt      , laddie%wmelt      )
    if (associated( laddie%entr      )) call deallocate_dist_shared( laddie%entr      , laddie%wentr      )
    if (associated( laddie%entr_dmin )) call deallocate_dist_shared( laddie%entr_dmin , laddie%wentr_dmin )
    if (associated( laddie%detr      )) call deallocate_dist_shared( laddie%detr      , laddie%wdetr      )
    if (associated( laddie%entr_tot  )) call deallocate_dist_shared( laddie%entr_tot  , laddie%wentr_tot  )
    if (associated( laddie%SGD       )) call deallocate_dist_shared( laddie%SGD       , laddie%wSGD       )
    call allocate_dist_shared( laddie%melt          , laddie%wmelt          , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [m s^-1]        Melting / freezing rate
    call allocate_dist_shared( laddie%entr          , laddie%wentr          , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [m s^-1]        Entrainment
    call allocate_dist_shared( laddie%entr_dmin     , laddie%wentr_dmin     , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [m s^-1]        Entrainment for D_min
    call allocate_dist_shared( laddie%detr          , laddie%wdetr          , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [m s^-1]        Detrainment
    call allocate_dist_shared( laddie%entr_tot      , laddie%wentr_tot      , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [m s^-1]        Total (net) entrainment
    call allocate_dist_shared( laddie%SGD           , laddie%wSGD           , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [m s^-1]        Subglacial discharge

    ! Horizontal fluxes
    if (associated( laddie%divQH )) call deallocate_dist_shared( laddie%divQH , laddie%wdivQH )
    if (associated( laddie%divQU )) call deallocate_dist_shared( laddie%divQU , laddie%wdivQU )
    if (associated( laddie%divQV )) call deallocate_dist_shared( laddie%divQV , laddie%wdivQV )
    if (associated( laddie%divQT )) call deallocate_dist_shared( laddie%divQT , laddie%wdivQT )
    if (associated( laddie%divQS )) call deallocate_dist_shared( laddie%divQS , laddie%wdivQS )
    call allocate_dist_shared( laddie%divQH         , laddie%wdivQH         , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [m^3 s^-1]      Divergence of layer thickness
    call allocate_dist_shared( laddie%divQU         , laddie%wdivQU         , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    ! [m^4 s^-2]      Divergence of momentum
    call allocate_dist_shared( laddie%divQV         , laddie%wdivQV         , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    ! [m^4 s^-2]
    call allocate_dist_shared( laddie%divQT         , laddie%wdivQT         , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [degC m^3 s^-1] Divergence of heat
    call allocate_dist_shared( laddie%divQS         , laddie%wdivQS         , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [PSU m^3 s^-1]  Divergence of salt

    ! Viscosities
    if (associated( laddie%viscU )) call deallocate_dist_shared( laddie%viscU , laddie%wviscU )
    if (associated( laddie%viscV )) call deallocate_dist_shared( laddie%viscV , laddie%wviscV )
    call allocate_dist_shared( laddie%viscU         , laddie%wviscU         , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    ! [m^2 s^-2]      Horizontal viscosity term
    call allocate_dist_shared( laddie%viscV         , laddie%wviscV         , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    ! [m^2 s^-2]

    ! Diffusivities
    if (associated( laddie%diffT )) call deallocate_dist_shared( laddie%diffT , laddie%wdiffT )
    if (associated( laddie%diffS )) call deallocate_dist_shared( laddie%diffS , laddie%wdiffS )
    call allocate_dist_shared( laddie%diffT         , laddie%wdiffT         , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [degC m s^-1]   Horizontal diffusivity of heat
    call allocate_dist_shared( laddie%diffS         , laddie%wdiffS         , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [PSU m s^-1]    Horizontal diffusivity of salt

    ! RHS terms
    if (associated( laddie%ddrho_amb_dx_b )) call deallocate_dist_shared( laddie%ddrho_amb_dx_b , laddie%wddrho_amb_dx_b )
    if (associated( laddie%ddrho_amb_dy_b )) call deallocate_dist_shared( laddie%ddrho_amb_dy_b , laddie%wddrho_amb_dy_b )
    if (associated( laddie%dH_dx_b        )) call deallocate_dist_shared( laddie%dH_dx_b        , laddie%wdH_dx_b        )
    if (associated( laddie%dH_dy_b        )) call deallocate_dist_shared( laddie%dH_dy_b        , laddie%wdH_dy_b        )
    if (associated( laddie%detr_b         )) call deallocate_dist_shared( laddie%detr_b         , laddie%wdetr_b         )
    call allocate_dist_shared( laddie%ddrho_amb_dx_b, laddie%wddrho_amb_dx_b, [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    ! [m^-1]          Horizontal derivative of buoyancy
    call allocate_dist_shared( laddie%ddrho_amb_dy_b, laddie%wddrho_amb_dy_b, [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    ! [m^-1]
    call allocate_dist_shared( laddie%dH_dx_b       , laddie%wdH_dx_b       , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    ! [m^-2]          Horizontal derivative of thickness
    call allocate_dist_shared( laddie%dH_dy_b       , laddie%wdH_dy_b       , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    ! [m^-2]
    call allocate_dist_shared( laddie%detr_b        , laddie%wdetr_b        , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    ! [m s^-1]        Detrainment on b grid

    ! Forward-Backward Runge-Kutta 3 scheme
    if (associated( laddie%Hstar )) call deallocate_dist_shared( laddie%Hstar , laddie%wHstar )
    call allocate_dist_shared( laddie%Hstar         , laddie%wHstar         , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! [m]               Intermediate layer thickness
    laddie%Hstar         ( mesh%pai_V%i1_nih  :mesh%pai_V%i2_nih  ) => laddie%Hstar

    ! Mapped variables
    if (associated( laddie%H_c     )) call deallocate_dist_shared( laddie%H_c     , laddie%wH_c     )
    if (associated( laddie%Hstar_b )) call deallocate_dist_shared( laddie%Hstar_b , laddie%wHstar_b )
    if (associated( laddie%Hstar_c )) call deallocate_dist_shared( laddie%Hstar_c , laddie%wHstar_c )
    call allocate_dist_shared( laddie%H_c           , laddie%wH_c           , [mesh%pai_E%i1_nih, mesh%pai_E%i2_nih])
    call allocate_dist_shared( laddie%Hstar_b       , laddie%wHstar_b       , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])
    call allocate_dist_shared( laddie%Hstar_c       , laddie%wHstar_c       , [mesh%pai_E%i1_nih, mesh%pai_E%i2_nih])

    ! Masks
    if (associated( laddie%mask_a    )) call deallocate_dist_shared( laddie%mask_a    , laddie%wmask_a    )
    if (associated( laddie%mask_gr_a )) call deallocate_dist_shared( laddie%mask_gr_a , laddie%wmask_gr_a )
    if (associated( laddie%mask_oc_a )) call deallocate_dist_shared( laddie%mask_oc_a , laddie%wmask_oc_a )
    if (associated( laddie%mask_b    )) call deallocate_dist_shared( laddie%mask_b    , laddie%wmask_b    )
    if (associated( laddie%mask_gl_b )) call deallocate_dist_shared( laddie%mask_gl_b , laddie%wmask_gl_b )
    if (associated( laddie%mask_cf_b )) call deallocate_dist_shared( laddie%mask_cf_b , laddie%wmask_cf_b )
    if (associated( laddie%mask_oc_b )) call deallocate_dist_shared( laddie%mask_oc_b , laddie%wmask_oc_b )
    call allocate_dist_shared( laddie%mask_a        , laddie%wmask_a        , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    !                 Mask on a-grid
    call allocate_dist_shared( laddie%mask_gr_a     , laddie%wmask_gr_a     , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    !                 Grounded mask on a-grid
    call allocate_dist_shared( laddie%mask_oc_a     , laddie%wmask_oc_a     , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    !                 Icefree ocean mask on a-grid
    call allocate_dist_shared( laddie%mask_b        , laddie%wmask_b        , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    !                 Mask on b-grid
    call allocate_dist_shared( laddie%mask_gl_b     , laddie%wmask_gl_b     , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    !                 Grounding line mask on b-grid
    call allocate_dist_shared( laddie%mask_cf_b     , laddie%wmask_cf_b     , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    !                 Calving front mask on b-grid
    call allocate_dist_shared( laddie%mask_oc_b     , laddie%wmask_oc_b     , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    !                 Icefree ocean mask on b-grid

    ! Domains
    if (associated( laddie%domain_a )) call deallocate_dist_shared( laddie%domain_a , laddie%wdomain_a )
    if (associated( laddie%domain_b )) call deallocate_dist_shared( laddie%domain_b , laddie%wdomain_b )
    call allocate_dist_shared( laddie%domain_a      , laddie%wdomain_a      , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])    ! []              Floating domain on a-grid
    call allocate_dist_shared( laddie%domain_b      , laddie%wdomain_b      , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])    ! []              Floating domain on b-grid

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_laddie_model

  SUBROUTINE allocate_laddie_timestep( mesh, npx)
    ! Allocate variables of the laddie model

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_laddie_timestep),             INTENT(INOUT) :: npx

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'allocate_laddie_timestep'

    ! Add routine to path
    CALL init_routine( routine_name)

    if (associated( npx%H   )) call deallocate_dist_shared( npx%H   , npx%wH   )
    if (associated( npx%H_b )) call deallocate_dist_shared( npx%H_b , npx%wH_b )
    if (associated( npx%H_c )) call deallocate_dist_shared( npx%H_c , npx%wH_c )
    if (associated( npx%U   )) call deallocate_dist_shared( npx%U   , npx%wU   )
    if (associated( npx%U_a )) call deallocate_dist_shared( npx%U_a , npx%wU_a )
    if (associated( npx%U_c )) call deallocate_dist_shared( npx%U_c , npx%wU_c )
    if (associated( npx%V   )) call deallocate_dist_shared( npx%V   , npx%wV   )
    if (associated( npx%V_a )) call deallocate_dist_shared( npx%V_a , npx%wV_a )
    if (associated( npx%V_c )) call deallocate_dist_shared( npx%V_c , npx%wV_c )
    if (associated( npx%T   )) call deallocate_dist_shared( npx%T   , npx%wT   )
    if (associated( npx%S   )) call deallocate_dist_shared( npx%S   , npx%wS   )

    call allocate_dist_shared( npx%H  , npx%wH  , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])   ! [m]             Layer thickness
    call allocate_dist_shared( npx%H_b, npx%wH_b, [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])   ! [m]             Layer thickness on b grid
    call allocate_dist_shared( npx%H_c, npx%wH_c, [mesh%pai_E%i1_nih, mesh%pai_E%i2_nih])   ! [m]             Layer thickness on c grid
    call allocate_dist_shared( npx%U  , npx%wU  , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])   ! [m s^-1]        2D velocity
    call allocate_dist_shared( npx%U_a, npx%wU_a, [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])   ! [m s^-1]        2D velocity on a grid
    call allocate_dist_shared( npx%U_c, npx%wU_c, [mesh%pai_E%i1_nih, mesh%pai_E%i2_nih])   ! [m s^-1]        2D velocity on b grid
    call allocate_dist_shared( npx%V  , npx%wV  , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])   ! [m s^-1]
    call allocate_dist_shared( npx%V_a, npx%wV_a, [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])   ! [m s^-1]
    call allocate_dist_shared( npx%V_c, npx%wV_c, [mesh%pai_E%i1_nih, mesh%pai_E%i2_nih])   ! [m s^-1]
    call allocate_dist_shared( npx%T  , npx%wT  , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])   ! [degC]          Temperature
    call allocate_dist_shared( npx%S  , npx%wS  , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])   ! [PSU]           Salinity

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_laddie_timestep

  subroutine allocate_laddie_forcing( mesh, forcing)
    ! Allocate variables of the laddie model

    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_laddie_forcing),              intent(inout) :: forcing

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'allocate_laddie_forcing'

    ! Add routine to path
    call init_routine( routine_name)

    ! Forcing
    call allocate_dist_shared( forcing%Hi                , forcing%wHi                , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])
    call allocate_dist_shared( forcing%Hs                , forcing%wHs                , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])
    call allocate_dist_shared( forcing%Hb                , forcing%wHb                , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])
    call allocate_dist_shared( forcing%Hib               , forcing%wHib               , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])
    call allocate_dist_shared( forcing%TAF               , forcing%wTAF               , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])
    call allocate_dist_shared( forcing%dHib_dx_b         , forcing%wdHib_dx_b         , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])
    call allocate_dist_shared( forcing%dHib_dy_b         , forcing%wdHib_dy_b         , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])
    call allocate_dist_shared( forcing%mask_icefree_land , forcing%wmask_icefree_land , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])
    call allocate_dist_shared( forcing%mask_icefree_ocean, forcing%wmask_icefree_ocean, [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])
    call allocate_dist_shared( forcing%mask_grounded_ice , forcing%wmask_grounded_ice , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])
    call allocate_dist_shared( forcing%mask_floating_ice , forcing%wmask_floating_ice , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])
    call allocate_dist_shared( forcing%mask_gl_fl        , forcing%wmask_gl_fl        , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])
    call allocate_dist_shared( forcing%mask_SGD          , forcing%wmask_SGD          , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])
    call allocate_dist_shared( forcing%mask              , forcing%wmask              , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih])
    call allocate_dist_shared( forcing%Ti                , forcing%wTi                , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih], [1, mesh%nz])
    call allocate_dist_shared( forcing%T_ocean           , forcing%wT_ocean           , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih], [1, C%nz_ocean])
    call allocate_dist_shared( forcing%S_ocean           , forcing%wS_ocean           , [mesh%pai_V%i1_nih, mesh%pai_V%i2_nih], [1, C%nz_ocean])
    call allocate_dist_shared( forcing%f_coriolis        , forcing%wf_coriolis        , [mesh%pai_Tri%i1_nih, mesh%pai_Tri%i2_nih])

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_laddie_forcing

END MODULE laddie_utilities

