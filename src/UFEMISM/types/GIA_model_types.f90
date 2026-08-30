module GIA_model_types

  use precisions, only: dp
  use grid_types, only: type_grid

  implicit none

  private

  public :: type_GIA_model, type_ELRA_model

  type type_GIA_model
    ! The GIA model data structure.

    ! The GIA model grid
    type(type_grid)                     :: grid

    ! Main data fields
    real(dp), dimension(:), allocatable :: relative_surface_load_mesh  ! [Pa] Surface load relative to the GIA-equilibrium reference geometry, on the UFEMISM mesh
    real(dp), dimension(:), allocatable :: relative_surface_load_grid  ! [Pa] Surface load relative to the GIA-equilibrium reference geometry, on the UFEMISM mesh

    ! Sub-models

    ! Timestepping
    real(dp)                            :: t_prev                      ! [yr] Time of the previous state
    real(dp)                            :: t_next                      ! [yr] Time of the next state
    real(dp), dimension(:), allocatable :: dHb_prev                    ! [m]  The previous state (bedrock deflection relative to the GIA-equilibrium reference geometry)
    real(dp), dimension(:), allocatable :: dHb_next                    ! [m]  The next state     (idem)

  end type type_GIA_model

  type type_ELRA_model
    ! The ELRA model data structure.

    real(dp), dimension(:  ), allocatable :: surface_load_mesh                ! [Pa]   Total surface load applied at each mesh point
    real(dp), dimension(:  ), allocatable :: dHb_grid_partial                 ! [m]    Partial bedrock displacement on the grid
    real(dp), dimension(:,:), allocatable :: dHb_grid_tot                     ! [m]    Total bedrock displacement on the grid
    real(dp), dimension(:,:), allocatable :: relaxation_time_grid             ! [yr]   Bedrock relaxation time on the grid
    real(dp), dimension(:,:), allocatable :: dHb_eq_grid                      ! [m]    Bedrock displacement at equilibrium (computed from relative surface load)
    real(dp), dimension(:,:), allocatable :: dHb_dt_grid                      ! [m/yr] Bedrock deformation rate on the grid
    real(dp), dimension(:  ), allocatable :: dHb_dt_grid_partial              ! [m/yr] Partial bedrock deformation rate on the grid
    real(dp), dimension(:  ), allocatable :: dHb_dt_mesh                      ! [m/yr] Bedrock deformation rate mapped on the mesh
    real(dp), dimension(:,:), allocatable :: relative_surface_load_grid_tot   ! [Pa]   Relative surface load (difference from equilibrium) on the full grid
    real(dp), dimension(:  ), allocatable :: surface_load_GIAeq               ! [Pa]   Reference surface load used for GIA equilibrium computation
    integer                               :: flex_prof_rad                    ! radius of the flexoral profile
    real(dp), dimension(:,:), allocatable :: flex_prof_grid                   ! radius of the flexoral profile on GIA grid

  end type type_ELRA_model

end module GIA_model_types
