module GIA_main

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: crash, init_routine, finalise_routine
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use ice_model_data, only: atype_ice_model_data
  use GIA_model_types, only: type_GIA_model, type_ELRA_model
  use region_types, only: type_model_region
  use grid_basic, only: setup_square_grid
  use reallocate_mod, only: reallocate_bounds
  use reference_geometry_types, only: type_reference_geometry
  use GIA_ELRA, only: run_ELRA_model, calculate_ELRA_bedrock_deformation_rate, initialise_ELRA_model, remap_ELRA_model
  use checksum_mod, only: checksum

  implicit none

  private

  public :: initialise_GIA_model, run_GIA_model, create_restart_file_GIA_model, &
    write_to_restart_file_GIA_model, remap_GIA_model

contains

  subroutine run_GIA_model( region)
    ! Calculate bedrock deformation at the desired time, and update
    ! predicted deformation if necessary

    ! In/output variables:
    type(type_model_region), intent(inout) :: region

    ! Local variables:
    character(len=*), parameter :: routine_name = 'run_GIA_model'
    real(dp)                    :: wt_prev, wt_next
    integer                     :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! If the desired time is beyond the time of the next modelled bedrock deformation,
    ! run the GIA model to calculate a new next modelled bedrock deformation.
    ! ============================================================================

    if (region%time == region%GIA%t_next) then
      ! Need to calculate new predicted bedrock deformation

      ! Store previous modelled bedrock deformation
      region%GIA%dHb_prev = region%GIA%dHb_next
      region%GIA%t_prev = region%GIA%t_next
      region%GIA%t_next = region%GIA%t_prev + C%dt_GIA

      ! Run the GIA model to calculate a new next modelled bedrock deformation
      select case (C%choice_GIA_model)
      case default
        CALL crash('unknown choice_GIA_model "' // trim( C%choice_GIA_model) // '"!')
      case ('none')
        ! No need to do anything
      case ('ELRA')
        call run_ELRA_model( region)
      end select

    elseif (region%time > region%GIA%t_next) then
      ! This should not be possible
      call crash('overshot the GIA time step')
    else
      ! We're within the current GIA prediction window
    end if

    ! Interpolate between previous and next modelled bedrock deformation
    ! to find the bedrock deformation and elevation at the desired time
    ! =================================================================

    ! Calculate time interpolation weights
    wt_prev = (region%GIA%t_next - region%time) / (region%GIA%t_next - region%GIA%t_prev)
    wt_next = 1._dp - wt_prev

    ! Interpolate modelled bedrock deformation to desired time
    do vi = region%mesh%vi1, region%mesh%vi2
      region%ice%dHb( vi) = wt_prev * region%GIA%dHb_prev( vi) + wt_next * region%GIA%dHb_next( vi)
    end do

    ! Calculate all other GIA quantities
    ! ==================================

    ! DO vi = region%mesh%vi1, region%mesh%vi2
    !   region%ice%dHb_dt( vi) = (region%GIA%dHb_next( vi) - region%GIA%dHb_prev( vi)) / C%dt_GIA
    ! END DO

    do vi = region%mesh%vi1, region%mesh%vi2
  	  region%ice%geom%Hb( vi) = region%refgeo_GIAeq%Hb( vi) + region%ice%dHb( vi)
	  end do

    call checksum( region%mesh%pai_V  , region%GIA%relative_surface_load_mesh, 'region%GIA%relative_surface_load_mesh')
    ! call checksum( region%GIA%grid%pai, region%GIA%relative_surface_load_grid, 'region%GIA%relative_surface_load_grid')
    call checksum( region%mesh%pai_V  , region%GIA%dHb_prev                  , 'region%GIA%dHb_prev')
    call checksum( region%mesh%pai_V  , region%GIA%dHb_next                  , 'region%GIA%dHb_next')
    call checksum( region%mesh%pai_V  , region%ice%dHb                       , 'region%ice%dHb')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_GIA_model

  subroutine initialise_GIA_model( mesh, GIA, region_name, refgeo_GIAeq, ELRA)

    ! In- and output variables
    type(type_mesh),               intent(in   ) :: mesh
    type(type_GIA_model),          intent(  out) :: GIA
    character(len=*),              intent(in   ) :: region_name
    type(type_reference_geometry), intent(in   ) :: refgeo_GIAeq
    type(type_ELRA_model),         intent(  out) :: ELRA

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'initialise_GIA_model'
    character(len=:), allocatable :: grid_name

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%primary)  write(*,"(A)") '   Initialising GIA model...'

    ! Create the square grid for the GIA model
    grid_name = 'square_grid_GIA_' // region_name
    call setup_square_grid( grid_name, mesh%xmin, mesh%xmax, mesh%ymin, mesh%ymax, &
       C%dx_GIA, GIA%grid, &
       lambda_M = mesh%lambda_M, phi_M = mesh%phi_M, beta_stereo = mesh%beta_stereo)

    ! Allocate memory for main variables

    ! In: relative surface load
    allocate( GIA%relative_surface_load_mesh( mesh%vi1:mesh%vi2)      , source = 0._dp)
    allocate( GIA%relative_surface_load_grid( GIA%grid%n1:GIA%grid%n2), source = 0._dp)

    ! Out: modelled bedrock deformation
    allocate( GIA%dHb_prev( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( GIA%dHb_next( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Model states for GIA model
    GIA%t_prev   = C%start_time_of_run
    GIA%t_next   = C%start_time_of_run

    ! Determine which GIA model to initialise
    select case (C%choice_GIA_model)
    case default
      call crash('unknown choice_GIA_model "' // trim( C%choice_GIA_model) // '"')
    case ('none')
      ! No need to do anything
    case ('ELRA')
      call initialise_ELRA_model( mesh, GIA%grid, ELRA, refgeo_GIAeq)
    end select

    call checksum( mesh%pai_V  , GIA%relative_surface_load_mesh, 'GIA%relative_surface_load_mesh')
    call checksum( GIA%grid%pai, GIA%relative_surface_load_grid, 'GIA%relative_surface_load_grid')
    call checksum( mesh%pai_V  , GIA%dHb_prev                  , 'GIA%dHb_prev')
    call checksum( mesh%pai_V  , GIA%dHb_next                  , 'GIA%dHb_next')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_GIA_model

  subroutine write_to_restart_file_GIA_model( mesh, GIA, region_name, time)

    ! In/output variables:
    type(type_mesh),      intent(in) :: mesh
    type(type_GIA_model), intent(in) :: GIA
    character(len=*),     intent(in) :: region_name
    real(dp),             intent(in) :: time

    ! Local variables:
    character(len=*), parameter :: routine_name = 'write_to_restart_file_GIA_model'

    ! Add routine to path
    call init_routine( routine_name)

    ! Write to the restart file of the chosen GIA model
    select case (C%choice_GIA_model)
    case default
      call crash('unknown choice_GIA_model "' // trim( C%choice_GIA_model) // '"')
    case ('none')
      ! No need to do anything
    case ('ELRA')
      ! No need to do anything
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file_GIA_model

  subroutine create_restart_file_GIA_model( mesh, GIA, region_name)

    ! In/output variables:
    type(type_mesh),      intent(in) :: mesh
    type(type_GIA_model), intent(in) :: GIA
    character(len=*),     intent(in) :: region_name

    ! Local variables:
    character(len=*), parameter :: routine_name = 'create_restart_file_GIA_model'

    ! Add routine to path
    call init_routine( routine_name)

    ! Create the restart file of the chosen GIA model
    select case (C%choice_GIA_model)
    case default
      call crash('unknown choice_GIA_model "' // trim( C%choice_GIA_model) // '"')
    case ('none')
      ! No need to do anything
    case ('ELRA')
      ! No need to do anything
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_GIA_model

  subroutine remap_GIA_model( mesh_old, mesh_new, GIA, refgeo_GIAeq, ELRA)

    ! In- and output variables
    type(type_mesh),               intent(in   ) :: mesh_old
    type(type_mesh),               intent(in   ) :: mesh_new
    type(type_GIA_model),          intent(inout) :: GIA
    type(type_reference_geometry), intent(in   ) :: refgeo_GIAeq
    type(type_ELRA_model),         intent(inout) :: ELRA

    ! Local variables:
    character(len=*), parameter :: routine_name = 'remap_GIA_model'

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%primary)  WRITE(*,"(A)") '    Remapping GIA model data to the new mesh...'

    ! Reallocate memory for main variables

    ! In: relative surface load
    call reallocate_bounds( GIA%relative_surface_load_mesh, mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( GIA%dHb_prev                  , mesh_new%vi1, mesh_new%vi2)
    call reallocate_bounds( GIA%dHb_next                  , mesh_new%vi1, mesh_new%vi2)

    ! Determine which GIA model to remap
    select case (C%choice_GIA_model)
    case default
      call crash('unknown choice_GIA_model "' // TRIM( C%choice_GIA_model) // '"')
    case ('none')
      ! No need to do anything
    case ('ELRA')
      call remap_ELRA_model( mesh_old, mesh_new, ELRA, refgeo_GIAeq, GIA%grid)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_GIA_model

end module GIA_main
