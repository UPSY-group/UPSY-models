module SMB_reconstructed

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use grid_basic, only: type_grid
  use ice_model_types, only: type_ice_model
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use mesh_ROI_polygons, only: calc_polygon_Patagonia
  use plane_geometry, only: is_in_polygon
  use mesh_data_smoothing, only: smooth_Gaussian
  use SMB_basic, only: atype_SMB_model, type_SMB_model_context_allocate, &
    type_SMB_model_context_initialise, type_SMB_model_context_run, type_SMB_model_context_remap

  implicit none

  private

  public :: type_SMB_model_reconstructed

  type, extends(atype_SMB_model) :: type_SMB_model_reconstructed

    contains

      procedure, public :: allocate_SMB_model   => allocate_SMB_model_reconstructed
      procedure, public :: initialise_SMB_model => initialise_SMB_model_reconstructed
      procedure, public :: run_SMB_model        => run_SMB_model_reconstructed_abs
      procedure, public :: remap_SMB_model      => remap_SMB_model_reconstructed

      procedure, private :: run_SMB_model_reconstructed

  end type type_SMB_model_reconstructed

contains

  subroutine allocate_SMB_model_reconstructed( self, context)

    ! In/output variables:
    class(type_SMB_model_reconstructed),   intent(inout) :: self
    type(type_SMB_model_context_allocate), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SMB_model_reconstructed'

    ! Add routine to path
    call init_routine( routine_name)

    call self%set_name('SMB_model_reconstructed')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_reconstructed

  subroutine initialise_SMB_model_reconstructed( self, context)

    ! In/output variables:
    class(type_SMB_model_reconstructed),     intent(inout) :: self
    type(type_SMB_model_context_initialise), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_reconstructed'

    ! Add routine to path
    call init_routine( routine_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_reconstructed

  subroutine run_SMB_model_reconstructed_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_reconstructed), intent(inout) :: self
    type(type_SMB_model_context_run),    intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_reconstructed_abs'

    ! Add routine to path
    call init_routine( routine_name)

    call self%run_SMB_model_reconstructed( self%mesh, context%grid_smooth, context%ice)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_reconstructed_abs

  subroutine remap_SMB_model_reconstructed( self, context)

    ! In/output variables:
    class(type_SMB_model_reconstructed), intent(inout) :: self
    type(type_SMB_model_context_remap),  intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_SMB_model_reconstructed'

    ! Add routine to path
    call init_routine( routine_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_reconstructed

  subroutine run_SMB_model_reconstructed( self, mesh, grid_smooth, ice)

    ! In/output variables:
    class(type_SMB_model_reconstructed), intent(inout) :: self
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_grid),                     intent(in   ) :: grid_smooth
    type(type_ice_model),                intent(in   ) :: ice

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'run_SMB_model_reconstructed'
    integer                                 :: vi
    real(dp), dimension(:,:  ), allocatable :: poly_ROI             ! Polygon defining reconstructed area
    real(dp), dimension(2)                  :: p                    ! Coordinates of a vertex
    real(dp), dimension(mesh%vi1:mesh%vi2)  :: SMB_smoothed         ! Smoothed SMB field
    real(dp)                                :: w_smooth             ! Weight of the smoothed SMB field
    real(dp), parameter                     :: r_smooth =  2.E4_dp  ! Radius used to smooth the SMB field
    real(dp), parameter                     :: Hs_ela   =  500._dp  ! Equilibrium line altitud: SMB becomes positive here
    real(dp), parameter                     :: Hs_tla   =  1500._dp ! Transitional line altitud: SMB reaches maximum here
    real(dp), parameter                     :: Hs_dla   =  2500._dp ! Desertification line altitude: SMB becomes zero here
    real(dp), parameter                     :: SMB_max  =  2._dp    ! Maximum SMB value allowed
    real(dp), parameter                     :: SMB_min  = -10._dp   ! Minimum SMB value allowed

    ! Add routine to path
    call init_routine( routine_name)

    if (.not. C%choice_regions_of_interest == 'Patagonia') THEN
      call crash('reconstructed SMB method only implemented for C%choice_regions_of_interest == Patagonia')
    end if

    ! Compute polygon for reconstruction
    call calc_polygon_Patagonia( poly_ROI)

    do vi = mesh%vi1, mesh%vi2

      ! Get x and y coordinates of this vertex
      p = mesh%V( vi,:)

      ! Check if point lies within our reconstruction polygon
      if (is_in_polygon(poly_ROI, p)) then
        ! If yes, check whether point lies above or below estimated transitional line altitude
        if (ice%Hs( vi) <= Hs_tla) then
          ! If below, SMB goes from 0 at the ELA to its estimated maximum at the TLA
          self%SMB( vi) = SMB_max * MAX( 0._dp, min( 1._dp, (ice%Hs( vi) - Hs_ela)/(Hs_tla - Hs_ela)))
        else
          ! If above, SMB goes from estimated maximum at the TLA to 0 at the DLA
          self%SMB( vi) = SMB_max * (1._dp - max( 0._dp, min( 1._dp, (ice%Hs( vi) - Hs_tla)/(Hs_dla - Hs_tla))))
        end if
      else
        ! If vertex lies outside of the reconstructed polygon, assume a negative
        ! SMB that counters (and then some) the flux convergence there
        self%SMB( vi) = min( 0._dp, max( SMB_min, ice%divQ( vi) - .5_dp))
      end if

    end do

    ! Smooth the reconstructed field
    SMB_smoothed = self%SMB( mesh%vi1:mesh%vi2)
    call smooth_Gaussian( mesh, grid_smooth, C%output_dir, SMB_smoothed, r_smooth)

    ! Only apply the smoothed field inside the reconstructed area
    ! to reduce the power of positive SMB there
    do vi = mesh%vi1, mesh%vi2
      ! Our vextex coordinates
      p = mesh%V( vi,:)
      ! Check if point lies inside polygon
      if (is_in_polygon(poly_ROI, p)) then
        ! Compute a weight based on Hs: the higher, the less smoothing
        w_smooth = max( 0._dp, min( 1._dp, ice%Hs( vi) / Hs_dla))
        ! Apply weighed smoothing
        self%SMB( vi) = w_smooth * self%SMB( vi) + (1._dp - w_smooth) * SMB_smoothed( vi)
      end if
    end do

    ! Smooth the field once more
    SMB_smoothed = self%SMB( mesh%vi1:mesh%vi2)
    call smooth_Gaussian( mesh, grid_smooth, C%output_dir, SMB_smoothed, r_smooth)

    ! Apply this second smoothing only outside of the reconstructed area
    ! to conserve the power of negative SMB there
    do vi = mesh%vi1, mesh%vi2
      p = mesh%V( vi,:)
      if (.not. is_in_polygon( poly_ROI, p)) self%SMB( vi) = SMB_smoothed( vi)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_reconstructed

end module SMB_reconstructed