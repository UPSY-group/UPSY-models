module subgrid_grounded_fractions_main
  !< Routines for calculating sub-grid grounded fractions

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use mpi_distributed_memory, only: gather_to_all
  use subgrid_grounded_fractions_bedrock_CDF
  use subgrid_grounded_fractions_bilin_TAF
   use ice_geometry_basics, only: thickness_above_floatation

  implicit none

  private

  public :: calc_grounded_fractions

contains

  subroutine calc_grounded_fractions( mesh, Hi, Hb, SL, dHb, fraction_gr, fraction_gr_b, mask_floating_ice, bedrock_cdf, bedrock_cdf_b)
    !< Calculate the sub-grid grounded-area fractions

    ! In- and output variables
    type(type_mesh),                                                     intent(in   ) :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2),                              intent(in   ) :: Hi
    real(dp), dimension(mesh%vi1:mesh%vi2),                              intent(in   ) :: Hb
    real(dp), dimension(mesh%vi1:mesh%vi2),                              intent(in   ) :: SL
    real(dp), dimension(mesh%vi1:mesh%vi2),                              intent(in   ) :: dHb
    logical,  dimension(mesh%vi1:mesh%vi2),                              intent(in   ) :: mask_floating_ice
    real(dp), dimension(mesh%vi1:mesh%vi2, C%subgrid_bedrock_cdf_nbins), intent(in   ) :: bedrock_cdf
    real(dp), dimension(mesh%ti1:mesh%ti2, C%subgrid_bedrock_cdf_nbins), intent(in   ) :: bedrock_cdf_b
    real(dp), dimension(mesh%vi1:mesh%vi2),                              intent(  out) :: fraction_gr
    real(dp), dimension(mesh%ti1:mesh%ti2),                              intent(  out) :: fraction_gr_b

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_grounded_fractions'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: TAF
    real(dp), dimension(mesh%vi1:mesh%vi2) :: fraction_gr_TAF_a
    real(dp), dimension(mesh%vi1:mesh%vi2) :: fraction_gr_CDF_a
    real(dp), dimension(mesh%ti1:mesh%ti2) :: fraction_gr_TAF_b
    real(dp), dimension(mesh%ti1:mesh%ti2) :: fraction_gr_CDF_b
    logical,  dimension(mesh%nV)           :: mask_floating_ice_tot
    integer                                :: ti, via, vib, vic, vi


    ! Add routine to path
    call init_routine( routine_name)


    do vi = mesh%vi1, mesh%vi2
      TAF = thickness_above_floatation( Hi( vi), Hb( vi), SL( vi))
    end do

    ! Use the specified way of calculating sub-grid grounded fractions
    select case (C%choice_subgrid_grounded_fraction)
    case default
      call crash('unknown choice_subgrid_grounded_fraction "' // &
        trim( C%choice_subgrid_grounded_fraction) // '"')
    case('bilin_interp_TAF')
      ! Bilinearly interpolate the thickness above floatation to calculate the grounded fractions

      call calc_grounded_fractions_bilin_interp_TAF_a( mesh, TAF, fraction_gr_TAF_a)
      call calc_grounded_fractions_bilin_interp_TAF_b( mesh, TAF, fraction_gr_TAF_b)

      fraction_gr   = fraction_gr_TAF_a
      fraction_gr_b = fraction_gr_TAF_b

    case ('bedrock_CDF')
      ! Use the sub-grid bedrock cumulative density functions to calculate the grounded fractions

      call calc_grounded_fractions_bedrock_CDF_a( mesh, Hi, SL, dHb, bedrock_cdf, fraction_gr_CDF_a)
      call calc_grounded_fractions_bedrock_CDF_b( mesh, Hi, SL, dHb, TAF, bedrock_cdf_b, fraction_gr_CDF_b)

      fraction_gr   = fraction_gr_CDF_a
      fraction_gr_b = fraction_gr_CDF_b

    case ('bilin_interp_TAF+bedrock_CDF')
      ! Use the TAF method at the grounding line, and the CDF method inland

      call calc_grounded_fractions_bilin_interp_TAF_a( mesh, TAF, fraction_gr_TAF_a)
      call calc_grounded_fractions_bilin_interp_TAF_b( mesh, TAF, fraction_gr_TAF_b)

      call calc_grounded_fractions_bedrock_CDF_a( mesh, Hi, SL, dHb, bedrock_cdf, fraction_gr_CDF_a)
      call calc_grounded_fractions_bedrock_CDF_b( mesh, Hi, SL, dHb,TAF, bedrock_cdf_b, fraction_gr_CDF_b)

      ! Gather global floating ice mask
      call gather_to_all( mask_floating_ice, mask_floating_ice_tot)

      ! a-grid (vertices): take the smallest value (used for basal melt?)
      fraction_gr = min( fraction_gr_TAF_a, fraction_gr_CDF_a)

      ! b-grid (triangles): take CDF inland, TAF at grounding line (used for basal friction)
      do ti = mesh%ti1, mesh%ti2

        ! The three vertices spanning triangle ti
        via = mesh%Tri( ti,1)
        vib = mesh%Tri( ti,2)
        vic = mesh%Tri( ti,3)

        if (mask_floating_ice_tot( via) .OR. mask_floating_ice_tot( vib) .OR. mask_floating_ice_tot( vic)) then
          ! At least one corner of this triangle is afloat; grounding line
          fraction_gr_b( ti) = fraction_gr_TAF_b( ti)
        else
          ! All three corners of the triangle are grounded: inland
          fraction_gr_b( ti) = fraction_gr_CDF_b( ti)
        end if

      end do

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grounded_fractions

end module subgrid_grounded_fractions_main