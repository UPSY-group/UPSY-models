module masks_mod
  !< Calculating masks (e.g. mask_ice, mask_shelf, etc.)

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mpi_distributed_memory, only: gather_to_all
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use ice_geometry_basics, only: is_floating
  use projections, only: oblique_sg_projection
  use plane_geometry, only: is_in_polygon
  use mesh_ROI_polygons
  use netcdf_io_main

  implicit none

  private

  public :: determine_masks, calc_mask_ROI, calc_mask_noice, calc_mask_noice_remove_Ellesmere, calc_mask_SGD

contains

  subroutine determine_masks( mesh, Hi, Hb, SL, &
    mask, mask_icefree_land, mask_icefree_ocean, mask_grounded_ice, mask_floating_ice, &
    mask_margin, mask_gl_fl, mask_gl_gr, mask_cf_gr, mask_cf_fl,mask_coastline)
    !< Determine the different masks

    ! In variables
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hi                      ! [m]       Ice thickness
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: Hb                      ! [m]       Bedrock elevation
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in   ) :: SL                      ! [m]       Water surface elevation

    ! Out variables
    logical,  dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_icefree_land       ! T: ice-free land , F: otherwise
    logical,  dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_icefree_ocean      ! T: ice-free ocean, F: otherwise
    logical,  dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_grounded_ice       ! T: grounded ice  , F: otherwise
    logical,  dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_floating_ice       ! T: floating ice  , F: otherwise
    logical,  dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_margin             ! T: ice next to ice-free, F: otherwise
    logical,  dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_gl_gr              ! T: grounded ice next to floating ice, F: otherwise
    logical,  dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_gl_fl              ! T: floating ice next to grounded ice, F: otherwise
    logical,  dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_cf_gr              ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
    logical,  dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_cf_fl              ! T: floating ice next to ice-free water (sea or lake), F: otherwise
    logical,  dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask_coastline          ! T: ice-free land next to ice-free ocean, F: otherwise
    integer,  dimension(mesh%vi1:mesh%vi2), intent(  out) :: mask

    ! Local variables:
    character(len=*), parameter :: routine_name = 'determine_masks'
    integer                     :: vi, ci, vj
    logical, dimension(mesh%nV) :: mask_icefree_land_tot
    logical, dimension(mesh%nV) :: mask_icefree_ocean_tot
    logical, dimension(mesh%nV) :: mask_grounded_ice_tot
    logical, dimension(mesh%nV) :: mask_floating_ice_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! === Basic masks ===
    ! ===================

    ! Initialise basic masks
    mask_icefree_land  = .false.
    mask_icefree_ocean = .false.
    mask_grounded_ice  = .false.
    mask_floating_ice  = .false.
    mask               = 0

    ! Calculate
    do vi = mesh%vi1, mesh%vi2

      if (is_floating( Hi( vi), Hb( vi), SL( vi))) then
        ! Ice thickness is below the floatation thickness; either floating ice, or ice-free ocean

        if (Hi( vi) > 0._dp) then
          ! Floating ice

          mask_floating_ice( vi) = .true.
          mask( vi) = C%type_floating_ice

        else
          ! Ice-free ocean

          mask_icefree_ocean( vi) = .true.
          mask( vi) = C%type_icefree_ocean

        end if

      else
        ! Ice thickness is above the floatation thickness; either grounded ice, or ice-free land

        if (Hi( vi) > 0._dp) then
          ! Grounded ice

          mask_grounded_ice( vi) = .true.
          mask( vi) = C%type_grounded_ice

        else
          ! Ice-free land

          mask_icefree_land( vi) = .true.
          mask( vi) = C%type_icefree_land

        end if

      end if

    end do

    ! === Transitional masks ===
    ! ==========================

    ! Gather basic masks to all processes
    call gather_to_all( mask_icefree_land , mask_icefree_land_tot )
    call gather_to_all( mask_icefree_ocean, mask_icefree_ocean_tot)
    call gather_to_all( mask_grounded_ice , mask_grounded_ice_tot )
    call gather_to_all( mask_floating_ice , mask_floating_ice_tot )

    ! Initialise transitional masks
    mask_margin    = .false.
    mask_gl_gr     = .false.
    mask_gl_fl     = .false.
    mask_cf_gr     = .false.
    mask_cf_fl     = .false.
    mask_coastline = .false.

    do vi = mesh%vi1, mesh%vi2

      ! Ice margin
      if (mask_grounded_ice_tot( vi) .OR. mask_floating_ice_tot( vi)) then
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          if (.not. (mask_grounded_ice_tot( vj) .OR. mask_floating_ice_tot( vj))) then
            mask_margin( vi) = .true.
            mask( vi) = C%type_margin
          end if
        end do
      end if

      ! Grounding line (grounded side)
      if (mask_grounded_ice_tot( vi)) then
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          if (mask_floating_ice_tot( vj)) then
            mask_gl_gr( vi) = .true.
            mask( vi) = C%type_groundingline_gr
          end if
        end do
      end if

      ! Grounding line (floating side)
      if (mask_floating_ice_tot( vi)) then
        do ci = 1, mesh%nC( vi)
          vj = mesh%C( vi,ci)
          if (mask_grounded_ice_tot( vj)) then
            mask_gl_fl( vi) = .true.
            mask( vi) = C%type_groundingline_fl
          end if
        end do
      end if

      ! Calving front (grounded)
      if (mask_grounded_ice_tot( vi)) then
        do ci = 1, mesh%nC(vi)
          vj = mesh%C( vi,ci)
          if (mask_icefree_ocean_tot( vj)) then
            mask_cf_gr( vi) = .true.
            mask( vi) = C%type_calvingfront_gr
          end if
        end do
      end if

      ! Calving front (floating)
      if (mask_floating_ice_tot( vi)) then
        do ci = 1, mesh%nC(vi)
          vj = mesh%C( vi,ci)
          if (mask_icefree_ocean_tot( vj)) then
            mask_cf_fl( vi) = .true.
            mask( vi) = C%type_calvingfront_fl
          end if
        end do
      end if

      ! Coastline
      if (mask_icefree_land_tot( vi)) then
        do ci = 1, mesh%nC(vi)
          vj = mesh%C( vi,ci)
          if (mask_icefree_ocean_tot( vj)) then
            mask_coastline( vi) = .true.
            mask( vi) = C%type_coastline
          end if
        end do
      end if

    end do ! do vi = mesh%vi1, mesh%vi2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine determine_masks

  subroutine calc_mask_ROI( mesh, ice, region_name)
    !< Calculate the ROI mask

    ! In/output variables:
    type(type_mesh),      intent(in   )           :: mesh
    type(type_ice_model), intent(inout)           :: ice
    character(len=3),     intent(in   )           :: region_name

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'calc_mask_ROI'
    character(len=256)                            :: all_names_ROI, name_ROI
    integer                                       :: vi, vj, ci, i, i_ROI
    real(dp), dimension(:,:  ), allocatable       :: poly_ROI
    real(dp), dimension(2)                        :: point

    ! Add routine to path
    call init_routine( routine_name)

    ! if no regions of interest are specified, do nothing
    if (C%choice_regions_of_interest == '') then
      call finalise_routine( routine_name)
      return
    end if

    all_names_ROI = C%choice_regions_of_interest
    i_ROI = 0

    do while (.true.)

      ! == Parse list of input ROIs
      ! ===========================

      ! Get the first region of interest from the list
      i = INDEX( all_names_ROI, '||')
      if (i == 0) then
        ! There is only one left in the list
        name_ROI = TRIM( all_names_ROI)
        all_names_ROI = ''
      else
        ! Get the first first one from the list and remove it
        name_ROI = all_names_ROI( 1:i-1)
        all_names_ROI = all_names_ROI( i+2:LEN_TRIM( all_names_ROI))
      end if

      ! == Check validity of requested ROIs
      ! ===================================

      ! Check if current region is indeed defined in the model
      select case (name_ROI)
        case ('')
          ! No region requested: don't need to do anything
          exit
        case ('PineIsland','Thwaites','Amery','RiiserLarsen','SipleCoast', 'LarsenC', &
              'TransMounts','DotsonCrosson', 'Franka_WAIS', 'Dotson_channel','Wilkes', &
              'Antarctic_Peninsula', 'Institute', &                                           ! Antarctica
              'Narsarsuaq','Nuuk','Jakobshavn','NGIS','Qaanaaq', &                            ! Greenland
              'Patagonia', &                                                                  ! Patagonia
              'CalvMIP_quarter')                                                              ! Idealised
          ! List of known regions of interest: these pass the test
        case default
          ! Region not found
          call crash('unknown region of interest "' // TRIM( name_ROI) // '"!')
      end select

      ! == Calculate ROIs
      ! =================

      ! Calculate the polygon describing the specified region of interest
      select case (region_name)
        case ('NAM')
          ! North america

          select case (name_ROI)
            case default
              ! Requested area not in this model domain; skip
              cycle
          end select

        case ('EAS')
          ! Eurasia

          select case (name_ROI)
            case default
              ! Requested area not in this model domain; skip
              cycle
          end select

        case ('GRL')
          ! Greenland

          select case (name_ROI)
            case ('Narsarsuaq')
              call calc_polygon_Narsarsuaq( poly_ROI)
            case ('Nuuk')
              call calc_polygon_Nuuk( poly_ROI)
            case ('Jakobshavn')
              call calc_polygon_Jakobshavn( poly_ROI)
            case ('NGIS')
              call calc_polygon_NGIS( poly_ROI)
            case ('Qaanaaq')
              call calc_polygon_Qaanaaq( poly_ROI)
            case default
              ! Requested area not in this model domain; skip
              cycle
          end select

        case ('ANT')

          select case (name_ROI)
            case ('PineIsland')
              call calc_polygon_Pine_Island_Glacier( poly_ROI)
            case ('Thwaites')
              call calc_polygon_Thwaites_Glacier( poly_ROI)
            case ('Amery')
              call calc_polygon_Amery_ice_shelf( poly_ROI)
            case ('RiiserLarsen')
              call calc_polygon_Riiser_Larsen_ice_shelf( poly_ROI)
            case ('SipleCoast')
              call calc_polygon_Siple_Coast( poly_ROI)
            case ('LarsenC')
              call calc_polygon_Larsen_ice_shelf( poly_ROI)
            case ('TransMounts')
              call calc_polygon_Transantarctic_Mountains( poly_ROI)
            case ('DotsonCrosson')
              call calc_polygon_DotsonCrosson_ice_shelf( poly_ROI)
            case ('Patagonia')
              call calc_polygon_Patagonia( poly_ROI)
            case ('CalvMIP_quarter')
              call calc_polygon_CalvMIP_quarter( poly_ROI)
            case ('Franka_WAIS')
              call calc_polygon_Franka_WAIS( poly_ROI)
            case ('Dotson_channel')
              call calc_polygon_Dotson_channel( poly_ROI)
            case ('Wilkes')
              call calc_polygon_Wilkes_basins( poly_ROI)
            case ('Antarctic_Peninsula')
              call calc_polygon_Antarctic_Peninsula( poly_ROI)
            case ('Institute')
              call calc_polygon_Institute_basin( poly_ROI)
            case default
              ! Requested area not in this model domain; skip
              cycle
          end select

        case default
          call crash('unknown region name "' // region_name // '"!')
      end select

      ! Check for each grid point whether it is located within the polygon of the ROI
      i_ROI = i_ROI+1
      do vi = mesh%vi1, mesh%vi2
        do ci = 1, mesh%nC(vi)
            vj = mesh%C( vi,ci)
            point = mesh%V( vj,:) ! Just to make sure it's in the right format
            if (is_in_polygon(poly_ROI, point)) then
              ice%mask_ROI(vi) = i_ROI
            end if
        end do
      end do ! do vi = mesh%vi1, mesh%vi2

      ! Clean up after yourself
      deallocate( poly_ROI)

    end do
    ice%nROI = i_ROI ! keep track of how many ROIs we actually have in the mask

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_mask_ROI

  subroutine calc_mask_noice( mesh, ice)
    !< Calculate the no-ice mask

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_mask_noice'
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise
    ! ==========

    ice%mask_noice = .false.

    ! domain-specific cases (mutually exclusive)
    ! ==========================================

    select case (C%choice_mask_noice)
    case default
      call crash('unknown choice_mask_noice "' // trim( C%choice_mask_noice) // '"')
      case ('none')
        ! Ice is (in principle) allowed everywhere

        ice%mask_noice = .false.

      case ('MISMIP_mod')
        ! Kill all ice when r > 900 km

        do vi = mesh%vi1, mesh%vi2
          if (NORM2( mesh%V( vi,:)) > 900E3_dp) then
            ice%mask_noice( vi) = .true.
          else
            ice%mask_noice( vi) = .false.
          end if
        end do

      case ('MISMIP+')
        ! Kill all ice when x > 640 km

        do vi = mesh%vi1, mesh%vi2
          if (mesh%V( vi,1) > 640E3_dp) then
            ice%mask_noice( vi) = .true.
          else
            ice%mask_noice( vi) = .false.
          end if
        end do

      case ('remove_Ellesmere')
        ! Prevent ice growth in the Ellesmere Island part of the Greenland domain

        call calc_mask_noice_remove_Ellesmere( mesh, ice%mask_noice)

      case ('Thule')
        ! Prevent ice growth in the Thule area
        do vi = mesh%vi1, mesh%vi2
          if (NORM2( mesh%V( vi,:)) > 750E3_dp) then
            ice%mask_noice( vi) = .true.
          else
            ice%mask_noice( vi) = .false.
          end if
        end do

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_mask_noice

  subroutine calc_mask_noice_remove_Ellesmere( mesh, mask_noice)
    !< Prevent ice growth in the Ellesmere Island part of the Greenland domain

    ! In- and output variables
    type(type_mesh),                        intent(in   ) :: mesh
    logical,  dimension(mesh%vi1:mesh%vi2), intent(inout) :: mask_noice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_mask_noice_remove_Ellesmere'
    integer                        :: vi
    real(dp), dimension(2)         :: pa_latlon, pb_latlon, pa, pb
    real(dp)                       :: xa, ya, xb, yb, yl_ab

    ! Add routine to path
    call init_routine( routine_name)

    ! The two endpoints in lat,lon
    pa_latlon = [76.74_dp, -74.79_dp]
    pb_latlon = [82.19_dp, -60.00_dp]

    ! The two endpoints in x,y
    call oblique_sg_projection( pa_latlon(2), pa_latlon(1), mesh%lambda_M, mesh%phi_M, mesh%beta_stereo, xa, ya)
    call oblique_sg_projection( pb_latlon(2), pb_latlon(1), mesh%lambda_M, mesh%phi_M, mesh%beta_stereo, xb, yb)

    pa = [xa,ya]
    pb = [xb,yb]

    do vi = mesh%vi1, mesh%vi2
      yl_ab = pa(2) + (mesh%V( vi,1) - pa(1)) * (pb(2)-pa(2)) / (pb(1)-pa(1))
      if (mesh%V( vi,2) > pa(2) .and. mesh%V( vi,2) > yl_ab .and. mesh%V( vi,1) < pb(1)) then
        mask_noice( vi) = .true.
      else
        mask_noice( vi) = .false.
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_mask_noice_remove_Ellesmere

  subroutine calc_mask_SGD( mesh, ice)
    !< Calculate the subglacial discharge mask (SGD)

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_mask_SGD'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_laddie_SGD)
      case ('none', 'read_transects')
        ! No SGD: don't need to do anything
      case ('idealised')
        call calc_mask_SGD_idealised(mesh, ice)
      case ('read_from_file')
        call calc_mask_SGD_from_file(mesh, ice)
      case default
        ! Choice not found
        call crash('unknown choice_laddie_SGD "' // TRIM( C%choice_laddie_SGD) // '"!')
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_mask_SGD

  subroutine calc_mask_SGD_idealised( mesh, ice)
    !< Calculate the subglacial discharge mask (SGD) - idealised cases for the MISMIPplus geometry

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_mask_SGD_idealised'
    integer                        :: vi
    real(dp)                       :: y_coord_channel

    ! Add routine to path
    call init_routine( routine_name)

    ! Set the y-coordinate of the subglacial channel based on idealised configuration.
    if (C%choice_laddie_SGD_idealised == 'MISMIPplus_PC') then
      y_coord_channel = 0._dp
    elseif (C%choice_laddie_SGD_idealised == 'MISMIPplus_PW') then
      y_coord_channel = 18000._dp
    elseif (C%choice_laddie_SGD_idealised == 'MISMIPplus_PE') then
      y_coord_channel = -18000._dp
    end if

    do vi = mesh%vi1, mesh%vi2
      ! Assign mask_SGD as true when vertice falls within a 5 km band centered on the channel.
      ! NOTE: This approach allows for some vertices along the transect at a given y-coordinate to be marked as false.
      !       Therefore, if the grounding line retreats to a position where no vertices are marked as true, no SGD will be applied.
      if (mesh%V( vi,2) < y_coord_channel + 2500._dp .and. mesh%V( vi,2) > y_coord_channel - 2500._dp) then
        ice%mask_SGD( vi) = .true.
      else
        ice%mask_SGD( vi) = .false.
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_mask_SGD_idealised

  subroutine calc_mask_SGD_from_file( mesh, ice)
    !< Calculate the subglacial discharge (SGD) mask by reading values from an external file.

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    type(type_ice_model), intent(inout) :: ice

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'calc_mask_SGD_from_file'
    integer                                :: vi
    REAL(dp), DIMENSION(mesh%vi1:mesh%vi2) :: temporary_mask_SGD

    ! Add routine to path
    call init_routine( routine_name)

    ! Read the SGD mask from a 2D input file. Values > 0 indicate presence of an SGD channel.
    call read_field_from_file_2D( C%filename_laddie_mask_SGD, 'mask_SGD', mesh, C%output_dir, temporary_mask_SGD)

    ! Assign mask_SGD as true where the read-in field is greater than zero.
    do vi = mesh%vi1, mesh%vi2
      if (temporary_mask_SGD( vi)>0) then
        ice%mask_SGD( vi) = .true.
      else
        ice%mask_SGD( vi) = .false.
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_mask_SGD_from_file

end module masks_mod
