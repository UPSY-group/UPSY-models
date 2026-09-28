module SMB_ITM_v2

  use mpi_basic, only: par
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use SMB_model_basic, only: atype_SMB_model
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use mpi_f08, only: MPI_WIN
  use ice_model_data, only: atype_ice_model_data
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use climate_model_types, only: type_climate_model, type_climate_model_snapshot
  use netcdf_io_main, only: read_field_from_file_2D, read_field_from_file_2D_monthly
  use parameters, only: freshwater_density, ice_density, T0, L_fusion, sec_per_year, R_gas, grav
  use climate_model_utilities, only: get_insolation_at_time
  use climate_realistic, only: initialise_insolation_forcing
  use reference_geometry_types, only: type_reference_geometry
  use grid_types, only: type_grid

  implicit none

  private

  public :: type_SMB_model_ITM_v2

! ===== Types =====
! =================

  type, extends(atype_SMB_model) :: type_SMB_model_ITM_v2
    !< Variables and functions that are specific to the ITM_v2 SMB model

      ! Main data fields
      real(dp), dimension(:,:), contiguous, pointer :: FirnAirContent   => null() !< [m] firn air content
      real(dp), dimension(:,:), contiguous, pointer :: Rainfall         => null() !< Monthly rainfall (m)
      real(dp), dimension(:,:), contiguous, pointer :: Snowfall         => null() !< Monthly snowfall (m)
      real(dp), dimension(:,:), contiguous, pointer :: SurfaceMelt      => null() !< Monthly melt (m)
      real(dp), dimension(:,:), contiguous, pointer :: Refreezing       => null() !< Monthly refreezing (m)
      real(dp), dimension(:,:), contiguous, pointer :: Runoff           => null() !< Monthly runoff (m)
      real(dp), dimension(:,:), contiguous, pointer :: Albedo           => null() !< Monthly albedo
      real(dp), dimension(:,:), contiguous, pointer :: SMB_monthly      => null() !< [m] Monthly SMB
      real(dp), dimension(:,:), contiguous, pointer :: Sublimation      => null() !< [m.w.e.] Monthly sublimation
      type(MPI_WIN) :: wFirnAirContent, wRainfall
      type(MPI_WIN) :: wSnowfall, wSurfaceMelt, wRefreezing
      type(MPI_WIN) :: wRunoff, wAlbedo, wSMB_monthly, wSublimation

      ! Tuning parameters for the ITM_v2 SMB model (different for each region, set from config)
      real(dp)  :: C_refr

      ! Ideally these parameters should not be region-dependent?
      real(dp)  :: albedo_water
      real(dp)  :: albedo_soil
      real(dp)  :: albedo_ice
      real(dp)  :: albedo_snow

    contains

      procedure, public :: allocate_SMB_model   => SMB_model_ITM_v2_allocate
      procedure, public :: deallocate_SMB_model => SMB_model_ITM_v2_deallocate
      procedure, public :: initialise_SMB_model => SMB_model_ITM_v2_initialise
      procedure, public :: run_SMB_model        => SMB_model_ITM_v2_run
      procedure, public :: remap_SMB_model      => SMB_model_ITM_v2_remap

      procedure, public :: get_SMB_model_name

      procedure, private :: initialise_ITM_v2_firn_from_file

  end type type_SMB_model_ITM_v2

contains

  subroutine SMB_model_ITM_v2_allocate( self)

    ! In/output variables:
    class(type_SMB_model_ITM_v2), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_ITM_v2_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is specific to the ITM_v2 SMB model

    call self%create_field( self%FirnAirContent, self%wFirnAirContent, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'FirnAirContent', &
      long_name = 'Air content in the firn layer', &
      units     = 'm')

    call self%create_field( self%Rainfall, self%wRainfall, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Rainfall', &
      long_name = 'Monthly rainfall', &
      units     = 'm')

    call self%create_field( self%Snowfall, self%wSnowfall, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Snowfall', &
      long_name = 'Monthly snowfall', &
      units     = 'm')

    call self%create_field( self%SurfaceMelt, self%wSurfaceMelt, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'SurfaceMelt', &
      long_name = 'Monthly melt', &
      units     = 'm')

    call self%create_field( self%Refreezing, self%wRefreezing, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Refreezing', &
      long_name = 'Monthly refreezing', &
      units     = 'm')

    call self%create_field( self%Runoff, self%wRunoff, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Runoff', &
      long_name = 'Monthly runoff', &
      units     = 'm')

    call self%create_field( self%Albedo, self%wAlbedo, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Albedo', &
      long_name = 'Monthly albedo', &
      units     = '-')

    call self%create_field( self%SMB_monthly, self%wSMB_monthly, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'SMB_monthly', &
      long_name = 'Monthly SMB', &
      units     = '-')

    call self%create_field( self%Sublimation, self%wSublimation, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Sublimation', &
      long_name = 'Monthly Sublimation', &
      units     = 'm.w.e.')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine SMB_model_ITM_v2_allocate

  subroutine SMB_model_ITM_v2_deallocate( self)

    ! In/output variables:
    class(type_SMB_model_ITM_v2), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_ITM_v2_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to SMB model ITM_v2

    nullify( self%FirnAirContent)
    nullify( self%Rainfall)
    nullify( self%Snowfall)
    nullify( self%SurfaceMelt)
    nullify( self%Refreezing)
    nullify( self%Runoff)
    nullify( self%Albedo)
    nullify( self%SMB_monthly)
    nullify( self%Sublimation)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine SMB_model_ITM_v2_deallocate

  subroutine SMB_model_ITM_v2_initialise( self, geom, refgeo_init, refgeo_PD)

    ! In/output variables
    class(type_SMB_model_ITM_v2),       intent(inout) :: self
    class(atype_ice_geometry_model_data), intent(in   ) :: geom
    type(type_reference_geometry),        intent(in   ) :: refgeo_init
    type(type_reference_geometry),        intent(in   ) :: refgeo_PD

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_ITM_v2_initialise'
    integer                     :: vi
    character(:), allocatable   :: choice_SMB_IMAUITM_init_firn

    ! Add routine to path
    call init_routine( routine_name)

    ! Initialise all the stuff that is specific to SMB model ITM_v2

    ! Determine which constants to use for this region
    select case (self%region_name())
    case default
      call crash('unknown self%region_name "' // self%region_name() // '"')
    case ('NAM')
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_NAM
      self%C_refr               = C%SMB_IMAUITM_C_refr_NAM
    case ('EAS')
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_EAS
      self%C_refr               = C%SMB_IMAUITM_C_refr_EAS
    case ('GRL')
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_GRL
      self%C_refr               = C%SMB_IMAUITM_C_refr_GRL
    case ('ANT')
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_ANT
      self%C_refr               = C%SMB_IMAUITM_C_refr_ANT
    end select

    ! Initialising albedo values
    self%albedo_water = C%SMB_IMAUITM_albedo_water
    self%albedo_soil  = C%SMB_IMAUITM_albedo_soil
    self%albedo_ice   = C%SMB_IMAUITM_albedo_ice
    self%albedo_snow  = C%SMB_IMAUITM_albedo_snow

    ! Initialise the firn layer
    select case (choice_SMB_IMAUITM_init_firn)
    case default
      call crash('unknown choice_SMB_IMAUITM_init_firn "' // trim( choice_SMB_IMAUITM_init_firn) // '"')
    case ('uniform')
      ! Initialise with a uniform firn layer over the ice sheet

      do vi = self%mesh%vi1, self%mesh%vi2
        self%FirnAirContent  ( vi,:) = C%SMB_ITM_initial_firn_air_content
        self%SurfaceMelt     ( vi,:) = 0._dp
        self%Snowfall        ( vi,:) = 0._dp
      end do

    case ('read_from_file')
      ! Initialise with the firn layer of a previous run
      call self%initialise_ITM_v2_firn_from_file( self%mesh, self%region_name())
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine SMB_model_ITM_v2_initialise

  subroutine initialise_ITM_v2_firn_from_file( self, mesh, region_name)
    !< Initialise the firn depth and meltpreviousyear data from a NetCDF file

    ! In/output variables
    class(type_SMB_model_ITM_v2), intent(inout) :: self
    type(type_mesh),                intent(in   ) :: mesh
    character(len=3),               intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_ITM_v2_firn_from_file'
    character(:), allocatable      :: filename_restart_firn
    real(dp)                       :: timeframe_restart_firn

    ! Add routine to path
    call init_routine( routine_name)

    ! Assume that SMB and geometry are read from the same restart file
    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name) // '"')
    case ('NAM')
      filename_restart_firn  = C%filename_firn_IMAUITM_NAM
      timeframe_restart_firn = C%timeframe_restart_firn_IMAUITM_NAM
    case ('EAS')
      filename_restart_firn  = C%filename_firn_IMAUITM_EAS
      timeframe_restart_firn = C%timeframe_restart_firn_IMAUITM_EAS
    case ('GRL')
      filename_restart_firn  = C%filename_firn_IMAUITM_GRL
      timeframe_restart_firn = C%timeframe_restart_firn_IMAUITM_GRL
    case ('ANT')
      filename_restart_firn  = C%filename_firn_IMAUITM_ANT
      timeframe_restart_firn = C%timeframe_restart_firn_IMAUITM_ANT
    end select

     ! Print to terminal
    if (par%primary)  write(*,"(A)") '   Initialising SMB-model firn layer from file "' &
      // UPSY%stru%colour_string( trim( filename_restart_firn),'light blue') // '"...'

    ! Read firn layer from then
    if (timeframe_restart_firn == 1E9_dp) THEN
      ! Assume the file has no time dimension
      call read_field_from_file_2D_monthly( filename_restart_firn, 'FirnAirContent', mesh, C%output_dir, self%FirnAirContent)
      call read_field_from_file_2D_monthly( filename_restart_firn, 'SurfaceMelt', mesh, C%output_dir, self%SurfaceMelt)
      call read_field_from_file_2D_monthly( filename_restart_firn, 'Snowfall', mesh, C%output_dir, self%Snowfall)
    else
      ! Assume the file has a time dimension, and read the specified timeframe
      call read_field_from_file_2D_monthly( filename_restart_firn, 'FirnAirContent', mesh, C%output_dir, self%FirnAirContent, time_to_read = timeframe_restart_firn)
      call read_field_from_file_2D_monthly( filename_restart_firn, 'SurfaceMelt', mesh, C%output_dir, self%SurfaceMelt, time_to_read = timeframe_restart_firn)
      call read_field_from_file_2D_monthly( filename_restart_firn, 'Snowfall', mesh, C%output_dir, self%Snowfall, time_to_read = timeframe_restart_firn)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ITM_v2_firn_from_file

  subroutine SMB_model_ITM_v2_run( self, time, ice, geom, climate, grid_smooth)

    ! NOTE: all the SMB components are in meters of water equivalent;
    !       the end result (SMB_monthly and SMB) are in meters of ice equivalent.

    ! In/output variables:
    class(type_SMB_model_ITM_v2),         intent(inout) :: self
    real(dp),                             intent(in   ) :: time
    class(atype_ice_model_data),          intent(in   ) :: ice
    class(atype_ice_geometry_model_data), intent(in   ) :: geom
    type(type_climate_model),             intent(inout) :: climate
    type(type_grid),                      intent(in   ) :: grid_smooth

    ! Local variables:
    character(len=*), parameter       :: routine_name = 'SMB_model_ITM_v2_run'
    integer                           :: vi
    integer                           :: m, mprev, i
    real(dp)                          :: snowfrac, surface_snow_density, temp_exponent
    real(dp)                          :: timeframe_init_insolation
    type(type_climate_model_snapshot) :: snapshot_dummy
    real(dp)                          :: Ec = 60000._dp ! [J/mol] Creep activation energy for densification
    real(dp)                          :: Eg = 42400._dp ! [J/mol] Grain growth activation energy
    real(dp)                          :: fac_temp
    real(dp)                          :: fac_scale_albedo = 0.5_dp ! [m] Exponential decay scale of albedo with fac
    real(dp)                          :: melt_scale_albedo = .30_dp ! [month m^-1] Linear scaling with melt
    real(dp)                          :: wind_speed     ! [m/s] Absolute wind speed
    real(dp)                          :: wind_sens      ! [?]
    real(dp)                          :: firn_temp      ! [K] Temperature in firn layer
    real(dp)                          :: albedo_ref
    real(dp)                          :: Qsurf          ! [W m^-2] Surface insolation accounting for cloudiness

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Run all the stuff that is specific to SMB model idealised

    do vi = self%mesh%vi1, self%mesh%vi2

      if (geom%mask_icefree_ocean( vi)) then
        ! Set everything to zero for ocean. No SMB allowed here.
        ! NOTE for advancing calving fronts, advection or extrapolation
        ! of the FirnAirContent should be added.
        self%SMB( vi) = 0._dp
        do m = 1, 12
          self%Albedo( vi, m) = self%albedo_water
          self%SurfaceMelt( vi, m) = 0._dp
          self%Snowfall( vi, m) = 0._dp
          self%Rainfall( vi, m) = 0._dp
          self%FirnAirContent( vi, m) = 0._dp
          self%Refreezing( vi, m) = 0._dp
          self%Runoff( vi, m) = 0._dp
          self%SMB_monthly( vi, m) = 0._dp
        end do

      else
        ! Potentially ice-covered cells, so compute ITM

        ! Compute monthly albedo and melt
        do m = 1, 12  ! Month loop

          ! Define the previous month index
          mprev = m - 1
          if (mprev==0) mprev = 12

          ! Determine the snow fraction based on an empirical fit to RACMO2.4p1 data
          ! The usage of a tanh-curve ensures convergence to 1 for low temperatures
          ! and 0 for high temperatures
          snowfrac = 0.5_dp * (1- tanh((climate%T2m( vi, m) - 274.46_dp) / 2.4736_dp))

          ! Exctract snowfall and rainfall from snowfraction and total precipitation
          self%Snowfall( vi, m) = climate%Precip( vi, m) *          snowfrac
          self%Rainfall( vi, m) = climate%Precip( vi, m) * (1._dp - snowfrac)

          if (geom%Hi( vi) > 0._dp) then

            ! Determine albedo based on melt in the previous month and snowfall this month
            ! Tuned to RACMO2.4p1 21st century ssp370 (Hofsteenge et al., 2026)
            albedo_ref = max(self%albedo_ice, self%albedo_snow - .91_dp * self%SurfaceMelt( vi, mprev))
            self%Albedo( vi,m) = .89_dp - (.89_dp - albedo_ref) * exp(-self%Snowfall( vi, m) / 30.9_dp)

            ! Determine ablation as a function of surface temperature 
            ! and albedo/insolation according following Bintanja et al. (2002)
            ! Retuned to RACMO2.4p1 data
            Qsurf = climate%Q_TOA( vi, m) * (0.72_dp + 3.92E-5_dp * geom%Hs( vi)) 

            self%SurfaceMelt( vi,m) = &
                C%SMB_ITM_C_melt_temp * max(0._dp, (climate%T2m( vi,m) - C%SMB_ITM_C_trans_temp))**2 &
                + C%SMB_ITM_C_melt_insol * max(0._dp, (1.0_dp - self%Albedo( vi,m)) * Qsurf &
                  * (climate%T2m( vi,m) - C%SMB_ITM_C_trans_temp))
          else
            ! Ice free land
            self%Albedo( vi, m) = self%albedo_soil
            self%SurfaceMelt( vi, m) = 0._dp
          end if

          ! Compute refreezing as the minimum value of 1) available liquid water
          ! and 2) available firn air content (= tipping bucket method?)

          if (geom%Hi( vi) > 0._dp) then
            self%Refreezing( vi, m) = min( &
              self%Rainfall( vi, m) + self%SurfaceMelt( vi, m), &
              self%FirnAirContent( vi, m))
          else
            ! Ice free land
            self%Refreezing( vi, m) = 0._dp
          end if

          ! Determine sublimation based on a fit to RACMO 2.4p1
          wind_sens = min( &
            max(0._dp, 4.0E-4_dp*(climate%T2m( vi, m) - 224._dp)), &
            -1.0E-2_dp * (climate%T2m( vi, m) - T0))

          wind_speed = hypot(climate%Wind_LR( vi, m), climate%Wind_DU( vi, m))
          self%Sublimation( vi, m) = wind_sens * wind_speed**2._dp * geom%Hs_slope( vi)

          ! Extract runoff and SMB
          self%Runoff( vi, m) = self%SurfaceMelt( vi, m) + self%Rainfall( vi, m) - self%Refreezing( vi, m)
          self%SMB_monthly( vi, m) = &
            self%Snowfall( vi, m) + self%Refreezing( vi, m) - self%SurfaceMelt( vi, m) - self%Sublimation( vi, m)

          ! Add this month's snow accumulation to next month's initial snow depth.
          if (geom%Hi( vi) > 0._dp) then

            ! Approximate surface snow density from Veldhuijzen et al. (2023)
            surface_snow_density = 376._dp + (sum(climate%T2m( vi, :))/12._dp - 235._dp) * 0.77_dp

            ! Determine firn layer temperature to ensure nonzero values during initialisation
            firn_temp = max(sum(climate%T2m( vi,:))/12._dp, ice%Ti( vi, 1))

            ! Compute the temperature-dependent exponent in Arthern et al. (2010),
            ! As used in Veldhuijzen et al. (2023)
            temp_exponent = exp(-Ec/(R_gas*climate%T2m( vi, m)) + Eg/(R_gas*firn_temp))

            ! Integrate firn air content over month:
            ! 1) Snowfall - melt adds a layer of firn at surface_snow_density, with all terms in mwe
            ! so added layer of firn = (S-M) * rho_fw/rho_ss (thickness of snowpack).
            ! Converted to FAC by multiplying with (rho_i-rho_ss)/rho_i (e.g., Kuipers Munnike et al., 2015)
            ! 2) Densification reduces FAC, following Arthern et al. (2010), and its rate scales with
            ! long-term accumulation (here taken as average monthly SMB, and with the temp_exponent

            self%FirnAirContent( vi, m) = max(0._dp, &
              self%FirnAirContent( vi, mprev) &
              + (self%Snowfall( vi, m) - self%SurfaceMelt( vi, m)) * freshwater_density / surface_snow_density &
                * (ice_density - surface_snow_density)/ice_density &
              - C%SMB_ITM_C_densification_rate * grav * sum(self%Snowfall( vi, :))/12._dp &
                * ice_density * temp_exponent * self%FirnAirContent( vi, mprev))

          else
            ! Ice free land
            self%FirnAirContent( vi, m) = 0._dp
          end if

        end do

        ! Integrate SMB over the full year and convert from water to ice equivalent
        self%SMB( vi) = sum(self%SMB_monthly( vi, :)) * freshwater_density / ice_density

      end if

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine SMB_model_ITM_v2_run

  subroutine SMB_model_ITM_v2_remap( self, mesh_new, time, refgeo_init, refgeo_PD, geom)

    ! In/output variables
    class(type_SMB_model_ITM_v2),        intent(inout) :: self
    type(type_mesh), target,               intent(in   ) :: mesh_new
    real(dp),                              intent(in   ) :: time
    type(type_reference_geometry), target, intent(in   ) :: refgeo_init, refgeo_PD
    class(atype_ice_geometry_model_data),  intent(in   ) :: geom

    ! Local variables:
    character(len=*), parameter :: routine_name = 'SMB_model_ITM_v2_remap'

    ! Add routine to path
    call init_routine( routine_name)

    ! Remap all the stuff that is specific to SMB model ITM_v2

    call self%remap_field( mesh_new, 'FirnAirContent'  , self%FirnAirContent   )
    call self%remap_field( mesh_new, 'Rainfall'        , self%Rainfall         )
    call self%remap_field( mesh_new, 'Snowfall'        , self%Snowfall         )
    call self%remap_field( mesh_new, 'SurfaceMelt'     , self%SurfaceMelt             )
    call self%remap_field( mesh_new, 'Refreezing'      , self%Refreezing       )
    call self%remap_field( mesh_new, 'Runoff'          , self%Runoff           )
    call self%remap_field( mesh_new, 'Albedo'          , self%Albedo           )
    call self%remap_field( mesh_new, 'SMB_monthly'     , self%SMB_monthly      )
    call self%remap_field( mesh_new, 'Sublimation'     , self%Sublimation      )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine SMB_model_ITM_v2_remap

  function get_SMB_model_name( self) result( SMB_model_name)
    class(type_SMB_model_ITM_v2), intent(in) :: self
    character(len=:), allocatable :: SMB_model_name
    SMB_model_name = 'ITM_v2'
  end function get_SMB_model_name

end module SMB_ITM_v2
