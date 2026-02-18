module SMB_IMAU_ITM

  use mpi_basic, only: par
  use UPSY_main, only: UPSY
  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use SMB_model_basic, only: atype_SMB_model, type_SMB_model_context_allocate, &
    type_SMB_model_context_initialise, type_SMB_model_context_run, &
    type_SMB_model_context_remap
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use mpi_f08, only: MPI_WIN
  use ice_model_types, only: type_ice_model
  use climate_model_types, only: type_climate_model
  use netcdf_io_main, only: read_field_from_file_2D, read_field_from_file_2D_monthly
  use parameters, only: ice_density, T0, L_fusion, sec_per_year

  implicit none

  private

  public :: type_SMB_model_IMAU_ITM

! ===== Types =====
! =================

  type, extends(atype_SMB_model) :: type_SMB_model_IMAU_ITM

      ! Main data fields
      real(dp), dimension(:  ), contiguous, pointer :: AlbedoSurf       => null() !< Surface albedo underneath the snow layer (water, rock or ice)
      real(dp), dimension(:  ), contiguous, pointer :: MeltPreviousYear => null() !< [m.w.e.] total melt in the previous year
      real(dp), dimension(:,:), contiguous, pointer :: FirnDepth        => null() !< [m] depth of the firn layer
      real(dp), dimension(:,:), contiguous, pointer :: Rainfall         => null() !< Monthly rainfall (m)
      real(dp), dimension(:,:), contiguous, pointer :: Snowfall         => null() !< Monthly snowfall (m)
      real(dp), dimension(:,:), contiguous, pointer :: AddedFirn        => null() !< Monthly added firn (m)
      real(dp), dimension(:,:), contiguous, pointer :: Melt             => null() !< Monthly melt (m)
      real(dp), dimension(:,:), contiguous, pointer :: Refreezing       => null() !< Monthly refreezing (m)
      real(dp), dimension(:  ), contiguous, pointer :: Refreezing_year  => null() !< Yearly  refreezing (m)
      real(dp), dimension(:,:), contiguous, pointer :: Runoff           => null() !< Monthly runoff (m)
      real(dp), dimension(:,:), contiguous, pointer :: Albedo           => null() !< Monthly albedo
      real(dp), dimension(:  ), contiguous, pointer :: Albedo_year      => null() !< Yearly albedo
      real(dp), dimension(:,:), contiguous, pointer :: SMB_monthly      => null() !< [m] Monthly SMB
      type(MPI_WIN) :: wAlbedoSurf, wMeltPreviousYear, wFirnDepth, wRainfall
      type(MPI_WIN) :: wSnowfall, wAddedFirn, wMelt, wRefreezing, wRefreezing_year
      type(MPI_WIN) :: wRunoff, wAlbedo, wAlbedo_year, wSMB_monthly

      ! Tuning parameters for the IMAU-ITM SMB model (different for each region, set from config)
      real(dp)  :: C_abl_constant
      real(dp)  :: C_abl_Ts
      real(dp)  :: C_abl_Q
      real(dp)  :: C_refr

      ! Ideally these parameters should not be region-dependent?
      real(dp)  :: albedo_water
      real(dp)  :: albedo_soil
      real(dp)  :: albedo_ice
      real(dp)  :: albedo_snow

    contains

      procedure, public :: allocate_SMB_model   => allocate_SMB_model_IMAU_ITM_abs
      procedure, public :: deallocate_SMB_model => deallocate_SMB_model_IMAU_ITM_abs
      procedure, public :: initialise_SMB_model => initialise_SMB_model_IMAU_ITM_abs
      procedure, public :: run_SMB_model        => run_SMB_model_IMAU_ITM_abs
      procedure, public :: remap_SMB_model      => remap_SMB_model_IMAU_ITM_abs

      procedure, private :: allocate_SMB_model_IMAU_ITM
      procedure, private :: initialise_SMB_model_IMAU_ITM
      procedure, private :: initialise_IMAU_ITM_firn_from_file
      procedure, private :: run_SMB_model_IMAU_ITM
      procedure, private :: remap_SMB_model_IMAU_ITM

  end type type_SMB_model_IMAU_ITM

contains

  subroutine allocate_SMB_model_IMAU_ITM_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_IMAU_ITM),                intent(inout) :: self
    type(type_SMB_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SMB_model_IMAU_ITM_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%allocate_SMB_model_IMAU_ITM

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_IMAU_ITM_abs

  subroutine deallocate_SMB_model_IMAU_ITM_abs( self)

    ! In/output variables:
    class(type_SMB_model_IMAU_ITM), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'deallocate_SMB_model_IMAU_ITM_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine deallocate_SMB_model_IMAU_ITM_abs

  subroutine initialise_SMB_model_IMAU_ITM_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_IMAU_ITM),                  intent(inout) :: self
    type(type_SMB_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_IMAU_ITM_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call self%initialise_SMB_model_IMAU_ITM( self%mesh, context%ice, self%region_name())

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_IMAU_ITM_abs

  subroutine run_SMB_model_IMAU_ITM_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_IMAU_ITM),           intent(inout) :: self
    type(type_SMB_model_context_run), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_IMAU_ITM_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call self%run_SMB_model_IMAU_ITM( self%mesh, context%ice, context%climate, self%region_name())

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_IMAU_ITM_abs

  subroutine remap_SMB_model_IMAU_ITM_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_IMAU_ITM),             intent(inout) :: self
    type(type_SMB_model_context_remap), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_SMB_model_IMAU_ITM_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%remap_SMB_model_IMAU_ITM( context%mesh_new)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_IMAU_ITM_abs



  subroutine allocate_SMB_model_IMAU_ITM( self)

    ! In/output variables
    class(type_SMB_model_IMAU_ITM), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SMB_model_IMAU_ITM'
    character(:), allocatable      :: choice_SMB_IMAU_ITM

    ! Add routine to path
    call init_routine( routine_name)

    call self%create_field( self%AlbedoSurf, self%wAlbedoSurf, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'AlbedoSurf', &
      long_name = 'Surface albedo underneath the snow layer (water, rock or ice)', &
      units     = '-')

    call self%create_field( self%MeltPreviousYear, self%wMeltPreviousYear, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'MeltPreviousYear', &
      long_name = 'Total melt in the previous year', &
      units     = 'm.w.e.')

    call self%create_field( self%FirnDepth, self%wFirnDepth, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'FirnDepth', &
      long_name = 'Depth of the firn layer', &
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

    call self%create_field( self%AddedFirn, self%wAddedFirn, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'AddedFirn', &
      long_name = 'Monthly added firn', &
      units     = 'm')

    call self%create_field( self%Melt, self%wMelt, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Melt', &
      long_name = 'Monthly melt', &
      units     = 'm')

    call self%create_field( self%Refreezing, self%wRefreezing, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Refreezing', &
      long_name = 'Monthly refreezing', &
      units     = 'm')

    call self%create_field( self%Refreezing_year, self%wRefreezing_year, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Refreezing_year', &
      long_name = 'Yearly refreezing', &
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

    call self%create_field( self%Albedo_year, self%wAlbedo_year, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Albedo_year', &
      long_name = 'Yearly albedo', &
      units     = '-')

    call self%create_field( self%SMB_monthly, self%wSMB_monthly, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'SMB_monthly', &
      long_name = 'Monthly SMB', &
      units     = '-')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_IMAU_ITM

  subroutine initialise_SMB_model_IMAU_ITM( self, mesh, ice, region_name)

    ! In/output variables
    class(type_SMB_model_IMAU_ITM), intent(inout) :: self
    type(type_mesh),                intent(in   ) :: mesh
    type(type_ice_model),           intent(in   ) :: ice
    character(len=3),               intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_IMAU_ITM'
    integer                        :: vi
    character(:), allocatable      :: choice_SMB_IMAUITM_init_firn

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine which constants to use for this region
    select case (region_name)
    case default
      call crash('unknown region_name "' // region_name // '"')
    case ('NAM')
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_NAM
      self%C_abl_constant       = C%SMB_IMAUITM_C_abl_constant_NAM
      self%C_abl_Ts             = C%SMB_IMAUITM_C_abl_Ts_NAM
      self%C_abl_Q              = C%SMB_IMAUITM_C_abl_Q_NAM
      self%C_refr               = C%SMB_IMAUITM_C_refr_NAM
    case ('EAS')
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_EAS
      self%C_abl_constant       = C%SMB_IMAUITM_C_abl_constant_EAS
      self%C_abl_Ts             = C%SMB_IMAUITM_C_abl_Ts_EAS
      self%C_abl_Q              = C%SMB_IMAUITM_C_abl_Q_EAS
      self%C_refr               = C%SMB_IMAUITM_C_refr_EAS
    case ('GRL')
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_GRL
      self%C_abl_constant       = C%SMB_IMAUITM_C_abl_constant_GRL
      self%C_abl_Ts             = C%SMB_IMAUITM_C_abl_Ts_GRL
      self%C_abl_Q              = C%SMB_IMAUITM_C_abl_Q_GRL
      self%C_refr               = C%SMB_IMAUITM_C_refr_GRL
    case ('ANT')
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_ANT
      self%C_abl_constant       = C%SMB_IMAUITM_C_abl_constant_ANT
      self%C_abl_Ts             = C%SMB_IMAUITM_C_abl_Ts_ANT
      self%C_abl_Q              = C%SMB_IMAUITM_C_abl_Q_ANT
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

      do vi = mesh%vi1, mesh%vi2
        if (ice%Hi( vi) > 0._dp) then
          self%FirnDepth       ( vi,:) = C%SMB_IMAUITM_initial_firn_thickness
          self%MeltPreviousYear( vi  ) = 0._dp
        else
          self%FirnDepth       ( vi,:) = 0._dp
          self%MeltPreviousYear( vi  ) = 0._dp
        end if
      end do

    case ('read_from_file')
      ! Initialise with the firn layer of a previous run
      call self%initialise_IMAU_ITM_firn_from_file( mesh, region_name)
    end select

    ! Initialise albedo
    do vi = mesh%vi1, mesh%vi2

      ! Background albedo
      if (ice%Hb( vi) < 0._dp) then
        self%AlbedoSurf( vi) = self%albedo_water
      else
        self%AlbedoSurf( vi) = self%albedo_soil
      end if
      if (ice%Hi( vi) > 0._dp) then
        self%AlbedoSurf( vi) = self%albedo_snow
      end if

      self%Albedo( vi,:) = self%AlbedoSurf( vi)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_IMAU_ITM

  subroutine initialise_IMAU_ITM_firn_from_file( self, mesh, region_name)
    !< Initialise the firn depth and meltpreviousyear data from a NetCDF file

    ! In/output variables
    class(type_SMB_model_IMAU_ITM), intent(inout) :: self
    type(type_mesh),                intent(in   ) :: mesh
    character(len=3),               intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_IMAU_ITM_firn_from_file'
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
      call read_field_from_file_2D_monthly( filename_restart_firn, 'FirnDepth', mesh, C%output_dir, self%FirnDepth)
      call read_field_from_file_2D( filename_restart_firn, 'MeltPreviousYear', mesh, C%output_dir, self%MeltPreviousYear)
    else
      ! Assume the file has a time dimension, and read the specified timeframe
      call read_field_from_file_2D_monthly( filename_restart_firn, 'FirnDepth', mesh, C%output_dir, self%FirnDepth, time_to_read = timeframe_restart_firn)
      call read_field_from_file_2D( filename_restart_firn, 'MeltPreviousYear', mesh, C%output_dir, self%MeltPreviousYear, time_to_read = timeframe_restart_firn)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_IMAU_ITM_firn_from_file

  subroutine run_SMB_model_IMAU_ITM( self, mesh, ice, climate, region_name)

    ! NOTE: all the SMB components are in meters of water equivalent;
    !       the end result (SMB_monthly and SMB) are in meters of ice equivalent.

    ! In/output variables:
    class(type_SMB_model_IMAU_ITM), intent(inout) :: self
    type(type_mesh),                intent(in   ) :: mesh
    type(type_ice_model),           intent(in   ) :: ice
    type(type_climate_model),       intent(in   ) :: climate
    character(len=3),               intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_IMAU_ITM'
    integer                        :: vi
    integer                        :: m, mprev
    real(dp)                       :: snowfrac, liquid_water, sup_imp_wat

    ! Add routine to call stack
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2

      ! Background albedo
      self%AlbedoSurf( vi) = self%albedo_soil
      if ((ice%mask_icefree_ocean( vi) .eqv. .true. .and. ice%mask_floating_ice( vi) .eqv. .false.) .or. ice%mask_noice( vi) .eqv. .true.) self%AlbedoSurf( vi) = self%albedo_water
      if (ice%mask_grounded_ice(   vi) .eqv. .true. .or.  ice%mask_floating_ice( vi) .eqv. .true.) self%AlbedoSurf( vi) = self%albedo_ice

      do m = 1, 12  ! Month loop

        mprev = m - 1
        if (mprev==0) mprev = 12

        self%Albedo( vi,m) = min( self%albedo_snow, max( self%AlbedoSurf( vi), self%albedo_snow - (self%albedo_snow - self%AlbedoSurf( vi))  * &
                             exp(-15._dp * self%FirnDepth( vi,mprev)) - 0.015_dp * self%MeltPreviousYear( vi)))
        if ((ice%mask_icefree_ocean( vi) .eqv. .true. .and. ice%mask_floating_ice( vi) .eqv. .false.) .or. ice%mask_noice( vi) .eqv. .true.) self%Albedo( vi,m) = self%albedo_water

        ! Determine ablation as a function of surface temperature and albedo/insolation according to Bintanja et al. (2002)
        self%Melt( vi,m) = max(0._dp, ( self%C_abl_Ts         * (climate%T2m( vi,m) - T0) + &
                                        self%C_abl_Q          * (1.0_dp - self%Albedo( vi,m)) * climate%Q_TOA( vi,m) - &
                                        self%C_abl_constant)  * sec_per_year / (L_fusion * 1000._dp * 12._dp))

        ! Determine accumulation with snow/rain fraction from Ohmura et al. (1999), liquid water content (rain and melt water) and snow depth

        ! NOTE: commented version is the old ANICE version, supposedly based on "physics" (which we cant check), but
        !       the new version was tuned to RACMO output and produced significantly better snow fractions...
        ! However there is still snowfall even if temperatures are at 300 K, which does not seem realistic.
        snowfrac = max(0._dp, min(1._dp, 0.5_dp   * (1 - atan((climate%T2m(vi,m) - T0) / 3.5_dp)  / 1.25664_dp))) ! ANICE "realistic" snow fractions
        !snowfrac = max(0._dp, min(1._dp, 0.725_dp * (1 - atan((climate%T2m( vi,m) - T0) / 5.95_dp) / 1.8566_dp))) ! IMAU-ICE "tuned" snow fractions

        self%Snowfall( vi,m) = climate%Precip( vi,m) *          snowfrac
        self%Rainfall( vi,m) = climate%Precip( vi,m) * (1._dp - snowfrac)

        ! Refreezing according to Janssens & Huybrechts (2000)
        ! The refreezing (=effective retention) is the minimum value of the amount of super imposed
        ! water and the available liquid water, with a maximum value of the total precipitation.
        ! (see also Huybrechts & de Wolde, 1999)

        ! Add this month's snow accumulation to next month's initial snow depth.
        self%AddedFirn( vi,m) = self%Snowfall( vi,m) - self%Melt( vi,m)
        self%FirnDepth( vi,m) = min( 10._dp, max( 0._dp, self%FirnDepth( vi,mprev) + self%AddedFirn( vi,m) ))

      end do

      ! Calculate refreezing for the whole year, divide equally over the 12 months, then calculate resulting runoff and SMB.
      ! This resolves the problem with refreezing where liquid water is mostly available in summer
      ! but "refreezing potential" mostly in winter, and there is no proper meltwater retention.
      sup_imp_wat  = self%C_refr * max( 0._dp, T0 - sum( climate%T2m( vi,:)) / 12._dp)
      liquid_water = sum( self%Rainfall( vi,:)) + sum( self%Melt( vi,:))

      ! Note: Refreezing is limited by the ability of the firn layer to store melt water. currently a ten meter firn layer can store
      ! 2.5 m of water. However, this is based on expert judgement, NOT empirical evidence.
      self%Refreezing_year( vi) = min( min( min( sup_imp_wat, liquid_water), sum(climate%Precip( vi,:))), 0.25_dp * sum( self%FirnDepth( vi,:) / 12._dp)) ! version from IMAU-ICE dev branch
      !SMB%Refreezing_year( vi) = min( min( sup_imp_wat, liquid_water), sum( climate%Precip( vi,:))) ! outdated version on main branch

      if (ice%mask_grounded_ice( vi)  .eqv. .false. .or. ice%mask_floating_ice( vi) .eqv. .false.) self%Refreezing_year( vi) = 0._dp
      if (ice%mask_icefree_ocean( vi) .eqv. .true.)                                                self%AddedFirn( vi,:)     = 0._dp ! Does it make sense to add firn over the ocean?!

      do m = 1, 12
        self%Refreezing(  vi,m) = self%Refreezing_year( vi) / 12._dp
        self%Runoff(      vi,m) = self%Melt( vi,m) + self%Rainfall( vi,m) - self%Refreezing( vi,m)
        self%SMB_monthly( vi,m) = self%Snowfall( vi,m) + self%Refreezing( vi,m) - self%Melt( vi,m)
      end do

      self%SMB( vi) = sum( self%SMB_monthly( vi,:))
      !if (ice%mask_icefree_ocean( vi) .eqv. .true.) SMB%SMB( vi) = 0._dp ! should we limit SMB over open ocean?

      ! Calculate total melt over this year, to be used for determining next year's albedo
      self%MeltPreviousYear( vi) = sum( self%Melt( vi,:))

    end do

    ! Convert final SMB from water to ice equivalent
    self%SMB_monthly( mesh%vi1:mesh%vi2,:) = self%SMB_monthly(  mesh%vi1:mesh%vi2,:) * 1000._dp / ice_density
    self%SMB(         mesh%vi1:mesh%vi2  ) = self%SMB(          mesh%vi1:mesh%vi2  ) * 1000._dp / ice_density

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_IMAU_ITM

  subroutine remap_SMB_model_IMAU_ITM( self, mesh_new)

    ! In/output variables
    class(type_SMB_model_IMAU_ITM), intent(inout) :: self
    type(type_mesh),                intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_SMB_model_IMAU_ITM'
    character(:), allocatable      :: choice_SMB_IMAU_ITM

    ! Add routine to path
    call init_routine( routine_name)

    call self%remap_field( mesh_new, 'AlbedoSurf'      , self%AlbedoSurf       )
    call self%remap_field( mesh_new, 'MeltPreviousYear', self%MeltPreviousYear )
    call self%remap_field( mesh_new, 'FirnDepth'       , self%FirnDepth        )
    call self%remap_field( mesh_new, 'Rainfall'        , self%Rainfall         )
    call self%remap_field( mesh_new, 'Snowfall'        , self%Snowfall         )
    call self%remap_field( mesh_new, 'AddedFirn'       , self%AddedFirn        )
    call self%remap_field( mesh_new, 'Melt'            , self%Melt             )
    call self%remap_field( mesh_new, 'Refreezing'      , self%Refreezing       )
    call self%remap_field( mesh_new, 'Refreezing_year' , self%Refreezing_year  )
    call self%remap_field( mesh_new, 'Runoff'          , self%Runoff           )
    call self%remap_field( mesh_new, 'Albedo'          , self%Albedo           )
    call self%remap_field( mesh_new, 'Albedo_year'     , self%Albedo_year      )
    call self%remap_field( mesh_new, 'SMB_monthly'     , self%SMB_monthly      )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_IMAU_ITM

end module SMB_IMAU_ITM