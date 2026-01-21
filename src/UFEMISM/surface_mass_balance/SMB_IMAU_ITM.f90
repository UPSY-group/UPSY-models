module SMB_IMAU_ITM

  ! IMAU-ITM SMB model

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use climate_model_types, only: type_climate_model
  use parameters, only: T0, L_fusion, sec_per_year, pi, ice_density
  use netcdf_io_main
  use global_forcing_types
  use mpi_f08, only: MPI_WIN
  use SMB_basic, only: atype_SMB_model, type_SMB_model_context_allocate, &
    type_SMB_model_context_initialise, type_SMB_model_context_run, type_SMB_model_context_remap
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension

  implicit none

  private

  public :: type_SMB_model_IMAU_ITM

  type, extends(atype_SMB_model) :: type_SMB_model_IMAU_ITM
    !< The IMAU Insolation-Temperature Model

      ! Main data fields
      real(dp), dimension(:  ), contiguous, pointer :: AlbedoSurf              ! Surface albedo underneath the snow layer (water, rock or ice)
      real(dp), dimension(:  ), contiguous, pointer :: MeltPreviousYear        ! [m.w.e.] total melt in the previous year
      real(dp), dimension(:,:), contiguous, pointer :: FirnDepth               ! [m] depth of the firn layer
      real(dp), dimension(:,:), contiguous, pointer :: Rainfall                ! Monthly rainfall (m)
      real(dp), dimension(:,:), contiguous, pointer :: Snowfall                ! Monthly snowfall (m)
      real(dp), dimension(:,:), contiguous, pointer :: AddedFirn               ! Monthly added firn (m)
      real(dp), dimension(:,:), contiguous, pointer :: Melt                    ! Monthly melt (m)
      real(dp), dimension(:,:), contiguous, pointer :: Refreezing              ! Monthly refreezing (m)
      real(dp), dimension(:  ), contiguous, pointer :: Refreezing_year         ! Yearly  refreezing (m)
      real(dp), dimension(:,:), contiguous, pointer :: Runoff                  ! Monthly runoff (m)
      real(dp), dimension(:,:), contiguous, pointer :: Albedo                  ! Monthly albedo
      real(dp), dimension(:  ), contiguous, pointer :: Albedo_year             ! Yearly albedo
      real(dp), dimension(:,:), contiguous, pointer :: SMB_monthly             ! [m] Monthly SMB
      type(MPI_WIN) :: wAlbedoSurf, wMeltPreviousYear, wFirnDepth, wRainfall
      type(MPI_WIN) :: wSnowfall, wAddedFirn, wMelt, wRefreezing, wRefreezing_year
      type(MPI_WIN) :: wRunoff, wAlbedo, wAlbedo_year, wSMB_monthly

      ! Tuning parameters for the IMAU-ITM SMB model (different for each region, set from config)
      real(dp) :: C_abl_constant
      real(dp) :: C_abl_Ts
      real(dp) :: C_abl_Q
      real(dp) :: C_refr

      ! Ideally these parameters should not be region-dependent?
      real(dp) :: albedo_water
      real(dp) :: albedo_soil
      real(dp) :: albedo_ice
      real(dp) :: albedo_snow

    contains

      procedure, public :: allocate_SMB_model   => allocate_SMB_model_IMAU_ITM
      procedure, public :: initialise_SMB_model => initialise_SMB_model_IMAU_ITM
      procedure, public :: run_SMB_model        => run_SMB_model_IMAU_ITM_abs
      procedure, public :: remap_SMB_model      => remap_SMB_model_IMAU_ITM

      procedure, private :: initialise_IMAUITM_firn_from_file
      procedure, private :: run_SMB_model_IMAU_ITM

  end type type_SMB_model_IMAU_ITM

contains

  subroutine allocate_SMB_model_IMAU_ITM( self, context)

    ! In/output variables:
    class(type_SMB_model_IMAU_ITM),        intent(inout) :: self
    type(type_SMB_model_context_allocate), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SMB_model_IMAU_ITM'

    ! Add routine to path
    call init_routine( routine_name)

    call self%set_name('SMB_model_IMAU_ITM')

    ! Create all model fields
    ! =======================

    call self%create_field( self%AlbedoSurf, self%wAlbedoSurf, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'AlbedoSurf', &
      long_name = 'background albedo', &
      units     = '-')

    call self%create_field( self%Rainfall, self%wRainfall, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Rainfall', &
      long_name = 'monthly rainfall', &
      units     = 'm')

    call self%create_field( self%Snowfall, self%wSnowfall, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Snowfall', &
      long_name = 'monthly snowfall', &
      units     = 'm')

    call self%create_field( self%AddedFirn, self%wAddedFirn, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'AddedFirn', &
      long_name = 'monthly added firn', &
      units     = 'm')

    call self%create_field( self%Melt, self%wMelt, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Melt', &
      long_name = 'monthly melt', &
      units     = 'm')

    call self%create_field( self%Refreezing, self%wRefreezing, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Refreezing', &
      long_name = 'monthly refreezing', &
      units     = 'm')

    call self%create_field( self%Refreezing_year, self%wRefreezing_year, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Refreezing_year', &
      long_name = 'annual refreezing', &
      units     = 'm')

    call self%create_field( self%Runoff, self%wRunoff, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Runoff', &
      long_name = 'monthly runoff', &
      units     = 'm')

    call self%create_field( self%Albedo, self%wAlbedo, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'Albedo', &
      long_name = 'monthly albedo', &
      units     = '-')

    call self%create_field( self%Albedo_year, self%wAlbedo_year, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Albedo_year', &
      long_name = 'annual albedo', &
      units     = '-')

    call self%create_field( self%SMB_monthly, self%wSMB_monthly, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'SMB_monthly', &
      long_name = 'monthly surface mass balance', &
      units     = 'm')

    call self%create_field( self%FirnDepth, self%wFirnDepth, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'FirnDepth', &
      long_name = 'monthly firn depth', &
      units     = 'm')

    call self%create_field( self%MeltPreviousYear, self%wMeltPreviousYear, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'MeltPreviousYear', &
      long_name = 'total melt in previous year', &
      units     = 'm')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_IMAU_ITM

  subroutine initialise_SMB_model_IMAU_ITM( self, context)

    ! In/output variables:
    class(type_SMB_model_IMAU_ITM),          intent(inout) :: self
    type(type_SMB_model_context_initialise), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_IMAU_ITM'
    character(:), allocatable      :: choice_SMB_IMAUITM_init_firn
    integer                        :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Determine which constants to use for this region
    select case (context%region_name)
    case default
      call crash('invalid region name "' // context%region_name // '"')
    case ('NAM')
      self%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_NAM
      self%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_NAM
      self%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_NAM
      self%C_refr         = C%SMB_IMAUITM_C_refr_NAM
    case ('EAS')
      self%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_EAS
      self%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_EAS
      self%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_EAS
      self%C_refr         = C%SMB_IMAUITM_C_refr_EAS
    case ('GRL')
      self%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_GRL
      self%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_GRL
      self%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_GRL
      self%C_refr         = C%SMB_IMAUITM_C_refr_GRL
    case ('ANT')
      self%C_abl_constant = C%SMB_IMAUITM_C_abl_constant_ANT
      self%C_abl_Ts       = C%SMB_IMAUITM_C_abl_Ts_ANT
      self%C_abl_Q        = C%SMB_IMAUITM_C_abl_Q_ANT
      self%C_refr         = C%SMB_IMAUITM_C_refr_ANT
    end select

    ! Initialising albedo values
    self%albedo_water = C%SMB_IMAUITM_albedo_water
    self%albedo_soil  = C%SMB_IMAUITM_albedo_soil
    self%albedo_ice   = C%SMB_IMAUITM_albedo_ice
    self%albedo_snow  = C%SMB_IMAUITM_albedo_snow

    ! Initialisation choice
    select case (context%region_name)
    case default
      call crash('invalid region name "' // context%region_name // '"')
    case ('NAM')
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_NAM
    case ('EAS')
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_EAS
    case ('GRL')
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_GRL
    case ('ANT')
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_ANT
    end select

    ! Initialise the firn layer
    select case (choice_SMB_IMAUITM_init_firn)
    case default
      call crash('unknown choice_SMB_IMAUITM_init_firn "' // trim( choice_SMB_IMAUITM_init_firn) // '"')
    case ('uniform')
      ! Initialise with a uniform firn layer over the ice sheet

      do vi = self%mesh%vi1, self%mesh%vi2
        if (context%ice%Hi( vi) > 0._dp) then
          self%FirnDepth(        vi,:) = C%SMB_IMAUITM_initial_firn_thickness
          self%MeltPreviousYear( vi  ) = 0._dp
        else
          self%FirnDepth(        vi,:) = 0._dp
          self%MeltPreviousYear( vi  ) = 0._dp
        end if
      end do

    case ('read_from_file')
      ! Initialise with the firn layer of a previous run
      call self%initialise_IMAUITM_firn_from_file( self%mesh, context%region_name)
    end select

    ! Initialise albedo
    do vi = self%mesh%vi1, self%mesh%vi2
      ! Background albedo
      if (context%ice%Hb( vi) < 0._dp) then
        self%AlbedoSurf( vi) = self%albedo_water
      else
        self%AlbedoSurf( vi) = self%albedo_soil
      end if
      if (context%ice%Hi( vi) > 0._dp) then
        self%AlbedoSurf(  vi) = self%albedo_snow
      end if
      self%Albedo( vi,:) = self%AlbedoSurf( vi)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_IMAU_ITM

  subroutine run_SMB_model_IMAU_ITM_abs( self, context)

    ! In/output variables:
    class(type_SMB_model_IMAU_ITM),   intent(inout) :: self
    type(type_SMB_model_context_run), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_IMAU_ITM_abs'

    ! Add routine to path
    call init_routine( routine_name)

    call self%run_SMB_model_IMAU_ITM( self%mesh, context%ice, context%climate)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_IMAU_ITM_abs

  subroutine remap_SMB_model_IMAU_ITM( self, context)

    ! In- and output variables
    class(type_SMB_model_IMAU_ITM),     intent(inout) :: self
    type(type_SMB_model_context_remap), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_SMB_model_IMAU_ITM'

    ! Add routine to path
    call init_routine( routine_name)

    call self%remap_field( context%mesh_new, 'AlbedoSurf'      , self%AlbedoSurf      )
    call self%remap_field( context%mesh_new, 'MeltPreviousYear', self%MeltPreviousYear)
    call self%remap_field( context%mesh_new, 'Refreezing_year' , self%Refreezing_year )
    call self%remap_field( context%mesh_new, 'Albedo_year'     , self%Albedo_year     )
    call self%remap_field( context%mesh_new, 'FirnDepth'       , self%FirnDepth       )
    call self%remap_field( context%mesh_new, 'Rainfall'        , self%Rainfall        )
    call self%remap_field( context%mesh_new, 'Snowfall'        , self%Snowfall        )
    call self%remap_field( context%mesh_new, 'AddedFirn'       , self%AddedFirn       )
    call self%remap_field( context%mesh_new, 'Melt'            , self%Melt            )
    call self%remap_field( context%mesh_new, 'Refreezing'      , self%Refreezing      )
    call self%remap_field( context%mesh_new, 'Runoff'          , self%Runoff          )
    call self%remap_field( context%mesh_new, 'Albedo'          , self%Albedo          )
    call self%remap_field( context%mesh_new, 'SMB_monthly'     , self%SMB_monthly     )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_IMAU_ITM



  subroutine run_SMB_model_IMAU_ITM( self, mesh, ice, climate)
    ! Run the IMAU-ITM SMB model.

    ! NOTE: all the SMB components are in meters of water equivalent;
    !       the end result (SMB and SMB_year) are in meters of ice equivalent.

    ! In- and output variables
    class(type_SMB_model_IMAU_ITM), intent(inout) :: self
    type(type_mesh),                intent(in)    :: mesh
    type(type_ice_model),           intent(in)    :: ice
    type(type_climate_model),       intent(in)    :: climate

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_IMAU_ITM'
    integer                        :: vi
    integer                        :: m,mprev
    real(dp)                       :: snowfrac, liquid_water, sup_imp_wat

    ! Add routine to path
    call init_routine( routine_name)

    DO vi = mesh%vi1, mesh%vi2

      ! Background albedo
      self%AlbedoSurf( vi) = self%albedo_soil
      IF ((ice%mask_icefree_ocean( vi) .eqv. .TRUE. .AND. ice%mask_floating_ice( vi) .eqv. .FALSE.) .OR. ice%mask_noice( vi) .eqv. .TRUE.) self%AlbedoSurf( vi) = self%albedo_water
      IF (ice%mask_grounded_ice(   vi) .eqv. .TRUE. .OR. ice%mask_floating_ice(  vi) .eqv. .TRUE.) self%AlbedoSurf( vi) = self%albedo_ice

      DO m = 1, 12

        mprev = m - 1
        IF (mprev==0) mprev = 12

        self%Albedo( vi,m) = MIN(self%albedo_snow, MAX( self%AlbedoSurf( vi), self%albedo_snow - (self%albedo_snow - self%AlbedoSurf( vi))  * &
                             EXP(-15._dp * self%FirnDepth( vi,mprev)) - 0.015_dp * self%MeltPreviousYear( vi)))
        IF ((ice%mask_icefree_ocean( vi) .eqv. .TRUE. .AND. ice%mask_floating_ice( vi) .eqv. .FALSE.) .OR. ice%mask_noice( vi) .eqv. .TRUE.) self%Albedo( vi,m) = self%albedo_water

        ! Determine ablation as a function of surface temperature and albedo/insolation according to Bintanja et al. (2002)
        self%Melt( vi,m) = MAX(0._dp, ( self%C_abl_Ts         * (climate%T2m( vi,m) - T0) + &
                                       self%C_abl_Q          * (1.0_dp - self%Albedo( vi,m)) * climate%Q_TOA( vi,m) - &
                                       self%C_abl_constant)  * sec_per_year / (L_fusion * 1000._dp * 12._dp))

        ! Determine accumulation with snow/rain fraction from Ohmura et al. (1999), liquid water content (rain and melt water) and snow depth

        ! NOTE: commented version is the old ANICE version, supposedly based on "physics" (which we cant check), but
        !       the new version was tuned to RACMO output and produced significantly better snow fractions...
        ! However there is still snowfall even if temperatures are at 300 K, which does not seem realistic.
        snowfrac = MAX(0._dp, MIN(1._dp, 0.5_dp   * (1 - ATAN((climate%T2m(vi,m) - T0) / 3.5_dp)  / 1.25664_dp))) ! ANICE "realistic" snow fractions
        !snowfrac = MAX(0._dp, MIN(1._dp, 0.725_dp * (1 - ATAN((climate%T2m( vi,m) - T0) / 5.95_dp) / 1.8566_dp))) ! IMAU-ICE "tuned" snow fractions

        self%Snowfall( vi,m) = climate%Precip( vi,m) *          snowfrac
        self%Rainfall( vi,m) = climate%Precip( vi,m) * (1._dp - snowfrac)

        ! Refreezing according to Janssens & Huybrechts (2000)
        ! The refreezing (=effective retention) is the minimum value of the amount of super imposed
        ! water and the available liquid water, with a maximum value of the total precipitation.
        ! (see also Huybrechts & de Wolde, 1999)

        ! Add this month's snow accumulation to next month's initial snow depth.
        self%AddedFirn( vi,m) = self%Snowfall( vi,m) - self%Melt( vi,m)
        self%FirnDepth( vi,m) = MIN(10._dp, MAX(0._dp, self%FirnDepth( vi,mprev) + self%AddedFirn( vi,m) ))

      END DO ! DO m = 1, 12

      ! Calculate refreezing for the whole year, divide equally over the 12 months, then calculate resulting runoff and SMB.
      ! This resolves the problem with refreezing where liquid water is mostly available in summer
      ! but "refreezing potential" mostly in winter, and there is no proper meltwater retention.
      sup_imp_wat  = self%C_refr * MAX(0._dp, T0 - SUM(climate%T2m( vi,:))/12._dp)
      liquid_water = SUM(self%Rainfall( vi,:)) + SUM(self%Melt( vi,:))

      ! Note: Refreezing is limited by the ability of the firn layer to store melt water. currently a ten meter firn layer can store
      ! 2.5 m of water. However, this is based on expert judgement, NOT empirical evidence.
      self%Refreezing_year( vi) = MIN( MIN( MIN( sup_imp_wat, liquid_water), SUM(climate%Precip( vi,:))), 0.25_dp * SUM(self%FirnDepth( vi,:)/12._dp)) ! version from IMAU-ICE dev branch
      !SMB%Refreezing_year( vi) = MIN( MIN( sup_imp_wat, liquid_water), SUM(climate%Precip( vi,:))) ! outdated version on main branch

      IF (ice%mask_grounded_ice( vi)  .eqv. .FALSE. .OR. ice%mask_floating_ice( vi) .eqv. .FALSE.) self%Refreezing_year( vi) = 0._dp
      IF (ice%mask_icefree_ocean( vi) .eqv. .TRUE.)                                                self%AddedFirn( vi,:)     = 0._dp ! Does it make sense to add firn over the ocean?!



      DO m = 1, 12
        self%Refreezing(  vi,m) = self%Refreezing_year( vi) / 12._dp
        self%Runoff(      vi,m) = self%Melt( vi,m) + self%Rainfall( vi,m) - self%Refreezing( vi,m)
        self%SMB_monthly( vi,m) = self%Snowfall( vi,m) + self%Refreezing( vi,m) - self%Melt( vi,m)
      END DO

      !IF (ice%mask_icefree_ocean( vi) .eqv. .TRUE.) SMB%SMB( vi) = 0._dp ! should we limit SMB over open ocean?

      ! Calculate total melt over this year, to be used for determining next year's albedo
      self%MeltPreviousYear( vi) = SUM(self%Melt( vi,:))

    END DO

    ! Convert final SMB from water to ice equivalent
    self%SMB_monthly( mesh%vi1:mesh%vi2,:) = self%SMB_monthly(  mesh%vi1:mesh%vi2,:) * 1000._dp / ice_density

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_IMAU_ITM

  subroutine initialise_IMAUITM_firn_from_file( self, mesh, region_name)
    ! If this is a restarted run, read the firn depth and meltpreviousyear data from the restart file

    ! In/output variables
    class(type_SMB_model_IMAU_ITM), intent(inout) :: self
    type(type_mesh),                intent(in   ) :: mesh
    character(len=3),               intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_IMAUITM_firn_from_file'
    character(:), allocatable      :: filename_restart_firn
    real(dp)                       :: timeframe_restart_firn

    ! Add routine to path
    call init_routine( routine_name)

    ! Assume that SMB and geometry are read from the same restart file
    select case (region_name)
    case default
      call crash('unknown region_name "' // trim( region_name) // '"')
    case ('NAM')
      filename_restart_firn = C%filename_firn_IMAUITM_NAM
      timeframe_restart_firn = C%timeframe_restart_firn_IMAUITM_NAM
    case ('EAS')
      filename_restart_firn = C%filename_firn_IMAUITM_EAS
      timeframe_restart_firn = C%timeframe_restart_firn_IMAUITM_EAS
    case ('GRL')
      filename_restart_firn = C%filename_firn_IMAUITM_GRL
      timeframe_restart_firn = C%timeframe_restart_firn_IMAUITM_GRL
    case ('ANT')
      filename_restart_firn = C%filename_firn_IMAUITM_ANT
      timeframe_restart_firn = C%timeframe_restart_firn_IMAUITM_ANT
    end select

     ! Print to terminal
    if (par%primary) write(*,"(A)") '   Initialising SMB-model firn layer from file "' // &
      colour_string( trim( filename_restart_firn),'light blue') // '"...'

    ! Read firn layer from file
    if (timeframe_restart_firn == 1E9_dp) then
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

  end subroutine initialise_IMAUITM_firn_from_file

end module SMB_IMAU_ITM