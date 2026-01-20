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
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use reallocate_dist_shared_mod, only: reallocate_dist_shared
  use mpi_f08, only: MPI_WIN
  use SMB_basic, only: atype_SMB_model

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

      procedure, public  :: init, run, remap
      procedure, private :: initialise_IMAUITM_firn_from_file

  end type type_SMB_model_IMAU_ITM

contains

  subroutine run( self, mesh, ice, climate)
    ! Run the IMAU-ITM SMB model.

    ! NOTE: all the SMB components are in meters of water equivalent;
    !       the end result (SMB and SMB_year) are in meters of ice equivalent.

    ! In- and output variables
    class(type_SMB_model_IMAU_ITM), intent(inout) :: self
    type(type_mesh),                intent(in)    :: mesh
    type(type_ice_model),           intent(in)    :: ice
    type(type_climate_model),       intent(in)    :: climate

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_IMAUITM'
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

  end subroutine run

  subroutine init( self, mesh, ice, region_name)

    ! In/output variables
    class(type_SMB_model_IMAU_ITM),    INTENT(INOUT) :: self
    TYPE(type_mesh),                   INTENT(IN)    :: mesh
    TYPE(type_ice_model),              INTENT(IN)    :: ice
    CHARACTER(LEN=3),                  INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_SMB_model_IMAUITM'
    INTEGER                                            :: vi
    CHARACTER(LEN=256)                                 :: choice_SMB_IMAUITM_init_firn

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which constants to use for this region
    IF     (region_name == 'NAM') THEN
      self%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_NAM
      self%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_NAM
      self%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_NAM
      self%C_refr                   = C%SMB_IMAUITM_C_refr_NAM
    ELSEIF (region_name == 'EAS') THEN
      self%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_EAS
      self%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_EAS
      self%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_EAS
      self%C_refr                   = C%SMB_IMAUITM_C_refr_EAS
    ELSEIF (region_name == 'GRL') THEN
      self%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_GRL
      self%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_GRL
      self%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_GRL
      self%C_refr                   = C%SMB_IMAUITM_C_refr_GRL
    ELSEIF (region_name == 'ANT') THEN
      self%C_abl_constant           = C%SMB_IMAUITM_C_abl_constant_ANT
      self%C_abl_Ts                 = C%SMB_IMAUITM_C_abl_Ts_ANT
      self%C_abl_Q                  = C%SMB_IMAUITM_C_abl_Q_ANT
      self%C_refr                   = C%SMB_IMAUITM_C_refr_ANT
    END IF

    ! Initialising albedo values
    self%albedo_water        = C%SMB_IMAUITM_albedo_water
    self%albedo_soil         = C%SMB_IMAUITM_albedo_soil
    self%albedo_ice          = C%SMB_IMAUITM_albedo_ice
    self%albedo_snow         = C%SMB_IMAUITM_albedo_snow

    ! Allocating necessary fields
    call allocate_dist_shared( self%SMB, self%wSMB, mesh%pai_V%n_nih)
    self%SMB( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => self%SMB

    call allocate_dist_shared( self%AlbedoSurf       , self%wAlbedoSurf      , mesh%pai_V%n_nih    )
    call allocate_dist_shared( self%Rainfall         , self%wRainfall        , mesh%pai_V%n_nih, 12)
    call allocate_dist_shared( self%Snowfall         , self%wSnowfall        , mesh%pai_V%n_nih, 12)
    call allocate_dist_shared( self%AddedFirn        , self%wAddedFirn       , mesh%pai_V%n_nih, 12)
    call allocate_dist_shared( self%Melt             , self%wMelt            , mesh%pai_V%n_nih, 12)
    call allocate_dist_shared( self%Refreezing       , self%wRefreezing      , mesh%pai_V%n_nih, 12)
    call allocate_dist_shared( self%Refreezing_year  , self%wRefreezing_year , mesh%pai_V%n_nih    )
    call allocate_dist_shared( self%Runoff           , self%wRunoff          , mesh%pai_V%n_nih, 12)
    call allocate_dist_shared( self%Albedo           , self%wAlbedo          , mesh%pai_V%n_nih, 12)
    call allocate_dist_shared( self%Albedo_year      , self%wAlbedo_year     , mesh%pai_V%n_nih    )
    call allocate_dist_shared( self%SMB_monthly      , self%wSMB_monthly     , mesh%pai_V%n_nih, 12)
    call allocate_dist_shared( self%FirnDepth        , self%wFirnDepth       , mesh%pai_V%n_nih, 12)
    call allocate_dist_shared( self%MeltPreviousYear , self%wMeltPreviousYear, mesh%pai_V%n_nih    )

    self%AlbedoSurf      ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih      ) => self%AlbedoSurf
    self%Rainfall        ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:12) => self%Rainfall
    self%Snowfall        ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:12) => self%Snowfall
    self%AddedFirn       ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:12) => self%AddedFirn
    self%Melt            ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:12) => self%Melt
    self%Refreezing      ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:12) => self%Refreezing
    self%Refreezing_year ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih      ) => self%Refreezing_year
    self%Runoff          ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:12) => self%Runoff
    self%Albedo          ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:12) => self%Albedo
    self%Albedo_year     ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih      ) => self%Albedo_year
    self%SMB_monthly     ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:12) => self%SMB_monthly
    self%FirnDepth       ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:12) => self%FirnDepth
    self%MeltPreviousYear( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih      ) => self%MeltPreviousYear

    ! Initialisation choice
    IF     (region_name == 'NAM') THEN
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_NAM
    ELSEIF (region_name == 'EAS') THEN
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_EAS
    ELSEIF (region_name == 'GRL') THEN
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_GRL
    ELSEIF (region_name == 'ANT') THEN
      choice_SMB_IMAUITM_init_firn = C%choice_SMB_IMAUITM_init_firn_ANT
    END IF

    ! Initialise the firn layer
    IF     (choice_SMB_IMAUITM_init_firn == 'uniform') THEN
      ! Initialise with a uniform firn layer over the ice sheet

      DO vi = mesh%vi1, mesh%vi2
        IF (ice%Hi( vi) > 0._dp) THEN
          self%FirnDepth(        vi,:) = C%SMB_IMAUITM_initial_firn_thickness
          self%MeltPreviousYear(   vi) = 0._dp
        ELSE
          self%FirnDepth(        vi,:) = 0._dp
          self%MeltPreviousYear(   vi) = 0._dp
        END IF
      END DO

    ELSEIF (choice_SMB_IMAUITM_init_firn == 'read_from_file') THEN
      ! Initialise with the firn layer of a previous run
      CALL self%initialise_IMAUITM_firn_from_file( mesh, region_name)
    ELSE
      CALL crash('unknown choice_SMB_IMAUITM_init_firn "' // TRIM( choice_SMB_IMAUITM_init_firn) // '"!')
    END IF

    ! Initialise albedo
    DO vi = mesh%vi1, mesh%vi2
      ! Background albedo
      IF (ice%Hb( vi) < 0._dp) THEN
        self%AlbedoSurf( vi) = self%albedo_water
      ELSE
        self%AlbedoSurf( vi) = self%albedo_soil
      END IF
      IF (ice%Hi( vi) > 0._dp) THEN
        self%AlbedoSurf(  vi) = self%albedo_snow
      END IF
      self%Albedo( vi,:) = self%AlbedoSurf( vi)
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)
  end subroutine init

  subroutine initialise_IMAUITM_firn_from_file( self, mesh, region_name)
    ! If this is a restarted run, read the firn depth and meltpreviousyear data from the restart file

    ! In/output variables
    class(type_SMB_model_IMAU_ITM),      INTENT(INOUT) :: self
    TYPE(type_mesh),                     INTENT(IN)    :: mesh
    CHARACTER(LEN=3),                    INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_IMAUITM_firn_from_file'
    CHARACTER(LEN=256)                                 :: filename_restart_firn
    REAL(dp)                                           :: timeframe_restart_firn
    !TYPE(type_restart_data)                            :: restart

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Assume that SMB and geometry are read from the same restart file
    SELECT CASE (region_name)
    CASE('NAM')
      filename_restart_firn = C%filename_firn_IMAUITM_NAM
      timeframe_restart_firn = C%timeframe_restart_firn_IMAUITM_NAM
    CASE('EAS')
      filename_restart_firn = C%filename_firn_IMAUITM_EAS
      timeframe_restart_firn = C%timeframe_restart_firn_IMAUITM_EAS
    CASE('GRL')
      filename_restart_firn = C%filename_firn_IMAUITM_GRL
      timeframe_restart_firn = C%timeframe_restart_firn_IMAUITM_GRL
    CASE('ANT')
      filename_restart_firn = C%filename_firn_IMAUITM_ANT
      timeframe_restart_firn = C%timeframe_restart_firn_IMAUITM_ANT
    CASE DEFAULT
        CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END SELECT

     ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising SMB-model firn layer from file "' // colour_string( TRIM(filename_restart_firn),'light blue') // '"...'


    ! Read firn layer from file
    IF (timeframe_restart_firn == 1E9_dp) THEN
      ! Assume the file has no time dimension
      CALL read_field_from_file_2D_monthly( filename_restart_firn, 'FirnDepth', mesh, C%output_dir, self%FirnDepth)
      CALL read_field_from_file_2D( filename_restart_firn, 'MeltPreviousYear', mesh, C%output_dir, self%MeltPreviousYear)
    ELSE
      ! Assume the file has a time dimension, and read the specified timeframe
      CALL read_field_from_file_2D_monthly( filename_restart_firn, 'FirnDepth', mesh, C%output_dir, self%FirnDepth, time_to_read = timeframe_restart_firn)
      CALL read_field_from_file_2D( filename_restart_firn, 'MeltPreviousYear', mesh, C%output_dir, self%MeltPreviousYear, time_to_read = timeframe_restart_firn)
    END IF


    ! Finalise routine path
    CALL finalise_routine( routine_name)

  end subroutine initialise_IMAUITM_firn_from_file

  subroutine remap( self, mesh_old, mesh_new)

    ! In- and output variables
    class(type_SMB_model_IMAU_ITM), intent(inout) :: self
    type(type_mesh),                intent(in   ) :: mesh_old
    type(type_mesh),                intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_SMB_model'

    ! Add routine to path
    call init_routine( routine_name)

    call reallocate_dist_shared( self%AlbedoSurf      , self%wAlbedoSurf      , mesh_new%pai_V%n_nih    )
    call reallocate_dist_shared( self%MeltPreviousYear, self%wMeltPreviousYear, mesh_new%pai_V%n_nih    )
    call reallocate_dist_shared( self%Refreezing_year , self%wRefreezing_year , mesh_new%pai_V%n_nih    )
    call reallocate_dist_shared( self%Albedo_year     , self%wAlbedo_year     , mesh_new%pai_V%n_nih    )
    call reallocate_dist_shared( self%FirnDepth       , self%wFirnDepth       , mesh_new%pai_V%n_nih, 12)
    call reallocate_dist_shared( self%Rainfall        , self%wRainfall        , mesh_new%pai_V%n_nih, 12)
    call reallocate_dist_shared( self%Snowfall        , self%wSnowfall        , mesh_new%pai_V%n_nih, 12)
    call reallocate_dist_shared( self%AddedFirn       , self%wAddedFirn       , mesh_new%pai_V%n_nih, 12)
    call reallocate_dist_shared( self%Melt            , self%wMelt            , mesh_new%pai_V%n_nih, 12)
    call reallocate_dist_shared( self%Refreezing      , self%wRefreezing      , mesh_new%pai_V%n_nih, 12)
    call reallocate_dist_shared( self%Runoff          , self%wRunoff          , mesh_new%pai_V%n_nih, 12)
    call reallocate_dist_shared( self%Albedo          , self%wAlbedo          , mesh_new%pai_V%n_nih, 12)
    call reallocate_dist_shared( self%SMB_monthly     , self%wSMB_monthly     , mesh_new%pai_V%n_nih, 12)

    self%AlbedoSurf      ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih      ) => self%AlbedoSurf
    self%Rainfall        ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:12) => self%Rainfall
    self%Snowfall        ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:12) => self%Snowfall
    self%AddedFirn       ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:12) => self%AddedFirn
    self%Melt            ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:12) => self%Melt
    self%Refreezing      ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:12) => self%Refreezing
    self%Refreezing_year ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih      ) => self%Refreezing_year
    self%Runoff          ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:12) => self%Runoff
    self%Albedo          ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:12) => self%Albedo
    self%Albedo_year     ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih      ) => self%Albedo_year
    self%SMB_monthly     ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:12) => self%SMB_monthly
    self%FirnDepth       ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:12) => self%FirnDepth
    self%MeltPreviousYear( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih      ) => self%MeltPreviousYear

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap

end module SMB_IMAU_ITM