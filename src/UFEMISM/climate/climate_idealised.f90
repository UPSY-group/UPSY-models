module climate_idealised

  ! Idealised climate models

! ===== Preamble =====
! ====================
  use precisions,                        only: dp
  use upsy_main,                         only: upsy
  use mpi_basic,                         only: par, sync
  use call_stack_and_comp_time_tracking, only: crash, init_routine, finalise_routine
  use model_configuration,               only: C
  use parameters
  use mesh_basic,                        only: type_mesh
  use ice_geometry_model_data,           only: atype_ice_geometry_model_data
  use climate_model_types,               only: type_climate_model

  implicit none

  private

  public :: type_climate_model_idealised

type, extends(atype_climate_model) :: type_climate_model_idealised
    contains

    procedure, public :: allocate_climate_model   => climate_model_idealised_allocate
    procedure, public :: deallocate_climate_model => climate_model_idealised_deallocate
    procedure, public :: initialise_climate_model => climate_model_idealised_initialise
    procedure, public :: run_climate_model        => climate_model_idealised_run
    procedure, public :: remap_climate_model      => climate_model_idealised_remap

    procedure, public :: get_climate_model_name

    procedure, private :: run_climate_model_idealised_EISMINT1

end type type_climate_model_idealised


contains
    ! ===== Main routines =====
    ! =========================

    subroutine climate_model_idealised_allocate( self)

    ! In/output variables:
    class(type_climate_model_idealised), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_idealised_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is specific to the idealised climate model

    ! Remove routine from call stack
    call finalise_routine( routine_name)

    end subroutine climate_model_idealised_allocate

    subroutine climate_model_idealised_deallocate( self)

    ! In/output variables:
    class(type_climate_model_idealised), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'climate_model_idealised_deallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Deallocate all the stuff that is specific to c.imate model prescribed


    ! Remove routine from call stack
    call finalise_routine( routine_name)

    end subroutine climate_model_idealised_deallocate

    subroutine climate_model_idealised_initialise(self, mesh, climate)
        ! Initialise the climate model
        !
        ! Use an idealised climate scheme
        ! In- and output variables
        class(type_climate_model_idealised),    intent(inout) :: self
        type(type_mesh),                        intent(in)    :: mesh
        type(type_climate_model),               intent(inout) :: climate

        ! Local variables:
        character(len=256), parameter                         :: routine_name = 'climate_model_idealised_initialise'

        ! Add routine to path
        call init_routine( routine_name)

        ! Print to terminal
        if (par%primary)  write(*,"(A)") '     Initialising idealised climate model "' // &
        UPSY%stru%colour_string( TRIM( C%choice_climate_model_idealised),'light blue') // '"...'

        ! Run the chosen idealised climate model
        if (C%choice_climate_model_idealised == 'EISMINT1_A' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_B' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_C' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_D' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_E' .OR. &
            C%choice_climate_model_idealised == 'EISMINT1_F') THEN
        ! No need to do anything
        else
        call crash('unknown choice_climate_model_idealised "' // TRIM( C%choice_climate_model_idealised) // '"')
        end if

        ! Finalise routine path
        call finalise_routine( routine_name)

    end subroutine climate_model_idealised_initialise

    subroutine climate_model_idealised_run(self, mesh, geom, climate, time)
        ! Calculate the climate
        !
        ! Use an idealised climate scheme

        ! In/output variables:
        class(type_climate_model_idealised),    intent(inout) :: self
        type(type_mesh),                        intent(in)    :: mesh
        class(atype_ice_geometry_model_data),   intent(in   ) :: geom
        type(type_climate_model),               intent(inout) :: climate
        real(dp),                               intent(in)    :: time

        ! Local variables:
        character(len=256), parameter                         :: routine_name = 'climate_model_idealised_run'

        ! Add routine to path
        call init_routine( routine_name)

        ! Run the chosen idealised climate model
        if     (C%choice_climate_model_idealised == 'EISMINT1_A' .OR. &
                C%choice_climate_model_idealised == 'EISMINT1_B' .OR. &
                C%choice_climate_model_idealised == 'EISMINT1_C' .OR. &
                C%choice_climate_model_idealised == 'EISMINT1_D' .OR. &
                C%choice_climate_model_idealised == 'EISMINT1_E' .OR. &
                C%choice_climate_model_idealised == 'EISMINT1_F') THEN
        call run_climate_model_idealised_EISMINT1( mesh, geom, climate, time)
        else
        call crash('unknown choice_climate_model_idealised "' // TRIM( C%choice_climate_model_idealised) // '"')
        end if

        ! Finalise routine path
        call finalise_routine( routine_name)

    end subroutine climate_model_idealised_run

    ! == EISMINT1
    ! ===========

    subroutine run_climate_model_idealised_EISMINT1(self, mesh, geom, climate, time)
        ! Temperature for the EISMINT1 experiments (Huybrechts et al., 1996)

        class(type_climate_model_idealised),  intent(inout) :: self
        type(type_mesh),                      intent(in)    :: mesh
        class(atype_ice_geometry_model_data), intent(in   ) :: geom
        type(type_climate_model),             intent(inout) :: climate
        real(dp),                             intent(in)    :: time

        ! Local variables:
        character(LEN=256), parameter                       :: routine_name = 'run_climate_model_idealised_EISMINT1'
        real(dp), parameter                                 :: x_summit = 0._dp      ! x-coordinate of ice divide [m]
        real(dp), parameter                                 :: y_summit = 0._dp      ! y-coordinate of ice divide [m]
        integer                                             :: vi
        real(dp)                                            :: x, y, d, T, dT

        ! Add routine to path
        call init_routine( routine_name)

        ! Set precipitation to zero - SMB is parameterised anyway...
        climate%Precip = 0._dp

        ! Baseline temperature
        if     (C%choice_climate_model_idealised == 'EISMINT1_A' .OR. &
                C%choice_climate_model_idealised == 'EISMINT1_B' .OR. &
                C%choice_climate_model_idealised == 'EISMINT1_C') THEN
        ! Moving margin: Huybrechts et al., Eq. 11

        do vi = mesh%vi1, mesh%vi2
            ! Calculate baseline temperature
            climate%T2m( vi,:) = 270._dp - 0.01_dp * geom%Hs( vi)
        end do

        elseif (C%choice_climate_model_idealised == 'EISMINT1_D' .OR. &
                C%choice_climate_model_idealised == 'EISMINT1_E' .OR. &
                C%choice_climate_model_idealised == 'EISMINT1_F') THEN
        ! Fixed margin: Huybrechts et al., Eq. 9

        do vi = mesh%vi1, mesh%vi2

            ! Calculate distance from ice divide (for fixed margin experiments, use square distance)
            x = mesh%V( vi,1)
            y = mesh%V( vi,2)
            d = MAX( ABS( x - x_summit), ABS( y - y_summit)) / 1E3_dp  ! [km]

            ! Calculate baseline temperature
            climate%T2m( vi,:) = 239._dp + (8.0E-08_dp * d**3)

        end do

        else
        call crash('unknown choice_climate_model_idealised"' // TRIM( C%choice_climate_model_idealised) // '"!')
        end if
        call sync

        ! Add temperature change for the glacial cycle experiments
        if     (C%choice_climate_model_idealised == 'EISMINT1_B' .OR. &
                C%choice_climate_model_idealised == 'EISMINT1_E') THEN
        ! 20,000-yr cyclicity

        T = 20E3_dp

        if (time > 0._dp) then
            dT = 10._dp * SIN( 2._dp * pi * time / T)
            climate%T2m = climate%T2m + dT
        end if

        elseif (C%choice_climate_model_idealised == 'EISMINT1_C' .OR. &
                C%choice_climate_model_idealised == 'EISMINT1_F') THEN
        ! 40,000-yr cyclicity

        T = 40E3_dp

        if (time > 0._dp) then
            dT = 10._dp * SIN( 2._dp * pi * time / T)
            climate%T2m = climate%T2m + dT
        end if

        end if
        call sync

        ! Finalise routine path
        call finalise_routine( routine_name)

    end subroutine run_climate_model_idealised_EISMINT1

    subroutine climate_model_idealised_remap( self, mesh_new, time, refgeo_init, refgeo_PD, geom)

        ! In/output variables:
        class(type_climate_model_idealised),   intent(inout) :: self
        type(type_mesh), target,               intent(in   ) :: mesh_new
        real(dp),                              intent(in   ) :: time
        type(type_reference_geometry), target, intent(in   ) :: refgeo_init, refgeo_PD
        class(atype_ice_geometry_model_data),  intent(in   ) :: geom

        ! Local variables:
        character(len=*), parameter :: routine_name = 'climate_model_idealised_remap'

        ! Add routine to call stack
        call init_routine( routine_name)

        ! Remap all the stuff that is specific to climate model idealised

        ! Re-initialise to read and remap the climatew from the input file again
        call self%initialise( geom, refgeo_init, refgeo_PD)

        ! Remove routine from call stack
        call finalise_routine( routine_name)

    end subroutine climate_model_idealised_remap

    function get_climate_model_name( self) result( climate_model_name)
        class(type_climate_model_idealised), intent(in) :: self
        character(len=:), allocatable :: climate_model_name
        SMB_model_name = 'idealised'
    end function get_climate_model_name

end module climate_idealised
