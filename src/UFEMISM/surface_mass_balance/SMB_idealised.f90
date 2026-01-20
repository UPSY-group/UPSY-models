module SMB_idealised

  ! Idealised SMB models

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine, colour_string
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use Halfar_SIA_solution, only: Halfar
  use SMB_basic, only: atype_SMB_model, type_SMB_model_context_allocate, type_SMB_model_context_initialise, &
    type_SMB_model_context_run, type_SMB_model_context_remap

  implicit none

  private

  public :: type_SMB_model_idealised

  type, extends(atype_SMB_model) :: type_SMB_model_idealised

    contains

      procedure, public :: allocate_SMB_model   => allocate_SMB_model_idealised
      procedure, public :: initialise_SMB_model => initialise_SMB_model_idealised
      procedure, public :: run_SMB_model        => run_SMB_model_idealised
      procedure, public :: remap_SMB_model      => remap_SMB_model_idealised

      procedure, private :: run_SMB_model_idealised_EISMINT1
      procedure, private :: run_SMB_model_idealised_Halfar_static

  end type type_SMB_model_idealised

contains

  subroutine allocate_SMB_model_idealised( self, context)

    ! In/output variables:
    class(type_SMB_model_idealised),       intent(inout) :: self
    type(type_SMB_model_context_allocate), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_SMB_model_idealised'

    ! Add routine to path
    call init_routine( routine_name)

    call self%set_name('SMB_model_idealised')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_SMB_model_idealised

  subroutine initialise_SMB_model_idealised( self, context)

    ! In/output variables:
    class(type_SMB_model_idealised),         intent(inout) :: self
    type(type_SMB_model_context_initialise), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_idealised'

    ! Add routine to path
    call init_routine( routine_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_idealised

  subroutine run_SMB_model_idealised( self, context)

    ! In/output variables:
    class(type_SMB_model_idealised),  intent(inout) :: self
    type(type_SMB_model_context_run), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_idealised'

    ! Add routine to path
    call init_routine( routine_name)

    ! Run the chosen idealised SMB model
    select case (C%choice_SMB_model_idealised)
    case default
      call crash('unknown choice_SMB_model_idealised "' // trim( C%choice_SMB_model_idealised) // '"')
    case ('EISMINT1_A', 'EISMINT1_B', 'EISMINT1_C', 'EISMINT1_D', 'EISMINT1_E', 'EISMINT1_F')
      call self%run_SMB_model_idealised_EISMINT1( self%mesh, context%time)
    case ('Halfar_static')
      call self%run_SMB_model_idealised_Halfar_static( self%mesh)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_idealised

  subroutine remap_SMB_model_idealised( self, context)

    ! In/output variables:
    class(type_SMB_model_idealised),    intent(inout) :: self
    type(type_SMB_model_context_remap), intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_SMB_model_idealised'

    ! Add routine to path
    call init_routine( routine_name)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_SMB_model_idealised

  subroutine run_SMB_model_idealised_EISMINT1( self, mesh, time)
    ! SMB for the EISMINT1 experiments (Huybrechts et al., 1996)

    ! In/output variables
    class(type_SMB_model_idealised), intent(inout) :: self
    type(type_mesh),                 intent(in   ) :: mesh
    real(dp),                        intent(in   ) :: time

    ! Local variables:
    character(len=256), parameter :: routine_name = 'run_SMB_model_idealised_EISMINT1'
    integer                       :: vi
    real(dp), parameter           :: s = 1E-2_dp    ! Mass balance change with distance from divide [m yr^-1 km^-1]
    real(dp)                      :: x, y, d, R_el, T

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_SMB_model_idealised)
    case default
      call crash('unknown choice_SMB_model_idealised "' // TRIM( C%choice_SMB_model_idealised) // '"')
    case ('EISMINT1_A')
      ! Moving margin, no cyclicity

      do vi = mesh%vi1, mesh%vi2

        ! Calculate distance from ice divide (for moving margin experiments, use Euclidean distance)
        x = mesh%V( vi,1)
        y = mesh%V( vi,2)
        d = SQRT( x**2 + y**2) / 1E3_dp  ! [km]

        ! Calculate distance from equilibrium line to ice divide
        R_el = 450._dp

        ! Calculate SMB (Huybrechts et al., Eq. 10)
        self%SMB( vi) = MIN( 0.5_dp, s * (R_el - d))

      end do

    case ('EISMINT1_B')
      ! Moving margin, 20,000-yr cyclicity

      T = 20E3_dp

      do vi = mesh%vi1, mesh%vi2

        ! Calculate distance from ice divide (for moving margin experiments, use Euclidean distance)
        x = mesh%V( vi,1)
        y = mesh%V( vi,2)
        d = SQRT( x**2 + y**2) / 1E3_dp  ! [km]

        ! Calculate distance from equilibrium line to ice divide (Huybrechts et al., Eq. 14)
        R_el = 450._dp + 100._dp * SIN( 2 * pi * time / T)

        ! Calculate SMB (Huybrechts et al., Eq. 10)
        self%SMB( vi) = MIN( 0.5_dp, s * (R_el - d))

      end do

    case ('EISMINT1_C')
      ! Moving margin, 40,000-yr cyclicity

      T = 40E3_dp

      do vi = mesh%vi1, mesh%vi2

        ! Calculate distance from ice divide (for moving margin experiments, use Euclidean distance)
        x = mesh%V( vi,1)
        y = mesh%V( vi,2)
        d = SQRT( x**2 + y**2) / 1E3_dp  ! [km]

        ! Calculate distance from equilibrium line to ice divide (Huybrechts et al., Eq. 14)
        R_el = 450._dp + 100._dp * SIN( 2._dp * pi * time / T)

        ! Calculate SMB (Huybrechts et al., Eq. 10)
        self%SMB( vi) = MIN( 0.5_dp, s * (R_el - d))

      end do

    case ('EISMINT1_D')
      ! Fixed margin, no cyclicity

      ! Calculate SMB (Huybrechts et al., Eq. 8)
      do vi = mesh%vi1, mesh%vi2
        self%SMB( vi) = 0.3_dp
      end do

    case ('EISMINT1_E')
      ! Fixed margin, 20,000-yr cyclicity

      T = 20E3_dp

      ! Calculate SMB (Huybrechts et al., Eq. 13)
      do vi = mesh%vi1, mesh%vi2
        self%SMB( vi) = 0.3_dp + 0.2_dp * SIN( 2._dp * pi * time / T)
      end do

    case ('EISMINT1_F')
      ! Fixed margin, 40,000-yr cyclicity

      T = 40E3_dp

      ! Calculate SMB (Huybrechts et al., Eq. 13)
      do vi = mesh%vi1, mesh%vi2
        self%SMB( vi) = 0.3_dp + 0.2_dp * SIN( 2._dp * pi * time / T)
      end do

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_idealised_EISMINT1

  subroutine run_SMB_model_idealised_Halfar_static( self, mesh)

    ! In/output variables
    class(type_SMB_model_idealised), intent(inout) :: self
    type(type_mesh),                 intent(in)    :: mesh

    ! Local variables:
    character(len=256), parameter :: routine_name = 'run_SMB_model_idealised_Halfar_static'
    integer                       :: vi

    ! Add routine to path
    call init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
      self%SMB( vi) = -1._dp * Halfar%dH_dt( C%uniform_Glens_flow_factor, C%Glens_flow_law_exponent, &
        C%refgeo_idealised_Halfar_H0, C%refgeo_idealised_Halfar_R0, &
        mesh%V( vi,1), mesh%V( vi,2), 0._dp)

      ! The analytical solution diverges to infinite dH/dt at the margin, limit this
      self%SMB( vi) = max( self%SMB( vi), -50._dp)
      if (sqrt( mesh%V( vi,1)**2 + mesh%V( vi,2)**2) > C%refgeo_idealised_Halfar_R0 - 1e-2_dp) then
        self%SMB( vi) = -50._dp
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_idealised_Halfar_static

end module SMB_idealised
