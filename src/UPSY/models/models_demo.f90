module models_demo

  use precisions, only: dp
  use parameters, only: pi, NaN
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use models_basic, only: atype_model
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_demo_model

  type, extends( atype_model) :: type_demo_model

    ! Some ice-model-esque data fields
    real(dp), dimension(:  ), contiguous, pointer :: H
    real(dp), dimension(:,:), contiguous, pointer :: u_3D
    real(dp), dimension(:,:), contiguous, pointer :: v_3D
    logical,  dimension(:  ), contiguous, pointer :: mask_ice
    real(dp), dimension(:,:), contiguous, pointer :: T2m
    type(MPI_WIN) :: wH, wu_3D, wv_3D, wmask_ice, wT2m

  contains

    procedure, public :: init => initialise_demo_model

  end type type_demo_model

contains

  subroutine initialise_demo_model( model, mesh, nz)

    ! In/output variables:
    class(type_demo_model), intent(  out) :: model
    type(type_mesh),        intent(in   ) :: mesh
    integer,                intent(in   ) :: nz

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_demo_model'
    integer                        :: vi,ti,k,m
    real(dp)                       :: x,y,cx,cy

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Set model metadata and grid
    call model%set_name('demo_model')
    call model%set_grid( mesh)

    ! Create model fields

    call model%create_field( model%H, model%wH, &
      mesh, Arakawa_grid%a(), &
      name      = 'H', &
      long_name = 'ice thickness', &
      units     = 'm')

    call model%create_field( model%u_3D, model%wu_3D, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = 'u_3D', &
      long_name = 'depth-dependent horizontal ice velocity in x-direction', &
      units     = 'm yr^-1')

    call model%create_field( model%v_3D, model%wv_3D, &
      mesh, Arakawa_grid%b(), third_dimension%ice_zeta( nz, 'regular'), &
      name      = 'v_3D', &
      long_name = 'depth-dependent horizontal ice velocity in y-direction', &
      units     = 'm yr^-1')

    call model%create_field( model%mask_ice, model%wmask_ice, &
      mesh, Arakawa_grid%a(), &
      name      = 'mask_ice', &
      long_name = 'ice mask', &
      units     = '-')

    call model%create_field( model%T2m, model%wT2m, &
      mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m', &
      long_name = 'Monthly 2-m air temperature', &
      units     = 'K')

    ! Initialise fields with simple analytical functions

    cx = mesh%xmax - mesh%xmin
    cy = mesh%ymax - mesh%ymin

    do vi = mesh%vi1, mesh%vi2

      x = mesh%V( vi,1)
      y = mesh%V( vi,2)

      model%H( vi) = max( 0._dp, cos( x * pi / cx) * cos( y * pi / cy) - 0.2_dp)
      model%mask_ice( vi) = model%H( vi) > 0._dp

      do m = 1, 12
        model%T2m( vi,m) = hypot( x,y) + sin( real(m,dp) * 2._dp * pi / 12._dp)
      end do

    end do

    do ti = mesh%ti1, mesh%ti2

      x = mesh%Trigc( ti,1)
      y = mesh%Trigc( ti,2)

      do k = 1, nz
        model%u_3D( ti,k) = max( 0._dp, cos( x * pi / cx) * cos( y * pi / cy) - 0.2_dp) + real( k,dp)
        model%v_3D( ti,k) =             sin( x * pi / cx) * sin( y * pi / cy)           + real( k,dp)
      end do

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_demo_model

end module models_demo
