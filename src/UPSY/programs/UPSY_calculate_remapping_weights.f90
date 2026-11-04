program UPSY_calculate_remapping_weights

  ! A program for calculating remapping weights between two
  ! meshes/grids. Call syntax:
  !
  ! mpiexec -n [n_cores] UPSY_calculate_remapping_weights [filename_grid_src] [Arakawa_grid_src] [filename_grid_dst] [Arakawa_grid_dst] [method] [filename_output]
  !
  ! n_cores          : number of CPU cores you want to run on
  ! filename_grid_src: path/to/source/grid_or_mesh.nc
  ! Arakawa_grid_src : Arakawa grid of source data (a, b, c, cx, cy)
  ! filename_grid_dst: path/to/destination/grid_or_mesh.nc
  ! Arakawa_grid_dst : Arakawa grid of destination data (a, b, c, cx, cy)
  ! method           : nearest_neighbour, linear, first_order_conservative, second_order_conservative
  ! filename_output  : path/to/output/file/containing/remapping_weights.nc

#include <petsc/finclude/petscsys.h>

  use basic_program_info, only: program_name
  use mpi_basic, only: par, sync
  use petsc, only: PetscErrorF, PetscFinalize
  use mpi_basic, only: initialise_parallelisation
  use parameters, only: initialise_constants
  use call_stack_and_comp_time_tracking, only: initialise_control_and_resource_tracker, init_routine, finalise_routine
  use crash_mod, only: crash
  use UPSY_main, only: UPSY
  use netcdf_io_main, only: inquire_xy_grid, inquire_lonlat_grid, inquire_mesh, &
    open_existing_netcdf_file_for_reading, setup_xy_grid_from_file, setup_lonlat_grid_from_file, &
    setup_mesh_from_file, close_netcdf_file
  use basic_grid_mod, only: atype_grid
  use xy_grid_basic, only: type_xy_grid
  use lonlat_grid_basic, only: type_lonlat_grid
  use mesh_basic, only: type_mesh
  use remapping_basic, only: type_map
  use Arakawa_grid_mod, only: Arakawa_grid

  implicit none

  integer             :: ierr
  character(len=1024) :: filename_grid_src, Arakawa_grid_src_str
  character(len=1024) :: filename_grid_dst, Arakawa_grid_dst_str
  character(len=1024) :: remapping_method_str
  character(len=1024) :: filename_output
  logical             :: file_exists
  class(atype_grid), allocatable :: grid_src, grid_dst
  type(type_map)      :: map

  program_name = 'UPSY_calculate_remapping_weights'

  ! Initialise MPI parallelisation and PETSc
  call initialise_parallelisation

  ! Initialise constants (pi, NaN, ...)
  call initialise_constants

  ! Initialise the control and resource tracker
  call initialise_control_and_resource_tracker

  ! Get the input arguments
  call get_input_arguments( filename_grid_src, Arakawa_grid_src_str, &
    filename_grid_dst, Arakawa_grid_dst_str, remapping_method_str, filename_output)

  ! Find out what kind of grids we're remapping between
  call read_grid( filename_grid_src, grid_src)
  call read_grid( filename_grid_dst, grid_dst)

  ! Set up the remapping context
  call map%ct%create( &
    grid_src, Arakawa_grid%create_from_string( Arakawa_grid_src_str), &
    grid_dst, Arakawa_grid%create_from_string( Arakawa_grid_dst_str), &
    UPSY%remapping%method%create_from_string( remapping_method_str))

  ! Create the remapping object
  call map%create()

  ! Write it to a NetCDF file
  call write_map_to_file( map, filename_output)

  ! Finalise PETSc and MPI parallelisation
  PetscCall( PetscFinalize( ierr))

contains

  subroutine get_input_arguments( filename_grid_src, Arakawa_grid_src_str, &
    filename_grid_dst, Arakawa_grid_dst_str, remapping_method_str, filename_output)

    ! In/output variables:
    character(len=*), intent(out) :: filename_grid_src
    character(len=*), intent(out) :: Arakawa_grid_src_str
    character(len=*), intent(out) :: filename_grid_dst
    character(len=*), intent(out) :: Arakawa_grid_dst_str
    character(len=*), intent(out) :: remapping_method_str
    character(len=*), intent(out) :: filename_output

    ! Local variables:
    character(len=*), parameter :: routine_name = 'get_input_arguments'

    ! Add routine to path
    call init_routine( routine_name)

    call safety_check_input_arguments()

    call getarg( 1, filename_grid_src   )
    call getarg( 2, Arakawa_grid_src_str)
    call getarg( 3, filename_grid_dst   )
    call getarg( 4, Arakawa_grid_dst_str)
    call getarg( 5, remapping_method_str)
    call getarg( 6, filename_output     )

    if (par%primary) then
      write(0,*) UPSY%stru%colour_string( trim( program_name) // ': calculating remapping weights.', 'green')
      write(0,*) '  Source      mesh grid   : ', UPSY%stru%colour_string( trim( filename_grid_src   ), 'light blue')
      write(0,*) '  Source      Arakawa grid: ', UPSY%stru%colour_string( trim( Arakawa_grid_src_str), 'light blue')
      write(0,*) '  Destination mesh grid   : ', UPSY%stru%colour_string( trim( filename_grid_dst   ), 'light blue')
      write(0,*) '  Destination Arakawa grid: ', UPSY%stru%colour_string( trim( Arakawa_grid_dst_str), 'light blue')
      write(0,*) '  Method                  : ', UPSY%stru%colour_string( trim( remapping_method_str), 'light blue')
      write(0,*) '  Output file             : ', UPSY%stru%colour_string( trim( filename_output     ), 'light blue')
    end if

    ! Check if the input files actually exist
    inquire( exist = file_exists, file = trim( filename_grid_src))
    if (.not. file_exists) call crash('file ' // trim( filename_grid_src) // ' not found')
    inquire( exist = file_exists, file = trim( filename_grid_dst))
    if (.not. file_exists) call crash('file ' // trim( filename_grid_dst) // ' not found')
    inquire( exist = file_exists, file = trim( filename_output))
    if (      file_exists) call crash('file ' // trim( filename_output) // ' already exists')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine get_input_arguments

  subroutine safety_check_input_arguments()

    ! Local variables:
    character(len=*), parameter :: routine_name = 'safety_check_input_arguments'

    ! Add routine to path
    call init_routine( routine_name)

    if (iargc() /= 6) then
      if (par%primary) then
        write(0,*) UPSY%stru%colour_string( trim( program_name) // 'Error: invalid call syntax for ' // trim( program_name) // '. Should be:', 'yellow')
        write(0,*) UPSY%stru%colour_string( '  mpiexec -n [n_cores] UPSY_calculate_remapping_weights [filename_grid_src] [Arakawa_grid_src] [filename_grid_dst] [Arakawa_grid_dst] [method] [filename_output]', 'yellow')
        write(0,*) UPSY%stru%colour_string( '    n_cores          : number of CPU cores you want to run on', 'yellow')
        write(0,*) UPSY%stru%colour_string( '    filename_grid_src: path/to/source/grid_or_mesh.nc', 'yellow')
        write(0,*) UPSY%stru%colour_string( '    Arakawa_grid_src : Arakawa grid of source data (a, b, c, cx, cy)', 'yellow')
        write(0,*) UPSY%stru%colour_string( '    filename_grid_dst: path/to/destination/grid_or_mesh.nc', 'yellow')
        write(0,*) UPSY%stru%colour_string( '    Arakawa_grid_dst : Arakawa grid of destination data (a, b, c, cx, cy)', 'yellow')
        write(0,*) UPSY%stru%colour_string( '    method           : nearest_neighbour, linear, first_order_conservative, second_order_conservative', 'yellow')
        write(0,*) UPSY%stru%colour_string( '    filename_output  : path/to/output/file/containing/remapping_weights.nc', 'yellow')
      end if
      call sync
      error stop
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine safety_check_input_arguments

  subroutine read_grid( filename, grid)

    ! In/output variables:
    character(len=*),               intent(in   ) :: filename
    class(atype_grid), allocatable, intent(inout) :: grid

    ! Local variables:
    character(len=*), parameter :: routine_name = 'read_grid'
    logical                     :: has_xy_grid, has_lonlat_grid, has_mesh
    integer                     :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! Inquire what kind of grid the file is defined on
    call inquire_xy_grid(     filename, has_xy_grid    )
    call inquire_lonlat_grid( filename, has_lonlat_grid)
    call inquire_mesh(        filename, has_mesh       )

    ! Safety
    if (has_xy_grid     .and. has_lonlat_grid) call crash('file ' // trim( filename) // ' contains both an x/y-grid and a lon/lat-grid')
    if (has_xy_grid     .and. has_mesh       ) call crash('file ' // trim( filename) // ' contains both an x/y-grid and a mesh')
    if (has_lonlat_grid .and. has_mesh       ) call crash('file ' // trim( filename) // ' contains both a lon/lat-grid and a mesh')

    ! Read the grid from the file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    if (has_xy_grid) then
      allocate( type_xy_grid :: grid)
      select type(xy_grid => grid)
      class is (type_xy_grid)
        call setup_xy_grid_from_file( filename, ncid, xy_grid)
      end select
    elseif (has_lonlat_grid) then
      allocate( type_lonlat_grid :: grid)
      select type(lonlat_grid => grid)
      class is (type_lonlat_grid)
        call setup_lonlat_grid_from_file( filename, ncid, lonlat_grid)
      end select
    elseif (has_mesh) then
      allocate( type_mesh :: grid)
      select type(mesh => grid)
      class is (type_mesh)
        call setup_mesh_from_file( filename, ncid, mesh)
      end select
    else
      call crash('could not find a valid grid in file ' // trim( filename))
    end if
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_grid

  subroutine write_map_to_file( map, filename_output)

    ! In/output variables
    type(type_map),   intent(in) :: map
    character(len=*), intent(in) :: filename_output

    ! Local variables
    character(len=*), parameter :: routine_name = 'write_map_to_file'

    ! Add routine to path
    call init_routine( routine_name)

    call crash('fixme!')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_map_to_file

end program UPSY_calculate_remapping_weights