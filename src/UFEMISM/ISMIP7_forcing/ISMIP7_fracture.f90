module ISMIP7_fracture

  ! The ISMIP7 protocol provides a single NetCDF file containing a mask
  ! for every year (1950-2299). The mask is integer-valued (as NetCDF
  ! doesn't support logicals), with 1 indicating fracture, implying
  ! that no (floating?) ice is allowed to exist there anymore.

  ! This setting is controlled by the config parameter 'do_apply_ISMIP7_fracture_mask_config'

  use precisions, only: dp
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: crash
  use models_basic, only: atype_model
  use mesh_types, only: type_mesh
  use mpi_f08, only: MPI_WIN
  use Arakawa_grid_mod, only: Arakawa_grid
  use netcdf_io_main, only: read_time_from_file, open_existing_netcdf_file_for_reading, &
    setup_xy_grid_from_file, close_netcdf_file, inquire_var, read_var_primary
  use grid_types, only: type_grid
  use remapping_types, only: type_map
  use remapping_grid_to_mesh_vertices, only: create_map_from_xy_grid_to_mesh_vertices
  use apply_maps, only: apply_map_xy_grid_to_mesh_2D
  use UPSY_main, only: UPSY
  use mpi_basic, only: par
  use mpi_distributed_memory_grid, only: distribute_gridded_data_from_primary
  use dist_to_hybrid_mod, only: dist_to_hybrid

  implicit none

  private

  public :: type_ISMIP7_fracture_model

  type, extends(atype_model) :: type_ISMIP7_fracture_model

      character(len=:), allocatable               :: filename
      real(dp), dimension(:), allocatable         :: time_in_file              ! List of timeframes in the NetCDF file

      type(type_grid)                             :: grid_raw                  ! The x/y-grid that the ISMIP7 folks provided the data on
      type(type_map)                              :: map                       ! Mapping object to remap data from the ISMIP7 grid to the UFEMISM mesh

      real(dp)                                    :: t0, t1                    ! Timestamps of enveloping timeframes
      real(dp), dimension(:), contiguous, pointer :: mask_dp_t0   => null()    ! Left  enveloping timeframe
      real(dp), dimension(:), contiguous, pointer :: mask_dp_t1   => null()    ! Right enveloping timeframe
      type(MPI_WIN) :: wmask_dp_t0, wmask_dp_t1

    contains

      procedure, public :: get_model_name
      procedure, public :: allocate   => ISMIP7_fracture_model_allocate
      procedure, public :: initialise => ISMIP7_fracture_model_initialise

      procedure, private :: initialise_remapping_object
      procedure, private :: update_timeframes
      procedure, private :: read_single_timeframe_from_netcdf

  end type type_ISMIP7_fracture_model

contains

  subroutine ISMIP7_fracture_model_allocate( self, region_name, mesh)

    ! In/output variables:
    class(type_ISMIP7_fracture_model), intent(inout) :: self
    character(len=*),                  intent(in   ) :: region_name
    type(type_mesh), target,           intent(in   ) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'ISMIP7_fracture_model_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Allocate all the stuff that is common to all models
    call self%allocate_model( region_name, mesh)

    ! Allocate all the stuff that is specific to the ISMIP7 mask-based hydrofracture model
    call self%create_field( self%mask_dp_t0, self%wmask_dp_t0, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'mask_dp_t0', &
      long_name = 'ISMIP7 hydrofracture mask at previous timeframe', &
      units     = '0-1', &
      remap_method = 'reallocate')

    call self%create_field( self%mask_dp_t1, self%wmask_dp_t1, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'mask_dp_t1', &
      long_name = 'ISMIP7 hydrofracture mask at next timeframe', &
      units     = '0-1', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ISMIP7_fracture_model_allocate

  subroutine ISMIP7_fracture_model_initialise( self)

    ! In/output variables:
    class(type_ISMIP7_fracture_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'ISMIP7_fracture_model_initialise'

    ! Add routine to call stack
    call init_routine( routine_name)

    self%filename = C%filename_ISMIP7_fracture_mask

    if (par%primary) write(0,*) '  Initialising ISMIP7 mask-based hydrofracture model from file "', &
      UPSY%stru%colour_string( trim( self%filename), 'light blue'), '"...'

    call self%initialise_remapping_object()

    call read_time_from_file( self%filename, self%time_in_file)

    self%t0 = C%start_time_of_run - 2._dp
    self%t1 = C%start_time_of_run - 1._dp
    call self%update_timeframes( C%start_time_of_run)

    call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ISMIP7_fracture_model_initialise

  subroutine initialise_remapping_object( self)

    ! In/output variables
    class(type_ISMIP7_fracture_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'initialise_remapping_object'
    integer                       :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! Read the grid from the first listed input file
    call open_existing_netcdf_file_for_reading( self%filename, ncid)
    call setup_xy_grid_from_file( self%filename, ncid, self%grid_raw)
    call close_netcdf_file( ncid)

    ! Calculate the remapping operator
    self%grid_raw%name = 'ISMIP7_input_grid_' // trim( self%name())
    call create_map_from_xy_grid_to_mesh_vertices( self%grid_raw, self%mesh, C%output_dir, self%map, '2nd_order_conservative')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_remapping_object

  subroutine update_timeframes( self, time)

    ! In/output variables:
    class(type_ISMIP7_fracture_model), intent(inout) :: self
    real(dp),                          intent(in   ) :: time

    ! Local variables:
    character(len=*), parameter :: routine_name = 'update_timeframes'
    integer                     :: ti0, ti1

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Determine which timeframes to read
    ti0 = 1
    ti1 = 2
    do while (self%time_in_file( ti1) < time .and. ti1 < size( self%time_in_file,1))
      ti0 = ti1
      ti1 = ti1 + 1
    end do

    self%t0 = self%time_in_file( ti0)
    self%t1 = self%time_in_file( ti1)

    call self%read_single_timeframe_from_netcdf( ti0, self%mask_dp_t0)
    call self%read_single_timeframe_from_netcdf( ti1, self%mask_dp_t1)

    call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine update_timeframes

  subroutine read_single_timeframe_from_netcdf( self, ti, mask_dp)

    ! In/output variables:
    class(type_ISMIP7_fracture_model), intent(in   ) :: self
    integer,                           intent(in   ) :: ti
    real(dp), dimension(:),            intent(inout) :: mask_dp

    ! Local variables
    character(len=*), parameter                      :: routine_name = 'read_single_timeframe_from_netcdf'
    integer,  dimension(:,:,:), allocatable          :: mask_grid_tot_with_time
    real(dp), dimension(:,:  ), allocatable          :: mask_grid_dp_tot
    integer                                          :: ncid, id_var
    real(dp), dimension(:    ), allocatable          :: mask_grid_dp_vec_partial
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2) :: mask_dp_dist

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) then
      write(0,'(A,I4,A,A,A)') '   Reading ISMIP7 hydrofracture mask for year ', floor( self%time_in_file( ti)), ' from file "', &
        UPSY%stru%colour_string( trim( self%filename), 'light blue'), '"...'
    end if

    ! Read raw gridded data to the primary
    if (par%primary) then
      allocate( mask_grid_tot_with_time( self%grid_raw%nx, self%grid_raw%ny, 1))
      allocate( mask_grid_dp_tot       ( self%grid_raw%nx, self%grid_raw%ny   ))
    else
      allocate( mask_grid_tot_with_time(0,0,0))
      allocate( mask_grid_dp_tot       (0,0  ))
    end if

    call open_existing_netcdf_file_for_reading( self%filename, ncid)
    call inquire_var( self%filename, ncid, 'mask', id_var)
    call read_var_primary( self%filename, ncid, id_var, mask_grid_tot_with_time)
    if (par%primary) mask_grid_dp_tot = real( mask_grid_tot_with_time( :,:,1), dp)
    deallocate( mask_grid_tot_with_time)
    call close_netcdf_file( ncid)

    ! Distribute gridded data to the processes
    allocate( mask_grid_dp_vec_partial( self%grid_raw%pai%i1: self%grid_raw%pai%i2))
    call distribute_gridded_data_from_primary( self%grid_raw, mask_grid_dp_vec_partial, mask_grid_dp_tot)
    deallocate( mask_grid_dp_tot)

    ! Remap data to mesh
    call apply_map_xy_grid_to_mesh_2D( self%grid_raw, self%mesh, self%map, mask_grid_dp_vec_partial, mask_dp_dist)
    call dist_to_hybrid( self%mesh%pai_V, mask_dp_dist, mask_dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_single_timeframe_from_netcdf

  function get_model_name( self) result( model_name)
    class(type_ISMIP7_fracture_model), intent(in) :: self
    character(len=:), allocatable :: model_name
    model_name = 'ISMIP7_fracture'
  end function get_model_name

end module ISMIP7_fracture
