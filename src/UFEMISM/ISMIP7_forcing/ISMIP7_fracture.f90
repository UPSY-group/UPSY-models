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
  use ice_geometry_model_data, only: atype_ice_geometry_model_data

  implicit none

  private

  public :: type_ISMIP7_fracture_model

  type, extends(atype_model) :: type_ISMIP7_fracture_model

      character(len=:), allocatable               :: filename
      real(dp), dimension(:), allocatable         :: time_in_file                    !       List of timeframes in the NetCDF file

      type(type_grid)                             :: grid_raw                        !       The x/y-grid that the ISMIP7 folks provided the data on
      type(type_map)                              :: map                             !       Mapping object to remap data from the ISMIP7 grid to the UFEMISM mesh

      real(dp)                                    :: t0, t1                          !       Timestamps of enveloping timeframes
      real(dp), dimension(:), contiguous, pointer :: mask_dp_t0         => null()    ! [0-1] Left  enveloping timeframe
      real(dp), dimension(:), contiguous, pointer :: mask_dp_t1         => null()    ! [0-1] Right enveloping timeframe
      type(MPI_WIN) :: wmask_dp_t0, wmask_dp_t1

      real(dp), dimension(:), contiguous, pointer :: max_ice_fraction   => null()    ! [0-1] Maximum allowed ice-covered fraction
      real(dp), dimension(:), contiguous, pointer :: Hi_max             => null()    ! [m]   Maximum allowed ice thickness
      real(dp), dimension(:), contiguous, pointer :: dHi_calved         => null()    ! [m]   Ice thickness lost to calving
      type(MPI_WIN) :: wmax_ice_fraction, wHi_max, wdHi_calved

    contains

      procedure, public :: get_model_name
      procedure, public :: allocate   => ISMIP7_fracture_model_allocate
      procedure, public :: initialise => ISMIP7_fracture_model_initialise
      procedure, public :: run        => ISMIP7_fracture_model_run

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

    call self%create_field( self%max_ice_fraction, self%wmax_ice_fraction, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'max_ice_fraction', &
      long_name = 'Maximum allowed ice-covered fraction', &
      units     = '0-1', &
      remap_method = 'reallocate')

    call self%create_field( self%Hi_max, self%wHi_max, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'Hi_max', &
      long_name = 'Maximum allowed ice thickness', &
      units     = 'm', &
      remap_method = 'reallocate')

    call self%create_field( self%dHi_calved, self%wdHi_calved, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'dHi_calved', &
      long_name = 'Ice thickness lost to calving', &
      units     = 'm', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ISMIP7_fracture_model_allocate

  subroutine ISMIP7_fracture_model_initialise( self)

    ! In/output variables:
    class(type_ISMIP7_fracture_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter :: routine_name = 'ISMIP7_fracture_model_initialise'

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

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ISMIP7_fracture_model_initialise

  subroutine ISMIP7_fracture_model_run( self, time, geom, Hi)

    ! In/output variables:
    class(type_ISMIP7_fracture_model),                                  intent(inout) :: self
    real(dp),                                                           intent(in   ) :: time
    class(atype_ice_geometry_model_data),                               intent(in   ) :: geom
    real(dp), dimension(self%mesh%vi1:self%mesh%vi2), intent(inout) :: Hi
    ! real(dp), dimension(self%mesh%pai_V%i1_nih:self%mesh%pai_V%i2_nih), intent(inout) :: Hi

    ! Local variables:
    character(len=*), parameter :: routine_name = 'ISMIP7_fracture_model_run'
    real(dp)                    :: wt0, wt1
    integer                     :: vi

    ! Add routine to call stack
    call init_routine( routine_name)

    if (time < self%t0 .or. time > self%t1) then
      call self%update_timeframes( time)
    end if

    ! Calculate time interpolation weights
    wt0 = (self%t1 - time) / (self%t1 - self%t0)
    wt0 = max( 0._dp, min( 1._dp, wt0 ))
    wt1 = 1._dp - wt0

    do vi = self%mesh%vi1, self%mesh%vi2

      if (.not. geom%mask_grounded_ice( vi)) then
        ! Hydrofracturing can only occur on shelves

        ! Interpolate timeframes in time to find the maximum allowed ice fraction
        ! Note that the provided mask is true (or actually 1) when hydrofracturing
        ! occurs, so the ice fraction is the inverse of that mask.
        self%max_ice_fraction( vi) = 1._dp - (wt0 * self%mask_dp_t0( vi) + wt1 * self%mask_dp_t1( vi))

        ! Calculate maximum allowed ice thickness
        self%Hi_max( vi) = geom%Hi_eff( vi) * self%max_ice_fraction( vi)

        ! Apply maximum allowed ice thickness
        if (Hi( vi) > self%Hi_max( vi)) then
          self%dHi_calved( vi) = Hi( vi) - self%Hi_max( vi)
          Hi( vi) = self%Hi_max( vi)
        else
          self%dHi_calved( vi) = 0._dp
        end if

      else
        self%max_ice_fraction( vi) = 1._dp
        self%Hi_max          ( vi) = Hi( vi)
        self%dHi_calved      ( vi) = 0._dp
      end if

    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ISMIP7_fracture_model_run

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
    call create_map_from_xy_grid_to_mesh_vertices( self%grid_raw, self%mesh, C%output_dir, self%map, '1st_order_conservative')

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
    integer                                          :: vi

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

    ! Read the raw mask data, which is supposed to be a logical, but since NetCDF doesn't support that,
    ! is read as an 8-bit (half-precision) integer, which only has values 0 and 1.
    call read_var_primary( self%filename, ncid, id_var, mask_grid_tot_with_time, start = [1,1,ti], count = [self%grid_raw%nx, self%grid_raw%ny, 1])
    ! Convert to floating point, so we can interpolate it to the mesh. The interpolated field
    ! essentially represents 'which fraction of each Voronoi cell overlaps with grid cells of value 1'.
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

    ! Limit mask values to [0,1]
    do vi = self%mesh%vi1, self%mesh%vi2
      mask_dp( vi) = max( 0._dp, min( 1._dp, mask_dp( vi)))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_single_timeframe_from_netcdf

  function get_model_name( self) result( model_name)
    class(type_ISMIP7_fracture_model), intent(in) :: self
    character(len=:), allocatable :: model_name
    model_name = 'ISMIP7_fracture'
  end function get_model_name

end module ISMIP7_fracture
