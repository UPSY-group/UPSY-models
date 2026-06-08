module ISMIP7_forcing_field_types

  use mpi_basic, only: par
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use basic_model_utilities, only: list_files_in_folder
  use model_configuration, only: C
  use crash_mod, only: crash
  use UPSY_main, only: UPSY
  use mesh_types, only: type_mesh
  use netcdf_io_main, only: read_time_from_file, read_field_from_file_2D, open_existing_netcdf_file_for_reading, &
    setup_xy_grid_from_file, close_netcdf_file, inquire_fill_value, inquire_var, read_var_primary
  use mpi_distributed_memory_grid, only: distribute_gridded_data_from_primary
  use smooth_gridded_data, only: extrapolate_fillvalue_Gaussian_grid
  use apply_maps, only: apply_map_xy_grid_to_mesh_2D, apply_map_xy_grid_to_mesh_3D
  use remapping_grid_to_mesh_vertices, only: create_map_from_xy_grid_to_mesh_vertices
  use parameters, only: NaN
  use grid_types, only: type_grid
  use remapping_types, only: type_map

  implicit none

  private

  public :: type_ISMIP7_forcing_field_monthly, type_ISMIP7_forcing_field_yearly

    ! Abstract type for monthly/yearly ISMIP7 forcing fields
    ! ======================================================

    type, abstract :: atype_ISMIP7_forcing_field
      !< Metadata of monthly/yearly ISMIP7 forcing fields

      character(len=1024)                            :: name          !           'tas', 'pr', 'acabf', 'dacabfdz', etc
      character(len=1024)                            :: foldername    !           Foldername that contains all files
      character(len=1024), dimension(:), allocatable :: filenames     !           Filenames

      real(dp), dimension(:), allocatable            :: timestamps    ! [years]   All time values in combined files

      integer                                        :: ti0 = -1      !           Index of timeframe before current time
      integer                                        :: ti1 = -1      !           Index of timeframe after current time

      type(type_grid)                                :: grid_raw      !           The x/y-grid that the ISMIP7 folks provided the data on
      type(type_map)                                 :: map           !           Mapping object to remap data from the ISMIP7 grid to the UFEMISM mesh

    contains

      procedure,                                       public  :: initialise                   => initialise_ISMIP7_forcing_field
      procedure,                                       public  :: update_and_interpolate       => update_and_interpolate_ISMIP7_forcing_field

      procedure,                                       private :: initialise_remapping_object
      procedure(allocate_timeframes_ifc),    deferred, private :: allocate_timeframes
      procedure(interpolate_timeframes_ifc), deferred, private :: interpolate_timeframes

    end type atype_ISMIP7_forcing_field

    ! Abstract interfaces for deferred procedures
    ! ===========================================

    abstract interface

      subroutine allocate_timeframes_ifc( self, mesh)
        import atype_ISMIP7_forcing_field, type_mesh
        class(atype_ISMIP7_forcing_field), intent(inout) :: self
        type(type_mesh),                   intent(in   ) :: mesh
      end subroutine allocate_timeframes_ifc

      subroutine interpolate_timeframes_ifc( self, mesh, time)
        import atype_ISMIP7_forcing_field, type_mesh, dp
        class(atype_ISMIP7_forcing_field), intent(inout) :: self
        type(type_mesh),                   intent(in   ) :: mesh
        real(dp),                          intent(in   ) :: time
      end subroutine interpolate_timeframes_ifc

    end interface

    ! Concrete types for monthly and yearly ISMIP7 forcing fields
    ! ===========================================================

    type, extends(atype_ISMIP7_forcing_field) :: type_ISMIP7_forcing_field_monthly
      !< Two enveloping timeframes and time-interpolated values of a single monthly ISMIP7 forcing field

      real(dp), dimension(:,:), allocatable :: val0          !< Values of timeframe before current time
      real(dp), dimension(:,:), allocatable :: val1          !< Values of timeframe after current time
      real(dp), dimension(:,:), allocatable :: val_interp    !< Interpolated values

    contains

      procedure, public  :: allocate_timeframes    => allocate_timeframes_monthly
      procedure, private :: interpolate_timeframes => interpolate_timeframes_monthly

    end type type_ISMIP7_forcing_field_monthly

    type, extends(atype_ISMIP7_forcing_field) :: type_ISMIP7_forcing_field_yearly
      !< Two enveloping timeframes and time-interpolated values of a single yearly ISMIP7 forcing field

      real(dp), dimension(:), allocatable :: val0          !< Values of timeframe before current time
      real(dp), dimension(:), allocatable :: val1          !< Values of timeframe after current time
      real(dp), dimension(:), allocatable :: val_interp    !< Interpolated values

    contains

      procedure, public  :: allocate_timeframes    => allocate_timeframes_yearly
      procedure, private :: interpolate_timeframes => interpolate_timeframes_yearly

    end type type_ISMIP7_forcing_field_yearly

  contains

    subroutine initialise_ISMIP7_forcing_field( self, ISMIP7_forcing_foldername, ISMIP7_forcing_version, mesh, field_name)

      ! In/output variables:
      class(atype_ISMIP7_forcing_field), intent(inout) :: self
      character(len=*),                  intent(in   ) :: ISMIP7_forcing_foldername
      character(len=*),                  intent(in   ) :: ISMIP7_forcing_version
      type(type_mesh),                   intent(in   ) :: mesh
      character(len=*),                  intent(in   ) :: field_name

      ! Local variables:
      character(len=*), parameter   :: routine_name = 'initialise_ISMIP7_forcing_field'
      character(len=:), allocatable :: filename

      ! Add routine to call stack
      call init_routine( routine_name)

      ! Define name
      self%name = field_name

      ! Get info from files
      call gather_fileinfo( ISMIP7_forcing_foldername, ISMIP7_forcing_version, &
        self%filenames, self%timestamps, self%name)

      ! Initialise the single remapping object that will be re-used for all
      ! NetCDF files of this field (which presumably are all defined on the same x/y-grid)
      call self%initialise_remapping_object( mesh)

      ! Allocate memory for the before, after, and interpolated timeframes
      call self%allocate_timeframes( mesh)

      ! Update and interpolate timeframes to the current model time
      call self%update_and_interpolate( mesh, C%start_time_of_run)

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine initialise_ISMIP7_forcing_field

    subroutine allocate_timeframes_monthly( self, mesh)

      ! In/output variables:
      class(type_ISMIP7_forcing_field_monthly), intent(inout) :: self
      type(type_mesh),                          intent(in   ) :: mesh

      ! Local variables:
      character(len=*), parameter :: routine_name = 'allocate_timeframes_monthly'

      ! Add routine to call stack
      call init_routine( routine_name)

      ! Deallocate if necessary
      if (allocated( self%val0       )) deallocate( self%val0      )
      if (allocated( self%val1       )) deallocate( self%val1      )
      if (allocated( self%val_interp )) deallocate( self%val_interp)

      ! Allocate memory for timeframes
      allocate (self%val0      ( mesh%vi1:mesh%vi2, 12), source = NaN)
      allocate (self%val1      ( mesh%vi1:mesh%vi2, 12), source = NaN)
      allocate (self%val_interp( mesh%vi1:mesh%vi2, 12), source = NaN)

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine allocate_timeframes_monthly

    subroutine allocate_timeframes_yearly( self, mesh)

      ! In/output variables:
      class(type_ISMIP7_forcing_field_yearly), intent(inout) :: self
      type(type_mesh),                         intent(in   ) :: mesh

      ! Local variables:
      character(len=*), parameter :: routine_name = 'allocate_timeframes_yearly'

      ! Add routine to call stack
      call init_routine( routine_name)

      ! Deallocate if necessary
      if (allocated( self%val0       )) deallocate( self%val0      )
      if (allocated( self%val1       )) deallocate( self%val1      )
      if (allocated( self%val_interp )) deallocate( self%val_interp)

      ! Allocate memory for timeframes
      allocate (self%val0      ( mesh%vi1:mesh%vi2), source = NaN)
      allocate (self%val1      ( mesh%vi1:mesh%vi2), source = NaN)
      allocate (self%val_interp( mesh%vi1:mesh%vi2), source = NaN)

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine allocate_timeframes_yearly

    subroutine initialise_remapping_object( self, mesh)

      ! In/output variables
      class(atype_ISMIP7_forcing_field), intent(inout) :: self
      type(type_mesh),                   intent(in   ) :: mesh

      ! Local variables:
      character(len=*), parameter   :: routine_name = 'initialise_remapping_object'
      character(len=:), allocatable :: filename
      integer                       :: ncid

      ! Add routine to path
      call init_routine( routine_name)

      ! Read the grid from the first listed input file
      filename = trim( self%filenames(1))
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call setup_xy_grid_from_file( filename, ncid, self%grid_raw)
      call close_netcdf_file( ncid)

      ! Calculate the remapping operator
      self%grid_raw%name = 'ISMIP7_input_grid_' // trim( self%name)
      call create_map_from_xy_grid_to_mesh_vertices( self%grid_raw, mesh, C%output_dir, self%map, '2nd_order_conservative')

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine initialise_remapping_object

    subroutine gather_fileinfo( ISMIP7_forcing_foldername, ISMIP7_forcing_version, filenames, timestamps, var_name)

      ! In/output variables
      character(len=*),                               intent(in   ) :: ISMIP7_forcing_foldername
      character(len=*),                               intent(in   ) :: ISMIP7_forcing_version
      character(len=1024), dimension(:), allocatable, intent(inout) :: filenames
      real(dp),            dimension(:), allocatable, intent(inout) :: timestamps
      character(len=*),                               intent(in   ) :: var_name

      ! Local variables:
      character(len=*), parameter                    :: routine_name = 'gather_fileinfo'
      character(len=:), allocatable                  :: foldername
      integer                                        :: i
      real(dp)                                       :: year

      ! Add routine to path
      call init_routine( routine_name)

      if (allocated( filenames )) deallocate( filenames)
      if (allocated( timestamps)) deallocate( timestamps)

      ! Construct foldernames
      foldername = trim( ISMIP7_forcing_foldername) // '/' // trim(var_name) // '/' // trim( ISMIP7_forcing_version)

      call list_files_in_folder( foldername, filenames, trim(var_name))
      if (size( filenames,1) == 0) call crash('could not find any valid NetCDF files in directory "' // trim( foldername) // '"')

      ! Read timestamps
      allocate( timestamps( size( filenames,1)))
      do i = 1, size( filenames,1)
        call read_year_from_netcdf_filename( filenames( i), year)
        ! Add half a year, so that we have the timestamp at the middle of the year rather than the start
        timestamps( i) = year + 0.5_dp
      end do

      ! Append foldername to filenames
      do i = 1, size( filenames,1)
        filenames( i) = trim( foldername) // '/' // trim( filenames( i))
      end do

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine gather_fileinfo

    subroutine read_year_from_netcdf_filename( filename, year)
      ! Rather than trying to read from the file's time dimension, which would mean
      ! dealing with whatever calendar it uses, just read it from the filename

      ! In/output variables
      character(len=*), intent(in   ) :: filename
      real(dp),         intent(  out) :: year

      ! Local variables
      character(len=*), parameter :: routine_name = 'read_year_from_netcdf_filename'
      character(len=4) :: year_str
      integer          :: year_int, stat

      ! Add routine to path
      call init_routine( routine_name)

      ! Since the files are called e.g. 'acabf_AIS_CESM2-WACCM_ssp585_SDBN1-8000m_v2_2015.nc',
      ! just assumed that the last few characters spell out the year

      year_str = filename( len_trim( filename) - 6 : len_trim( filename) - 3)
      year_int = UPSY%stru%str2int( year_str, stat)
      if (stat /= 0) call crash('could not read year from filename "' // trim( filename) // '"')

      year = real( year_int,dp)

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine read_year_from_netcdf_filename

    subroutine update_and_interpolate_ISMIP7_forcing_field( self, mesh, time)

      ! In/output variables:
      class(atype_ISMIP7_forcing_field), intent(inout) :: self
      type(type_mesh),                   intent(in   ) :: mesh
      real(dp),                          intent(in   ) :: time

      ! Local variables
      character(len=*), parameter :: routine_name = 'update_and_interpolate_ISMIP7_forcing_field'
      integer                     :: ti0_old, ti1_old

      ! Add routine to path
      call init_routine( routine_name)

      ! Get current bracket indices
      ti0_old = self%ti0
      ti1_old = self%ti1

      ! Update the indices of time slices before and after current time
      call update_bracket_indices( self%timestamps, self%ti0, self%ti1, time)

      ! Update timeframes if necessary
      if (self%ti0 /= ti0_old) then
        select type (f => self)
          class default
            call crash('invalid extentsion of atype_ISMIP7_forcing_field')
          class is (type_ISMIP7_forcing_field_monthly)
            call read_single_timeframe_from_netcdf_monthly( f, mesh, f%ti0, f%val0)
          class is (type_ISMIP7_forcing_field_yearly)
            call read_single_timeframe_from_netcdf_yearly( f, mesh, f%ti0, f%val0)
        end select
      end if

      if (self%ti1 /= ti1_old) then
        select type (f => self)
          class default
            call crash('invalid extentsion of atype_ISMIP7_forcing_field')
          class is (type_ISMIP7_forcing_field_monthly)
            call read_single_timeframe_from_netcdf_monthly( f, mesh, f%ti1, f%val1)
          class is (type_ISMIP7_forcing_field_yearly)
            call read_single_timeframe_from_netcdf_yearly( f, mesh, f%ti1, f%val1)
        end select
      end if

      ! Interpolate between timeframes
      call self%interpolate_timeframes( mesh, time)

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine update_and_interpolate_ISMIP7_forcing_field

    subroutine update_bracket_indices( timestamps, ti0, ti1, time)

      ! In/output variables:
      real(dp), dimension(:), intent(in   ) :: timestamps
      integer,                intent(inout) :: ti0
      integer,                intent(inout) :: ti1
      real(dp),               intent(in   ) :: time

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'update_bracket_indices'
      integer                        :: i, n

      ! Add routine to call stack
      call init_routine( routine_name)

      ! Determine total numer of timeframes available for this field
      n = size(timestamps)

      if (time <= timestamps(1)) then
        ! Model time before first available time value, return first two indices
        ti0 = 1
        ti1 = 2

      elseif (time >= timestamps(n)) then
        ! Model time after last available time value, return last two indices
        ti0 = n-1
        ti1 = n

      else
        ! Model time within array, return bracketing indices
        do i = 1, n
          if (timestamps(i) <= time) then
            ti0 = i
            ti1 = i+1
          end if
        end do

      end if

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine update_bracket_indices

    subroutine read_single_timeframe_from_netcdf_monthly( self, mesh, ti, val)

      ! In/output variables:
      type(type_ISMIP7_forcing_field_monthly),       intent(in   ) :: self
      type(type_mesh),                               intent(in   ) :: mesh
      integer,                                       intent(in   ) :: ti
      real(dp), dimension( mesh%vi1:mesh%vi2, 1:12), intent(inout) :: val

      ! Local variables
      character(len=*), parameter             :: routine_name = 'update_single_timeframe_monthly'
      character(len=:), allocatable           :: filename
      real(dp), dimension(:,:,:), allocatable :: d_grid_tot
      integer                                 :: ncid, id_var
      real(dp)                                :: fill_value
      real(dp), dimension(:,:  ), allocatable :: d_grid_vec_partial
      real(dp)                                :: sigma

      ! Add routine to path
      call init_routine( routine_name)

      filename = trim( self%filenames( ti))

      if (par%primary) then
        write(0,*) '   Reading ISMIP7 monthly forcing from file: ', &
          UPSY%stru%colour_string( trim( filename), 'light blue')
      end if

      ! Read raw gridded data to the primary
      if (par%primary) then
        allocate( d_grid_tot( self%grid_raw%nx, self%grid_raw%ny, 12), source = NaN)
      else
        allocate( d_grid_tot(0,0,0))
      end if

      call open_existing_netcdf_file_for_reading( filename, ncid)
      call inquire_fill_value( filename, ncid, self%name, fill_value)
      call inquire_var( filename, ncid, self%name, id_var)
      call read_var_primary( filename, ncid, id_var, d_grid_tot)
      call close_netcdf_file( ncid)

      ! Distribute gridded data to the processes
      allocate( d_grid_vec_partial( self%grid_raw%pai%i1: self%grid_raw%pai%i2, 12))
      call distribute_gridded_data_from_primary( self%grid_raw, d_grid_vec_partial, d_grid_tot)
      deallocate( d_grid_tot)

      ! Extrapolate data into fill_value cells
      sigma = self%grid_raw%dx * 2._dp
      call extrapolate_fillvalue_Gaussian_grid( self%grid_raw, d_grid_vec_partial, sigma, fill_value)

      ! Remap data to mesh
      call apply_map_xy_grid_to_mesh_3D( self%grid_raw, mesh, self%map, d_grid_vec_partial, val)

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine read_single_timeframe_from_netcdf_monthly

    subroutine read_single_timeframe_from_netcdf_yearly( self, mesh, ti, val)

      ! In/output variables:
      type(type_ISMIP7_forcing_field_yearly),  intent(in   ) :: self
      type(type_mesh),                         intent(in   ) :: mesh
      integer,                                 intent(in   ) :: ti
      real(dp), dimension( mesh%vi1:mesh%vi2), intent(inout) :: val

      ! Local variables
      character(len=*), parameter             :: routine_name = 'update_single_timeframe_yearly'
      character(len=:), allocatable           :: filename
      real(dp), dimension(:,:,:), allocatable :: d_grid_tot_with_time
      real(dp), dimension(:,:  ), allocatable :: d_grid_tot
      integer                                 :: ncid, id_var
      real(dp)                                :: fill_value
      real(dp), dimension(:    ), allocatable :: d_grid_vec_partial
      real(dp)                                :: sigma

      ! Add routine to path
      call init_routine( routine_name)

      filename = trim(self%filenames( ti))

      if (par%primary) then
        write(0,*) '   Reading ISMIP7 yearly climate forcing from file: ', &
          UPSY%stru%colour_string( trim( filename), 'light blue')
      end if

      ! Read raw gridded data to the primary
      if (par%primary) then
        allocate( d_grid_tot_with_time( self%grid_raw%nx, self%grid_raw%ny, 1), source = NaN)
        allocate( d_grid_tot          ( self%grid_raw%nx, self%grid_raw%ny   ), source = NaN)
      else
        allocate( d_grid_tot_with_time(0,0,0))
        allocate( d_grid_tot          (0,0  ))
      end if

      call open_existing_netcdf_file_for_reading( filename, ncid)
      call inquire_fill_value( filename, ncid, self%name, fill_value)
      call inquire_var( filename, ncid, self%name, id_var)
      call read_var_primary( filename, ncid, id_var, d_grid_tot_with_time)
      if (par%primary) d_grid_tot = d_grid_tot_with_time( :,:,1)
      deallocate( d_grid_tot_with_time)
      call close_netcdf_file( ncid)

      ! Distribute gridded data to the processes
      allocate( d_grid_vec_partial( self%grid_raw%pai%i1: self%grid_raw%pai%i2))
      call distribute_gridded_data_from_primary( self%grid_raw, d_grid_vec_partial, d_grid_tot)
      deallocate( d_grid_tot)

      ! Extrapolate data into fill_value cells
      sigma = self%grid_raw%dx * 2._dp
      call extrapolate_fillvalue_Gaussian_grid( self%grid_raw, d_grid_vec_partial, sigma, fill_value)

      ! Remap data to mesh
      call apply_map_xy_grid_to_mesh_2D( self%grid_raw, mesh, self%map, d_grid_vec_partial, val)

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine read_single_timeframe_from_netcdf_yearly

    subroutine interpolate_timeframes_monthly( self, mesh, time)
      class(type_ISMIP7_forcing_field_monthly), intent(inout) :: self
      type(type_mesh),                          intent(in   ) :: mesh
      real(dp),                                 intent(in   ) :: time
      real(dp) :: w0, w1
      call calc_interpolation_weights( self%timestamps( self%ti0), self%timestamps( self%ti1), time, w0, w1)
      self%val_interp( mesh%vi1:mesh%vi2, :) = w0 * self%val0( mesh%vi1:mesh%vi2, :) + w1 * self%val1( mesh%vi1:mesh%vi2, :)
    end subroutine interpolate_timeframes_monthly

    subroutine interpolate_timeframes_yearly( self, mesh, time)
      class(type_ISMIP7_forcing_field_yearly), intent(inout) :: self
      type(type_mesh),                         intent(in   ) :: mesh
      real(dp),                                intent(in   ) :: time
      real(dp) :: w0, w1
      call calc_interpolation_weights( self%timestamps( self%ti0), self%timestamps( self%ti1), time, w0, w1)
      self%val_interp( mesh%vi1:mesh%vi2) = w0 * self%val0( mesh%vi1:mesh%vi2) + w1 * self%val1( mesh%vi1:mesh%vi2)
    end subroutine interpolate_timeframes_yearly

    subroutine calc_interpolation_weights( timestamp0, timestamp1, time, w0, w1)
      real(dp), intent(in   ) :: timestamp0
      real(dp), intent(in   ) :: timestamp1
      real(dp), intent(in   ) :: time
      real(dp), intent(  out) :: w0
      real(dp), intent(  out) :: w1
      ! Calculate linear interpolation weights
      if (time < timestamp0) then
        ! Model time before bracket times, put full weight on the first timeframe
        w0 = 1._dp
      elseif (time > timestamp1) then
        ! Model time after bracket times, put full weight on the last timeframe
        w0 = 0._dp
      else
        ! Model time between bracket times, determine interpolated weight
        w0 = (timestamp1 - time) / (timestamp1 - timestamp0)
      end if
      w1 = 1._dp - w0
    end subroutine calc_interpolation_weights

  end module ISMIP7_forcing_field_types