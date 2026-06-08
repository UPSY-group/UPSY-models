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
    setup_xy_grid_from_file, close_netcdf_file
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

      procedure(initialise_ISMIP7_forcing_field_ifc ), deferred :: initialise
      procedure(update_ISMIP7_forcing_field_ifc     ), deferred :: update
      procedure(interpolate_ISMIP7_forcing_field_ifc), deferred :: interpolate

      procedure, private :: initialise_remapping_object

    end type atype_ISMIP7_forcing_field

    ! Abstract interfaces for deferred procedures
    ! ===========================================

    abstract interface

      subroutine initialise_ISMIP7_forcing_field_ifc( self, ISMIP7_forcing_foldername, ISMIP7_forcing_version, mesh, field_name)
        import atype_ISMIP7_forcing_field, type_mesh
        class(atype_ISMIP7_forcing_field), intent(inout) :: self
        character(len=*),                  intent(in   ) :: ISMIP7_forcing_foldername
        character(len=*),                  intent(in   ) :: ISMIP7_forcing_version
        type(type_mesh),                   intent(in   ) :: mesh
        character(len=*),                  intent(in   ) :: field_name
      end subroutine initialise_ISMIP7_forcing_field_ifc

      subroutine update_ISMIP7_forcing_field_ifc( self, mesh, time)
        import atype_ISMIP7_forcing_field, type_mesh, dp
        class(atype_ISMIP7_forcing_field), intent(inout) :: self
        type(type_mesh),                   intent(in   ) :: mesh
        real(dp),                          intent(in   ) :: time
      end subroutine update_ISMIP7_forcing_field_ifc

      subroutine interpolate_ISMIP7_forcing_field_ifc( self, mesh, time)
        import atype_ISMIP7_forcing_field, type_mesh, dp
        class(atype_ISMIP7_forcing_field), intent(inout) :: self
        type(type_mesh),                   intent(in   ) :: mesh
        real(dp),                          intent(in   ) :: time
      end subroutine interpolate_ISMIP7_forcing_field_ifc

    end interface

    ! Concrete types for monthly and yearly ISMIP7 forcing fields
    ! ===========================================================

    type, extends(atype_ISMIP7_forcing_field) :: type_ISMIP7_forcing_field_monthly
      !< Two enveloping timeframes and time-interpolated values of a single monthly ISMIP7 forcing field

      real(dp), dimension(:,:), allocatable :: val0          !< Values of timeframe before current time
      real(dp), dimension(:,:), allocatable :: val1          !< Values of timeframe after current time
      real(dp), dimension(:,:), allocatable :: val_interp    !< Interpolated values

    contains

      procedure :: initialise  => initialise_monthly
      procedure :: update      => update_timeframes_monthly
      procedure :: interpolate => interpolate_single_field_monthly

    end type type_ISMIP7_forcing_field_monthly

    type, extends(atype_ISMIP7_forcing_field) :: type_ISMIP7_forcing_field_yearly
      !< Two enveloping timeframes and time-interpolated values of a single yearly ISMIP7 forcing field

      real(dp), dimension(:), allocatable :: val0          !< Values of timeframe before current time
      real(dp), dimension(:), allocatable :: val1          !< Values of timeframe after current time
      real(dp), dimension(:), allocatable :: val_interp    !< Interpolated values

    contains

      procedure :: initialise  => initialise_yearly
      procedure :: update      => update_timeframes_yearly
      procedure :: interpolate => interpolate_single_field_yearly

    end type type_ISMIP7_forcing_field_yearly

  contains

    subroutine initialise_monthly( self, ISMIP7_forcing_foldername, ISMIP7_forcing_version, mesh, field_name)

      ! In/output variables:
      class(type_ISMIP7_forcing_field_monthly), intent(inout) :: self
      character(len=*),                         intent(in   ) :: ISMIP7_forcing_foldername
      character(len=*),                         intent(in   ) :: ISMIP7_forcing_version
      type(type_mesh),                          intent(in   ) :: mesh
      character(len=*),                         intent(in   ) :: field_name

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'initialise_monthly'
      character(len=1024)            :: filename

      ! Add routine to call stack
      call init_routine( routine_name)

      ! Define name
      self%name = field_name

      ! Get info from files
      call gather_fileinfo( ISMIP7_forcing_foldername, ISMIP7_forcing_version, &
        self%filenames, self%timestamps, self%name)

      call self%initialise_remapping_object( mesh)

      ! Deallocate if necessary
      if (allocated( self%val0       )) deallocate( self%val0      )
      if (allocated( self%val1       )) deallocate( self%val1      )
      if (allocated( self%val_interp )) deallocate( self%val_interp)

      ! Allocate memory for timeframes
      allocate (self%val0      ( mesh%vi1:mesh%vi2, 12), source = NaN)
      allocate (self%val1      ( mesh%vi1:mesh%vi2, 12), source = NaN)
      allocate (self%val_interp( mesh%vi1:mesh%vi2, 12), source = NaN)

      ! Update timeframes to the current model time
      call self%update( mesh, C%start_time_of_run)

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine initialise_monthly

    subroutine initialise_yearly( self, ISMIP7_forcing_foldername, ISMIP7_forcing_version, mesh, field_name)

      ! In/output variables:
      class(type_ISMIP7_forcing_field_yearly), intent(inout) :: self
      character(len=*),                        intent(in   ) :: ISMIP7_forcing_foldername
      character(len=*),                        intent(in   ) :: ISMIP7_forcing_version
      type(type_mesh),                         intent(in   ) :: mesh
      character(len=*),                        intent(in   ) :: field_name

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'initialise_yearly'
      character(len=1024)            :: filename

      ! Add routine to call stack
      call init_routine( routine_name)

      ! Define name
      self%name = field_name

      ! Get info from files
      call gather_fileinfo( ISMIP7_forcing_foldername, ISMIP7_forcing_version, &
        self%filenames, self%timestamps, self%name)

      call self%initialise_remapping_object( mesh)

      ! Deallocate if necessary
      if (allocated( self%val0       )) deallocate( self%val0      )
      if (allocated( self%val1       )) deallocate( self%val1      )
      if (allocated( self%val_interp )) deallocate( self%val_interp)

      ! Allocate memory for timeframes
      allocate (self%val0      ( mesh%vi1:mesh%vi2), source = NaN)
      allocate (self%val1      ( mesh%vi1:mesh%vi2), source = NaN)
      allocate (self%val_interp( mesh%vi1:mesh%vi2), source = NaN)

      ! Update timeframes to the current model time
      call self%update( mesh, C%start_time_of_run)

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine initialise_yearly

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

    subroutine update_timeframes_monthly( self, mesh, time)

      ! In/output variables:
      class(type_ISMIP7_forcing_field_monthly), intent(inout) :: self
      type(type_mesh),                          intent(in   ) :: mesh
      real(dp),                                 intent(in   ) :: time

      ! Local variables
      character(len=*), parameter :: routine_name = 'update_timeframes_monthly'
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
        call update_single_timeframe_monthly( mesh, self, self%ti0, self%val0)
      end if

      if (self%ti1 /= ti1_old) then
        call update_single_timeframe_monthly( mesh, self, self%ti1, self%val1)
      end if

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine update_timeframes_monthly

    subroutine update_timeframes_yearly( self, mesh, time)

      ! In/output variables:
      class(type_ISMIP7_forcing_field_yearly), intent(inout) :: self
      type(type_mesh),                         intent(in   ) :: mesh
      real(dp),                                intent(in   ) :: time

      ! Local variables
      character(len=*), parameter :: routine_name = 'update_timeframes_yearly'
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
        call update_single_timeframe_yearly( mesh, self, self%ti0, self%val0)
      end if

      if (self%ti1 /= ti1_old) then
        call update_single_timeframe_yearly( mesh, self, self%ti1, self%val1)
      end if

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine update_timeframes_yearly

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

    subroutine update_single_timeframe_monthly( mesh, field, ti, val)

      ! In/output variables:
      type(type_mesh),                               intent(in   ) :: mesh
      type(type_ISMIP7_forcing_field_monthly),       intent(in   ) :: field
      integer,                                       intent(in   ) :: ti
      real(dp), dimension( mesh%vi1:mesh%vi2, 1:12), intent(inout) :: val

      ! Local variables
      character(len=*), parameter             :: routine_name = 'update_single_timeframe_monthly'
      character(len=:), allocatable           :: filename
      real(dp), dimension(:), allocatable     :: time_from_file
      integer                                 :: mi
      real(dp), dimension( mesh%vi1:mesh%vi2) :: d_month

      ! Add routine to path
      call init_routine( routine_name)

      filename = trim(field%filenames( ti))

      ! Read time dimension from file
      call read_time_from_file( filename, time_from_file)
      if (size( time_from_file,1) /= 12) call crash('file "' // trim( filename) // '" doesnt have 12 months')

      if (par%primary) then
        write(0,*) '   Reading ISMIP7 monthly climate forcing from file: ', &
          UPSY%stru%colour_string( trim( filename), 'light blue')
      end if

      ! Read all 12 months individually with extrapolation
      do mi = 1, 12
        call read_field_from_file_2D( filename, trim(field%name), mesh, C%output_dir, d_month, &
          time_to_read = time_from_file( mi), extrapolate_fillvalues = .true.)
        val( mesh%vi1:mesh%vi2, mi) = d_month( mesh%vi1:mesh%vi2)
      end do

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine update_single_timeframe_monthly

    subroutine update_single_timeframe_yearly( mesh, field, ti, val)

      ! In/output variables:
      type(type_mesh),                          intent(in   ) :: mesh
      type(type_ISMIP7_forcing_field_yearly),   intent(in   ) :: field
      integer,                                  intent(in   ) :: ti
      real(dp), dimension( mesh%vi1:mesh%vi2),  intent(inout) :: val

      ! Local variables
      character(len=*), parameter             :: routine_name = 'update_single_timeframe_yearly'
      character(len=:), allocatable           :: filename
      real(dp), dimension(:), allocatable     :: time_from_file

      ! Add routine to path
      call init_routine( routine_name)

      filename = trim(field%filenames( ti))

      ! Read time dimension from file
      call read_time_from_file( filename, time_from_file)
      if (size( time_from_file,1) /= 1) call crash('file "' // trim( filename) // '" doesnt have 1 timeframe')

      if (par%primary) then
        write(0,*) '   Reading ISMIP7 yearly climate forcing from file: ', &
          UPSY%stru%colour_string( trim( filename), 'light blue')
      end if

      ! Read with extrapolation
      call read_field_from_file_2D( filename, trim(field%name), mesh, C%output_dir, val, &
        time_to_read = time_from_file( 1), extrapolate_fillvalues = .true.)

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine update_single_timeframe_yearly

    subroutine interpolate_single_field_monthly( self, mesh, time)

      ! In/output variables:
      class(type_ISMIP7_forcing_field_monthly), intent(inout) :: self
      type(type_mesh),                          intent(in   ) :: mesh
      real(dp),                                 intent(in   ) :: time

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'interpolate_single_field_monthly'
      real(dp)                       :: w0, w1
      real(dp)                       :: days

      ! Add routine to call stack
      call init_routine( routine_name)

      call get_interpolation_weights( self%timestamps( self%ti0), self%timestamps( self%ti1), time, w0, w1)

      ! Apply interpolation to values
      self%val_interp( mesh%vi1:mesh%vi2, :) = w0 * self%val0( mesh%vi1:mesh%vi2, :) + w1 * self%val1( mesh%vi1:mesh%vi2, :)

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine interpolate_single_field_monthly

    subroutine interpolate_single_field_yearly( self, mesh, time)

      ! In/output variables:
      class(type_ISMIP7_forcing_field_yearly), intent(inout) :: self
      type(type_mesh),                         intent(in   ) :: mesh
      real(dp),                                intent(in   ) :: time

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'interpolate_single_field_yearly'
      real(dp)                       :: w0, w1
      real(dp)                       :: days

      ! Add routine to call stack
      call init_routine( routine_name)

      call get_interpolation_weights( self%timestamps( self%ti0), self%timestamps( self%ti1), time, w0, w1)

      ! Apply interpolation to values
      self%val_interp( mesh%vi1:mesh%vi2) = w0 * self%val0( mesh%vi1:mesh%vi2) + w1 * self%val1( mesh%vi1:mesh%vi2)

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine interpolate_single_field_yearly

    subroutine get_interpolation_weights( timestamp0, timestamp1, time, w0, w1)

      ! In/output variables:
      real(dp), intent(in   ) :: timestamp0
      real(dp), intent(in   ) :: timestamp1
      real(dp), intent(in   ) :: time
      real(dp), intent(  out) :: w0
      real(dp), intent(  out) :: w1

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'get_interpolation_weights'

      ! Add routine to call stack
      call init_routine( routine_name)

      ! Get weights
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

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine get_interpolation_weights

  end module ISMIP7_forcing_field_types