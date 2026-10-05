module TIPMIP_forcing_field_types

  use mpi_basic, only: par
  use mpi_f08, only: MPI_WIN, MPI_BCAST, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use basic_model_utilities, only: list_files_in_folder
  use model_configuration, only: C
  use crash_mod, only: crash, warning
  use UPSY_main, only: UPSY
  use mesh_types, only: type_mesh
  use netcdf_io_main
  use mpi_distributed_memory_grid, only: distribute_gridded_data_from_primary
  use smooth_gridded_data, only: extrapolate_fillvalue_Gaussian_grid
  use apply_maps, only: apply_map_xy_grid_to_mesh_2D, apply_map_xy_grid_to_mesh_3D
  use remapping_grid_to_mesh_vertices, only: create_map_from_xy_grid_to_mesh_vertices
  use parameters, only: NaN, freshwater_density, sec_per_year, ice_density
  use grid_types, only: type_grid
  use remapping_types, only: type_map
  use dist_to_hybrid_mod, only: dist_to_hybrid
  use models_basic, only: atype_model
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use calendar, only: convert_month_to_days, convert_time_to_days

  implicit none

  private

  public :: type_TIPMIP_forcing_field

    ! Type for TIPMIP forcing fields
    ! ======================================================

    type :: type_TIPMIP_forcing_field
      !< Metadata of TIPMIP forcing fields

      character(len=1024)                            :: name          !           'tas_anomaly', 'pr_ratio'
      character(len=1024)                            :: foldername    !           Foldername that contains all files
      character(len=1024), dimension(:), allocatable :: filenames     !           Filenames

      real(dp), dimension(:), allocatable            :: timestamps    ! [years]   All year values in combined files
      integer, dimension(:), allocatable             :: fileindices   !           Index of file for each year

      real(dp)                                       :: y0            !           Year before current time
      real(dp)                                       :: y1            !           Year after current time

      type(type_grid)                                :: grid_raw      !           The x/y-grid that the TIPMIP folks provided the data on
      type(type_map)                                 :: map           !           Mapping object to remap data from the TIPMIP grid to the UFEMISM mesh

    !< Variables and functions that are specific to the TIPMIP forcing field
      real(dp), dimension(:,:), contiguous, pointer :: val0       => null()   !< Values of timeframe before current time
      real(dp), dimension(:,:), contiguous, pointer :: val1       => null()   !< Values of timeframe after current time
      real(dp), dimension(:,:), contiguous, pointer :: val_interp => null()   !< Time-interpolated values
      type(MPI_WIN) :: wval0, wval1, wval_interp

    contains

      procedure, public  :: allocate => allocate_TIPMIP_forcing_field
      procedure, public  :: initialise => initialise_TIPMIP_forcing_field
      procedure, public  :: update_and_interpolate

      procedure, private :: initialise_remapping_object
      procedure, private :: interpolate_timeframes
      procedure, private :: gather_fileinfo
      procedure, private :: read_single_timeframe_from_netcdf
      procedure, private :: update_bracket_years

    end type type_TIPMIP_forcing_field

  contains

    subroutine allocate_TIPMIP_forcing_field( self, model, name, long_name, units)

      ! In/output variables:
      class(type_TIPMIP_forcing_field), intent(inout) :: self
      class(atype_model),              intent(inout) :: model
      character(len=*),                intent(in   ) :: name, long_name, units

      ! Local variables:
      character(len=*), parameter   :: routine_name = 'allocate_TIPMIP_forcing_field'

      ! Add routine to call stack
      call init_routine( routine_name)

      self%name = name

      ! Create model fields for all three timeframes
      call model%create_field( self%val0, self%wval0, &
        model%mesh, Arakawa_grid%a(), third_dimension%month(), &
        name      = trim( name) // '_val0', &
        long_name = trim( long_name) // ' - timeframe 0', &
        units     = trim( units), &
        remap_method = 'reallocate')

      call model%create_field( self%val1, self%wval1, &
        model%mesh, Arakawa_grid%a(), third_dimension%month(), &
        name      = trim( name) // '_val1', &
        long_name = trim( long_name) // ' - timeframe 1', &
        units     = trim( units), &
        remap_method = 'reallocate')

      call model%create_field( self%val_interp, self%wval_interp, &
        model%mesh, Arakawa_grid%a(), third_dimension%month(), &
        name      = trim( name) // '_val_interp', &
        long_name = trim( long_name) // ' - time-interpolated', &
        units     = trim( units), &
        remap_method = 'reallocate')

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine allocate_TIPMIP_forcing_field

    subroutine initialise_TIPMIP_forcing_field( self, TIPMIP_forcing_foldername, mesh)

      ! In/output variables:
      class(type_TIPMIP_forcing_field), intent(inout) :: self
      character(len=*),                intent(in   ) :: TIPMIP_forcing_foldername
      type(type_mesh),                 intent(in   ) :: mesh

      ! Local variables:
      character(len=*), parameter   :: routine_name = 'initialise_TIPMIP_forcing_field'
      character(len=:), allocatable :: filename

      ! Add routine to call stack
      call init_routine( routine_name)

      self%foldername = trim( TIPMIP_forcing_foldername) // '/' // trim(self%name)

      ! Get info from files
      call self%gather_fileinfo()

      ! Initialise the single remapping object that will be re-used for all
      ! NetCDF files of this field (which presumably are all defined on the same x/y-grid)
      call self%initialise_remapping_object( mesh)

      ! Update and interpolate timeframes to the current model time
      call self%update_and_interpolate( mesh, C%start_time_of_run)

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine initialise_TIPMIP_forcing_field

    subroutine initialise_remapping_object( self, mesh)

      ! In/output variables
      class(type_TIPMIP_forcing_field), intent(inout) :: self
      type(type_mesh),                 intent(in   ) :: mesh

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
      self%grid_raw%name = 'TIPMIP_input_grid_' // trim( self%name)
      call create_map_from_xy_grid_to_mesh_vertices( self%grid_raw, mesh, C%output_dir, self%map, '2nd_order_conservative')

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine initialise_remapping_object

    subroutine gather_fileinfo( self)

      ! In/output variables
      class(type_TIPMIP_forcing_field), intent(inout) :: self

      ! Local variables:
      character(len=*), parameter          :: routine_name = 'gather_fileinfo'
      character(len=1024)                  :: filename
      integer                              :: i, ncid, id_dim_time, nt, id_var_time
      integer                              :: ierr
      real(dp), dimension(:), allocatable  :: time_tmp
      real(dp)                             :: y_start, y_end
      integer                              :: cnt, t
      real(dp), dimension(1000)            :: time_buff
      integer, dimension(1000)             :: fi_buff

      ! Add routine to path
      call init_routine( routine_name)

      if (allocated( self%filenames  )) deallocate( self%filenames)
      if (allocated( self%timestamps )) deallocate( self%timestamps)
      if (allocated( self%fileindices)) deallocate( self%fileindices)

      ! Extract filenames
      call list_files_in_folder( trim(self%foldername), self%filenames, trim(self%name))
      if (size( self%filenames,1) == 0) call crash('could not find any valid NetCDF files in directory "' // trim( self%foldername) // '"')

      ! Get years contained in combined files, and associated fileindices

      ! Initialise counter for timeframes
      cnt = 1

      do i = 1, size(self%filenames)

        ! Construct the full filename
        filename = trim(self%foldername) // '/' // trim(self%filenames( i))

        ! Open file and extract time variable
        call open_existing_netcdf_file_for_reading( filename, ncid)
        call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)
        call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

        ! Copy time variable into a temporary array and close netcdf
        allocate( time_tmp( nt))
        call read_var_primary( filename, ncid, id_var_time, time_tmp)
        call close_netcdf_file( ncid)

        ! Copy temporary time array to all processes
        call MPI_BCAST( time_tmp(:), nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        ! Get first and last years. By disallowing a residual and converting to integer, these result in full years
        call convert_time_to_days( time_tmp(1) , y_start, calendar='noleap', refyear=0, allow_residual=.false.)
        call convert_time_to_days( time_tmp(nt), y_end  , calendar='noleap', refyear=0, allow_residual=.false.)

        deallocate( time_tmp)

        do t = int(y_start), int(y_end)
          ! Add time value to buffer
          time_buff( cnt) = t
          ! Add file index to buffer
          fi_buff( cnt) = i
          ! Increase count
          cnt = cnt + 1
        end do

      end do

      ! Allocate arrays for timestamps and file indices, combined for all files in this folder
      allocate( self%timestamps ( cnt-1))
      allocate( self%fileindices( cnt-1))
      
      ! Copy all valid values from buffer
      self%timestamps  = time_buff( 1:cnt-1)
      self%fileindices = fi_buff  ( 1:cnt-1)

      ! TODO write out to check

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine gather_fileinfo

    subroutine update_and_interpolate( self, mesh, time)

      ! In/output variables:
      class(type_TIPMIP_forcing_field), intent(inout) :: self
      type(type_mesh),                 intent(in   ) :: mesh
      real(dp),                        intent(in   ) :: time

      ! Local variables
      character(len=*), parameter :: routine_name = 'update_and_interpolate'
      real(dp)                    :: y0_old, y1_old
      real(dp), parameter         :: eps = 1e-8_dp

      ! Add routine to path
      call init_routine( routine_name)

      ! Get current bracket years
      y0_old = self%y0
      y1_old = self%y1

      ! Update the calendar years before and after the current time
      call self%update_bracket_years( time)

      ! Update timeframes if necessary
      if (self%y0 /= y0_old) then
        if (abs(self%y0 - y1_old) < eps) then
          ! Copy data from old timeframe 1
          self%val0( mesh%vi1:mesh%vi2,:) = self%val1( mesh%vi1:mesh%vi2,:)
        else
          ! Read new timeframe from NetCDF
          call read_single_timeframe_from_netcdf( self, mesh, self%y0, self%val0)
        end if
      end if

      if (self%y1 /= y1_old) then
        call read_single_timeframe_from_netcdf( self, mesh, self%y1, self%val1)
      end if

      ! Interpolate between timeframes
      call self%interpolate_timeframes( mesh, time)

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine update_and_interpolate

    subroutine update_bracket_years( self, time)

      ! In/output variables:
      class(type_TIPMIP_forcing_field), intent(inout) :: self
      real(dp),                        intent(in   ) :: time

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'update_bracket_years'
      integer                        :: i, n

      ! Add routine to call stack
      call init_routine( routine_name)

      ! Determine total numer of timeframes available for this field
      n = size(self%timestamps)

      if (time <= self%timestamps(1)) then
        ! Model time before first available time value, return first two values
        self%y0 = self%timestamps(1)
        self%y1 = self%timestamps(2)

      elseif (time >= self%timestamps(n)) then
        ! Model time after last available time value, return last two indices
        self%y0 = self%timestamps(n-1)
        self%y1 = self%timestamps(n)

      else
        ! Extract bracket years directly from time
        self%y0 = floor(time - 0.5_dp)
        self%y1 = self%y0 + 1._dp

      end if

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine update_bracket_years

    subroutine read_single_timeframe_from_netcdf( self, mesh, y, val)

      ! In/output variables:
      class(type_TIPMIP_forcing_field), intent(in   ) :: self
      type(type_mesh),                  intent(in   ) :: mesh
      real(dp),                         intent(in   ) :: y
      real(dp), dimension(:,:),         intent(inout) :: val

      ! Local variables
      character(len=*), parameter               :: routine_name = 'read_single_timeframe_from_netcdf'
      character(len=:), allocatable             :: filename
      real(dp), dimension(:,:,:), allocatable   :: d_grid_tot
      integer                                   :: ncid, id_var
      real(dp)                                  :: fill_value
      real(dp), dimension(:,:,:), allocatable   :: d_grid_with_time
      real(dp), dimension(:,:  ), allocatable   :: d_grid_vec_partial
      real(dp)                                  :: sigma
      real(dp), dimension(mesh%vi1:mesh%vi2,12) :: d_dist
      real(dp), parameter                       :: eps = 1e-8_dp
      integer                                   :: i, idx, m, ti
      real(dp)                                  :: time_to_read

      ! Add routine to path
      call init_routine( routine_name)

      ! Get index and filename containing requested year y
      idx = -1
      do i = 1, size(self%timestamps)
        if (abs( y - self%timestamps( i)) < eps) then
          idx = self%fileindices( i)
        end if
      end do

      filename = trim(self%foldername) // '/' // trim(self%filenames( idx))

      if (par%primary) then
        write(0,*) '   Reading TIPMIP forcing from file: ', &
          UPSY%stru%colour_string( trim( filename), 'light blue')
      end if

      ! Read raw gridded data to the primary
      if (par%primary) then
        allocate( d_grid_with_time( self%grid_raw%nx, self%grid_raw%ny, 1))
        allocate( d_grid_tot( self%grid_raw%nx, self%grid_raw%ny, 12), source = NaN)
      else
        allocate( d_grid_with_time(0,0,0))
        allocate( d_grid_tot(0,0,0))
      end if

      ! Read data from file
      call open_existing_netcdf_file_for_reading( filename, ncid)
      call inquire_fill_value( filename, ncid, self%name, fill_value)
      call inquire_var( filename, ncid, self%name, id_var)
      do m = 1, 12
        call convert_month_to_days( int(y), m, time_to_read)
        call find_timeframe( filename, ncid, time_to_read, ti)
        call read_var_primary( filename, ncid, id_var, d_grid_with_time, start = (/1, 1, ti /), &
          count = (/ self%grid_raw%nx, self%grid_raw%ny, 1/))
        if (par%primary) d_grid_tot( :, :, m) = d_grid_with_time( :, :, 1)
      end do
      call close_netcdf_file( ncid)

      deallocate( d_grid_with_time)

      ! Distribute gridded data to the processes
      allocate( d_grid_vec_partial( self%grid_raw%pai%i1: self%grid_raw%pai%i2, 12))
      call distribute_gridded_data_from_primary( self%grid_raw, d_grid_vec_partial, d_grid_tot)
      deallocate( d_grid_tot)

      ! Extrapolate data into fill_value cells
      sigma = self%grid_raw%dx * 2._dp
      call extrapolate_fillvalue_Gaussian_grid( self%grid_raw, d_grid_vec_partial, sigma, fill_value)

      ! Remap data to mesh
      call apply_map_xy_grid_to_mesh_3D( self%grid_raw, mesh, self%map, d_grid_vec_partial, d_dist)
      call dist_to_hybrid( mesh%pai_V, 12, d_dist, val)

      ! Finalise routine path
      call finalise_routine( routine_name)

    end subroutine read_single_timeframe_from_netcdf

    subroutine interpolate_timeframes( self, mesh, time)
      class(type_TIPMIP_forcing_field),  intent(inout) :: self
      type(type_mesh),                  intent(in   ) :: mesh
      real(dp),                         intent(in   ) :: time
      real(dp) :: w0, w1

      ! Calculate interpolation based on center years (hence the +0.5)
      call calc_interpolation_weights( self%y0 + 0.5_dp, self%y1 + 0.5_dp, time, w0, w1)
      self%val_interp( mesh%vi1:mesh%vi2, :) = w0 * self%val0( mesh%vi1:mesh%vi2, :) + w1 * self%val1( mesh%vi1:mesh%vi2, :)
    end subroutine interpolate_timeframes

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

  end module TIPMIP_forcing_field_types
