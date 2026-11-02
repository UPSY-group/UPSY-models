module climate_retreat_mask
  !< Prescribed ice-shelf retreat from static or time-dependent NetCDF masks
  !
  ! The mask file contains a variable "mask" with values in [0,1], on an x/y-grid,
  ! a lon/lat-grid, or a mesh. A time-dependent mask has a time axis in years.
  ! Between two frames the mask is interpolated linearly in time; outside the
  ! time axis the nearest frame is held. Retreat is applied wherever the mask
  ! exceeds retreat_mask_threshold, either as removal of floating ice
  ! (shelf_collapse_type = 'calving') or as prescribed shelf melt
  ! (shelf_collapse_type = 'BMB'). Missing values are not allowed in the file.

  use mpi_f08, only: MPI_BCAST, MPI_INTEGER, MPI_DOUBLE_PRECISION, MPI_CHARACTER, MPI_COMM_WORLD, &
    MPI_ALLREDUCE, MPI_IN_PLACE, MPI_SUM
  use netcdf, only: NF90_NOERR, NF90_CHAR, NF90_GLOBAL, NF90_INQUIRE_ATTRIBUTE, NF90_GET_ATT, NF90_MAX_VAR_DIMS
  use precisions, only: dp
  use mpi_basic, only: par
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash, warning
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_geometry_model_data, only: atype_ice_geometry_model_data
  use climate_model_types, only: type_climate_model
  use netcdf_io_main
  use remapping_main, only: map_from_mesh_to_mesh_2D
  use apply_maps, only: clear_all_maps_involving_this_mesh
  use checksum_mod, only: checksum
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

  implicit none

  private

  public :: initialise_climate_retreat_mask, run_climate_retreat_mask, remap_climate_retreat_mask
  public :: retreat_mask_threshold, retreat_mask_BMB_shelf, resolve_region_filename

  ! Retreat is applied where the (time-interpolated) mask exceeds this value
  real(dp), parameter :: retreat_mask_threshold = 0.01_dp
  ! [m/yr] Shelf melt prescribed where the mask is active (shelf_collapse_type = 'BMB')
  real(dp), parameter :: retreat_mask_BMB_shelf = -400._dp

contains

  subroutine initialise_climate_retreat_mask( mesh, geom, climate, region_name)
    !< Validate the mask file and set up the retreat mask on the model mesh

    ! In/output variables:
    type(type_mesh),                      intent(in   ) :: mesh
    class(atype_ice_geometry_model_data), intent(in   ) :: geom
    type(type_climate_model),             intent(inout) :: climate
    character(len=3),                     intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_climate_retreat_mask'
    character(len=1024)            :: filename
    integer                        :: n

    ! Add routine to path
    call init_routine( routine_name)

    filename = C%ISMIP_future_shelf_collapse_forcing_filename

    select case (C%shelf_collapse_type)
    case default
      call crash('unknown shelf_collapse_type "' // trim( C%shelf_collapse_type) // '"')
    case ('calving', 'BMB')
    end select
    if (C%retreat_mask_applied_only_to_open_ocean .and. .not. C%retreat_mask_without_time) &
      call crash('retreat_mask_applied_only_to_open_ocean requires a static retreat mask')

    if (par%primary) write(0,'(A)') '     Initialising retreat mask from "' // trim( filename) // '"...'

    ! Reject missing or out-of-range values before they reach the spatial remapping
    call check_retreat_mask_values( filename)

    ! Time axis of a time-dependent mask
    if (.not. C%retreat_mask_without_time) then
      if (allocated( climate%retreat%times)) deallocate( climate%retreat%times)
      call read_time_from_file( filename, climate%retreat%times)
      n = size( climate%retreat%times)
      if (n < 1) call crash('retreat mask file "' // trim( filename) // '" has an empty time axis')
      if (any( .not. ieee_is_finite( climate%retreat%times))) &
        call crash('retreat mask file "' // trim( filename) // '" has non-finite times')
      if (n > 1) then
        if (any( climate%retreat%times( 2:n) <= climate%retreat%times( 1:n-1))) &
          call crash('retreat mask times in file "' // trim( filename) // '" must be strictly increasing')
      end if
      call check_retreat_mask_time_units( filename)
      call report_retreat_mask_time_axis( climate%retreat%times)
    end if

    ! Allocate memory
    call allocate_retreat_mask( mesh, climate)
    climate%retreat%region_name = region_name

    ! Open-ocean reference, fixed for the entire experiment
    if (C%retreat_mask_applied_only_to_open_ocean) then
      if (C%retreat_mask_open_ocean_reference_filename /= '') then
        if (count( [C%do_NAM, C%do_EAS, C%do_GRL, C%do_ANT]) > 1 .and. &
          index( C%retreat_mask_open_ocean_reference_filename, '{region}') == 0) &
          call crash('with more than one model region, retreat_mask_open_ocean_reference_filename_config ' // &
            'must contain the token {region}')
        climate%retreat%open_ocean_reference_filename = &
          resolve_region_filename( C%retreat_mask_open_ocean_reference_filename, region_name)
        call read_open_ocean_reference( mesh, climate)
      else
        climate%retreat%open_ocean( mesh%vi1:mesh%vi2) = geom%mask_icefree_ocean( mesh%vi1:mesh%vi2)
        climate%retreat%open_ocean_reference_filename = trim( C%output_dir)
        n = len_trim( climate%retreat%open_ocean_reference_filename)
        if (n > 0) then
          if (climate%retreat%open_ocean_reference_filename( n:n) /= '/') &
            climate%retreat%open_ocean_reference_filename = trim( climate%retreat%open_ocean_reference_filename) // '/'
        end if
        climate%retreat%open_ocean_reference_filename = trim( climate%retreat%open_ocean_reference_filename) // &
          'retreat_mask_open_ocean_reference_' // region_name // '.nc'
        call write_open_ocean_reference( mesh, climate)
      end if
      if (par%primary) write(0,'(A)') '     Retreat-mask open-ocean reference: "' // &
        trim( climate%retreat%open_ocean_reference_filename) // '"'
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_climate_retreat_mask

  subroutine run_climate_retreat_mask( mesh, climate, time)
    !< Update the retreat mask to the current model time

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    type(type_climate_model), intent(inout) :: climate
    real(dp),                 intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_climate_retreat_mask'
    character(len=1024)            :: filename
    integer                        :: i0, i1, n, vi
    real(dp)                       :: t, w1

    ! Add routine to path
    call init_routine( routine_name)

    filename = C%ISMIP_future_shelf_collapse_forcing_filename

    if (C%retreat_mask_without_time) then
      ! Static mask: read once per mesh

      if (.not. climate%retreat%loaded) then
        call read_field_from_file_2D( filename, 'mask', mesh, C%output_dir, climate%retreat%mask)
        call limit_retreat_mask( mesh, climate%retreat%mask)
        if (C%retreat_mask_applied_only_to_open_ocean) then
          do vi = mesh%vi1, mesh%vi2
            if (.not. climate%retreat%open_ocean( vi)) climate%retreat%mask( vi) = 0._dp
          end do
        end if
        climate%retreat%loaded = .true.
      end if

    else
      ! Time-dependent mask: interpolate linearly between the two enveloping frames of the
      ! file, holding the nearest frame outside the time axis

      n = size( climate%retreat%times)
      t = max( climate%retreat%times( 1), min( time, climate%retreat%times( n)))

      if (.not. climate%retreat%loaded .or. t < climate%retreat%t0 .or. t > climate%retreat%t1) then

        ! Find the enveloping frames
        i1 = 1
        do while (i1 < n)
          if (climate%retreat%times( i1) >= t) exit
          i1 = i1 + 1
        end do
        i0 = max( 1, i1-1)

        ! Frame 0: reuse the previous frame 1 when stepping forward, else read it
        if (climate%retreat%loaded .and. climate%retreat%times( i0) == climate%retreat%t1) then
          climate%retreat%mask0 = climate%retreat%mask1
        else
          call read_field_from_file_2D( filename, 'mask', mesh, C%output_dir, climate%retreat%mask0, &
            time_to_read = climate%retreat%times( i0))
          call limit_retreat_mask( mesh, climate%retreat%mask0)
        end if

        ! Frame 1
        if (i1 == i0) then
          climate%retreat%mask1 = climate%retreat%mask0
        else
          call read_field_from_file_2D( filename, 'mask', mesh, C%output_dir, climate%retreat%mask1, &
            time_to_read = climate%retreat%times( i1))
          call limit_retreat_mask( mesh, climate%retreat%mask1)
        end if

        climate%retreat%t0 = climate%retreat%times( i0)
        climate%retreat%t1 = climate%retreat%times( i1)
        climate%retreat%loaded = .true.

      end if

      w1 = 0._dp
      if (climate%retreat%t1 > climate%retreat%t0) &
        w1 = (t - climate%retreat%t0) / (climate%retreat%t1 - climate%retreat%t0)
      climate%retreat%mask = (1._dp - w1) * climate%retreat%mask0 + w1 * climate%retreat%mask1

    end if

    call checksum( mesh%pai_V, climate%retreat%mask, 'climate%retreat%mask')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_climate_retreat_mask

  subroutine remap_climate_retreat_mask( mesh_new, climate, time)
    !< Set up the retreat mask on a new mesh

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh_new
    type(type_climate_model), intent(inout) :: climate
    real(dp),                 intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_climate_retreat_mask'

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_retreat_mask( mesh_new, climate)

    ! Map the open-ocean reference from its original mesh, not from the previous model mesh
    if (C%retreat_mask_applied_only_to_open_ocean) call read_open_ocean_reference( mesh_new, climate)

    ! Reload immediately: the ice dynamics run before the climate in the next model step
    call run_climate_retreat_mask( mesh_new, climate, time)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_climate_retreat_mask

! ===== Utilities =====
! =====================

  subroutine allocate_retreat_mask( mesh, climate)
    !< (Re)allocate the retreat mask fields on the model mesh

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    type(type_climate_model), intent(inout) :: climate

    if (allocated( climate%retreat%mask0)) deallocate( climate%retreat%mask0)
    if (allocated( climate%retreat%mask1)) deallocate( climate%retreat%mask1)
    if (allocated( climate%retreat%mask )) deallocate( climate%retreat%mask )
    allocate( climate%retreat%mask0( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( climate%retreat%mask1( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( climate%retreat%mask ( mesh%vi1:mesh%vi2), source = 0._dp)

    if (C%retreat_mask_applied_only_to_open_ocean) then
      if (allocated( climate%retreat%open_ocean)) deallocate( climate%retreat%open_ocean)
      allocate( climate%retreat%open_ocean( mesh%vi1:mesh%vi2), source = .false.)
    end if

    climate%retreat%loaded = .false.

  end subroutine allocate_retreat_mask

  subroutine limit_retreat_mask( mesh, mask)
    !< Remove the small overshoots of the conservative remapping from a validated mask

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: mask

    if (any( .not. ieee_is_finite( mask))) call crash('non-finite retreat mask after remapping')
    mask = max( 0._dp, min( 1._dp, mask))

  end subroutine limit_retreat_mask

  subroutine check_retreat_mask_values( filename)
    !< Check that all values of the mask variable are finite, within [0,1], and not
    !< equal to a declared missing value
    !
    ! NetCDF missing values (_FillValue, missing_value, NaN) would otherwise enter the
    ! conservative remapping, and a large positive sentinel would activate retreat. A
    ! declared missing value within [0,1] cannot be told apart from data, so any value
    ! equal to it is rejected as well.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename

    ! Local variables:
    character(len=1024), parameter           :: routine_name = 'check_retreat_mask_values'
    real(dp), parameter                      :: tol = 1e-6_dp
    integer                                  :: ncid, id_var, id_dim_time, var_type, ndims_of_var, di, ti, nt, ierr
    integer, dimension( NF90_MAX_VAR_DIMS)   :: dims_of_var
    integer, dimension(3)                    :: n
    logical                                  :: has_time
    real(dp), dimension(:    ), allocatable  :: d1
    real(dp), dimension(:,:  ), allocatable  :: d2
    real(dp), dimension(:,:,:), allocatable  :: d3
    integer, dimension(3)                    :: n_bad          ! (non-finite, outside [0,1], equal to a missing value)
    real(dp), dimension(2)                   :: d_range
    real(dp), dimension(:), allocatable      :: missing_values
    real(dp)                                 :: missing_value_found

    ! Add routine to path
    call init_routine( routine_name)

    call open_existing_netcdf_file_for_reading( filename, ncid)

    call inquire_var_multopt( filename, ncid, 'mask', id_var, var_type = var_type, &
      ndims_of_var = ndims_of_var, dims_of_var = dims_of_var)
    if (id_var == -1) call crash('no variable "mask" in retreat mask file "' // trim( filename) // '"')

    ! Determine whether the last dimension is time
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time)
    has_time = id_dim_time /= -1 .and. ndims_of_var > 0
    if (has_time) has_time = dims_of_var( ndims_of_var) == id_dim_time
    if (C%retreat_mask_without_time .and. has_time) call crash('retreat mask in file "' // &
      trim( filename) // '" has a time dimension, but retreat_mask_without_time = .true.')
    if (.not. C%retreat_mask_without_time .and. .not. has_time) call crash('retreat mask in file "' // &
      trim( filename) // '" has no time dimension, but retreat_mask_without_time = .false.')

    n = 1
    if (ndims_of_var > 3) call crash('retreat mask in file "' // trim( filename) // '" has too many dimensions')
    do di = 1, ndims_of_var
      call inquire_dim_info( filename, ncid, dims_of_var( di), dim_length = n( di))
    end do

    ! Declared missing values (_FillValue and missing_value, which may be a vector)
    allocate( missing_values( 0))
    if (par%primary) then
      call add_missing_values( '_FillValue')
      call add_missing_values( 'missing_value')
    end if

    ! Read all values frame by frame on the primary and count invalid values
    n_bad   = 0
    d_range = [huge( 1._dp), -huge( 1._dp)]
    missing_value_found = 0._dp
    nt = 1
    if (has_time) nt = n( ndims_of_var)

    do ti = 1, nt
      select case (ndims_of_var)
      case default
        call crash('unsupported layout of retreat mask in file "' // trim( filename) // '"')
      case (1)
        allocate( d1( n(1)))
        call read_var_primary( filename, ncid, id_var, d1)
        if (par%primary) call count_invalid( reshape( d1, [size( d1)]))
        deallocate( d1)
      case (2)
        if (has_time) then
          allocate( d2( n(1), 1))
          call read_var_primary( filename, ncid, id_var, d2, start = [1, ti], count = [n(1), 1])
        else
          allocate( d2( n(1), n(2)))
          call read_var_primary( filename, ncid, id_var, d2)
        end if
        if (par%primary) call count_invalid( reshape( d2, [size( d2)]))
        deallocate( d2)
      case (3)
        if (.not. has_time) call crash('unsupported layout of retreat mask in file "' // trim( filename) // '"')
        allocate( d3( n(1), n(2), 1))
        call read_var_primary( filename, ncid, id_var, d3, start = [1, 1, ti], count = [n(1), n(2), 1])
        if (par%primary) call count_invalid( reshape( d3, [size( d3)]))
        deallocate( d3)
      end select
    end do

    call close_netcdf_file( ncid)

    call MPI_BCAST( n_bad              , 3, MPI_INTEGER         , 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( d_range            , 2, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( missing_value_found, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    if (sum( n_bad) > 0) call crash('retreat mask in file "' // trim( filename) // '" has {int_01} non-finite ' // &
      'values, {int_02} values outside [0,1] (range {dp_01} to {dp_02}), and {int_03} values equal to a ' // &
      'declared _FillValue or missing_value (e.g. {dp_03}). Missing values are not supported: set them ' // &
      'explicitly to 0 (no retreat) or 1 (retreat), and remove missing-value attributes that equal valid values.', &
      int_01 = n_bad( 1), int_02 = n_bad( 2), int_03 = n_bad( 3), dp_01 = d_range( 1), dp_02 = d_range( 2), &
      dp_03 = missing_value_found)

    ! Finalise routine path
    call finalise_routine( routine_name)

  contains

    subroutine add_missing_values( att_name)
      character(len=*), intent(in) :: att_name
      integer                             :: att_type, att_len
      real(dp), dimension(:), allocatable :: values
      if (NF90_INQUIRE_ATTRIBUTE( ncid, id_var, att_name, xtype = att_type, len = att_len) /= NF90_NOERR) return
      if (att_type == NF90_CHAR) call crash('attribute ' // att_name // ' of the retreat mask in file "' // &
        trim( filename) // '" is not numeric')
      allocate( values( att_len))
      if (NF90_GET_ATT( ncid, id_var, att_name, values) /= NF90_NOERR) call crash('could not read attribute ' // &
        att_name // ' of the retreat mask in file "' // trim( filename) // '"')
      missing_values = [missing_values, pack( values, ieee_is_finite( values))]
    end subroutine add_missing_values

    subroutine count_invalid( d)
      real(dp), dimension(:), intent(in) :: d
      integer :: i
      do i = 1, size( d)
        if (.not. ieee_is_finite( d( i))) then
          n_bad( 1) = n_bad( 1) + 1
        else
          d_range( 1) = min( d_range( 1), d( i))
          d_range( 2) = max( d_range( 2), d( i))
          if (d( i) < -tol .or. d( i) > 1._dp + tol) n_bad( 2) = n_bad( 2) + 1
          if (any( missing_values == d( i))) then
            n_bad( 3) = n_bad( 3) + 1
            missing_value_found = d( i)
          end if
        end if
      end do
    end subroutine count_invalid

  end subroutine check_retreat_mask_values

  subroutine check_retreat_mask_time_units( filename)
    !< Check that the time axis of the mask is given in years
    !
    ! Units relative to a reference date (e.g. "days since 1850-01-01") are not converted, so they are rejected.

    ! In/output variables:
    character(len=*), intent(in   ) :: filename

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'check_retreat_mask_time_units'
    integer                        :: ncid, id_var_time, att_type, att_len, units_status, ierr, i
    character(len=256)             :: units

    ! Add routine to path
    call init_routine( routine_name)

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)

    ! units_status: 0 = no units, 1 = years, 2 = relative to a reference date, 3 = other
    units = ''
    units_status = 0
    if (par%primary) then
      if (NF90_INQUIRE_ATTRIBUTE( ncid, id_var_time, 'units', xtype = att_type, len = att_len) == NF90_NOERR) then
        units_status = 3
        if (att_type == NF90_CHAR .and. att_len <= len( units)) then
          if (NF90_GET_ATT( ncid, id_var_time, 'units', units) == NF90_NOERR) then
            do i = 1, len_trim( units)
              if (units( i:i) >= 'A' .and. units( i:i) <= 'Z') units( i:i) = achar( iachar( units( i:i)) + 32)
              if (iachar( units( i:i)) == 0) units( i:i) = ' '
            end do
            units = adjustl( units)
            select case (trim( units))
            case ('years', 'year', 'yr', 'yrs', 'a')
              units_status = 1
            case default
              if (index( units, ' since ') > 0) units_status = 2
            end select
          end if
        end if
      end if
    end if

    call close_netcdf_file( ncid)

    call MPI_BCAST( units_status, 1        , MPI_INTEGER  , 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST( units       , len( units), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

    select case (units_status)
    case (0)
      call warning('time variable of retreat mask file "' // trim( filename) // &
        '" has no units attribute; its values are used as model time in years')
    case (1)
      ! Years, as used for the model time
    case (2)
      call crash('time units "' // trim( units) // '" of retreat mask file "' // trim( filename) // &
        '" are relative to a reference date; provide the time axis in years of model time')
    case default
      call crash('unsupported time units "' // trim( units) // '" of retreat mask file "' // &
        trim( filename) // '"; expected years')
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_retreat_mask_time_units

  subroutine report_retreat_mask_time_axis( times)
    !< Report the time axis and where frames are held or undersampled

    ! In/output variables:
    real(dp), dimension(:), intent(in   ) :: times

    ! Local variables:
    integer :: n

    n = size( times)

    if (par%primary) write(0,'(A,I0,A,F0.3,A,F0.3,A)') '     Retreat mask has ', n, &
      ' time frame(s) from t = ', times( 1), ' to t = ', times( n), ' yr'

    if (C%end_time_of_run < times( 1) .or. C%start_time_of_run > times( n)) then
      call warning('the run lies entirely outside the retreat mask time axis ({dp_01} to {dp_02} yr); ' // &
        'the nearest frame is held throughout', dp_01 = times( 1), dp_02 = times( n))
    elseif (C%start_time_of_run < times( 1) .or. C%end_time_of_run > times( n)) then
      call warning('the run extends beyond the retreat mask time axis ({dp_01} to {dp_02} yr); ' // &
        'the nearest frame is held there', dp_01 = times( 1), dp_02 = times( n))
    end if

    if (n > 1 .and. C%do_asynchronous_climate) then
      if (C%dt_climate > minval( times( 2:n) - times( 1:n-1))) &
        call warning('dt_climate ({dp_01} yr) exceeds the smallest spacing of the retreat mask frames; ' // &
          'the mask is only updated at climate time steps', dp_01 = C%dt_climate)
    end if

  end subroutine report_retreat_mask_time_axis

  subroutine write_open_ocean_reference( mesh, climate)
    !< Write the open-ocean reference on the current mesh, so that later mesh
    !< updates and restarted run segments can use exactly the same reference

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    type(type_climate_model), intent(in   ) :: climate

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'write_open_ocean_reference'
    character(len=1024)                    :: filename
    integer                                :: ncid
    real(dp), dimension(mesh%vi1:mesh%vi2) :: d

    ! Add routine to path
    call init_routine( routine_name)

    filename = climate%retreat%open_ocean_reference_filename

    d = 0._dp
    where (climate%retreat%open_ocean) d = 1._dp

    call create_new_netcdf_file_for_writing( filename, ncid)
    call add_attribute_char( filename, ncid, NF90_GLOBAL, 'region', climate%retreat%region_name)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call add_time_dimension_to_file( filename, ncid)
    call add_field_mesh_dp_2D( filename, ncid, 'open_ocean_reference', &
      long_name = 'Ice-free ocean in the initial geometry (retreat-mask eligibility)', units = '-')
    call write_time_to_file( filename, ncid, C%start_time_of_run)
    call write_to_field_multopt_mesh_dp_2D( mesh, filename, ncid, 'open_ocean_reference', d)
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_open_ocean_reference

  subroutine read_open_ocean_reference( mesh, climate)
    !< Map the open-ocean reference from the mesh it was written on to the model mesh
    !
    ! The reference must belong to the same model region and projection, and cover the
    ! model domain; its mesh and resolution may differ. Its values must be exactly 0 or 1.

    ! In/output variables:
    type(type_mesh),          intent(in   ) :: mesh
    type(type_climate_model), intent(inout) :: climate

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'read_open_ocean_reference'
    real(dp), parameter                    :: tol_angle = 1e-6_dp
    character(len=1024)                    :: filename
    character(len=16)                      :: region_ref
    integer                                :: ncid, id_dim_time, nt, vi, n_invalid, ierr
    type(type_mesh)                        :: mesh_ref
    real(dp), dimension(:), allocatable    :: d_ref
    real(dp), dimension(mesh%vi1:mesh%vi2) :: d
    real(dp)                               :: tol_dist

    ! Add routine to path
    call init_routine( routine_name)

    filename = climate%retreat%open_ocean_reference_filename

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_mesh_from_file( filename, ncid, mesh_ref)
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)
    region_ref = ''
    if (par%primary) then
      if (NF90_GET_ATT( ncid, NF90_GLOBAL, 'region', region_ref) /= NF90_NOERR) region_ref = ''
    end if
    call MPI_BCAST( region_ref, len( region_ref), MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    call close_netcdf_file( ncid)

    if (id_dim_time == -1 .or. nt /= 1) call crash('open-ocean reference file "' // trim( filename) // &
      '" must contain exactly one time frame')

    ! Region and projection must match; the mesh itself may differ
    if (trim( region_ref) /= climate%retreat%region_name) call crash('open-ocean reference file "' // &
      trim( filename) // '" belongs to region "' // trim( region_ref) // '", not to region "' // &
      climate%retreat%region_name // '"')
    if (abs( mesh_ref%lambda_M    - mesh%lambda_M   ) > tol_angle .or. &
        abs( mesh_ref%phi_M       - mesh%phi_M      ) > tol_angle .or. &
        abs( mesh_ref%beta_stereo - mesh%beta_stereo) > tol_angle) &
      call crash('open-ocean reference file "' // trim( filename) // '" uses a different map projection')
    tol_dist = 1e-9_dp * (mesh%xmax - mesh%xmin)
    if (mesh%xmin < mesh_ref%xmin - tol_dist .or. mesh%xmax > mesh_ref%xmax + tol_dist .or. &
        mesh%ymin < mesh_ref%ymin - tol_dist .or. mesh%ymax > mesh_ref%ymax + tol_dist) &
      call crash('open-ocean reference file "' // trim( filename) // '" does not cover the model domain')

    ! Give the reference mesh a unique name, so that no mapping object of an unrelated mesh
    ! with the same name (e.g. the first mesh of a restarted run) is reused
    mesh_ref%name = 'retreat_mask_open_ocean_reference_mesh'

    allocate( d_ref( mesh_ref%vi1:mesh_ref%vi2))
    call read_field_from_mesh_file_dp_2D( filename, 'open_ocean_reference', d_ref, time_to_read = C%start_time_of_run)

    ! The reference must be binary before it is mapped
    n_invalid = 0
    do vi = mesh_ref%vi1, mesh_ref%vi2
      if (.not. (d_ref( vi) == 0._dp .or. d_ref( vi) == 1._dp)) n_invalid = n_invalid + 1
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, n_invalid, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
    if (n_invalid > 0) call crash('open-ocean reference file "' // trim( filename) // &
      '" has {int_01} values that are not exactly 0 or 1', int_01 = n_invalid)

    ! Nearest-neighbour mapping keeps the reference binary
    call map_from_mesh_to_mesh_2D( mesh_ref, mesh, C%output_dir, d_ref, d, method = 'nearest_neighbour')
    climate%retreat%open_ocean( mesh%vi1:mesh%vi2) = d > 0.5_dp

    ! Release the mapping object of the reference mesh
    call clear_all_maps_involving_this_mesh( mesh_ref)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine read_open_ocean_reference

  function resolve_region_filename( filename_template, region_name) result( filename)
    !< Replace the token {region} in a file name by the region name

    ! In/output variables:
    character(len=*), intent(in) :: filename_template
    character(len=3), intent(in) :: region_name
    character(len=1024)          :: filename

    ! Local variables:
    integer :: i

    filename = filename_template
    i = index( filename, '{region}')
    do while (i > 0)
      filename = filename( 1:i-1) // region_name // filename( i+8:)
      i = index( filename, '{region}')
    end do

  end function resolve_region_filename

end module climate_retreat_mask
