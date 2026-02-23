submodule( fields_registry) fields_registry_submod_io

  use fields_dimensions, only: type_third_dimension
  use fields_basic, only: atype_field_2D, atype_field_3D

contains

  subroutine write_to_netcdf( self, filename, ncid)
    ! NOTE: assumes the NetCDF file already exists, is open, and
    !       already contains the grid/mesh

    ! In/output variables:
    class(type_fields_registry), intent(in) :: self
    character(len=*),            intent(in) :: filename
    integer,                     intent(in) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_netcdf'
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    do i = 1, self%n
      call self%items(i)%p%write_to_netcdf( filename, ncid)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf

  subroutine read_from_netcdf( self, filename, ncid)
    ! NOTE: assumes the NetCDF file is defined on
    !       the same grid/mesh as the fields

    ! In/output variables:
    class(type_fields_registry), intent(inout) :: self
    character(len=*),            intent(in   ) :: filename
    integer,                     intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter                        :: routine_name = 'read_from_netcdf'
    integer                                               :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    do i = 1, self%n
      call self%items(i)%p%read_from_netcdf( filename, ncid)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf

end submodule fields_registry_submod_io