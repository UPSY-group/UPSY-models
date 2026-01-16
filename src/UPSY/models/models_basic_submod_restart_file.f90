submodule( models_basic) models_basic_submod_restart_file

use netcdf_io_main

contains

  subroutine write_to_restart_file( model, output_dir, filename)

    ! In/output variables:
    class(atype_model),                  intent(in   ) :: model
    character(len=*),                    intent(in   ) :: output_dir
    character(:), allocatable, optional, intent(  out) :: filename

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_restart_file'
    character(:), allocatable      :: filename_loc
    integer                        :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    filename_loc = trim(  output_dir) // '/restart_file_' // trim( model%name()) // '.nc'
    if (present( filename)) filename = filename_loc

    call create_new_netcdf_file_for_writing( filename_loc, ncid)

    ! Set up grid/mesh
    select type (grid => model%grid_val)
    class default
      call crash('invalid model%grid type - programming error')
    class is (type_grid)
      call setup_xy_grid_in_netcdf_file( filename_loc, ncid, grid)
    class is (type_mesh)
      call setup_mesh_in_netcdf_file( filename_loc, ncid, grid)
    end select

    call model%flds_reg%write_to_netcdf( filename_loc, ncid)

    call close_netcdf_file( ncid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file

  subroutine read_from_restart_file( model, filename)

    ! In/output variables:
    class(atype_model), intent(inout) :: model
    character(len=*),   intent(in   ) :: filename

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'read_from_restart_file'
    integer                        :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    call open_existing_netcdf_file_for_reading( filename, ncid)
    call model%flds_reg%read_from_netcdf( filename, ncid)
    call close_netcdf_file( ncid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_restart_file

end submodule models_basic_submod_restart_file
