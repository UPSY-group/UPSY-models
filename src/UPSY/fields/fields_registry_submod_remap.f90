submodule( fields_registry) fields_registry_submod_remap

contains

  subroutine remap( flds_reg, mesh_new)

    ! In/output variables:
    class(type_fields_registry), intent(inout) :: flds_reg
    type(type_mesh),             intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'type_fields_registry_remap'
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    do i = 1, flds_reg%n
      call flds_reg%items(i)%p%remap( mesh_new)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap

end submodule fields_registry_submod_remap