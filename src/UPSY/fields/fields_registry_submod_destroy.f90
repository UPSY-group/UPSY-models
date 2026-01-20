submodule( fields_registry) fields_registry_submod_destroy

contains

  subroutine destroy( flds_reg)

    ! In/output variables:
    class(type_fields_registry), intent(inout) :: flds_reg

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'type_fields_registry_destroy'
    integer                         :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    do i = 1, flds_reg%n

      ! Deallocate the shared memory of the actual array
      select type (f => flds_reg%items(i)%p)
      class default
        call crash('invalid field class')
      class is (type_field_logical_2D)
        call deallocate_dist_shared( f%d_nih, f%w)
      class is (type_field_logical_3D)
        call deallocate_dist_shared( f%d_nih, f%w)
      class is (type_field_int_2D)
        call deallocate_dist_shared( f%d_nih, f%w)
      class is (type_field_int_3D)
        call deallocate_dist_shared( f%d_nih, f%w)
      class is (type_field_dp_2D)
        call deallocate_dist_shared( f%d_nih, f%w)
      class is (type_field_dp_3D)
        call deallocate_dist_shared( f%d_nih, f%w)
      end select

      ! Deallocate the field
      deallocate( flds_reg%items(i)%p)

    end do

    deallocate( flds_reg%items)
    flds_reg%n     = 0
    flds_reg%n_max = 0

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine destroy

end submodule fields_registry_submod_destroy