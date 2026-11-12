submodule( fields_registry) fields_registry_submod_checksum

contains

  subroutine checksum_fields_registry( self)
    class(type_fields_registry), intent(in) :: self
    integer :: i
    do i = 1, self%n
      call self%items(i)%p%checksum()
    end do
  end subroutine checksum_fields_registry

end submodule fields_registry_submod_checksum