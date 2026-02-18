submodule( fields_registry) fields_registry_submod_reallocate_field

contains

  subroutine reallocate_field_logical_2D( self, mesh_new, field_name, d_nih)
    class(type_fields_registry),                intent(inout) :: self
    character(len=*),                           intent(in   ) :: field_name
    type(type_mesh),                            intent(in   ) :: mesh_new
    logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    call self%items( self%find( field_name))%p%reallocate( mesh_new, d_nih)
  end subroutine reallocate_field_logical_2D

  subroutine reallocate_field_logical_3D( self, mesh_new, field_name, d_nih)
    class(type_fields_registry),                  intent(inout) :: self
    character(len=*),                             intent(in   ) :: field_name
    type(type_mesh),                              intent(in   ) :: mesh_new
    logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    call self%items( self%find( field_name))%p%reallocate( mesh_new, d_nih)
  end subroutine reallocate_field_logical_3D

  subroutine reallocate_field_int_2D( self, mesh_new, field_name, d_nih)
    class(type_fields_registry),                intent(inout) :: self
    character(len=*),                           intent(in   ) :: field_name
    type(type_mesh),                            intent(in   ) :: mesh_new
    integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    call self%items( self%find( field_name))%p%reallocate( mesh_new, d_nih)
  end subroutine reallocate_field_int_2D

  subroutine reallocate_field_int_3D( self, mesh_new, field_name, d_nih)
    class(type_fields_registry),                  intent(inout) :: self
    character(len=*),                             intent(in   ) :: field_name
    type(type_mesh),                              intent(in   ) :: mesh_new
    integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    call self%items( self%find( field_name))%p%reallocate( mesh_new, d_nih)
  end subroutine reallocate_field_int_3D

  subroutine reallocate_field_dp_2D( self, mesh_new, field_name, d_nih)
    class(type_fields_registry),                 intent(inout) :: self
    character(len=*),                            intent(in   ) :: field_name
    type(type_mesh),                             intent(in   ) :: mesh_new
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
    call self%items( self%find( field_name))%p%reallocate( mesh_new, d_nih)
  end subroutine reallocate_field_dp_2D

  subroutine reallocate_field_dp_3D( self, mesh_new, field_name, d_nih)
    class(type_fields_registry),                   intent(inout) :: self
    character(len=*),                              intent(in   ) :: field_name
    type(type_mesh),                               intent(in   ) :: mesh_new
    real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    call self%items( self%find( field_name))%p%reallocate( mesh_new, d_nih)
  end subroutine reallocate_field_dp_3D

end submodule fields_registry_submod_reallocate_field