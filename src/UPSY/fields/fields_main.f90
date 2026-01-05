module fields_main

  use fields_basic, only: &
    atype_field, atype_field_2D, atype_field_3D, &
    type_field_logical_2D, type_field_int_2D, type_field_dp_2D, &
    type_field_logical_3D, type_field_int_3D, type_field_dp_3D
  use fields_registry, only: type_fields_registry

  implicit none

  private

  public :: &
    atype_field, atype_field_2D, atype_field_3D, &
    type_field_logical_2D, type_field_int_2D, type_field_dp_2D, &
    type_field_logical_3D, type_field_int_3D, type_field_dp_3D
  public :: type_fields_registry

end module fields_main