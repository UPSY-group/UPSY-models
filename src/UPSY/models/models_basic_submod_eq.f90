submodule( models_basic) models_basic_submod_eq

contains

  function test_model_equality( model1, model2) result( res)

    ! In/output variables:
    class(atype_model), intent(in) :: model1, model2
    logical                        :: res

    res = model1%is_name( model2%name()) .and. model1%mesh%name == model2%mesh%name
    if (.not. res) return
    res = res .and. model1%flds_reg == model2%flds_reg

  end function test_model_equality

end submodule models_basic_submod_eq
