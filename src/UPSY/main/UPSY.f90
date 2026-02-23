module UPSY_main
  !< The main UPSY module that should be imported into models built with it.

use models_basic, only: atype_model, atype_model_context_allocate, &
  atype_model_context_initialise, atype_model_context_run, atype_model_context_remap
use string_module, only: type_string_utilities

implicit none

private

public :: UPSY
public :: atype_model, atype_model_context_allocate, &
  atype_model_context_initialise, atype_model_context_run, atype_model_context_remap

type type_UPSY
    !< The complete UPSY toolkit
    private
    type(type_string_utilities), public :: stru
  contains
end type type_UPSY

type(type_UPSY) :: UPSY

end module UPSY_main