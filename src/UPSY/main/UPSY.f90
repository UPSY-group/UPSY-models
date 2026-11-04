module UPSY_main
  !< The main UPSY module that should be imported into models built with it.

use models_basic, only: atype_model
use string_module, only: type_string_utilities
use remapping_basic, only: type_remapping

implicit none

private

public :: UPSY
public :: atype_model

type type_UPSY
    !< The complete UPSY toolkit
    private
    type(type_string_utilities), public :: stru
    type(type_remapping),        public :: remapping
  contains
end type type_UPSY

type(type_UPSY) :: UPSY

end module UPSY_main