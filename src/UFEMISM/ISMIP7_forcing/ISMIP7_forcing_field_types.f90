module ISMIP7_forcing_field_types

  use precisions, only: dp

  implicit none

  private

  public :: type_ISMIP7_forcing_field_monthly, type_ISMIP7_forcing_field_yearly

  type, abstract :: type_ISMIP7_forcing_field
    !< Metadata of monthly/yearly ISMIP7 forcing fields

    character(len=1024)                            :: name          !           'tas', 'pr', etc
    character(len=1024)                            :: foldername    !           Foldername that contains all files
    character(len=1024), dimension(:), allocatable :: filenames     !           Filenames

    real(dp), dimension(:), allocatable            :: timestamps    ! [years]   All time values in combined files

    integer                                        :: ti0 = -1      !           Time index before current time
    integer                                        :: ti1 = -1      !           Time index after current time

  end type type_ISMIP7_forcing_field

  type, extends(type_ISMIP7_forcing_field) :: type_ISMIP7_forcing_field_monthly
    !< Two enveloping timeframes and time-interpolated values of a single monthly ISMIP7 forcing field

    real(dp), dimension(:,:), allocatable          :: val0          !           Values at timeslice before current time
    real(dp), dimension(:,:), allocatable          :: val1          !           Values at timeslice after current time
    real(dp), dimension(:,:), allocatable          :: val_interp    !           Interpolated values

  end type type_ISMIP7_forcing_field_monthly

  type, extends(type_ISMIP7_forcing_field) :: type_ISMIP7_forcing_field_yearly
    !< Two enveloping timeframes and time-interpolated values of a single yearly ISMIP7 forcing field

    real(dp), dimension(:), allocatable            :: val0          !           Values at timeslice before current time
    real(dp), dimension(:), allocatable            :: val1          !           Values at timeslice after current time
    real(dp), dimension(:), allocatable            :: val_interp    !           Interpolated values

  end type type_ISMIP7_forcing_field_yearly

  end module ISMIP7_forcing_field_types