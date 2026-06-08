module ISMIP7_forcing_field_types

  use precisions, only: dp

  implicit none

  private

  public :: type_climate_field_ISMIP7_monthly, type_climate_field_ISMIP7_yearly

  type type_climate_field_ISMIP7_monthly
    ! Data and metadata of monthly tas/pr fields and anomalies

    character(len=1024)                            :: name          !           'tas', 'pr', etc
    character(len=1024)                            :: foldername    !           Foldername that contains all files
    character(len=1024), dimension(:), allocatable :: filenames     !           Filenames

    real(dp), dimension(:), allocatable            :: timestamps    ! [years]   All time values in combined files

    real(dp), dimension(:,:), allocatable          :: val0          !           Values at timeslice before current time
    real(dp), dimension(:,:), allocatable          :: val1          !           Values at timeslice after current time
    real(dp), dimension(:,:), allocatable          :: val_interp    !           Interpolated values

    integer                                        :: ti0 = -1      !           Time index before current time
    integer                                        :: ti1 = -1      !           Time index after current time

  end type type_climate_field_ISMIP7_monthly

  type type_climate_field_ISMIP7_yearly
    ! Data and metadata of yearly fields (dtsdz)

    character(len=1024)                            :: name          !           'dtsdz', etc
    character(len=1024)                            :: foldername    !           Foldername that contains all files
    character(len=1024), dimension(:), allocatable :: filenames     !           Filenames

    real(dp), dimension(:), allocatable            :: timestamps    ! [years]   All time values in combined files

    real(dp), dimension(:), allocatable            :: val0          !           Values at timeslice before current time
    real(dp), dimension(:), allocatable            :: val1          !           Values at timeslice after current time
    real(dp), dimension(:), allocatable            :: val_interp    !           Interpolated values

    integer                                        :: ti0 = -1      !           Time index before current time
    integer                                        :: ti1 = -1      !           Time index after current time

  end type type_climate_field_ISMIP7_yearly

  end module ISMIP7_forcing_field_types