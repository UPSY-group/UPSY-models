submodule (fields_basic) fields_basic_submod_print_info

  use mpi_f08, only: MPI_GATHER, MPI_INTEGER, MPI_COMM_WORLD

contains

  subroutine print_info( self)

    ! In/output variables:
    class(atype_field), intent(in) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'print_field_info'
    class(*), pointer              :: field_grid
    character(len=1024)            :: field_grid_type
    character(len=1024)            :: field_grid_name
    type(type_Arakawa_grid)        :: field_Arakawa_grid
    character(len=1024)            :: field_Arakawa_grid_name
    character(len=1024)            :: field_dimension
    character(len=1024)            :: field_data_type
    type(type_third_dimension)     :: field_third_dimension
    integer, dimension(0:par%n-1)  :: lbs, ubs
    integer, dimension(0:par%n-1)  :: lbs1, ubs1, lbs2, ubs2
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Grid
    field_grid              => self%grid()
    field_Arakawa_grid      =  self%Arakawa_grid()
    field_Arakawa_grid_name =  field_Arakawa_grid%str()

    select type (grid => field_grid)
    class default
      call crash('invalid field%grid class')
    class is (type_grid)
      field_grid_type = 'x/y-grid'
      field_grid_name = grid%name
    class is (type_mesh)
      field_grid_type = 'mesh'
      field_grid_name = grid%name
    end select

    ! Third dimension
    select type (f => self)
    class default
      call crash('invalid field class')
    class is (atype_field_2D)
      field_dimension = '2-D'
    class is (atype_field_3D)
      field_dimension = '3-D'
      field_third_dimension = f%third_dimension()
    end select

    ! Data type
    select type( f => self)
    class default
      call crash('invalid field class')
    class is (type_field_logical_2D)
      field_data_type = 'logical'
    class is (type_field_logical_3D)
      field_data_type = 'logical'
    class is (type_field_int_2D)
      field_data_type = 'int'
    class is (type_field_int_3D)
      field_data_type = 'int'
    class is (type_field_dp_2D)
      field_data_type = 'dp'
    class is (type_field_dp_3D)
      field_data_type = 'dp'
    end select

    ! Bounds
    select case (field_dimension)
    case default
      call crash('invalid field_dimension')
    case ('2-D')
      call gather_field_bounds_2D( self, lbs, ubs)
    case ('3-D')
      call gather_field_bounds_3D( self, lbs1, ubs1, lbs2, ubs2)
    end select

    if (par%primary) then
      write(0,*) ''
      write(0,*) '     Field: ', colour_string( trim( self%name()),'light blue')

      write(0,*) '       Long name: ', trim( self%long_name())
      write(0,*) '       Parent   : ', trim( field_grid_type), ' "', &
        trim( field_grid_name), '" (', trim( field_Arakawa_grid_name), '-grid)'

      select case (field_dimension)
      case default
        call crash('invalid field_dimension')
      case ('2-D')
        write(0,*) '       Type     : 2-D, ', trim( field_data_type)
      case ('3-D')
        write(0,*) '       Type     : 3-D (', trim( field_third_dimension%name), '), ', &
          trim( field_data_type)
      end select

      write(0,*) '       Units    : [', trim( self%units()), ']'

      select case (field_dimension)
      case default
        call crash('invalid field_dimension')
      case ('2-D')
        do i = 0, par%n-1
          if (i == 0) then
            write(0,*) '       Bounds   : process ', i, ' owns [', lbs(i), ' - ', ubs(i), ']'
          else
            write(0,*) '                          ', i, '      [', lbs(i), ' - ', ubs(i), ']'
          end if
        end do
      case ('3-D')
        do i = 0, par%n-1
          if (i == 0) then
            write(0,*) '       Bounds   : process ', i, ' owns [', lbs1(i), ' - ', ubs1(i), ']', &
              ', [', lbs2(i), ' - ', ubs2(i), ']'
          else
            write(0,*) '                          ', i, '      [', lbs1(i), ' - ', ubs1(i), ']', &
              ', [', lbs2(i), ' - ', ubs2(i), ']'
          end if
        end do
      end select

    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine print_info

  function field_lbound( self, dim) result( lb)

    ! In/output variables:
    class(atype_field), intent(in) :: self
    integer,            intent(in) :: dim
    integer                        :: lb

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'field_lbound'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f => self)
    class default
      call crash('invalid field type')
    class is (type_field_logical_2D)
      lb = lbound( f%d_nih, dim)
    class is (type_field_int_2D)
      lb = lbound( f%d_nih, dim)
    class is (type_field_dp_2D)
      lb = lbound( f%d_nih, dim)
    class is (type_field_logical_3D)
      lb = lbound( f%d_nih, dim)
    class is (type_field_int_3D)
      lb = lbound( f%d_nih, dim)
    class is (type_field_dp_3D)
      lb = lbound( f%d_nih, dim)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function field_lbound

  function field_ubound( self, dim) result( ub)

    ! In/output variables:
    class(atype_field), intent(in) :: self
    integer,            intent(in) :: dim
    integer                        :: ub

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'field_ubound'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f => self)
    class default
      call crash('invalid field type')
    class is (type_field_logical_2D)
      ub = ubound( f%d_nih, dim)
    class is (type_field_int_2D)
      ub = ubound( f%d_nih, dim)
    class is (type_field_dp_2D)
      ub = ubound( f%d_nih, dim)
    class is (type_field_logical_3D)
      ub = ubound( f%d_nih, dim)
    class is (type_field_int_3D)
      ub = ubound( f%d_nih, dim)
    class is (type_field_dp_3D)
      ub = ubound( f%d_nih, dim)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function field_ubound

  subroutine gather_field_bounds_2D( field, lbs, ubs)

    ! In/output variables:
    class(atype_field),            intent(in   ) :: field
    integer, dimension(0:par%n-1), intent(  out) :: lbs, ubs

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'gather_field_bounds_2D'
    integer                        :: lb, ub, ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    lb = field%lbound( 1)
    ub = field%ubound( 1)

    call MPI_GATHER( lb, 1, MPI_INTEGER, lbs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER( ub, 1, MPI_INTEGER, ubs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine gather_field_bounds_2D

  subroutine gather_field_bounds_3D( field, lbs1, ubs1, lbs2, ubs2)

    ! In/output variables:
    class(atype_field),            intent(in   ) :: field
    integer, dimension(0:par%n-1), intent(  out) :: lbs1, ubs1, lbs2, ubs2

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'gather_field_bounds_3D'
    integer                        :: lb1, ub1, lb2, ub2, ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    lb1 = field%lbound( 1)
    ub1 = field%ubound( 1)
    lb2 = field%lbound( 2)
    ub2 = field%ubound( 2)

    call MPI_GATHER( lb1, 1, MPI_INTEGER, lbs1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER( ub1, 1, MPI_INTEGER, ubs1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER( lb2, 1, MPI_INTEGER, lbs2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER( ub2, 1, MPI_INTEGER, ubs2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine gather_field_bounds_3D

end submodule fields_basic_submod_print_info