submodule (fields_basic) fields_basic_submod_print_info

  use mpi_f08, only: MPI_GATHER, MPI_INTEGER, MPI_COMM_WORLD

contains

  subroutine print_info( field)

    ! In/output variables:
    class(atype_field), intent(in) :: field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'print_field_info'
    character(len=1024)            :: field_name
    character(len=1024)            :: field_long_name
    character(len=1024)            :: field_dimension
    character(len=1024)            :: field_parent_grid_type
    character(len=1024)            :: field_parent_grid_name
    type(type_Arakawa_grid)        :: field_parent_Arakawa_grid
    character(len=1024)            :: field_parent_Arakawa_grid_str
    type(type_third_dimension)     :: field_third_dimension
    character(len=1024)            :: field_third_dimension_name
    character(len=1024)            :: field_units
    integer, dimension(0:par%n-1)  :: lbs, ubs
    integer, dimension(0:par%n-1)  :: lbs1, ubs1, lbs2, ubs2
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    field_name                    = field%name()
    field_long_name               = field%long_name()
    field_units                   = field%units()
    field_parent_Arakawa_grid     = field%Arakawa_grid()
    field_parent_Arakawa_grid_str = field_parent_Arakawa_grid%str()

    ! Field dimensions
    field_dimension            = ''
    field_third_dimension_name = ''

    select type (p => field)
    class default
      call crash('invalid field type')
    class is (type_field_grid_logical_2D)
      field_dimension = '2-D'
      call gather_field_bounds_2D( field, lbs, ubs)
    class is (type_field_grid_logical_3D)
      field_dimension = '3-D'
      call gather_field_bounds_3D( field, lbs1, ubs1, lbs2, ubs2)
      field_third_dimension = p%third_dimension()
      field_third_dimension_name = field_third_dimension%name
    class is (type_field_grid_int_2D)
      field_dimension = '2-D'
      call gather_field_bounds_2D( field, lbs, ubs)
    class is (type_field_grid_int_3D)
      field_dimension = '3-D'
      call gather_field_bounds_3D( field, lbs1, ubs1, lbs2, ubs2)
      field_third_dimension = p%third_dimension()
      field_third_dimension_name = field_third_dimension%name
    class is (type_field_grid_dp_2D)
      field_dimension = '2-D'
      call gather_field_bounds_2D( field, lbs, ubs)
    class is (type_field_grid_dp_3D)
      field_dimension = '3-D'
      call gather_field_bounds_3D( field, lbs1, ubs1, lbs2, ubs2)
      field_third_dimension = p%third_dimension()
      field_third_dimension_name = field_third_dimension%name
    class is (type_field_mesh_logical_2D)
      field_dimension = '2-D'
      call gather_field_bounds_2D( field, lbs, ubs)
    class is (type_field_mesh_logical_3D)
      field_dimension = '3-D'
      call gather_field_bounds_3D( field, lbs1, ubs1, lbs2, ubs2)
      field_third_dimension = p%third_dimension()
      field_third_dimension_name = field_third_dimension%name
    class is (type_field_mesh_int_2D)
      field_dimension = '2-D'
      call gather_field_bounds_2D( field, lbs, ubs)
    class is (type_field_mesh_int_3D)
      field_dimension = '3-D'
      call gather_field_bounds_3D( field, lbs1, ubs1, lbs2, ubs2)
      field_third_dimension = p%third_dimension()
      field_third_dimension_name = field_third_dimension%name
    class is (type_field_mesh_dp_2D)
      field_dimension = '2-D'
      call gather_field_bounds_2D( field, lbs, ubs)
    class is (type_field_mesh_dp_3D)
      field_dimension = '3-D'
      call gather_field_bounds_3D( field, lbs1, ubs1, lbs2, ubs2)
      field_third_dimension = p%third_dimension()
      field_third_dimension_name = field_third_dimension%name
    end select

    ! Parent grid can be either an x/y-grid or a mesh
    select type (p => field)
    class default
      call crash('invalid field type')
    class is (atype_field_grid)
      field_parent_grid_type = 'grid'
      field_parent_grid_name = p%parent%name
      class is (atype_field_mesh)
      field_parent_grid_type = 'mesh'
      field_parent_grid_name = p%parent%name
    end select

    if (par%primary) then
      write(0,*) '     Field: ', colour_string( trim( field_name),'light blue')

      write(0,*) '       Long name: ', trim( field_long_name)
      write(0,*) '       Parent   : ', trim( field_parent_grid_type), ' "', &
      trim( field_parent_grid_name), '" (', trim( field_parent_Arakawa_grid_str), '-grid)'
      write(0,*) '       Units    : [', trim( field_units), ']'

      select case (field_dimension)
      case default
        call crash('invalid field_dimension')
      case ('2-D')
        do i = 0, par%n-1
          if (i == 0) then
            write(0,*) '       Dimension: 2-D - Process ', i, ' owns [', lbs(i), ' - ', ubs(i), ']'
          else
            write(0,*) '                                ', i, '      [', lbs(i), ' - ', ubs(i), ']'
          end if
        end do
      case ('3-D')
        do i = 0, par%n-1
          if (i == 0) then
            write(0,*) '       Dimension: 3-D - Process ', i, ' owns [', lbs1(i), ' - ', ubs1(i), ']', &
              ', [', lbs2(i), ' - ', ubs2(i), ']'
          else
            write(0,*) '                                ', i, '      [', lbs1(i), ' - ', ubs1(i), ']', &
              ', [', lbs2(i), ' - ', ubs2(i), ']'
          end if
        end do
      end select

    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine print_info

  subroutine gather_field_bounds_2D( field, lbs, ubs)

    ! In/output variables:
    class(atype_field),            intent(in   ) :: field
    integer, dimension(0:par%n-1), intent(  out) :: lbs, ubs

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'gather_field_bounds'
    integer                        :: lb, ub, ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (p => field)
    class default
      call crash('invalid field type')
    class is (type_field_grid_logical_2D)
      lb = lbound( p%d,1)
      ub = ubound( p%d,1)
    class is (type_field_grid_int_2D)
      lb = lbound( p%d,1)
      ub = ubound( p%d,1)
    class is (type_field_grid_dp_2D)
      lb = lbound( p%d,1)
      ub = ubound( p%d,1)
    class is (type_field_mesh_logical_2D)
      lb = lbound( p%d,1)
      ub = ubound( p%d,1)
    class is (type_field_mesh_int_2D)
      lb = lbound( p%d,1)
      ub = ubound( p%d,1)
    class is (type_field_mesh_dp_2D)
      lb = lbound( p%d,1)
      ub = ubound( p%d,1)
    end select

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

    select type (p => field)
    class default
      call crash('invalid field type')
    class is (type_field_grid_logical_3D)
      lb1 = lbound( p%d,1)
      ub1 = ubound( p%d,1)
      lb2 = lbound( p%d,2)
      ub2 = ubound( p%d,2)
    class is (type_field_grid_int_3D)
      lb1 = lbound( p%d,1)
      ub1 = ubound( p%d,1)
      lb2 = lbound( p%d,2)
      ub2 = ubound( p%d,2)
    class is (type_field_grid_dp_3D)
      lb1 = lbound( p%d,1)
      ub1 = ubound( p%d,1)
      lb2 = lbound( p%d,2)
      ub2 = ubound( p%d,2)
    class is (type_field_mesh_logical_3D)
      lb1 = lbound( p%d,1)
      ub1 = ubound( p%d,1)
      lb2 = lbound( p%d,2)
      ub2 = ubound( p%d,2)
    class is (type_field_mesh_int_3D)
      lb1 = lbound( p%d,1)
      ub1 = ubound( p%d,1)
      lb2 = lbound( p%d,2)
      ub2 = ubound( p%d,2)
    class is (type_field_mesh_dp_3D)
      lb1 = lbound( p%d,1)
      ub1 = ubound( p%d,1)
      lb2 = lbound( p%d,2)
      ub2 = ubound( p%d,2)
    end select

    call MPI_GATHER( lb1, 1, MPI_INTEGER, lbs1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER( ub1, 1, MPI_INTEGER, ubs1, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER( lb2, 1, MPI_INTEGER, lbs2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_GATHER( ub2, 1, MPI_INTEGER, ubs2, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine gather_field_bounds_3D

end submodule fields_basic_submod_print_info