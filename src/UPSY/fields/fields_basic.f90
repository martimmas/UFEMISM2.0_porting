module fields_basic

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid, type_Arakawa_grid
  use fields_dimensions, only: third_dimension, type_third_dimension
  use tests_main, only: test_tol_grid, test_tol_mesh
  use mpi_f08, only: MPI_WIN, MPI_GATHER, MPI_INTEGER, MPI_COMM_WORLD

  implicit none

  private

  public :: &
    atype_field, atype_field_grid, atype_field_mesh, &
    type_field_grid_logical_2D, type_field_grid_logical_3D, &
    type_field_grid_int_2D, type_field_grid_int_3D, &
    type_field_grid_dp_2D, type_field_grid_dp_3D, &
    type_field_mesh_logical_2D, type_field_mesh_logical_3D, &
    type_field_mesh_int_2D, type_field_mesh_int_3D, &
    type_field_mesh_dp_2D, type_field_mesh_dp_3D

  ! Abstract basic field type and mesh/grid-based field types
  ! =========================================================

  type, abstract :: atype_field

    ! Metadata
    character(len=1024) :: name
    character(len=1024) :: long_name
    character(len=1024) :: units

    ! Parent grid
    ! type(type_grid), pointer :: parent    ! Defined in extended types atype_field_grid and atype_field_mesh!
    ! type(type_mesh), pointer :: parent
    type(type_Arakawa_grid)  :: Arakawa_grid

    ! Pointer to array containing the actual field data
    ! logical, pointer :: d(:)              ! Defined in extended types type_field_grid_XXX, etc.
    type(MPI_WIN) :: w

  contains

    procedure :: print_field_info
    procedure :: is_parent
    procedure :: is_Arakawa_grid
    procedure :: is_third_dimension

  end type atype_field

  type, abstract, extends(atype_field) :: atype_field_grid
    type(type_grid), pointer :: parent
  end type atype_field_grid

  type, abstract, extends(atype_field) :: atype_field_mesh
    type(type_mesh), pointer :: parent
  end type atype_field_mesh

  ! Concrete field types
  ! ====================

  ! Grid-based fields

  type, extends(atype_field_grid) :: type_field_grid_logical_2D
    logical, pointer :: d(:)
  end type type_field_grid_logical_2D

  type, extends(atype_field_grid) :: type_field_grid_logical_3D
    type(type_third_dimension) :: third_dimension
    logical, pointer :: d(:,:)
  end type type_field_grid_logical_3D

  type, extends(atype_field_grid) :: type_field_grid_int_2D
    integer, pointer :: d(:)
  end type type_field_grid_int_2D

  type, extends(atype_field_grid) :: type_field_grid_int_3D
    type(type_third_dimension) :: third_dimension
    integer, pointer :: d(:,:)
  end type type_field_grid_int_3D

  type, extends(atype_field_grid) :: type_field_grid_dp_2D
    real(dp), pointer :: d(:)
  end type type_field_grid_dp_2D

  type, extends(atype_field_grid) :: type_field_grid_dp_3D
    type(type_third_dimension) :: third_dimension
    real(dp), pointer :: d(:,:)
  end type type_field_grid_dp_3D

  ! Mesh-based fields

  type, extends(atype_field_mesh) :: type_field_mesh_logical_2D
    logical, pointer :: d(:)
  end type type_field_mesh_logical_2D

  type, extends(atype_field_mesh) :: type_field_mesh_logical_3D
    type(type_third_dimension) :: third_dimension
    logical, pointer :: d(:,:)
  end type type_field_mesh_logical_3D

  type, extends(atype_field_mesh) :: type_field_mesh_int_2D
    integer, pointer :: d(:)
  end type type_field_mesh_int_2D

  type, extends(atype_field_mesh) :: type_field_mesh_int_3D
    type(type_third_dimension) :: third_dimension
    integer, pointer :: d(:,:)
  end type type_field_mesh_int_3D

  type, extends(atype_field_mesh) :: type_field_mesh_dp_2D
    real(dp), pointer :: d(:)
  end type type_field_mesh_dp_2D

  type, extends(atype_field_mesh) :: type_field_mesh_dp_3D
    type(type_third_dimension) :: third_dimension
    real(dp), pointer :: d(:,:)
  end type type_field_mesh_dp_3D

contains

  subroutine print_field_info( field)

    ! In/output variables:
    class(atype_field), intent(in) :: field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'print_field_info'
    character(len=1024)            :: field_name
    character(len=1024)            :: field_long_name
    character(len=1024)            :: field_dimension
    character(len=1024)            :: field_parent_grid_type
    character(len=1024)            :: field_parent_grid_name
    character(len=1024)            :: field_parent_Arakawa_grid
    character(len=1024)            :: field_third_dimension_name
    character(len=1024)            :: field_units
    integer, dimension(0:par%n-1)  :: lbs, ubs
    integer, dimension(0:par%n-1)  :: lbs1, ubs1, lbs2, ubs2
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    field_name                = field%name
    field_long_name           = field%long_name
    field_units               = field%units
    field_parent_Arakawa_grid = field%Arakawa_grid%str()

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
      field_third_dimension_name = p%third_dimension%name
    class is (type_field_grid_int_2D)
      field_dimension = '2-D'
      call gather_field_bounds_2D( field, lbs, ubs)
    class is (type_field_grid_int_3D)
      field_dimension = '3-D'
      call gather_field_bounds_3D( field, lbs1, ubs1, lbs2, ubs2)
      field_third_dimension_name = p%third_dimension%name
    class is (type_field_grid_dp_2D)
      field_dimension = '2-D'
      call gather_field_bounds_2D( field, lbs, ubs)
    class is (type_field_grid_dp_3D)
      field_dimension = '3-D'
      call gather_field_bounds_3D( field, lbs1, ubs1, lbs2, ubs2)
      field_third_dimension_name = p%third_dimension%name
    class is (type_field_mesh_logical_2D)
      field_dimension = '2-D'
      call gather_field_bounds_2D( field, lbs, ubs)
    class is (type_field_mesh_logical_3D)
      field_dimension = '3-D'
      call gather_field_bounds_3D( field, lbs1, ubs1, lbs2, ubs2)
      field_third_dimension_name = p%third_dimension%name
    class is (type_field_mesh_int_2D)
      field_dimension = '2-D'
      call gather_field_bounds_2D( field, lbs, ubs)
    class is (type_field_mesh_int_3D)
      field_dimension = '3-D'
      call gather_field_bounds_3D( field, lbs1, ubs1, lbs2, ubs2)
      field_third_dimension_name = p%third_dimension%name
    class is (type_field_mesh_dp_2D)
      field_dimension = '2-D'
      call gather_field_bounds_2D( field, lbs, ubs)
    class is (type_field_mesh_dp_3D)
      field_dimension = '3-D'
      call gather_field_bounds_3D( field, lbs1, ubs1, lbs2, ubs2)
      field_third_dimension_name = p%third_dimension%name
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
      trim( field_parent_grid_name), '" (', trim( field_parent_Arakawa_grid), '-grid)'
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

  end subroutine print_field_info

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

  function is_parent( field, parent) result( res)

    ! In/output variables:
    class(atype_field), intent(in) :: field
    class(*),           intent(in) :: parent
    logical                        :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_parent'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => parent)
    class is (type_grid)
      res = is_parent_grid( field, p)
    class is (type_mesh)
      res = is_parent_mesh( field, p)
    class default
      res = .false.
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_parent

  function is_parent_grid( field, grid) result( res)

    ! In/output variables:
    class(atype_field), intent(in) :: field
    type(type_grid),    intent(in) :: grid
    logical                        :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_parent_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => field)
    class is (atype_field_grid)
      res = test_tol_grid( p%parent, grid, p%parent%tol_dist)
    class default
      res = .false.
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_parent_grid

  function is_parent_mesh( field, mesh) result( res)

    ! In/output variables:
    class(atype_field), intent(in) :: field
    type(type_mesh),    intent(in) :: mesh
    logical                        :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_parent_mesh'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => field)
    class is (atype_field_mesh)
      res = test_tol_mesh( p%parent, mesh, p%parent%tol_dist)
    class default
      res = .false.
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_parent_mesh

  function is_Arakawa_grid( field, field_Arakawa_grid) result( res)

    ! In/output variables:
    class(atype_field),      intent(in) :: field
    type(type_Arakawa_grid), intent(in) :: field_Arakawa_grid
    logical                             :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_Arakawa_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => field)
    class is (atype_field_grid)
      res = p%Arakawa_grid == field_Arakawa_grid
    class is (atype_field_mesh)
      res = p%Arakawa_grid == field_Arakawa_grid
    class default
      res = .false.
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_Arakawa_grid

  function is_third_dimension( field, field_third_dimension) result( res)

    ! In/output variables:
    class(atype_field),         intent(in) :: field
    type(type_third_dimension), intent(in) :: field_third_dimension
    logical                                :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_third_dimension'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => field)
    class default
      call crash('class/type not implemented')
    class is (type_field_grid_logical_2D)
      res = .false.
    class is (type_field_grid_logical_3D)
      res = p%third_dimension == field_third_dimension
    class is (type_field_grid_int_2D)
      res = .false.
    class is (type_field_grid_int_3D)
      res = p%third_dimension == field_third_dimension
    class is (type_field_grid_dp_2D)
      res = .false.
    class is (type_field_grid_dp_3D)
      res = p%third_dimension == field_third_dimension
    class is (type_field_mesh_logical_2D)
      res = .false.
    class is (type_field_mesh_logical_3D)
      res = p%third_dimension == field_third_dimension
    class is (type_field_mesh_int_2D)
      res = .false.
    class is (type_field_mesh_int_3D)
      res = p%third_dimension == field_third_dimension
    class is (type_field_mesh_dp_2D)
      res = .false.
    class is (type_field_mesh_dp_3D)
      res = p%third_dimension == field_third_dimension
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_third_dimension

end module fields_basic