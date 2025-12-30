module fields_basic

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid, type_Arakawa_grid
  use fields_dimensions, only: third_dimension, type_third_dimension
  use tests_main, only: test_tol_grid, test_tol_mesh
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: &
    atype_field, atype_field_grid, atype_field_mesh, &
    type_field_grid_logical_2D, type_field_grid_logical_3D, &
    type_field_grid_int_2D, type_field_grid_int_3D, &
    type_field_grid_dp_2D, type_field_grid_dp_3D, &
    type_field_mesh_logical_2D, type_field_mesh_logical_3D, &
    type_field_mesh_int_2D, type_field_mesh_int_3D, &
    type_field_mesh_dp_2D, type_field_mesh_dp_3D, &
    type_field_collection, add_initialised_field_to_collection, find_field_by_name

  ! Abstract base type and mesh/grid-based types
  ! ============================================

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

    procedure(print_field_info_ifc), deferred :: print_field_info
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
  contains
    procedure :: print_field_info  => print_field_grid_info_logical_2D
  end type type_field_grid_logical_2D

  type, extends(atype_field_grid) :: type_field_grid_logical_3D
    type(type_third_dimension) :: third_dimension
    logical, pointer :: d(:,:)
  contains
    procedure :: print_field_info  => print_field_grid_info_logical_3D
  end type type_field_grid_logical_3D

  type, extends(atype_field_grid) :: type_field_grid_int_2D
    integer, pointer :: d(:)
  contains
    procedure :: print_field_info  => print_field_grid_info_int_2D
  end type type_field_grid_int_2D

  type, extends(atype_field_grid) :: type_field_grid_int_3D
    type(type_third_dimension) :: third_dimension
    integer, pointer :: d(:,:)
  contains
    procedure :: print_field_info  => print_field_grid_info_int_3D
  end type type_field_grid_int_3D

  type, extends(atype_field_grid) :: type_field_grid_dp_2D
    real(dp), pointer :: d(:)
  contains
    procedure :: print_field_info  => print_field_grid_info_dp_2D
  end type type_field_grid_dp_2D

  type, extends(atype_field_grid) :: type_field_grid_dp_3D
    type(type_third_dimension) :: third_dimension
    real(dp), pointer :: d(:,:)
  contains
    procedure :: print_field_info  => print_field_grid_info_dp_3D
  end type type_field_grid_dp_3D

  ! Mesh-based fields

  type, extends(atype_field_mesh) :: type_field_mesh_logical_2D
    logical, pointer :: d(:)
  contains
    procedure :: print_field_info  => print_field_mesh_info_logical_2D
  end type type_field_mesh_logical_2D

  type, extends(atype_field_mesh) :: type_field_mesh_logical_3D
    type(type_third_dimension) :: third_dimension
    logical, pointer :: d(:,:)
  contains
    procedure :: print_field_info  => print_field_mesh_info_logical_3D
  end type type_field_mesh_logical_3D

  type, extends(atype_field_mesh) :: type_field_mesh_int_2D
    integer, pointer :: d(:)
  contains
    procedure :: print_field_info  => print_field_mesh_info_int_2D
  end type type_field_mesh_int_2D

  type, extends(atype_field_mesh) :: type_field_mesh_int_3D
    type(type_third_dimension) :: third_dimension
    integer, pointer :: d(:,:)
  contains
    procedure :: print_field_info  => print_field_mesh_info_int_3D
  end type type_field_mesh_int_3D

  type, extends(atype_field_mesh) :: type_field_mesh_dp_2D
    real(dp), pointer :: d(:)
  contains
    procedure :: print_field_info  => print_field_mesh_info_dp_2D
  end type type_field_mesh_dp_2D

  type, extends(atype_field_mesh) :: type_field_mesh_dp_3D
    type(type_third_dimension) :: third_dimension
    real(dp), pointer :: d(:,:)
  contains
    procedure :: print_field_info  => print_field_mesh_info_dp_3D
  end type type_field_mesh_dp_3D

  ! Wrapper so we can have a mixed-type array
  ! =========================================

  type :: type_field_box
     class(atype_field), allocatable :: p
  end type type_field_box

  ! Field collection
  ! ================

  type :: type_field_collection
     type(type_field_box), allocatable :: items(:)
     integer                           :: n     = 0
     integer                           :: n_max = 0
   contains
     procedure :: add_initialised_field_to_collection
     procedure :: print_fields_info
     procedure :: find_field_by_name
  end type type_field_collection

  ! Interfaces for deferred procedures
  ! ==================================

  abstract interface

    subroutine print_field_info_ifc( self)
      import :: atype_field
      class(atype_field), intent(in) :: self
    end subroutine print_field_info_ifc

  end interface

  ! Interfaces to procedures defined in submodules
  ! ==============================================

  ! print_field_info
  interface

    module subroutine print_field_grid_info_logical_2D( self)
      class(type_field_grid_logical_2D), intent(in) :: self
    end subroutine print_field_grid_info_logical_2D

    module subroutine print_field_grid_info_logical_3D( self)
      class(type_field_grid_logical_3D), intent(in) :: self
    end subroutine print_field_grid_info_logical_3D

    module subroutine print_field_grid_info_int_2D( self)
      class(type_field_grid_int_2D), intent(in) :: self
    end subroutine print_field_grid_info_int_2D

    module subroutine print_field_grid_info_int_3D( self)
      class(type_field_grid_int_3D), intent(in) :: self
    end subroutine print_field_grid_info_int_3D

    module subroutine print_field_grid_info_dp_2D( self)
      class(type_field_grid_dp_2D), intent(in) :: self
    end subroutine print_field_grid_info_dp_2D

    module subroutine print_field_grid_info_dp_3D( self)
      class(type_field_grid_dp_3D), intent(in) :: self
    end subroutine print_field_grid_info_dp_3D

    module subroutine print_field_mesh_info_logical_2D( self)
      class(type_field_mesh_logical_2D), intent(in) :: self
    end subroutine print_field_mesh_info_logical_2D

    module subroutine print_field_mesh_info_logical_3D( self)
      class(type_field_mesh_logical_3D), intent(in) :: self
    end subroutine print_field_mesh_info_logical_3D

    module subroutine print_field_mesh_info_int_2D( self)
      class(type_field_mesh_int_2D), intent(in) :: self
    end subroutine print_field_mesh_info_int_2D

    module subroutine print_field_mesh_info_int_3D( self)
      class(type_field_mesh_int_3D), intent(in) :: self
    end subroutine print_field_mesh_info_int_3D

    module subroutine print_field_mesh_info_dp_2D( self)
      class(type_field_mesh_dp_2D), intent(in) :: self
    end subroutine print_field_mesh_info_dp_2D

    module subroutine print_field_mesh_info_dp_3D( self)
      class(type_field_mesh_dp_3D), intent(in) :: self
    end subroutine print_field_mesh_info_dp_3D

  end interface

contains

  subroutine add_initialised_field_to_collection( self, field)

    ! In/output variables:
    class(type_field_collection), intent(inout) :: self
    class(atype_field),           intent(in   ) :: field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_initialised_field_to_collection'
    logical                        :: is_in_use
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    ! If the collection is not alocated yet, allocate it.
    if (.not. allocated( self%items)) then
      self%n     = 1
      self%n_max = 1
      allocate( self%items(1))
      allocate( self%items(1)%p, source = field)
      call finalise_routine( routine_name)
      return
    end if

    ! If the collection is full, extend it
    if (self%n == self%n_max) then
      self%n_max = self%n_max * 2
      call extend_field_collection( self%items, self%n, self%n_max)
    end if

    ! Check that this name is not already in use
    is_in_use = .false.
    do i = 1, self%n
      is_in_use = is_in_use .or. self%items(i)%p%name == field%name
    end do
    if (is_in_use) call crash('a field of name "' // trim( field%name) // '" already exists')

    ! Add field to collection
    self%n = self%n+1
    allocate( self%items( self%n)%p, source = field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine add_initialised_field_to_collection

  subroutine extend_field_collection( arr, n_keep, n_max_new)

    ! In/output variables:
    type(type_field_box), allocatable, intent(inout) :: arr(:)
    integer,                           intent(in   ) :: n_keep, n_max_new

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'extend_field_collection'
    type(type_field_box), allocatable :: tmp(:)
    integer                           :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    allocate( tmp( n_keep))
    do i = 1, n_keep
      allocate( tmp(i)%p, source = arr(i)%p)
    end do

    deallocate( arr)
    allocate( arr( n_max_new))

    do i = 1, n_keep
      allocate( arr(i)%p, source = tmp(i)%p)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine extend_field_collection

  subroutine print_fields_info( self)

    ! In/output variables:
    class(type_field_collection), intent(in) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'print_fields_info'
    integer                        :: i

    if (.not. allocated( self%items)) return

    ! Add routine to call stack
    call init_routine( routine_name)

    do i = 1, self%n
       call self%items(i)%p%print_field_info
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine print_fields_info

  function find_field_by_name( self, name) result(i)

    ! In/output variables:
    class(type_field_collection), intent(in) :: self
    character(len=*),             intent(in) :: name
    integer                                  :: i

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'find_field_by_name'
    logical                        :: found_it
    integer                        :: ii

    ! Add routine to call stack
    call init_routine( routine_name)

    if (.not. allocated(self%items)) then
      call crash('field collection not allocated')
    end if

    i = 0
    found_it = .false.
    do ii = 1, size( self%items)
      if (allocated( self%items(ii)%p)) then
        if (self%items(ii)%p%name == name) then
          found_it = .true.
          i = ii
          exit
        end if
      end if
    end do

    if (.not. found_it) then
      call crash('could not find field "' // trim( name) // '" in field collection')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function find_field_by_name

  function is_parent( self, parent) result( res)

    ! In/output variables:
    class(atype_field), intent(in) :: self
    class(*),           intent(in) :: parent
    logical                        :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_parent'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => parent)
    class is (type_grid)
      res = is_parent_grid( self, p)
    class is (type_mesh)
      res = is_parent_mesh( self, p)
    class default
      res = .false.
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_parent

  function is_parent_grid( self, grid) result( res)

    ! In/output variables:
    class(atype_field), intent(in) :: self
    type(type_grid),    intent(in) :: grid
    logical                        :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_parent_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => self)
    class is (atype_field_grid)
      res = test_tol_grid( p%parent, grid, p%parent%tol_dist)
    class default
      res = .false.
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_parent_grid

  function is_parent_mesh( self, mesh) result( res)

    ! In/output variables:
    class(atype_field), intent(in) :: self
    type(type_mesh),    intent(in) :: mesh
    logical                        :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_parent_mesh'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => self)
    class is (atype_field_mesh)
      res = test_tol_mesh( p%parent, mesh, p%parent%tol_dist)
    class default
      res = .false.
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_parent_mesh

  function is_Arakawa_grid( self, field_Arakawa_grid) result( res)

    ! In/output variables:
    class(atype_field),      intent(in) :: self
    type(type_Arakawa_grid), intent(in) :: field_Arakawa_grid
    logical                             :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_Arakawa_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => self)
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

  function is_third_dimension( self, field_third_dimension) result( res)

    ! In/output variables:
    class(atype_field),         intent(in) :: self
    type(type_third_dimension), intent(in) :: field_third_dimension
    logical                                :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_third_dimension'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => self)
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