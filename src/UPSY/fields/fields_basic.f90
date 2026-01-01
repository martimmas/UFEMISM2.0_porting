module fields_basic

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, colour_string
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid, type_Arakawa_grid
  use fields_dimensions, only: third_dimension, type_third_dimension
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
    type_field_mesh_dp_2D, type_field_mesh_dp_3D

  ! Abstract basic field type
  ! =========================

  type, abstract :: atype_field

    ! Metadata
    character(len=1024), private :: name_val
    character(len=1024), private :: long_name_val
    character(len=1024), private :: units_val

    ! Parent grid
    ! type(type_grid), pointer :: parent    ! Defined in extended types atype_field_grid and atype_field_mesh!
    ! type(type_mesh), pointer :: parent
    type(type_Arakawa_grid),    private :: Arakawa_grid_val
    type(type_third_dimension), private :: third_dimension_val

    ! Pointer to array containing the actual field data
    ! logical, pointer :: d(:)              ! Defined in extended types type_field_grid_XXX, etc.
    type(MPI_WIN), public :: w

  contains

    procedure, public :: set_metadata
    procedure, public :: name      => get_name
    procedure, public :: long_name => get_long_name
    procedure, public :: units     => get_units

    procedure, public :: set_parent_Arakawa_grid
    procedure, public :: Arakawa_grid => get_parent_Arakawa_grid

    procedure, public :: set_third_dimension
    procedure, public :: third_dimension => get_third_dimension

    procedure, public :: print_info
    procedure, public :: is_parent

  end type atype_field

  ! Interfaces to type-bound procedures defined in submodules
  interface

    module subroutine set_metadata( field, name, long_name, units)
      class(atype_field), intent(inout) :: field
      character(len=*),   intent(in   ) :: name
      character(len=*),   intent(in   ) :: long_name
      character(len=*),   intent(in   ) :: units
    end subroutine set_metadata

    module function get_name( field) result( name)
      class(atype_field), intent(in) :: field
      character(:), allocatable      :: name
    end function get_name

    module function get_long_name( field) result( long_name)
      class(atype_field), intent(in) :: field
      character(:), allocatable      :: long_name
    end function get_long_name

    module function get_units( field) result( units)
      class(atype_field), intent(in) :: field
      character(:), allocatable      :: units
    end function get_units

    module subroutine print_info( field)
      class(atype_field), intent(in) :: field
    end subroutine print_info

    module function is_parent( field, parent) result( res)
      class(atype_field), intent(in) :: field
      class(*),           intent(in) :: parent
      logical                        :: res
    end function is_parent

    module subroutine set_parent_Arakawa_grid( field, field_Arakawa_grid)
      class(atype_field),      intent(inout) :: field
      type(type_Arakawa_grid), intent(in   ) :: field_Arakawa_grid
    end subroutine set_parent_Arakawa_grid

    module function get_parent_Arakawa_grid( field) result( field_Arakawa_grid)
      class(atype_field), intent(in) :: field
      type(type_Arakawa_grid)        :: field_Arakawa_grid
    end function get_parent_Arakawa_grid

    module subroutine set_third_dimension( field, field_third_dimension)
      class(atype_field),         intent(inout) :: field
      type(type_third_dimension), intent(in   ) :: field_third_dimension
    end subroutine set_third_dimension

    module function get_third_dimension( field) result( field_third_dimension)
      class(atype_field), intent(in) :: field
      type(type_third_dimension)     :: field_third_dimension
    end function get_third_dimension

  end interface

  ! Extedned abstract field type for grid- and mesh-based fields
  ! ============================================================

  type, abstract, extends(atype_field) :: atype_field_grid
    type(type_grid), pointer, public :: parent
  contains
    procedure, public :: set_parent_grid
  end type atype_field_grid

  type, abstract, extends(atype_field) :: atype_field_mesh
    type(type_mesh), pointer, public :: parent
  contains
    procedure, public :: set_parent_mesh
  end type atype_field_mesh

  ! Interfaces to type-bound procedures defined in submodules
  interface

    module subroutine set_parent_grid( field, grid)
      class(atype_field_grid), intent(inout) :: field
      type(type_grid), target, intent(in   ) :: grid
    end subroutine set_parent_grid

    module subroutine set_parent_mesh( field, mesh)
      class(atype_field_mesh), intent(inout) :: field
      type(type_mesh), target, intent(in   ) :: mesh
    end subroutine set_parent_mesh

  end interface

  ! Concrete field types
  ! ====================

  ! Grid-based fields

  type, extends(atype_field_grid) :: type_field_grid_logical_2D
    logical, pointer, public :: d(:)
  end type type_field_grid_logical_2D

  type, extends(atype_field_grid) :: type_field_grid_logical_3D
    logical, pointer, public :: d(:,:)
  end type type_field_grid_logical_3D

  type, extends(atype_field_grid) :: type_field_grid_int_2D
    integer, pointer, public :: d(:)
  end type type_field_grid_int_2D

  type, extends(atype_field_grid) :: type_field_grid_int_3D
    integer, pointer, public :: d(:,:)
  end type type_field_grid_int_3D

  type, extends(atype_field_grid) :: type_field_grid_dp_2D
    real(dp), pointer, public :: d(:)
  end type type_field_grid_dp_2D

  type, extends(atype_field_grid) :: type_field_grid_dp_3D
    real(dp), pointer, public :: d(:,:)
  end type type_field_grid_dp_3D

  ! Mesh-based fields

  type, extends(atype_field_mesh) :: type_field_mesh_logical_2D
    logical, pointer, public :: d(:)
  end type type_field_mesh_logical_2D

  type, extends(atype_field_mesh) :: type_field_mesh_logical_3D
    logical, pointer, public :: d(:,:)
  end type type_field_mesh_logical_3D

  type, extends(atype_field_mesh) :: type_field_mesh_int_2D
    integer, pointer, public :: d(:)
  end type type_field_mesh_int_2D

  type, extends(atype_field_mesh) :: type_field_mesh_int_3D
    integer, pointer, public :: d(:,:)
  end type type_field_mesh_int_3D

  type, extends(atype_field_mesh) :: type_field_mesh_dp_2D
    real(dp), pointer, public :: d(:)
  end type type_field_mesh_dp_2D

  type, extends(atype_field_mesh) :: type_field_mesh_dp_3D
    real(dp), pointer, public :: d(:,:)
  end type type_field_mesh_dp_3D

end module fields_basic