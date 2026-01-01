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

  ! Abstract basic field type
  ! =========================

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

    procedure :: print_info
    procedure :: is_parent

  end type atype_field

  ! Interfaces to type-bound procedures defined in submodules
  interface
    module subroutine print_info( field)
      class(atype_field), intent(in) :: field
    end subroutine print_info
    module function is_parent( field, parent) result( res)
      class(atype_field), intent(in) :: field
      class(*),           intent(in) :: parent
      logical                        :: res
    end function is_parent
  end interface

  ! Extedned abstract field type for grid- and mesh-based fields
  ! ============================================================

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

end module fields_basic