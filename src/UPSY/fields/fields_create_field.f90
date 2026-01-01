module fields_create_field

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: type_Arakawa_grid
  use fields_dimensions, only: type_third_dimension
  use mpi_f08, only: MPI_WIN
  use fields_basic, only: &
    type_field_grid_logical_2D, type_field_grid_logical_3D, &
    type_field_grid_int_2D, type_field_grid_int_3D, &
    type_field_grid_dp_2D, type_field_grid_dp_3D, &
    type_field_mesh_logical_2D, type_field_mesh_logical_3D, &
    type_field_mesh_int_2D, type_field_mesh_int_3D, &
    type_field_mesh_dp_2D, type_field_mesh_dp_3D
  use fields_init_field, only: init_field
  use fields_field_collection, only: type_field_collection

  implicit none

  private

  public :: create_field

  interface create_field
    procedure :: create_field_grid_logical_2D
    procedure :: create_field_grid_logical_3D
    procedure :: create_field_grid_int_2D
    procedure :: create_field_grid_int_3D
    procedure :: create_field_grid_dp_2D
    procedure :: create_field_grid_dp_3D
    procedure :: create_field_mesh_logical_2D
    procedure :: create_field_mesh_logical_3D
    procedure :: create_field_mesh_int_2D
    procedure :: create_field_mesh_int_3D
    procedure :: create_field_mesh_dp_2D
    procedure :: create_field_mesh_dp_3D
  end interface create_field

contains

  ! Grid-based fields
  ! =================

  subroutine create_field_grid_logical_2D( bof, d, w, grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_collection),                intent(inout) :: bof
    logical, dimension(:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                              intent(inout) :: w
    type(type_grid), target,                    intent(in   ) :: grid
    type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
    character(len=*), optional,                 intent(in   ) :: name
    character(len=*), optional,                 intent(in   ) :: long_name
    character(len=*), optional,                 intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'create_field_grid_logical_2D'
    type(type_field_grid_logical_2D), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_non_optional_optional_arguments( name, long_name, units)
    call init_field( field, d, w, grid, &
      field_Arakawa_grid, name, long_name, units)
    call bof%add_initialised_field_to_collection( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_grid_logical_2D

  subroutine create_field_grid_logical_3D( bof, d, w, grid, field_third_dimension, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_collection),                  intent(inout) :: bof
    logical, dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                intent(inout) :: w
    type(type_grid), target,                      intent(in   ) :: grid
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    character(len=*), optional,                   intent(in   ) :: name
    character(len=*), optional,                   intent(in   ) :: long_name
    character(len=*), optional,                   intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'create_field_grid_logical_3D'
    type(type_field_grid_logical_3D), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_non_optional_optional_arguments( name, long_name, units)
    call init_field( field, d, w, grid, field_third_dimension, &
      field_Arakawa_grid, name, long_name, units)
    call bof%add_initialised_field_to_collection( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_grid_logical_3D

  subroutine create_field_grid_int_2D( bof, d, w, grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_collection),                intent(inout) :: bof
    integer, dimension(:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                              intent(inout) :: w
    type(type_grid), target,                    intent(in   ) :: grid
    type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
    character(len=*), optional,                 intent(in   ) :: name
    character(len=*), optional,                 intent(in   ) :: long_name
    character(len=*), optional,                 intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'create_field_grid_int_2D'
    type(type_field_grid_int_2D), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_non_optional_optional_arguments( name, long_name, units)
    call init_field( field, d, w, grid, &
      field_Arakawa_grid, name, long_name, units)
    call bof%add_initialised_field_to_collection( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_grid_int_2D

  subroutine create_field_grid_int_3D( bof, d, w, grid, field_third_dimension, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_collection),                  intent(inout) :: bof
    integer, dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                intent(inout) :: w
    type(type_grid), target,                      intent(in   ) :: grid
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    character(len=*), optional,                   intent(in   ) :: name
    character(len=*), optional,                   intent(in   ) :: long_name
    character(len=*), optional,                   intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'create_field_grid_int_3D'
    type(type_field_grid_int_3D), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_non_optional_optional_arguments( name, long_name, units)
    call init_field( field, d, w, grid, field_third_dimension, &
      field_Arakawa_grid, name, long_name, units)
    call bof%add_initialised_field_to_collection( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_grid_int_3D

  subroutine create_field_grid_dp_2D( bof, d, w, grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_collection),                 intent(inout) :: bof
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                               intent(inout) :: w
    type(type_grid), target,                     intent(in   ) :: grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*), optional,                  intent(in   ) :: name
    character(len=*), optional,                  intent(in   ) :: long_name
    character(len=*), optional,                  intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'create_field_grid_dp_2D'
    type(type_field_grid_dp_2D), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_non_optional_optional_arguments( name, long_name, units)
    call init_field( field, d, w, grid, &
      field_Arakawa_grid, name, long_name, units)
    call bof%add_initialised_field_to_collection( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_grid_dp_2D

  subroutine create_field_grid_dp_3D( bof, d, w, grid, field_third_dimension, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_collection),                   intent(inout) :: bof
    real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                 intent(inout) :: w
    type(type_grid), target,                       intent(in   ) :: grid
    type(type_third_dimension),                    intent(in   ) :: field_third_dimension
    type(type_Arakawa_grid),                       intent(in   ) :: field_Arakawa_grid
    character(len=*), optional,                    intent(in   ) :: name
    character(len=*), optional,                    intent(in   ) :: long_name
    character(len=*), optional,                    intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'create_field_grid_dp_3D'
    type(type_field_grid_dp_3D), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_non_optional_optional_arguments( name, long_name, units)
    call init_field( field, d, w, grid, field_third_dimension, &
      field_Arakawa_grid, name, long_name, units)
    call bof%add_initialised_field_to_collection( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_grid_dp_3D

  ! Mesh-based fields
  ! =================

  subroutine create_field_mesh_logical_2D( bof, d, w, mesh, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_collection),                intent(inout) :: bof
    logical, dimension(:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                              intent(inout) :: w
    type(type_mesh), target,                    intent(in   ) :: mesh
    type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
    character(len=*), optional,                 intent(in   ) :: name
    character(len=*), optional,                 intent(in   ) :: long_name
    character(len=*), optional,                 intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'create_field_mesh_logical_2D'
    type(type_field_mesh_logical_2D), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_non_optional_optional_arguments( name, long_name, units)
    call init_field( field, d, w, mesh, &
      field_Arakawa_grid, name, long_name, units)
    call bof%add_initialised_field_to_collection( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_mesh_logical_2D

  subroutine create_field_mesh_logical_3D( bof, d, w, mesh, field_third_dimension, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_collection),                  intent(inout) :: bof
    logical, dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                intent(inout) :: w
    type(type_mesh), target,                      intent(in   ) :: mesh
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    character(len=*), optional,                   intent(in   ) :: name
    character(len=*), optional,                   intent(in   ) :: long_name
    character(len=*), optional,                   intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'create_field_mesh_logical_3D'
    type(type_field_mesh_logical_3D), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_non_optional_optional_arguments( name, long_name, units)
    call init_field( field, d, w, mesh, field_third_dimension, &
      field_Arakawa_grid, name, long_name, units)
    call bof%add_initialised_field_to_collection( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_mesh_logical_3D

  subroutine create_field_mesh_int_2D( bof, d, w, mesh, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_collection),                intent(inout) :: bof
    integer, dimension(:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                              intent(inout) :: w
    type(type_mesh), target,                    intent(in   ) :: mesh
    type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
    character(len=*), optional,                 intent(in   ) :: name
    character(len=*), optional,                 intent(in   ) :: long_name
    character(len=*), optional,                 intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'create_field_mesh_int_2D'
    type(type_field_mesh_int_2D), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_non_optional_optional_arguments( name, long_name, units)
    call init_field( field, d, w, mesh, &
      field_Arakawa_grid, name, long_name, units)
    call bof%add_initialised_field_to_collection( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_mesh_int_2D

  subroutine create_field_mesh_int_3D( bof, d, w, mesh, field_third_dimension, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_collection),                  intent(inout) :: bof
    integer, dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                intent(inout) :: w
    type(type_mesh), target,                      intent(in   ) :: mesh
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    character(len=*), optional,                   intent(in   ) :: name
    character(len=*), optional,                   intent(in   ) :: long_name
    character(len=*), optional,                   intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'create_field_mesh_int_3D'
    type(type_field_mesh_int_3D), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_non_optional_optional_arguments( name, long_name, units)
    call init_field( field, d, w, mesh, field_third_dimension, &
      field_Arakawa_grid, name, long_name, units)
    call bof%add_initialised_field_to_collection( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_mesh_int_3D

  subroutine create_field_mesh_dp_2D( bof, d, w, mesh, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_collection),                 intent(inout) :: bof
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                               intent(inout) :: w
    type(type_mesh), target,                     intent(in   ) :: mesh
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*), optional,                  intent(in   ) :: name
    character(len=*), optional,                  intent(in   ) :: long_name
    character(len=*), optional,                  intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'create_field_mesh_dp_2D'
    type(type_field_mesh_dp_2D), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_non_optional_optional_arguments( name, long_name, units)
    call init_field( field, d, w, mesh, &
      field_Arakawa_grid, name, long_name, units)
    call bof%add_initialised_field_to_collection( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_mesh_dp_2D

  subroutine create_field_mesh_dp_3D( bof, d, w, mesh, field_third_dimension, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_collection),                   intent(inout) :: bof
    real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                 intent(inout) :: w
    type(type_mesh), target,                       intent(in   ) :: mesh
    type(type_third_dimension),                    intent(in   ) :: field_third_dimension
    type(type_Arakawa_grid),                       intent(in   ) :: field_Arakawa_grid
    character(len=*), optional,                    intent(in   ) :: name
    character(len=*), optional,                    intent(in   ) :: long_name
    character(len=*), optional,                    intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter   :: routine_name = 'create_field_mesh_dp_3D'
    type(type_field_mesh_dp_3D), allocatable :: field

    ! Add routine to call stack
    call init_routine( routine_name)

    call check_non_optional_optional_arguments( name, long_name, units)
    call init_field( field, d, w, mesh, field_third_dimension, &
      field_Arakawa_grid, name, long_name, units)
    call bof%add_initialised_field_to_collection( field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_field_mesh_dp_3D

  ! Support functions
  ! =================

  subroutine check_non_optional_optional_arguments( name, long_name, units)

    character(len=*), optional, intent(in) :: name
    character(len=*), optional, intent(in) :: long_name
    character(len=*), optional, intent(in) :: units

    ! The optional arguments are actually non-optional; it's just that
    ! this approach results in more readable code on the call site.
    if (.not. present( name     )) call crash('input variable "name" not provided')
    if (.not. present( long_name)) call crash('input variable "long_name" not provided')
    if (.not. present( units    )) call crash('input variable "units" not provided')

  end subroutine check_non_optional_optional_arguments

end module fields_create_field