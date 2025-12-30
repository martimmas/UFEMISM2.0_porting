module fields_init_field

  use precisions, only: dp
  use parameters, only: NaN
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: type_Arakawa_grid, Arakawa_grid
  use fields_dimensions, only: type_third_dimension
  use mpi_f08, only: MPI_WIN
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use fields_basic, only: atype_field, atype_field_grid, atype_field_mesh, &
    type_field_grid_logical_2D, type_field_grid_logical_3D, &
    type_field_grid_int_2D, type_field_grid_int_3D, &
    type_field_grid_dp_2D, type_field_grid_dp_3D, &
    type_field_mesh_logical_2D, type_field_mesh_logical_3D, &
    type_field_mesh_int_2D, type_field_mesh_int_3D, &
    type_field_mesh_dp_2D, type_field_mesh_dp_3D

  implicit none

  private

  public :: init_field

  interface init_field
    procedure :: init_field_grid_logical_2D
    procedure :: init_field_grid_logical_3D
    procedure :: init_field_grid_int_2D
    procedure :: init_field_grid_int_3D
    procedure :: init_field_grid_dp_2D
    procedure :: init_field_grid_dp_3D
    procedure :: init_field_mesh_logical_2D
    procedure :: init_field_mesh_logical_3D
    procedure :: init_field_mesh_int_2D
    procedure :: init_field_mesh_int_3D
    procedure :: init_field_mesh_dp_2D
    procedure :: init_field_mesh_dp_3D
  end interface init_field

contains

  ! Grid-based fields
  ! =================

  subroutine init_field_grid_logical_2D( field, d, w, grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_grid_logical_2D),           intent(  out) :: field
    logical, dimension(:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                              intent(inout) :: w
    type(type_grid), target,                    intent(in   ) :: grid
    type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
    character(len=*),                           intent(in   ) :: name
    character(len=*),                           intent(in   ) :: long_name
    character(len=*),                           intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'init_field_grid_logical_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call set_field_metadata( field, name, long_name, units)
    call set_field_parent_grid( field, grid, field_Arakawa_grid)

    ! Allocate memory for field data array and bind field pointer to it
    if (field_Arakawa_grid == Arakawa_grid%a()) then
      call allocate_dist_shared( d, w, grid%n_loc)
      d      ( grid%n1: grid%n2) => d
      field%d( grid%n1: grid%n2) => d
    else
      call crash('staggered square grids are not supported')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine init_field_grid_logical_2D

  subroutine init_field_grid_logical_3D( field, d, w, grid, field_third_dimension, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_grid_logical_3D),             intent(  out) :: field
    logical, dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                intent(inout) :: w
    type(type_grid), target,                      intent(in   ) :: grid
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    character(len=*),                             intent(in   ) :: name
    character(len=*),                             intent(in   ) :: long_name
    character(len=*),                             intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'init_field_grid_logical_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call set_field_metadata( field, name, long_name, units)
    call set_field_parent_grid( field, grid, field_Arakawa_grid)

    field%third_dimension = field_third_dimension

    ! Allocate memory for field data array and bind field pointer to it
    if (field_Arakawa_grid == Arakawa_grid%a()) then
      call allocate_dist_shared( d, w, grid%n_loc, field_third_dimension%n)
      d      ( grid%n1: grid%n2, 1: field_third_dimension%n) => d
      field%d( grid%n1: grid%n2, 1: field_third_dimension%n) => d
    else
      call crash('staggered square grids are not supported')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine init_field_grid_logical_3D

  subroutine init_field_grid_int_2D( field, d, w, grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_grid_int_2D),               intent(  out) :: field
    integer, dimension(:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                              intent(inout) :: w
    type(type_grid), target,                    intent(in   ) :: grid
    type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
    character(len=*),                           intent(in   ) :: name
    character(len=*),                           intent(in   ) :: long_name
    character(len=*),                           intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'init_field_grid_int_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call set_field_metadata( field, name, long_name, units)
    call set_field_parent_grid( field, grid, field_Arakawa_grid)

    ! Allocate memory for field data array and bind field pointer to it
    if (field_Arakawa_grid == Arakawa_grid%a()) then
      call allocate_dist_shared( d, w, grid%n_loc)
      d      ( grid%n1: grid%n2) => d
      field%d( grid%n1: grid%n2) => d
    else
      call crash('staggered square grids are not supported')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine init_field_grid_int_2D

  subroutine init_field_grid_int_3D( field, d, w, grid, field_third_dimension, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_grid_int_3D),                 intent(  out) :: field
    integer, dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                intent(inout) :: w
    type(type_grid), target,                      intent(in   ) :: grid
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    character(len=*),                             intent(in   ) :: name
    character(len=*),                             intent(in   ) :: long_name
    character(len=*),                             intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'init_field_grid_int_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call set_field_metadata( field, name, long_name, units)
    call set_field_parent_grid( field, grid, field_Arakawa_grid)

    field%third_dimension = field_third_dimension

    ! Allocate memory for field data array and bind field pointer to it
    if (field_Arakawa_grid == Arakawa_grid%a()) then
      call allocate_dist_shared( d, w, grid%n_loc, field_third_dimension%n)
      d      ( grid%n1: grid%n2, 1: field_third_dimension%n) => d
      field%d( grid%n1: grid%n2, 1: field_third_dimension%n) => d
    else
      call crash('staggered square grids are not supported')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine init_field_grid_int_3D

  subroutine init_field_grid_dp_2D( field, d, w, grid, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_grid_dp_2D),                 intent(  out) :: field
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                               intent(inout) :: w
    type(type_grid), target,                     intent(in   ) :: grid
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                            intent(in   ) :: name
    character(len=*),                            intent(in   ) :: long_name
    character(len=*),                            intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'init_field_grid_dp_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call set_field_metadata( field, name, long_name, units)
    call set_field_parent_grid( field, grid, field_Arakawa_grid)

    ! Allocate memory for field data array and bind field pointer to it
    if (field_Arakawa_grid == Arakawa_grid%a()) then
      call allocate_dist_shared( d, w, grid%n_loc)
      d      ( grid%n1: grid%n2) => d
      field%d( grid%n1: grid%n2) => d
    else
      call crash('staggered square grids are not supported')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine init_field_grid_dp_2D

  subroutine init_field_grid_dp_3D( field, d, w, grid, field_third_dimension, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_grid_dp_3D),                   intent(  out) :: field
    real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                 intent(inout) :: w
    type(type_grid), target,                       intent(in   ) :: grid
    type(type_third_dimension),                    intent(in   ) :: field_third_dimension
    type(type_Arakawa_grid),                       intent(in   ) :: field_Arakawa_grid
    character(len=*),                              intent(in   ) :: name
    character(len=*),                              intent(in   ) :: long_name
    character(len=*),                              intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'init_field_grid_dp_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call set_field_metadata( field, name, long_name, units)
    call set_field_parent_grid( field, grid, field_Arakawa_grid)

    field%third_dimension = field_third_dimension

    ! Allocate memory for field data array and bind field pointer to it
    if (field_Arakawa_grid == Arakawa_grid%a()) then
      call allocate_dist_shared( d, w, grid%n_loc, field_third_dimension%n)
      d      ( grid%n1: grid%n2, 1: field_third_dimension%n) => d
      field%d( grid%n1: grid%n2, 1: field_third_dimension%n) => d
    else
      call crash('staggered square grids are not supported')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine init_field_grid_dp_3D

  ! Mesh-based fields
  ! =================

  subroutine init_field_mesh_logical_2D( field, d, w, mesh, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_mesh_logical_2D),           intent(  out) :: field
    logical, dimension(:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                              intent(inout) :: w
    type(type_mesh), target,                    intent(in   ) :: mesh
    type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
    character(len=*),                           intent(in   ) :: name
    character(len=*),                           intent(in   ) :: long_name
    character(len=*),                           intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'init_field_mesh_logical_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call set_field_metadata( field, name, long_name, units)
    call set_field_parent_mesh( field, mesh, field_Arakawa_grid)

    ! Allocate memory for field data array and bind field pointer to it
    if (field_Arakawa_grid == Arakawa_grid%a()) then
      call allocate_dist_shared( d, w, mesh%nV_loc)
      d      ( mesh%vi1: mesh%vi2) => d
      field%d( mesh%vi1: mesh%vi2) => d
    elseif (field_Arakawa_grid == Arakawa_grid%b()) then
      call allocate_dist_shared( d, w, mesh%nTri_loc)
      d      ( mesh%ti1: mesh%ti2) => d
      field%d( mesh%ti1: mesh%ti2) => d
    elseif (field_Arakawa_grid == Arakawa_grid%c()) then
      call allocate_dist_shared( d, w, mesh%nE_loc)
      d      ( mesh%ei1: mesh%ei2) => d
      field%d( mesh%ei1: mesh%ei2) => d
    else
      call crash('invalid Arakawa grid')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine init_field_mesh_logical_2D

  subroutine init_field_mesh_logical_3D( field, d, w, mesh, field_third_dimension, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_mesh_logical_3D),             intent(  out) :: field
    logical, dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                intent(inout) :: w
    type(type_mesh), target,                      intent(in   ) :: mesh
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    character(len=*),                             intent(in   ) :: name
    character(len=*),                             intent(in   ) :: long_name
    character(len=*),                             intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'init_field_mesh_logical_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call set_field_metadata( field, name, long_name, units)
    call set_field_parent_mesh( field, mesh, field_Arakawa_grid)

    field%third_dimension = field_third_dimension

    ! Allocate memory for field data array and bind field pointer to it
    if (field_Arakawa_grid == Arakawa_grid%a()) then
      call allocate_dist_shared( d, w, mesh%nV_loc, field_third_dimension%n)
      d      ( mesh%vi1: mesh%vi2, 1: field_third_dimension%n) => d
      field%d( mesh%vi1: mesh%vi2, 1: field_third_dimension%n) => d
    elseif (field_Arakawa_grid == Arakawa_grid%b()) then
      call allocate_dist_shared( d, w, mesh%nTri_loc, field_third_dimension%n)
      d      ( mesh%ti1: mesh%ti2, 1: field_third_dimension%n) => d
      field%d( mesh%ti1: mesh%ti2, 1: field_third_dimension%n) => d
    elseif (field_Arakawa_grid == Arakawa_grid%c()) then
      call allocate_dist_shared( d, w, mesh%nE_loc, field_third_dimension%n)
      d      ( mesh%ei1: mesh%ei2, 1: field_third_dimension%n) => d
      field%d( mesh%ei1: mesh%ei2, 1: field_third_dimension%n) => d
    else
      call crash('invalid Arakawa grid')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine init_field_mesh_logical_3D

  subroutine init_field_mesh_int_2D( field, d, w, mesh, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_mesh_int_2D),               intent(  out) :: field
    integer, dimension(:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                              intent(inout) :: w
    type(type_mesh), target,                    intent(in   ) :: mesh
    type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
    character(len=*),                           intent(in   ) :: name
    character(len=*),                           intent(in   ) :: long_name
    character(len=*),                           intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'init_field_mesh_int_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call set_field_metadata( field, name, long_name, units)
    call set_field_parent_mesh( field, mesh, field_Arakawa_grid)

    ! Allocate memory for field data array and bind field pointer to it
    if (field_Arakawa_grid == Arakawa_grid%a()) then
      call allocate_dist_shared( d, w, mesh%nV_loc)
      d      ( mesh%vi1: mesh%vi2) => d
      field%d( mesh%vi1: mesh%vi2) => d
    elseif (field_Arakawa_grid == Arakawa_grid%b()) then
      call allocate_dist_shared( d, w, mesh%nTri_loc)
      d      ( mesh%ti1: mesh%ti2) => d
      field%d( mesh%ti1: mesh%ti2) => d
    elseif (field_Arakawa_grid == Arakawa_grid%c()) then
      call allocate_dist_shared( d, w, mesh%nE_loc)
      d      ( mesh%ei1: mesh%ei2) => d
      field%d( mesh%ei1: mesh%ei2) => d
    else
      call crash('invalid Arakawa grid')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine init_field_mesh_int_2D

  subroutine init_field_mesh_int_3D( field, d, w, mesh, field_third_dimension, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_mesh_int_3D),                 intent(  out) :: field
    integer, dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                intent(inout) :: w
    type(type_mesh), target,                      intent(in   ) :: mesh
    type(type_third_dimension),                   intent(in   ) :: field_third_dimension
    type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
    character(len=*),                             intent(in   ) :: name
    character(len=*),                             intent(in   ) :: long_name
    character(len=*),                             intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'init_field_mesh_int_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call set_field_metadata( field, name, long_name, units)
    call set_field_parent_mesh( field, mesh, field_Arakawa_grid)

    field%third_dimension = field_third_dimension

    ! Allocate memory for field data array and bind field pointer to it
    if (field_Arakawa_grid == Arakawa_grid%a()) then
      call allocate_dist_shared( d, w, mesh%nV_loc, field_third_dimension%n)
      d      ( mesh%vi1: mesh%vi2, 1: field_third_dimension%n) => d
      field%d( mesh%vi1: mesh%vi2, 1: field_third_dimension%n) => d
    elseif (field_Arakawa_grid == Arakawa_grid%b()) then
      call allocate_dist_shared( d, w, mesh%nTri_loc, field_third_dimension%n)
      d      ( mesh%ti1: mesh%ti2, 1: field_third_dimension%n) => d
      field%d( mesh%ti1: mesh%ti2, 1: field_third_dimension%n) => d
    elseif (field_Arakawa_grid == Arakawa_grid%c()) then
      call allocate_dist_shared( d, w, mesh%nE_loc, field_third_dimension%n)
      d      ( mesh%ei1: mesh%ei2, 1: field_third_dimension%n) => d
      field%d( mesh%ei1: mesh%ei2, 1: field_third_dimension%n) => d
    else
      call crash('invalid Arakawa grid')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine init_field_mesh_int_3D

  subroutine init_field_mesh_dp_2D( field, d, w, mesh, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_mesh_dp_2D),                 intent(  out) :: field
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                               intent(inout) :: w
    type(type_mesh), target,                     intent(in   ) :: mesh
    type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
    character(len=*),                            intent(in   ) :: name
    character(len=*),                            intent(in   ) :: long_name
    character(len=*),                            intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'init_field_mesh_dp_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call set_field_metadata( field, name, long_name, units)
    call set_field_parent_mesh( field, mesh, field_Arakawa_grid)

    ! Allocate memory for field data array and bind field pointer to it
    if (field_Arakawa_grid == Arakawa_grid%a()) then
      call allocate_dist_shared( d, w, mesh%nV_loc)
      d      ( mesh%vi1: mesh%vi2) => d
      field%d( mesh%vi1: mesh%vi2) => d
    elseif (field_Arakawa_grid == Arakawa_grid%b()) then
      call allocate_dist_shared( d, w, mesh%nTri_loc)
      d      ( mesh%ti1: mesh%ti2) => d
      field%d( mesh%ti1: mesh%ti2) => d
    elseif (field_Arakawa_grid == Arakawa_grid%c()) then
      call allocate_dist_shared( d, w, mesh%nE_loc)
      d      ( mesh%ei1: mesh%ei2) => d
      field%d( mesh%ei1: mesh%ei2) => d
    else
      call crash('invalid Arakawa grid')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine init_field_mesh_dp_2D

  subroutine init_field_mesh_dp_3D( field, d, w, mesh, field_third_dimension, &
    field_Arakawa_grid, name, long_name, units)

    ! In/output variables:
    type(type_field_mesh_dp_3D),                   intent(  out) :: field
    real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d
    type(MPI_WIN),                                 intent(inout) :: w
    type(type_mesh), target,                       intent(in   ) :: mesh
    type(type_third_dimension),                    intent(in   ) :: field_third_dimension
    type(type_Arakawa_grid),                       intent(in   ) :: field_Arakawa_grid
    character(len=*),                              intent(in   ) :: name
    character(len=*),                              intent(in   ) :: long_name
    character(len=*),                              intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'init_field_mesh_dp_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call set_field_metadata( field, name, long_name, units)
    call set_field_parent_mesh( field, mesh, field_Arakawa_grid)

    field%third_dimension = field_third_dimension

    ! Allocate memory for field data array and bind field pointer to it
    if (field_Arakawa_grid == Arakawa_grid%a()) then
      call allocate_dist_shared( d, w, mesh%nV_loc, field_third_dimension%n)
      d      ( mesh%vi1: mesh%vi2, 1: field_third_dimension%n) => d
      field%d( mesh%vi1: mesh%vi2, 1: field_third_dimension%n) => d
    elseif (field_Arakawa_grid == Arakawa_grid%b()) then
      call allocate_dist_shared( d, w, mesh%nTri_loc, field_third_dimension%n)
      d      ( mesh%ti1: mesh%ti2, 1: field_third_dimension%n) => d
      field%d( mesh%ti1: mesh%ti2, 1: field_third_dimension%n) => d
    elseif (field_Arakawa_grid == Arakawa_grid%c()) then
      call allocate_dist_shared( d, w, mesh%nE_loc, field_third_dimension%n)
      d      ( mesh%ei1: mesh%ei2, 1: field_third_dimension%n) => d
      field%d( mesh%ei1: mesh%ei2, 1: field_third_dimension%n) => d
    else
      call crash('invalid Arakawa grid')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine init_field_mesh_dp_3D

  ! Support functions
  ! =================

  subroutine set_field_metadata( field, name, long_name, units)

    ! In/output variables:
    class(atype_field), intent(inout) :: field
    character(len=*),   intent(in   ) :: name
    character(len=*),   intent(in   ) :: long_name
    character(len=*),   intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_field_metadata'

    ! Add routine to call stack
    call init_routine( routine_name)

    field%name      = name
    field%long_name = long_name
    field%units     = units

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_field_metadata

  subroutine set_field_parent_grid( field, grid, field_Arakawa_grid)

    ! In/output variables:
    class(atype_field_grid), intent(inout) :: field
    type(type_grid), target, intent(in   ) :: grid
    type(type_Arakawa_grid), intent(in   ) :: field_Arakawa_grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_field_parent_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    field%parent       => grid
    field%Arakawa_grid = field_Arakawa_grid

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_field_parent_grid

  subroutine set_field_parent_mesh( field, mesh, field_Arakawa_grid)

    ! In/output variables:
    class(atype_field_mesh), intent(inout) :: field
    type(type_mesh), target, intent(in   ) :: mesh
    type(type_Arakawa_grid), intent(in   ) :: field_Arakawa_grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_field_parent_mesh'

    ! Add routine to call stack
    call init_routine( routine_name)

    field%parent       => mesh
    field%Arakawa_grid = field_Arakawa_grid

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_field_parent_mesh

end module fields_init_field