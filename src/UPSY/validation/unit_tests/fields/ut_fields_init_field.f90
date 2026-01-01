module ut_fields_init_field

  use precisions, only: dp
  use mpi_basic, only: par
  use parameters, only: pi
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use fields_main
  use ut_basic
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use grid_basic, only: setup_square_grid
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_dimensions, only: third_dimension
  use mpi_f08, only: MPI_WIN
  use mesh_memory, only: allocate_mesh_primary, crop_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_secondary, only: calc_all_secondary_mesh_data

  implicit none

  private

  public :: test_init_field

contains

  subroutine test_init_field( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_init_field'
    character(len=1024), parameter :: test_name_local = 'init_field'
    character(len=1024)            :: test_name
    real(dp)                       :: alpha_min, res_max
    real(dp), parameter            :: xmin = -1._dp
    real(dp), parameter            :: xmax =  1._dp
    real(dp), parameter            :: ymin = -1._dp
    real(dp), parameter            :: ymax =  1._dp
    type(type_mesh)                :: mesh

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call allocate_mesh_primary( mesh, 'test_mesh', 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 23.2_dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call crop_mesh_primary( mesh)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)

    call test_init_field_grid_logical  ( test_name)
    call test_init_field_grid_int      ( test_name)
    call test_init_field_grid_dp       ( test_name)
    call test_init_field_mesh_logical_a( test_name, mesh)
    call test_init_field_mesh_logical_b( test_name, mesh)
    call test_init_field_mesh_logical_c( test_name, mesh)
    call test_init_field_mesh_int_a    ( test_name, mesh)
    call test_init_field_mesh_int_b    ( test_name, mesh)
    call test_init_field_mesh_int_c    ( test_name, mesh)
    call test_init_field_mesh_dp_a     ( test_name, mesh)
    call test_init_field_mesh_dp_b     ( test_name, mesh)
    call test_init_field_mesh_dp_c     ( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_init_field

  subroutine test_init_field_grid_logical( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_init_field_grid_logical'
    character(len=1024), parameter                :: test_name_local = 'grid/logical'
    character(len=1024)                           :: test_name
    type(type_grid), target                       :: grid
    character(len=1024)                           :: name, long_name, units
    integer                                       :: nz
    type(type_field_grid_logical_2D), allocatable :: field_grid_2D
    type(type_field_grid_logical_3D), allocatable :: field_grid_3D_zeta
    type(type_field_grid_logical_3D), allocatable :: field_grid_3D_month
    type(type_field_grid_logical_3D), allocatable :: field_grid_3D_ocean
    logical, dimension(:  ), contiguous, pointer  :: d_grid_2D
    logical, dimension(:,:), contiguous, pointer  :: d_grid_3D_zeta
    logical, dimension(:,:), contiguous, pointer  :: d_grid_3D_month
    logical, dimension(:,:), contiguous, pointer  :: d_grid_3D_ocean
    type(MPI_WIN)                                 :: wd_grid_2D
    type(MPI_WIN)                                 :: wd_grid_3D_zeta
    type(MPI_WIN)                                 :: wd_grid_3D_month
    type(MPI_WIN)                                 :: wd_grid_3D_ocean
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call setup_square_grid( 'dummy_grid', 0._dp, 1._dp, 0._dp, 1._dp, 0.1_dp, grid)

    ! 2-D

    name      = 'd_grid_2D'
    long_name = 'd_grid_2D_long_name'
    units     = 'd_grid_2D_units'

    call init_field( field_grid_2D, d_grid_2D, wd_grid_2D, &
      grid, Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_grid_2D,1)
    ub1_a = ubound( d_grid_2D,1)

    lb1_f = lbound( field_grid_2D%d,1)
    ub1_f = ubound( field_grid_2D%d,1)

    d_grid_2D( grid%n1+1) = .true.

    call unit_test( (&
      field_grid_2D%name()      == name .and. &
      field_grid_2D%long_name() == long_name .and. &
      field_grid_2D%units()     == units .and. &
      field_grid_2D%is_parent( grid) .and. &
      field_grid_2D%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      field_grid_2D%d( grid%n1+1) .eqv. .true.), &
      trim( test_name) // '_2D')

    ! 3-D (ice zeta)

    name      = 'd_grid_3D_zeta'
    long_name = 'd_grid_3D_zeta_long_name'
    units     = 'd_grid_3D_zeta_units'
    nz        = 10

    call init_field( field_grid_3D_zeta, d_grid_3D_zeta, wd_grid_3D_zeta, &
      grid, third_dimension%ice_zeta( nz), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_grid_3D_zeta,1)
    ub1_a = ubound( d_grid_3D_zeta,1)

    lb1_f = lbound( field_grid_3D_zeta%d,1)
    ub1_f = ubound( field_grid_3D_zeta%d,1)

    lb2_a = lbound( d_grid_3D_zeta,2)
    ub2_a = ubound( d_grid_3D_zeta,2)

    lb2_f = lbound( field_grid_3D_zeta%d,2)
    ub2_f = ubound( field_grid_3D_zeta%d,2)

    d_grid_3D_zeta( grid%n1+1, 3) = .true.

    call unit_test( (&
      field_grid_3D_zeta%name()      == name .and. &
      field_grid_3D_zeta%long_name() == long_name .and. &
      field_grid_3D_zeta%units()     == units .and. &
      field_grid_3D_zeta%is_parent( grid) .and. &
      field_grid_3D_zeta%is_parent( third_dimension%ice_zeta( nz)) .and. &
      field_grid_3D_zeta%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_grid_3D_zeta%d( grid%n1+1, 3) .eqv. .true.), &
      trim( test_name) // '_3D_zeta')

    ! 3-D (month)

    name      = 'd_grid_3D_month'
    long_name = 'd_grid_3D_month_long_name'
    units     = 'd_grid_3D_month_units'

    call init_field( field_grid_3D_month, d_grid_3D_month, wd_grid_3D_month, &
      grid, third_dimension%month(), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_grid_3D_month,1)
    ub1_a = ubound( d_grid_3D_month,1)

    lb1_f = lbound( field_grid_3D_month%d,1)
    ub1_f = ubound( field_grid_3D_month%d,1)

    lb2_a = lbound( d_grid_3D_month,2)
    ub2_a = ubound( d_grid_3D_month,2)

    lb2_f = lbound( field_grid_3D_month%d,2)
    ub2_f = ubound( field_grid_3D_month%d,2)

    d_grid_3D_month( grid%n1+1, 3) = .true.

    call unit_test( (&
      field_grid_3D_month%name()      == name .and. &
      field_grid_3D_month%long_name() == long_name .and. &
      field_grid_3D_month%units()     == units .and. &
      field_grid_3D_month%is_parent( grid) .and. &
      field_grid_3D_month%is_parent( third_dimension%month()) .and. &
      field_grid_3D_month%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      field_grid_3D_month%d( grid%n1+1, 3) .eqv. .true.), &
      trim( test_name) // '_3D_month')

    ! 3-D (ocean depth)

    name      = 'd_grid_3D_ocean'
    long_name = 'd_grid_3D_ocean_long_name'
    units     = 'd_grid_3D_ocean_units'
    nz        = 20

    call init_field( field_grid_3D_ocean, d_grid_3D_ocean, wd_grid_3D_ocean, &
      grid, third_dimension%ocean_depth( nz), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_grid_3D_ocean,1)
    ub1_a = ubound( d_grid_3D_ocean,1)

    lb1_f = lbound( field_grid_3D_ocean%d,1)
    ub1_f = ubound( field_grid_3D_ocean%d,1)

    lb2_a = lbound( d_grid_3D_ocean,2)
    ub2_a = ubound( d_grid_3D_ocean,2)

    lb2_f = lbound( field_grid_3D_ocean%d,2)
    ub2_f = ubound( field_grid_3D_ocean%d,2)

    d_grid_3D_ocean( grid%n1+1, 3) = .true.

    call unit_test( (&
      field_grid_3D_ocean%name()      == name .and. &
      field_grid_3D_ocean%long_name() == long_name .and. &
      field_grid_3D_ocean%units()     == units .and. &
      field_grid_3D_ocean%is_parent( grid) .and. &
      field_grid_3D_ocean%is_parent( third_dimension%ocean_depth( nz)) .and. &
      field_grid_3D_ocean%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_grid_3D_ocean%d( grid%n1+1, 3) .eqv. .true.), &
      trim( test_name) // '_3D_ocean')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_init_field_grid_logical

  subroutine test_init_field_grid_int( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_init_field_grid_int'
    character(len=1024), parameter               :: test_name_local = 'grid/int'
    character(len=1024)                          :: test_name
    type(type_grid), target                      :: grid
    character(len=1024)                          :: name, long_name, units
    integer                                      :: nz
    type(type_field_grid_int_2D), allocatable    :: field_grid_2D
    type(type_field_grid_int_3D), allocatable    :: field_grid_3D_zeta
    type(type_field_grid_int_3D), allocatable    :: field_grid_3D_month
    type(type_field_grid_int_3D), allocatable    :: field_grid_3D_ocean
    integer, dimension(:  ), contiguous, pointer :: d_grid_2D
    integer, dimension(:,:), contiguous, pointer :: d_grid_3D_zeta
    integer, dimension(:,:), contiguous, pointer :: d_grid_3D_month
    integer, dimension(:,:), contiguous, pointer :: d_grid_3D_ocean
    type(MPI_WIN)                                :: wd_grid_2D
    type(MPI_WIN)                                :: wd_grid_3D_zeta
    type(MPI_WIN)                                :: wd_grid_3D_month
    type(MPI_WIN)                                :: wd_grid_3D_ocean
    integer                                      :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                      :: lb1_f, ub1_f, lb2_f, ub2_f
    integer, parameter                           :: test_val = 1337

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call setup_square_grid( 'dummy_grid', 0._dp, 1._dp, 0._dp, 1._dp, 0.1_dp, grid)

    ! 2-D

    name      = 'd_grid_2D'
    long_name = 'd_grid_2D_long_name'
    units     = 'd_grid_2D_units'

    call init_field( field_grid_2D, d_grid_2D, wd_grid_2D, &
      grid, Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_grid_2D,1)
    ub1_a = ubound( d_grid_2D,1)

    lb1_f = lbound( field_grid_2D%d,1)
    ub1_f = ubound( field_grid_2D%d,1)

    d_grid_2D( grid%n1+1) = test_val

    call unit_test( (&
      field_grid_2D%name()      == name .and. &
      field_grid_2D%long_name() == long_name .and. &
      field_grid_2D%units()     == units .and. &
      field_grid_2D%is_parent( grid) .and. &
      field_grid_2D%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      field_grid_2D%d( grid%n1+1) == test_val), &
      trim( test_name) // '_2D')

    ! 3-D (ice zeta)

    name      = 'd_grid_3D_zeta'
    long_name = 'd_grid_3D_zeta_long_name'
    units     = 'd_grid_3D_zeta_units'
    nz        = 10

    call init_field( field_grid_3D_zeta, d_grid_3D_zeta, wd_grid_3D_zeta, &
      grid, third_dimension%ice_zeta( nz), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_grid_3D_zeta,1)
    ub1_a = ubound( d_grid_3D_zeta,1)

    lb1_f = lbound( field_grid_3D_zeta%d,1)
    ub1_f = ubound( field_grid_3D_zeta%d,1)

    lb2_a = lbound( d_grid_3D_zeta,2)
    ub2_a = ubound( d_grid_3D_zeta,2)

    lb2_f = lbound( field_grid_3D_zeta%d,2)
    ub2_f = ubound( field_grid_3D_zeta%d,2)

    d_grid_3D_zeta( grid%n1+1, 3) = test_val

    call unit_test( (&
      field_grid_3D_zeta%name()      == name .and. &
      field_grid_3D_zeta%long_name() == long_name .and. &
      field_grid_3D_zeta%units()     == units .and. &
      field_grid_3D_zeta%is_parent( grid) .and. &
      field_grid_3D_zeta%is_parent( third_dimension%ice_zeta( nz)) .and. &
      field_grid_3D_zeta%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_grid_3D_zeta%d( grid%n1+1, 3) == test_val), &
      trim( test_name) // '_3D_zeta')

    ! 3-D (month)

    name      = 'd_grid_3D_month'
    long_name = 'd_grid_3D_month_long_name'
    units     = 'd_grid_3D_month_units'

    call init_field( field_grid_3D_month, d_grid_3D_month, wd_grid_3D_month, &
      grid, third_dimension%month(), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_grid_3D_month,1)
    ub1_a = ubound( d_grid_3D_month,1)

    lb1_f = lbound( field_grid_3D_month%d,1)
    ub1_f = ubound( field_grid_3D_month%d,1)

    lb2_a = lbound( d_grid_3D_month,2)
    ub2_a = ubound( d_grid_3D_month,2)

    lb2_f = lbound( field_grid_3D_month%d,2)
    ub2_f = ubound( field_grid_3D_month%d,2)

    d_grid_3D_month( grid%n1+1, 3) = test_val

    call unit_test( (&
      field_grid_3D_month%name()      == name .and. &
      field_grid_3D_month%long_name() == long_name .and. &
      field_grid_3D_month%units()     == units .and. &
      field_grid_3D_month%is_parent( grid) .and. &
      field_grid_3D_month%is_parent( third_dimension%month()) .and. &
      field_grid_3D_month%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      field_grid_3D_month%d( grid%n1+1, 3) == test_val), &
      trim( test_name) // '_3D_month')

    ! 3-D (ocean depth)

    name      = 'd_grid_3D_ocean'
    long_name = 'd_grid_3D_ocean_long_name'
    units     = 'd_grid_3D_ocean_units'
    nz        = 20

    call init_field( field_grid_3D_ocean, d_grid_3D_ocean, wd_grid_3D_ocean, &
      grid, third_dimension%ocean_depth( nz), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_grid_3D_ocean,1)
    ub1_a = ubound( d_grid_3D_ocean,1)

    lb1_f = lbound( field_grid_3D_ocean%d,1)
    ub1_f = ubound( field_grid_3D_ocean%d,1)

    lb2_a = lbound( d_grid_3D_ocean,2)
    ub2_a = ubound( d_grid_3D_ocean,2)

    lb2_f = lbound( field_grid_3D_ocean%d,2)
    ub2_f = ubound( field_grid_3D_ocean%d,2)

    d_grid_3D_ocean( grid%n1+1, 3) = test_val

    call unit_test( (&
      field_grid_3D_ocean%name()      == name .and. &
      field_grid_3D_ocean%long_name() == long_name .and. &
      field_grid_3D_ocean%units()     == units .and. &
      field_grid_3D_ocean%is_parent( grid) .and. &
      field_grid_3D_ocean%is_parent( third_dimension%ocean_depth( nz)) .and. &
      field_grid_3D_ocean%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_grid_3D_ocean%d( grid%n1+1, 3) == test_val), &
      trim( test_name) // '_3D_ocean')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_init_field_grid_int

  subroutine test_init_field_grid_dp( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_init_field_grid_dp'
    character(len=1024), parameter                :: test_name_local = 'grid/dp'
    character(len=1024)                           :: test_name
    type(type_grid), target                       :: grid
    character(len=1024)                           :: name, long_name, units
    integer                                       :: nz
    type(type_field_grid_dp_2D), allocatable      :: field_grid_2D
    type(type_field_grid_dp_3D), allocatable      :: field_grid_3D_zeta
    type(type_field_grid_dp_3D), allocatable      :: field_grid_3D_month
    type(type_field_grid_dp_3D), allocatable      :: field_grid_3D_ocean
    real(dp), dimension(:  ), contiguous, pointer :: d_grid_2D
    real(dp), dimension(:,:), contiguous, pointer :: d_grid_3D_zeta
    real(dp), dimension(:,:), contiguous, pointer :: d_grid_3D_month
    real(dp), dimension(:,:), contiguous, pointer :: d_grid_3D_ocean
    type(MPI_WIN)                                 :: wd_grid_2D
    type(MPI_WIN)                                 :: wd_grid_3D_zeta
    type(MPI_WIN)                                 :: wd_grid_3D_month
    type(MPI_WIN)                                 :: wd_grid_3D_ocean
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f
    real(dp), parameter                           :: test_val = 13.37_dp

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call setup_square_grid( 'dummy_grid', 0._dp, 1._dp, 0._dp, 1._dp, 0.1_dp, grid)

    ! 2-D

    name      = 'd_grid_2D'
    long_name = 'd_grid_2D_long_name'
    units     = 'd_grid_2D_units'

    call init_field( field_grid_2D, d_grid_2D, wd_grid_2D, &
      grid, Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_grid_2D,1)
    ub1_a = ubound( d_grid_2D,1)

    lb1_f = lbound( field_grid_2D%d,1)
    ub1_f = ubound( field_grid_2D%d,1)

    d_grid_2D( grid%n1+1) = test_val

    call unit_test( (&
      field_grid_2D%name()      == name .and. &
      field_grid_2D%long_name() == long_name .and. &
      field_grid_2D%units()     == units .and. &
      field_grid_2D%is_parent( grid) .and. &
      field_grid_2D%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      field_grid_2D%d( grid%n1+1) == test_val), &
      trim( test_name) // '_2D')

    ! 3-D (ice zeta)

    name      = 'd_grid_3D_zeta'
    long_name = 'd_grid_3D_zeta_long_name'
    units     = 'd_grid_3D_zeta_units'
    nz        = 10

    call init_field( field_grid_3D_zeta, d_grid_3D_zeta, wd_grid_3D_zeta, &
      grid, third_dimension%ice_zeta( nz), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_grid_3D_zeta,1)
    ub1_a = ubound( d_grid_3D_zeta,1)

    lb1_f = lbound( field_grid_3D_zeta%d,1)
    ub1_f = ubound( field_grid_3D_zeta%d,1)

    lb2_a = lbound( d_grid_3D_zeta,2)
    ub2_a = ubound( d_grid_3D_zeta,2)

    lb2_f = lbound( field_grid_3D_zeta%d,2)
    ub2_f = ubound( field_grid_3D_zeta%d,2)

    d_grid_3D_zeta( grid%n1+1, 3) = test_val

    call unit_test( (&
      field_grid_3D_zeta%name()      == name .and. &
      field_grid_3D_zeta%long_name() == long_name .and. &
      field_grid_3D_zeta%units()     == units .and. &
      field_grid_3D_zeta%is_parent( grid) .and. &
      field_grid_3D_zeta%is_parent( third_dimension%ice_zeta( nz)) .and. &
      field_grid_3D_zeta%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_grid_3D_zeta%d( grid%n1+1, 3) == test_val), &
      trim( test_name) // '_3D_zeta')

    ! 3-D (month)

    name      = 'd_grid_3D_month'
    long_name = 'd_grid_3D_month_long_name'
    units     = 'd_grid_3D_month_units'

    call init_field( field_grid_3D_month, d_grid_3D_month, wd_grid_3D_month, &
      grid, third_dimension%month(), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_grid_3D_month,1)
    ub1_a = ubound( d_grid_3D_month,1)

    lb1_f = lbound( field_grid_3D_month%d,1)
    ub1_f = ubound( field_grid_3D_month%d,1)

    lb2_a = lbound( d_grid_3D_month,2)
    ub2_a = ubound( d_grid_3D_month,2)

    lb2_f = lbound( field_grid_3D_month%d,2)
    ub2_f = ubound( field_grid_3D_month%d,2)

    d_grid_3D_month( grid%n1+1, 3) = test_val

    call unit_test( (&
      field_grid_3D_month%name()      == name .and. &
      field_grid_3D_month%long_name() == long_name .and. &
      field_grid_3D_month%units()     == units .and. &
      field_grid_3D_month%is_parent( grid) .and. &
      field_grid_3D_month%is_parent( third_dimension%month()) .and. &
      field_grid_3D_month%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      field_grid_3D_month%d( grid%n1+1, 3) == test_val), &
      trim( test_name) // '_3D_month')

    ! 3-D (ocean depth)

    name      = 'd_grid_3D_ocean'
    long_name = 'd_grid_3D_ocean_long_name'
    units     = 'd_grid_3D_ocean_units'
    nz        = 20

    call init_field( field_grid_3D_ocean, d_grid_3D_ocean, wd_grid_3D_ocean, &
      grid, third_dimension%ocean_depth( nz), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_grid_3D_ocean,1)
    ub1_a = ubound( d_grid_3D_ocean,1)

    lb1_f = lbound( field_grid_3D_ocean%d,1)
    ub1_f = ubound( field_grid_3D_ocean%d,1)

    lb2_a = lbound( d_grid_3D_ocean,2)
    ub2_a = ubound( d_grid_3D_ocean,2)

    lb2_f = lbound( field_grid_3D_ocean%d,2)
    ub2_f = ubound( field_grid_3D_ocean%d,2)

    d_grid_3D_ocean( grid%n1+1, 3) = test_val

    call unit_test( (&
      field_grid_3D_ocean%name()      == name .and. &
      field_grid_3D_ocean%long_name() == long_name .and. &
      field_grid_3D_ocean%units()     == units .and. &
      field_grid_3D_ocean%is_parent( grid) .and. &
      field_grid_3D_ocean%is_parent( third_dimension%ocean_depth( nz)) .and. &
      field_grid_3D_ocean%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_grid_3D_ocean%d( grid%n1+1, 3) == test_val), &
      trim( test_name) // '_3D_ocean')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_init_field_grid_dp

  subroutine test_init_field_mesh_logical_a( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),        intent(in) :: test_name_parent
    type(type_mesh), target, intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_init_field_mesh_logical_a'
    character(len=1024), parameter                :: test_name_local = 'mesh/a/logical'
    character(len=1024)                           :: test_name
    character(len=1024)                           :: name, long_name, units
    integer                                       :: nz
    type(type_field_mesh_logical_2D), allocatable :: field_mesh_2D
    type(type_field_mesh_logical_3D), allocatable :: field_mesh_3D_zeta
    type(type_field_mesh_logical_3D), allocatable :: field_mesh_3D_month
    type(type_field_mesh_logical_3D), allocatable :: field_mesh_3D_ocean
    logical, dimension(:  ), contiguous, pointer  :: d_mesh_2D
    logical, dimension(:,:), contiguous, pointer  :: d_mesh_3D_zeta
    logical, dimension(:,:), contiguous, pointer  :: d_mesh_3D_month
    logical, dimension(:,:), contiguous, pointer  :: d_mesh_3D_ocean
    type(MPI_WIN)                                 :: wd_mesh_2D
    type(MPI_WIN)                                 :: wd_mesh_3D_zeta
    type(MPI_WIN)                                 :: wd_mesh_3D_month
    type(MPI_WIN)                                 :: wd_mesh_3D_ocean
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! 2-D, a-mesh

    name      = 'd_mesh_2D'
    long_name = 'd_mesh_2D_long_name'
    units     = 'd_mesh_2D_units'

    call init_field( field_mesh_2D, d_mesh_2D, wd_mesh_2D, &
      mesh, Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_mesh_2D,1)
    ub1_a = ubound( d_mesh_2D,1)

    lb1_f = lbound( field_mesh_2D%d,1)
    ub1_f = ubound( field_mesh_2D%d,1)

    d_mesh_2D( mesh%vi1+1) = .true.

    call unit_test( (&
      field_mesh_2D%name()      == name .and. &
      field_mesh_2D%long_name() == long_name .and. &
      field_mesh_2D%units()     == units .and. &
      field_mesh_2D%is_parent( mesh) .and. &
      field_mesh_2D%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == mesh%vi1 .and. &
      ub1_a == mesh%vi2 .and. &
      lb1_f == mesh%vi1 .and. &
      ub1_f == mesh%vi2 .and. &
      field_mesh_2D%d( mesh%vi1+1) .eqv. .true.), &
      trim( test_name) // '_2D')

    ! 3-D (ice zeta), a-mesh

    name      = 'd_mesh_3D_zeta'
    long_name = 'd_mesh_3D_zeta_long_name'
    units     = 'd_mesh_3D_zeta_units'
    nz        = 10

    call init_field( field_mesh_3D_zeta, d_mesh_3D_zeta, wd_mesh_3D_zeta, &
      mesh, third_dimension%ice_zeta( nz), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_zeta,1)
    ub1_a = ubound( d_mesh_3D_zeta,1)

    lb1_f = lbound( field_mesh_3D_zeta%d,1)
    ub1_f = ubound( field_mesh_3D_zeta%d,1)

    lb2_a = lbound( d_mesh_3D_zeta,2)
    ub2_a = ubound( d_mesh_3D_zeta,2)

    lb2_f = lbound( field_mesh_3D_zeta%d,2)
    ub2_f = ubound( field_mesh_3D_zeta%d,2)

    d_mesh_3D_zeta( mesh%vi1+1, 3) = .true.

    call unit_test( (&
      field_mesh_3D_zeta%name()      == name .and. &
      field_mesh_3D_zeta%long_name() == long_name .and. &
      field_mesh_3D_zeta%units()     == units .and. &
      field_mesh_3D_zeta%is_parent( mesh) .and. &
      field_mesh_3D_zeta%is_parent( third_dimension%ice_zeta( nz)) .and. &
      field_mesh_3D_zeta%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == mesh%vi1 .and. &
      ub1_a == mesh%vi2 .and. &
      lb1_f == mesh%vi1 .and. &
      ub1_f == mesh%vi2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_zeta%d( mesh%vi1+1, 3) .eqv. .true.), &
      trim( test_name) // '_3D_zeta')

    ! 3-D (month), a-mesh

    name      = 'd_mesh_3D_month'
    long_name = 'd_mesh_3D_month_long_name'
    units     = 'd_mesh_3D_month_units'

    call init_field( field_mesh_3D_month, d_mesh_3D_month, wd_mesh_3D_month, &
      mesh, third_dimension%month(), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_month,1)
    ub1_a = ubound( d_mesh_3D_month,1)

    lb1_f = lbound( field_mesh_3D_month%d,1)
    ub1_f = ubound( field_mesh_3D_month%d,1)

    lb2_a = lbound( d_mesh_3D_month,2)
    ub2_a = ubound( d_mesh_3D_month,2)

    lb2_f = lbound( field_mesh_3D_month%d,2)
    ub2_f = ubound( field_mesh_3D_month%d,2)

    d_mesh_3D_month( mesh%vi1+1, 3) = .true.

    call unit_test( (&
      field_mesh_3D_month%name()      == name .and. &
      field_mesh_3D_month%long_name() == long_name .and. &
      field_mesh_3D_month%units()     == units .and. &
      field_mesh_3D_month%is_parent( mesh) .and. &
      field_mesh_3D_month%is_parent( third_dimension%month()) .and. &
      field_mesh_3D_month%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == mesh%vi1 .and. &
      ub1_a == mesh%vi2 .and. &
      lb1_f == mesh%vi1 .and. &
      ub1_f == mesh%vi2 .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      field_mesh_3D_month%d( mesh%vi1+1, 3) .eqv. .true.), &
      trim( test_name) // '_3D_month')

    ! 3-D (ocean depth), a-mesh

    name      = 'd_mesh_3D_ocean'
    long_name = 'd_mesh_3D_ocean_long_name'
    units     = 'd_mesh_3D_ocean_units'
    nz        = 20

    call init_field( field_mesh_3D_ocean, d_mesh_3D_ocean, wd_mesh_3D_ocean, &
      mesh, third_dimension%ocean_depth( nz), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_ocean,1)
    ub1_a = ubound( d_mesh_3D_ocean,1)

    lb1_f = lbound( field_mesh_3D_ocean%d,1)
    ub1_f = ubound( field_mesh_3D_ocean%d,1)

    lb2_a = lbound( d_mesh_3D_ocean,2)
    ub2_a = ubound( d_mesh_3D_ocean,2)

    lb2_f = lbound( field_mesh_3D_ocean%d,2)
    ub2_f = ubound( field_mesh_3D_ocean%d,2)

    d_mesh_3D_ocean( mesh%vi1+1, 3) = .true.

    call unit_test( (&
      field_mesh_3D_ocean%name()      == name .and. &
      field_mesh_3D_ocean%long_name() == long_name .and. &
      field_mesh_3D_ocean%units()     == units .and. &
      field_mesh_3D_ocean%is_parent( mesh) .and. &
      field_mesh_3D_ocean%is_parent( third_dimension%ocean_depth( nz)) .and. &
      field_mesh_3D_ocean%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == mesh%vi1 .and. &
      ub1_a == mesh%vi2 .and. &
      lb1_f == mesh%vi1 .and. &
      ub1_f == mesh%vi2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_ocean%d( mesh%vi1+1, 3) .eqv. .true.), &
      trim( test_name) // '_3D_ocean')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_init_field_mesh_logical_a

  subroutine test_init_field_mesh_logical_b( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),        intent(in) :: test_name_parent
    type(type_mesh), target, intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_init_field_mesh_logical_b'
    character(len=1024), parameter                :: test_name_local = 'mesh/b/logical'
    character(len=1024)                           :: test_name
    character(len=1024)                           :: name, long_name, units
    integer                                       :: nz
    type(type_field_mesh_logical_2D), allocatable :: field_mesh_2D
    type(type_field_mesh_logical_3D), allocatable :: field_mesh_3D_zeta
    type(type_field_mesh_logical_3D), allocatable :: field_mesh_3D_month
    type(type_field_mesh_logical_3D), allocatable :: field_mesh_3D_ocean
    logical, dimension(:  ), contiguous, pointer  :: d_mesh_2D
    logical, dimension(:,:), contiguous, pointer  :: d_mesh_3D_zeta
    logical, dimension(:,:), contiguous, pointer  :: d_mesh_3D_month
    logical, dimension(:,:), contiguous, pointer  :: d_mesh_3D_ocean
    type(MPI_WIN)                                 :: wd_mesh_2D
    type(MPI_WIN)                                 :: wd_mesh_3D_zeta
    type(MPI_WIN)                                 :: wd_mesh_3D_month
    type(MPI_WIN)                                 :: wd_mesh_3D_ocean
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! 2-D, a-mesh

    name      = 'd_mesh_2D'
    long_name = 'd_mesh_2D_long_name'
    units     = 'd_mesh_2D_units'

    call init_field( field_mesh_2D, d_mesh_2D, wd_mesh_2D, &
      mesh, Arakawa_grid%b(), name, long_name, units)

    lb1_a = lbound( d_mesh_2D,1)
    ub1_a = ubound( d_mesh_2D,1)

    lb1_f = lbound( field_mesh_2D%d,1)
    ub1_f = ubound( field_mesh_2D%d,1)

    d_mesh_2D( mesh%ti1+1) = .true.

    call unit_test( (&
      field_mesh_2D%name()      == name .and. &
      field_mesh_2D%long_name() == long_name .and. &
      field_mesh_2D%units()     == units .and. &
      field_mesh_2D%is_parent( mesh) .and. &
      field_mesh_2D%is_parent( Arakawa_grid%b()) .and. &
      lb1_a == mesh%ti1 .and. &
      ub1_a == mesh%ti2 .and. &
      lb1_f == mesh%ti1 .and. &
      ub1_f == mesh%ti2 .and. &
      field_mesh_2D%d( mesh%ti1+1) .eqv. .true.), &
      trim( test_name) // '_2D')

    ! 3-D (ice zeta), a-mesh

    name      = 'd_mesh_3D_zeta'
    long_name = 'd_mesh_3D_zeta_long_name'
    units     = 'd_mesh_3D_zeta_units'
    nz        = 10

    call init_field( field_mesh_3D_zeta, d_mesh_3D_zeta, wd_mesh_3D_zeta, &
      mesh, third_dimension%ice_zeta( nz), Arakawa_grid%b(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_zeta,1)
    ub1_a = ubound( d_mesh_3D_zeta,1)

    lb1_f = lbound( field_mesh_3D_zeta%d,1)
    ub1_f = ubound( field_mesh_3D_zeta%d,1)

    lb2_a = lbound( d_mesh_3D_zeta,2)
    ub2_a = ubound( d_mesh_3D_zeta,2)

    lb2_f = lbound( field_mesh_3D_zeta%d,2)
    ub2_f = ubound( field_mesh_3D_zeta%d,2)

    d_mesh_3D_zeta( mesh%ti1+1, 3) = .true.

    call unit_test( (&
      field_mesh_3D_zeta%name()      == name .and. &
      field_mesh_3D_zeta%long_name() == long_name .and. &
      field_mesh_3D_zeta%units()     == units .and. &
      field_mesh_3D_zeta%is_parent( mesh) .and. &
      field_mesh_3D_zeta%is_parent( third_dimension%ice_zeta( nz)) .and. &
      field_mesh_3D_zeta%is_parent( Arakawa_grid%b()) .and. &
      lb1_a == mesh%ti1 .and. &
      ub1_a == mesh%ti2 .and. &
      lb1_f == mesh%ti1 .and. &
      ub1_f == mesh%ti2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_zeta%d( mesh%ti1+1, 3) .eqv. .true.), &
      trim( test_name) // '_3D_zeta')

    ! 3-D (month), a-mesh

    name      = 'd_mesh_3D_month'
    long_name = 'd_mesh_3D_month_long_name'
    units     = 'd_mesh_3D_month_units'

    call init_field( field_mesh_3D_month, d_mesh_3D_month, wd_mesh_3D_month, &
      mesh, third_dimension%month(), Arakawa_grid%b(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_month,1)
    ub1_a = ubound( d_mesh_3D_month,1)

    lb1_f = lbound( field_mesh_3D_month%d,1)
    ub1_f = ubound( field_mesh_3D_month%d,1)

    lb2_a = lbound( d_mesh_3D_month,2)
    ub2_a = ubound( d_mesh_3D_month,2)

    lb2_f = lbound( field_mesh_3D_month%d,2)
    ub2_f = ubound( field_mesh_3D_month%d,2)

    d_mesh_3D_month( mesh%ti1+1, 3) = .true.

    call unit_test( (&
      field_mesh_3D_month%name()      == name .and. &
      field_mesh_3D_month%long_name() == long_name .and. &
      field_mesh_3D_month%units()     == units .and. &
      field_mesh_3D_month%is_parent( mesh) .and. &
      field_mesh_3D_month%is_parent( third_dimension%month()) .and. &
      field_mesh_3D_month%is_parent( Arakawa_grid%b()) .and. &
      lb1_a == mesh%ti1 .and. &
      ub1_a == mesh%ti2 .and. &
      lb1_f == mesh%ti1 .and. &
      ub1_f == mesh%ti2 .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      field_mesh_3D_month%d( mesh%ti1+1, 3) .eqv. .true.), &
      trim( test_name) // '_3D_month')

    ! 3-D (ocean depth), a-mesh

    name      = 'd_mesh_3D_ocean'
    long_name = 'd_mesh_3D_ocean_long_name'
    units     = 'd_mesh_3D_ocean_units'
    nz        = 20

    call init_field( field_mesh_3D_ocean, d_mesh_3D_ocean, wd_mesh_3D_ocean, &
      mesh, third_dimension%ocean_depth( nz), Arakawa_grid%b(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_ocean,1)
    ub1_a = ubound( d_mesh_3D_ocean,1)

    lb1_f = lbound( field_mesh_3D_ocean%d,1)
    ub1_f = ubound( field_mesh_3D_ocean%d,1)

    lb2_a = lbound( d_mesh_3D_ocean,2)
    ub2_a = ubound( d_mesh_3D_ocean,2)

    lb2_f = lbound( field_mesh_3D_ocean%d,2)
    ub2_f = ubound( field_mesh_3D_ocean%d,2)

    d_mesh_3D_ocean( mesh%ti1+1, 3) = .true.

    call unit_test( (&
      field_mesh_3D_ocean%name()      == name .and. &
      field_mesh_3D_ocean%long_name() == long_name .and. &
      field_mesh_3D_ocean%units()     == units .and. &
      field_mesh_3D_ocean%is_parent( mesh) .and. &
      field_mesh_3D_ocean%is_parent( third_dimension%ocean_depth( nz)) .and. &
      field_mesh_3D_ocean%is_parent( Arakawa_grid%b()) .and. &
      lb1_a == mesh%ti1 .and. &
      ub1_a == mesh%ti2 .and. &
      lb1_f == mesh%ti1 .and. &
      ub1_f == mesh%ti2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_ocean%d( mesh%ti1+1, 3) .eqv. .true.), &
      trim( test_name) // '_3D_ocean')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_init_field_mesh_logical_b

  subroutine test_init_field_mesh_logical_c( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),        intent(in) :: test_name_parent
    type(type_mesh), target, intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_init_field_mesh_logical_c'
    character(len=1024), parameter                :: test_name_local = 'mesh/c/logical'
    character(len=1024)                           :: test_name
    character(len=1024)                           :: name, long_name, units
    integer                                       :: nz
    type(type_field_mesh_logical_2D), allocatable :: field_mesh_2D
    type(type_field_mesh_logical_3D), allocatable :: field_mesh_3D_zeta
    type(type_field_mesh_logical_3D), allocatable :: field_mesh_3D_month
    type(type_field_mesh_logical_3D), allocatable :: field_mesh_3D_ocean
    logical, dimension(:  ), contiguous, pointer  :: d_mesh_2D
    logical, dimension(:,:), contiguous, pointer  :: d_mesh_3D_zeta
    logical, dimension(:,:), contiguous, pointer  :: d_mesh_3D_month
    logical, dimension(:,:), contiguous, pointer  :: d_mesh_3D_ocean
    type(MPI_WIN)                                 :: wd_mesh_2D
    type(MPI_WIN)                                 :: wd_mesh_3D_zeta
    type(MPI_WIN)                                 :: wd_mesh_3D_month
    type(MPI_WIN)                                 :: wd_mesh_3D_ocean
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! 2-D, a-mesh

    name      = 'd_mesh_2D'
    long_name = 'd_mesh_2D_long_name'
    units     = 'd_mesh_2D_units'

    call init_field( field_mesh_2D, d_mesh_2D, wd_mesh_2D, &
      mesh, Arakawa_grid%c(), name, long_name, units)

    lb1_a = lbound( d_mesh_2D,1)
    ub1_a = ubound( d_mesh_2D,1)

    lb1_f = lbound( field_mesh_2D%d,1)
    ub1_f = ubound( field_mesh_2D%d,1)

    d_mesh_2D( mesh%ei1+1) = .true.

    call unit_test( (&
      field_mesh_2D%name()      == name .and. &
      field_mesh_2D%long_name() == long_name .and. &
      field_mesh_2D%units()     == units .and. &
      field_mesh_2D%is_parent( mesh) .and. &
      field_mesh_2D%is_parent( Arakawa_grid%c()) .and. &
      lb1_a == mesh%ei1 .and. &
      ub1_a == mesh%ei2 .and. &
      lb1_f == mesh%ei1 .and. &
      ub1_f == mesh%ei2 .and. &
      field_mesh_2D%d( mesh%ei1+1) .eqv. .true.), &
      trim( test_name) // '_2D')

    ! 3-D (ice zeta), a-mesh

    name      = 'd_mesh_3D_zeta'
    long_name = 'd_mesh_3D_zeta_long_name'
    units     = 'd_mesh_3D_zeta_units'
    nz        = 10

    call init_field( field_mesh_3D_zeta, d_mesh_3D_zeta, wd_mesh_3D_zeta, &
      mesh, third_dimension%ice_zeta( nz), Arakawa_grid%c(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_zeta,1)
    ub1_a = ubound( d_mesh_3D_zeta,1)

    lb1_f = lbound( field_mesh_3D_zeta%d,1)
    ub1_f = ubound( field_mesh_3D_zeta%d,1)

    lb2_a = lbound( d_mesh_3D_zeta,2)
    ub2_a = ubound( d_mesh_3D_zeta,2)

    lb2_f = lbound( field_mesh_3D_zeta%d,2)
    ub2_f = ubound( field_mesh_3D_zeta%d,2)

    d_mesh_3D_zeta( mesh%ei1+1, 3) = .true.

    call unit_test( (&
      field_mesh_3D_zeta%name()      == name .and. &
      field_mesh_3D_zeta%long_name() == long_name .and. &
      field_mesh_3D_zeta%units()     == units .and. &
      field_mesh_3D_zeta%is_parent( mesh) .and. &
      field_mesh_3D_zeta%is_parent( third_dimension%ice_zeta( nz)) .and. &
      field_mesh_3D_zeta%is_parent( Arakawa_grid%c()) .and. &
      lb1_a == mesh%ei1 .and. &
      ub1_a == mesh%ei2 .and. &
      lb1_f == mesh%ei1 .and. &
      ub1_f == mesh%ei2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_zeta%d( mesh%ei1+1, 3) .eqv. .true.), &
      trim( test_name) // '_3D_zeta')

    ! 3-D (month), a-mesh

    name      = 'd_mesh_3D_month'
    long_name = 'd_mesh_3D_month_long_name'
    units     = 'd_mesh_3D_month_units'

    call init_field( field_mesh_3D_month, d_mesh_3D_month, wd_mesh_3D_month, &
      mesh, third_dimension%month(), Arakawa_grid%c(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_month,1)
    ub1_a = ubound( d_mesh_3D_month,1)

    lb1_f = lbound( field_mesh_3D_month%d,1)
    ub1_f = ubound( field_mesh_3D_month%d,1)

    lb2_a = lbound( d_mesh_3D_month,2)
    ub2_a = ubound( d_mesh_3D_month,2)

    lb2_f = lbound( field_mesh_3D_month%d,2)
    ub2_f = ubound( field_mesh_3D_month%d,2)

    d_mesh_3D_month( mesh%ei1+1, 3) = .true.

    call unit_test( (&
      field_mesh_3D_month%name()      == name .and. &
      field_mesh_3D_month%long_name() == long_name .and. &
      field_mesh_3D_month%units()     == units .and. &
      field_mesh_3D_month%is_parent( mesh) .and. &
      field_mesh_3D_month%is_parent( third_dimension%month()) .and. &
      field_mesh_3D_month%is_parent( Arakawa_grid%c()) .and. &
      lb1_a == mesh%ei1 .and. &
      ub1_a == mesh%ei2 .and. &
      lb1_f == mesh%ei1 .and. &
      ub1_f == mesh%ei2 .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      field_mesh_3D_month%d( mesh%ei1+1, 3) .eqv. .true.), &
      trim( test_name) // '_3D_month')

    ! 3-D (ocean depth), a-mesh

    name      = 'd_mesh_3D_ocean'
    long_name = 'd_mesh_3D_ocean_long_name'
    units     = 'd_mesh_3D_ocean_units'
    nz        = 20

    call init_field( field_mesh_3D_ocean, d_mesh_3D_ocean, wd_mesh_3D_ocean, &
      mesh, third_dimension%ocean_depth( nz), Arakawa_grid%c(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_ocean,1)
    ub1_a = ubound( d_mesh_3D_ocean,1)

    lb1_f = lbound( field_mesh_3D_ocean%d,1)
    ub1_f = ubound( field_mesh_3D_ocean%d,1)

    lb2_a = lbound( d_mesh_3D_ocean,2)
    ub2_a = ubound( d_mesh_3D_ocean,2)

    lb2_f = lbound( field_mesh_3D_ocean%d,2)
    ub2_f = ubound( field_mesh_3D_ocean%d,2)

    d_mesh_3D_ocean( mesh%ei1+1, 3) = .true.

    call unit_test( (&
      field_mesh_3D_ocean%name()      == name .and. &
      field_mesh_3D_ocean%long_name() == long_name .and. &
      field_mesh_3D_ocean%units()     == units .and. &
      field_mesh_3D_ocean%is_parent( mesh) .and. &
      field_mesh_3D_ocean%is_parent( third_dimension%ocean_depth( nz)) .and. &
      field_mesh_3D_ocean%is_parent( Arakawa_grid%c()) .and. &
      lb1_a == mesh%ei1 .and. &
      ub1_a == mesh%ei2 .and. &
      lb1_f == mesh%ei1 .and. &
      ub1_f == mesh%ei2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_ocean%d( mesh%ei1+1, 3) .eqv. .true.), &
      trim( test_name) // '_3D_ocean')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_init_field_mesh_logical_c

  subroutine test_init_field_mesh_int_a( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),        intent(in) :: test_name_parent
    type(type_mesh), target, intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_init_field_mesh_int_a'
    character(len=1024), parameter               :: test_name_local = 'mesh/a/int'
    character(len=1024)                          :: test_name
    character(len=1024)                          :: name, long_name, units
    integer                                      :: nz
    type(type_field_mesh_int_2D), allocatable    :: field_mesh_2D
    type(type_field_mesh_int_3D), allocatable    :: field_mesh_3D_zeta
    type(type_field_mesh_int_3D), allocatable    :: field_mesh_3D_month
    type(type_field_mesh_int_3D), allocatable    :: field_mesh_3D_ocean
    integer, dimension(:  ), contiguous, pointer :: d_mesh_2D
    integer, dimension(:,:), contiguous, pointer :: d_mesh_3D_zeta
    integer, dimension(:,:), contiguous, pointer :: d_mesh_3D_month
    integer, dimension(:,:), contiguous, pointer :: d_mesh_3D_ocean
    type(MPI_WIN)                                :: wd_mesh_2D
    type(MPI_WIN)                                :: wd_mesh_3D_zeta
    type(MPI_WIN)                                :: wd_mesh_3D_month
    type(MPI_WIN)                                :: wd_mesh_3D_ocean
    integer                                      :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                      :: lb1_f, ub1_f, lb2_f, ub2_f
    integer, parameter                           :: test_val = 1337

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! 2-D, a-mesh

    name      = 'd_mesh_2D'
    long_name = 'd_mesh_2D_long_name'
    units     = 'd_mesh_2D_units'

    call init_field( field_mesh_2D, d_mesh_2D, wd_mesh_2D, &
      mesh, Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_mesh_2D,1)
    ub1_a = ubound( d_mesh_2D,1)

    lb1_f = lbound( field_mesh_2D%d,1)
    ub1_f = ubound( field_mesh_2D%d,1)

    d_mesh_2D( mesh%vi1+1) = test_val

    call unit_test( (&
      field_mesh_2D%name()      == name .and. &
      field_mesh_2D%long_name() == long_name .and. &
      field_mesh_2D%units()     == units .and. &
      field_mesh_2D%is_parent( mesh) .and. &
      field_mesh_2D%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == mesh%vi1 .and. &
      ub1_a == mesh%vi2 .and. &
      lb1_f == mesh%vi1 .and. &
      ub1_f == mesh%vi2 .and. &
      field_mesh_2D%d( mesh%vi1+1) == test_val), &
      trim( test_name) // '_2D')

    ! 3-D (ice zeta), a-mesh

    name      = 'd_mesh_3D_zeta'
    long_name = 'd_mesh_3D_zeta_long_name'
    units     = 'd_mesh_3D_zeta_units'
    nz        = 10

    call init_field( field_mesh_3D_zeta, d_mesh_3D_zeta, wd_mesh_3D_zeta, &
      mesh, third_dimension%ice_zeta( nz), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_zeta,1)
    ub1_a = ubound( d_mesh_3D_zeta,1)

    lb1_f = lbound( field_mesh_3D_zeta%d,1)
    ub1_f = ubound( field_mesh_3D_zeta%d,1)

    lb2_a = lbound( d_mesh_3D_zeta,2)
    ub2_a = ubound( d_mesh_3D_zeta,2)

    lb2_f = lbound( field_mesh_3D_zeta%d,2)
    ub2_f = ubound( field_mesh_3D_zeta%d,2)

    d_mesh_3D_zeta( mesh%vi1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_zeta%name()      == name .and. &
      field_mesh_3D_zeta%long_name() == long_name .and. &
      field_mesh_3D_zeta%units()     == units .and. &
      field_mesh_3D_zeta%is_parent( mesh) .and. &
      field_mesh_3D_zeta%is_parent( third_dimension%ice_zeta( nz)) .and. &
      field_mesh_3D_zeta%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == mesh%vi1 .and. &
      ub1_a == mesh%vi2 .and. &
      lb1_f == mesh%vi1 .and. &
      ub1_f == mesh%vi2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_zeta%d( mesh%vi1+1, 3) == test_val), &
      trim( test_name) // '_3D_zeta')

    ! 3-D (month), a-mesh

    name      = 'd_mesh_3D_month'
    long_name = 'd_mesh_3D_month_long_name'
    units     = 'd_mesh_3D_month_units'

    call init_field( field_mesh_3D_month, d_mesh_3D_month, wd_mesh_3D_month, &
      mesh, third_dimension%month(), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_month,1)
    ub1_a = ubound( d_mesh_3D_month,1)

    lb1_f = lbound( field_mesh_3D_month%d,1)
    ub1_f = ubound( field_mesh_3D_month%d,1)

    lb2_a = lbound( d_mesh_3D_month,2)
    ub2_a = ubound( d_mesh_3D_month,2)

    lb2_f = lbound( field_mesh_3D_month%d,2)
    ub2_f = ubound( field_mesh_3D_month%d,2)

    d_mesh_3D_month( mesh%vi1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_month%name()      == name .and. &
      field_mesh_3D_month%long_name() == long_name .and. &
      field_mesh_3D_month%units()     == units .and. &
      field_mesh_3D_month%is_parent( mesh) .and. &
      field_mesh_3D_month%is_parent( third_dimension%month()) .and. &
      field_mesh_3D_month%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == mesh%vi1 .and. &
      ub1_a == mesh%vi2 .and. &
      lb1_f == mesh%vi1 .and. &
      ub1_f == mesh%vi2 .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      field_mesh_3D_month%d( mesh%vi1+1, 3) == test_val), &
      trim( test_name) // '_3D_month')

    ! 3-D (ocean depth), a-mesh

    name      = 'd_mesh_3D_ocean'
    long_name = 'd_mesh_3D_ocean_long_name'
    units     = 'd_mesh_3D_ocean_units'
    nz        = 20

    call init_field( field_mesh_3D_ocean, d_mesh_3D_ocean, wd_mesh_3D_ocean, &
      mesh, third_dimension%ocean_depth( nz), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_ocean,1)
    ub1_a = ubound( d_mesh_3D_ocean,1)

    lb1_f = lbound( field_mesh_3D_ocean%d,1)
    ub1_f = ubound( field_mesh_3D_ocean%d,1)

    lb2_a = lbound( d_mesh_3D_ocean,2)
    ub2_a = ubound( d_mesh_3D_ocean,2)

    lb2_f = lbound( field_mesh_3D_ocean%d,2)
    ub2_f = ubound( field_mesh_3D_ocean%d,2)

    d_mesh_3D_ocean( mesh%vi1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_ocean%name()      == name .and. &
      field_mesh_3D_ocean%long_name() == long_name .and. &
      field_mesh_3D_ocean%units()     == units .and. &
      field_mesh_3D_ocean%is_parent( mesh) .and. &
      field_mesh_3D_ocean%is_parent( third_dimension%ocean_depth( nz)) .and. &
      field_mesh_3D_ocean%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == mesh%vi1 .and. &
      ub1_a == mesh%vi2 .and. &
      lb1_f == mesh%vi1 .and. &
      ub1_f == mesh%vi2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_ocean%d( mesh%vi1+1, 3) == test_val), &
      trim( test_name) // '_3D_ocean')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_init_field_mesh_int_a

  subroutine test_init_field_mesh_int_b( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),        intent(in) :: test_name_parent
    type(type_mesh), target, intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_init_field_mesh_int_b'
    character(len=1024), parameter               :: test_name_local = 'mesh/b/int'
    character(len=1024)                          :: test_name
    character(len=1024)                          :: name, long_name, units
    integer                                      :: nz
    type(type_field_mesh_int_2D), allocatable    :: field_mesh_2D
    type(type_field_mesh_int_3D), allocatable    :: field_mesh_3D_zeta
    type(type_field_mesh_int_3D), allocatable    :: field_mesh_3D_month
    type(type_field_mesh_int_3D), allocatable    :: field_mesh_3D_ocean
    integer, dimension(:  ), contiguous, pointer :: d_mesh_2D
    integer, dimension(:,:), contiguous, pointer :: d_mesh_3D_zeta
    integer, dimension(:,:), contiguous, pointer :: d_mesh_3D_month
    integer, dimension(:,:), contiguous, pointer :: d_mesh_3D_ocean
    type(MPI_WIN)                                :: wd_mesh_2D
    type(MPI_WIN)                                :: wd_mesh_3D_zeta
    type(MPI_WIN)                                :: wd_mesh_3D_month
    type(MPI_WIN)                                :: wd_mesh_3D_ocean
    integer                                      :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                      :: lb1_f, ub1_f, lb2_f, ub2_f
    integer, parameter                           :: test_val = 1337

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! 2-D, a-mesh

    name      = 'd_mesh_2D'
    long_name = 'd_mesh_2D_long_name'
    units     = 'd_mesh_2D_units'

    call init_field( field_mesh_2D, d_mesh_2D, wd_mesh_2D, &
      mesh, Arakawa_grid%b(), name, long_name, units)

    lb1_a = lbound( d_mesh_2D,1)
    ub1_a = ubound( d_mesh_2D,1)

    lb1_f = lbound( field_mesh_2D%d,1)
    ub1_f = ubound( field_mesh_2D%d,1)

    d_mesh_2D( mesh%ti1+1) = test_val

    call unit_test( (&
      field_mesh_2D%name()      == name .and. &
      field_mesh_2D%long_name() == long_name .and. &
      field_mesh_2D%units()     == units .and. &
      field_mesh_2D%is_parent( mesh) .and. &
      field_mesh_2D%is_parent( Arakawa_grid%b()) .and. &
      lb1_a == mesh%ti1 .and. &
      ub1_a == mesh%ti2 .and. &
      lb1_f == mesh%ti1 .and. &
      ub1_f == mesh%ti2 .and. &
      field_mesh_2D%d( mesh%ti1+1) == test_val), &
      trim( test_name) // '_2D')

    ! 3-D (ice zeta), a-mesh

    name      = 'd_mesh_3D_zeta'
    long_name = 'd_mesh_3D_zeta_long_name'
    units     = 'd_mesh_3D_zeta_units'
    nz        = 10

    call init_field( field_mesh_3D_zeta, d_mesh_3D_zeta, wd_mesh_3D_zeta, &
      mesh, third_dimension%ice_zeta( nz), Arakawa_grid%b(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_zeta,1)
    ub1_a = ubound( d_mesh_3D_zeta,1)

    lb1_f = lbound( field_mesh_3D_zeta%d,1)
    ub1_f = ubound( field_mesh_3D_zeta%d,1)

    lb2_a = lbound( d_mesh_3D_zeta,2)
    ub2_a = ubound( d_mesh_3D_zeta,2)

    lb2_f = lbound( field_mesh_3D_zeta%d,2)
    ub2_f = ubound( field_mesh_3D_zeta%d,2)

    d_mesh_3D_zeta( mesh%ti1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_zeta%name()      == name .and. &
      field_mesh_3D_zeta%long_name() == long_name .and. &
      field_mesh_3D_zeta%units()     == units .and. &
      field_mesh_3D_zeta%is_parent( mesh) .and. &
      field_mesh_3D_zeta%is_parent( third_dimension%ice_zeta( nz)) .and. &
      field_mesh_3D_zeta%is_parent( Arakawa_grid%b()) .and. &
      lb1_a == mesh%ti1 .and. &
      ub1_a == mesh%ti2 .and. &
      lb1_f == mesh%ti1 .and. &
      ub1_f == mesh%ti2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_zeta%d( mesh%ti1+1, 3) == test_val), &
      trim( test_name) // '_3D_zeta')

    ! 3-D (month), a-mesh

    name      = 'd_mesh_3D_month'
    long_name = 'd_mesh_3D_month_long_name'
    units     = 'd_mesh_3D_month_units'

    call init_field( field_mesh_3D_month, d_mesh_3D_month, wd_mesh_3D_month, &
      mesh, third_dimension%month(), Arakawa_grid%b(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_month,1)
    ub1_a = ubound( d_mesh_3D_month,1)

    lb1_f = lbound( field_mesh_3D_month%d,1)
    ub1_f = ubound( field_mesh_3D_month%d,1)

    lb2_a = lbound( d_mesh_3D_month,2)
    ub2_a = ubound( d_mesh_3D_month,2)

    lb2_f = lbound( field_mesh_3D_month%d,2)
    ub2_f = ubound( field_mesh_3D_month%d,2)

    d_mesh_3D_month( mesh%ti1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_month%name()      == name .and. &
      field_mesh_3D_month%long_name() == long_name .and. &
      field_mesh_3D_month%units()     == units .and. &
      field_mesh_3D_month%is_parent( mesh) .and. &
      field_mesh_3D_month%is_parent( third_dimension%month()) .and. &
      field_mesh_3D_month%is_parent( Arakawa_grid%b()) .and. &
      lb1_a == mesh%ti1 .and. &
      ub1_a == mesh%ti2 .and. &
      lb1_f == mesh%ti1 .and. &
      ub1_f == mesh%ti2 .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      field_mesh_3D_month%d( mesh%ti1+1, 3) == test_val), &
      trim( test_name) // '_3D_month')

    ! 3-D (ocean depth), a-mesh

    name      = 'd_mesh_3D_ocean'
    long_name = 'd_mesh_3D_ocean_long_name'
    units     = 'd_mesh_3D_ocean_units'
    nz        = 20

    call init_field( field_mesh_3D_ocean, d_mesh_3D_ocean, wd_mesh_3D_ocean, &
      mesh, third_dimension%ocean_depth( nz), Arakawa_grid%b(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_ocean,1)
    ub1_a = ubound( d_mesh_3D_ocean,1)

    lb1_f = lbound( field_mesh_3D_ocean%d,1)
    ub1_f = ubound( field_mesh_3D_ocean%d,1)

    lb2_a = lbound( d_mesh_3D_ocean,2)
    ub2_a = ubound( d_mesh_3D_ocean,2)

    lb2_f = lbound( field_mesh_3D_ocean%d,2)
    ub2_f = ubound( field_mesh_3D_ocean%d,2)

    d_mesh_3D_ocean( mesh%ti1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_ocean%name()      == name .and. &
      field_mesh_3D_ocean%long_name() == long_name .and. &
      field_mesh_3D_ocean%units()     == units .and. &
      field_mesh_3D_ocean%is_parent( mesh) .and. &
      field_mesh_3D_ocean%is_parent( third_dimension%ocean_depth( nz)) .and. &
      field_mesh_3D_ocean%is_parent( Arakawa_grid%b()) .and. &
      lb1_a == mesh%ti1 .and. &
      ub1_a == mesh%ti2 .and. &
      lb1_f == mesh%ti1 .and. &
      ub1_f == mesh%ti2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_ocean%d( mesh%ti1+1, 3) == test_val), &
      trim( test_name) // '_3D_ocean')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_init_field_mesh_int_b

  subroutine test_init_field_mesh_int_c( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),        intent(in) :: test_name_parent
    type(type_mesh), target, intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_init_field_mesh_int_c'
    character(len=1024), parameter               :: test_name_local = 'mesh/c/int'
    character(len=1024)                          :: test_name
    character(len=1024)                          :: name, long_name, units
    integer                                      :: nz
    type(type_field_mesh_int_2D), allocatable    :: field_mesh_2D
    type(type_field_mesh_int_3D), allocatable    :: field_mesh_3D_zeta
    type(type_field_mesh_int_3D), allocatable    :: field_mesh_3D_month
    type(type_field_mesh_int_3D), allocatable    :: field_mesh_3D_ocean
    integer, dimension(:  ), contiguous, pointer :: d_mesh_2D
    integer, dimension(:,:), contiguous, pointer :: d_mesh_3D_zeta
    integer, dimension(:,:), contiguous, pointer :: d_mesh_3D_month
    integer, dimension(:,:), contiguous, pointer :: d_mesh_3D_ocean
    type(MPI_WIN)                                :: wd_mesh_2D
    type(MPI_WIN)                                :: wd_mesh_3D_zeta
    type(MPI_WIN)                                :: wd_mesh_3D_month
    type(MPI_WIN)                                :: wd_mesh_3D_ocean
    integer                                      :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                      :: lb1_f, ub1_f, lb2_f, ub2_f
    integer, parameter                           :: test_val = 1337

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! 2-D, a-mesh

    name      = 'd_mesh_2D'
    long_name = 'd_mesh_2D_long_name'
    units     = 'd_mesh_2D_units'

    call init_field( field_mesh_2D, d_mesh_2D, wd_mesh_2D, &
      mesh, Arakawa_grid%c(), name, long_name, units)

    lb1_a = lbound( d_mesh_2D,1)
    ub1_a = ubound( d_mesh_2D,1)

    lb1_f = lbound( field_mesh_2D%d,1)
    ub1_f = ubound( field_mesh_2D%d,1)

    d_mesh_2D( mesh%ei1+1) = test_val

    call unit_test( (&
      field_mesh_2D%name()      == name .and. &
      field_mesh_2D%long_name() == long_name .and. &
      field_mesh_2D%units()     == units .and. &
      field_mesh_2D%is_parent( mesh) .and. &
      field_mesh_2D%is_parent( Arakawa_grid%c()) .and. &
      lb1_a == mesh%ei1 .and. &
      ub1_a == mesh%ei2 .and. &
      lb1_f == mesh%ei1 .and. &
      ub1_f == mesh%ei2 .and. &
      field_mesh_2D%d( mesh%ei1+1) == test_val), &
      trim( test_name) // '_2D')

    ! 3-D (ice zeta), a-mesh

    name      = 'd_mesh_3D_zeta'
    long_name = 'd_mesh_3D_zeta_long_name'
    units     = 'd_mesh_3D_zeta_units'
    nz        = 10

    call init_field( field_mesh_3D_zeta, d_mesh_3D_zeta, wd_mesh_3D_zeta, &
      mesh, third_dimension%ice_zeta( nz), Arakawa_grid%c(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_zeta,1)
    ub1_a = ubound( d_mesh_3D_zeta,1)

    lb1_f = lbound( field_mesh_3D_zeta%d,1)
    ub1_f = ubound( field_mesh_3D_zeta%d,1)

    lb2_a = lbound( d_mesh_3D_zeta,2)
    ub2_a = ubound( d_mesh_3D_zeta,2)

    lb2_f = lbound( field_mesh_3D_zeta%d,2)
    ub2_f = ubound( field_mesh_3D_zeta%d,2)

    d_mesh_3D_zeta( mesh%ei1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_zeta%name()      == name .and. &
      field_mesh_3D_zeta%long_name() == long_name .and. &
      field_mesh_3D_zeta%units()     == units .and. &
      field_mesh_3D_zeta%is_parent( mesh) .and. &
      field_mesh_3D_zeta%is_parent( third_dimension%ice_zeta( nz)) .and. &
      field_mesh_3D_zeta%is_parent( Arakawa_grid%c()) .and. &
      lb1_a == mesh%ei1 .and. &
      ub1_a == mesh%ei2 .and. &
      lb1_f == mesh%ei1 .and. &
      ub1_f == mesh%ei2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_zeta%d( mesh%ei1+1, 3) == test_val), &
      trim( test_name) // '_3D_zeta')

    ! 3-D (month), a-mesh

    name      = 'd_mesh_3D_month'
    long_name = 'd_mesh_3D_month_long_name'
    units     = 'd_mesh_3D_month_units'

    call init_field( field_mesh_3D_month, d_mesh_3D_month, wd_mesh_3D_month, &
      mesh, third_dimension%month(), Arakawa_grid%c(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_month,1)
    ub1_a = ubound( d_mesh_3D_month,1)

    lb1_f = lbound( field_mesh_3D_month%d,1)
    ub1_f = ubound( field_mesh_3D_month%d,1)

    lb2_a = lbound( d_mesh_3D_month,2)
    ub2_a = ubound( d_mesh_3D_month,2)

    lb2_f = lbound( field_mesh_3D_month%d,2)
    ub2_f = ubound( field_mesh_3D_month%d,2)

    d_mesh_3D_month( mesh%ei1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_month%name()      == name .and. &
      field_mesh_3D_month%long_name() == long_name .and. &
      field_mesh_3D_month%units()     == units .and. &
      field_mesh_3D_month%is_parent( mesh) .and. &
      field_mesh_3D_month%is_parent( third_dimension%month()) .and. &
      field_mesh_3D_month%is_parent( Arakawa_grid%c()) .and. &
      lb1_a == mesh%ei1 .and. &
      ub1_a == mesh%ei2 .and. &
      lb1_f == mesh%ei1 .and. &
      ub1_f == mesh%ei2 .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      field_mesh_3D_month%d( mesh%ei1+1, 3) == test_val), &
      trim( test_name) // '_3D_month')

    ! 3-D (ocean depth), a-mesh

    name      = 'd_mesh_3D_ocean'
    long_name = 'd_mesh_3D_ocean_long_name'
    units     = 'd_mesh_3D_ocean_units'
    nz        = 20

    call init_field( field_mesh_3D_ocean, d_mesh_3D_ocean, wd_mesh_3D_ocean, &
      mesh, third_dimension%ocean_depth( nz), Arakawa_grid%c(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_ocean,1)
    ub1_a = ubound( d_mesh_3D_ocean,1)

    lb1_f = lbound( field_mesh_3D_ocean%d,1)
    ub1_f = ubound( field_mesh_3D_ocean%d,1)

    lb2_a = lbound( d_mesh_3D_ocean,2)
    ub2_a = ubound( d_mesh_3D_ocean,2)

    lb2_f = lbound( field_mesh_3D_ocean%d,2)
    ub2_f = ubound( field_mesh_3D_ocean%d,2)

    d_mesh_3D_ocean( mesh%ei1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_ocean%name()      == name .and. &
      field_mesh_3D_ocean%long_name() == long_name .and. &
      field_mesh_3D_ocean%units()     == units .and. &
      field_mesh_3D_ocean%is_parent( mesh) .and. &
      field_mesh_3D_ocean%is_parent( third_dimension%ocean_depth( nz)) .and. &
      field_mesh_3D_ocean%is_parent( Arakawa_grid%c()) .and. &
      lb1_a == mesh%ei1 .and. &
      ub1_a == mesh%ei2 .and. &
      lb1_f == mesh%ei1 .and. &
      ub1_f == mesh%ei2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_ocean%d( mesh%ei1+1, 3) == test_val), &
      trim( test_name) // '_3D_ocean')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_init_field_mesh_int_c

  subroutine test_init_field_mesh_dp_a( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),        intent(in) :: test_name_parent
    type(type_mesh), target, intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_init_field_mesh_dp_a'
    character(len=1024), parameter                :: test_name_local = 'mesh/a/dp'
    character(len=1024)                           :: test_name
    character(len=1024)                           :: name, long_name, units
    integer                                       :: nz
    type(type_field_mesh_dp_2D), allocatable      :: field_mesh_2D
    type(type_field_mesh_dp_3D), allocatable      :: field_mesh_3D_zeta
    type(type_field_mesh_dp_3D), allocatable      :: field_mesh_3D_month
    type(type_field_mesh_dp_3D), allocatable      :: field_mesh_3D_ocean
    real(dp), dimension(:  ), contiguous, pointer :: d_mesh_2D
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_3D_zeta
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_3D_month
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_3D_ocean
    type(MPI_WIN)                                 :: wd_mesh_2D
    type(MPI_WIN)                                 :: wd_mesh_3D_zeta
    type(MPI_WIN)                                 :: wd_mesh_3D_month
    type(MPI_WIN)                                 :: wd_mesh_3D_ocean
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f
    real(dp), parameter                           :: test_val = 13.37_dp

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! 2-D, a-mesh

    name      = 'd_mesh_2D'
    long_name = 'd_mesh_2D_long_name'
    units     = 'd_mesh_2D_units'

    call init_field( field_mesh_2D, d_mesh_2D, wd_mesh_2D, &
      mesh, Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_mesh_2D,1)
    ub1_a = ubound( d_mesh_2D,1)

    lb1_f = lbound( field_mesh_2D%d,1)
    ub1_f = ubound( field_mesh_2D%d,1)

    d_mesh_2D( mesh%vi1+1) = test_val

    call unit_test( (&
      field_mesh_2D%name()      == name .and. &
      field_mesh_2D%long_name() == long_name .and. &
      field_mesh_2D%units()     == units .and. &
      field_mesh_2D%is_parent( mesh) .and. &
      field_mesh_2D%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == mesh%vi1 .and. &
      ub1_a == mesh%vi2 .and. &
      lb1_f == mesh%vi1 .and. &
      ub1_f == mesh%vi2 .and. &
      field_mesh_2D%d( mesh%vi1+1) == test_val), &
      trim( test_name) // '_2D')

    ! 3-D (ice zeta), a-mesh

    name      = 'd_mesh_3D_zeta'
    long_name = 'd_mesh_3D_zeta_long_name'
    units     = 'd_mesh_3D_zeta_units'
    nz        = 10

    call init_field( field_mesh_3D_zeta, d_mesh_3D_zeta, wd_mesh_3D_zeta, &
      mesh, third_dimension%ice_zeta( nz), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_zeta,1)
    ub1_a = ubound( d_mesh_3D_zeta,1)

    lb1_f = lbound( field_mesh_3D_zeta%d,1)
    ub1_f = ubound( field_mesh_3D_zeta%d,1)

    lb2_a = lbound( d_mesh_3D_zeta,2)
    ub2_a = ubound( d_mesh_3D_zeta,2)

    lb2_f = lbound( field_mesh_3D_zeta%d,2)
    ub2_f = ubound( field_mesh_3D_zeta%d,2)

    d_mesh_3D_zeta( mesh%vi1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_zeta%name()      == name .and. &
      field_mesh_3D_zeta%long_name() == long_name .and. &
      field_mesh_3D_zeta%units()     == units .and. &
      field_mesh_3D_zeta%is_parent( mesh) .and. &
      field_mesh_3D_zeta%is_parent( third_dimension%ice_zeta( nz)) .and. &
      field_mesh_3D_zeta%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == mesh%vi1 .and. &
      ub1_a == mesh%vi2 .and. &
      lb1_f == mesh%vi1 .and. &
      ub1_f == mesh%vi2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_zeta%d( mesh%vi1+1, 3) == test_val), &
      trim( test_name) // '_3D_zeta')

    ! 3-D (month), a-mesh

    name      = 'd_mesh_3D_month'
    long_name = 'd_mesh_3D_month_long_name'
    units     = 'd_mesh_3D_month_units'

    call init_field( field_mesh_3D_month, d_mesh_3D_month, wd_mesh_3D_month, &
      mesh, third_dimension%month(), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_month,1)
    ub1_a = ubound( d_mesh_3D_month,1)

    lb1_f = lbound( field_mesh_3D_month%d,1)
    ub1_f = ubound( field_mesh_3D_month%d,1)

    lb2_a = lbound( d_mesh_3D_month,2)
    ub2_a = ubound( d_mesh_3D_month,2)

    lb2_f = lbound( field_mesh_3D_month%d,2)
    ub2_f = ubound( field_mesh_3D_month%d,2)

    d_mesh_3D_month( mesh%vi1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_month%name()      == name .and. &
      field_mesh_3D_month%long_name() == long_name .and. &
      field_mesh_3D_month%units()     == units .and. &
      field_mesh_3D_month%is_parent( mesh) .and. &
      field_mesh_3D_month%is_parent( third_dimension%month()) .and. &
      field_mesh_3D_month%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == mesh%vi1 .and. &
      ub1_a == mesh%vi2 .and. &
      lb1_f == mesh%vi1 .and. &
      ub1_f == mesh%vi2 .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      field_mesh_3D_month%d( mesh%vi1+1, 3) == test_val), &
      trim( test_name) // '_3D_month')

    ! 3-D (ocean depth), a-mesh

    name      = 'd_mesh_3D_ocean'
    long_name = 'd_mesh_3D_ocean_long_name'
    units     = 'd_mesh_3D_ocean_units'
    nz        = 20

    call init_field( field_mesh_3D_ocean, d_mesh_3D_ocean, wd_mesh_3D_ocean, &
      mesh, third_dimension%ocean_depth( nz), Arakawa_grid%a(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_ocean,1)
    ub1_a = ubound( d_mesh_3D_ocean,1)

    lb1_f = lbound( field_mesh_3D_ocean%d,1)
    ub1_f = ubound( field_mesh_3D_ocean%d,1)

    lb2_a = lbound( d_mesh_3D_ocean,2)
    ub2_a = ubound( d_mesh_3D_ocean,2)

    lb2_f = lbound( field_mesh_3D_ocean%d,2)
    ub2_f = ubound( field_mesh_3D_ocean%d,2)

    d_mesh_3D_ocean( mesh%vi1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_ocean%name()      == name .and. &
      field_mesh_3D_ocean%long_name() == long_name .and. &
      field_mesh_3D_ocean%units()     == units .and. &
      field_mesh_3D_ocean%is_parent( mesh) .and. &
      field_mesh_3D_ocean%is_parent( third_dimension%ocean_depth( nz)) .and. &
      field_mesh_3D_ocean%is_parent( Arakawa_grid%a()) .and. &
      lb1_a == mesh%vi1 .and. &
      ub1_a == mesh%vi2 .and. &
      lb1_f == mesh%vi1 .and. &
      ub1_f == mesh%vi2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_ocean%d( mesh%vi1+1, 3) == test_val), &
      trim( test_name) // '_3D_ocean')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_init_field_mesh_dp_a

  subroutine test_init_field_mesh_dp_b( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),        intent(in) :: test_name_parent
    type(type_mesh), target, intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_init_field_mesh_dp_b'
    character(len=1024), parameter                :: test_name_local = 'mesh/b/dp'
    character(len=1024)                           :: test_name
    character(len=1024)                           :: name, long_name, units
    integer                                       :: nz
    type(type_field_mesh_dp_2D), allocatable      :: field_mesh_2D
    type(type_field_mesh_dp_3D), allocatable      :: field_mesh_3D_zeta
    type(type_field_mesh_dp_3D), allocatable      :: field_mesh_3D_month
    type(type_field_mesh_dp_3D), allocatable      :: field_mesh_3D_ocean
    real(dp), dimension(:  ), contiguous, pointer :: d_mesh_2D
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_3D_zeta
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_3D_month
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_3D_ocean
    type(MPI_WIN)                                 :: wd_mesh_2D
    type(MPI_WIN)                                 :: wd_mesh_3D_zeta
    type(MPI_WIN)                                 :: wd_mesh_3D_month
    type(MPI_WIN)                                 :: wd_mesh_3D_ocean
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f
    real(dp), parameter                           :: test_val = 13.37_dp

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! 2-D, a-mesh

    name      = 'd_mesh_2D'
    long_name = 'd_mesh_2D_long_name'
    units     = 'd_mesh_2D_units'

    call init_field( field_mesh_2D, d_mesh_2D, wd_mesh_2D, &
      mesh, Arakawa_grid%b(), name, long_name, units)

    lb1_a = lbound( d_mesh_2D,1)
    ub1_a = ubound( d_mesh_2D,1)

    lb1_f = lbound( field_mesh_2D%d,1)
    ub1_f = ubound( field_mesh_2D%d,1)

    d_mesh_2D( mesh%ti1+1) = test_val

    call unit_test( (&
      field_mesh_2D%name()      == name .and. &
      field_mesh_2D%long_name() == long_name .and. &
      field_mesh_2D%units()     == units .and. &
      field_mesh_2D%is_parent( mesh) .and. &
      field_mesh_2D%is_parent( Arakawa_grid%b()) .and. &
      lb1_a == mesh%ti1 .and. &
      ub1_a == mesh%ti2 .and. &
      lb1_f == mesh%ti1 .and. &
      ub1_f == mesh%ti2 .and. &
      field_mesh_2D%d( mesh%ti1+1) == test_val), &
      trim( test_name) // '_2D')

    ! 3-D (ice zeta), a-mesh

    name      = 'd_mesh_3D_zeta'
    long_name = 'd_mesh_3D_zeta_long_name'
    units     = 'd_mesh_3D_zeta_units'
    nz        = 10

    call init_field( field_mesh_3D_zeta, d_mesh_3D_zeta, wd_mesh_3D_zeta, &
      mesh, third_dimension%ice_zeta( nz), Arakawa_grid%b(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_zeta,1)
    ub1_a = ubound( d_mesh_3D_zeta,1)

    lb1_f = lbound( field_mesh_3D_zeta%d,1)
    ub1_f = ubound( field_mesh_3D_zeta%d,1)

    lb2_a = lbound( d_mesh_3D_zeta,2)
    ub2_a = ubound( d_mesh_3D_zeta,2)

    lb2_f = lbound( field_mesh_3D_zeta%d,2)
    ub2_f = ubound( field_mesh_3D_zeta%d,2)

    d_mesh_3D_zeta( mesh%ti1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_zeta%name()      == name .and. &
      field_mesh_3D_zeta%long_name() == long_name .and. &
      field_mesh_3D_zeta%units()     == units .and. &
      field_mesh_3D_zeta%is_parent( mesh) .and. &
      field_mesh_3D_zeta%is_parent( third_dimension%ice_zeta( nz)) .and. &
      field_mesh_3D_zeta%is_parent( Arakawa_grid%b()) .and. &
      lb1_a == mesh%ti1 .and. &
      ub1_a == mesh%ti2 .and. &
      lb1_f == mesh%ti1 .and. &
      ub1_f == mesh%ti2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_zeta%d( mesh%ti1+1, 3) == test_val), &
      trim( test_name) // '_3D_zeta')

    ! 3-D (month), a-mesh

    name      = 'd_mesh_3D_month'
    long_name = 'd_mesh_3D_month_long_name'
    units     = 'd_mesh_3D_month_units'

    call init_field( field_mesh_3D_month, d_mesh_3D_month, wd_mesh_3D_month, &
      mesh, third_dimension%month(), Arakawa_grid%b(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_month,1)
    ub1_a = ubound( d_mesh_3D_month,1)

    lb1_f = lbound( field_mesh_3D_month%d,1)
    ub1_f = ubound( field_mesh_3D_month%d,1)

    lb2_a = lbound( d_mesh_3D_month,2)
    ub2_a = ubound( d_mesh_3D_month,2)

    lb2_f = lbound( field_mesh_3D_month%d,2)
    ub2_f = ubound( field_mesh_3D_month%d,2)

    d_mesh_3D_month( mesh%ti1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_month%name()      == name .and. &
      field_mesh_3D_month%long_name() == long_name .and. &
      field_mesh_3D_month%units()     == units .and. &
      field_mesh_3D_month%is_parent( mesh) .and. &
      field_mesh_3D_month%is_parent( third_dimension%month()) .and. &
      field_mesh_3D_month%is_parent( Arakawa_grid%b()) .and. &
      lb1_a == mesh%ti1 .and. &
      ub1_a == mesh%ti2 .and. &
      lb1_f == mesh%ti1 .and. &
      ub1_f == mesh%ti2 .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      field_mesh_3D_month%d( mesh%ti1+1, 3) == test_val), &
      trim( test_name) // '_3D_month')

    ! 3-D (ocean depth), a-mesh

    name      = 'd_mesh_3D_ocean'
    long_name = 'd_mesh_3D_ocean_long_name'
    units     = 'd_mesh_3D_ocean_units'
    nz        = 20

    call init_field( field_mesh_3D_ocean, d_mesh_3D_ocean, wd_mesh_3D_ocean, &
      mesh, third_dimension%ocean_depth( nz), Arakawa_grid%b(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_ocean,1)
    ub1_a = ubound( d_mesh_3D_ocean,1)

    lb1_f = lbound( field_mesh_3D_ocean%d,1)
    ub1_f = ubound( field_mesh_3D_ocean%d,1)

    lb2_a = lbound( d_mesh_3D_ocean,2)
    ub2_a = ubound( d_mesh_3D_ocean,2)

    lb2_f = lbound( field_mesh_3D_ocean%d,2)
    ub2_f = ubound( field_mesh_3D_ocean%d,2)

    d_mesh_3D_ocean( mesh%ti1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_ocean%name()      == name .and. &
      field_mesh_3D_ocean%long_name() == long_name .and. &
      field_mesh_3D_ocean%units()     == units .and. &
      field_mesh_3D_ocean%is_parent( mesh) .and. &
      field_mesh_3D_ocean%is_parent( third_dimension%ocean_depth( nz)) .and. &
      field_mesh_3D_ocean%is_parent( Arakawa_grid%b()) .and. &
      lb1_a == mesh%ti1 .and. &
      ub1_a == mesh%ti2 .and. &
      lb1_f == mesh%ti1 .and. &
      ub1_f == mesh%ti2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_ocean%d( mesh%ti1+1, 3) == test_val), &
      trim( test_name) // '_3D_ocean')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_init_field_mesh_dp_b

  subroutine test_init_field_mesh_dp_c( test_name_parent, mesh)

    ! In/output variables:
    character(len=*),        intent(in) :: test_name_parent
    type(type_mesh), target, intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_init_field_mesh_dp_c'
    character(len=1024), parameter                :: test_name_local = 'mesh/c/dp'
    character(len=1024)                           :: test_name
    character(len=1024)                           :: name, long_name, units
    integer                                       :: nz
    type(type_field_mesh_dp_2D), allocatable      :: field_mesh_2D
    type(type_field_mesh_dp_3D), allocatable      :: field_mesh_3D_zeta
    type(type_field_mesh_dp_3D), allocatable      :: field_mesh_3D_month
    type(type_field_mesh_dp_3D), allocatable      :: field_mesh_3D_ocean
    real(dp), dimension(:  ), contiguous, pointer :: d_mesh_2D
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_3D_zeta
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_3D_month
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_3D_ocean
    type(MPI_WIN)                                 :: wd_mesh_2D
    type(MPI_WIN)                                 :: wd_mesh_3D_zeta
    type(MPI_WIN)                                 :: wd_mesh_3D_month
    type(MPI_WIN)                                 :: wd_mesh_3D_ocean
    integer                                       :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                       :: lb1_f, ub1_f, lb2_f, ub2_f
    real(dp), parameter                           :: test_val = 13.37_dp

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! 2-D, a-mesh

    name      = 'd_mesh_2D'
    long_name = 'd_mesh_2D_long_name'
    units     = 'd_mesh_2D_units'

    call init_field( field_mesh_2D, d_mesh_2D, wd_mesh_2D, &
      mesh, Arakawa_grid%c(), name, long_name, units)

    lb1_a = lbound( d_mesh_2D,1)
    ub1_a = ubound( d_mesh_2D,1)

    lb1_f = lbound( field_mesh_2D%d,1)
    ub1_f = ubound( field_mesh_2D%d,1)

    d_mesh_2D( mesh%ei1+1) = test_val

    call unit_test( (&
      field_mesh_2D%name()      == name .and. &
      field_mesh_2D%long_name() == long_name .and. &
      field_mesh_2D%units()     == units .and. &
      field_mesh_2D%is_parent( mesh) .and. &
      field_mesh_2D%is_parent( Arakawa_grid%c()) .and. &
      lb1_a == mesh%ei1 .and. &
      ub1_a == mesh%ei2 .and. &
      lb1_f == mesh%ei1 .and. &
      ub1_f == mesh%ei2 .and. &
      field_mesh_2D%d( mesh%ei1+1) == test_val), &
      trim( test_name) // '_2D')

    ! 3-D (ice zeta), a-mesh

    name      = 'd_mesh_3D_zeta'
    long_name = 'd_mesh_3D_zeta_long_name'
    units     = 'd_mesh_3D_zeta_units'
    nz        = 10

    call init_field( field_mesh_3D_zeta, d_mesh_3D_zeta, wd_mesh_3D_zeta, &
      mesh, third_dimension%ice_zeta( nz), Arakawa_grid%c(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_zeta,1)
    ub1_a = ubound( d_mesh_3D_zeta,1)

    lb1_f = lbound( field_mesh_3D_zeta%d,1)
    ub1_f = ubound( field_mesh_3D_zeta%d,1)

    lb2_a = lbound( d_mesh_3D_zeta,2)
    ub2_a = ubound( d_mesh_3D_zeta,2)

    lb2_f = lbound( field_mesh_3D_zeta%d,2)
    ub2_f = ubound( field_mesh_3D_zeta%d,2)

    d_mesh_3D_zeta( mesh%ei1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_zeta%name()      == name .and. &
      field_mesh_3D_zeta%long_name() == long_name .and. &
      field_mesh_3D_zeta%units()     == units .and. &
      field_mesh_3D_zeta%is_parent( mesh) .and. &
      field_mesh_3D_zeta%is_parent( third_dimension%ice_zeta( nz)) .and. &
      field_mesh_3D_zeta%is_parent( Arakawa_grid%c()) .and. &
      lb1_a == mesh%ei1 .and. &
      ub1_a == mesh%ei2 .and. &
      lb1_f == mesh%ei1 .and. &
      ub1_f == mesh%ei2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_zeta%d( mesh%ei1+1, 3) == test_val), &
      trim( test_name) // '_3D_zeta')

    ! 3-D (month), a-mesh

    name      = 'd_mesh_3D_month'
    long_name = 'd_mesh_3D_month_long_name'
    units     = 'd_mesh_3D_month_units'

    call init_field( field_mesh_3D_month, d_mesh_3D_month, wd_mesh_3D_month, &
      mesh, third_dimension%month(), Arakawa_grid%c(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_month,1)
    ub1_a = ubound( d_mesh_3D_month,1)

    lb1_f = lbound( field_mesh_3D_month%d,1)
    ub1_f = ubound( field_mesh_3D_month%d,1)

    lb2_a = lbound( d_mesh_3D_month,2)
    ub2_a = ubound( d_mesh_3D_month,2)

    lb2_f = lbound( field_mesh_3D_month%d,2)
    ub2_f = ubound( field_mesh_3D_month%d,2)

    d_mesh_3D_month( mesh%ei1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_month%name()      == name .and. &
      field_mesh_3D_month%long_name() == long_name .and. &
      field_mesh_3D_month%units()     == units .and. &
      field_mesh_3D_month%is_parent( mesh) .and. &
      field_mesh_3D_month%is_parent( third_dimension%month()) .and. &
      field_mesh_3D_month%is_parent( Arakawa_grid%c()) .and. &
      lb1_a == mesh%ei1 .and. &
      ub1_a == mesh%ei2 .and. &
      lb1_f == mesh%ei1 .and. &
      ub1_f == mesh%ei2 .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      field_mesh_3D_month%d( mesh%ei1+1, 3) == test_val), &
      trim( test_name) // '_3D_month')

    ! 3-D (ocean depth), a-mesh

    name      = 'd_mesh_3D_ocean'
    long_name = 'd_mesh_3D_ocean_long_name'
    units     = 'd_mesh_3D_ocean_units'
    nz        = 20

    call init_field( field_mesh_3D_ocean, d_mesh_3D_ocean, wd_mesh_3D_ocean, &
      mesh, third_dimension%ocean_depth( nz), Arakawa_grid%c(), name, long_name, units)

    lb1_a = lbound( d_mesh_3D_ocean,1)
    ub1_a = ubound( d_mesh_3D_ocean,1)

    lb1_f = lbound( field_mesh_3D_ocean%d,1)
    ub1_f = ubound( field_mesh_3D_ocean%d,1)

    lb2_a = lbound( d_mesh_3D_ocean,2)
    ub2_a = ubound( d_mesh_3D_ocean,2)

    lb2_f = lbound( field_mesh_3D_ocean%d,2)
    ub2_f = ubound( field_mesh_3D_ocean%d,2)

    d_mesh_3D_ocean( mesh%ei1+1, 3) = test_val

    call unit_test( (&
      field_mesh_3D_ocean%name()      == name .and. &
      field_mesh_3D_ocean%long_name() == long_name .and. &
      field_mesh_3D_ocean%units()     == units .and. &
      field_mesh_3D_ocean%is_parent( mesh) .and. &
      field_mesh_3D_ocean%is_parent( third_dimension%ocean_depth( nz)) .and. &
      field_mesh_3D_ocean%is_parent( Arakawa_grid%c()) .and. &
      lb1_a == mesh%ei1 .and. &
      ub1_a == mesh%ei2 .and. &
      lb1_f == mesh%ei1 .and. &
      ub1_f == mesh%ei2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      field_mesh_3D_ocean%d( mesh%ei1+1, 3) == test_val), &
      trim( test_name) // '_3D_ocean')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_init_field_mesh_dp_c

end module ut_fields_init_field