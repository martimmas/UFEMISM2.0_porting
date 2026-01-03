module ut_fields_create_field

  use precisions, only: dp
  use mpi_basic, only: par
  use parameters, only: pi
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use fields_main, only: &
    atype_field, atype_field_2D, atype_field_3D, &
    type_field_logical_2D, type_field_int_2D, type_field_dp_2D, &
    type_field_logical_3D, type_field_int_3D, type_field_dp_3D, &
    type_fields_registry
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

  public :: test_create_field

contains

  subroutine test_create_field( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_create_field'
    character(len=1024), parameter :: test_name_local = 'create_field'
    character(len=1024)            :: test_name
    type(type_grid)                :: grid
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

    call setup_square_grid( 'dummy_grid', 0._dp, 1._dp, 0._dp, 1._dp, 0.1_dp, grid)

    call allocate_mesh_primary( mesh, 'dummy_mesh', 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 23.2_dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call crop_mesh_primary( mesh)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)

    call test_create_field_grid_logical  ( test_name, grid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_create_field

  subroutine test_create_field_grid_logical( test_name_parent, grid)

    ! In/output variables:
    ! In/output variables:
    character(len=*),           intent(in   ) :: test_name_parent
    type(type_grid), target,    intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_create_field_grid_logical'
    character(len=1024), parameter               :: test_name_local = 'grid/logical'
    character(len=1024)                          :: test_name
    type(type_fields_registry)                   :: flds_reg
    character(len=1024)                          :: name, long_name, units
    integer                                      :: nz
    logical, dimension(:  ), contiguous, pointer :: d_2D
    logical, dimension(:,:), contiguous, pointer :: d_3D_zeta
    logical, dimension(:,:), contiguous, pointer :: d_3D_month
    logical, dimension(:,:), contiguous, pointer :: d_3D_ocean
    type(MPI_WIN)                                :: wd_2D
    type(MPI_WIN)                                :: wd_3D_zeta
    type(MPI_WIN)                                :: wd_3D_month
    type(MPI_WIN)                                :: wd_3D_ocean
    integer                                      :: i
    integer                                      :: lb1_a, ub1_a, lb2_a, ub2_a
    integer                                      :: lb1_f, ub1_f, lb2_f, ub2_f
    logical, parameter                           :: test_val = .true.

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! 2-D
    ! ===

    name      = 'd_grid_logical_2D'
    long_name = 'd_grid_logical_2D_long_name'
    units     = 'd_grid_logical_2D_units'

    call flds_reg%create_field( d_2D, wd_2D, &
      grid, Arakawa_grid%a(), &
      name      = name, &
      long_name = long_name, &
      units     = units)

    i = flds_reg%find( name)

    lb1_a = lbound( d_2D,1)
    ub1_a = ubound( d_2D,1)

    lb1_f = flds_reg%items(i)%p%lbound( 1)
    ub1_f = flds_reg%items(i)%p%ubound( 1)

    select type (f => flds_reg%items(i)%p)
    class default
      call crash('unexpected field type')
    class is (type_field_logical_2D)
      f%d( grid%n1+1) = test_val
    end select

    call unit_test( (&
      flds_reg%items(i)%p%name()      == name .and. &
      flds_reg%items(i)%p%long_name() == long_name .and. &
      flds_reg%items(i)%p%units()     == units .and. &
      flds_reg%items(i)%p%is_grid( grid) .and. &
      flds_reg%items(i)%p%is_Arakawa_grid( Arakawa_grid%a()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      d_2D( grid%n1+1) .eqv. test_val), &
      trim( test_name) // '_2D')

    ! 3-D (zeta)
    ! ==========

    name      = 'd_grid_logical_3D_zeta'
    long_name = 'd_grid_logical_3D_zeta_long_name'
    units     = 'd_grid_logical_3D_zeta_units'
    nz        = 10

    call flds_reg%create_field( d_3D_zeta, wd_3D_zeta, &
      grid, Arakawa_grid%a(), third_dimension%ice_zeta( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units)

    i = flds_reg%find( name)

    lb1_a = lbound( d_3D_zeta,1)
    ub1_a = ubound( d_3D_zeta,1)

    lb2_a = lbound( d_3D_zeta,2)
    ub2_a = ubound( d_3D_zeta,2)

    lb1_f = flds_reg%items(i)%p%lbound( 1)
    ub1_f = flds_reg%items(i)%p%ubound( 1)

    lb2_f = flds_reg%items(i)%p%lbound( 2)
    ub2_f = flds_reg%items(i)%p%ubound( 2)

    select type (f => flds_reg%items(i)%p)
    class default
      call crash('unexpected field type')
    class is (type_field_logical_3D)
      f%d( grid%n1+1,3) = test_val
    end select

    call unit_test( (&
      flds_reg%items(i)%p%name()      == name .and. &
      flds_reg%items(i)%p%long_name() == long_name .and. &
      flds_reg%items(i)%p%units()     == units .and. &
      flds_reg%items(i)%p%is_grid( grid) .and. &
      flds_reg%items(i)%p%is_Arakawa_grid( Arakawa_grid%a()) .and. &
      flds_reg%items(i)%p%is_third_dimension( third_dimension%ice_zeta( nz)) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      d_3D_zeta( grid%n1+1,3) .eqv. test_val), &
      trim( test_name) // '_3D_zeta')

    ! 3-D (month)
    ! ===========

    name      = 'd_grid_logical_3D_month'
    long_name = 'd_grid_logical_3D_month_long_name'
    units     = 'd_grid_logical_3D_month_units'

    call flds_reg%create_field( d_3D_month, wd_3D_month, &
      grid, Arakawa_grid%a(), third_dimension%month(), &
      name      = name, &
      long_name = long_name, &
      units     = units)

    i = flds_reg%find( name)

    lb1_a = lbound( d_3D_month,1)
    ub1_a = ubound( d_3D_month,1)

    lb2_a = lbound( d_3D_month,2)
    ub2_a = ubound( d_3D_month,2)

    lb1_f = flds_reg%items(i)%p%lbound( 1)
    ub1_f = flds_reg%items(i)%p%ubound( 1)

    lb2_f = flds_reg%items(i)%p%lbound( 2)
    ub2_f = flds_reg%items(i)%p%ubound( 2)

    select type (f => flds_reg%items(i)%p)
    class default
      call crash('unexpected field type')
    class is (type_field_logical_3D)
      f%d( grid%n1+1,3) = test_val
    end select

    call unit_test( (&
      flds_reg%items(i)%p%name()      == name .and. &
      flds_reg%items(i)%p%long_name() == long_name .and. &
      flds_reg%items(i)%p%units()     == units .and. &
      flds_reg%items(i)%p%is_grid( grid) .and. &
      flds_reg%items(i)%p%is_Arakawa_grid( Arakawa_grid%a()) .and. &
      flds_reg%items(i)%p%is_third_dimension( third_dimension%month()) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      lb2_a == 1  .and. &
      ub2_a == 12 .and. &
      lb2_f == 1  .and. &
      ub2_f == 12 .and. &
      d_3D_month( grid%n1+1,3) .eqv. test_val), &
      trim( test_name) // '_3D_month')

    ! 3-D (ocean)
    ! ==========

    name      = 'd_grid_logical_3D_ocean'
    long_name = 'd_grid_logical_3D_ocean_long_name'
    units     = 'd_grid_logical_3D_ocean_units'
    nz        = 20

    call flds_reg%create_field( d_3D_ocean, wd_3D_ocean, &
      grid, Arakawa_grid%a(), third_dimension%ocean_depth( nz), &
      name      = name, &
      long_name = long_name, &
      units     = units)

    i = flds_reg%find( name)

    lb1_a = lbound( d_3D_ocean,1)
    ub1_a = ubound( d_3D_ocean,1)

    lb2_a = lbound( d_3D_ocean,2)
    ub2_a = ubound( d_3D_ocean,2)

    lb1_f = flds_reg%items(i)%p%lbound( 1)
    ub1_f = flds_reg%items(i)%p%ubound( 1)

    lb2_f = flds_reg%items(i)%p%lbound( 2)
    ub2_f = flds_reg%items(i)%p%ubound( 2)

    select type (f => flds_reg%items(i)%p)
    class default
      call crash('unexpected field type')
    class is (type_field_logical_3D)
      f%d( grid%n1+1,3) = test_val
    end select

    call unit_test( (&
      flds_reg%items(i)%p%name()      == name .and. &
      flds_reg%items(i)%p%long_name() == long_name .and. &
      flds_reg%items(i)%p%units()     == units .and. &
      flds_reg%items(i)%p%is_grid( grid) .and. &
      flds_reg%items(i)%p%is_Arakawa_grid( Arakawa_grid%a()) .and. &
      flds_reg%items(i)%p%is_third_dimension( third_dimension%ocean_depth( nz)) .and. &
      lb1_a == grid%n1 .and. &
      ub1_a == grid%n2 .and. &
      lb1_f == grid%n1 .and. &
      ub1_f == grid%n2 .and. &
      lb2_a == 1  .and. &
      ub2_a == nz .and. &
      lb2_f == 1  .and. &
      ub2_f == nz .and. &
      d_3D_ocean( grid%n1+1,3) .eqv. test_val), &
      trim( test_name) // '_3D_ocean')

    ! DENK DROM
    call flds_reg%print_info

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_create_field_grid_logical

end module ut_fields_create_field