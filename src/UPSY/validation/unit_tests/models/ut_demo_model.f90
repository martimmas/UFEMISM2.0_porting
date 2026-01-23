module ut_demo_model

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning
  use mesh_types, only: type_mesh
  use demo_model_a, only: type_demo_model_a
  use demo_model_b, only: type_demo_model_b
  use ut_basic, only: unit_test, foldername_unit_tests_output

  implicit none

  private

  public :: test_demo_model

contains

  subroutine test_demo_model( test_name_parent, mesh1, mesh2)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh1, mesh2

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_demo_model'
    character(len=1024), parameter :: test_name_local = 'demo_model'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_demo_model_a( test_name, mesh1, mesh2)
    call test_demo_model_b( test_name, mesh1, mesh2)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_demo_model

  subroutine test_demo_model_a( test_name_parent, mesh1, mesh2)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh1, mesh2

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_demo_model_a'
    character(len=1024), parameter :: test_name_local = 'a'
    character(len=1024)            :: test_name
    integer, parameter             :: nz = 10
    real(dp), parameter            :: H0 = 0.1_dp
    real(dp), parameter            :: till_friction_angle_uniform = 0.1_dp
    real(dp), parameter            :: beta_sq_uniform = 1e4_dp
    real(dp), parameter            :: H_new = 0.2_dp
    real(dp), parameter            :: dH = 1._dp
    type(type_demo_model_a)        :: a1, a2
    character(:), allocatable      :: filename_a1

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate the demo model and test if that worked
    call a1%allocate( a1%ct_allocate( 'demo_model_a1', 'aaa', mesh1, nz))
    call unit_test( ( &
      a1%name()        == 'demo_model_a1' .and. &
      a1%region_name() == 'aaa' .and. &
      a1%mesh%name     == mesh1%name .and. &
      size( a1%H   ,1) == mesh1%pai_V%n_nih .and. &
      size( a1%u_3D,1) == mesh1%pai_Tri%n_nih .and. &
      size( a1%u_3D,2) == nz &
      ), trim( test_name) // '/allocate')

    ! Initialise the demo model and test if that worked
    call a1%initialise( a1%ct_initialise( H0, till_friction_angle_uniform, beta_sq_uniform))
    call unit_test( (&
      minval( a1%H) == H0 .and. &
      minval( a1%till_friction_angle) == till_friction_angle_uniform .and. &
      maxval( a1%till_friction_angle) == till_friction_angle_uniform &
      ), trim( test_name) // '/initialise')

    ! Run the demo model and test if that worked
    call a1%run( a1%ct_run( H_new, dH))
    call unit_test( (&
      minval( a1%H) == H0 + dH &
      ), trim( test_name) // '/run')

    ! Remap the demo model and test if that worked
    call a1%remap( a1%ct_remap( mesh2))
    call unit_test( ( &
      a1%mesh%name     == mesh2%name .and. &
      size( a1%H   ,1) == mesh2%pai_V%n_nih .and. &
      size( a1%u_3D,1) == mesh2%pai_Tri%n_nih .and. &
      size( a1%u_3D,2) == nz &
      ), trim( test_name) // '/remap')

    ! Write the demo model to a restart file,
    ! then allocate a new demo model and initialise it from that
    ! restart file; test if all of that worked
    call a1%write_to_restart_file( foldername_unit_tests_output, filename_a1)
    call a2%allocate( a2%ct_allocate( 'demo_model_a2', 'aaa', mesh2, nz))
    call a2%read_from_restart_file( filename_a1)
    call unit_test( a1 == a2, trim( test_name) // '/restart')

    ! Deallocate the demo model and test if that worked
    call a1%deallocate
    call unit_test( ( &
      a1%name() == 'empty_model' .and. &
      a1%region_name() == '!!!' .and. &
      .not. associated( a1%mesh) .and. &
      .not. associated( a1%H) .and. &
      .not. associated( a1%u_3D) .and. &
      .not. associated( a1%v_3D) .and. &
      .not. associated( a1%mask_ice) .and. &
      .not. associated( a1%T2m) .and. &
      .not. associated( a1%till_friction_angle) &
      ), trim( test_name) // '/deallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_demo_model_a

  subroutine test_demo_model_b( test_name_parent, mesh1, mesh2)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh1, mesh2

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_demo_model_b'
    character(len=1024), parameter :: test_name_local = 'b'
    character(len=1024)            :: test_name
    integer, parameter             :: nz = 10
    real(dp), parameter            :: H0 = 0.1_dp
    real(dp), parameter            :: till_friction_angle_uniform = 0.1_dp
    real(dp), parameter            :: beta_sq_uniform = 1e4_dp
    real(dp), parameter            :: H_new = 0.2_dp
    real(dp), parameter            :: dH = 1._dp
    type(type_demo_model_b)        :: b1, b2
    character(:), allocatable      :: filename_b1

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate the demo model and test if that worked
    call b1%allocate( b1%ct_allocate( 'demo_model_b1', 'aaa', mesh1, nz))
    call unit_test( ( &
      b1%name()        == 'demo_model_b1' .and. &
      b1%region_name() == 'aaa' .and. &
      b1%mesh%name     == mesh1%name .and. &
      size( b1%H   ,1) == mesh1%pai_V%n_nih .and. &
      size( b1%u_3D,1) == mesh1%pai_Tri%n_nih .and. &
      size( b1%u_3D,2) == nz &
      ), trim( test_name) // '/allocate')

    ! Initialise the demo model and test if that worked
    call b1%initialise( b1%ct_initialise( H0, till_friction_angle_uniform, beta_sq_uniform))
    call unit_test( (&
      minval( b1%H) == H0 .and. &
      minval( b1%beta_sq) == beta_sq_uniform .and. &
      maxval( b1%beta_sq) == beta_sq_uniform &
      ), trim( test_name) // '/initialise')

    ! Run the demo model and test if that worked
    call b1%run( b1%ct_run( H_new, dH))
    call unit_test( (&
      minval( b1%H) == H_new .and. &
      maxval( b1%H) == H_new &
      ), trim( test_name) // '/run')

    ! Remap the demo model and test if that worked
    call b1%remap( b1%ct_remap( mesh2))
    call unit_test( ( &
      b1%mesh%name     == mesh2%name .and. &
      size( b1%H   ,1) == mesh2%pai_V%n_nih .and. &
      size( b1%u_3D,1) == mesh2%pai_Tri%n_nih .and. &
      size( b1%u_3D,2) == nz &
      ), trim( test_name) // '/remap')

    ! Write the demo model to a restart file,
    ! then allocate a new demo model and initialise it from that
    ! restart file; test if all of that worked
    call b1%write_to_restart_file( foldername_unit_tests_output, filename_b1)
    call b2%allocate( b2%ct_allocate( 'demo_model_b2', 'aaa', mesh2, nz))
    call b2%read_from_restart_file( filename_b1)
    call unit_test( b1 == b2, trim( test_name) // '/restart')

    ! Deallocate the demo model and test if that worked
    call b1%deallocate
    call unit_test( ( &
      b1%name() == 'empty_model' .and. &
      b1%region_name() == '!!!' .and. &
      .not. associated( b1%mesh) .and. &
      .not. associated( b1%H) .and. &
      .not. associated( b1%u_3D) .and. &
      .not. associated( b1%v_3D) .and. &
      .not. associated( b1%mask_ice) .and. &
      .not. associated( b1%T2m) .and. &
      .not. associated( b1%beta_sq) &
      ), trim( test_name) // '/deallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_demo_model_b

end module ut_demo_model