module ut_demo_model

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use demo_model_a, only: type_demo_model_a
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

    call unit_test( .true., test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_demo_model

end module ut_demo_model