module ut_models

  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  ! use ut_models_io, only: test_models_io

  implicit none

  private

  public :: unit_tests_models_main

contains

  subroutine unit_tests_models_main( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_models_main'
    character(len=1024), parameter :: test_name_local = 'models'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! call test_models_io( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_models_main

end module ut_models