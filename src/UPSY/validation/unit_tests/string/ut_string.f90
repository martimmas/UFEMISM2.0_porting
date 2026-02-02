module ut_string

  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use precisions, only: dp
  use ut_basic, only: unit_test
  use string_module, only: separate_strings_by_double_vertical_bars

  implicit none

  private

  public :: unit_tests_string_main

contains

  subroutine unit_tests_string_main( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'unit_tests_string_main'
    character(len=*), parameter   :: test_name_local = 'string'
    character(len=:), allocatable :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_separate_strings_by_double_vertical_bars( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_string_main

  subroutine test_separate_strings_by_double_vertical_bars( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter                    :: routine_name = 'test_separate_strings_by_double_vertical_bars'
    character(len=*), parameter                    :: test_name_local = 'separate_strings_by_double_vertical_bars'
    character(len=:), allocatable                  :: test_name
    character(len=:), allocatable                  :: str
    character(len=1024), dimension(:), allocatable :: strs

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Some test cases
    str = 'john||william||dede||a||ladida,di+da&  di|da'
    call separate_strings_by_double_vertical_bars( str, strs)
    call unit_test( &
      size(strs,1) == 5 .and. &
      strs(1)=='john' .and. &
      strs(2)=='william' .and. &
      strs(3)=='dede' .and. &
      strs(4)=='a' .and. &
      strs(5)=='ladida,di+da&  di|da', test_name // '/multiple_strings')
    deallocate( str, strs)

    str = 'pete'
    call separate_strings_by_double_vertical_bars( str, strs)
    call unit_test( &
      size(strs,1) == 1 .and. &
      strs(1)=='pete', test_name // '/single_string')
    deallocate( str, strs)

    str = ''
    call separate_strings_by_double_vertical_bars( str, strs)
    call unit_test( size(strs,1) == 0, test_name // '/empty_string')
    deallocate( str, strs)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_separate_strings_by_double_vertical_bars

end module ut_string