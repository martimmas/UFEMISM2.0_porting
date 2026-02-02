module ut_string

  use UPSY_main, only: UPSY
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use precisions, only: dp
  use ut_basic, only: unit_test
  use mpi_basic, only: par
  use parameters, only: pi, NaN

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
    call test_colour_string                           ( test_name)
    call test_insert_val_into_string_int              ( test_name)
    call test_insert_val_into_string_dp               ( test_name)

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
    call UPSY%stru%separate_strings_by_double_vertical_bars( str, strs)
    call unit_test( &
      size(strs,1) == 5 .and. &
      strs(1)=='john' .and. &
      strs(2)=='william' .and. &
      strs(3)=='dede' .and. &
      strs(4)=='a' .and. &
      strs(5)=='ladida,di+da&  di|da', test_name // '/multiple_strings')
    deallocate( str, strs)

    str = 'pete'
    call UPSY%stru%separate_strings_by_double_vertical_bars( str, strs)
    call unit_test( &
      size(strs,1) == 1 .and. &
      strs(1)=='pete', test_name // '/single_string')
    deallocate( str, strs)

    str = ''
    call UPSY%stru%separate_strings_by_double_vertical_bars( str, strs)
    call unit_test( size(strs,1) == 0, test_name // '/empty_string')
    deallocate( str, strs)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_separate_strings_by_double_vertical_bars

  subroutine test_colour_string( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'test_colour_string'
    character(len=*), parameter   :: test_name_local = 'colour_string'
    character(len=:), allocatable :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be gray', 'gray')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be red', 'red')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be green', 'green')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be yellow', 'yellow')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be blue', 'blue')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be pink', 'pink')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be light blue', 'light blue')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  This string should be white', 'anything else')

    call unit_test( .true., test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_colour_string

  subroutine test_insert_val_into_string_int( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'test_insert_val_into_string_int'
    character(len=*), parameter   :: test_name_local = 'insert_val_into_string_int'
    character(len=:), allocatable :: test_name
    character(len=:), allocatable :: str

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    str = 'Test: {int_01}' // &
      ', {int_02}' // &
      ', john' // &
      ', {int_04}' // &
      ', {int_05}' // &
      ', a_very_long_and_nonsensical_marker_should_work_too&&*@#(#(#))' // &
      ', {int_07}' // &
      ', {int_08}' // &
      ', {int_09}'

    str = UPSY%stru%insert_val_into_string_int( str, '{int_01}', 1)
    str = UPSY%stru%insert_val_into_string_int( str, '{int_02}', 2)
    str = UPSY%stru%insert_val_into_string_int( str, 'john', 3)
    str = UPSY%stru%insert_val_into_string_int( str, '{int_04}', 4)
    str = UPSY%stru%insert_val_into_string_int( str, '{int_05}', 5)
    str = UPSY%stru%insert_val_into_string_int( str, 'a_very_long_and_nonsensical_marker_should_work_too&&*@#(#(#))', 6)
    str = UPSY%stru%insert_val_into_string_int( str, '{int_07}', 7)
    str = UPSY%stru%insert_val_into_string_int( str, '{int_08}', 8)
    str = UPSY%stru%insert_val_into_string_int( str, '{int_09}', 9)

    call unit_test( str == 'Test: 1, 2, 3, 4, 5, 6, 7, 8, 9', test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_insert_val_into_string_int

  subroutine test_insert_val_into_string_dp( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'test_insert_val_into_string_dp'
    character(len=*), parameter   :: test_name_local = 'insert_val_into_string_dp'
    character(len=:), allocatable :: test_name
    character(len=:), allocatable :: str

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    str = 'Test: {int_01}' // &
      ', {dp_02}' // &
      ', john' // &
      ', a_very_long_and_nonsensical_marker_should_work_too&&*@#(#(#))'

    str = UPSY%stru%insert_val_into_string_dp( str, '{int_01}', 1._dp)
    str = UPSY%stru%insert_val_into_string_dp( str, '{dp_02}', pi)
    str = UPSY%stru%insert_val_into_string_dp( str, 'john', 3e12_dp)
    str = UPSY%stru%insert_val_into_string_dp( str, 'a_very_long_and_nonsensical_marker_should_work_too&&*@#(#(#))', NaN)

    call unit_test( str == 'Test: 0.10000000000000E+01, 0.31415926535898E+01, 0.30000000000000E+13,                  NaN', test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_insert_val_into_string_dp

end module ut_string