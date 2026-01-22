module ut_models_io

  use precisions, only: dp
  use parameters, only: pi
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary, crop_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_secondary, only: calc_all_secondary_mesh_data
  ! use models_demo, only: type_demo_model
  use ut_basic, only: unit_test, foldername_unit_tests_output

  implicit none

  private

  ! public :: test_models_io

contains

  ! subroutine test_models_io( test_name_parent)

  !   ! In/output variables:
  !   character(len=*), intent(in) :: test_name_parent

  !   ! Local variables:
  !   character(len=1024), parameter :: routine_name = 'test_models_io'
  !   character(len=1024), parameter :: test_name_local = 'io'
  !   character(len=1024)            :: test_name
  !   real(dp)                       :: alpha_min, res_max
  !   real(dp), parameter            :: xmin = -1._dp
  !   real(dp), parameter            :: xmax =  1._dp
  !   real(dp), parameter            :: ymin = -1._dp
  !   real(dp), parameter            :: ymax =  1._dp
  !   type(type_mesh)                :: mesh

  !   ! Add routine to call stack
  !   call init_routine( routine_name)

  !   ! Add test name to list
  !   test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  !   call allocate_mesh_primary( mesh, 'dummy_mesh', 100, 200)
  !   call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

  !   ! Refine the test mesh
  !   alpha_min = 25._dp * pi / 180._dp
  !   res_max = pi / 23.2_dp
  !   call refine_mesh_uniform( mesh, res_max, alpha_min)
  !   call crop_mesh_primary( mesh)
  !   call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)

  !   ! Run all unit tests
  !   call test_write_to_read_from_restart_file( test_name, mesh)

  !   ! Remove routine from call stack
  !   call finalise_routine( routine_name)

  ! end subroutine test_models_io

  ! subroutine test_write_to_read_from_restart_file( test_name_parent, mesh)

  !   ! In/output variables:
  !   character(len=*), intent(in) :: test_name_parent
  !   type(type_mesh),  intent(in) :: mesh

  !   ! Local variables:
  !   character(len=1024), parameter :: routine_name = 'test_write_to_read_from_restart_file'
  !   character(len=1024), parameter :: test_name_local = 'write_to_read_from_restart_file'
  !   character(len=1024)            :: test_name
  !   type(type_demo_model)          :: demo_model1, demo_model2
  !   integer, parameter             :: nz = 15
  !   character(:), allocatable      :: filename

  !   ! Add routine to call stack
  !   call init_routine( routine_name)

  !   ! Add test name to list
  !   test_name = trim( test_name_parent) // '/' // trim( test_name_local)

  !   call demo_model1%init( mesh, nz)
  !   call demo_model1%write_to_restart_file( foldername_unit_tests_output, filename)

  !   call demo_model2%init( mesh, nz)
  !   call demo_model2%read_from_restart_file( filename)

  !   call unit_test( demo_model1 == demo_model2, test_name)

  !   ! Remove routine from call stack
  !   call finalise_routine( routine_name)

  ! end subroutine test_write_to_read_from_restart_file

end module ut_models_io