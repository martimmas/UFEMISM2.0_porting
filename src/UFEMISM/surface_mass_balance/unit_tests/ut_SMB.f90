module ut_SMB

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use ut_basic, only: unit_test
  use tests_main, only: test_tol
  use model_configuration, only: C
  use SMB_model, only: atype_SMB_model, create_SMB_model
  use parameters, only: pi
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary, crop_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_disc_calc_matrix_operators_2D, only: calc_all_matrix_operators_mesh
  use ice_model_types, only: type_ice_model
  use climate_model_types, only: type_climate_model
  use grid_types, only: type_grid
  use Halfar_SIA_solution, only: Halfar
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, MPI_WIN
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use deallocate_dist_shared_mod, only: deallocate_dist_shared
  use netcdf_io_main, only: create_new_netcdf_file_for_writing, &
    setup_mesh_in_netcdf_file, add_month_dimension_to_file, &
    add_field_mesh_dp_2D_notime, write_to_field_multopt_mesh_dp_2D_notime, &
    close_netcdf_file

  implicit none

  private

  public :: unit_tests_SMB_main

contains

  subroutine unit_tests_SMB_main( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_SMB_main'
    character(len=1024), parameter :: test_name_local = 'SMB'
    character(len=1024)            :: test_name
    real(dp)                       :: alpha_min, res_max
    real(dp), parameter            :: xmin = -750e3_dp
    real(dp), parameter            :: xmax =  750e3_dp
    real(dp), parameter            :: ymin = -750e3_dp
    real(dp), parameter            :: ymax =  750e3_dp
    type(type_mesh)                :: mesh

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Create a simple test mesh
    call allocate_mesh_primary( mesh, 'dummy_mesh_1', 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = 50e3_dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call crop_mesh_primary( mesh)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)
    call calc_all_matrix_operators_mesh( mesh)

    ! Run all unit tests
    call test_SMB_idealised ( test_name, mesh)
    call test_SMB_prescribed( test_name, mesh)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_SMB_main

  subroutine test_SMB_idealised( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_SMB_idealised'
    character(len=1024), parameter      :: test_name_local = 'idealised'
    character(len=1024)                 :: test_name
    class(atype_SMB_model), allocatable :: SMB
    type(type_ice_model)     , target   :: ice
    type(type_climate_model) , target   :: climate
    type(type_grid)          , target   :: grid_smooth
    real(dp)                            :: time
    integer                             :: vi
    real(dp)                            :: dH_dt
    logical                             :: test_result
    integer                             :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Create idealised SMB model
    call create_SMB_model( SMB, 'idealised')
    call SMB%allocate  ( SMB%ct_allocate( 'SMB_idealised', 'ANT', mesh))
    call SMB%initialise( SMB%ct_initialise( ice))

    ! Run idealised SMB model for static Halfar solution
    C%choice_SMB_model_idealised = 'Halfar_static'
    C%refgeo_idealised_Halfar_H0 = 3000._dp
    C%refgeo_idealised_Halfar_R0 = 500e3_dp
    C%uniform_Glens_flow_factor  = 1e-16_dp
    time = 0._dp
    call SMB%run( SMB%ct_run( time, ice, climate, grid_smooth))

    ! Verify that it worked
    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      if (hypot( mesh%V( vi,1), mesh%V( vi,2)) < C%refgeo_idealised_Halfar_R0 * 0.9_dp) then
        ! Because there are some exceptions for the ice margin and ice-free area...
        dH_dt = Halfar%dH_dt( C%uniform_Glens_flow_factor, C%Glens_flow_law_exponent, &
          C%refgeo_idealised_Halfar_H0, C%refgeo_idealised_Halfar_R0, &
          mesh%V( vi,1), mesh%V( vi,2), time)
        test_result = test_result .and. test_tol( SMB%SMB(vi), -dH_dt, 1e-5_dp)
      end if
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call SMB%deallocate

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_SMB_idealised

  subroutine test_SMB_prescribed( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_SMB_prescribed'
    character(len=1024), parameter              :: test_name_local = 'prescribed'
    character(len=1024)                         :: test_name
    real(dp), dimension(:), contiguous, pointer :: SMB_ref => null()
    type(MPI_WIN)                               :: wSMB_ref
    integer                                     :: vi
    character(:), allocatable                   :: filename
    integer                                     :: ncid
    class(atype_SMB_model), allocatable         :: SMB
    type(type_ice_model)     , target           :: ice
    type(type_climate_model) , target           :: climate
    type(type_grid)          , target           :: grid_smooth
    real(dp)                                    :: time
    logical                                     :: test_result
    integer                                     :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Set up reference SMB and T2m fields
    call allocate_dist_shared( SMB_ref, wSMB_ref, mesh%pai_V%n_nih)
    do vi = mesh%vi1, mesh%vi2
      SMB_ref( vi) = hypot( mesh%V( vi,1), mesh%V( vi,2))
    end do

    ! Write these to a temporary NetCDF file
    filename = trim( C%output_dir) // '/test_SMB_prescribed_SMB_ref.nc'
    call create_new_netcdf_file_for_writing( filename, ncid)
    call setup_mesh_in_netcdf_file( filename, ncid, mesh)
    call add_field_mesh_dp_2D_notime( filename, ncid, 'SMB')
    call write_to_field_multopt_mesh_dp_2D_notime( mesh, filename, ncid, 'SMB', SMB_ref)
    call close_netcdf_file( ncid)

    ! Create prescribed SMB model
    C%choice_SMB_prescribed_ANT    = 'SMB_no_time'
    C%filename_SMB_prescribed_ANT  = filename
    C%timeframe_SMB_prescribed_ANT = 1e9_dp

    call create_SMB_model( SMB, 'prescribed')
    call SMB%allocate  ( SMB%ct_allocate( 'SMB_prescribed', 'ANT', mesh))
    call SMB%initialise( SMB%ct_initialise( ice))

    ! Verify that it worked
    test_result = .true.
    do vi = mesh%vi1, mesh%vi2
      test_result = test_result .and. test_tol( SMB%SMB(vi), SMB_ref( vi), 1e-5_dp)
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, test_name)

    ! Clean up after yourself
    call SMB%deallocate
    call deallocate_dist_shared( SMB_ref, wSMB_ref)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_SMB_prescribed

end module ut_SMB