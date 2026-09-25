module ut_climate_retreat_mask

  ! Unit tests for the prescribed ice-shelf retreat mask

  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD
  use netcdf, only: NF90_CREATE, NF90_NETCDF4, NF90_CLOBBER, NF90_DEF_DIM, NF90_DEF_VAR, NF90_PUT_ATT, &
    NF90_ENDDEF, NF90_PUT_VAR, NF90_CLOSE, NF90_DOUBLE, NF90_FLOAT, NF90_INT, NF90_INT64, NF90_UNLIMITED, &
    NF90_NOERR, NF90_STRERROR
  use precisions, only: dp, int8
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use ut_basic, only: unit_test, foldername_unit_tests_output
  use tests_main, only: test_tol
  use model_configuration, only: C
  use parameters, only: pi
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary, crop_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use mesh_disc_calc_matrix_operators_2D, only: calc_all_matrix_operators_mesh
  use ice_geometry_model_basic, only: type_ice_geometry_model
  use climate_model_types, only: type_climate_model
  use netcdf_io_main, only: read_time_from_file
  use climate_retreat_mask, only: initialise_climate_retreat_mask, run_climate_retreat_mask, &
    remap_climate_retreat_mask, resolve_region_filename

  implicit none

  private

  public :: unit_tests_climate_retreat_mask_main

  ! Irregular time axis and spatially uniform mask values of the test frames
  real(dp), dimension(3), parameter :: times_test  = [2000._dp, 2007._dp, 2020._dp]
  real(dp), dimension(3), parameter :: values_test = [0._dp, 0.35_dp, 1._dp]

contains

  subroutine unit_tests_climate_retreat_mask_main( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'unit_tests_climate_retreat_mask_main'
    character(len=1024), parameter :: test_name_local = 'climate_retreat_mask'
    character(len=1024)            :: test_name
    real(dp)                       :: alpha_min, res_max
    real(dp), parameter            :: xmin = -500e3_dp
    real(dp), parameter            :: xmax =  500e3_dp
    real(dp), parameter            :: ymin = -500e3_dp
    real(dp), parameter            :: ymax =  500e3_dp
    type(type_mesh)                :: mesh

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Create a simple test mesh
    call allocate_mesh_primary( mesh, 'retreat_mask_test_mesh', 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)
    alpha_min = 25._dp * pi / 180._dp
    res_max = 50e3_dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call crop_mesh_primary( mesh)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)
    call calc_all_matrix_operators_mesh( mesh)

    ! Run all unit tests
    call test_time_axis_types  ( test_name)
    call test_transient_mask   ( test_name, mesh)
    call test_single_frame_mask( test_name, mesh)
    call test_open_ocean_reference( test_name, mesh)
    call test_region_filename     ( test_name)

    ! Leave the retreat mask switched off for any later tests
    C%do_use_ISMIP_future_shelf_collapse_forcing = .false.

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine unit_tests_climate_retreat_mask_main

  subroutine test_time_axis_types( test_name_parent)
    ! Time axes of all supported types must arrive identically on every process

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'test_time_axis_types'
    character(len=1024)                 :: test_name, filename
    character(len=6), dimension(4)      :: type_names = ['double', 'float ', 'int   ', 'int64 ']
    integer,          dimension(4)      :: types
    real(dp), dimension(:), allocatable :: times
    logical                             :: test_result
    integer                             :: i, ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    types = [NF90_DOUBLE, NF90_FLOAT, NF90_INT, NF90_INT64]

    do i = 1, size( types)
      test_name = trim( test_name_parent) // '/read_time_' // trim( type_names( i))
      filename = trim( foldername_unit_tests_output) // '/retreat_mask_time_' // trim( type_names( i)) // '.nc'
      if (par%primary) call write_mask_file( filename, types( i), times_test, values_test)
      call sync

      call read_time_from_file( filename, times)
      test_result = size( times) == size( times_test)
      if (test_result) test_result = all( times == times_test)
      call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
      call unit_test( test_result, test_name)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_time_axis_types

  subroutine test_transient_mask( test_name_parent, mesh)
    ! Linear interpolation between irregular frames, cached frames, endpoint holding, and remapping

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'test_transient_mask'
    character(len=1024)                    :: test_name
    type(type_ice_geometry_model)          :: geom
    type(type_climate_model)               :: climate
    real(dp), dimension(7), parameter      :: t_run  = [1990._dp, 2000._dp, 2003.5_dp, 2007._dp, 2013.5_dp, 2020._dp, 2050._dp]
    real(dp), dimension(7), parameter      :: expect = [0._dp, 0._dp, 0.175_dp, 0.35_dp, 0.675_dp, 1._dp, 1._dp]
    logical                                :: test_result
    integer                                :: i, ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    test_name = trim( test_name_parent) // '/transient'

    call geom%allocate( 'ANT', mesh)
    call set_retreat_config( trim( foldername_unit_tests_output) // '/retreat_mask_time_double.nc', .false., .false.)
    call initialise_climate_retreat_mask( mesh, geom, climate, 'ANT')

    test_result = .true.
    do i = 1, size( t_run)
      call run_climate_retreat_mask( mesh, climate, t_run( i))
      test_result = test_result .and. all_close( climate%retreat%mask, expect( i))
      ! Frames bracket the time within the file interval
      if (t_run( i) == 2003.5_dp) test_result = test_result .and. &
        climate%retreat%t0 == 2000._dp .and. climate%retreat%t1 == 2007._dp
      if (t_run( i) == 2013.5_dp) test_result = test_result .and. &
        climate%retreat%t0 == 2007._dp .and. climate%retreat%t1 == 2020._dp
    end do
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '/interpolation_and_holding')

    ! Remapping reloads both frames
    call remap_climate_retreat_mask( mesh, climate, 2013.5_dp)
    test_result = all_close( climate%retreat%mask, 0.675_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '/remap')

    call geom%deallocate()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_transient_mask

  subroutine test_single_frame_mask( test_name_parent, mesh)
    ! A time axis with a single frame is held at all times

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_single_frame_mask'
    character(len=1024)            :: filename
    type(type_ice_geometry_model)  :: geom
    type(type_climate_model)       :: climate
    logical                        :: test_result
    integer                        :: ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    filename = trim( foldername_unit_tests_output) // '/retreat_mask_single_frame.nc'
    if (par%primary) call write_mask_file( filename, NF90_DOUBLE, [2007._dp], [0.35_dp])
    call sync

    call geom%allocate( 'ANT', mesh)
    call set_retreat_config( filename, .false., .false.)
    call initialise_climate_retreat_mask( mesh, geom, climate, 'ANT')

    call run_climate_retreat_mask( mesh, climate, 1990._dp)
    test_result = all_close( climate%retreat%mask, 0.35_dp)
    call run_climate_retreat_mask( mesh, climate, 2050._dp)
    test_result = test_result .and. all_close( climate%retreat%mask, 0.35_dp)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name_parent) // '/single_frame')

    call geom%deallocate()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_single_frame_mask

  subroutine test_open_ocean_reference( test_name_parent, mesh)
    ! The open-ocean reference is derived once, written, and reused after remapping and
    ! in a new run segment, independently of the geometry at that time

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_open_ocean_reference'
    character(len=1024)            :: test_name, filename, reference_filename
    type(type_ice_geometry_model)  :: geom
    type(type_climate_model)       :: climate, climate_restart
    logical                        :: test_result
    integer                        :: vi, ierr

    ! Add routine to call stack
    call init_routine( routine_name)

    test_name = trim( test_name_parent) // '/open_ocean_reference'

    filename = trim( foldername_unit_tests_output) // '/retreat_mask_static.nc'
    if (par%primary) call write_mask_file( filename, -1, [0._dp], [1._dp])
    call sync

    ! Initial geometry: ice-free ocean for x > 0
    call geom%allocate( 'ANT', mesh)
    do vi = mesh%vi1, mesh%vi2
      geom%mask_icefree_ocean( vi) = mesh%V( vi,1) > 0._dp
    end do

    call set_retreat_config( filename, .true., .true.)
    call initialise_climate_retreat_mask( mesh, geom, climate, 'ANT')
    call run_climate_retreat_mask( mesh, climate, 0._dp)
    test_result = mask_matches_reference( mesh, climate)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '/derived')

    ! The geometry changes, but remapping reads the written reference
    reference_filename = climate%retreat%open_ocean_reference_filename
    geom%mask_icefree_ocean( mesh%vi1:mesh%vi2) = .true.
    call remap_climate_retreat_mask( mesh, climate, 1._dp)
    test_result = mask_matches_reference( mesh, climate)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '/remap')

    ! A new run segment given the written reference ignores its own initial geometry
    C%retreat_mask_open_ocean_reference_filename = reference_filename
    call initialise_climate_retreat_mask( mesh, geom, climate_restart, 'ANT')
    call run_climate_retreat_mask( mesh, climate_restart, 1._dp)
    test_result = mask_matches_reference( mesh, climate_restart)
    call MPI_ALLREDUCE( MPI_IN_PLACE, test_result, 1, MPI_LOGICAL, MPI_LAND, MPI_COMM_WORLD, ierr)
    call unit_test( test_result, trim( test_name) // '/restart')

    C%retreat_mask_open_ocean_reference_filename = ''
    call geom%deallocate()

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_open_ocean_reference

  subroutine test_region_filename( test_name_parent)
    ! The token {region} in the reference file name selects the file of each model region

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_region_filename'
    logical                        :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    test_result = &
      resolve_region_filename( 'run/reference_{region}.nc', 'GRL') == 'run/reference_GRL.nc' .and. &
      resolve_region_filename( '{region}/reference_{region}.nc', 'ANT') == 'ANT/reference_ANT.nc' .and. &
      resolve_region_filename( 'run/reference_ANT.nc', 'GRL') == 'run/reference_ANT.nc'
    call unit_test( test_result, trim( test_name_parent) // '/region_filename')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_region_filename

! ===== Utilities =====
! =====================

  subroutine set_retreat_config( filename, static, open_ocean)

    ! In/output variables:
    character(len=*), intent(in) :: filename
    logical,          intent(in) :: static, open_ocean

    C%do_use_ISMIP_future_shelf_collapse_forcing   = .true.
    C%ISMIP_future_shelf_collapse_forcing_filename = filename
    C%shelf_collapse_type                          = 'calving'
    C%retreat_mask_without_time                    = static
    C%retreat_mask_applied_only_to_open_ocean      = open_ocean
    C%start_time_of_run                            = 2000._dp
    C%end_time_of_run                              = 2020._dp

  end subroutine set_retreat_config

  function all_close( mask, value) result( res)

    ! In/output variables:
    real(dp), dimension(:), intent(in) :: mask
    real(dp),               intent(in) :: value
    logical                            :: res

    ! Local variables:
    integer :: i

    res = .true.
    do i = 1, size( mask)
      res = res .and. test_tol( mask( i), value, 1e-12_dp)
    end do

  end function all_close

  function mask_matches_reference( mesh, climate) result( res)
    ! A uniform static mask of one is restricted to the ice-free ocean for x > 0

    ! In/output variables:
    type(type_mesh),          intent(in) :: mesh
    type(type_climate_model), intent(in) :: climate
    logical                              :: res

    ! Local variables:
    integer :: vi

    res = .true.
    do vi = mesh%vi1, mesh%vi2
      if (mesh%V( vi,1) > 0._dp) then
        res = res .and. climate%retreat%open_ocean( vi) .and. test_tol( climate%retreat%mask( vi), 1._dp, 1e-12_dp)
      else
        res = res .and. (.not. climate%retreat%open_ocean( vi)) .and. climate%retreat%mask( vi) == 0._dp
      end if
    end do

  end function mask_matches_reference

  subroutine write_mask_file( filename, time_type, times, values)
    ! Write a mask file on an x/y-grid covering the test mesh, with spatially uniform
    ! values per frame. A negative time_type writes a static mask without time.

    ! In/output variables:
    character(len=*),       intent(in) :: filename
    integer,                intent(in) :: time_type
    real(dp), dimension(:), intent(in) :: times, values

    ! Local variables:
    integer, parameter                          :: n = 25
    integer                                     :: ncid, id_x, id_y, id_t, id_vx, id_vy, id_vt, id_m, i
    real(dp), dimension(n)                      :: x
    real(dp), dimension(:,:,:), allocatable     :: d

    do i = 1, n
      x( i) = -600e3_dp + real( i-1, dp) * 50e3_dp
    end do

    call check_nc( NF90_CREATE( trim( filename), ior( NF90_NETCDF4, NF90_CLOBBER), ncid), filename)
    call check_nc( NF90_DEF_DIM( ncid, 'x', n, id_x), filename)
    call check_nc( NF90_DEF_DIM( ncid, 'y', n, id_y), filename)
    call check_nc( NF90_DEF_VAR( ncid, 'x', NF90_DOUBLE, [id_x], id_vx), filename)
    call check_nc( NF90_DEF_VAR( ncid, 'y', NF90_DOUBLE, [id_y], id_vy), filename)
    if (time_type < 0) then
      call check_nc( NF90_DEF_VAR( ncid, 'mask', NF90_DOUBLE, [id_x, id_y], id_m), filename)
    else
      call check_nc( NF90_DEF_DIM( ncid, 'time', NF90_UNLIMITED, id_t), filename)
      call check_nc( NF90_DEF_VAR( ncid, 'time', time_type, [id_t], id_vt), filename)
      call check_nc( NF90_PUT_ATT( ncid, id_vt, 'units', 'years'), filename)
      call check_nc( NF90_DEF_VAR( ncid, 'mask', NF90_DOUBLE, [id_x, id_y, id_t], id_m), filename)
    end if
    call check_nc( NF90_ENDDEF( ncid), filename)

    call check_nc( NF90_PUT_VAR( ncid, id_vx, x), filename)
    call check_nc( NF90_PUT_VAR( ncid, id_vy, x), filename)
    allocate( d( n, n, size( values)))
    do i = 1, size( values)
      d( :,:,i) = values( i)
    end do
    if (time_type < 0) then
      call check_nc( NF90_PUT_VAR( ncid, id_m, d( :,:,1)), filename)
    else
      select case (time_type)
      case (NF90_INT)
        call check_nc( NF90_PUT_VAR( ncid, id_vt, nint( times)), filename)
      case (NF90_INT64)
        call check_nc( NF90_PUT_VAR( ncid, id_vt, int( times, int8)), filename)
      case default
        call check_nc( NF90_PUT_VAR( ncid, id_vt, times), filename)
      end select
      call check_nc( NF90_PUT_VAR( ncid, id_m, d), filename)
    end if
    call check_nc( NF90_CLOSE( ncid), filename)

  end subroutine write_mask_file

  subroutine check_nc( stat, filename)
    ! Stop the unit tests if writing a test input file fails

    ! In/output variables:
    integer,          intent(in) :: stat
    character(len=*), intent(in) :: filename

    if (stat /= NF90_NOERR) call crash('writing test file "' // trim( filename) // '" failed: ' // &
      trim( NF90_STRERROR( stat)))

  end subroutine check_nc

end module ut_climate_retreat_mask
