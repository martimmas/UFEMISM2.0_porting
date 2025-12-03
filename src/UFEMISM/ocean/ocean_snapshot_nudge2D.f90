module ocean_snapshot_nudge2D

  use precisions, only: dp
  use parameters, only: NaN
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning, colour_string
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use ocean_model_types, only: type_ocean_model_snapshot_nudge2D
  use netcdf_io_main
  use remapping_main

  implicit none

  private

  public :: initialise_ocean_model_snapshot_nudge2D, run_ocean_model_snapshot_nudge2D

contains

  subroutine run_ocean_model_snapshot_nudge2D( mesh, ice, snapshot_nudge2D, region_name)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ice_model),                    intent(in   ) :: ice
    type(type_ocean_model_snapshot_nudge2D), intent(inout) :: snapshot_nudge2D
    character(len=3),                        intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'run_ocean_model_snapshot_nudge2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine run_ocean_model_snapshot_nudge2D

  subroutine initialise_ocean_model_snapshot_nudge2D( mesh, snapshot_nudge2D, region_name)

    ! In/output variables:
    type(type_mesh),                         intent(in   ) :: mesh
    type(type_ocean_model_snapshot_nudge2D), intent(inout) :: snapshot_nudge2D
    character(len=3),                        intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'initialise_ocean_model_snapshot_nudge2D'
    character(len=1024)                   :: filename
    integer                               :: ncid
    real(dp), dimension(:,:), allocatable :: T_ref_partial_raw_layers
    real(dp), dimension(:,:), allocatable :: S_ref_partial_raw_layers

    ! Add routine to call stack
    call init_routine( routine_name)

    select case (region_name)
      case ('NAM')
        filename = C%filename_ocean_snapshot_NAM
      case ('EAS')
        filename = C%filename_ocean_snapshot_EAS
      case ('GRL')
        filename = C%filename_ocean_snapshot_GRL
      case ('ANT')
        filename = C%filename_ocean_snapshot_ANT
      case default
        call crash('unknown region_name "' // region_name // '"')
    end select

    ! Allocate memory for meshed ocean data
    allocate( snapshot_nudge2D%T_ref       ( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate( snapshot_nudge2D%S_ref       ( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate( snapshot_nudge2D%deltaT_nudge( mesh%vi1:mesh%vi2            ), source = NaN)
    allocate( snapshot_nudge2D%T           ( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)
    allocate( snapshot_nudge2D%S           ( mesh%vi1:mesh%vi2, C%nz_ocean), source = NaN)

    ! Set up the grid from the file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call setup_xy_grid_from_file( filename, ncid, snapshot_nudge2D%grid_ref)
    call setup_depth_from_file( filename, ncid, snapshot_nudge2D%ndepth_ref, snapshot_nudge2D%depth_ref)
    call close_netcdf_file( ncid)

    ! Allocate memory for gridded 3-D ocean snapshot data
    allocate( snapshot_nudge2D%T_ref_grid( snapshot_nudge2D%grid_ref%n1: snapshot_nudge2D%grid_ref%n2, snapshot_nudge2D%ndepth_ref), source = NaN)
    allocate( snapshot_nudge2D%S_ref_grid( snapshot_nudge2D%grid_ref%n1: snapshot_nudge2D%grid_ref%n2, snapshot_nudge2D%ndepth_ref), source = NaN)
    allocate( snapshot_nudge2D%T_grid    ( snapshot_nudge2D%grid_ref%n1: snapshot_nudge2D%grid_ref%n2, snapshot_nudge2D%ndepth_ref), source = NaN)
    allocate( snapshot_nudge2D%S_grid    ( snapshot_nudge2D%grid_ref%n1: snapshot_nudge2D%grid_ref%n2, snapshot_nudge2D%ndepth_ref), source = NaN)

    ! Read gridded 3-D ocean snapshot data
    call read_field_from_xy_file_dp_3D_ocean( filename, field_name_options_T_ocean, snapshot_nudge2D%T_ref_grid)
    call read_field_from_xy_file_dp_3D_ocean( filename, field_name_options_S_ocean, snapshot_nudge2D%S_ref_grid)

    ! Allocate memory for 3-D ocean snapshot data on the mesh, but with the original vertical layers
    allocate( T_ref_partial_raw_layers( mesh%vi1:mesh%vi2, snapshot_nudge2D%ndepth_ref))
    allocate( S_ref_partial_raw_layers( mesh%vi1:mesh%vi2, snapshot_nudge2D%ndepth_ref))

    ! Remap data horizontally
    call map_from_xy_grid_to_mesh_3D( snapshot_nudge2D%grid_ref, mesh, trim( C%output_dir), snapshot_nudge2D%T_ref_grid, T_ref_partial_raw_layers)
    call map_from_xy_grid_to_mesh_3D( snapshot_nudge2D%grid_ref, mesh, trim( C%output_dir), snapshot_nudge2D%S_ref_grid, S_ref_partial_raw_layers)

    ! Remap data vertically
    call map_from_vertical_to_vertical_2D_ocean( mesh, snapshot_nudge2D%depth_ref, C%z_ocean, T_ref_partial_raw_layers, snapshot_nudge2D%T_ref)
    call map_from_vertical_to_vertical_2D_ocean( mesh, snapshot_nudge2D%depth_ref, C%z_ocean, S_ref_partial_raw_layers, snapshot_nudge2D%S_ref)

    call create_ocean_model_snapshot_nudge2D_output_file( snapshot_nudge2D)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model_snapshot_nudge2D

  subroutine create_ocean_model_snapshot_nudge2D_output_file( snapshot_nudge2D)

    ! In/output variables:
    type(type_ocean_model_snapshot_nudge2D), intent(inout) :: snapshot_nudge2D

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'create_ocean_model_snapshot_nudge2D_output_file'
    character(len=1024)                   :: filename
    integer                               :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    filename = trim( C%output_dir) // 'ocean_snapshot_nudged.nc'
    snapshot_nudge2D%output_filename = filename

    ! Print to terminal
    if (par%primary) write(0,'(a)') '     Creating ocean snapshot+nudge2D output file "' // &
      colour_string( trim( filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( filename, ncid)

    ! Set up the grid in the file
    call setup_xy_grid_in_netcdf_file( filename, ncid, snapshot_nudge2D%grid_ref)

    ! Add a time dimension to the file
    call add_time_dimension_to_file( filename, ncid)

    ! Add a depth dimension to the file
    call add_depth_dimension_to_file( filename, ncid, snapshot_nudge2D%depth_ref)

    ! Add the data fields to the file
    call add_field_grid_dp_3D_ocean( filename, ncid, 't_an' , long_name = 'Ocean temperature', units = 'degrees C')
    call add_field_grid_dp_3D_ocean( filename, ncid, 's_an' , long_name = 'Ocean salinity', units = 'PSU')

    ! Close the file
    call close_netcdf_file( ncid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine create_ocean_model_snapshot_nudge2D_output_file

end module ocean_snapshot_nudge2D
