module LADDIE_main_model

  ! The main stand-alone model of LADDIE

! ===== Preamble =====
! ====================

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_IN_PLACE, MPI_INTEGER, MPI_SUM, MPI_WTIME
  use precisions, only: dp
  use mpi_basic, only: par, sync 
  use control_resources_and_error_messaging, only: happy, warning, crash, init_routine, finalise_routine, &
    colour_string, str2int, int2str, insert_val_into_string_dp
  use model_configuration                                    , ONLY: C
  use parameters
  use reference_geometry_types, only: type_reference_geometry
  use laddie_model_types, only: type_laddie_model
  use laddie_forcing_types, only: type_laddie_forcing
  use laddie_main, only: run_laddie_model
  use mesh_types, only: type_mesh
  use ocean_main, only: initialise_ocean_model, run_ocean_model
  use netcdf_io_main
  use mesh_creation_main, only: create_mesh_from_gridded_geometry, create_mesh_from_meshed_geometry, write_mesh_success
  use grid_basic, only: setup_square_grid
  use mesh_output_files, only: create_main_regional_output_file_mesh, write_to_main_regional_output_file_mesh
  use grid_output_files, only: create_main_regional_output_file_grid, write_to_main_regional_output_file_grid, &
    create_main_regional_output_file_grid_ROI, write_to_main_regional_output_file_grid_ROI
  use scalar_output_files, only: create_scalar_regional_output_file, buffer_scalar_output, write_to_scalar_regional_output_file
  use mesh_ROI_polygons
  use apply_maps, only: clear_all_maps_involving_this_mesh
  use mesh_memory, only: deallocate_mesh

  implicit none

contains

! ===== Main routine =====
! ========================

  subroutine run_laddie_standalone( laddie, t_end, forcing)
    ! Integrate the model until t_end

    implicit none

    ! In/output variables:
    type(type_laddie_model)                            , intent(inout) :: laddie
    real(dp)                                           , intent(in)    :: t_end
    type(type_laddie_forcing)                          , intent(in)    :: forcing

    ! Local variables:
    character(len=256)                                                 :: routine_name

    ! Add routine to path
    routine_name = 'run_laddie_standalone'  
    call init_routine( routine_name)

    if (par%primary) write (0,'(A)') ' Running laddie standalone '

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_laddie_standalone

! ===== Model initialisation =====
! ================================

  subroutine initialise_laddie_standalone( laddie, forcing)
    ! Initialise the model

    implicit none

    ! In/output variables:
    type(type_laddie_model)                            , intent(inout) :: laddie
    type(type_laddie_forcing)                          , intent(inout) :: forcing

    ! Local variables:
    character(len=256), parameter                                      :: routine_name = 'initialise_laddie_standalone'

    ! Add routine to path
    call init_routine( routine_name)

    if (par%primary) write (0,'(A)') ' Initialising laddie standalone '

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_laddie_standalone

end module LADDIE_main_model
