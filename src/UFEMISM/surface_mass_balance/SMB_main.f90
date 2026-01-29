MODULE SMB_main

  ! The main SMB model module.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE grid_basic                                             , ONLY: type_grid
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model
  USE SMB_model_types                                        , ONLY: type_SMB_model
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  use mesh_ROI_polygons, only: calc_polygon_Patagonia
  use plane_geometry, only: is_in_polygon
  use mesh_data_smoothing, only: smooth_Gaussian
  use netcdf_io_main

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_SMB_model( mesh, grid_smooth, ice, climate, SMB, region_name, time)
    ! Calculate the surface mass balance

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_grid),            target,     INTENT(IN)    :: grid_smooth
    TYPE(type_ice_model),       target,     INTENT(IN)    :: ice
    TYPE(type_climate_model),   target,     INTENT(IN)    :: climate
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if we need to calculate a new SMB
    IF (C%do_asynchronous_SMB) THEN
      ! Asynchronous coupling: do not calculate a new SMB in
      ! every model loop, but only at its own separate time step

      ! Check if this is the next SMB time step
      IF (time == SMB%t_next) THEN
        ! Go on to calculate a new SMB
        SMB%t_next = time + C%dt_SMB
      ELSEIF (time > SMB%t_next) THEN
        ! This should not be possible
        CALL crash('overshot the SMB time step')
      ELSE
        ! It is not yet time to calculate a new SMB
        CALL finalise_routine( routine_name)
        RETURN
      END IF

    ELSE ! IF (C%do_asynchronous_SMB) THEN
      ! Synchronous coupling: calculate a new SMB in every model loop
      SMB%t_next = time + C%dt_SMB
    END IF

    ! Determine which SMB model to run for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Run the chosen SMB model
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        SMB%SMB = C%uniform_SMB
      CASE ('idealised')
        call SMB%idealised%run( SMB%idealised%ct_run( time, ice, climate, grid_smooth))
        SMB%SMB( mesh%vi1:mesh%vi2) = SMB%idealised%SMB( mesh%vi1:mesh%vi2)
      CASE ('prescribed')
        call SMB%prescribed%run( SMB%prescribed%ct_run( time, ice, climate, grid_smooth))
        SMB%SMB( mesh%vi1:mesh%vi2) = SMB%prescribed%SMB( mesh%vi1:mesh%vi2)
      CASE ('reconstructed')
        call SMB%reconstructed%run( SMB%reconstructed%ct_run( time, ice, climate, grid_smooth))
        SMB%SMB( mesh%vi1:mesh%vi2) = SMB%reconstructed%SMB( mesh%vi1:mesh%vi2)
      CASE ('snapshot_plus_anomalies')
        call SMB%snapshot_plus_anomalies%run( SMB%snapshot_plus_anomalies%ct_run( time, ice, climate, grid_smooth))
        SMB%SMB( mesh%vi1:mesh%vi2) = SMB%snapshot_plus_anomalies%SMB( mesh%vi1:mesh%vi2)
      CASE ('IMAU-ITM')
        call SMB%IMAUITM%run( SMB%IMAUITM%ct_run( time, ice, climate, grid_smooth))
        SMB%SMB( mesh%vi1:mesh%vi2) = SMB%IMAUITM%SMB( mesh%vi1:mesh%vi2)
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model

  SUBROUTINE initialise_SMB_model( mesh, ice, SMB, region_name)
    ! Initialise the SMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),      target,      INTENT(IN)    :: ice
    TYPE(type_SMB_model),                   INTENT(OUT)   :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising surface mass balance model...'

    ! Determine which SMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Allocate memory for main variables
    ALLOCATE( SMB%SMB( mesh%vi1:mesh%vi2))
    SMB%SMB = 0._dp

    ! Set time of next calculation to start time
    SMB%t_next = C%start_time_of_run

    ! Determine which SMB model to initialise
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        SMB%SMB = C%uniform_SMB
      CASE ('idealised')
        call SMB%idealised%allocate  ( SMB%idealised%ct_allocate( 'SMB_idealised', region_name, mesh))
        call SMB%idealised%initialise( SMB%idealised%ct_initialise( ice))
      CASE ('prescribed')
        call SMB%prescribed%allocate  ( SMB%prescribed%ct_allocate( 'SMB_prescribed', region_name, mesh))
        call SMB%prescribed%initialise( SMB%prescribed%ct_initialise( ice))
      CASE ('reconstructed')
        call SMB%reconstructed%allocate  ( SMB%reconstructed%ct_allocate( 'SMB_reconstructed', region_name, mesh))
        call SMB%reconstructed%initialise( SMB%reconstructed%ct_initialise( ice))
      CASE ('snapshot_plus_anomalies')
        call SMB%snapshot_plus_anomalies%allocate  ( SMB%snapshot_plus_anomalies%ct_allocate( 'SMB_snapshot_plus_anomalies', region_name, mesh))
        call SMB%snapshot_plus_anomalies%initialise( SMB%snapshot_plus_anomalies%ct_initialise( ice))
      CASE ('IMAU-ITM')
        call SMB%IMAUITM%allocate  ( SMB%IMAUITM%ct_allocate( 'SMB_IMAU_ITM', region_name, mesh))
        call SMB%IMAUITM%initialise( SMB%IMAUITM%ct_initialise( ice))
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SMB_model

  SUBROUTINE write_to_restart_file_SMB_model( mesh, SMB, region_name, time)
    ! Write to the restart file for the SMB model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                   INTENT(IN)    :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'write_to_restart_file_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which SMB model to use for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Write to the restart file of the chosen SMB model
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('prescribed')
        ! No need to do anything
      CASE ('reconstructed')
        ! No need to do anything
      CASE ('IMAU-ITM')
        call write_to_restart_file_SMB_model_region(mesh, SMB, region_name, time)
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE write_to_restart_file_SMB_model

  SUBROUTINE write_to_restart_file_SMB_model_region(mesh, SMB, region_name, time)
    ! Write to the restart NetCDF file for the SMB model

    ! In/output variables:
    TYPE(type_mesh),          INTENT(IN) :: mesh
    TYPE(type_SMB_model),     INTENT(IN) :: SMB
    CHARACTER(LEN=3),         INTENT(IN) :: region_name
    REAL(dp),                 INTENT(IN) :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER        :: routine_name = 'write_to_restart_file_SMB_model_region'
    INTEGER                              :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Print to terminal
    IF (par%primary) WRITE(0,'(A)') '   Writing to SMB restart file "' // &
      colour_string( TRIM( SMB%restart_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    CALL open_existing_netcdf_file_for_writing( SMB%restart_filename, ncid)

    ! Write the time to the file
    CALL write_time_to_file( SMB%restart_filename, ncid, time)
    ! month dimension is already written when adding it to file

    ! ! Write the SMB fields to the file
    ! TODO: do we need to check if IMAUITM is being used before writing these files?
    CALL write_to_field_multopt_mesh_dp_2D_monthly( mesh, SMB%restart_filename, ncid, 'SMB_monthly', SMB%IMAUITM%SMB_monthly)
    CALL write_to_field_multopt_mesh_dp_2D_monthly( mesh, SMB%restart_filename, ncid, 'FirnDepth', SMB%IMAUITM%FirnDepth)
    CALL write_to_field_multopt_mesh_dp_2D(         mesh, SMB%restart_filename, ncid, 'MeltPreviousYear', SMB%IMAUITM%MeltPreviousYear)
    CALL write_to_field_multopt_mesh_dp_2D(         mesh, SMB%restart_filename, ncid, 'SMB', SMB%SMB)

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)
  END SUBROUTINE write_to_restart_file_SMB_model_region

  SUBROUTINE create_restart_file_SMB_model( mesh, SMB, region_name)
    ! Create the restart file for the SMB model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_file_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine which SMB model to use for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Create the restart file of the chosen SMB model
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        ! No need to do anything
      CASE ('prescribed')
        ! No need to do anything
      CASE ('reconstructed')
        ! No need to do anything
      CASE ('IMAU-ITM')
        call create_restart_file_SMB_model_region(mesh, SMB, region_name)
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_SMB_model

  SUBROUTINE create_restart_file_SMB_model_region(mesh, SMB, region_name)
    ! Create the restart file for the SMB model

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_SMB_model),                   INTENT(INOUT) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'create_restart_file_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model
    CHARACTER(LEN=256)                                    :: filename_base
    INTEGER                                               :: ncid

    ! Add routine to path
    CALL init_routine( routine_name)

    ! If no NetCDF output should be created, do nothing
    IF (.NOT. C%do_create_netcdf_output) THEN
      CALL finalise_routine( routine_name)
      RETURN
    END IF

    ! Set the filename
    filename_base = TRIM( C%output_dir) // 'restart_SMB_' // region_name
    CALL generate_filename_XXXXXdotnc( filename_base, SMB%restart_filename)

    ! Print to terminal
    IF (par%primary) WRITE(0,'(A)') '   Creating SMB model restart file "' // &
      colour_string( TRIM( SMB%restart_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    CALL create_new_netcdf_file_for_writing( SMB%restart_filename, ncid)

    ! Set up the mesh in the file
    CALL setup_mesh_in_netcdf_file( SMB%restart_filename, ncid, mesh)

    ! Add a time dimension to the file
    CALL add_month_dimension_to_file( SMB%restart_filename, ncid)
    CALL add_time_dimension_to_file( SMB%restart_filename, ncid)


    ! Add the data fields to the file
    CALL add_field_mesh_dp_2D( SMB%restart_filename, ncid, 'SMB',                 long_name = 'Surface mass balance',            units = 'm/yr')
    CALL add_field_mesh_dp_2D_monthly( SMB%restart_filename, ncid, 'SMB_monthly', long_name = 'monthly Surface mass balance',    units = 'm/yr')
    CALL add_field_mesh_dp_2D_monthly( SMB%restart_filename, ncid, 'FirnDepth',   long_name = 'Firn Depth',                      units = 'm')
    CALL add_field_mesh_dp_2D( SMB%restart_filename, ncid, 'MeltPreviousYear',    long_name = 'Total melt in the previous year', units = 'm w.e.')

    ! Close the file
    CALL close_netcdf_file( ncid)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE create_restart_file_SMB_model_region

  SUBROUTINE remap_SMB_model( mesh_old, mesh_new, SMB, time, region_name)
    ! Remap the SMB model

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_old
    TYPE(type_mesh),                        INTENT(IN)    :: mesh_new
    TYPE(type_SMB_model),                   INTENT(inout) :: SMB
    real(dp),                               intent(in)    :: time
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'remap_SMB_model'
    CHARACTER(LEN=256)                                    :: choice_SMB_model

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '    Remapping surface mass balance model data to the new mesh...'

    ! Determine which SMB model to initialise for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_model = C%choice_SMB_model_NAM
      CASE ('EAS')
        choice_SMB_model = C%choice_SMB_model_EAS
      CASE ('GRL')
        choice_SMB_model = C%choice_SMB_model_GRL
      CASE ('ANT')
        choice_SMB_model = C%choice_SMB_model_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // region_name // '"')
    END SELECT

    ! Reallocate memory for main variables
    CALL reallocate_bounds( SMB%SMB, mesh_new%vi1, mesh_new%vi2)

    ! Determine which SMB model to initialise
    SELECT CASE (choice_SMB_model)
      CASE ('uniform')
        ! No need to do anything
      CASE ('idealised')
        call SMB%idealised%remap( SMB%idealised%ct_remap( mesh_new, time))
      CASE ('prescribed')
        call SMB%prescribed%remap( SMB%prescribed%ct_remap( mesh_new, time))
      CASE ('reconstructed')
        call SMB%reconstructed%remap( SMB%reconstructed%ct_remap( mesh_new, time))
      CASE ('snapshot_plus_anomalies')
        call SMB%snapshot_plus_anomalies%remap( SMB%snapshot_plus_anomalies%ct_remap( mesh_new, time))
      CASE ('IMAU-ITM')
        call SMB%IMAUITM%remap( SMB%IMAUITM%ct_remap( mesh_new, time))
      CASE DEFAULT
        CALL crash('unknown choice_SMB_model "' // TRIM( choice_SMB_model) // '"')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE remap_SMB_model

END MODULE SMB_main
