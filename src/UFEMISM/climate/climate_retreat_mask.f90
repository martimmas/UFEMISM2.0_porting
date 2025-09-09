module climate_retreat_mask

! ===== Preamble =====
! ====================

  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, sync
  use control_resources_and_error_messaging                  , only: crash, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C
  use parameters
  use mpi_f08, only: MPI_ALLREDUCE, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, MPI_INTEGER
  use mesh_types                                             , only: type_mesh
  use ice_model_types                                        , only: type_ice_model
  use grid_types                                             , only: type_grid
  use climate_model_types                                    , only: type_climate_model, type_climate_model_matrix, type_climate_model_snapshot
  use global_forcing_types                                   , only: type_global_forcing
  use SMB_model_types, only: type_SMB_model
  use climate_realistic                                      , only: initialise_climate_model_realistic, initialise_insolation_forcing, remap_snapshot
  use reallocate_mod                                         , only: reallocate_bounds
  use netcdf_io_main
  use mesh_data_smoothing, only: smooth_Gaussian
  use SMB_IMAU_ITM, only: run_SMB_model_IMAUITM, initialise_SMB_model_IMAUITM
  use climate_matrix_utilities, only: allocate_climate_snapshot, read_climate_snapshot, adapt_precip_CC, adapt_precip_Roe, get_insolation_at_time
  use assertions_basic, only: assert

 implicit none

  private

  public :: run_climate_model_matrix
  public :: initialise_climate_matrix
  public :: remap_climate_matrix_model

contains

! == ISMIP-style Antarctica future phase (SMB + aSMB + ST + aST) forcing
! ======================================================================

  SUBROUTINE run_climate_retreat_mask( mesh, climate_matrix, time, ice)
    ! Run the regional climate model
    !
    ! Use the ISMIP-style (SMB + aSMB + ST + aST) yearly forcing

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                    INTENT(INOUT) :: mesh
    TYPE(type_climate_matrix_regional), INTENT(INOUT) :: climate_matrix
    REAL(dp),                           INTENT(IN)    :: time
    TYPE(type_ice_model),               INTENT(IN)    :: ice

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                           :: routine_name = 'run_climate_model_ISMIP_style_future'
    REAL(dp)                                                :: wt0, wt1
    INTEGER                                                 :: vi

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time < climate_matrix%ISMIP_style%t0 .OR. time > climate_matrix%ISMIP_style%t1) THEN

      ! Find and read the two global time frames
      CALL sync
      CALL update_ISMIP_style_future_timeframes( mesh, climate_matrix, time)

    END IF ! IF (time >= climate_matrix%SMB_direct%t0 .AND. time <= climate_matrix%SMB_direct%t1) THEN

    ! Interpolate the two timeframes in time
    wt0 = (climate_matrix%ISMIP_style%t1 - time) / (climate_matrix%ISMIP_style%t1 - climate_matrix%ISMIP_style%t0)
    wt1 = 1._dp - wt0

    DO vi = mesh%vi1, mesh%vi2

      climate_matrix%ISMIP_style%aSMB(   vi) = (wt0 * climate_matrix%ISMIP_style%aSMB0(   vi)) + &
                                               (wt1 * climate_matrix%ISMIP_style%aSMB1(   vi))
      climate_matrix%ISMIP_style%aST(    vi) = (wt0 * climate_matrix%ISMIP_style%aST0(    vi)) + &
                                               (wt1 * climate_matrix%ISMIP_style%aST1(    vi))

    END DO
    CALL sync

    ! Apply the anomaly to calculate the applied SMB and temperature
    DO vi = mesh%vi1, mesh%vi2

      climate_matrix%ISMIP_style%SMB( vi) = climate_matrix%ISMIP_style%SMB_ref( vi) + &
                                            climate_matrix%ISMIP_style%aSMB(    vi) * sec_per_year / ice_density

      climate_matrix%ISMIP_style%ST(  vi) = climate_matrix%ISMIP_style%ST_ref(  vi) + &
                                            climate_matrix%ISMIP_style%aST(     vi)

    END DO
    CALL sync

    ! Set the final values in the "applied" climate snapshot
    DO vi = mesh%vi1, mesh%vi2

      climate_matrix%applied%T2m( vi,:) = climate_matrix%ISMIP_style%ST( vi)

    END DO
    CALL sync


    ! == Shelf collapse forcing ==
    ! ============================

    IF (C%do_use_ISMIP_future_shelf_collapse_forcing) THEN

      ! Interpolate the two timeframes in time
      wt0 = (climate_matrix%ISMIP_style%shelf_collapse_mask_t1 - time) / (climate_matrix%ISMIP_style%shelf_collapse_mask_t1 - climate_matrix%ISMIP_style%shelf_collapse_mask_t0)
      wt1 = 1._dp - wt0

      DO vi = mesh%vi1, mesh%vi2

        climate_matrix%ISMIP_style%shelf_collapse_mask( vi) = (wt0 * climate_matrix%ISMIP_style%shelf_collapse_mask0( vi)) + &
                                                              (wt1 * climate_matrix%ISMIP_style%shelf_collapse_mask1( vi))

      END DO
      CALL sync

    END IF ! IF (C%do_use_ISMIP_future_shelf_collapse_forcing) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_ISMIP_style_future

  SUBROUTINE update_climate_retreat_mask( mesh, climate_matrix, time)
    ! Update the two timeframes of the ISMIP-style future yearly (SMB + aSMB + ST + aST) forcing data

    USE netcdf_input_module, ONLY: read_field_from_file_2D

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                          INTENT(INOUT) :: mesh
    TYPE(type_climate_matrix_regional),       INTENT(INOUT) :: climate_matrix
    REAL(dp),                                 INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                           :: routine_name = 'update_ISMIP_style_future_timeframes'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Read timeframes from file
    CALL read_field_from_file_2D( C%ISMIP_future_atmosphere_anomaly_filename, 'smb_anomaly', mesh, climate_matrix%ISMIP_style%aSMB0, 'ANT', time_to_read = REAL( FLOOR( time),dp)        )
    CALL read_field_from_file_2D( C%ISMIP_future_atmosphere_anomaly_filename, 'smb_anomaly', mesh, climate_matrix%ISMIP_style%aSMB1, 'ANT', time_to_read = REAL( FLOOR( time),dp) + 1._dp)
    CALL read_field_from_file_2D( C%ISMIP_future_atmosphere_anomaly_filename, 'ts_anomaly' , mesh, climate_matrix%ISMIP_style%aST0 , 'ANT', time_to_read = REAL( FLOOR( time),dp)        )
    CALL read_field_from_file_2D( C%ISMIP_future_atmosphere_anomaly_filename, 'ts_anomaly' , mesh, climate_matrix%ISMIP_style%aST1 , 'ANT', time_to_read = REAL( FLOOR( time),dp) + 1._dp)


    ! == Shelf collapse forcing ==
    ! ============================

    IF (C%do_use_ISMIP_future_shelf_collapse_forcing) THEN

      ! Read timeframes from file
      CALL read_field_from_file_2D( C%ISMIP_future_shelf_collapse_forcing_filename, 'mask', mesh, climate_matrix%ISMIP_style%shelf_collapse_mask0, 'ANT', time_to_read = REAL( FLOOR( time),dp)        )
      CALL read_field_from_file_2D( C%ISMIP_future_shelf_collapse_forcing_filename, 'mask', mesh, climate_matrix%ISMIP_style%shelf_collapse_mask1, 'ANT', time_to_read = REAL( FLOOR( time),dp) + 1._dp)

    END IF ! IF (C%do_use_ISMIP_future_shelf_collapse_forcing) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 13)

  END SUBROUTINE update_ISMIP_style_future_timeframes

  SUBROUTINE initialise_climate_retreat_mask( mesh, climate_matrix)
    ! Initialise the regional climate model
    !
    ! Use the ISMIP-style (SMB + aSMB + ST + aST) yearly forcing

    USE netcdf_input_module, ONLY: read_field_from_file_2D

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                    INTENT(INOUT) :: mesh
    TYPE(type_climate_matrix_regional), INTENT(INOUT) :: climate_matrix

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                           :: routine_name = 'initialise_climate_model_ISMIP_style_future'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Allocate memory for baseline temperature and SMB
    CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%SMB_ref, climate_matrix%ISMIP_style%wSMB_ref)
    CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%ST_ref , climate_matrix%ISMIP_style%wST_ref )

    ! Read baseline temperature and SMB
    CALL read_field_from_file_2D( C%ISMIP_future_atmosphere_baseline_filename, 'SMB', mesh, climate_matrix%ISMIP_style%SMB_ref, 'ANT')
    CALL read_field_from_file_2D( C%ISMIP_future_atmosphere_baseline_filename, 'ST' , mesh, climate_matrix%ISMIP_style%ST_ref , 'ANT')

    ! Allocate memory for the timestamps of the two timeframes
    CALL allocate_shared_dp_0D( climate_matrix%ISMIP_style%t0, climate_matrix%ISMIP_style%wt0)
    CALL allocate_shared_dp_0D( climate_matrix%ISMIP_style%t1, climate_matrix%ISMIP_style%wt1)

    IF (par%master) THEN
      ! Give impossible values to timeframes, so that the first call to run_climate_model_ISMIP_style
      ! is guaranteed to first read two new timeframes from the NetCDF file
      climate_matrix%ISMIP_style%t0 = C%start_time_of_run - 100._dp
      climate_matrix%ISMIP_style%t1 = C%start_time_of_run - 90._dp
    END IF ! IF (par%master) THEN
    CALL sync

    ! Allocate memory for the two timeframes
    CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%aSMB0  , climate_matrix%ISMIP_style%waSMB0  )
    CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%aST0   , climate_matrix%ISMIP_style%waST0   )

    CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%aSMB1  , climate_matrix%ISMIP_style%waSMB1  )
    CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%aST1   , climate_matrix%ISMIP_style%waST1   )

    ! Allocate memory for the time-interpolated values of aSMB, dSMBdz, ST, and dSTdz
    CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%aSMB   , climate_matrix%ISMIP_style%waSMB   )
    CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%aST    , climate_matrix%ISMIP_style%waST    )

    ! Allocate memory for the applied values of SMB and ST (i.e. after applying the anomaly and elevation correction)
    CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%SMB    , climate_matrix%ISMIP_style%wSMB    )
    CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%ST     , climate_matrix%ISMIP_style%wST     )

    ! Lastly, allocate memory for the "applied" snapshot
    CALL allocate_climate_snapshot_regional( mesh, climate_matrix%applied, name = 'applied')


    ! == Shelf collapse forcing ==
    ! ============================

    IF (C%do_use_ISMIP_future_shelf_collapse_forcing) THEN

      ! Allocate memory
      CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%shelf_collapse_mask0, climate_matrix%ISMIP_style%wshelf_collapse_mask0)
      CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%shelf_collapse_mask1, climate_matrix%ISMIP_style%wshelf_collapse_mask1)
      CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%shelf_collapse_mask , climate_matrix%ISMIP_style%wshelf_collapse_mask )

      ! Allocate memory for the timestamps of the two timeframes
      CALL allocate_shared_dp_0D( climate_matrix%ISMIP_style%shelf_collapse_mask_t0, climate_matrix%ISMIP_style%wshelf_collapse_mask_t0)
      CALL allocate_shared_dp_0D( climate_matrix%ISMIP_style%shelf_collapse_mask_t1, climate_matrix%ISMIP_style%wshelf_collapse_mask_t1)

      IF (par%master) THEN
        ! Give impossible values to timeframes, so that the first call to run_climate_model_ISMIP_style
        ! is guaranteed to first read two new timeframes from the NetCDF file
        climate_matrix%ISMIP_style%shelf_collapse_mask_t0 = C%start_time_of_run - 100._dp
        climate_matrix%ISMIP_style%shelf_collapse_mask_t1 = C%start_time_of_run - 90._dp
      END IF ! IF (par%master) THEN
      CALL sync

    END IF ! IF (C%do_use_ISMIP_future_shelf_collapse_forcing) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name, n_extra_windows_expected = 1)

  END SUBROUTINE initialise_climate_model_ISMIP_style_future









  ! NOTE: this bit was in the old BMB module

    ! implement ISMIP6 shelf collapse forcing as a very high melt rate, just as in ABUMIP
    IF (C%do_use_ISMIP_future_shelf_collapse_forcing) THEN
      DO vi = mesh%vi1, mesh%vi2
        IF (shelf_collapse_mask( vi) > 0.01_dp) THEN
          BMB%BMB_shelf( vi) = -400._dp
        END IF
      END DO
    END IF

end module climate_retreat_mask