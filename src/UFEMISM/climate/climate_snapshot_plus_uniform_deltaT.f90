MODULE climate_snapshot_plus_uniform_deltaT

  ! Snapshot + dT climate model

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string, warning, insert_val_into_string_int,insert_val_into_string_dp
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model, type_climate_model_snapshot
  USE global_forcing_types                                   , ONLY: type_global_forcing
  use climate_realistic                                      , only: initialise_climate_model_realistic, initialise_insolation_forcing, apply_geometry_downscaling_corrections
  USE global_forcings_main
  USE netcdf_io_main
  USE netcdf_basic
  use mpi_distributed_memory, only: distribute_from_primary
  use climate_matrix_utilities, only: allocate_climate_snapshot, read_climate_snapshot,

  IMPLICIT NONE

  private

  public :: run_climate_model_snapshot_plus_uniform_deltaT
  public :: initialise_climate_model_snapshot_plus_uniform_deltaT


CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_climate_model_snapshot_plus_uniform_deltaT( mesh, ice, climate, forcing, time)
    ! Calculate the climate
    !
    ! Use a snapshot plus a prescribed uniform deltaT

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    TYPE(type_global_forcing),              INTENT(IN)    :: forcing
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_climate_model_snapshot_plus_uniform_deltaT'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen snapshot+deltaT climate model

    ! Update temperature and precipitation fields based on the mismatch between
    ! the ice sheet surface elevation in the forcing climate and the model's ice sheet surface elevation
    CALL apply_geometry_downscaling_corrections( mesh, ice, climate)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_snapshot_plus_uniform_deltaT

  SUBROUTINE initialise_climate_model_snapshot_plus_uniform_deltaT( mesh, ice, climate, forcing, region_name)
    ! Initialise the climate model
    !
    ! Use a realistic climate scheme

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    TYPE(type_global_forcing),              INTENT(IN)    :: forcing
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_climate_model_snapshot_plus_uniform_deltaT'
    CHARACTER(LEN=256)                                    :: filename_climate_snapshot
    LOGICAL                                               :: do_lapse_rates
    REAL(dp)                                              :: timeframe_init_insolation

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '     Initialising realistic climate model "' // &
      colour_string( TRIM( C%choice_climate_model_realistic),'light blue') // '"...'

    ! Run the chosen realistic climate model
    climate%snapshot%has_insolation = .FALSE.
    IF (C%choice_climate_model_realistic == 'snapshot') THEN
      ! Read single-time data from external file

      ! Determine which climate model to initialise for this region
      IF     (region_name == 'NAM') THEN
        filename_climate_snapshot                        = C%filename_climate_snapshot_unif_dT_NAM
        climate%snapshot_deltaT%precip_CC_correction     = C%precip_CC_correction_NAM
        climate%snapshot_deltaT%snapshot%lapse_rate_temp = C%lapse_rate_temp_NAM
        IF (C%choice_SMB_model_NAM == 'IMAU-ITM') THEN
           climate%snapshot%has_insolation = .TRUE.
        END IF
      ELSEIF (region_name == 'EAS') THEN
        filename_climate_snapshot                        = C%filename_climate_snapshot_unif_dT_EAS
        climate%snapshot_deltaT%precip_CC_correction     = C%precip_CC_correction_EAS
        climate%snapshot_deltaT%snapshot%lapse_rate_temp = C%lapse_rate_temp_EAS
        IF (C%choice_SMB_model_EAS == 'IMAU-ITM') THEN
           climate%snapshot%has_insolation = .TRUE.
        END IF
      ELSEIF (region_name == 'GRL') THEN
        filename_climate_snapshot                        = C%filename_climate_snapshot_unif_dT_GRL
        climate%snapshot_deltaT%precip_CC_correction     = C%precip_CC_correction_GRL
        climate%snapshot_deltaT%snapshot%lapse_rate_temp = C%lapse_rate_temp_GRL
        IF (C%choice_SMB_model_GRL == 'IMAU-ITM') THEN
           climate%snapshot%has_insolation = .TRUE.
        END IF
      ELSEIF (region_name == 'ANT') THEN
        filename_climate_snapshot                        = C%filename_climate_snapshot_unif_dT_ANT
        climate%snapshot_deltaT%precip_CC_correction.    = C%precip_CC_correction_ANT
        climate%snapshot_deltaT%snapshot%lapse_rate_temp = C%lapse_rate_temp_ANT
        IF (C%choice_SMB_model_ANT == 'IMAU-ITM') THEN
           climate%snapshot%has_insolation = .TRUE.
        END IF
      ELSE
        CALL crash('unknown region_name "' // region_name // '"')
      END IF

      CALL read_field_from_file_2D(         filename_climate_snapshot, 'Hs'    , mesh, C%output_dir, climate%snapshot%Hs)
      CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'T2m'   , mesh, C%output_dir, climate%T2m)
      CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'Precip', mesh, C%output_dir, climate%Precip)


      ! Adding deltaT to the temperature field (uniform in space and time)
      do vi = mesh%vi1, mesh%vi2
        do m = 1, 12
            climate%T2m( vi, m) = climate%snapshot_deltaT%snapshot%T2m + climate%snapshot_deltaT%deltaT
        end do
      end do

      ! apply corrections (increase in Precip due to deltaT, plus downscaling correction)
      call apply_precipitation_CC_correction(mesh, climate)
      call apply_geometry_downscaling_corrections( mesh, ice, climate)
      

      ! Initialises the insolation (if needed)
      IF (climate%snapshot%has_insolation) THEN
        IF (C%choice_insolation_forcing == 'none') THEN
          CALL crash('Chosen climate or SMB model cannot be used with choice_insolation_forcing = "none"!')
        ELSE
          CALL initialise_insolation_forcing( climate%snapshot, mesh)
          IF (C%start_time_of_run < 0._dp) THEN
            timeframe_init_insolation = C%start_time_of_run
          ELSE
            timeframe_init_insolation = 0._dp
          END IF
          CALL get_insolation_at_time( mesh, timeframe_init_insolation, climate%snapshot)
        END IF
      END IF

    ELSE
      CALL crash('unknown choice_climate_model_realistic "' // TRIM( C%choice_climate_model_realistic) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_snapshot_plus_uniform_deltaT

  SUBROUTINE apply_precipitation_CC_correction(mesh, climate)
  ! Applies a simple Clausius-Clapeyron correction to temperature based on the prescribed deltaT

    IMPLICIT NONE

    TYPE(type_mesh),                       INTENT(IN)    :: mesh
    TYPE(type_climate_model),              INTENT(INOUT) :: climate

    ! Local Variables
    CHARACTER(LEN=256), PARAMETER                        :: routine_name = 'apply_precipitation_CC_correction'
    INTEGER                                              :: vi, m
    
    ! Add routine to path
    CALL init_routine( routine_name)

    do vi = mesh%vi1, mesh%vi2
        do m = 1, 12
            ! Precip(\Delta T) = Precip(PD) \times 1.068^{\Delta T}
            climate%Precip( vi, m) = climate%snapshot_plus_uniform_deltaT%snapshot%Precip( vi, m) * climate%snapshot_plus_uniform_deltaT%precip_CC_correction**climate%snapshot_plus_uniform_deltaT%deltaT
        end do
    end do

    ! Finalise routine path
    CALL finalise_routine( routine_name)
  END SUBROUTINE apply_precipitation_CC_correction

  SUBROUTINE apply_geometry_downscaling_corrections( mesh, ice, climate)
    ! Applies the lapse rate corrections for temperature and precipitation
    ! to correct for the mismatch between T and P at the forcing's ice surface elevation and the model's ice surface elevation

    IMPLICIT NONE

    TYPE(type_mesh),                       INTENT(IN)    :: mesh
    TYPE(type_ice_model),                  INTENT(IN)    :: ice
    TYPE(type_climate_model),              INTENT(INOUT) :: climate

    ! Local Variables
    CHARACTER(LEN=256), PARAMETER                        :: routine_name = 'apply_geometry_downscaling_corrections'
    INTEGER                                              :: vi, m
    REAL(dp)                                             :: deltaH, deltaT, deltaP
    REAL(dp), DIMENSION(:,:), ALLOCATABLE                :: T_inv, T_inv_ref

    ! Add routine to path
    CALL init_routine( routine_name)

    allocate( T_inv     (mesh%vi1:mesh%vi2, 12))
    allocate( T_inv_ref (mesh%vi1:mesh%vi2, 12))


    do vi = mesh%vi1, mesh%vi2
    ! we only apply corrections where it is not open ocean
    if (ice%mask_icefree_ocean( vi) .eqv. .FALSE.) then
        deltaT  = (ice%Hs( vi) - climate%snapshot_deltaT%snapshot%Hs( vi)) * (-1._dp * abs(climate%snapshot_deltaT%snapshot%lapse_rate_temp))
        do m = 1, 12
        ! Do corrections - based on Eq. 11 of Albrecht et al. (2020; TC) for PISM
        climate%T2m( vi, m)    = climate%T2m( vi, m)    + deltaT

        ! Calculate inversion-layer temperatures
        T_inv_ref( vi, m) = 88.9_dp + 0.67_dp *  climate%T2m( vi, m)
        T_inv(     vi, m) = 88.9_dp + 0.67_dp * (climate%T2m( vi, m) - climate%snapshot_deltaT%snapshot%lapse_rate_temp * (ice%Hs( vi) - climate%snapshot_deltaT%snapshot%Hs( vi)))
        ! Correct precipitation based on a simple Clausius-Clapeyron method (Jouzel & Merlivat, 1984; Huybrechts, 2002)
        ! Same as implemented in IMAU-ICE
        climate%Precip( vi, m) = climate%Precip( vi, m) * (T_inv_ref( vi, m) / T_inv( vi, m))**2 * EXP(22.47_dp * (T0 / T_inv_ref( vi, m) - T0 / T_inv( vi, m)))

        end do ! m
    end if
    end do ! vi

    deallocate(T_inv)
    deallocate(T_inv_ref)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE apply_geometry_downscaling_corrections

  


END MODULE climate_snapshot_plus_uniform_deltaT
