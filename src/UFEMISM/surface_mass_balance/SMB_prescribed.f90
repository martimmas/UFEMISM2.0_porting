MODULE SMB_prescribed

  ! Prescribed SMB forcing

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, warning, happy, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model
  use netcdf_io_main

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_SMB_model_prescribed( mesh, ice, SMB, region_name, time)
    ! Calculate the surface mass balance
    !
    ! Prescribed SMB forcing

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_SMB_model_prescribed'
    CHARACTER(LEN=256)                                    :: choice_SMB_prescribed

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine the type of prescribed SMB forcing for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_prescribed  = C%choice_SMB_prescribed_NAM
      CASE ('EAS')
        choice_SMB_prescribed  = C%choice_SMB_prescribed_EAS
      CASE ('GRL')
        choice_SMB_prescribed  = C%choice_SMB_prescribed_GRL
      CASE ('ANT')
        choice_SMB_prescribed  = C%choice_SMB_prescribed_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END SELECT

    ! Run the chosen type of prescribed SMB forcing
    SELECT CASE (choice_SMB_prescribed)
      CASE ('SMB_no_time')
        ! SMB only, no time
        CALL run_SMB_model_prescribed_notime( mesh, SMB)
      CASE DEFAULT
        CALL crash('unknown choice_SMB_prescribed "' // TRIM( choice_SMB_prescribed) // '"!')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model_prescribed

  SUBROUTINE initialise_SMB_model_prescribed( mesh, SMB, region_name)
    ! Initialise the SMB model
    !
    ! Prescribed SMB forcing

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_SMB_model_prescribed'
    CHARACTER(LEN=256)                                    :: choice_SMB_prescribed

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine the type of prescribed SMB forcing for this region
    SELECT CASE (region_name)
      CASE ('NAM')
        choice_SMB_prescribed  = C%choice_SMB_prescribed_NAM
      CASE ('EAS')
        choice_SMB_prescribed  = C%choice_SMB_prescribed_EAS
      CASE ('GRL')
        choice_SMB_prescribed  = C%choice_SMB_prescribed_GRL
      CASE ('ANT')
        choice_SMB_prescribed  = C%choice_SMB_prescribed_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END SELECT

    ! Initialised the chosen type of prescribed SMB forcing
    SELECT CASE (choice_SMB_prescribed)
      CASE ('SMB_no_time')
        CALL initialise_SMB_model_prescribed_notime( mesh, SMB, region_name)
      CASE DEFAULT
        CALL crash('unknown choice_SMB_prescribed "' // TRIM( choice_SMB_prescribed) // '"!')
    END SELECT

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SMB_model_prescribed

  ! == SMB only, no time
  ! ====================

  SUBROUTINE run_SMB_model_prescribed_notime( mesh, SMB)
    ! Calculate the surface mass balance
    !
    ! Prescribed SMB forcing
    !
    ! SMB only, no time

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: SMB

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_SMB_model_prescribed_notime'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! No need to do anything, as the SMB was already read during initialisation

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_SMB_model_prescribed_notime

  SUBROUTINE initialise_SMB_model_prescribed_notime( mesh, SMB, region_name)
    ! Initialise the SMB model
    !
    ! Prescribe SMB from a file without a time dimension
    !
    ! SMB only, no time

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: SMB
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_SMB_model_prescribed_notime'
    CHARACTER(LEN=256)                                    :: filename_SMB_prescribed
    REAL(dp)                                              :: timeframe_SMB_prescribed

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Determine filename for this model region
    SELECT CASE (region_name)
      CASE ('NAM')
        filename_SMB_prescribed  = C%filename_SMB_prescribed_NAM
        timeframe_SMB_prescribed = C%timeframe_SMB_prescribed_NAM
      CASE ('EAS')
        filename_SMB_prescribed  = C%filename_SMB_prescribed_EAS
        timeframe_SMB_prescribed = C%timeframe_SMB_prescribed_EAS
      CASE ('GRL')
        filename_SMB_prescribed  = C%filename_SMB_prescribed_GRL
        timeframe_SMB_prescribed = C%timeframe_SMB_prescribed_GRL
      CASE ('ANT')
        filename_SMB_prescribed  = C%filename_SMB_prescribed_ANT
        timeframe_SMB_prescribed = C%timeframe_SMB_prescribed_ANT
      CASE DEFAULT
        CALL crash('unknown region_name "' // TRIM( region_name) // '"!')
    END SELECT

    ! Print to terminal
    IF (par%primary)  WRITE(*,"(A)") '   Initialising SMB from file "' // colour_string( TRIM( filename_SMB_prescribed),'light blue') // '"...'

    ! Read SMB from file
    IF (timeframe_SMB_prescribed == 1E9_dp) THEN
      ! Assume the file has no time dimension
      CALL read_field_from_file_2D( filename_SMB_prescribed, 'SMB||surface_mass_balance||', mesh, C%output_dir, SMB)
    ELSE
      ! Assume the file has a time dimension, and read the specified timeframe
      CALL read_field_from_file_2D( filename_SMB_prescribed, 'SMB||surface_mass_balance||', mesh, C%output_dir, SMB, time_to_read = timeframe_SMB_prescribed)
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_SMB_model_prescribed_notime

END MODULE SMB_prescribed
