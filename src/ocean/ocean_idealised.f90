MODULE ocean_idealised

  ! Idealised ocean models

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE ocean_model_types                                      , ONLY: type_ocean_model

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_ocean_model_idealised( mesh, ice, ocean)
    ! Calculate the ocean
    !
    ! Use an idealised ocean scheme

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_ocean_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen idealised ocean model
    IF (C%choice_ocean_model_idealised == 'ISOMIP_WARM' .OR. &
        C%choice_ocean_model_idealised == 'ISOMIP_COLD') THEN
      ! No need to do anything 
    ELSE
      CALL crash('unknown choice_ocean_model_idealised "' // TRIM( C%choice_ocean_model_idealised) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_ocean_model_idealised

  SUBROUTINE initialise_ocean_model_idealised( mesh, ocean)
    ! Initialise the ocean model
    !
    ! Use an idealised ocean scheme

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ocean_model),                 INTENT(INOUT) :: ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_ocean_model_idealised'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '     Initialising idealised ocean model "' // &
      colour_string( TRIM( C%choice_ocean_model_idealised),'light blue') // '"...'

    ! Run the chosen idealised ocean model
    IF (C%choice_ocean_model_idealised == 'ISOMIP_WARM' .OR. &
        C%choice_ocean_model_idealised == 'ISOMIP_COLD') THEN
      CALL initialise_ocean_model_idealised_ISOMIP( mesh, ocean)
    ELSE
      CALL crash('unknown choice_ocean_model_idealised "' // TRIM( C%choice_ocean_model_idealised) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ocean_model_idealised

  ! == ISOMIP ==
  ! ============

  SUBROUTINE initialise_ocean_model_idealised_ISOMIP( mesh, ocean)
    ! 

    IMPLICIT NONE

    TYPE(type_mesh),                      INTENT(IN)    :: mesh
    TYPE(type_ocean_model),               INTENT(INOUT) :: ocean

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                       :: routine_name = 'initialise_ocean_model_idealised_ISOMIP'
    INTEGER                                             :: vi
    INTEGER                                             :: k
    REAL(dp), PARAMETER                                 :: z1 = 720_dp  ! [m] Reference depth
    REAL(dp), PARAMETER                                 :: T0 = -1.9_dp ! [degC] Surface temperature
    REAL(dp)                                            :: T1           ! [degC] Reference temperature
    REAL(dp), PARAMETER                                 :: S0 = 33.8_dp ! [PSU]  Surface salinity
    REAL(dp)                                            :: S1           ! [PSU]  Reference salinity

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Define scenario-dependent parameters
    IF (C%choice_ocean_model_idealised == 'ISOMIP_WARM') THEN
      T1 = 1.0_dp
      S1 = 34.7_dp
    ELSE IF (C%choice_ocean_model_idealised == 'ISOMIP_COLD') THEN
      T1 = -1.9_dp
      S1 = 34.55_dp
    END IF
    
    DO vi = mesh%vi1, mesh%vi2
      DO k = 1, C%nz_ocean
        ocean%T( vi, k) = T0 + (T1-T0)*C%z_ocean( k)/z1      
        ocean%S( vi, k) = S0 + (S1-S0)*C%z_ocean( k)/z1      
      END DO
    END DO

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_ocean_model_idealised_ISOMIP

END MODULE ocean_idealised
