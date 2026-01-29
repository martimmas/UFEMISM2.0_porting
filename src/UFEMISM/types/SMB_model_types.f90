MODULE SMB_model_types

  ! The different data types used in the SMB modules

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  use SMB_idealised, only: type_SMB_model_idealised
  use SMB_prescribed, only: type_SMB_model_prescribed
  use SMB_reconstructed, only: type_SMB_model_reconstructed
  use SMB_snapshot_plus_anomalies, only: type_SMB_model_snp_p_anml
  use SMB_IMAU_ITM, only: type_SMB_model_IMAU_ITM

  IMPLICIT NONE

  TYPE type_SMB_model
    ! The SMB model data structure.

    ! Main data fields
    REAL(dp), DIMENSION(:    ), ALLOCATABLE      :: SMB                       ! Yearly  SMB (m)

    ! Sub-models
    REAL(dp), DIMENSION(:    ), ALLOCATABLE      :: SMB_correction            ! [m.i.e./yr] Surface mass balance
    type(type_SMB_model_idealised)               :: idealised
    type(type_SMB_model_prescribed)              :: prescribed
    type(type_SMB_model_reconstructed)           :: reconstructed
    TYPE(type_SMB_model_IMAU_ITM)                :: IMAUITM
    type(type_SMB_model_snp_p_anml)              :: snapshot_plus_anomalies

    ! Timestepping
    REAL(dp)                                     :: t_next

    ! Metadata
    CHARACTER(LEN=256)                           :: restart_filename          ! Name for generated restart file

  END TYPE type_SMB_model

CONTAINS

END MODULE SMB_model_types