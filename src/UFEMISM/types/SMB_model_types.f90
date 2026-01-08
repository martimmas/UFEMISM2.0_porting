module SMB_model_types

  use precisions, only: dp

  implicit none

  private

  public :: type_SMB_model, type_SMB_model_IMAU_ITM, type_SMB_model_snapshot_plus_anomalies

  type type_SMB_model_IMAU_ITM
    !< The IMAU-ITM SMB model data structure

    ! Main data fields
    real(dp), dimension(:  ), allocatable :: AlbedoSurf              ! Surface albedo underneath the snow layer (water, rock or ice)
    real(dp), dimension(:  ), allocatable :: MeltPreviousYear        ! [m.w.e.] total melt in the previous year
    real(dp), dimension(:,:), allocatable :: FirnDepth               ! [m] depth of the firn layer
    real(dp), dimension(:,:), allocatable :: Rainfall                ! Monthly rainfall (m)
    real(dp), dimension(:,:), allocatable :: Snowfall                ! Monthly snowfall (m)
    real(dp), dimension(:,:), allocatable :: AddedFirn               ! Monthly added firn (m)
    real(dp), dimension(:,:), allocatable :: Melt                    ! Monthly melt (m)
    real(dp), dimension(:,:), allocatable :: Refreezing              ! Monthly refreezing (m)
    real(dp), dimension(:  ), allocatable :: Refreezing_year         ! Yearly  refreezing (m)
    real(dp), dimension(:,:), allocatable :: Runoff                  ! Monthly runoff (m)
    real(dp), dimension(:,:), allocatable :: Albedo                  ! Monthly albedo
    real(dp), dimension(:  ), allocatable :: Albedo_year             ! Yearly albedo
    real(dp), dimension(:,:), allocatable :: SMB_monthly             ! [m] Monthly SMB

    ! Tuning parameters for the IMAU-ITM SMB model (different for each region, set from config)
    real(dp)  :: C_abl_constant
    real(dp)  :: C_abl_Ts
    real(dp)  :: C_abl_Q
    real(dp)  :: C_refr

    ! Ideally these parameters should not be region-dependent?
    real(dp)  :: albedo_water
    real(dp)  :: albedo_soil
    real(dp)  :: albedo_ice
    real(dp)  :: albedo_snow

  end type type_SMB_model_IMAU_ITM

  type type_SMB_model_snapshot_plus_anomalies

    ! Baseline climate
    real(dp), dimension(:,:), allocatable :: T2m_baseline
    real(dp), dimension(:  ), allocatable :: SMB_baseline

    ! Two anomaly timeframes enveloping the current model time
    real(dp)                              :: anomaly_t0
    real(dp), dimension(:  ), allocatable :: T2m_anomaly_0
    real(dp), dimension(:  ), allocatable :: SMB_anomaly_0

    real(dp)                              :: anomaly_t1
    real(dp), dimension(:  ), allocatable :: T2m_anomaly_1
    real(dp), dimension(:  ), allocatable :: SMB_anomaly_1

    ! Time-weighted anomaly
    real(dp), dimension(:  ), allocatable :: T2m_anomaly
    real(dp), dimension(:  ), allocatable :: SMB_anomaly

    ! Applied climate
    real(dp), dimension(:,:), allocatable :: T2m    ! = baseline + anomaly
    real(dp), dimension(:  ), allocatable :: SMB

  end type type_SMB_model_snapshot_plus_anomalies

  type type_SMB_model
    ! The SMB model data structure.

    ! Main data fields
    real(dp), dimension(:    ), allocatable      :: SMB                       ! Yearly  SMB (m)

    ! Sub-models
    real(dp), dimension(:    ), allocatable      :: SMB_correction            ! [m.i.e./yr] Surface mass balance
    TYPE(type_SMB_model_IMAU_ITM)                :: IMAUITM
    type(type_SMB_model_snapshot_plus_anomalies) :: snapshot_plus_anomalies

    ! Timestepping
    real(dp)                                     :: t_next

    ! Metadata
    character(:), allocatable                    :: restart_filename          ! Name for generated restart file

  end type type_SMB_model

end module SMB_model_types