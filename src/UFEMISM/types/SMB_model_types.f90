module SMB_model_types

  use precisions, only: dp
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_SMB_model_snapshot_plus_anomalies

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

end module SMB_model_types