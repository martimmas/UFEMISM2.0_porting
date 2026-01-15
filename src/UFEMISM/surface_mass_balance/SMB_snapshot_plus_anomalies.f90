module SMB_snapshot_plus_anomalies

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mesh_types, only: type_mesh
  use climate_model_types, only: type_climate_model
  use parameters, only: NaN
  use model_configuration, only: C
  use netcdf_io_main
  use mpi_f08, only: MPI_WIN, MPI_BCAST, MPI_DOUBLE_PRECISION, MPI_COMM_WORLD
  use allocate_dist_shared_mod, only: allocate_dist_shared

  implicit none

  private

  public :: &
    type_SMB_model_snapshot_plus_anomalies, &
    run_climate_model_SMB_snapshot_plus_anomalies, &
    initialise_SMB_model_snapshot_plus_anomalies, &
    run_SMB_model_snapshot_plus_anomalies

  type type_SMB_model_snapshot_plus_anomalies

    ! Baseline climate
    real(dp), dimension(:,:), contiguous, pointer :: T2m_baseline
    real(dp), dimension(:  ), contiguous, pointer :: SMB_baseline
    type(MPI_WIN) :: wT2m_baseline, wSMB_baseline

    ! Two anomaly timeframes enveloping the current model time
    real(dp)                              :: anomaly_t0
    real(dp), dimension(:  ), contiguous, pointer :: T2m_anomaly_0
    real(dp), dimension(:  ), contiguous, pointer :: SMB_anomaly_0
    type(MPI_WIN) :: wT2m_anomaly_0, wSMB_anomaly_0

    real(dp)                              :: anomaly_t1
    real(dp), dimension(:  ), contiguous, pointer :: T2m_anomaly_1
    real(dp), dimension(:  ), contiguous, pointer :: SMB_anomaly_1
    type(MPI_WIN) :: wT2m_anomaly_1, wSMB_anomaly_1

    ! Time-weighted anomaly
    real(dp), dimension(:  ), contiguous, pointer :: T2m_anomaly
    real(dp), dimension(:  ), contiguous, pointer :: SMB_anomaly
    type(MPI_WIN) :: wT2m_anomaly, wSMB_anomaly

    ! Applied climate
    real(dp), dimension(:,:), contiguous, pointer :: T2m    ! = baseline + anomaly
    real(dp), dimension(:  ), contiguous, pointer :: SMB
    type(MPI_WIN) :: wT2m, wSMB

  end type type_SMB_model_snapshot_plus_anomalies

contains

  subroutine run_climate_model_SMB_snapshot_plus_anomalies( mesh, climate, snapshot_plus_anomalies, time)

    ! In/output variables:
    type(type_mesh),                              intent(in   ) :: mesh
    type(type_climate_model),                     intent(inout) :: climate
    type(type_SMB_model_snapshot_plus_anomalies), intent(inout) :: snapshot_plus_anomalies
    real(dp),                                     intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_climate_model_SMB_snapshot_plus_anomalies'
    real(dp)                       :: w0, w1
    integer                        :: m

    ! Add routine to path
    call init_routine( routine_name)

    ! If the current model time falls outside the enveloping window
    ! of the two timeframes that have been read, update them
    if (time < snapshot_plus_anomalies%anomaly_t0 .or. &
        time > snapshot_plus_anomalies%anomaly_t1) then
      call update_timeframes( mesh, snapshot_plus_anomalies, time)
    end if

    ! Interpolate between the two timeframes to find the applied anomaly
    w0 = (snapshot_plus_anomalies%anomaly_t1 - time) / &
         (snapshot_plus_anomalies%anomaly_t1 - snapshot_plus_anomalies%anomaly_t0)
    w1 = 1._dp - w0

    snapshot_plus_anomalies%T2m_anomaly = &
      w0 * snapshot_plus_anomalies%T2m_anomaly_0 + &
      w1 * snapshot_plus_anomalies%T2m_anomaly_1

    ! Add anomaly to snapshot to find the applied temperature
    do m = 1, 12
      snapshot_plus_anomalies%T2m( mesh%vi1:mesh%vi2,m) = &
        snapshot_plus_anomalies%T2m_baseline( mesh%vi1:mesh%vi2,m) + &
        snapshot_plus_anomalies%T2m_anomaly ( mesh%vi1:mesh%vi2  )
    end do

    ! Copy to climate model
    climate%T2m( mesh%vi1:mesh%vi2,:) = snapshot_plus_anomalies%T2m( mesh%vi1:mesh%vi2,:)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_climate_model_SMB_snapshot_plus_anomalies

  subroutine run_SMB_model_snapshot_plus_anomalies( mesh, snapshot_plus_anomalies, time)

    ! In/output variables:
    type(type_mesh),                              intent(in   ) :: mesh
    type(type_SMB_model_snapshot_plus_anomalies), intent(inout) :: snapshot_plus_anomalies
    real(dp),                                     intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'run_SMB_model_snapshot_plus_anomalies'
    real(dp)                       :: w0, w1

    ! Add routine to path
    call init_routine( routine_name)

    ! If the current model time falls outside the enveloping window
    ! of the two timeframes that have been read, update them
    if (time < snapshot_plus_anomalies%anomaly_t0 .or. &
        time > snapshot_plus_anomalies%anomaly_t1) then
      call update_timeframes( mesh, snapshot_plus_anomalies, time)
    end if

    ! Interpolate between the two timeframes to find the applied anomaly
    w0 = (snapshot_plus_anomalies%anomaly_t1 - time) / &
         (snapshot_plus_anomalies%anomaly_t1 - snapshot_plus_anomalies%anomaly_t0)
    w1 = 1._dp - w0

    snapshot_plus_anomalies%SMB_anomaly( mesh%vi1: mesh%vi2) = &
      w0 * snapshot_plus_anomalies%SMB_anomaly_0( mesh%vi1: mesh%vi2) + &
      w1 * snapshot_plus_anomalies%SMB_anomaly_1( mesh%vi1: mesh%vi2)

    ! Add anomaly to snapshot to find the applied SMB
    snapshot_plus_anomalies%SMB( mesh%vi1: mesh%vi2) = &
      snapshot_plus_anomalies%SMB_baseline( mesh%vi1: mesh%vi2) + &
      snapshot_plus_anomalies%SMB_anomaly( mesh%vi1: mesh%vi2)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_SMB_model_snapshot_plus_anomalies

  subroutine initialise_SMB_model_snapshot_plus_anomalies( mesh, snapshot_plus_anomalies)

    ! In/output variables:
    type(type_mesh),                              intent(in   ) :: mesh
    type(type_SMB_model_snapshot_plus_anomalies), intent(inout) :: snapshot_plus_anomalies

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_SMB_model_snapshot_plus_anomalies'

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate memory

    ! Baseline climate
    call allocate_dist_shared( snapshot_plus_anomalies%T2m_baseline, snapshot_plus_anomalies%wT2m_baseline, mesh%pai_V%n_nih, 12)
    call allocate_dist_shared( snapshot_plus_anomalies%SMB_baseline, snapshot_plus_anomalies%wSMB_baseline, mesh%pai_V%n_nih    )
    snapshot_plus_anomalies%T2m_baseline( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:12) => snapshot_plus_anomalies%T2m_baseline
    snapshot_plus_anomalies%SMB_baseline( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih      ) => snapshot_plus_anomalies%SMB_baseline

    ! Two anomaly snapshots enveloping the current model time
    call allocate_dist_shared( snapshot_plus_anomalies%T2m_anomaly_0, snapshot_plus_anomalies%wT2m_anomaly_0, mesh%pai_V%n_nih)
    call allocate_dist_shared( snapshot_plus_anomalies%SMB_anomaly_0, snapshot_plus_anomalies%wSMB_anomaly_0, mesh%pai_V%n_nih)
    call allocate_dist_shared( snapshot_plus_anomalies%T2m_anomaly_1, snapshot_plus_anomalies%wT2m_anomaly_1, mesh%pai_V%n_nih)
    call allocate_dist_shared( snapshot_plus_anomalies%SMB_anomaly_1, snapshot_plus_anomalies%wSMB_anomaly_1, mesh%pai_V%n_nih)
    snapshot_plus_anomalies%T2m_anomaly_0( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => snapshot_plus_anomalies%T2m_anomaly_0
    snapshot_plus_anomalies%SMB_anomaly_0( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => snapshot_plus_anomalies%SMB_anomaly_0
    snapshot_plus_anomalies%T2m_anomaly_1( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => snapshot_plus_anomalies%T2m_anomaly_1
    snapshot_plus_anomalies%SMB_anomaly_1( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => snapshot_plus_anomalies%SMB_anomaly_1

    ! Time-weighted anomaly
    call allocate_dist_shared( snapshot_plus_anomalies%T2m_anomaly, snapshot_plus_anomalies%wT2m_anomaly, mesh%pai_V%n_nih)
    call allocate_dist_shared( snapshot_plus_anomalies%SMB_anomaly, snapshot_plus_anomalies%wSMB_anomaly, mesh%pai_V%n_nih)
    snapshot_plus_anomalies%T2m_anomaly( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => snapshot_plus_anomalies%T2m_anomaly
    snapshot_plus_anomalies%SMB_anomaly( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => snapshot_plus_anomalies%SMB_anomaly

    ! Applied climate
    call allocate_dist_shared( snapshot_plus_anomalies%T2m, snapshot_plus_anomalies%wT2m, mesh%pai_V%n_nih, 12)
    call allocate_dist_shared( snapshot_plus_anomalies%SMB, snapshot_plus_anomalies%wSMB, mesh%pai_V%n_nih    )
    snapshot_plus_anomalies%T2m( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:12) => snapshot_plus_anomalies%T2m
    snapshot_plus_anomalies%SMB( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih      ) => snapshot_plus_anomalies%SMB

    ! Read baseline snapshot
    call read_field_from_file_2D_monthly( C%SMB_snp_p_anml_filename_snapshot_T2m, 'T2m', &
      mesh, C%output_dir, snapshot_plus_anomalies%T2m_baseline)
    call read_field_from_file_2D( C%SMB_snp_p_anml_filename_snapshot_SMB, 'SMB', &
      mesh, C%output_dir, snapshot_plus_anomalies%SMB_baseline)

    ! Initialise anomaly timeframes
    snapshot_plus_anomalies%anomaly_t0 = C%start_time_of_run - 200._dp
    snapshot_plus_anomalies%anomaly_t1 = C%start_time_of_run - 100._dp
    call update_timeframes( mesh, snapshot_plus_anomalies, C%start_time_of_run)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_SMB_model_snapshot_plus_anomalies

  subroutine update_timeframes( mesh, snapshot_plus_anomalies, time)

    ! In/output variables:
    type(type_mesh),                              intent(in   ) :: mesh
    type(type_SMB_model_snapshot_plus_anomalies), intent(inout) :: snapshot_plus_anomalies
    real(dp),                                     intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'update_timeframes'
    character(len=1024)                 :: filename
    integer                             :: ncid, id_dim_time, nt, id_var_time, ierr
    real(dp), dimension(:), allocatable :: time_from_file
    integer                             :: ti0, ti1

    ! Add routine to path
    call init_routine( routine_name)

    filename = trim( C%SMB_snp_p_anml_filename_anomalies)

    ! Read time variable from the file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    call check_time( filename, ncid)
    call inquire_dim_multopt( filename, ncid, field_name_options_time, id_dim_time, dim_length = nt)
    call inquire_var_multopt( filename, ncid, field_name_options_time, id_var_time)
    allocate( time_from_file( nt))
    call read_var_primary( filename, ncid, id_var_time, time_from_file)
    call MPI_BCAST( time_from_file(:), nt, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call close_netcdf_file( ncid)

    ! Find the two timeframes
    if (time < time_from_file( 1)) then
      if (par%primary) call warning('model time before start of anomaly file; using first timeframe')
      ti0 = 1
      ti1 = 2
    elseif (time > time_from_file( size( time_from_file,1))) then
      if (par%primary) call warning('model time beyond end of anomaly file; using last timeframe')
      ti0 = size( time_from_file,1) - 1
      ti1 = size( time_from_file,1)
    else
      ti0 = 1
      ti1 = 2
      do while (time_from_file( ti1) < time .and. ti1 < size( time_from_file,1))
        ti0 = ti1
        ti1 = ti1 + 1
      end do
    end if

    snapshot_plus_anomalies%anomaly_t0 = time_from_file( ti0)
    snapshot_plus_anomalies%anomaly_t1 = time_from_file( ti1)

    ! Read the two timeframes
    call read_field_from_file_2D( filename, 'T2m_anomaly', &
      mesh, C%output_dir, snapshot_plus_anomalies%T2m_anomaly_0, &
      time_to_read = snapshot_plus_anomalies%anomaly_t0)
    call read_field_from_file_2D( filename, 'T2m_anomaly', &
      mesh, C%output_dir, snapshot_plus_anomalies%T2m_anomaly_1, &
      time_to_read = snapshot_plus_anomalies%anomaly_t1)
    call read_field_from_file_2D( filename, 'SMB_anomaly', &
      mesh, C%output_dir, snapshot_plus_anomalies%SMB_anomaly_0, &
      time_to_read = snapshot_plus_anomalies%anomaly_t0)
    call read_field_from_file_2D( filename, 'SMB_anomaly', &
      mesh, C%output_dir, snapshot_plus_anomalies%SMB_anomaly_1, &
      time_to_read = snapshot_plus_anomalies%anomaly_t1)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_timeframes

end module SMB_snapshot_plus_anomalies
