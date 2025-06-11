module ocean_realistic

  ! Realistic ocean models

! ===== Preamble =====
! ====================

  use precisions                                             , only: dp
  use mpi_basic                                              , only: par, sync
  use control_resources_and_error_messaging                  , only: crash, init_routine, finalise_routine, colour_string
  use model_configuration                                    , only: C
  use parameters
  use mesh_types                                             , only: type_mesh
  use ice_model_types                                        , only: type_ice_model
  use ocean_model_types                                      , only: type_ocean_model
  use netcdf_io_main
  use ocean_extrapolation                                    , only: extrapolate_ocean_forcing

  implicit none

contains

! ===== Main routines =====
! =========================

  subroutine run_ocean_model_realistic( mesh, ice, ocean, time)
    ! Calculate the ocean
    !
    ! Use an realistic ocean scheme

    ! In/output variables:
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(in)    :: ice
    type(type_ocean_model),                 intent(inout) :: ocean
    real(dp),                               intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'run_ocean_model_realistic'

    ! Add routine to path
    call init_routine( routine_name)

    ! Run the chosen realistic ocean model
    select case (C%choice_ocean_model_realistic)
      case ('snapshot')

        ! Apply extrapolation method if required
        select case (C%choice_ocean_extrapolation_method)
          case('initialisation')
            ! nothing to do here 
          case default
            call crash('unknown choice_ocean_extrapolation_method "' // trim( C%choice_ocean_extrapolation_method) // '"')
        end select

      case ('transient')
        call run_ocean_model_transient(mesh, ocean, time)

      case default
        call crash('unknown choice_ocean_model_realistic "' // trim( C%choice_ocean_model_realistic) // '"')
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_ocean_model_realistic

  subroutine initialise_ocean_model_realistic( mesh, ice, ocean, region_name, start_time_of_run)
    ! Initialise the ocean model
    !
    ! Use an realistic ocean scheme

    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(in)    :: ice
    type(type_ocean_model),                 intent(inout) :: ocean
    character(len=3),                       intent(in)    :: region_name
    real(dp),                               intent(in)    :: start_time_of_run

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'initialise_ocean_model_realistic'
    character(len=256)                                    :: filename_ocean_snapshot

    ! Add routine to path
    call init_routine( routine_name)

    ! Print to terminal
    if (par%primary)  write(*,"(A)") '     Initialising realistic ocean model "' // &
      colour_string( trim( C%choice_ocean_model_realistic),'light blue') // '"...'

    ! Run the chosen realistic ocean model
    select case (C%choice_ocean_model_realistic)
      case ('snapshot')
        ! Read single-time data from external file

        select case (region_name)
          case ('NAM')
            filename_ocean_snapshot = C%filename_ocean_snapshot_NAM
          case ('EAS')
            filename_ocean_snapshot = C%filename_ocean_snapshot_EAS
          case ('GRL')
            filename_ocean_snapshot = C%filename_ocean_snapshot_GRL
          case ('ANT')
            filename_ocean_snapshot = C%filename_ocean_snapshot_ANT
          case default
            call crash('unknown region_name "' // region_name // '"')
        end select

        ! Fill in  main variables
        call read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_T_ocean, mesh, ocean%T)
        call read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_S_ocean, mesh, ocean%S)

        ! Apply extrapolation method if required
        select case (C%choice_ocean_extrapolation_method)
          case('initialisation')
            call extrapolate_ocean_forcing( mesh, ice, ocean%T)
            call extrapolate_ocean_forcing( mesh, ice, ocean%S) 
          case default
            call crash('unknown choice_ocean_extrapolation_method "' // trim( C%choice_ocean_extrapolation_method) // '"')
        end select

      case ('transient')  
        select case (C%choice_ocean_model_transient)
          case ('deltaT')

            ! We need the snapshot and the dT to apply to it
            select case (region_name)
              case ('NAM')
                filename_ocean_snapshot = C%filename_ocean_snapshot_NAM
                filename_ocean_dT       = C%filename_ocean_dT_NAM
              case ('EAS')
                filename_ocean_snapshot = C%filename_ocean_snapshot_EAS
                filename_ocean_dT       = C%filename_ocean_dT_EAS
              case ('GRL')
                filename_ocean_snapshot = C%filename_ocean_snapshot_GRL
                filename_ocean_dT       = C%filename_ocean_dT_GRL
              case ('ANT')
                filename_ocean_snapshot = C%filename_ocean_snapshot_ANT
                filename_ocean_dT       = C%filename_ocean_dT_ANT
              case default
                call crash('unknown region_name "' // region_name // '"')
            end select

            ! Fill in  main variables
            call read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_T_ocean,  mesh, ocean%transient%T0)
            call read_field_from_file_3D_ocean( filename_ocean_snapshot, field_name_options_S_ocean,  mesh, ocean%transient%S0)

            ! Allocating timeframe variables; the series itself is allocated in the read function below
            allocate(ocean%transient%dT_t0)
            allocate(ocean%transient%dT_t1)
            allocate(ocean%transient%dT_at_t0)
            allocate(ocean%transient%dT_at_t1)
            allocate( ocean%transient%T0( mesh%vi1:mesh%vi2,C%nz_ocean))
            allocate( ocean%transient%S0( mesh%vi1:mesh%vi2,C%nz_ocean))
            ocean%transient%T0 = 0._dp
            ocean%transient%S0 = 0._dp
            
            call read_field_from_series_file(   filename_ocean_dT,       field_name_options_dT_ocean, ocean%transient%dT_series, ocean%transient%dT_series_time)
            call update_dT_timeframes_from_curve(ocean, start_time_of_run)

            ! Apply extrapolation method if required
            select case (C%choice_ocean_extrapolation_method)
              case('initialisation')
                call extrapolate_ocean_forcing( mesh, ice, ocean%transient%T0)
                call extrapolate_ocean_forcing( mesh, ice, ocean%transient%S0) 
              case default
                call crash('unknown choice_ocean_extrapolation_method "' // trim( C%choice_ocean_extrapolation_method) // '"')
            end select

          case default
            call crash('unknown choice_ocean_model_transient "' // trim( C%choice_ocean_model_transient) // '"')
      end select
        

      case default
        call crash('unknown choice_ocean_model_realistic "' // trim( C%choice_ocean_model_realistic) // '"')
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_ocean_model_realistic

  SUBROUTINE update_dT_timeframes_from_curve(ocean, time)
    ! Update the ocean dT timeframes so we can interpolate between two points in the ocean dT curve

    IMPLICIT NONE

    TYPE(type_ocean_model),         INTENT(INOUT)   :: ocean
    REAL(dp),                       INTENT(IN)      :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_dT_timeframes_from_curve'
    INTEGER                                            :: ti0, ti1, tii, ncid, nt
    CHARACTER(LEN=256)                                 :: str
    REAL(dp)                                           :: dt_min 

    ! Add routine to path
    CALL init_routine( routine_name)

    
    ! Update sea level
    ! Find timeframe closest to desired time
    nt = size(ocean%transient%dT_series_time)
    if (ocean%transient%dT_series_time( 1) > time) then
      ! Desired time beyond lower limit
      ti0 = 1
    elseif (ocean%transient%dT_series_time( nt) < time) then
      ! Desired time beyond upper limit
      ti0 = nt
    else
      ! Desired time is within the file time
      dt_min = huge( 1._dp)
      do tii = 1, nt
        if (abs( ocean%transient%dT_series_time( tii) - time) < dt_min) then
          ti0 = tii
          dt_min = abs( ocean%transient%dT_series_time( tii) - time)
        end if
      end do
      if (dt_min > 0._dp) then
        !call warning('desired timeframe at t = {dp_01} not present in sea level record; reading data from closest match at t = {dp_02} instead!', &
        !  dp_01 = time, dp_02 = forcing%sea_level_time( ti0))
      end if
    end if
      
    
    ocean%transient%dT_t0    = ocean%transient%dT_series_time(ti0)
    ocean%transient%dT_at_t0 = ocean%transient%dT_series(ti0)
      
    ! if the desired time is after t0, we take one record after for t1
    if (time >= ocean%transient%dT_t0) then
      if (ti0 == size(ocean%transient%dT_series_time)) then
        IF (par%primary) WRITE(0,*) 'desired timeframe is at or beyond the last record. Using last available value for both timeframes...'
        ocean%transient%dT_t1    = ocean%transient%dT_series_time(ti0)
        ocean%transient%dT_at_t0 = ocean%transient%dT_series(ti0)
      else
        ocean%transient%dT_t1    = ocean%transient%dT_series_time(ti0+1)
        ocean%transient%dT_at_t1 = ocean%transient%dT_series(ti0+1)
      end if
    else
      ! otherwise we read one record before for t0, and that record is t1
      if (ti0 == 1) then
        IF (par%primary) WRITE(0,*) 'desired timeframe is at or before the first record. Using first available value for both timeframes...'
        ocean%transient%dT_t1    = ocean%transient%dT_series_time(ti0)
        ocean%transient%dT_at_t1 = ocean%transient%dT_series(ti0)
      else
        ocean%transient%dT_t1    = ocean%transient%dT_series_time(ti0)
        ocean%transient%dT_at_t1 = ocean%transient%dT_series(ti0)
        ocean%transient%dT_t0    = ocean%transient%dT_series_time(ti0-1)
        ocean%transient%dT_at_t0 = ocean%transient%dT_series(ti0-1)
      end if
    end if

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_sealevel_timeframes_from_curve

  subroutine run_ocean_model_transient(mesh, ocean, time)
  ! Apply the dT anomalies to the snapshot ocean

    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ocean_model),                 intent(inout) :: ocean
    real(dp),                               intent(in)    :: time

    ! Local variables
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'run_ocean_model_transient'
    REAL(dp)                                           :: wt0,wt1, dT_at_time

    ! Add routine to path
    CALL init_routine( routine_name)


    if (ocean%transient%dT_t1 == ocean%transient%dT_t0) then
      wt0 = 0._dp
      wt1 = 1._dp
    else
      if (time > ocean%transient%dT_t1) then
        wt0 = 0._dp
      elseif (time < ocean%transient%dT_t0) then
        wt0 = 1._dp
      else
        wt0 = (ocean%transient%dT_t1 - time) / (ocean%transient%dT_t1 - ocean%transient%dT_t0)
      end if
      wt1 = 1._dp - wt0
    end if

    ! Interpolate the two timeframes - constant dT over the entire region
    dT_at_time = wt0 * ocean%transient%dT_at_t0 + wt1 * ocean%transient%dT_at_t1

    do vi = mesh%vi1, mesh%vi2
      do z = 1, C%nz_ocean
        ocean%T(vi, z) = ocean%transient%T0(vi, z) + dT_at_time 
      end do
    end do


    ! Finalise routine path
    CALL finalise_routine( routine_name)

  end subroutine run_ocean_model_transient

end module ocean_realistic
