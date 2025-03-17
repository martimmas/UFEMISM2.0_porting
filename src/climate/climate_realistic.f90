MODULE climate_realistic

  ! Realistic climate models

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string, warning, insert_val_into_string_int,insert_val_into_string_dp
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE ice_model_types                                        , ONLY: type_ice_model
  USE climate_model_types                                    , ONLY: type_climate_model, type_global_forcing
  use netcdf_io_main

  IMPLICIT NONE

CONTAINS

! ===== Main routines =====
! =========================

  SUBROUTINE run_climate_model_realistic( mesh, ice, climate, forcing, time)
    ! Calculate the climate
    !
    ! Use an realistic climate scheme

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_ice_model),                   INTENT(IN)    :: ice
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    TYPE(type_global_forcing),              INTENT(INOUT) :: forcing
    REAL(dp),                               INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'run_climate_model_realistic'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Run the chosen realistic climate model
    IF     (C%choice_climate_model_realistic == 'snapshot') THEN
      ! Do nothing
    ELSEIF (C%choice_climate_model_realistic == 'climate_matrix') THEN
      ! This is probably where we will update insolation, CO2, etc...
      CALL crash('choice_climate_model_realistic climate_matrix not implemented yet!"')
      CALL get_insolation_at_time( mesh, time, forcing, climate%Q_TOA)
      !CALL get_climate_at_time( mesh, time, forcing, climate)
    ELSE
      CALL crash('unknown choice_climate_model_realistic "' // TRIM( C%choice_climate_model_realistic) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE run_climate_model_realistic

  SUBROUTINE initialise_climate_model_realistic( mesh, climate, forcing, region_name)
    ! Initialise the climate model
    !
    ! Use a realistic climate scheme

    IMPLICIT NONE

    ! In- and output variables
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_climate_model),               INTENT(INOUT) :: climate
    TYPE(type_global_forcing),              INTENT(OUT)   :: forcing
    CHARACTER(LEN=3),                       INTENT(IN)    :: region_name

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'initialise_climate_model_realistic'
    CHARACTER(LEN=256)                                    :: filename_climate_snapshot
    REAL(dp)                                              :: timeframe_init_insolation

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Print to terminal
    IF (par%master)  WRITE(*,"(A)") '     Initialising realistic climate model "' // &
      colour_string( TRIM( C%choice_climate_model_realistic),'light blue') // '"...'

    ! Run the chosen realistic climate model
    IF (C%choice_climate_model_realistic == 'snapshot') THEN
      ! Read single-time data from external file

      ! Determine which climate model to initialise for this region
      IF     (region_name == 'NAM') THEN
        filename_climate_snapshot = C%filename_climate_snapshot_NAM
      ELSEIF (region_name == 'EAS') THEN
        filename_climate_snapshot = C%filename_climate_snapshot_EAS
      ELSEIF (region_name == 'GRL') THEN
        filename_climate_snapshot = C%filename_climate_snapshot_GRL
      ELSEIF (region_name == 'ANT') THEN
        filename_climate_snapshot = C%filename_climate_snapshot_ANT
      ELSE
        CALL crash('unknown region_name "' // region_name // '"')
      END IF

      CALL read_field_from_file_2D( filename_climate_snapshot, 'Hs', mesh, climate%Hs)
      CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'T2m', mesh, climate%T2m)
      CALL read_field_from_file_2D_monthly( filename_climate_snapshot, 'Precip', mesh, climate%Precip)
      allocate(climate%lambda( mesh%vi1:mesh%vi2))
      allocate(climate%Q_TOA(  mesh%vi1:mesh%vi2,12))
      allocate(climate%Albedo( mesh%vi1:mesh%vi2,12))
      allocate(climate%I_abs(  mesh%vi1:mesh%vi2))

      IF (par%master)  WRITE(*,"(A)") '     Initialising global forcings...'
      CALL initialise_global_forcings( mesh, forcing)

      ! If the simulation is properly set up with times in [ka], we just get the absolute value of the initial time
      ! TODO: what is the standard? time in [ka] or in "[a]"
      IF (C%start_time_of_run < 0._dp) THEN
        timeframe_init_insolation = C%start_time_of_run
      ELSE
        timeframe_init_insolation = 0._dp
      END IF
      IF (par%master)  WRITE(*,"(A)") '     Calling getting insolation at time...'
      CALL get_insolation_at_time( mesh, timeframe_init_insolation, forcing, climate%Q_TOA) ! TODO: check logic

    ELSE
      CALL crash('unknown choice_climate_model_realistic "' // TRIM( C%choice_climate_model_realistic) // '"')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_climate_model_realistic

  SUBROUTINE initialise_global_forcings( mesh, forcing)
    ! initialise the forcing structure to get d18O, CO2, insolation, etc...

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                   INTENT(IN)    :: mesh
    TYPE(type_global_forcing),         INTENT(OUT)   :: forcing

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                    :: routine_name = 'initialise_global_forcings'

    ! Add routine to path
    CALL init_routine( routine_name)

    ! TODO: checks with what exactly we need to load here to know which global forcings need to be read
    ! e.g., insolation, CO2, d18O, etc...
    ! read and load the insolation data
    CALL initialise_insolation_forcing( forcing, mesh)

    ! CO2 record - not yet implemented
    
    ! d18O record - not yet implemented

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_global_forcings

  ! == Insolation
  SUBROUTINE initialise_insolation_forcing( forcing, mesh)
    ! initialise the insolation series in the forcing structure

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_global_forcing),           INTENT(OUT) :: forcing
    TYPE(type_mesh),                      INTENT(IN) :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                    :: routine_name = 'initialise_insolation_forcing'
    CHARACTER(LEN=256)                               :: str
    INTEGER                                          :: ncid
    REAL(dp)                                         :: closest_t0

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_insolation_forcing == 'none') THEN
      ! No insolation included, likely because we're running an idealised-geometry experiment
    ELSEIF (C%choice_insolation_forcing == 'static' .OR. &
            C%choice_insolation_forcing == 'realistic') THEN
      
      ! Initialise insolation
      ! The times at which we have insolation fields from Laskar, between which we'll interpolate
      ! to find the insolation at model time (assuming that t0 <= model_time <= t1)
      IF (par%master)   WRITE(0,*) ' Initialising insolation data from ', TRIM(C%filename_insolation), '...'

      ! Memory allocation
      ALLOCATE( forcing%ins_t0)
      ALLOCATE( forcing%ins_t1)
      ALLOCATE( forcing%ins_ti0)
      ALLOCATE( forcing%ins_ti1)
      ALLOCATE( forcing%ins_nlat)
      ALLOCATE( forcing%ins_nlon)
      ALLOCATE(forcing%ins_lat            (   forcing%ins_nlat))
      ALLOCATE(forcing%ins_Q_TOA0         (mesh%vi1:mesh%vi2,12))
      ALLOCATE(forcing%ins_Q_TOA1         (mesh%vi1:mesh%vi2,12))
      forcing%ins_t0     = C%start_time_of_run
      forcing%ins_nlat   = 181
      forcing%ins_nlon   = 360
      forcing%ins_lat    = 0._dp
      forcing%ins_Q_TOA0 = 0._dp
      forcing%ins_Q_TOA1 = 0._dp
      
      ! Read the fields at ins_t0
      call read_field_from_file_0D( C%filename_insolation, field_name_options_time, closest_t0, time_to_read = forcing%ins_t0)

      if (par%master) WRITE(0,*) '     Reading Q_TOA0...'
      call read_field_from_file_1D_monthly( C%filename_insolation, field_name_options_insolation, mesh, forcing%ins_Q_TOA0, time_to_read = forcing%ins_t0)
      
      ! if the start time is after the closest t0, we read one record after for t1
      if (C%start_time_of_run >= closest_t0) then
        if (par%master) WRITE(0,*) '     start time is after closest ins_t0, reading one step further...'
        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, forcing%ins_t1, time_to_read = C%start_time_of_run+1000._dp)
      else
        ! otherwise we read one record before for t1
        if (par%master) WRITE(0,*) '     start time is before closest ins_t0, reading one step earlier...'
        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, forcing%ins_t1, time_to_read = C%start_time_of_run-1000._dp)
      end if

      call read_field_from_file_1D_monthly( C%filename_insolation, field_name_options_insolation, mesh, forcing%ins_Q_TOA1, time_to_read = forcing%ins_t1)

      IF (C%start_time_of_run < forcing%ins_t0) THEN
        CALL warning(' Model time starts before start of insolation record; the model will crash (lol) or use the first available record, which might be awfully wrong')
      END IF
      IF (C%end_time_of_run > forcing%ins_t1) THEN
        CALL warning(' Model time will reach beyond end of insolation record; constant extrapolation will be used in that case!')
      END IF

    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_insolation_forcing

  SUBROUTINE get_insolation_at_time( mesh, time, forcing, Q_TOA)
    ! Get monthly insolation at time t on the regional grid

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                        INTENT(IN)    :: mesh
    TYPE(type_global_forcing),              INTENT(INOUT) :: forcing
    REAL(dp),                               INTENT(IN)    :: time
    REAL(dp), DIMENSION(:,:), ALLOCATABLE,  INTENT(OUT)   :: Q_TOA

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                    :: routine_name = 'get_insolation_at_time'
    REAL(dp)                                         :: time_applied
    INTEGER                                          :: vi,m !,ilat_l,ilat_u
    REAL(dp)                                         :: wt0, wt1!, wlat_l, wlat_u ! not necessary?
    ! REAL(dp), DIMENSION(:  ), ALLOCATABLE            ::  Q_TOA_int ! not necessary?
    ! INTEGER                                          :: wQ_TOA_int ! not necessary?

    ! Add routine to path
    CALL init_routine( routine_name)

    time_applied = 0._dp

    ! Safety
    IF     (C%choice_insolation_forcing == 'none') THEN
      CALL crash('insolation should not be used when choice_insolation_forcing = "none"!')
    ELSEIF (C%choice_insolation_forcing == 'static') THEN
      time_applied = C%static_insolation_time
    ELSEIF (C%choice_insolation_forcing == 'realistic') THEN
      time_applied = time
    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    IF (time_applied < forcing%ins_t0 .OR. time_applied > forcing%ins_t1) THEN
      IF (par%master)  WRITE(0,*) '   Model time is out of the current insolation timeframes. Updating timeframes...'
      CALL update_insolation_timeframes_from_file( forcing, time_applied, mesh)
    END IF

    ! TODO: Is it necessary to allocate shared memory for timeframe-interpolated lat-month-only insolation?
    ALLOCATE(Q_TOA             (mesh%vi1:mesh%vi2,12))

    ! Calculate timeframe interpolation weights
    wt0 = (forcing%ins_t1 - time_applied) / (forcing%ins_t1 - forcing%ins_t0)
    wt1 = 1._dp - wt0

    ! Interpolate the two timeframes
    do vi = mesh%vi1, mesh%vi2
      do m = 1, 12
        Q_TOA(vi, m) = wt0 * forcing%ins_Q_TOA0(vi, m) + wt1 * forcing%ins_Q_TOA1(vi, m)
      end do
    end do

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE get_insolation_at_time

  SUBROUTINE update_insolation_timeframes_from_file( forcing, time, mesh)
    ! Read the NetCDF file containing the insolation forcing data. Only read the time frames enveloping the current
    ! coupling timestep to save on memory usage. Only done by master.

    ! NOTE: assumes time in forcing file is in kyr

    IMPLICIT NONE

    TYPE(type_mesh),                      INTENT(IN)   :: mesh
    TYPE(type_global_forcing),         INTENT(INOUT)   :: forcing
    REAL(dp),                             INTENT(IN)   :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_insolation_timeframes_from_file'
    INTEGER                                            :: ti0, ti1, ncid
    CHARACTER(LEN=256)                                 :: str

    ! Add routine to path
    CALL init_routine( routine_name)

    IF     (C%choice_insolation_forcing == 'none') THEN
      CALL crash('insolation should not be used when choice_insolation_forcing = "none"!')
    ELSEIF (C%choice_insolation_forcing == 'static' .OR. &
            C%choice_insolation_forcing == 'realistic') THEN

      ! Update insolation
      ! Find time indices to be read
      IF (par%master) THEN

        call read_field_from_file_0D( C%filename_insolation, field_name_options_time, forcing%ins_t0, time_to_read = time)
        call read_field_from_file_1D_monthly( C%filename_insolation, field_name_options_insolation, mesh, forcing%ins_Q_TOA0, time_to_read = forcing%ins_t0)
        
        ! if the desired time is after t0, we read one record after for t1
        if (time >= forcing%ins_t0) then
          call read_field_from_file_0D( C%filename_insolation, field_name_options_time, forcing%ins_t1, time_to_read = time+1000._dp)
        else
        ! otherwise we read one record before for t1
          call read_field_from_file_0D( C%filename_insolation, field_name_options_time, forcing%ins_t1, time_to_read = time-1000._dp)
        end if

      END IF ! IF (par%master) THEN

      call read_field_from_file_1D_monthly( C%filename_insolation, field_name_options_insolation, mesh, forcing%ins_Q_TOA1, time_to_read = forcing%ins_t1)

    ELSE
      CALL crash('unknown choice_insolation_forcing "' // TRIM( C%choice_insolation_forcing) // '"!')
    END IF

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_insolation_timeframes_from_file
  
! == Prescribed CO2 record
  SUBROUTINE update_CO2_at_model_time( time)
    ! Interpolate the data in forcing%CO2 to find the value at the queried time.
    ! If time lies outside the range of forcing%CO2_time, return the first/last value
    !
    ! NOTE: assumes time is listed in yr BP, so LGM would be -21000.0, and 0.0 corresponds to January 1st 1900.
    !
    ! NOTE: calculates average value over the preceding 30 years. For paleo this doesn't matter
    !       in the least, but for the historical period this makes everything more smooth.

    IMPLICIT NONE

    ! In/output variables:
    REAL(dp),                            INTENT(IN)    :: time

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'update_CO2_at_model_time'
    INTEGER                                            :: ti1, ti2, til, tiu
    REAL(dp)                                           :: a, b, tl, tu, intCO2, dintCO2
    REAL(dp), PARAMETER                                :: dt_smooth = 60._dp

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      ! Observed CO2 is needed for these forcing methods.
    ELSE
      CALL crash('should only be called when choice_forcing_method = "CO2_direct"!')
    END IF

    IF (par%master) THEN

      IF     (time < MINVAL( forcing%CO2_time)) THEN
        ! Model time before start of CO2 record; using constant extrapolation
        forcing%CO2_obs = forcing%CO2_record( 1)
      ELSEIF (time > MAXVAL( forcing%CO2_time)) THEN
        ! Model time beyond end of CO2 record; using constant extrapolation
        forcing%CO2_obs = forcing%CO2_record( C%CO2_record_length)
      ELSE

        ! Find range of raw time frames enveloping model time
        ti1 = 1
        DO WHILE (forcing%CO2_time( ti1) < time - dt_smooth .AND. ti1 < C%CO2_record_length)
          ti1 = ti1 + 1
        END DO
        ti1 = MAX( 1, ti1 - 1)

        ti2 = 2
        DO WHILE (forcing%CO2_time( ti2) < time             .AND. ti2 < C%CO2_record_length)
          ti2 = ti2 + 1
        END DO

        ! Calculate conservatively-remapped time-averaged CO2
        intCO2 = 0._dp
        DO til = ti1, ti2 - 1
          tiu = til + 1

          ! Linear interpolation between til and tiu: CO2( t) = a + b*t
          b = (forcing%CO2_record( tiu) - forcing%CO2_record( til)) / (forcing%CO2_time( tiu) - forcing%CO2_time( til))
          a = forcing%CO2_record( til) - b*forcing%CO2_time( til)

          ! Window of overlap between [til,tiu] and [t - dt_smooth, t]
          tl = MAX( forcing%CO2_time( til), time - dt_smooth)
          tu = MIN( forcing%CO2_time( tiu), time            )
          dintCO2 = (tu - tl) * (a + b * (tl + tu) / 2._dp)
          intCO2 = intCO2 + dintCO2
        END DO
        forcing%CO2_obs = intCO2 / dt_smooth

      END IF

    END IF
    CALL sync

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE update_CO2_at_model_time
  SUBROUTINE initialise_CO2_record
    ! Read the CO2 record specified in C%filename_CO2_record. Assumes this is an ASCII text file with at least two columns (time in kyr and CO2 in ppmv)
    ! and the number of rows being equal to C%CO2_record_length

    ! NOTE: assumes time is listed in kyr BP (so LGM would be -21.0); converts to yr after reading!

    IMPLICIT NONE

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'initialise_CO2_record'
    INTEGER                                            :: i,ios

    ! Add routine to path
    CALL init_routine( routine_name)

    ! Safety
    IF     (C%choice_forcing_method == 'CO2_direct') THEN
      ! Observed CO2 is needed for these forcing methods.
    ELSE
      CALL crash('should only be called when choice_forcing_method = "CO2_direct"!')
    END IF

    ! Allocate shared memory to take the data
    allocate( forcing%CO2_time(   C%CO2_record_length))
    allocate( forcing%CO2_record( C%CO2_record_length))
    CALL allocate_shared_dp_1D( C%CO2_record_length, forcing%CO2_time,   forcing%wCO2_time  )
    CALL allocate_shared_dp_1D( C%CO2_record_length, forcing%CO2_record, forcing%wCO2_record)
    CALL allocate_shared_dp_0D(                      forcing%CO2_obs,    forcing%wCO2_obs   )

    ! Read CO2 record (time and values) from specified text file
    IF (par%master) THEN

      WRITE(0,*) ' Reading CO2 record from ', TRIM(C%filename_CO2_record), '...'
! check this?!
      OPEN(   UNIT = 1337, FILE=C%filename_CO2_record, ACTION='READ')

      DO i = 1, C%CO2_record_length
        READ( UNIT = 1337, FMT=*, IOSTAT=ios) forcing%CO2_time( i), forcing%CO2_record( i)
        IF (ios /= 0) THEN
          CALL crash('length of text file "' // TRIM(C%filename_CO2_record) // '" does not match C%CO2_record_length!')
        END IF
      END DO

      CLOSE( UNIT  = 1337)

      IF (C%start_time_of_run/1000._dp < forcing%CO2_time(1)) THEN
         CALL warning(' Model time starts before start of CO2 record; constant extrapolation will be used in that case!')
      END IF
      IF (C%end_time_of_run/1000._dp > forcing%CO2_time(C%CO2_record_length)) THEN
         CALL warning(' Model time will reach beyond end of CO2 record; constant extrapolation will be used in that case!')
      END IF

      ! Convert from kyr to yr
      forcing%CO2_time = forcing%CO2_time * 1000._dp

    END IF ! IF (par%master)
    CALL sync

    ! Set the value for the current (starting) model time
    CALL update_CO2_at_model_time( C%start_time_of_run)

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE initialise_CO2_record

END MODULE climate_realistic
