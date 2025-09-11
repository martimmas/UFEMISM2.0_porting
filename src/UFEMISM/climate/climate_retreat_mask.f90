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

  public :: run_climate_retreat_mask
  public :: update_climate_retreat_mask
  public :: initialise_climate_retreat_mask

contains

! == ISMIP-style Antarctica future retreat mask
! ======================================================================

  subroutine run_climate_retreat_mask( mesh, climate, time, ice)

    ! In/output variables:
    type(type_mesh),                    intent(in)    :: mesh
    type(type_climate),                 intent(inout) :: climate
    real(dp),                           intent(in)    :: time
    type(type_ice_model),               intent(in)    :: ice

    ! Local variables:
    character(len=256), parameter                           :: routine_name = 'run_climate_retreat_mask'
    integer                                                 :: vi

    ! Add routine to path
    call init_routine( routine_name)

    ! Check if the requested time is enveloped by the two timeframes;
    ! if not, read the two relevant timeframes from the NetCDF file
    if (time < climate%ISMIP_style%shelf_collapse_mask_t0 .OR. time > climate%ISMIP_style%shelf_collapse_mask_t1) then

      ! Find and read the two global time frames
      !call sync
      ! this needs to be fixed
      call update_ISMIP_style_future_timeframes( mesh, climate, time)

    end if ! IF (time >= climate_matrix%SMB_direct%t0 .AND. time <= climate_matrix%SMB_direct%t1) THEN

    ! Interpolate the two timeframes in time
    wt0 = (climate%ISMIP_style%shelf_collapse_mask_t1 - time) / (climate%ISMIP_style%shelf_collapse_mask_t1 - climate%ISMIP_style%shelf_collapse_mask_t0)
    wt1 = 1._dp - wt0

    do vi = mesh%vi1, mesh%vi2

      climate%ISMIP_style%shelf_collapse_mask( vi) = (wt0 * climate%ISMIP_style%shelf_collapse_mask0( vi)) + &
                                                        (wt1 * climate%ISMIP_style%shelf_collapse_mask1( vi))

    end do
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_climate_retreat_mask

  subroutine update_ISMIP_style_future_timeframes( mesh, climate, time)
    ! Update the two timeframes of the ISMIP-style forcing data

    ! In/output variables:
    type(type_mesh),                          intent(in   ) :: mesh
    type(type_climate),                       intent(inout) :: climate
    real(dp),                                 intent(in)    :: time

    ! Local variables:
    character(len=256), parameter                           :: routine_name = 'update_ISMIP_style_future_timeframes'

    ! Add routine to path
    call init_routine( routine_name)

    ! Read timeframes from file
    call read_field_from_file_2D( C%ISMIP_future_shelf_collapse_forcing_filename, 'mask', mesh, climate_matrix%ISMIP_style%shelf_collapse_mask0, 'ANT', time_to_read = REAL( FLOOR( time),dp)        )
    call read_field_from_file_2D( C%ISMIP_future_shelf_collapse_forcing_filename, 'mask', mesh, climate_matrix%ISMIP_style%shelf_collapse_mask1, 'ANT', time_to_read = REAL( FLOOR( time),dp) + 1._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine update_ISMIP_style_future_timeframes

! start with the initialise subroutine to understand what I need to restructure
  subroutine initialise_climate_retreat_mask( mesh, climate)
    ! Use the ISMIP-style forcing

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_climate_ISMIP_style),     intent(inout) :: climate

    ! Local variables:
    character(len=256), parameter                           :: routine_name = 'initialise_climate_retreat_mask'

    ! Add routine to path
    call init_routine( routine_name)

! I think all the code above is not needed... and even the IF is not needed anymore as I will call this part as something else from climate main...
      ! Allocate memory
    allocate( climate%ISMIP_style%shelf_collapse_mask0( mesh%vi1:mesh%vi2))
    allocate( climate%ISMIP_style%shelf_collapse_mask1( mesh%vi1:mesh%vi2))
    allocate( climate%ISMIP_style%shelf_collapse_mask( mesh%vi1:mesh%vi2))

      !CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%shelf_collapse_mask0, climate_matrix%ISMIP_style%wshelf_collapse_mask0)
      !CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%shelf_collapse_mask1, climate_matrix%ISMIP_style%wshelf_collapse_mask1)
      !CALL allocate_shared_dp_1D( mesh%nV, climate_matrix%ISMIP_style%shelf_collapse_mask , climate_matrix%ISMIP_style%wshelf_collapse_mask )

! what about this one? I think is not needed? as it is like a real value
      ! Allocate memory for the timestamps of the two timeframes
    allocate( climate%ISMIP_style%shelf_collapse_mask_t0)
    allocate( climate%ISMIP_style%shelf_collapse_mask_t1)
      !CALL allocate_shared_dp_0D( climate%ISMIP_style%shelf_collapse_mask_t0, climate_matrix%ISMIP_style%wshelf_collapse_mask_t0)
      !CALL allocate_shared_dp_0D( climate%ISMIP_style%shelf_collapse_mask_t1, climate_matrix%ISMIP_style%wshelf_collapse_mask_t1)
! I THINK NEXT LINE IS NOT NEEDED...
      if (par%master) then
        ! Give impossible values to timeframes, so that the first call to run_climate_model_ISMIP_style
        ! is guaranteed to first read two new timeframes from the NetCDF file
        climate%ISMIP_style%shelf_collapse_mask_t0 = C%start_time_of_run - 100._dp
        climate%ISMIP_style%shelf_collapse_mask_t1 = C%start_time_of_run - 90._dp
      end if ! IF (par%master) THEN
      call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_climate_retreat_mask









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