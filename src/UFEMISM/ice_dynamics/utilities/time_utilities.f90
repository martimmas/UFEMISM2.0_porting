module time_utilities
  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string, warning, crash
  use model_configuration, only: C
  use parameters, only: sec_per_year

  implicit none

  function days_since_ISMIP_basetime( time) result( ndays)
    ! Calculate the number of days since ISMIP basetime
    ! Assume basetime equals t = 0

    ! In/output variables:
    real(dp),                            intent(in)    :: time
    real(dp)                                           :: ndays

    ! Local variables:
    ! character(len=1024), parameter :: routine_name = 'days_since_ISMIP_basetime'

    ! Add routine to path
    ! call init_routine( routine_name)

    ndays = time * 360._dp

    ! Finalise routine path
    ! call finalise_routine( routine_name)

  end function days_since_ISMIP_basetime

end module time_utilities