module ocean_snapshot_plus_unif_dT

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
  use series_utilities
  use ocean_realistic
  use ocean_GlacialIndex

  implicit none

contains

! ===== Main routines =====
! =========================
subroutine initialise_ocean_model_snapshot_plus_unif_dT(mesh, ice, ocean, region_name)
  ! Initialise the ocean model
    !
    ! Using a snapshot plus a static, spatially uniform deltaT

    ! In- and output variables
    type(type_mesh),                        intent(in)    :: mesh
    type(type_ice_model),                   intent(in)    :: ice
    type(type_ocean_model),                 intent(inout) :: ocean
    character(len=3),                       intent(in)    :: region_name

    ! Local variables:
    character(len=256), parameter                         :: routine_name = 'initialise_ocean_model_snapshot_plus_unif_dT'
    integer                                               :: vi, z

    ! Add routine to path
    call init_routine( routine_name)

    ! Read single-time data from external file
    select case (region_name)
      case ('NAM')
        ocean%deltaT%dT = C%ocean_uniform_deltaT_NAM
      case ('EAS')
        ocean%deltaT%dT = C%ocean_uniform_deltaT_EAS
      case ('GRL')
        ocean%deltaT%dT = C%ocean_uniform_deltaT_GRL
      case ('ANT')
        ocean%deltaT%dT = C%ocean_uniform_deltaT_ANT
      case default
        call crash('unknown region_name "' // region_name // '"')
    end select

    call initialise_ocean_model_snapshot(mesh, ice, ocean, region_name)
    ! adds the deltaT
    do vi = mesh%vi1, mesh%vi2
    do z = 1, C%nz_ocean
      ocean%T(vi, z) = ocean%T(vi, z) + ocean%deltaT%dT
    end do
    end do

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  end subroutine initialise_ocean_model_snapshot_plus_unif_dT

end module ocean_snapshot_plus_unif_dT