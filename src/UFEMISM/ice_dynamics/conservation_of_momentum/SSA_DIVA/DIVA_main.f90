module DIVA_main

  ! Routines for calculating ice velocities using the Depth-Integrated Viscosity Approximation (DIVA)

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model, type_ice_velocity_solver_DIVA
  use bed_roughness_model_types, only: type_bed_roughness_model
  use DIVA_infinite_slab, only: initialise_DIVA_solver_infinite_slab, solve_DIVA_infinite_slab, &
    remap_DIVA_solver_infinite_slab, initialise_DIVA_velocities_from_file_infinite_slab, &
    create_restart_file_DIVA_infinite_slab, write_to_restart_file_DIVA_infinite_slab

  implicit none

  private

  public :: initialise_DIVA_solver, solve_DIVA, remap_DIVA_solver, &
    create_restart_file_DIVA, write_to_restart_file_DIVA

contains

  ! == Main routines

  subroutine initialise_DIVA_solver( mesh, DIVA, region_name)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(  out) :: DIVA
    character(len=3),                    intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_DIVA_solver'
    character(len=256)             :: choice_initial_velocity

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%BC_ice_front)
    case default
      call crash('unknown BC_ice_front "' // trim( C%BC_ice_front) // '"')
    case ('infinite_slab')
      call initialise_DIVA_solver_infinite_slab( mesh, DIVA, region_name)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_DIVA_solver

  subroutine solve_DIVA( mesh, ice, bed_roughness, DIVA, n_visc_its, n_Axb_its, &
    BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    !< Calculate ice velocities by solving the Depth-Integrated Viscosity Approximation

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_model),                intent(inout) :: ice
    type(type_bed_roughness_model),      intent(in   ) :: bed_roughness
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA
    integer,                             intent(  out) :: n_visc_its               ! Number of non-linear viscosity iterations
    integer,                             intent(  out) :: n_Axb_its                ! Number of iterations in iterative solver for linearised momentum balance
    integer,  dimension(:), optional,    intent(in   ) :: BC_prescr_mask_b         ! Mask of triangles where velocity is prescribed
    real(dp), dimension(:), optional,    intent(in   ) :: BC_prescr_u_b            ! Prescribed velocities in the x-direction
    real(dp), dimension(:), optional,    intent(in   ) :: BC_prescr_v_b            ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'solve_DIVA'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%BC_ice_front)
    case default
      call crash('unknown BC_ice_front "' // trim( C%BC_ice_front) // '"')
    case ('infinite_slab')
      call solve_DIVA_infinite_slab( mesh, ice, bed_roughness, DIVA, n_visc_its, n_Axb_its, &
        BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_DIVA

  subroutine remap_DIVA_solver( mesh_old, mesh_new, DIVA)

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh_old
    type(type_mesh),                     intent(in   ) :: mesh_new
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_DIVA_solver'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%BC_ice_front)
    case default
      call crash('unknown BC_ice_front "' // trim( C%BC_ice_front) // '"')
    case ('infinite_slab')
      call remap_DIVA_solver_infinite_slab( mesh_old, mesh_new, DIVA)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remap_DIVA_solver

  ! == Initialisation

  subroutine initialise_DIVA_velocities_from_file( mesh, DIVA, region_name)
    !< Initialise the velocities for the DIVA solver from an external NetCDF file

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA
    character(len=3),                    intent(in   ) :: region_name

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_DIVA_velocities_from_file'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%BC_ice_front)
    case default
      call crash('unknown BC_ice_front "' // trim( C%BC_ice_front) // '"')
    case ('infinite_slab')
      call initialise_DIVA_velocities_from_file_infinite_slab( mesh, DIVA, region_name)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_DIVA_velocities_from_file

  ! == Restart NetCDF files

  subroutine write_to_restart_file_DIVA( mesh, DIVA, time)
    ! Write to the restart NetCDF file for the DIVA solver

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(in   ) :: DIVA
    real(dp),                            intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_restart_file_DIVA'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%BC_ice_front)
    case default
      call crash('unknown BC_ice_front "' // trim( C%BC_ice_front) // '"')
    case ('infinite_slab')
      call write_to_restart_file_DIVA_infinite_slab( mesh, DIVA, time)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file_DIVA

  subroutine create_restart_file_DIVA( mesh, DIVA)
    ! Create a restart NetCDF file for the DIVA solver
    ! Includes generation of the procedural filename (e.g. "restart_DIVA_00001.nc")

    ! In/output variables:
    type(type_mesh),                     intent(in   ) :: mesh
    type(type_ice_velocity_solver_DIVA), intent(inout) :: DIVA

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_restart_file_DIVA'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%BC_ice_front)
    case default
      call crash('unknown BC_ice_front "' // trim( C%BC_ice_front) // '"')
    case ('infinite_slab')
      call create_restart_file_DIVA_infinite_slab( mesh, DIVA)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_restart_file_DIVA

end module DIVA_main
