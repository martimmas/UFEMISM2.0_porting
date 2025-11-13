module bed_roughness_nudging_H_dHdt_Akker

  use precisions, only: dp
  use control_resources_and_error_messaging, only: warning, crash, init_routine, finalise_routine
  use model_configuration, only: C
  use parameters
  use mesh_types, only: type_mesh
  use grid_basic, only: type_grid
  use ice_model_types, only: type_ice_model
  use reference_geometry_types, only: type_reference_geometry
  use bed_roughness_model_types, only: type_bed_roughness_model, type_bed_roughness_nudging_model_H_dHdt_Akker
  use mesh_utilities, only: extrapolate_Gaussian
  use nudging_utilities, only: calc_nudging_vs_extrapolation_masks

  implicit none

  private

  public :: initialise_bed_roughness_nudging_H_dHdt_Akker, run_bed_roughness_nudging_H_dHdt_Akker

contains

  subroutine run_bed_roughness_nudging_H_dHdt_Akker( mesh, grid_smooth, ice, &
    target_geometry, bed_roughness_prev, bed_roughness_next, nudge)
    ! Run the bed roughness nuding model based on local values of H and dH/dt (based on Akker et al. 2025)

    ! In/output variables:
    type(type_mesh),                                     intent(in   ) :: mesh
    type(type_grid),                                     intent(in   ) :: grid_smooth
    type(type_ice_model),                                intent(in   ) :: ice
    type(type_reference_geometry),                       intent(in   ) :: target_geometry
    real(dp), dimension(mesh%vi1:mesh%vi2),              intent(in   ) :: bed_roughness_prev
    real(dp), dimension(mesh%vi1:mesh%vi2),              intent(  out) :: bed_roughness_next
    type(type_bed_roughness_nudging_model_H_dHdt_Akker), intent(inout) :: nudge

    ! Local variables:
    character(len=256), parameter          :: routine_name = 'run_bed_roughness_nudging_H_dHdt_Akker'
    integer                                :: vi
    real(dp)                               :: w_Hb, phi_fric_relax, bed_roughness_target

    ! Add routine to path
    call init_routine( routine_name)

    call calc_nudging_vs_extrapolation_masks( mesh, ice, &
      nudge%mask_calc_dCdt_from_nudging, &
      nudge%mask_calc_dCdt_from_extrapolation, &
      nudge%mask_extrapolation)

    nudge%dC_dt = 0._dp

    do vi = mesh%vi1, mesh%vi2

      if (nudge%mask_calc_dCdt_from_nudging( vi)) then

        ! Bedrock elevation weight
        w_Hb = (ice%Hb( vi) - C%bednudge_H_dHdt_Akker_lowerHb) / (C%bednudge_H_dHdt_Akker_upperHb - C%bednudge_H_dHdt_Akker_lowerHb)
        ! Limit weight to [0,1]
        w_Hb = MIN( 1._dp, MAX( 0._dp, w_Hb))
        ! Compute target relaxation friction based on bedrock elevation
        bed_roughness_target = (1._dp - w_Hb) * C%bednudge_H_dHdt_Akker_minCr + w_Hb * C%bednudge_H_dHdt_Akker_maxCr

        ! Tim's big equation
        nudge%dC_dt( vi) = -bed_roughness_prev( vi) * (&
            (ice%Hs( vi) - target_geometry%Hs( vi)) / (C%bednudge_H_dHdt_Akker_H0 * C%bednudge_H_dHdt_Akker_tau) &
          + (2._dp / C%bednudge_H_dHdt_Akker_H0 * ice%dHs_dt( vi)) &
           - C%bednudge_H_dHdt_Akker_r / C%bednudge_H_dHdt_Akker_tau * log( bed_roughness_prev( vi) / bed_roughness_target))

      end if

    end do

    ! Perform the extrapolation - mask: 2 -> use as seed; 1 -> extrapolate; 0 -> ignore
    call extrapolate_Gaussian( mesh, nudge%mask_extrapolation, nudge%dC_dt, 50e3_dp)

    ! Calculate predicted bed roughness at t+dt
    bed_roughness_next = max( C%generic_bed_roughness_min, min( C%generic_bed_roughness_max, &
      bed_roughness_prev + C%bed_roughness_nudging_dt * nudge%dC_dt ))

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine run_bed_roughness_nudging_H_dHdt_Akker

  subroutine initialise_bed_roughness_nudging_H_dHdt_Akker( mesh, nudge)
    ! Initialise the bed roughness nudging model based on local values of H and dH/dt (based on Akker et al. 2025)

    ! In/output variables:
    type(type_mesh),                                    intent(in   ) :: mesh
    type(type_bed_roughness_nudging_model_H_dHdt_Akker), intent(inout) :: nudge

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_bed_roughness_nudging_H_dHdt_Akker'

    ! Add routine to path
    call init_routine( routine_name)

    ! Nudging masks
    allocate( nudge%mask_calc_dCdt_from_nudging      ( mesh%vi1:mesh%vi2), source = .false.)
    allocate( nudge%mask_calc_dCdt_from_extrapolation( mesh%vi1:mesh%vi2), source = .false.)
    allocate( nudge%mask_extrapolation               ( mesh%vi1:mesh%vi2), source = 0)

    ! Intermediate terms
    !allocate( nudge%C       ( mesh%vi1:mesh%vi2), source = 0._dp)
    !allocate( nudge%Laplac_C( mesh%vi1:mesh%vi2), source = 0._dp)
    allocate( nudge%dC_dt   ( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_bed_roughness_nudging_H_dHdt_Akker

end module bed_roughness_nudging_H_dHdt_Akker