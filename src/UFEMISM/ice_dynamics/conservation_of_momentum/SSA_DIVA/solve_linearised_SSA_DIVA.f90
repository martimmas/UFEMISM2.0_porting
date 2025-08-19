module solve_linearised_SSA_DIVA

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph_pair
  use solve_linearised_SSA_DIVA_infinite_slab, only: solve_SSA_DIVA_linearised_infinite_slab
  use solve_linearised_SSA_DIVA_ocean_pressure, only: solve_SSA_DIVA_linearised_ocean_pressure

  implicit none

  private

  public :: solve_SSA_DIVA_linearised

contains

  subroutine solve_SSA_DIVA_linearised( mesh, graphs, u_b, v_b, &
    Hi_a, Hb_a, SL_a, &
    N_b, dN_dx_b, dN_dy_b, &
    basal_friction_coefficient_b, tau_dx_b, tau_dy_b, u_b_prev, v_b_prev, &
    PETSc_rtol, PETSc_abstol, n_Axb_its, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    !< Solve the linearised SSA

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_graph_pair),                  intent(in   ) :: graphs
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: u_b, v_b
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(inout) :: Hi_a, Hb_a, SL_a
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: N_b, dN_dx_b, dN_dy_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: basal_friction_coefficient_b
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: tau_dx_b, tau_dy_b
    real(dp), dimension(mesh%nTri),         intent(inout) :: u_b_prev, v_b_prev
    real(dp),                               intent(in   ) :: PETSc_rtol, PETSc_abstol
    integer,                                intent(  out) :: n_Axb_its             ! Number of iterations used in the iterative solver
    integer,  dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_mask_b      ! Mask of triangles where velocity is prescribed
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_u_b         ! Prescribed velocities in the x-direction
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: BC_prescr_v_b         ! Prescribed velocities in the y-direction

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'solve_SSA_DIVA_linearised'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%BC_ice_front)
    case default
      call crash('unknown BC_ice_front "' // trim( C%BC_ice_front) // '"')
    case ('infinite_slab')
      call solve_SSA_DIVA_linearised_infinite_slab( mesh, u_b, v_b, N_b, dN_dx_b, dN_dy_b, &
        basal_friction_coefficient_b, tau_dx_b, tau_dy_b, u_b_prev, v_b_prev, &
        PETSc_rtol, PETSc_abstol, n_Axb_its, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    case ('ocean_pressure')
      call solve_SSA_DIVA_linearised_ocean_pressure( mesh, graphs, u_b, v_b, &
        Hi_a, Hb_a, SL_a, &
        N_b, dN_dx_b, dN_dy_b, &
        basal_friction_coefficient_b, tau_dx_b, tau_dy_b, u_b_prev, v_b_prev, &
        PETSc_rtol, PETSc_abstol, n_Axb_its, BC_prescr_mask_b, BC_prescr_u_b, BC_prescr_v_b)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSA_DIVA_linearised

end module solve_linearised_SSA_DIVA
