module solve_linearised_SSA_DIVA_ocean_pressure

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph_pair
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use CSR_matrix_basics, only: allocate_matrix_CSR_dist, add_entry_CSR_dist, read_single_row_CSR_dist, &
    finalise_matrix_CSR_dist
  use mesh_utilities, only: find_ti_copy_ISMIP_HOM_periodic, find_ti_copy_SSA_icestream_infinite
  use mpi_distributed_memory, only: gather_to_all
  use petsc_basic, only: solve_matrix_equation_CSR_PETSc
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared, &
    gather_dist_shared_to_all
  use mpi_f08, only: MPI_WIN
  use mesh_graph_mapping, only: map_mesh_triangles_to_graph, map_graph_to_mesh_triangles, &
    map_mesh_vertices_to_graph_ghost_nodes
  use ice_geometry_basics, only: height_of_water_column_at_ice_front
  use parameters, only: ice_density, seawater_density, grav

  use netcdf_io_main

  implicit none

  private

  public :: solve_SSA_DIVA_linearised_ocean_pressure

contains

  subroutine solve_SSA_DIVA_linearised_ocean_pressure( mesh, graphs, u_b, v_b, &
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
    character(len=1024), parameter      :: routine_name = 'solve_SSA_DIVA_linearised_ocean_pressure'

    real(dp), dimension(:), pointer     :: u_b_g     => null()
    real(dp), dimension(:), pointer     :: v_b_g     => null()
    type(MPI_WIN)                       :: wu_b_g, wv_b_g
    real(dp), dimension(:), pointer     :: Hi_b_ggh  => null()
    real(dp), dimension(:), pointer     :: Hb_b_ggh  => null()
    real(dp), dimension(:), pointer     :: SL_b_ggh  => null()
    real(dp), dimension(:), pointer     :: Ho_b_ggh  => null()
    type(MPI_WIN)                       :: wHi_b_ggh, wHb_b_ggh, wSL_b_ggh, wHo_b_ggh
    real(dp), dimension(:), pointer     :: N_b_g     => null()
    real(dp), dimension(:), pointer     :: dN_dx_b_g => null()
    real(dp), dimension(:), pointer     :: dN_dy_b_g => null()
    type(MPI_WIN)                       :: wN_b_g, wdN_dx_b_g, wdN_dy_b_g
    real(dp), dimension(:), pointer     :: basal_friction_coefficient_b_g => null()
    real(dp), dimension(:), pointer     :: tau_dx_b_g => null()
    real(dp), dimension(:), pointer     :: tau_dy_b_g => null()
    type(MPI_WIN)                       :: wbasal_friction_coefficient_b_g, wtau_dx_b_g, wtau_dy_b_g

    integer                             :: ncols, ncols_loc, nrows, nrows_loc, nnz_est_proc
    type(type_sparse_matrix_CSR_dp)     :: A_CSR
    real(dp), dimension(:), allocatable :: bb_b_g
    real(dp), dimension(:), allocatable :: uv_buv_g
    integer                             :: row_niuv, ni, uv

    ! Add routine to path
    call init_routine( routine_name)

    ! Store the previous solution
    call gather_to_all( u_b, u_b_prev)
    call gather_to_all( v_b, v_b_prev)

    ! Allocate hybrid distrobuted/shared memory
    call allocate_dist_shared( u_b_g                         , wu_b_g                         , graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( v_b_g                         , wv_b_g                         , graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( Hi_b_ggh                      , wHi_b_ggh                      , graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( Hb_b_ggh                      , wHb_b_ggh                      , graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( SL_b_ggh                      , wSL_b_ggh                      , graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( Ho_b_ggh                      , wHo_b_ggh                      , graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( N_b_g                         , wN_b_g                         , graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( dN_dx_b_g                     , wdN_dx_b_g                     , graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( dN_dy_b_g                     , wdN_dy_b_g                     , graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( basal_friction_coefficient_b_g, wbasal_friction_coefficient_b_g, graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( tau_dx_b_g                    , wtau_dx_b_g                    , graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( tau_dy_b_g                    , wtau_dy_b_g                    , graphs%graph_b%pai%n_nih)

    ! Map ice model data from the mesh triangles to the b-graph
    call map_mesh_triangles_to_graph           ( mesh, u_b                         , graphs%graph_b, u_b_g                         )
    call map_mesh_triangles_to_graph           ( mesh, v_b                         , graphs%graph_b, v_b_g                         )
    call map_mesh_vertices_to_graph_ghost_nodes( mesh, Hi_a                        , graphs%graph_b, Hi_b_ggh                      )
    call map_mesh_vertices_to_graph_ghost_nodes( mesh, Hb_a                        , graphs%graph_b, Hb_b_ggh                      )
    call map_mesh_vertices_to_graph_ghost_nodes( mesh, SL_a                        , graphs%graph_b, SL_b_ggh                      )
    call map_mesh_triangles_to_graph           ( mesh, N_b                         , graphs%graph_b, N_b_g                         )
    call map_mesh_triangles_to_graph           ( mesh, dN_dx_b                     , graphs%graph_b, dN_dx_b_g                     )
    call map_mesh_triangles_to_graph           ( mesh, dN_dy_b                     , graphs%graph_b, dN_dy_b_g                     )
    call map_mesh_triangles_to_graph           ( mesh, basal_friction_coefficient_b, graphs%graph_b, basal_friction_coefficient_b_g)
    call map_mesh_triangles_to_graph           ( mesh, tau_dx_b                    , graphs%graph_b, tau_dx_b_g                    )
    call map_mesh_triangles_to_graph           ( mesh, tau_dy_b                    , graphs%graph_b, tau_dy_b_g                    )

    ! Calculate height of the water column in contact with the ice front
    do ni = graphs%graph_b%ni1, graphs%graph_b%ni2
      if (graphs%graph_b%is_ghost( ni)) then
        Ho_b_ggh( ni) = height_of_water_column_at_ice_front( Hi_b_ggh( ni), Hb_b_ggh( ni), SL_b_ggh( ni))
      end if
    end do

    ! == Initialise the stiffness matrix using the native UFEMISM CSR-matrix format
    ! =============================================================================

    ! Matrix size
    ncols           = graphs%graph_b%n     * 2      ! from
    ncols_loc       = graphs%graph_b%n_loc * 2
    nrows           = graphs%graph_b%n     * 2      ! to
    nrows_loc       = graphs%graph_b%n_loc * 2
    nnz_est_proc    = graphs%M2_ddx_b_b%nnz * 4

    call allocate_matrix_CSR_dist( A_CSR, nrows, ncols, nrows_loc, ncols_loc, nnz_est_proc)

    ! Allocate memory for the load vector and the solution
    allocate( bb_b_g(   graphs%graph_b%ni1*2-1: graphs%graph_b%ni2*2))
    allocate( uv_buv_g( graphs%graph_b%ni1*2-1: graphs%graph_b%ni2*2))

    ! Fill in the current velocity solution
    do ni = graphs%graph_b%ni1, graphs%graph_b%ni2

      ! u
      row_niuv = graphs%biuv2n( ni,1)
      uv_buv_g( row_niuv) = u_b_g( ni)

      ! v
      row_niuv = graphs%biuv2n( ni,2)
      uv_buv_g( row_niuv) = v_b_g( ni)

    end do

    ! == Construct the stiffness matrix for the linearised SSA
    ! ========================================================

    do row_niuv = A_CSR%i1, A_CSR%i2

      ni = graphs%n2biuv( row_niuv,1)
      uv = graphs%n2biuv( row_niuv,2)

      ! if (BC_prescr_mask_b( ti) == 1) then
      !   ! Dirichlet boundary condition; velocities are prescribed for this triangle

    !     ! Stiffness matrix: diagonal element set to 1
    !     call add_entry_CSR_dist( A_CSR, row_tiuv, row_tiuv, 1._dp)

    !     ! Load vector: prescribed velocity
    !     if     (uv == 1) then
    !       bb( row_tiuv) = BC_prescr_u_b( ti)
    !     elseif (uv == 2) then
    !       bb( row_tiuv) = BC_prescr_v_b( ti)
    !     else
    !       call crash('uv can only be 1 or 2!')
    !     end if

      if (graphs%graph_b%is_ghost( ni)) then
        ! Ice margin: apply ocean pressure boundary condition

        call calc_SSA_DIVA_stiffness_matrix_row_BC_ice_front( graphs, &
          N_b_g, Hi_b_ggh, Ho_b_ggh, A_CSR, bb_b_g, row_niuv)

      else
        ! No boundary conditions apply; solve the SSA

        if (C%do_include_SSADIVA_crossterms) then
          ! Calculate matrix coefficients for the full SSA/DIVA
          call calc_SSA_DIVA_stiffness_matrix_row_free( graphs, N_b_g, dN_dx_b_g, dN_dy_b_g, &
            basal_friction_coefficient_b_g, tau_dx_b_g, tau_dy_b_g, A_CSR, bb_b_g, row_niuv)
        else
          ! Calculate matrix coefficients for the SSA sans the gradients of the effective viscosity (the "cross-terms")
          call calc_SSA_DIVA_sans_stiffness_matrix_row_free( graphs, N_b_g, &
            basal_friction_coefficient_b_g, tau_dx_b_g, tau_dy_b_g, A_CSR, bb_b_g, row_niuv)
        end if

      end if

    end do

    call finalise_matrix_CSR_dist( A_CSR)

    ! ! == Solve the matrix equation
    ! ! ============================

    ! Use PETSc to solve the matrix equation
    call solve_matrix_equation_CSR_PETSc( A_CSR, bb_b_g, uv_buv_g, PETSc_rtol, PETSc_abstol, &
      n_Axb_its)

    ! Disentangle the u and v components of the velocity solution
    do ni = graphs%graph_b%ni1, graphs%graph_b%ni2

      ! u
      row_niuv = graphs%biuv2n( ni,1)
      u_b_g( ni) = uv_buv_g( row_niuv)

      ! v
      row_niuv = graphs%biuv2n( ni,2)
      v_b_g( ni) = uv_buv_g( row_niuv)

    end do

    ! Map velocity solution from the b-graph to the mesh triangles
    call map_graph_to_mesh_triangles( graphs%graph_b, u_b_g, mesh, u_b)
    call map_graph_to_mesh_triangles( graphs%graph_b, v_b_g, mesh, v_b)

    ! Clean up after yourself
    call deallocate_dist_shared( u_b_g                         , wu_b_g                         )
    call deallocate_dist_shared( v_b_g                         , wv_b_g                         )
    call deallocate_dist_shared( Hi_b_ggh                      , wHi_b_ggh                      )
    call deallocate_dist_shared( Hb_b_ggh                      , wHb_b_ggh                      )
    call deallocate_dist_shared( SL_b_ggh                      , wSL_b_ggh                      )
    call deallocate_dist_shared( Ho_b_ggh                      , wHo_b_ggh                      )
    call deallocate_dist_shared( N_b_g                         , wN_b_g                         )
    call deallocate_dist_shared( dN_dx_b_g                     , wdN_dx_b_g                     )
    call deallocate_dist_shared( dN_dy_b_g                     , wdN_dy_b_g                     )
    call deallocate_dist_shared( basal_friction_coefficient_b_g, wbasal_friction_coefficient_b_g)
    call deallocate_dist_shared( tau_dx_b_g                    , wtau_dx_b_g                    )
    call deallocate_dist_shared( tau_dy_b_g                    , wtau_dy_b_g                    )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine solve_SSA_DIVA_linearised_ocean_pressure

  subroutine calc_SSA_DIVA_stiffness_matrix_row_free( graphs, N_b_g, dN_dx_b_g, dN_dy_b_g, &
    basal_friction_coefficient_b_g, tau_dx_b_g, tau_dy_b_g, A_CSR, bb_b_g, row_niuv)
    !< Add coefficients to this matrix row to represent the linearised SSA

    ! The SSA reads;
    !
    !   d/dx [ 2 N ( 2 du/dx + dv/dy )] + d/dy [ N ( du/dy + dv/dx)] - beta_b u = -tau_dx
    !
    !   d/dy [ 2 N ( 2 dv/dy + du/dx )] + d/dx [ N ( dv/dx + du/dy)] - beta_b v = -tau_dy
    !
    ! Using the chain rule, this expands to read:
    !
    !   4 N d2u/dx2 + 4 dN/dx du/dx + 2 N d2v/dxdy + 2 dN/dx dv/dy + ...
    !     N d2u/dy2 +   dN/dy du/dy +   N d2v/dxdy +   dN/dy dv/dx - beta_b u = -tau_dx
    !
    !   4 N d2v/dy2 + 4 dN/dy dv/dy + 2 N d2u/dxdy + 2 dN/dy du/dx + ...
    !     N d2v/dx2 +   dN/dx dv/dx +   N d2u/dxdy +   dN/dx du/dy - beta_b v = -tau_dy
    !
    ! Rearranging to gather the terms involving u and v gives:
    !
    !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_b u + ...
    !   3 N d2v/dxdy + 2 dN/dx dv/dy +             dN/dy dv/dx = -tau_dx
    !
    !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_b v + ...
    !   3 N d2u/dxdy + 2 dN/dy du/dx +             dN/dx du/dy = -tau_dy

    ! In/output variables:
    type(type_graph_pair),                                                    intent(in   ) :: graphs
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: N_b_g, dN_dx_b_g, dN_dy_b_g
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: basal_friction_coefficient_b_g
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: tau_dx_b_g, tau_dy_b_g
    type(type_sparse_matrix_CSR_dp),                                          intent(inout) :: A_CSR
    real(dp), dimension(graphs%graph_b%ni1*2-1: graphs%graph_b%ni2*2),        intent(inout) :: bb_b_g
    integer,                                                                  intent(in   ) :: row_niuv

    ! Local variables:
    integer                             :: ni, uv
    real(dp)                            :: N, dN_dx, dN_dy, basal_friction_coefficient, tau_dx, tau_dy
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_ddx_val
    real(dp), dimension(:), allocatable :: single_row_ddy_val
    real(dp), dimension(:), allocatable :: single_row_d2dx2_val
    real(dp), dimension(:), allocatable :: single_row_d2dxdy_val
    real(dp), dimension(:), allocatable :: single_row_d2dy2_val
    integer                             :: single_row_nnz
    real(dp)                            :: Au, Av
    integer                             :: k, nj, col_nju, col_njv

    ! Relevant indices for this graph node
    ni = graphs%n2biuv( row_niuv,1)
    uv = graphs%n2biuv( row_niuv,2)

    ! N, dN/dx, dN/dy, basal_friction_coefficient_b, tau_dx, and tau_dy on this graph node
    N                          = N_b_g                         ( ni)
    dN_dx                      = dN_dx_b_g                     ( ni)
    dN_dy                      = dN_dy_b_g                     ( ni)
    basal_friction_coefficient = basal_friction_coefficient_b_g( ni)
    tau_dx                     = tau_dx_b_g                    ( ni)
    tau_dy                     = tau_dy_b_g                    ( ni)

    ! allocate memory for single matrix rows
    allocate( single_row_ind(        graphs%graph_a%nC_mem*2))
    allocate( single_row_ddx_val(    graphs%graph_a%nC_mem*2))
    allocate( single_row_ddy_val(    graphs%graph_a%nC_mem*2))
    allocate( single_row_d2dx2_val(  graphs%graph_a%nC_mem*2))
    allocate( single_row_d2dxdy_val( graphs%graph_a%nC_mem*2))
    allocate( single_row_d2dy2_val(  graphs%graph_a%nC_mem*2))

    ! Read coefficients of the operator matrices
    call read_single_row_CSR_dist( graphs%M2_ddx_b_b   , ni, single_row_ind, single_row_ddx_val   , single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_ddy_b_b   , ni, single_row_ind, single_row_ddy_val   , single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_d2dx2_b_b , ni, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_d2dxdy_b_b, ni, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_d2dy2_b_b , ni, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    if (uv == 1) then
      ! x-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        nj      = single_row_ind( k)
        col_nju = graphs%biuv2n( nj,1)
        col_njv = graphs%biuv2n( nj,2)

        !   4 N d2u/dx2  + 4 dN/dx du/dx + N d2u/dy2 + dN/dy du/dy - beta_b u + ...
        !   3 N d2v/dxdy + 2 dN/dx dv/dy +             dN/dy dv/dx = -tau_dx

        ! Combine the mesh operators
        Au = 4._dp * N     * single_row_d2dx2_val(  k) + &  ! 4  N    d2u/dx2
             4._dp * dN_dx * single_row_ddx_val(    k) + &  ! 4 dN/dx du/dx
                     N     * single_row_d2dy2_val(  k) + &  !    N    d2u/dy2
                     dN_dy * single_row_ddy_val(    k)      !   dN/dy du/dy
        if (nj == ni) Au = Au - basal_friction_coefficient  ! - beta_b u

        Av = 3._dp * N     * single_row_d2dxdy_val( k) + &  ! 3  N    d2v/dxdy
             2._dp * dN_dx * single_row_ddy_val(    k) + &  ! 2 dN/dx dv/dy
                     dN_dy * single_row_ddx_val(    k)      !   dN/dy dv/dx

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_niuv, col_nju, Au)
        call add_entry_CSR_dist( A_CSR, row_niuv, col_njv, Av)

      end do

      ! Load vector
      bb_b_g( row_niuv) = -tau_dx

    elseif (uv == 2) then
      ! y-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        nj      = single_row_ind( k)
        col_nju = graphs%biuv2n( nj,1)
        col_njv = graphs%biuv2n( nj,2)

        !   4 N d2v/dy2  + 4 dN/dy dv/dy + N d2v/dx2 + dN/dx dv/dx - beta_b v + ...
        !   3 N d2u/dxdy + 2 dN/dy du/dx +             dN/dx du/dy = -tau_dy

        ! Combine the mesh operators
        Av = 4._dp * N     * single_row_d2dy2_val(  k) + &  ! 4  N    d2v/dy2
             4._dp * dN_dy * single_row_ddy_val(    k) + &  ! 4 dN/dy dv/dy
                     N     * single_row_d2dx2_val(  k) + &  !    N    d2v/dx2
                     dN_dx * single_row_ddx_val(    k)      !   dN/dx dv/dx
        if (nj == ni) Av = Av - basal_friction_coefficient  ! - beta_b v

        Au = 3._dp * N     * single_row_d2dxdy_val( k) + &  ! 3  N    d2u/dxdy
             2._dp * dN_dy * single_row_ddx_val(    k) + &  ! 2 dN/dy du/dx
                     dN_dx * single_row_ddy_val(    k)      !   dN/dx du/dy

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_niuv, col_nju, Au)
        call add_entry_CSR_dist( A_CSR, row_niuv, col_njv, Av)

      end do

      ! Load vector
      bb_b_g( row_niuv) = -tau_dy

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_SSA_DIVA_stiffness_matrix_row_free

  subroutine calc_SSA_DIVA_sans_stiffness_matrix_row_free( graphs, N_b_g, basal_friction_coefficient_b_g, &
    tau_dx_b_g, tau_dy_b_g, A_CSR, bb_b_g, row_niuv)
    !< Add coefficients to this matrix row to represent the linearised SSA
    !< sans the gradients of the effective viscosity (the "cross-terms")

    ! The SSA reads;
    !
    !   d/dx [ 2 N ( 2 du/dx + dv/dy )] + d/dy [ N ( du/dy + dv/dx)] - beta_b u = -tau_dx
    !
    !   d/dy [ 2 N ( 2 dv/dy + du/dx )] + d/dx [ N ( dv/dx + du/dy)] - beta_b v = -tau_dy
    !
    ! Using the chain rule, this expands to read:
    !
    !   4 N d2u/dx2 + 4 dN/dx du/dx + 2 N d2v/dxdy + 2 dN/dx dv/dy + ...
    !     N d2u/dy2 +   dN/dy du/dy +   N d2v/dxdy +   dN/dy dv/dx - beta_b u = -tau_dx
    !
    !   4 N d2v/dy2 + 4 dN/dy dv/dy + 2 N d2u/dxdy + 2 dN/dy du/dx + ...
    !     N d2v/dx2 +   dN/dx dv/dx +   N d2u/dxdy +   dN/dx du/dy - beta_b v = -tau_dy
    !
    ! The "sans" approximation neglects the gradients dN/dx, dN/dy of N:
    !
    !   4 N d2u/dx2 + N d2u/dy2 + 3 N d2v/dxdy - beta_b u = -tau_dx
    !   4 N d2v/dy2 + N d2v/dx2 + 3 N d2u/dxdy - beta_b v = -tau_dy
    !
    ! Dividing both sides by N yields:
    !
    !   4 d2u/dx2 + d2u/dy2 + 3 d2v/dxdy - beta_b u / N = -tau_dx / N
    !   4 d2v/dy2 + d2v/dx2 + 3 d2u/dxdy - beta_b v / N = -tau_dy / N
    !
    ! Note that there is no clear mathematical or physical reason why this should be allowed.
    ! However, while I (Tijn Berends, 2023) have found a few cases where there are noticeable
    ! differences in the solutions (e.g. ISMIP-HOM experiments with high strain rates),
    ! most of the time the difference with respect to the full SSA/DIVA is very small.
    ! The "sans" option makes the solver quite a lot more stable and therefore faster.
    ! Someone really ought to perform some proper experiments to determine whether or not
    ! this should be the default.

    ! In/output variables:
    type(type_graph_pair),                                                    intent(in   ) :: graphs
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: N_b_g
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: basal_friction_coefficient_b_g
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: tau_dx_b_g, tau_dy_b_g
    type(type_sparse_matrix_CSR_dp),                                          intent(inout) :: A_CSR
    real(dp), dimension(graphs%graph_b%ni1*2-1: graphs%graph_b%ni2*2),        intent(inout) :: bb_b_g
    integer,                                                                  intent(in   ) :: row_niuv

    ! Local variables:
    integer                             :: ni, uv
    real(dp)                            :: N, basal_friction_coefficient, tau_dx, tau_dy
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_ddx_val
    real(dp), dimension(:), allocatable :: single_row_ddy_val
    real(dp), dimension(:), allocatable :: single_row_d2dx2_val
    real(dp), dimension(:), allocatable :: single_row_d2dxdy_val
    real(dp), dimension(:), allocatable :: single_row_d2dy2_val
    integer                             :: single_row_nnz
    real(dp)                            :: Au, Av
    integer                             :: k, nj, col_nju, col_njv

    ! Relevant indices for this triangle
    ni = graphs%n2biuv( row_niuv,1)
    uv = graphs%n2biuv( row_niuv,2)

    ! N, beta_b, tau_dx, and tau_dy on this graph node
    N                          = N_b_g                         ( ni)
    basal_friction_coefficient = basal_friction_coefficient_b_g( ni)
    tau_dx                     = tau_dx_b_g                    ( ni)
    tau_dy                     = tau_dy_b_g                    ( ni)

    ! allocate memory for single matrix rows
    allocate( single_row_ind(        graphs%graph_a%nC_mem*2))
    allocate( single_row_ddx_val(    graphs%graph_a%nC_mem*2))
    allocate( single_row_ddy_val(    graphs%graph_a%nC_mem*2))
    allocate( single_row_d2dx2_val(  graphs%graph_a%nC_mem*2))
    allocate( single_row_d2dxdy_val( graphs%graph_a%nC_mem*2))
    allocate( single_row_d2dy2_val(  graphs%graph_a%nC_mem*2))

    ! Read coefficients of the operator matrices
    call read_single_row_CSR_dist( graphs%M2_ddx_b_b   , ni, single_row_ind, single_row_ddx_val   , single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_ddy_b_b   , ni, single_row_ind, single_row_ddy_val   , single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_d2dx2_b_b , ni, single_row_ind, single_row_d2dx2_val , single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_d2dxdy_b_b, ni, single_row_ind, single_row_d2dxdy_val, single_row_nnz)
    call read_single_row_CSR_dist( graphs%M2_d2dy2_b_b , ni, single_row_ind, single_row_d2dy2_val , single_row_nnz)

    if (uv == 1) then
      ! x-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        nj      = single_row_ind( k)
        col_nju = graphs%biuv2n( nj,1)
        col_njv = graphs%biuv2n( nj,2)

        !   4 d2u/dx2 + d2u/dy2 + 3 d2v/dxdy - beta_b u / N = -tau_dx / N

        ! Combine the mesh operators
        Au = 4._dp * single_row_d2dx2_val(  k) + &             ! 4 d2u/dx2
                     single_row_d2dy2_val(  k)                 !   d2u/dy2
        if (nj == ni) Au = Au - basal_friction_coefficient / N ! - beta_b u / N

        Av = 3._dp * single_row_d2dxdy_val( k)                 ! 3 d2v/dxdy

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_niuv, col_nju, Au)
        call add_entry_CSR_dist( A_CSR, row_niuv, col_njv, Av)

      end do

      ! Load vector
      bb_b_g( row_niuv) = -tau_dx / N

    elseif (uv == 2) then
      ! y-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        nj      = single_row_ind( k)
        col_nju = graphs%biuv2n( nj,1)
        col_njv = graphs%biuv2n( nj,2)

        !   4 d2v/dy2 + d2v/dx2 + 3 d2u/dxdy - beta_b v / N = -tau_dy / N

        ! Combine the mesh operators
        Av = 4._dp * single_row_d2dy2_val(  k) + &             ! 4 d2v/dy2
                     single_row_d2dx2_val(  k)                 !   d2v/dx2
        if (nj == ni) Av = Av - basal_friction_coefficient / N ! - beta_b v / N

        Au = 3._dp * single_row_d2dxdy_val( k)                 ! 3 d2u/dxdy

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_niuv, col_nju, Au)
        call add_entry_CSR_dist( A_CSR, row_niuv, col_njv, Av)

      end do

      ! Load vector
      bb_b_g( row_niuv) = -tau_dy / N

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_SSA_DIVA_sans_stiffness_matrix_row_free

  subroutine calc_SSA_DIVA_stiffness_matrix_row_BC_ice_front( graphs, N_b_g, Hi_b_ggh, Ho_b_ggh, &
    A_CSR, bb_b_g, row_niuv)
    ! The ocean-pressure boundary condition to the SSA at the ice front reads;
    !
    !   2 N ( 2 du/dx + dv/dy ) n_x + N ( du/dy + dv/dx) n_y = (1/2 rho_i g H^2 - 1/2 rho_sw g Ho^2) n_x
    !
    !   2 N ( 2 dv/dy + du/dx ) n_y + N ( dv/dx + du/dy) n_x = (1/2 rho_i g H^2 - 1/2 rho_sw g Ho^2) n_y
    !
    ! (See also Robinson et al., 2020, Eq. 19)
    !
    ! Rearranging to gather the terms involving u and v gives:
    !
    !   4 N n_x du/dx + N n_y du/dy + ...
    !   2 N n_x dv/dy + N n_y dv/dx = (1/2 rho_i g H^2 - 1/2 rho_sw g Ho^2) n_x
    !
    !   4 N n_y dv/dy + N n_x dv/dx + ...
    !   2 N n_y du/dx + N n_x du/dy = (1/2 rho_i g H^2 - 1/2 rho_sw g Ho^2) n_y
    !
    ! The height of the water column in contact with the ice front is defined as:
    !
    !   Ho = min( max( SL - Hb, 0), rho_i / rho_sw H)

    ! In/output variables:
    type(type_graph_pair),                                                    intent(in   ) :: graphs
    real(dp), dimension(graphs%graph_b%pai%i1_nih:graphs%graph_b%pai%i2_nih), intent(in   ) :: N_b_g, Hi_b_ggh, Ho_b_ggh
    type(type_sparse_matrix_CSR_dp),                                          intent(inout) :: A_CSR
    real(dp), dimension(A_CSR%i1:A_CSR%i2),                                   intent(inout) :: bb_b_g
    integer,                                                                  intent(in   ) :: row_niuv

    ! Local variables:
    integer                             :: ni, uv
    real(dp)                            :: N, Hi, Ho, n_x, n_y
    integer,  dimension(:), allocatable :: single_row_ind
    real(dp), dimension(:), allocatable :: single_row_ddx_val
    real(dp), dimension(:), allocatable :: single_row_ddy_val
    integer                             :: single_row_nnz
    real(dp)                            :: Au, Av
    integer                             :: k, nj, col_nju, col_njv

    ! Relevant indices for this graph node
    ni = graphs%n2biuv( row_niuv,1)
    uv = graphs%n2biuv( row_niuv,2)

    ! N, Hi, Ho, and outward unit normal vector on this graph node
    N   = N_b_g( ni)
    Hi  = Hi_b_ggh( ni)
    Ho  = Ho_b_ggh( ni)
    n_x = graphs%graph_b%ghost_nhat( ni,1)
    n_y = graphs%graph_b%ghost_nhat( ni,2)

    ! Allocate memory for single matrix rows
    allocate( single_row_ind(     graphs%graph_a%nC_mem*2))
    allocate( single_row_ddx_val( graphs%graph_a%nC_mem*2))
    allocate( single_row_ddy_val( graphs%graph_a%nC_mem*2))

    ! Read coefficients of the operator matrices
    call read_single_row_CSR_dist( graphs%M_ddx_ghost_b_b, ni, single_row_ind, single_row_ddx_val, single_row_nnz)
    call read_single_row_CSR_dist( graphs%M_ddy_ghost_b_b, ni, single_row_ind, single_row_ddy_val, single_row_nnz)

    if (uv == 1) then
      ! x-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        nj      = single_row_ind( k)
        col_nju = graphs%biuv2n( nj,1)
        col_njv = graphs%biuv2n( nj,2)

        !   4 N n_x du/dx + N n_y du/dy + ...
        !   2 N n_x dv/dy + N n_y dv/dx = (1/2 rho_i g H^2 - 1/2 rho_sw g Ho^2) n_x

        ! Combine the mesh operators
        Au = 4._dp * N * n_x * single_row_ddx_val( k) + &  ! 4 N n_x du/dx
                     N * n_y * single_row_ddy_val( k)      !   N n_y du/dy

        Av = 2._dp * N * n_x * single_row_ddy_val( k) + &  ! 2 N n_x dv/dy
                     N * n_y * single_row_ddx_val( k)      !   N n_y dv/dx

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_niuv, col_nju, Au)
        call add_entry_CSR_dist( A_CSR, row_niuv, col_njv, Av)

      end do

      ! Load vector: b = (1/2 rho_i g H^2 - 1/2 rho_sw g Ho^2) n_x
      bb_b_g( row_niuv) = (&
          0.5_dp * ice_density      * grav * Hi**2 &
        - 0.5_dp * seawater_density * grav * Ho**2) * n_x

    elseif (uv == 2) then
      ! y-component

      do k = 1, single_row_nnz

        ! Relevant indices for this neighbouring triangle
        nj      = single_row_ind( k)
        col_nju = graphs%biuv2n( nj,1)
        col_njv = graphs%biuv2n( nj,2)

        !   4 N n_y dv/dy + N n_x dv/dx + ...
        !   2 N n_y du/dx + N n_x du/dy = (1/2 rho_i g H^2 - 1/2 rho_sw g Ho^2) n_y

        ! Combine the mesh operators
        Av = 4._dp * N * n_y * single_row_ddy_val( k) + &  ! 4 N n_y dv/dy
                     N * n_x * single_row_ddx_val( k)      !   N n_x dv/dx

        Au = 2._dp * N * n_y * single_row_ddx_val( k) + &  ! 2 N n_y du/dx
                     N * n_x * single_row_ddy_val( k)      !   N n_x du/dy

        ! Add coefficients to the stiffness matrix
        call add_entry_CSR_dist( A_CSR, row_niuv, col_nju, Au)
        call add_entry_CSR_dist( A_CSR, row_niuv, col_njv, Av)

      end do

      ! Load vector: b = (1/2 rho_i g H^2 - 1/2 rho_sw g Ho^2) n_x
      bb_b_g( row_niuv) = (&
          0.5_dp * ice_density      * grav * Hi**2 &
        - 0.5_dp * seawater_density * grav * Ho**2) * n_y

    else
      call crash('uv can only be 1 or 2!')
    end if

  end subroutine calc_SSA_DIVA_stiffness_matrix_row_BC_ice_front

end module solve_linearised_SSA_DIVA_ocean_pressure
