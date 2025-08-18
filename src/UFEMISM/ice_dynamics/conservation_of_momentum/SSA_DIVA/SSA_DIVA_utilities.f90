module SSA_DIVA_utilities

  use mpi_f08, only: MPI_COMM_WORLD, MPI_ALLREDUCE, MPI_DOUBLE_PRECISION, MPI_IN_PLACE, MPI_SUM
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use model_configuration, only: C
  use parameters, only: ice_density, grav
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph_pair
  use ice_model_types, only: type_ice_model
  use mesh_disc_apply_operators, only: map_a_b_2D, ddx_a_b_2D, ddy_a_b_2D, ddx_b_a_2D, ddy_b_a_2D
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared
  use mpi_f08, only: MPI_WIN
  use mesh_graph_mapping, only: map_mesh_vertices_to_graph, map_graph_to_mesh_triangles, &
    map_mesh_triangles_to_graph, map_graph_to_mesh_vertices
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_1D_wrapper

  implicit none

  private

  public :: calc_driving_stress, calc_horizontal_strain_rates, relax_viscosity_iterations, &
    apply_velocity_limits, calc_L2_norm_uv

contains

  subroutine calc_driving_stress( mesh, graphs, ice, tau_dx_b, tau_dy_b)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_graph_pair),                  intent(in   ) :: graphs
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(  out) :: tau_dx_b, tau_dy_b

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'calc_driving_stress'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%BC_ice_front)
    case default
      call crash('unknown BC_ice_front "' // trim( C%BC_ice_front) // '"')
    case ('infinite_slab')
      call calc_driving_stress_infinite_slab( mesh, ice, tau_dx_b, tau_dy_b)
    case ('ocean_pressure')
      call calc_driving_stress_ocean_pressure( mesh, graphs, ice, tau_dx_b, tau_dy_b)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_driving_stress

  subroutine calc_driving_stress_infinite_slab( mesh, ice, tau_dx_b, tau_dy_b)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(  out) :: tau_dx_b, tau_dy_b

    ! Local variables:
    character(len=1024), parameter      :: routine_name = 'calc_driving_stress_infinite_slab'
    real(dp), dimension(:), allocatable :: Hi_b
    real(dp), dimension(:), allocatable :: dHs_dx_b
    real(dp), dimension(:), allocatable :: dHs_dy_b
    integer                             :: ti

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate shared memory
    allocate( Hi_b(     mesh%ti1:mesh%ti2))
    allocate( dHs_dx_b( mesh%ti1:mesh%ti2))
    allocate( dHs_dy_b( mesh%ti1:mesh%ti2))

    ! Calculate Hi, dHs/dx, and dHs/dy on the b-grid
    call map_a_b_2D( mesh, ice%Hi, Hi_b    )
    call ddx_a_b_2D( mesh, ice%Hs, dHs_dx_b)
    call ddy_a_b_2D( mesh, ice%Hs, dHs_dy_b)

    ! Calculate the driving stress
    do ti = mesh%ti1, mesh%ti2
      tau_dx_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dx_b( ti)
      tau_dy_b( ti) = -ice_density * grav * Hi_b( ti) * dHs_dy_b( ti)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_driving_stress_infinite_slab

  subroutine calc_driving_stress_ocean_pressure( mesh, graphs, ice, tau_dx_b, tau_dy_b)

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_graph_pair),                  intent(in   ) :: graphs
    type(type_ice_model),                   intent(in   ) :: ice
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(  out) :: tau_dx_b, tau_dy_b

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'calc_driving_stress_ocean_pressure'
    real(dp), dimension(:), pointer :: Hi_a_g => null()
    real(dp), dimension(:), pointer :: Hs_a_g => null()
    type(MPI_WIN)                   :: wHi_a_g, wHs_a_g
    real(dp), dimension(:), pointer :: Hi_b_g => null()
    real(dp), dimension(:), pointer :: dHs_dx_b_g => null()
    real(dp), dimension(:), pointer :: dHs_dy_b_g => null()
    type(MPI_WIN)                   :: wHi_b_g, wdHs_dx_b_g, wdHs_dy_b_g
    real(dp), dimension(:), pointer :: tau_dx_b_g => null()
    real(dp), dimension(:), pointer :: tau_dy_b_g => null()
    type(MPI_WIN)                   :: wtau_dx_b_g, wtau_dy_b_g
    integer                         :: ni

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( Hi_a_g    , wHi_a_g    , graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( Hs_a_g    , wHs_a_g    , graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( Hi_b_g    , wHi_b_g    , graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( dHs_dx_b_g, wdHs_dx_b_g, graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( dHs_dy_b_g, wdHs_dy_b_g, graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( tau_dx_b_g, wtau_dx_b_g, graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( tau_dy_b_g, wtau_dy_b_g, graphs%graph_b%pai%n_nih)

    ! Map ice model data from the mesh vertices to the a-graph
    call map_mesh_vertices_to_graph( mesh, ice%Hi, graphs%graph_a, Hi_a_g)
    call map_mesh_vertices_to_graph( mesh, ice%Hs, graphs%graph_a, Hs_a_g)

    ! Calculate ice thickness and surface slopes on the b-graph
    call multiply_CSR_matrix_with_vector_1D_wrapper( graphs%M_map_a_b, &
      graphs%graph_a%pai, Hi_a_g, graphs%graph_b%pai, Hi_b_g, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = graphs%graph_a%buffer1_g_nih, buffer_yy_nih = graphs%graph_b%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( graphs%M_ddx_a_b, &
      graphs%graph_a%pai, Hs_a_g, graphs%graph_b%pai, dHs_dx_b_g, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = graphs%graph_a%buffer1_g_nih, buffer_yy_nih = graphs%graph_b%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( graphs%M_ddy_a_b, &
      graphs%graph_a%pai, Hs_a_g, graphs%graph_b%pai, dHs_dy_b_g, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = graphs%graph_a%buffer1_g_nih, buffer_yy_nih = graphs%graph_b%buffer2_g_nih)

    ! Calculate the driving stress
    do ni = graphs%graph_b%ni1, graphs%graph_b%ni2
      tau_dx_b_g( ni) = -ice_density * grav * Hi_b_g( ni) * dHs_dx_b_g( ni)
      tau_dy_b_g( ni) = -ice_density * grav * Hi_b_g( ni) * dHs_dy_b_g( ni)
    end do

    ! Map driving stress from the b-graph to the mesh triangles
    call map_graph_to_mesh_triangles( graphs%graph_b, tau_dx_b_g, mesh, tau_dx_b)
    call map_graph_to_mesh_triangles( graphs%graph_b, tau_dy_b_g, mesh, tau_dy_b)

    ! Clean up after yourself
    call deallocate_dist_shared( Hi_a_g    , wHi_a_g    )
    call deallocate_dist_shared( Hs_a_g    , wHs_a_g    )
    call deallocate_dist_shared( Hi_b_g    , wHi_b_g    )
    call deallocate_dist_shared( dHs_dx_b_g, wdHs_dx_b_g)
    call deallocate_dist_shared( dHs_dy_b_g, wdHs_dy_b_g)
    call deallocate_dist_shared( tau_dx_b_g, wtau_dx_b_g)
    call deallocate_dist_shared( tau_dy_b_g, wtau_dy_b_g)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_driving_stress_ocean_pressure

  subroutine calc_horizontal_strain_rates( mesh, graphs, u_b, v_b, du_dx_a, du_dy_a, dv_dx_a, dv_dy_a)
    !< Calculate the vertically averaged horizontal strain rates

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_graph_pair),                  intent(in   ) :: graphs
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: u_b, v_b
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: du_dx_a, du_dy_a, dv_dx_a, dv_dy_a

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_horizontal_strain_rates'

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%BC_ice_front)
    case default
      call crash('unknown BC_ice_front "' // trim( C%BC_ice_front) // '"')
    case ('infinite_slab')
      call calc_horizontal_strain_rates_infinite_slab( mesh, u_b, v_b, &
        du_dx_a, du_dy_a, dv_dx_a, dv_dy_a)
    case ('ocean_pressure')
      call calc_horizontal_strain_rates_ocean_pressure( mesh, graphs, u_b, v_b, &
        du_dx_a, du_dy_a, dv_dx_a, dv_dy_a)
    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_horizontal_strain_rates

  subroutine calc_horizontal_strain_rates_infinite_slab( mesh, u_b, v_b, &
    du_dx_a, du_dy_a, dv_dx_a, dv_dy_a)
    !< Calculate the vertically averaged horizontal strain rates

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: u_b, v_b
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: du_dx_a, du_dy_a, dv_dx_a, dv_dy_a

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_horizontal_strain_rates_infinite_slab'

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate the strain rates
    call ddx_b_a_2D( mesh, u_b, du_dx_a)
    call ddy_b_a_2D( mesh, u_b, du_dy_a)
    call ddx_b_a_2D( mesh, v_b, dv_dx_a)
    call ddy_b_a_2D( mesh, v_b, dv_dy_a)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_horizontal_strain_rates_infinite_slab

  subroutine calc_horizontal_strain_rates_ocean_pressure( mesh, graphs, u_b, v_b, &
    du_dx_a, du_dy_a, dv_dx_a, dv_dy_a)
    !< Calculate the vertically averaged horizontal strain rates

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    type(type_graph_pair),                  intent(in   ) :: graphs
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: u_b, v_b
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(  out) :: du_dx_a, du_dy_a, dv_dx_a, dv_dy_a

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_horizontal_strain_rates_ocean_pressure'
    real(dp), dimension(:), pointer :: u_b_g => null()
    real(dp), dimension(:), pointer :: v_b_g => null()
    type(MPI_WIN)                   :: wu_b_g, wv_b_g
    real(dp), dimension(:), pointer :: du_dx_a_g => null()
    real(dp), dimension(:), pointer :: du_dy_a_g => null()
    real(dp), dimension(:), pointer :: dv_dx_a_g => null()
    real(dp), dimension(:), pointer :: dv_dy_a_g => null()
    type(MPI_WIN)                   :: wdu_dx_a_g, wdu_dy_a_g, wdv_dx_a_g, wdv_dy_a_g

    ! Add routine to path
    call init_routine( routine_name)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( u_b_g    , wu_b_g    , graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( v_b_g    , wv_b_g    , graphs%graph_b%pai%n_nih)
    call allocate_dist_shared( du_dx_a_g, wdu_dx_a_g, graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( du_dy_a_g, wdu_dy_a_g, graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( dv_dx_a_g, wdv_dx_a_g, graphs%graph_a%pai%n_nih)
    call allocate_dist_shared( dv_dy_a_g, wdv_dy_a_g, graphs%graph_a%pai%n_nih)

    ! Map ice model data from the mesh triangles to the b-graph
    call map_mesh_triangles_to_graph( mesh, u_b, graphs%graph_b, u_b_g)
    call map_mesh_triangles_to_graph( mesh, v_b, graphs%graph_b, v_b_g)

    ! ! Calculate the strain rates
    call multiply_CSR_matrix_with_vector_1D_wrapper( graphs%M_ddx_b_a, &
      graphs%graph_b%pai, u_b_g, graphs%graph_a%pai, du_dx_a_g, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = graphs%graph_b%buffer1_g_nih, buffer_yy_nih = graphs%graph_a%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( graphs%M_ddy_b_a, &
      graphs%graph_b%pai, u_b_g, graphs%graph_a%pai, du_dy_a_g, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = graphs%graph_b%buffer1_g_nih, buffer_yy_nih = graphs%graph_a%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( graphs%M_ddx_b_a, &
      graphs%graph_b%pai, v_b_g, graphs%graph_a%pai, dv_dx_a_g, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = graphs%graph_b%buffer1_g_nih, buffer_yy_nih = graphs%graph_a%buffer2_g_nih)

    call multiply_CSR_matrix_with_vector_1D_wrapper( graphs%M_ddy_b_a, &
      graphs%graph_b%pai, v_b_g, graphs%graph_a%pai, dv_dy_a_g, &
      xx_is_hybrid = .true., yy_is_hybrid = .true., &
      buffer_xx_nih = graphs%graph_b%buffer1_g_nih, buffer_yy_nih = graphs%graph_a%buffer2_g_nih)

    ! Map strain rates from the a-graph to the mesh vertices
    call map_graph_to_mesh_vertices( graphs%graph_a, du_dx_a_g, mesh, du_dx_a)
    call map_graph_to_mesh_vertices( graphs%graph_a, du_dy_a_g, mesh, du_dy_a)
    call map_graph_to_mesh_vertices( graphs%graph_a, dv_dx_a_g, mesh, dv_dx_a)
    call map_graph_to_mesh_vertices( graphs%graph_a, dv_dy_a_g, mesh, dv_dy_a)

    ! Clean up after yourself
    call deallocate_dist_shared( u_b_g    , wu_b_g    )
    call deallocate_dist_shared( v_b_g    , wv_b_g    )
    call deallocate_dist_shared( du_dx_a_g, wdu_dx_a_g)
    call deallocate_dist_shared( du_dy_a_g, wdu_dy_a_g)
    call deallocate_dist_shared( dv_dx_a_g, wdv_dx_a_g)
    call deallocate_dist_shared( dv_dy_a_g, wdv_dy_a_g)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_horizontal_strain_rates_ocean_pressure

  subroutine relax_viscosity_iterations( mesh, u_b, v_b, u_b_prev, v_b_prev, visc_it_relax)
    !< Reduce the change between velocity solutions

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: u_b, v_b
    real(dp), dimension(mesh%nTri),         intent(in   ) :: u_b_prev, v_b_prev
    real(dp),                               intent(in   ) :: visc_it_relax

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'relax_viscosity_iterations'
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2
      u_b( ti) = (visc_it_relax * u_b( ti)) + ((1._dp - visc_it_relax) * u_b_prev( ti))
      v_b( ti) = (visc_it_relax * v_b( ti)) + ((1._dp - visc_it_relax) * v_b_prev( ti))
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine relax_viscosity_iterations

  subroutine apply_velocity_limits( mesh, u_b, v_b)
    !< Limit velocities for improved stability

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(inout) :: u_b, v_b

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'apply_velocity_limits'
    integer                        :: ti
    real(dp)                       :: uabs

    ! Add routine to path
    call init_routine( routine_name)

    do ti = mesh%ti1, mesh%ti2

      ! Calculate absolute speed
      uabs = sqrt( u_b( ti)**2 + v_b( ti)**2)

      ! Reduce velocities if necessary
      if (uabs > C%vel_max) then
        u_b( ti) = u_b( ti) * C%vel_max / uabs
        v_b( ti) = v_b( ti) * C%vel_max / uabs
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_velocity_limits

  subroutine calc_L2_norm_uv( mesh, u_b, v_b, u_b_prev, v_b_prev, L2_uv)
    !< Calculate the L2-norm of the two consecutive velocity solutions

    ! In/output variables:
    type(type_mesh),                        intent(in   ) :: mesh
    real(dp), dimension(mesh%ti1:mesh%ti2), intent(in   ) :: u_b, v_b
    real(dp), dimension(mesh%nTri),         intent(in   ) :: u_b_prev, v_b_prev
    real(dp),                               intent(  out) :: L2_uv

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_visc_iter_UV_resid'
    integer                        :: ierr
    integer                        :: ti
    real(dp)                       :: res1, res2

    ! Add routine to path
    call init_routine( routine_name)

    res1 = 0._dp
    res2 = 0._dp

    do ti = mesh%ti1, mesh%ti2

      res1 = res1 + (u_b( ti) - u_b_prev( ti))**2
      res1 = res1 + (v_b( ti) - v_b_prev( ti))**2

      res2 = res2 + (u_b( ti) + u_b_prev( ti))**2
      res2 = res2 + (v_b( ti) + v_b_prev( ti))**2

    end do

    ! Combine results from all processes
    call MPI_ALLREDUCE( MPI_IN_PLACE, res1, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
    call MPI_ALLREDUCE( MPI_IN_PLACE, res2, 1, MPI_doUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

    ! Calculate L2-norm
    L2_uv = 2._dp * res1 / max( res2, 1E-8_dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_L2_norm_uv

end module SSA_DIVA_utilities
