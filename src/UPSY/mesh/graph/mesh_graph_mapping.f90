module mesh_graph_mapping

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph
  use mpi_f08, only: MPI_WIN
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use deallocate_dist_shared_mod, only: deallocate_dist_shared
  use mpi_distributed_memory, only: gather_to_all
  use mpi_distributed_shared_memory, only: gather_dist_shared_to_all
  use mpi_basic, only: par, sync

  implicit none

  private

  public :: map_mesh_vertices_to_graph

  interface map_mesh_vertices_to_graph
    procedure :: map_mesh_vertices_to_graph_logical_2D
    procedure :: map_mesh_vertices_to_graph_logical_3D
    procedure :: map_mesh_vertices_to_graph_int_2D
    procedure :: map_mesh_vertices_to_graph_int_3D
    procedure :: map_mesh_vertices_to_graph_dp_2D
    procedure :: map_mesh_vertices_to_graph_dp_3D
  end interface map_mesh_vertices_to_graph

contains

  subroutine map_mesh_vertices_to_graph_logical_2D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                                       intent(in   ) :: mesh
    logical, dimension(:),                                 intent(in   ) :: d_mesh
    type(type_graph),                                      intent(in   ) :: graph
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_vertices_to_graph_logical_2D'
    logical, dimension(:), pointer :: d_mesh_tot => null()
    type(MPI_WIN)                  :: wd_mesh_tot
    integer                        :: ni, vi

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_mesh_tot, wd_mesh_tot, mesh%nV)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_V%n_loc) then
      call gather_to_all( d_mesh, d_mesh_tot)
    elseif (size( d_mesh,1) == mesh%pai_V%n_nih) then
      call gather_dist_shared_to_all( mesh%pai_V, d_mesh, d_mesh_tot)
    else
      call crash('invalid size for d_mesh')
    end if

    do ni = graph%pai%i1, graph%pai%i2
      vi = graph%ni2mi( ni)
      d_graph_nih( ni) = d_mesh_tot( vi)
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_tot, wd_mesh_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_vertices_to_graph_logical_2D

  subroutine map_mesh_vertices_to_graph_logical_3D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    logical, dimension(:,:), target, intent(in   ) :: d_mesh
    type(type_graph),                intent(in   ) :: graph
    logical, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_vertices_to_graph_logical_3D'
    integer                        :: k
    logical, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_mesh_vertices_to_graph_logical_2D( mesh, d_mesh_k, graph, d_graph_nih_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_vertices_to_graph_logical_3D

  subroutine map_mesh_vertices_to_graph_int_2D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                                       intent(in   ) :: mesh
    integer, dimension(:),                                 intent(in   ) :: d_mesh
    type(type_graph),                                      intent(in   ) :: graph
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_vertices_to_graph_int_2D'
    integer, dimension(:), pointer :: d_mesh_tot => null()
    type(MPI_WIN)                  :: wd_mesh_tot
    integer                        :: ni, vi

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_mesh_tot, wd_mesh_tot, mesh%nV)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_V%n_loc) then
      call gather_to_all( d_mesh, d_mesh_tot)
    elseif (size( d_mesh,1) == mesh%pai_V%n_nih) then
      call gather_dist_shared_to_all( mesh%pai_V, d_mesh, d_mesh_tot)
    else
      call crash('invalid size for d_mesh')
    end if

    do ni = graph%pai%i1, graph%pai%i2
      vi = graph%ni2mi( ni)
      d_graph_nih( ni) = d_mesh_tot( vi)
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_tot, wd_mesh_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_vertices_to_graph_int_2D

  subroutine map_mesh_vertices_to_graph_int_3D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    integer, dimension(:,:), target, intent(in   ) :: d_mesh
    type(type_graph),                intent(in   ) :: graph
    integer, dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),     target, intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_mesh_vertices_to_graph_int_3D'
    integer                        :: k
    integer, dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_mesh_vertices_to_graph_int_2D( mesh, d_mesh_k, graph, d_graph_nih_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_vertices_to_graph_int_3D

  subroutine map_mesh_vertices_to_graph_dp_2D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                                        intent(in   ) :: mesh
    real(dp), dimension(:),                                 intent(in   ) :: d_mesh
    type(type_graph),                                       intent(in   ) :: graph
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih), intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_mesh_vertices_to_graph_dp_2D'
    real(dp), dimension(:), pointer :: d_mesh_tot => null()
    type(MPI_WIN)                   :: wd_mesh_tot
    integer                         :: ni, vi

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_dist_shared( d_mesh_tot, wd_mesh_tot, mesh%nV)

    ! Check if the mesh data is distributed or hybrid distributed/shared
    if (size( d_mesh,1) == mesh%pai_V%n_loc) then
      call gather_to_all( d_mesh, d_mesh_tot)
    elseif (size( d_mesh,1) == mesh%pai_V%n_nih) then
      call gather_dist_shared_to_all( mesh%pai_V, d_mesh, d_mesh_tot)
    else
      call crash('invalid size for d_mesh')
    end if

    do ni = graph%pai%i1, graph%pai%i2
      vi = graph%ni2mi( ni)
      d_graph_nih( ni) = d_mesh_tot( vi)
    end do

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_tot, wd_mesh_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_vertices_to_graph_dp_2D

  subroutine map_mesh_vertices_to_graph_dp_3D( mesh, d_mesh, graph, d_graph_nih)

    ! In/output variables:
    type(type_mesh),                  intent(in   ) :: mesh
    real(dp), dimension(:,:), target, intent(in   ) :: d_mesh
    type(type_graph),                 intent(in   ) :: graph
    real(dp), dimension(graph%pai%i1_nih:graph%pai%i2_nih, &
      1:size(d_mesh,2)),      target, intent(  out) :: d_graph_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'map_mesh_vertices_to_graph_dp_3D'
    integer                         :: k
    real(dp), dimension(:), pointer :: d_mesh_k, d_graph_nih_k

    ! Add routine to path
    call init_routine( routine_name)

    do k = 1, size( d_mesh,2)
      d_mesh_k      => d_mesh     ( :,k)
      d_graph_nih_k => d_graph_nih( :,k)
      call map_mesh_vertices_to_graph_dp_2D( mesh, d_mesh_k, graph, d_graph_nih_k)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine map_mesh_vertices_to_graph_dp_3D

end module mesh_graph_mapping