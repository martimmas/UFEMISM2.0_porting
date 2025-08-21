module graph_memory

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use graph_types, only: type_graph
  use reallocate_mod, only: reallocate
  use mpi_distributed_shared_memory, only: deallocate_dist_shared

  implicit none

  private

  public :: allocate_graph_primary, crop_graph_primary, deallocate_graph

contains

  subroutine allocate_graph_primary( graph, n, n_mesh, nC_mem)

    ! In/output variables:
    type(type_graph), intent(inout) :: graph
    integer,          intent(in   ) :: n, n_mesh, nC_mem

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'allocate_graph_primary'

    ! Add routine to path
    call init_routine( routine_name)

    graph%n      = 0
    graph%nC_mem = nC_mem

    ! Mapping between parent mesh and graph
    allocate( graph%ni2mi     ( n        ), source = 0)
    allocate( graph%mi2ni     ( n_mesh   ), source = 0)

    ! Node coordinates and connectivity
    allocate( graph%V         ( n, 2     ), source = 0._dp)
    allocate( graph%nC        ( n        ), source = 0)
    allocate( graph%C         ( n, nC_mem), source = 0)

    ! Ghost nodes
    allocate( graph%is_ghost  ( n        ), source = .false.)
    allocate( graph%V_ghost_BC( n, 2     ), source = 0._dp)
    allocate( graph%ghost_nhat( n, 2     ), source = 0._dp)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine allocate_graph_primary

  subroutine crop_graph_primary( graph)

    ! In/output variables:
    type(type_graph), intent(inout) :: graph

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'crop_graph_primary'

    ! Add routine to path
    call init_routine( routine_name)

    ! Mapping between parent mesh and graph
    call reallocate( graph%ni2mi, graph%n)

    ! Node coordinates and connectivity
    call reallocate( graph%V         , graph%n, 2           )
    call reallocate( graph%nC        , graph%n              )
    call reallocate( graph%C         , graph%n, graph%nC_mem)

    ! Ghost nodes
    call reallocate( graph%is_ghost  , graph%n              )
    call reallocate( graph%V_ghost_BC, graph%n, 2           )
    call reallocate( graph%ghost_nhat, graph%n, 2           )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine crop_graph_primary

  subroutine deallocate_graph( graph)

    ! In/output variables:
    type(type_graph), intent(inout) :: graph

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'deallocate_graph'

    ! Add routine to path
    call init_routine( routine_name)

    ! Mapping between parent mesh and graph
    deallocate( graph%ni2mi)
    deallocate( graph%mi2ni)

    ! Node coordinates and connectivity
    deallocate( graph%V )
    deallocate( graph%nC)
    deallocate( graph%C )

    ! Ghost nodes
    deallocate( graph%is_ghost  )
    deallocate( graph%V_ghost_BC)
    deallocate( graph%ghost_nhat)

    ! Parallelisation ranges
    deallocate( graph%owning_process)
    deallocate( graph%owning_node   )

    call deallocate_dist_shared( graph%buffer1_g_nih, graph%wbuffer1_g_nih)
    call deallocate_dist_shared( graph%buffer2_g_nih, graph%wbuffer2_g_nih)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_graph

end module graph_memory