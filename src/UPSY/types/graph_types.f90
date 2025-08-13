module graph_types

  use precisions, only: dp
  use CSR_sparse_matrix_type, only: type_sparse_matrix_CSR_dp
  use parallel_array_info_type, only: type_par_arr_info
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_graph

  type type_graph

    ! Metadata
    character(len=1024)                     :: parent_mesh_name              !< Name of this graph's parent mesh
    integer                                 :: n                             !< Total number of nodes in the graph
    integer                                 :: nn                            !< Number of non-ghost nodes in the graph
    integer                                 :: ng                            !< Number of ghost nodes in the graph
    integer                                 :: nC_mem                        !< Maximum allowed number of connections per node

    ! Mapping between parent mesh and graph
    integer,  dimension(:    ), allocatable :: mi2ni                         !<     Graph node (ni) to mesh node (mi) translation table
    integer,  dimension(:    ), allocatable :: ni2mi                         !<     Graph node (ni) to mesh node (mi) translation table

    ! Node coordinates and connectivity
    real(dp), dimension(:,:  ), allocatable :: V                             !< [m] Node coordinates
    integer,  dimension(:    ), allocatable :: nC                            !<     Number  of other nodes each node is connected to
    integer,  dimension(:,:  ), allocatable :: C                             !<     Indices of other nodes each node is connected to

    ! Ghost nodes
    logical,  dimension(:    ), allocatable :: is_ghost                      !<     Whether or not a node is a ghost node
    real(dp), dimension(:,:  ), allocatable :: ghost_nhat                    !<     Unit normal vector at each ghost node

    ! Parallelisation ranges
    integer                                 :: ni1, ni2, n_loc               !< Each process "owns" n_loc nodes ni1 - ni2, so that n_loc = ni2 + 1 - ni1

    integer,  dimension(:    ), allocatable :: owning_process                !< Which process owns each node
    integer,  dimension(:    ), allocatable :: owning_node                   !< Which HPC node owns each graph node (yes, this is confusing...)

    type(type_par_arr_info)                 :: pai                           !< Parallelisation info for graph-based fields

    real(dp), dimension(:    ), contiguous, pointer :: buffer1_d_g_nih  => null()     !< Pre-allocated buffer memory on the graph
    real(dp), dimension(:    ), contiguous, pointer :: buffer2_d_g_nih  => null()     !< Pre-allocated buffer memory on the graph
    ! real(dp), dimension(:,:  ), contiguous, pointer :: buffer1_d_gk_nih => null()     !< Pre-allocated buffer memory on the graph
    ! real(dp), dimension(:,:  ), contiguous, pointer :: buffer2_d_gk_nih => null()     !< Pre-allocated buffer memory on the graph
    type(MPI_WIN) :: wbuffer1_d_g_nih, wbuffer2_d_g_nih     !< MPI window to pre-allocated buffer memory on the graph
    ! type(MPI_WIN) :: wbuffer1_d_g_nih, wbuffer2_d_g_nih, wbuffer1_d_gk_nih, wbuffer2_d_gk_nih     !< MPI window to pre-allocated buffer memory on the graph

  end type type_graph

end module graph_types