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
    character(len=1024)                   :: parent_mesh_name     !< Name of this graph's parent mesh
    integer                               :: n                    !< Total number of nodes in the graph
    integer                               :: nn                   !< Number of non-ghost nodes in the graph
    integer                               :: ng                   !< Number of ghost nodes in the graph
    integer                               :: nC_mem               !< Maximum allowed number of connections per node

    ! Mapping between parent mesh and graph
    integer,  dimension(:  ), allocatable :: mi2ni                !<     Graph node (ni) to mesh node (mi) translation table
    integer,  dimension(:  ), allocatable :: ni2mi                !<     Graph node (ni) to mesh node (mi) translation table

    ! Node coordinates and connectivity
    real(dp), dimension(:,:), allocatable :: V                    !< [m] Node coordinates
    integer,  dimension(:  ), allocatable :: nC                   !<     Number  of other nodes each node is connected to
    integer,  dimension(:,:), allocatable :: C                    !<     Indices of other nodes each node is connected to

    ! Ghost nodes
    logical,  dimension(:  ), allocatable :: is_ghost             !<     Whether or not a node is a ghost node
    real(dp), dimension(:,:), allocatable :: ghost_nhat           !<     Unit normal vector at each ghost node

  end type type_graph

end module graph_types
