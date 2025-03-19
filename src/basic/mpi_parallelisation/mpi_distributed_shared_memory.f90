module mpi_distributed_shared_memory

  use allocate_dist_shared_mod, only: allocate_dist_shared
  use deallocate_dist_shared_mod, only: deallocate_dist_shared
  use gather_dist_shared_to_primary_mod, only: gather_dist_shared_to_primary
  use gather_dist_shared_to_all_mod, only: gather_dist_shared_to_all
  use distribute_dist_shared_from_primary_mod, only: distribute_dist_shared_from_primary

  implicit none

  private

  public :: allocate_dist_shared, deallocate_dist_shared, gather_dist_shared_to_primary, &
    gather_dist_shared_to_all, distribute_dist_shared_from_primary

contains

end module mpi_distributed_shared_memory
