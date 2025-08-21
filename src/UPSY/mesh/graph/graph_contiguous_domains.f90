module graph_contiguous_domains

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use graph_types, only: type_graph
  use sorting, only: quick_n_dirty_sort

  implicit none

  private

  public :: enforce_contiguous_process_domains_graph

contains

  !> Shuffle nodes so that the nodes owned by each process form
  !> a contiguous domain with short borders, to minimise the data volumes for halo exchanges.
  subroutine enforce_contiguous_process_domains_graph( graph)

    ! In/output variables:
    type(type_graph), intent(inout) :: graph

    ! Local variables:
    character(len=256), parameter              :: routine_name = 'enforce_contiguous_process_domains_graph'
    real(dp), dimension(graph%n)               :: xx
    integer,  dimension(graph%n)               :: ni_new2ni_old, ni_old2ni_new
    integer                                    :: ni_old, ni_new, ci, nj_old, nj_new, ni, mi
    real(dp), dimension(graph%n,2)             :: V_old
    integer,  dimension(graph%n)               :: nC_old
    integer,  dimension(graph%n, graph%nC_mem) :: C_old
    integer,  dimension(graph%n)               :: ni_old2mi
    integer,  dimension(size(graph%mi2ni))     :: mi2ni_old
    logical,  dimension(graph%n)               :: is_ghost_old
    real(dp), dimension(graph%n,2)             :: V_ghost_BC_old
    real(dp), dimension(graph%n,2)             :: ghost_nhat_old

    ! Add routine to path
    call init_routine( routine_name)

    ! Sort nodes by x-coordinate
    xx = graph%V( :,1)
    call quick_n_dirty_sort( xx, ni_new2ni_old)

    ! Calculate translation table in the opposite direction
    do ni_new = 1, graph%n
      ni_old = ni_new2ni_old( ni_new)
      ni_old2ni_new( ni_old) = ni_new
    end do

    ! Shuffle node data: V, nC, C
    ! ===========================

    ni_old2mi      = graph%ni2mi
    mi2ni_old      = graph%mi2ni
    V_old          = graph%V
    nC_old         = graph%nC
    C_old          = graph%C
    is_ghost_old   = graph%is_ghost
    V_ghost_BC_old = graph%V_ghost_BC
    ghost_nhat_old = graph%ghost_nhat

    graph%ni2mi      = 0
    graph%mi2ni      = 0
    graph%V          = 0._dp
    graph%nC         = 0
    graph%C          = 0
    graph%is_ghost   = .false.
    graph%V_ghost_BC = 0._dp
    graph%ghost_nhat = 0._dp

    do ni_new = 1, graph%n

      ! This new node corresponds to this old node
      ni_old = ni_new2ni_old( ni_new)

      ! ni2mi, mi2ni
      mi = ni_old2mi( ni_old)
      graph%ni2mi( ni_new) = mi
      if (.not. is_ghost_old( ni_old)) graph%mi2ni( mi) = ni_new

      ! V
      graph%V( ni_new,:) = V_old( ni_old,:)

      ! nC
      graph%nC( ni_new) = nC_old( ni_old)

      ! C
      do ci = 1, nC_old( ni_old)
        nj_old = C_old( ni_old, ci)
        nj_new = ni_old2ni_new( nj_old)
        graph%C( ni_new,ci) = nj_new
      end do

      ! is_ghost
      graph%is_ghost( ni_new) = is_ghost_old( ni_old)

      ! V_ghost_BC
      graph%V_ghost_BC( ni_new,:) = V_ghost_BC_old( ni_old,:)

      ! ghost_nhat
      graph%ghost_nhat( ni_new,:) = ghost_nhat_old( ni_old,:)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine enforce_contiguous_process_domains_graph

end module graph_contiguous_domains