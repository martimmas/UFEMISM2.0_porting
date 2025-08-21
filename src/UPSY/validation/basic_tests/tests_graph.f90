module tests_graph

  ! Basic tests for graphs.

  use precisions, only: dp
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mpi_basic, only: par

  implicit none

  private

  public :: test_graph_is_self_consistent

contains

  function test_graph_is_self_consistent( mesh, graph) result( isso)

    ! In/output variables:
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_is_self_consistent'

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.

    isso = isso .and. test_graph_matches_mesh( mesh, graph)
    isso = isso .and. test_graph_connectivity_is_self_consistent( graph)
    isso = isso .and. test_graph_to_mesh_connectivity_is_self_consistent( graph)
    isso = isso .and. test_graph_nodes_are_sorted( graph)
    isso = isso .and. test_graph_ghost_nodes( mesh, graph)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_is_self_consistent

  function test_graph_matches_mesh( mesh, graph) result( isso)

    ! In/output variables:
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_matches_mesh'

    ! Add routine to path
    call init_routine( routine_name)

    if (graph%parent_mesh_name == trim( mesh%name) // '_vertices') then
      isso = test_graph_matches_mesh_vertices( mesh, graph)
    elseif (graph%parent_mesh_name == trim( mesh%name) // '_triangles') then
      isso = test_graph_matches_mesh_triangles( mesh, graph)
    else
      call crash('mesh is not this graphs parent mesh')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_matches_mesh

  function test_graph_matches_mesh_vertices( mesh, graph) result( isso)

    ! In/output variables:
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_matches_mesh'
    integer                        :: ni, vi

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.
    do ni = 1, graph%n
      if (.not. graph%is_ghost( ni)) then
        vi = graph%ni2mi( ni)
        isso = isso .and. norm2( graph%V( ni,:) - mesh%V( vi,:)) < mesh%tol_dist
      else
        ! An a-graph should not have any ghost nodes at all
        isso = .false.
      end if
    end do

    if (par%primary .and. .not. isso) call warning('graph node coordinates do not match corresponding mesh vertices')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_matches_mesh_vertices

  function test_graph_matches_mesh_triangles( mesh, graph) result( isso)

    ! In/output variables:
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_matches_mesh_triangles'
    integer                        :: ni, ti

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.
    do ni = 1, graph%n
      if (.not. graph%is_ghost( ni)) then
        ti = graph%ni2mi( ni)
        isso = isso .and. norm2( graph%V( ni,:) - mesh%TriGC( ti,:)) < mesh%tol_dist
      end if
    end do

    if (par%primary .and. .not. isso) call warning('graph node coordinates do not match corresponding mesh triangles')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_matches_mesh_triangles

  function test_graph_connectivity_is_self_consistent( graph) result( isso)

    ! In/output variables:
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_connectivity_is_self_consistent'
    integer                        :: ni, ci, nj, cj, nk
    logical                        :: found_reverse

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.

    do ni = 1, graph%n
      do ci = 1, graph%nC( ni)
        nj = graph%C( ni,ci)
        found_reverse = .false.
        do cj = 1, graph%nC( nj)
          nk = graph%C( nj,cj)
          if (nk == ni) then
            found_reverse = .true.
            exit
          end if
        end do
        isso = isso .and. found_reverse
      end do
    end do

    if (par%primary .and. .not. isso) call warning('graph connectivity is not self-consistent')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_connectivity_is_self_consistent

  function test_graph_to_mesh_connectivity_is_self_consistent( graph) result( isso)

    ! In/output variables:
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_to_mesh_connectivity_is_self_consistent'
    integer                        :: ni, ci, nj, cj, nk, mi
    logical                        :: found_reverse

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.

    do mi = 1, size( graph%mi2ni)
      ni = graph%mi2ni( mi)
      if (ni > 0) isso = isso .and. graph%ni2mi( ni) == mi
    end do

    do ni = 1, graph%n
      if (.not. graph%is_ghost( ni)) then
        mi = graph%ni2mi( ni)
        isso = isso .and. graph%mi2ni( mi) == ni
      end if
    end do

    if (par%primary .and. .not. isso) call warning('graph-to-mesh connectivity is not self-consistent')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_to_mesh_connectivity_is_self_consistent

  function test_graph_ghost_nodes( mesh, graph) result( isso)

    ! In/output variables:
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_ghost_nodes'

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.
    isso = isso .and. test_graph_ghost_nodes_connectivity( mesh, graph)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_ghost_nodes

  function test_graph_ghost_nodes_connectivity( mesh, graph) result( isso)

    ! In/output variables:
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_ghost_nodes_connectivity'
    integer                        :: ni

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.
    do ni = 1, graph%n
      if (graph%is_ghost( ni)) then
        isso = isso .and. graph%nC( ni) == 1
      end if
    end do

    if (par%primary .and. .not. isso) call warning('graph ghost nodes do not all have only a single connection')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_ghost_nodes_connectivity

  function test_graph_nodes_are_sorted( graph) result( isso)

    ! In/output variables:
    type(type_graph), intent(in) :: graph
    logical                      :: isso

    ! Local variables:
    character(len=256), parameter :: routine_name = 'test_graph_nodes_are_sorted'
    integer                       :: ni

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.
    do ni = 2, graph%n
      isso = isso .and. graph%V( ni,1) >= graph%V( ni-1,1)
    end do

    if (par%primary .and. .not. isso) call warning('graph nodes are not sorted in the x-dimension')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_nodes_are_sorted

end module tests_graph
