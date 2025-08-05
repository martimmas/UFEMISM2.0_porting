module delete_vertices

  ! Delete a vertex from the mesh and update the Delaunay triangulation accordingly.

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use assertions_basic, only: assert
  use switch_array_elements, only: switch_rows
  use switch_vertices_triangles, only: switch_vertices, switch_triangles
  use mesh_utilities, only: update_triangle_circumcenter
  use flip_triangles, only: add_triangle_pairs_around_triangle_to_Delaunay_check_stack, &
    flip_triangles_until_Delaunay
  use plane_geometry, only: cross2

  implicit none

  private

  public :: delete_vertex

contains

  subroutine delete_vertex( mesh, vi_kill, &
    vi_new2vi_old, vi_old2vi_new, ti_new2ti_old, ti_old2ti_new)
    ! Delete vertex vi from the mesh

    ! In/output variables:
    type(type_mesh),                    intent(inout) :: mesh
    integer,                            intent(in   ) :: vi_kill
    integer, dimension(:), allocatable, intent(  out) :: vi_new2vi_old, vi_old2vi_new
    integer, dimension(:), allocatable, intent(  out) :: ti_new2ti_old, ti_old2ti_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex'

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (mesh%VBI( vi_kill) > 0) call crash('delete_vertex only works on free vertices for now')

    if (mesh%nC( vi_kill) == 3) then
      call crash('delete_vertex_nCeq3 not implemented yet')
      ! call delete_vertex_nCeq3( mesh, vi_kill, &
      !   vi_new2vi_old, vi_old2vi_new, ti_new2ti_old, ti_old2ti_new)
    else
      call delete_vertex_nCge4( mesh, vi_kill, &
        vi_new2vi_old, vi_old2vi_new, ti_new2ti_old, ti_old2ti_new)
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex

  ! Delete a vertex which has 4 or more neighbours
  ! ==============================================

  subroutine delete_vertex_nCge4( mesh, vi_kill, &
    vi_new2vi_old, vi_old2vi_new, ti_new2ti_old, ti_old2ti_new)
    ! Before, local geometry looks like this:
    !
    !                tj_opp()
    !
    !  vj_opp(nvj_opp) ------ vj_opp(1)
    !           /  \ ti_opp() /  \
    !     tj1  /    \        /    \  tj4
    !         / ti1  \      /  ti4 \
    !        /        \    /        \
    !   vj_clock ---- vi_kill ---- vj_anti
    !         \          |          /
    !           \  ti2   |   ti3  /
    !             \      |      /
    !        tj2    \    |    /    tj3
    !                 vj_focus
    !
    ! After, it looks like this:
    !
    !                tj_opp()
    !
    !  vj_opp(nvj_opp) ------ vj_opp(1)
    !           /  \ ti_opp() /  \
    !     tj1  /    |        |    \  tj4
    !         /     \        /     \
    !        /  ti2  |      |  ti3  \
    !   vj_clock     \      /      vj_anti
    !         \       |    |        /
    !           \     \    /      /
    !             \    |  |     /
    !        tj2    \   \ /   /    tj3
    !                 vj_focus

    ! In/output variables:
    type(type_mesh),                    intent(inout) :: mesh
    integer,                            intent(in   ) :: vi_kill
    integer, dimension(:), allocatable, intent(  out) :: vi_new2vi_old, vi_old2vi_new
    integer, dimension(:), allocatable, intent(  out) :: ti_new2ti_old, ti_old2ti_new

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'delete_vertex_nCge4'
    integer                            :: vi, ti
    integer                            :: vj_clock, vj_focus, vj_anti
    integer, dimension(:), allocatable :: vj_opp
    integer                            :: ti1, ti2, ti3, ti4
    integer, dimension(:), allocatable :: ti_opp
    integer                            :: tj1, tj2, tj3, tj4
    integer, dimension(:), allocatable :: tj_opp
    integer                            :: vii, tii

    ! Add routine to path
    call init_routine( routine_name)

    allocate( vi_new2vi_old( mesh%nV))
    allocate( vi_old2vi_new( mesh%nV))
    do vi = 1, mesh%nV
      vi_new2vi_old( vi) = vi
      vi_old2vi_new( vi) = vi
    end do

    allocate( ti_new2ti_old( mesh%nTri))
    allocate( ti_old2ti_new( mesh%nTri))
    do ti = 1, mesh%nTri
      ti_new2ti_old( ti) = ti
      ti_old2ti_new( ti) = ti
    end do

    ! Determine the local geometry
    call delete_vertex_nCge4_local_geometry( mesh, vi_kill, &
      vj_clock, vj_focus, vj_anti, vj_opp, &
      ti1, ti2, ti3, ti4, ti_opp, tj1, tj2, tj3, tj4, tj_opp)

    call delete_vertex_nCge4_nC_C(       mesh, vi_kill, vj_clock, vj_focus, vj_anti, vj_opp, ti1, ti2, ti3, ti4, ti_opp, tj1, tj2, tj3, tj4, tj_opp)
    call delete_vertex_nCge4_niTri_iTri( mesh, vi_kill, vj_clock, vj_focus, vj_anti, vj_opp, ti1, ti2, ti3, ti4, ti_opp, tj1, tj2, tj3, tj4, tj_opp)
    call delete_vertex_nCge4_Tri(        mesh, vi_kill, vj_clock, vj_focus, vj_anti, vj_opp, ti1, ti2, ti3, ti4, ti_opp, tj1, tj2, tj3, tj4, tj_opp)
    call delete_vertex_nCge4_TriC(       mesh, vi_kill, vj_clock, vj_focus, vj_anti, vj_opp, ti1, ti2, ti3, ti4, ti_opp, tj1, tj2, tj3, tj4, tj_opp)

    call delete_vertex_V  ( mesh, vi_kill, vi_new2vi_old, vi_old2vi_new)
    call delete_vertex_Tri( mesh, ti1    , ti_new2ti_old, ti_old2ti_new)
    call delete_vertex_Tri( mesh, ti4    , ti_new2ti_old, ti_old2ti_new)

    vj_clock = vi_old2vi_new( vj_clock)
    vj_focus = vi_old2vi_new( vj_focus)
    vj_anti  = vi_old2vi_new( vj_anti )
    do vii = 1, size( vj_opp,1)
      vj_opp( vii) = vi_old2vi_new( vj_opp( vii))
    end do

    ti2 = ti_old2ti_new( ti2)
    ti3 = ti_old2ti_new( ti3)
    do tii = 1, size( ti_opp,1)
      ti_opp( tii) = ti_old2ti_new( ti_opp( tii))
    end do

    call update_triangle_circumcenter( mesh, ti2)
    call update_triangle_circumcenter( mesh, ti3)
    do tii = 1, size( ti_opp,1)
      call update_triangle_circumcenter( mesh, ti_opp( tii))
    end do

    mesh%check_Delaunay_map    = .false.
    mesh%check_Delaunay_stack  = 0
    mesh%check_Delaunay_stackN = 0

    call add_triangle_pairs_around_triangle_to_Delaunay_check_stack( mesh, ti2)
    call add_triangle_pairs_around_triangle_to_Delaunay_check_stack( mesh, ti3)
    do tii = 1, size( ti_opp,1)
      call add_triangle_pairs_around_triangle_to_Delaunay_check_stack( mesh, ti_opp( tii))
    end do

    ! Iteratively flip triangle pairs until the local Delaunay
    ! criterion is satisfied everywhere
    call flip_triangles_until_Delaunay( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_nCge4

  subroutine delete_vertex_nCge4_local_geometry( mesh, vi_kill, &
    vj_clock, vj_focus, vj_anti, vj_opp, &
    ti1, ti2, ti3, ti4, ti_opp, tj1, tj2, tj3, tj4, tj_opp)

    ! Find the local geometry:
    !
    !                tj_opp()
    !
    !  vj_opp(nvj_opp) ------ vj_opp(1)
    !           /  \ ti_opp() /  \
    !     tj1  /    \        /    \  tj4
    !         / ti1  \      /  ti4 \
    !        /        \    /        \
    !   vj_clock ---- vi_kill ---- vj_anti
    !         \          |          /
    !           \  ti2   |   ti3  /
    !             \      |      /
    !        tj2    \    |    /    tj3
    !                 vj_focus

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    integer,                            intent(in   ) :: vi_kill
    integer,                            intent(  out) :: vj_clock, vj_focus, vj_anti
    integer, dimension(:), allocatable, intent(  out) :: vj_opp
    integer,                            intent(  out) :: ti1, ti2, ti3, ti4
    integer, dimension(:), allocatable, intent(  out) :: ti_opp
    integer,                            intent(  out) :: tj1, tj2, tj3, tj4
    integer, dimension(:), allocatable, intent(  out) :: tj_opp

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_nCge4_local_geometry'
    logical                        :: found_it
    integer                        :: iti2, ti, n1, n2, iti1, iti3, iti4, n, tj, tii

    ! Add routine to path
    call init_routine( routine_name)

    ! Vertices

    ! Find vj_focus, which must be at a convex corner
    call find_vj_focus_et_al_nCge4( mesh, vi_kill, &
      vj_clock, vj_focus, vj_anti, vj_opp)

    ! Inside triangles
    found_it = .false.
    do iti2 = 1, mesh%niTri( vi_kill)
      ti = mesh%iTri( vi_kill,iti2)
      do n1 = 1, 3
        n2 = n1 + 1
        if (n2 == 4) n2 = 1

        if (mesh%Tri( ti,n1) == vj_focus .and. mesh%Tri( ti,n2) == vi_kill) then

          iti1 = iti2 - 1
          if (iti1 == 0) iti1 = mesh%niTri( vi_kill)
          ti1 = mesh%iTri( vi_kill,iti1)

          ti2 = mesh%iTri( vi_kill,iti2)

          iti3 = iti2 + 1
          if (iti3 > mesh%niTri( vi_kill)) iti3 = 1
          ti3 = mesh%iTri( vi_kill,iti3)

          iti4 = iti3 + 1
          if (iti4 > mesh%niTri( vi_kill)) iti4 = 1
          ti4 = mesh%iTri( vi_kill,iti4)

          found_it = .true.
          exit

        end if

      end do
      if (found_it) exit
    end do

    if (iti4 > iti1) then
      ! [iti1,iti2,iti3,iti4] form a contiguous block within iTri(vi_kill,:)
      ti_opp = [mesh%iTri( vi_kill, iti4+1:mesh%niTri( vi_kill)), mesh%iTri( vi_kill, 1:iti1-1)]
    else
      ! The block is non-contiguous
      ti_opp = mesh%iTri( vi_kill, iti4+1: iti1-1)
    end if

    ! Outside triangles
    tj1 = 0
    do n = 1, 3
      tj = mesh%TriC( ti1,n)
      if (tj /= ti1 .and. tj /= ti2 .and. tj /= ti3 .and. tj /= ti4 .and. .not. any( tj == ti_opp)) then
        tj1 = tj
        exit
      end if
    end do

    tj2 = 0
    do n = 1, 3
      tj = mesh%TriC( ti2,n)
      if (tj /= ti1 .and. tj /= ti2 .and. tj /= ti3 .and. tj /= ti4 .and. .not. any( tj == ti_opp)) then
        tj2 = tj
        exit
      end if
    end do

    tj3 = 0
    do n = 1, 3
      tj = mesh%TriC( ti3,n)
      if (tj /= ti1 .and. tj /= ti2 .and. tj /= ti3 .and. tj /= ti4 .and. .not. any( tj == ti_opp)) then
        tj3 = tj
        exit
      end if
    end do

    tj4 = 0
    do n = 1, 3
      tj = mesh%TriC( ti4,n)
      if (tj /= ti1 .and. tj /= ti2 .and. tj /= ti3 .and. tj /= ti4 .and. .not. any( tj == ti_opp)) then
        tj4 = tj
        exit
      end if
    end do

    allocate( tj_opp( size( ti_opp,1)), source = 0)
    do tii = 1, size( ti_opp,1)
      ti = ti_opp( tii)
      do n = 1, 3
        tj = mesh%TriC( ti,n)
        if (tj /= ti1 .and. tj /= ti2 .and. tj /= ti3 .and. tj /= ti4 .and. .not. any( tj == ti_opp)) then
          tj_opp( tii) = tj
          exit
        end if
      end do
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_nCge4_local_geometry

  subroutine find_vj_focus_et_al_nCge4( mesh, vi_kill, &
    vj_clock, vj_focus, vj_anti, vj_opp)
    !< Find a valid vj_focus (from where you can "see" all other neighbours of vi_kill)
    !< ...which is not guaranteed, as they might span a non-convex polygon

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    integer,                            intent(in   ) :: vi_kill
    integer,                            intent(  out) :: vj_clock, vj_focus, vj_anti
    integer, dimension(:), allocatable, intent(  out) :: vj_opp

    ! Local variables:
    integer, dimension(:), allocatable :: CC
    integer                            :: nit

    CC = mesh%C( vi_kill, 1: mesh%nC( vi_kill))

    nit = 0
    do while (nit <= mesh%nC( vi_kill))
      nit = nit + 1

      CC = [CC( 2: mesh%nC( vi_kill)), CC( 1)];

      vj_clock = CC( 1);
      vj_focus = CC( 2);
      vj_anti  = CC( 3);
      vj_opp   = CC( 4: mesh%nC( vi_kill));

      if (is_valid_order_vj_focus_nCge4( mesh, vj_clock, vj_focus, vj_anti, vj_opp)) then
        exit
      end if
    end do

  end subroutine find_vj_focus_et_al_nCge4

  function is_valid_order_vj_focus_nCge4( mesh, vj_clock, vj_focus, vj_anti, vj_opp) result( isso)
    ! All new triangles resulting from this ordering must have their
    ! vertices ordered counter-clockwise

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    integer,               intent(in   ) :: vj_clock, vj_focus, vj_anti
    integer, dimension(:), intent(in   ) :: vj_opp
    logical                              :: isso

    ! Local variables:
    integer, dimension(3) :: tri
    integer               :: vii

    isso = .true.

    ! ti3
    tri = [vj_focus, vj_anti, vj_opp( 1)]
    isso = isso .and. is_ordered_counterclockwise( mesh, tri)

    ! ti_opp
    do vii = 1, size( vj_opp,1)-1
      tri = [vj_focus, vj_opp( vii), vj_opp( vii+1)]
      isso = isso .and. is_ordered_counterclockwise( mesh, tri)
    end do

    ! ti3
    tri = [vj_focus, vj_opp( size( vj_opp,1)), vj_clock]
    isso = isso .and. is_ordered_counterclockwise( mesh, tri)

  end function is_valid_order_vj_focus_nCge4

  pure function is_ordered_counterclockwise( mesh, tri) result( isso)
    !< Check if the vertices in a hypothetical triangle are ordered counter-clockwise

    ! In/output variables:
    type(type_mesh),       intent(in   ) :: mesh
    integer, dimension(3), intent(in   ) :: tri
    logical                              :: isso

    ! Local arguments:
    real(dp), dimension(2) :: pa, pb, pc

    pa = mesh%V( tri(1),:)
    pb = mesh%V( tri(2),:)
    pc = mesh%V( tri(3),:)

    isso = cross2( pb-pa, pc-pa) > 0._dp

  end function is_ordered_counterclockwise

  subroutine delete_vertex_nCge4_nC_C( mesh, vi_kill, &
    vj_clock, vj_focus, vj_anti, vj_opp, &
    ti1, ti2, ti3, ti4, ti_opp, tj1, tj2, tj3, tj4, tj_opp)

    ! In/output variables:
    type(type_mesh),       intent(inout) :: mesh
    integer,               intent(in   ) :: vi_kill
    integer,               intent(in   ) :: vj_clock, vj_focus, vj_anti
    integer, dimension(:), intent(in   ) :: vj_opp
    integer,               intent(in   ) :: ti1, ti2, ti3, ti4
    integer, dimension(:), intent(in   ) :: ti_opp
    integer,               intent(in   ) :: tj1, tj2, tj3, tj4
    integer, dimension(:), intent(in   ) :: tj_opp

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_nCge4_nC_C'
    integer                        :: vii

    ! Add routine to path
    call init_routine( routine_name)

    ! vj_focus: replace connection to vi_kill by vj_opp
    call replace_vj_in_C_vi_with_vks( mesh, vj_focus, vi_kill, vj_opp)

    ! vj_clock: remove connection to vi_kill
    call remove_vj_in_C_vi( mesh, vj_clock, vi_kill)

    ! vj_anti: remove connection to vi_kill
    call remove_vj_in_C_vi( mesh, vj_anti, vi_kill)

    ! vj_opp: replace connection to vi_kill by vj_focus
    do vii = 1, size( vj_opp,1)
      call replace_vj_in_C_vi_with_vks( mesh, vj_opp( vii), vi_kill, [vj_focus])
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_nCge4_nC_C

  subroutine delete_vertex_nCge4_niTri_iTri( mesh, vi_kill, &
    vj_clock, vj_focus, vj_anti, vj_opp, &
    ti1, ti2, ti3, ti4, ti_opp, tj1, tj2, tj3, tj4, tj_opp)

    ! In/output variables:
    type(type_mesh),       intent(inout) :: mesh
    integer,               intent(in   ) :: vi_kill
    integer,               intent(  out) :: vj_clock, vj_focus, vj_anti
    integer, dimension(:), intent(  out) :: vj_opp
    integer,               intent(  out) :: ti1, ti2, ti3, ti4
    integer, dimension(:), intent(  out) :: ti_opp
    integer,               intent(  out) :: tj1, tj2, tj3, tj4
    integer, dimension(:), intent(  out) :: tj_opp

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_nCge4_niTri_iTri'

    ! Add routine to path
    call init_routine( routine_name)

    ! vj_focus: add ti_opp(:) after ti3
    call add_tjs_in_iTri_vi_after_ti( mesh, vj_focus, ti3, ti_opp)

    ! vj_clock: remove ti1 from iTri
    call remove_ti_in_iTri_vi( mesh, vj_clock, ti1)

    ! vj_anti: remove ti4 from iTri
    call remove_ti_in_iTri_vi( mesh, vj_anti, ti4)

    ! vj_opp: replace ti1 with ti2, and ti4 with ti3
    call replace_ti_in_iTri_vi_with_tj( mesh, vj_opp( size( vj_opp,1)), ti1, ti2)
    call replace_ti_in_iTri_vi_with_tj( mesh, vj_opp( 1              ), ti4, ti3)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_nCge4_niTri_iTri

  subroutine delete_vertex_nCge4_Tri( mesh, vi_kill, &
    vj_clock, vj_focus, vj_anti, vj_opp, &
    ti1, ti2, ti3, ti4, ti_opp, tj1, tj2, tj3, tj4, tj_opp)

    ! In/output variables:
    type(type_mesh),       intent(inout) :: mesh
    integer,               intent(in   ) :: vi_kill
    integer,               intent(in   ) :: vj_clock, vj_focus, vj_anti
    integer, dimension(:), intent(in   ) :: vj_opp
    integer,               intent(in   ) :: ti1, ti2, ti3, ti4
    integer, dimension(:), intent(in   ) :: ti_opp
    integer,               intent(in   ) :: tj1, tj2, tj3, tj4
    integer, dimension(:), intent(in   ) :: tj_opp

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_nCge4_Tri'
    integer                        :: tii

    ! Add routine to path
    call init_routine( routine_name)

    ! ti1: nothing changes (will be deleted)

    ! ti2: replace vi_kill with vj_opp( end)
    call replace_vi_in_Tri_ti_with_vj( mesh, ti2, vi_kill, vj_opp( size( vj_opp,1)))

    ! ti3: replace vi_kill with vj_opp( 1)
    call replace_vi_in_Tri_ti_with_vj( mesh, ti3, vi_kill, vj_opp( 1))

    ! ti4: nothing changes (will be deleted)

    ! ti_opp: replace vi_kill with vj_focus
    do tii = 1, size( ti_opp,1)
      call replace_vi_in_Tri_ti_with_vj( mesh, ti_opp( tii), vi_kill, vj_focus)
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_nCge4_Tri

  subroutine delete_vertex_nCge4_TriC( mesh, vi_kill, &
    vj_clock, vj_focus, vj_anti, vj_opp, &
    ti1, ti2, ti3, ti4, ti_opp, tj1, tj2, tj3, tj4, tj_opp)

    ! In/output variables:
    type(type_mesh),       intent(inout) :: mesh
    integer,               intent(in   ) :: vi_kill
    integer,               intent(in   ) :: vj_clock, vj_focus, vj_anti
    integer, dimension(:), intent(in   ) :: vj_opp
    integer,               intent(in   ) :: ti1, ti2, ti3, ti4
    integer, dimension(:), intent(in   ) :: ti_opp
    integer,               intent(in   ) :: tj1, tj2, tj3, tj4
    integer, dimension(:), intent(in   ) :: tj_opp

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_nCge4_TriC'

    ! Add routine to path
    call init_routine( routine_name)

    ! ti1: nothing changes (will be deleted)

    ! ti2: replace ti1 in TriC with tj1, and ti3 with ti_opp( end)
    call replace_tj_in_TriC_ti_with_tk( mesh, ti2, ti1, tj1)
    if (size( ti_opp,1) >= 1) then
      call replace_tj_in_TriC_ti_with_tk( mesh, ti2, ti3, ti_opp( size( ti_opp,1)))
    end if

    ! ti3: replace ti4 in TriC with tj4, and ti2 with ti_opp( 1)
    call replace_tj_in_TriC_ti_with_tk( mesh, ti3, ti4, tj4)
    if (size( ti_opp,1) >= 1) then
      call replace_tj_in_TriC_ti_with_tk( mesh, ti3, ti2, ti_opp( 1))
    end if

    ! ti4: nothing changes (will be deleted)

    ! ti_opp( 1): replace ti4 in TriC with ti3
    if (size( ti_opp,1) >= 1) then
      call replace_tj_in_TriC_ti_with_tk( mesh, ti_opp( 1), ti4, ti3)
    end if

    ! ti_opp( end): replace ti1 in TriC with ti2
    if (size( ti_opp,1) >= 1) then
      call replace_tj_in_TriC_ti_with_tk( mesh, ti_opp( size( ti_opp,1)), ti1, ti2)
    end if

    ! tj1: replace ti1 in TriC with ti2
    call replace_tj_in_TriC_ti_with_tk( mesh, tj1, ti1, ti2)

    ! tj2: nothing changes
    ! tj3: nothing changes

    ! tj4: replace ti4 in TriC with ti3
    call replace_tj_in_TriC_ti_with_tk( mesh, tj4, ti4, ti3)

    ! tj_opp: nothing changes

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_nCge4_TriC

  ! General operations on primary mesh data
  ! =======================================

  subroutine replace_vj_in_C_vi_with_vks( mesh, vi, vj, vks)
    ! Replace vj in C(vi,:) with vk(:)

    ! In/output variables:
    type(type_mesh),       intent(inout) :: mesh
    integer,               intent(in   ) :: vi, vj
    integer, dimension(:), intent(in   ) :: vks

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'replace_vj_in_C_vi_with_vks'
    logical                        :: found_it
    integer                        :: ci

    ! Add routine to path
    call init_routine( routine_name)

    found_it = .false.
    do ci = 1, mesh%nC( vi)
      if (mesh%C( vi,ci) == vj) then
        found_it = .true.
        mesh%C( vi,1:mesh%nC( vi) - 1 + size( vks,1)) = &
          [mesh%C( vi,1:ci-1), vks, mesh%C( vi,ci+1:mesh%nC( vi))]
        mesh%nC( vi) = mesh%nC( vi) - 1 + size( vks,1)
        exit
      end if
    end do
    if (.not. found_it) call crash('couldnt find vj in C(vi,:)')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine replace_vj_in_C_vi_with_vks

  subroutine remove_vj_in_C_vi( mesh, vi, vj)
    ! Remove vj in C(vi,:)

    ! In/output variables:
    type(type_mesh),       intent(inout) :: mesh
    integer,               intent(in   ) :: vi, vj

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remove_vj_in_C_vi'
    logical                        :: found_it
    integer                        :: ci

    ! Add routine to path
    call init_routine( routine_name)

    found_it = .false.
    do ci = 1, mesh%nC( vi)
      if (mesh%C( vi,ci) == vj) then
        found_it = .true.
        mesh%C( vi,1:mesh%nC( vi)) = &
          [mesh%C( vi,1:ci-1), mesh%C( vi,ci+1:mesh%nC( vi)), 0]
        mesh%nC( vi) = mesh%nC( vi) - 1
        exit
      end if
    end do
    if (.not. found_it) call crash('couldnt find vj in C(vi,:)')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remove_vj_in_C_vi

  subroutine remove_ti_in_iTri_vi( mesh, vi, ti)
    ! Remove ti in iTri( vi,:)

    ! In/output variables:
    type(type_mesh),       intent(inout) :: mesh
    integer,               intent(in   ) :: vi, ti

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remove_ti_in_iTri_vi'
    logical                        :: found_it
    integer                        :: iti

    ! Add routine to path
    call init_routine( routine_name)

    found_it = .false.
    do iti = 1, mesh%niTri( vi)
      if (mesh%iTri( vi,iti) == ti) then
        found_it = .true.
        mesh%iTri( vi,1:mesh%niTri( vi)) = &
          [mesh%iTri( vi,1:iti-1), mesh%iTri( vi,iti+1:mesh%niTri( vi)), 0]
        mesh%niTri( vi) = mesh%niTri( vi) - 1
        exit
      end if
    end do
    if (.not. found_it) call crash('couldnt find ti in iTri(vi,:)')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine remove_ti_in_iTri_vi

  subroutine replace_ti_in_iTri_vi_with_tj( mesh, vi, ti, tj)
    ! Replace ti in iTri(vi,:) with tj

    ! In/output variables:
    type(type_mesh),       intent(inout) :: mesh
    integer,               intent(in   ) :: vi, ti, tj

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'replace_ti_in_iTri_vi_with_tj'
    logical                        :: found_it
    integer                        :: iti

    ! Add routine to path
    call init_routine( routine_name)

    found_it = .false.
    do iti = 1, mesh%niTri( vi)
      if (mesh%iTri( vi,iti) == ti) then
        found_it = .true.
        mesh%iTri( vi,iti) = tj
        exit
      end if
    end do
    if (.not. found_it) call crash('couldnt find ti in iTri(vi,:)')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine replace_ti_in_iTri_vi_with_tj

  subroutine add_tjs_in_iTri_vi_after_ti( mesh, vi, ti, tjs)
    ! Add tjs in iTri( vi,:) after tj

    ! In/output variables:
    type(type_mesh),       intent(inout) :: mesh
    integer,               intent(in   ) :: vi, ti
    integer, dimension(:), intent(in   ) :: tjs

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_tjs_in_iTri_vi_after_ti'
    logical                        :: found_it
    integer                        :: iti

    ! Add routine to path
    call init_routine( routine_name)

    found_it = .false.
    do iti = 1, mesh%niTri( vi)
      if (mesh%iTri( vi,iti) == ti) then
        found_it = .true.
        mesh%iTri( vi,1:mesh%niTri( vi) + size( tjs,1)) = &
          [mesh%iTri( vi,1:iti), tjs, mesh%iTri( vi,iti+1:mesh%niTri( vi))]
        mesh%niTri( vi) = mesh%niTri( vi) + size( tjs,1)
        exit
      end if
    end do
    if (.not. found_it) call crash('couldnt find ti in iTri(vi,:)')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_tjs_in_iTri_vi_after_ti

  subroutine replace_vi_in_Tri_ti_with_vj( mesh, ti, vi, vj)
    ! Replace vi in Tri(ti,:) with tj

    ! In/output variables:
    type(type_mesh),       intent(inout) :: mesh
    integer,               intent(in   ) :: ti, vi, vj

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'replace_vi_in_Tri_ti_with_vj'
    logical                        :: found_it
    integer                        :: n

    ! Add routine to path
    call init_routine( routine_name)

    found_it = .false.
    do n = 1, 3
      if (mesh%Tri( ti,n) == vi) then
        found_it = .true.
        mesh%Tri( ti,n) = vj
        exit
      end if
    end do
    if (.not. found_it) call crash('couldnt find vi in Tri(ti,:)')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine replace_vi_in_Tri_ti_with_vj

  subroutine replace_tj_in_TriC_ti_with_tk( mesh, ti, tj, tk)
    ! Replace tj in TriC(ti,:) with tk

    ! In/output variables:
    type(type_mesh),       intent(inout) :: mesh
    integer,               intent(in   ) :: ti, tj, tk

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'replace_tj_in_TriC_ti_with_tk'
    logical                        :: found_it
    integer                        :: n

    ! Add routine to path
    call init_routine( routine_name)

    found_it = .false.
    do n = 1, 3
      if (mesh%TriC( ti,n) == tj) then
        found_it = .true.
        mesh%TriC( ti,n) = tk
        exit
      end if
    end do
    if (.not. found_it) call crash('couldnt find tj in TriC(ti,:)')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine replace_tj_in_TriC_ti_with_tk

  subroutine delete_vertex_V( mesh, vi_kill, vi_new2vi_old, vi_old2vi_new)

    ! In/output variables:
    type(type_mesh),                    intent(inout) :: mesh
    integer,                            intent(in   ) :: vi_kill
    integer, dimension(:), allocatable, intent(inout) :: vi_new2vi_old, vi_old2vi_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_V'
    integer                        :: n

    ! Add routine to path
    call init_routine( routine_name)

    call switch_vertices( mesh, vi_kill, mesh%nV)

    call switch_rows( vi_new2vi_old, vi_kill, mesh%nV)
    call switch_rows( vi_old2vi_new, vi_kill, mesh%nV)
    vi_new2vi_old = vi_new2vi_old( 1: mesh%nV)
    vi_old2vi_new( vi_kill) = 0

    mesh%nV    = mesh%nV - 1
    mesh%V     = mesh%V(     1:mesh%nV,:)
    mesh%nC    = mesh%nC(    1:mesh%nV  )
    mesh%C     = mesh%C(     1:mesh%nV,:)
    mesh%niTri = mesh%niTri( 1:mesh%nV  )
    mesh%iTri  = mesh%iTri(  1:mesh%nV,:)
    mesh%VBI   = mesh%VBI(   1:mesh%nV  )

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_V

  subroutine delete_vertex_Tri( mesh, ti_kill, ti_new2ti_old, ti_old2ti_new)

    ! In/output variables:
    type(type_mesh),                    intent(inout) :: mesh
    integer,                            intent(in   ) :: ti_kill
    integer, dimension(:), allocatable, intent(inout) :: ti_new2ti_old, ti_old2ti_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'delete_vertex_Tri'
    logical                        :: found_it
    integer                        :: n

    ! Add routine to path
    call init_routine( routine_name)

    call switch_triangles( mesh, ti_kill, mesh%nTri)

    call switch_rows( ti_new2ti_old, ti_kill, mesh%nTri)
    call switch_rows( ti_old2ti_new, ti_kill, mesh%nTri)
    ti_new2ti_old = ti_new2ti_old( 1: mesh%nTri)
    ti_old2ti_new( ti_kill) = 0

    mesh%nTri  = mesh%nTri - 1
    mesh%Tri   = mesh%Tri(   1:mesh%nTri,:)
    mesh%Tricc = mesh%Tricc( 1:mesh%nTri,:)
    mesh%TriC  = mesh%TriC(  1:mesh%nTri,:)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertex_Tri

end module delete_vertices
