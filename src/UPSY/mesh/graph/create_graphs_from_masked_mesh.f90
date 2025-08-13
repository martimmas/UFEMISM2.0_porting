module create_graphs_from_masked_mesh

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph
  use mpi_distributed_memory, only: gather_to_all
  use mpi_distributed_shared_memory, only: allocate_dist_shared
  use plane_geometry, only: mirror_p_across_qr
  use assertions_basic, only: assert

  implicit none

  private

  public :: create_graph_from_masked_mesh_a, create_graph_from_masked_mesh_b, &
    test_graph_connectivity_is_self_consistent

  contains

  subroutine create_graph_from_masked_mesh_a( mesh, mask_a, graph)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_a
    type(type_graph),                      intent(  out) :: graph

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_graph_from_masked_mesh_a'
    logical, dimension(1:mesh%nV  ) :: mask_a_tot
    integer                         :: ni, vi, ci, vj, nj

    ! Add routine to path
    call init_routine( routine_name)

    call gather_to_all( mask_a, mask_a_tot)

    ! Set graph metadata
    graph%parent_mesh_name = mesh%name
    graph%n                = count( mask_a_tot)
    graph%nn               = 0
    graph%ng               = 0
    graph%nC_mem           = maxval( mesh%nC)

    ! Allocate memory
    allocate( graph%mi2ni     ( mesh%nV              ), source = 0)
    allocate( graph%ni2mi     ( graph%n              ), source = 0)

    allocate( graph%V         ( graph%n, 2           ), source = 0._dp)
    allocate( graph%nC        ( graph%n              ), source = 0)
    allocate( graph%C         ( graph%n, graph%nC_mem), source = 0)

    allocate( graph%is_ghost  ( graph%n              ), source = .false.)
    allocate( graph%ghost_nhat( graph%n, 2           ), source = 0._dp)

    ! Create vertex-to-node mapping
    ni = 0
    do vi = 1, mesh%nV
      if (mask_a_tot( vi)) then
        ni = ni + 1
        graph%mi2ni( vi) = ni
        graph%ni2mi( ni) = vi
      end if
    end do

    ! Construct reduced graph
    do vi = 1, mesh%nV
      if (mask_a_tot( vi)) then

        ni = graph%mi2ni( vi)

        ! Add this masked vertex
        graph%V ( ni,:) = mesh%V( vi,:)

        ! Add its connectivity
        graph%nC( ni) = 0
        do ci = 1, mesh%nC( vi)

          vj = mesh%C( vi,ci)

          if (mask_a_tot( vj)) then
            ! This connection points to a masked vertex

            nj = graph%mi2ni( vj)
            graph%nC( ni) = graph%nC( ni) + 1
            graph%C( ni, graph%nC( ni)) = nj

          end if
        end do

      end if
    end do

#if (DO_ASSERTIONS)
    call assert(test_graph_connectivity_is_self_consistent( graph), 'inconsistent graph connectivity')
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_graph_from_masked_mesh_a

  subroutine create_graph_from_masked_mesh_b( mesh, mask_a, graph)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_a
    type(type_graph),                      intent(  out) :: graph

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_graph_from_masked_mesh_b'
    logical, dimension(1:mesh%nV  ) :: mask_a_tot
    logical, dimension(1:mesh%nTri) :: mask_b_tot
    integer                         :: n_mask_b
    logical, dimension(1:mesh%nE  ) :: mask_boundary_c_tot
    integer                         :: n_boundary_c
    integer                         :: ni, ti, ng, n, ei, tj, nj, til, tir, vi, vj, vi1, vi2
    real(dp), dimension(2)          :: p, q, r, s, normal_vector

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate b-grid and c-grid masks
    call gather_to_all( mask_a, mask_a_tot)
    call calc_masks_b_c( mesh, mask_a_tot, mask_b_tot, n_mask_b, mask_boundary_c_tot, n_boundary_c)

    ! Set graph metadata
    graph%parent_mesh_name = mesh%name
    graph%n                = n_mask_b + n_boundary_c
    graph%nn               = n_mask_b
    graph%ng               = n_boundary_c
    graph%nC_mem           = 3

    ! Allocate memory
    allocate( graph%mi2ni     ( mesh%nTri            ), source = 0)
    allocate( graph%ni2mi     ( graph%n              ), source = 0)

    allocate( graph%V         ( graph%n, 2           ), source = 0._dp)
    allocate( graph%nC        ( graph%n              ), source = 0)
    allocate( graph%C         ( graph%n, graph%nC_mem), source = 0)

    allocate( graph%is_ghost  ( graph%n              ), source = .false.)
    allocate( graph%ghost_nhat( graph%n, 2           ), source = 0._dp)

    ! Create triangle-to-node mapping
    ni = 0
    do ti = 1, mesh%nTri
      if (mask_b_tot( ti)) then
        ni = ni + 1
        graph%mi2ni( ti) = ni
        graph%ni2mi( ni) = ti
      end if
    end do

    ! Construct reduced graph
    ng = ni
    do ti = 1, mesh%nTri
      if (mask_b_tot( ti)) then

        ni = graph%mi2ni( ti)

        ! Add this masked triangle
        graph%V ( ni,:) = mesh%TriGC( ti,:)
        graph%nC( ni  ) = 3

        ! Add its connectivity
        do n = 1, 3

          ei = mesh%TriE( ti,n)
          tj = mesh%TriC( ti,n)

          if (mask_b_tot( tj)) then
            ! This connection points to a masked triangle

            nj = graph%mi2ni( tj)
            graph%C( ni,n) = nj

          else
            ! This connection points to a ghost node
            ng = ng + 1
            graph%C( ni,n) = ng

            ! Add the ghost node as well
            graph%is_ghost( ng  ) = .true.
            graph%ni2mi   ( ng  ) = ei
            graph%nC      ( ng  ) = 1
            graph%C       ( ng,1) = ni

            til = mesh%ETri( ei,1)
            tir = mesh%ETri( ei,2)
            vi  = mesh%EV( ei,1)
            vj  = mesh%EV( ei,2)

            ! Define vi1 and vi2 such that the line from vi1 to vi2 has
            ! the masked mesh interior on its left-hand side
            if (ti == til) then
              vi1 = vi
              vi2 = vj
            elseif (ti == tir) then
              vi1 = vj
              vi2 = vi
            else
              call crash('inconsistency in EV/ETri')
            end if

            ! Calculate coordinates of the ghost node
            p = mesh%TriGC( ti,:)
            q = mesh%V( vi1,:)
            r = mesh%V( vi2,:)
            s = mirror_p_across_qr( p, q, r)
            graph%V( ng,:) = s

            ! Calculate unit normal vector
            normal_vector = graph%V( ng,:) - mesh%TriGC( ti,:)
            graph%ghost_nhat( ng,:) = normal_vector / norm2( normal_vector)

          end if
        end do

      end if
    end do

#if (DO_ASSERTIONS)
    call assert(test_graph_connectivity_is_self_consistent( graph), 'inconsistent graph connectivity')
#endif

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_graph_from_masked_mesh_b

  subroutine calc_masks_b_c( mesh, mask_a_tot, mask_b_tot, n_mask_b, mask_boundary_c_tot, n_boundary_c)

    ! In/output variables:
    type(type_mesh),               intent(in   ) :: mesh
    logical, dimension(mesh%nV  ), intent(in   ) :: mask_a_tot
    logical, dimension(mesh%nTri), intent(  out) :: mask_b_tot
    integer,                       intent(  out) :: n_mask_b
    logical, dimension(mesh%nE  ), intent(  out) :: mask_boundary_c_tot
    integer,                       intent(  out) :: n_boundary_c

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'calc_masks_b_c'
    integer                         :: ti, n_a, n, vi, ei, tj

    ! Add routine to path
    call init_routine( routine_name)

    ! b mask
    mask_b_tot = .false.
    n_mask_b = 0
    do ti = 1, mesh%nTri
      n_a = 0
      do n = 1, 3
        vi = mesh%Tri( ti,n)
        if (mask_a_tot( vi)) n_a = n_a + 1
      end do
      if (n_a == 3) then
        mask_b_tot( ti) = .true.
        n_mask_b = n_mask_b + 1
      end if
    end do

    ! c boundary mask
    mask_boundary_c_tot = .false.
    n_boundary_c = 0
    do ei = 1, mesh%nE
      if (mesh%EBI( ei) == 0) then
        ti = mesh%ETri( ei,1)
        tj = mesh%ETri( ei,2)
        if (xor( mask_b_tot( ti), mask_b_tot( tj))) then
          mask_boundary_c_tot( ei) = .true.
          n_boundary_c = n_boundary_c + 1
        end if
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_masks_b_c

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

    ! Finalise routine path
    call finalise_routine( routine_name)

  end function test_graph_connectivity_is_self_consistent

end module create_graphs_from_masked_mesh