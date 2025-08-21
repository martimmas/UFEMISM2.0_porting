module create_graphs_from_masked_mesh

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph
  use graph_memory, only: allocate_graph_primary, crop_graph_primary
  use mpi_distributed_memory, only: gather_to_all
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared
  use plane_geometry, only: mirror_p_across_qr
  use assertions_basic, only: assert
  use graph_contiguous_domains, only: enforce_contiguous_process_domains_graph
  use graph_parallelisation, only: setup_graph_parallelisation
  use tests_main, only: test_graph_is_self_consistent
  use parameters, only: NaN
  use mpi_basic, only: par

  implicit none

  private

  public :: create_graph_from_masked_mesh_a, create_graph_from_masked_mesh_b

contains

  subroutine create_graph_from_masked_mesh_a( mesh, mask_a, nz, graph)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_a
    integer,                               intent(in   ) :: nz
    type(type_graph),                      intent(  out) :: graph

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_graph_from_masked_mesh_a'
    logical, dimension(1:mesh%nV  ) :: mask_a_tot
    integer                         :: ni, vi, ci, vj, nj

    ! Add routine to path
    call init_routine( routine_name)

    call gather_to_all( mask_a, mask_a_tot)

    ! Allocate graph and set metadata
    call allocate_graph_primary( graph, mesh%nV, mesh%nV, maxval( mesh%nC))

    graph%parent_mesh_name = trim( mesh%name) // '_vertices'
    graph%xmin             = mesh%xmin
    graph%xmax             = mesh%xmax
    graph%ymin             = mesh%ymin
    graph%ymax             = mesh%ymax

    ! Create vertex-to-regular-node mapping
    graph%n = 0
    do vi = 1, mesh%nV
      if (mask_a_tot( vi)) then
        graph%n = graph%n + 1
        graph%mi2ni( vi     ) = graph%n
        graph%ni2mi( graph%n) = vi
      end if
    end do

    ! Construct regular nodes and their interconnectivity
    do ni = 1, graph%n

      ! This regular graph node corresponds to a masked vertex
      vi = graph%ni2mi( ni)

      graph%V         ( ni,:) = mesh%V( vi,:)
      graph%is_ghost  ( ni  ) = .false.
      graph%V_ghost_BC( ni,:) = NaN
      graph%ghost_nhat( ni,:) = NaN

      ! Add connections to other regular nodes
      graph%nC( ni) = 0
      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        if (mask_a_tot( vj)) then
          ! This connection points to another regular node
          nj = graph%mi2ni( vj)
          graph%nC( ni) = graph%nC( ni) + 1
          graph%C( ni, graph%nC( ni)) = nj
        end if
      end do

    end do

    ! Finalisation
    call crop_graph_primary( graph)
    call enforce_contiguous_process_domains_graph( graph)
    call setup_graph_parallelisation( graph, nz)

#if (DO_ASSERTIONS)
    call assert(test_graph_is_self_consistent( mesh, graph), 'inconsistent graph connectivity')
#endif

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 4)

  end subroutine create_graph_from_masked_mesh_a

  subroutine create_graph_from_masked_mesh_b( mesh, mask_a, nz, graph)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_a
    integer,                               intent(in   ) :: nz
    type(type_graph),                      intent(  out) :: graph

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'create_graph_from_masked_mesh_b'
    logical, dimension(1:mesh%nTri) :: mask_b_tot
    integer                         :: ni, ti, ng, n, ei, tj, nj, til, tir, vi, vj, vi1, vi2
    real(dp), dimension(2)          :: p, q, r, s, normal_vector
    integer                         :: n_unmasked_neighbours

    ! Add routine to path
    call init_routine( routine_name)

    ! Calculate b-grid mask
    call calc_mask_b( mesh, mask_a, mask_b_tot)

    ! Allocate graph and set metadata
    call allocate_graph_primary( graph, mesh%nTri*3, mesh%nTri, 3)

    graph%parent_mesh_name = trim( mesh%name) // '_triangles'
    graph%xmin             = mesh%xmin
    graph%xmax             = mesh%xmax
    graph%ymin             = mesh%ymin
    graph%ymax             = mesh%ymax

    ! Create triangle-to-regular-node mapping
    graph%n = 0
    do ti = 1, mesh%nTri
      if (mask_b_tot( ti)) then
        graph%n = graph%n + 1
        graph%mi2ni( ti     ) = graph%n
        graph%ni2mi( graph%n) = ti
      end if
    end do

    ! Construct regular nodes and their interconnectivity
    do ni = 1, graph%n

      ! This regular graph node corresponds to a masked triangle
      ti = graph%ni2mi( ni)

      graph%V         ( ni,:) = mesh%TriGC( ti,:)
      graph%is_ghost  ( ni  ) = .false.
      graph%V_ghost_BC( ni,:) = NaN
      graph%ghost_nhat( ni,:) = NaN

      ! Add connections to other regular nodes
      graph%nC( ni) = 0
      do n = 1, 3
        tj = mesh%TriC( ti,n)
        nj = graph%mi2ni( tj)
        if (nj > 0) then
          ! This connection points to another regular node
          graph%nC( ni) = graph%nC( ni) + 1
          graph%C( ni, graph%nC( ni)) = nj
        end if
      end do

    end do

    ! Add ghost nodes
    do ti = 1, mesh%nTri
      if (mask_b_tot( ti)) then

        n_unmasked_neighbours = 0
        do n = 1, 3
          tj = mesh%TriC( ti,n)
          if (.not. mask_b_tot( tj)) then
            n_unmasked_neighbours = n_unmasked_neighbours + 1
          end if
        end do

        if (n_unmasked_neighbours == 1) then
          call add_ghost_node_1_unmasked_neighbour( mesh, graph, mask_b_tot, ti)
        elseif (n_unmasked_neighbours == 2) then
          call add_ghost_node_2_unmasked_neighbours( mesh, graph, mask_b_tot, ti)
        end if

      end if
    end do

    ! Finalisation
    call crop_graph_primary( graph)
    call enforce_contiguous_process_domains_graph( graph)
    call setup_graph_parallelisation( graph, nz)

#if (DO_ASSERTIONS)
    call assert(test_graph_is_self_consistent( mesh, graph), 'inconsistent graph connectivity')
#endif

    ! Finalise routine path
    call finalise_routine( routine_name, n_extra_MPI_windows_expected = 4)

  end subroutine create_graph_from_masked_mesh_b

  subroutine calc_mask_b( mesh, mask_a, mask_b_tot)

    ! In/output variables:
    type(type_mesh),                       intent(in   ) :: mesh
    logical, dimension(mesh%vi1:mesh%vi2), intent(in   ) :: mask_a
    logical, dimension(mesh%nTri),         intent(  out) :: mask_b_tot

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calc_mask_b'
    logical, dimension(mesh%nV  )  :: mask_a_tot
    integer                        :: ti, n_a, n, vi

    ! Add routine to path
    call init_routine( routine_name)

    call gather_to_all( mask_a, mask_a_tot)

    ! b mask
    mask_b_tot = .false.
    do ti = 1, mesh%nTri
      n_a = 0
      do n = 1, 3
        vi = mesh%Tri( ti,n)
        if (mask_a_tot( vi)) n_a = n_a + 1
      end do
      if (n_a == 3) then
        mask_b_tot( ti) = .true.
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_mask_b

  subroutine add_ghost_node_1_unmasked_neighbour( mesh, graph, mask_b_tot, ti)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_graph),                intent(inout) :: graph
    logical, dimension(1:mesh%nTri), intent(in   ) :: mask_b_tot
    integer,                         intent(in   ) :: ti

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_ghost_node_1_unmasked_neighbour'
    integer                        :: ni_reg, via, vib, vic, tj_unmasked, tjb, tjc, ni_ghost
    real(dp), dimension(2)         :: V_reg, q, r, V_ghost, V_ghost_BC, ghost_nhat

    ! Add routine to path
    call init_routine( routine_name)

    ! The regular node adjacent to the new ghost node
    ni_reg = graph%mi2ni( ti)
    V_reg  = graph%V( ni_reg,:)

    call add_ghost_node_1_unmasked_neighbour_local_geometry( mesh, mask_b_tot, ti, &
      via, vib, vic, tj_unmasked, tjb, tjc)

    ! Define the ghost node as the mirror image of the regular node across the line [vib,vic]
    q = mesh%V( vib,:)
    r = mesh%V( vic,:)
    V_ghost = mirror_p_across_qr( V_reg, q, r)

    ! Use this ghost node to calculate boundary conditions at the domain margin,
    ! i.e. at the midpoint between the ghost node and the adjacent regular node
    V_ghost_BC = (V_ghost + graph%V( ni_reg,:)) / 2._dp

    ! Calculate the outward unit normal vector at this ghost node
    ghost_nhat = (V_ghost - V_reg) / norm2( V_ghost - V_reg)

    ! Add the ghost node to the graph
    graph%n = graph%n + 1
    ni_ghost = graph%n

    graph%V         ( ni_ghost,:) = V_ghost
    graph%nC        ( ni_ghost  ) = 1
    graph%C         ( ni_ghost,1) = ni_reg
    graph%is_ghost  ( ni_ghost  ) = .true.
    graph%V_ghost_BC( ni_ghost,:) = V_ghost_BC
    graph%ghost_nhat( ni_ghost,:) = ghost_nhat

    ! Add the reverse connection to the regular node
    graph%nC( ni_reg) = graph%nC( ni_reg) + 1
    graph%C( ni_reg, graph%nC( ni_reg)) = ni_ghost

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_ghost_node_1_unmasked_neighbour

  subroutine add_ghost_node_1_unmasked_neighbour_local_geometry( mesh, mask_b_tot, ti, &
    via, vib, vic, tj_unmasked, tjb, tjc)
    ! The local geometry looks like this:
    !             \   /
    !          --- vic ---
    !             /   \
    !            /     \
    !      tjb  /       \  tj_unmasked
    !          /   ti    \
    !     \   /           \   /
    !  --- via ----------- vib ---
    !     /   \           /   \
    !              tjc

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    logical, dimension(1:mesh%nTri), intent(in   ) :: mask_b_tot
    integer,                         intent(in   ) :: ti
    integer,                         intent(  out) :: via, vib, vic, tj_unmasked, tjb, tjc

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_ghost_node_1_unmasked_neighbour_local_geometry'
    integer                        :: n1, n2, n3, tj

    ! Add routine to path
    call init_routine( routine_name)

    do n1 = 1, 3

      n2 = n1 + 1
      if (n2 == 4) n2 = 1
      n3 = n2 + 1
      if (n3 == 4) n3 = 1

      tj = mesh%TriC( ti,n1)

      if (.not. mask_b_tot( tj)) then
        via         = mesh%Tri ( ti,n1)
        vib         = mesh%Tri ( ti,n2)
        vic         = mesh%Tri ( ti,n3)
        tj_unmasked = mesh%TriC( ti,n1)
        tjb         = mesh%TriC( ti,n2)
        tjc         = mesh%TriC( ti,n3)
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_ghost_node_1_unmasked_neighbour_local_geometry

  subroutine add_ghost_node_2_unmasked_neighbours( mesh, graph, mask_b_tot, ti)

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_graph),                intent(inout) :: graph
    logical, dimension(1:mesh%nTri), intent(in   ) :: mask_b_tot
    integer,                         intent(in   ) :: ti

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_ghost_node_2_unmasked_neighbours'
    integer                        :: ni_reg, via, vib, vic, tj_masked, tjb, tjc, ni_ghost
    real(dp), dimension(2)         :: V_reg, V_ghost, V_ghost_BC, ghost_nhat

    ! Add routine to path
    call init_routine( routine_name)

    ! The regular node adjacent to the new ghost node
    ni_reg = graph%mi2ni( ti)
    V_reg  = graph%V( ni_reg,:)

    call add_ghost_node_2_unmasked_neighbours_local_geometry( mesh, mask_b_tot, ti, &
      via, vib, vic, tj_masked, tjb, tjc)

    V_ghost = mesh%TriGC( ti,:) + (mesh%TriGC( ti,:) - mesh%TriGC( tj_masked,:))

    ! Use this ghost node to calculate boundary conditions at the domain margin,
    ! i.e. at the midpoint between the ghost node and the adjacent regular node
    V_ghost_BC = (V_ghost + graph%V( ni_reg,:)) / 2._dp

    ! Calculate the outward unit normal vector at this ghost node
    ghost_nhat = (V_ghost - V_reg) / norm2( V_ghost - V_reg)

    ! Add the ghost node to the graph
    graph%n = graph%n + 1
    ni_ghost = graph%n

    graph%V         ( ni_ghost,:) = V_ghost
    graph%nC        ( ni_ghost  ) = 1
    graph%C         ( ni_ghost,1) = ni_reg
    graph%is_ghost  ( ni_ghost  ) = .true.
    graph%V_ghost_BC( ni_ghost,:) = V_ghost_BC
    graph%ghost_nhat( ni_ghost,:) = ghost_nhat

    ! Add the reverse connection to the regular node
    graph%nC( ni_reg) = graph%nC( ni_reg) + 1
    graph%C( ni_reg, graph%nC( ni_reg)) = ni_ghost

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_ghost_node_2_unmasked_neighbours

  subroutine add_ghost_node_2_unmasked_neighbours_local_geometry( mesh, mask_b_tot, ti, &
    via, vib, vic, tj_masked, tjb, tjc)
    ! The local geometry looks like this:
    !             \   /
    !          --- vic ---
    !             /   \
    !            /     \
    !      tjb  /       \  tj_masked
    !          /   ti    \
    !     \   /           \   /
    !  --- via ----------- vib ---
    !     /   \           /   \
    !              tjc

    ! In/output variables:
    type(type_mesh),                 intent(in   ) :: mesh
    logical, dimension(1:mesh%nTri), intent(in   ) :: mask_b_tot
    integer,                         intent(in   ) :: ti
    integer,                         intent(  out) :: via, vib, vic, tj_masked, tjb, tjc

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_ghost_node_2_unmasked_neighbours_local_geometry'
    integer                        :: n1, n2, n3, tj

    ! Add routine to path
    call init_routine( routine_name)

    do n1 = 1, 3

      n2 = n1 + 1
      if (n2 == 4) n2 = 1
      n3 = n2 + 1
      if (n3 == 4) n3 = 1

      tj = mesh%TriC( ti,n1)

      if (mask_b_tot( tj)) then
        via       = mesh%Tri ( ti,n1)
        vib       = mesh%Tri ( ti,n2)
        vic       = mesh%Tri ( ti,n3)
        tj_masked = mesh%TriC( ti,n1)
        tjb       = mesh%TriC( ti,n2)
        tjc       = mesh%TriC( ti,n3)
      end if

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_ghost_node_2_unmasked_neighbours_local_geometry

end module create_graphs_from_masked_mesh