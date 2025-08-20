module ut_mesh_graphs_mapping

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mpi_basic, only: par
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph
  use create_graphs_from_masked_mesh, only: create_graph_from_masked_mesh_a, create_graph_from_masked_mesh_b
  use mpi_f08, only: MPI_WIN
  use mpi_distributed_memory, only: gather_to_all
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use deallocate_dist_shared_mod, only: deallocate_dist_shared
  use mesh_graph_mapping, only: map_mesh_vertices_to_graph, map_graph_to_mesh_vertices, &
    map_mesh_triangles_to_graph, map_graph_to_mesh_triangles

  implicit none

  private

  public :: test_mesh_graph_mapping

contains

  subroutine test_mesh_graph_mapping( test_name_parent, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_mesh_graph_mapping'
    character(len=1024), parameter     :: test_name_local = 'mapping'
    character(len=1024)                :: test_name
    logical, dimension(:), allocatable :: mask_a
    integer                            :: nz
    integer                            :: vi
    type(type_graph)                   :: graph_a, graph_b

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Define a masked set of vertices
    allocate( mask_a( mesh%vi1:mesh%vi2), source = .false.)
    do vi = mesh%vi1, mesh%vi2
      if (hypot( mesh%V( vi,1), mesh%V( vi,2)) < mesh%xmax * 0.5_dp) mask_a( vi) = .true.
    end do

    ! Create graphs from the masked vertices and triangles
    nz = 12
    call create_graph_from_masked_mesh_a( mesh, mask_a, nz, graph_a)
    call create_graph_from_masked_mesh_b( mesh, mask_a, nz, graph_b)

    call test_mesh_vertices_to_graph_mapping ( test_name, mesh, graph_a)
    call test_graph_to_mesh_vertices_mapping ( test_name, mesh, graph_a, mask_a)
    call test_mesh_triangles_to_graph_mapping( test_name, mesh, graph_b)
    call test_graph_to_mesh_triangles_mapping( test_name, mesh, graph_b)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_graph_mapping

  ! Mesh vertices to graph
  ! ======================

  subroutine test_mesh_vertices_to_graph_mapping( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_mesh_vertices_to_graph_mapping'
    character(len=1024), parameter :: test_name_local = 'mesh_vertices_to_graph'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_mesh_vertices_to_graph_mapping_logical_2D( test_name, mesh, graph)
    call test_mesh_vertices_to_graph_mapping_logical_3D( test_name, mesh, graph)
    call test_mesh_vertices_to_graph_mapping_int_2D    ( test_name, mesh, graph)
    call test_mesh_vertices_to_graph_mapping_int_3D    ( test_name, mesh, graph)
    call test_mesh_vertices_to_graph_mapping_dp_2D     ( test_name, mesh, graph)
    call test_mesh_vertices_to_graph_mapping_dp_3D     ( test_name, mesh, graph)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping

  subroutine test_mesh_vertices_to_graph_mapping_logical_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_mesh_vertices_to_graph_mapping_logical_2D'
    character(len=1024), parameter             :: test_name_local = 'logical_2D'
    character(len=1024)                        :: test_name
    logical, dimension(:), allocatable         :: d_mesh_loc_ex
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    logical, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    logical, dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                    :: vi, ni, d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex     ( mesh%vi1 :mesh%vi2 ), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_V%n_nih)
    d_mesh_nih_ex( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      d = ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2)))
      d_mesh_loc_ex( vi) = modulo( d,2) == 0
      d_mesh_nih_ex( vi) = modulo( d,2) == 0
    end do
    do ni = graph%ni1, graph%ni2
      d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)))
      d_graph_nih_ex( ni) = modulo( d,2) == 0
    end do

    ! Map data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)
    call unit_test( test_eqv( d_graph_nih_ex, d_graph_nih_mapped), trim( test_name) // '/dist_hybrid')

    call map_mesh_vertices_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)
    call unit_test( test_eqv( d_graph_nih_ex, d_graph_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping_logical_2D

  subroutine test_mesh_vertices_to_graph_mapping_logical_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_mesh_vertices_to_graph_mapping_logical_3D'
    character(len=1024), parameter               :: test_name_local = 'logical_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    logical, dimension(:,:), allocatable         :: d_mesh_loc_ex
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                      :: vi, ni, d, k

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%vi1 :mesh%vi2, nz), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_V%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih,1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih,1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih,1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      do k = 1, nz
        d = ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))) + k
        d_mesh_loc_ex( vi,k) = modulo( d,2) == 0
        d_mesh_nih_ex( vi,k) = modulo( d,2) == 0
      end do
    end do
    do ni = graph%ni1, graph%ni2
      do k = 1, nz
        d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))) + k
        d_graph_nih_ex( ni,k) = modulo( d,2) == 0
      end do
    end do

    ! Map data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)
    call unit_test( test_eqv( d_graph_nih_ex, d_graph_nih_mapped), trim( test_name) // '/dist_hybrid')

    call map_mesh_vertices_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)
    call unit_test( test_eqv( d_graph_nih_ex, d_graph_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping_logical_3D

  subroutine test_mesh_vertices_to_graph_mapping_int_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_mesh_vertices_to_graph_mapping_int_2D'
    character(len=1024), parameter             :: test_name_local = 'int_2D'
    character(len=1024)                        :: test_name
    integer, dimension(:), allocatable         :: d_mesh_loc_ex
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    integer, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    integer, dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                    :: vi, ni, d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex     ( mesh%vi1 :mesh%vi2 ), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_V%n_nih)
    d_mesh_nih_ex( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      d = ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2)))
      d_mesh_loc_ex( vi) = d
      d_mesh_nih_ex( vi) = d
    end do
    do ni = graph%ni1, graph%ni2
      d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)))
      d_graph_nih_ex( ni) = d
    end do

    ! Map data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)
    call unit_test( test_eq( d_graph_nih_ex, d_graph_nih_mapped), trim( test_name) // '/dist_hybrid')

    call map_mesh_vertices_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)
    call unit_test( test_eq( d_graph_nih_ex, d_graph_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping_int_2D

  subroutine test_mesh_vertices_to_graph_mapping_int_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_mesh_vertices_to_graph_mapping_int_3D'
    character(len=1024), parameter               :: test_name_local = 'int_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    integer, dimension(:,:), allocatable         :: d_mesh_loc_ex
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                      :: vi, ni, d, k

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%vi1 :mesh%vi2, nz), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_V%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih,1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih,1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih,1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      do k = 1, nz
        d = ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))) + k
        d_mesh_loc_ex( vi,k) = d
        d_mesh_nih_ex( vi,k) = d
      end do
    end do
    do ni = graph%ni1, graph%ni2
      do k = 1, nz
        d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))) + k
        d_graph_nih_ex( ni,k) = d
      end do
    end do

    ! Map data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)
    call unit_test( test_eq( d_graph_nih_ex, d_graph_nih_mapped), trim( test_name) // '/dist_hybrid')

    call map_mesh_vertices_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)
    call unit_test( test_eq( d_graph_nih_ex, d_graph_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping_int_3D

  subroutine test_mesh_vertices_to_graph_mapping_dp_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_mesh_vertices_to_graph_mapping_dp_2D'
    character(len=1024), parameter              :: test_name_local = 'dp_2D'
    character(len=1024)                         :: test_name
    real(dp), dimension(:), allocatable         :: d_mesh_loc_ex
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                               :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                     :: vi, ni
    real(dp)                                    :: d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex     ( mesh%vi1 :mesh%vi2 ), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_V%n_nih)
    d_mesh_nih_ex( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      d = 1000._dp + hypot( mesh%V( vi,1), mesh%V( vi,2))
      d_mesh_loc_ex( vi) = d
      d_mesh_nih_ex( vi) = d
    end do
    do ni = graph%ni1, graph%ni2
      d = 1000._dp + hypot( graph%V( ni,1), graph%V( ni,2))
      d_graph_nih_ex( ni) = d
    end do

    ! Map data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)
    call unit_test( test_tol( d_graph_nih_ex, d_graph_nih_mapped, 1e-12_dp), trim( test_name) // '/dist_hybrid')

    call map_mesh_vertices_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)
    call unit_test( test_tol( d_graph_nih_ex, d_graph_nih_mapped, 1e-12_dp), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping_dp_2D

  subroutine test_mesh_vertices_to_graph_mapping_dp_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_mesh_vertices_to_graph_mapping_dp_3D'
    character(len=1024), parameter                :: test_name_local = 'dp_3D'
    character(len=1024)                           :: test_name
    integer                                       :: nz
    real(dp), dimension(:,:), allocatable         :: d_mesh_loc_ex
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                 :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                       :: vi, ni, k
    real(dp)                                      :: d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%vi1 :mesh%vi2, nz), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_V%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih,1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih,1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih,1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      do k = 1, nz
        d = 1000._dp + hypot( mesh%V( vi,1), mesh%V( vi,2)) + real( k,dp)
        d_mesh_loc_ex( vi,k) = d
        d_mesh_nih_ex( vi,k) = d
      end do
    end do
    do ni = graph%ni1, graph%ni2
      do k = 1, nz
        d = 1000._dp + hypot( graph%V( ni,1), graph%V( ni,2)) + real( k,dp)
        d_graph_nih_ex( ni,k) = d
      end do
    end do

    ! Map data from the mesh to the graph
    call map_mesh_vertices_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)
    call unit_test( test_eq( d_graph_nih_ex, d_graph_nih_mapped), trim( test_name) // '/dist_hybrid')

    call map_mesh_vertices_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)
    call unit_test( test_eq( d_graph_nih_ex, d_graph_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_vertices_to_graph_mapping_dp_3D

  ! Graph to mesh vertices
  ! ======================

  subroutine test_graph_to_mesh_vertices_mapping( test_name_parent, mesh, graph, mask_a)

    ! In/output variables:
    character(len=*),      intent(in) :: test_name_parent
    type(type_mesh),       intent(in) :: mesh
    type(type_graph),      intent(in) :: graph
    logical, dimension(:), intent(in) :: mask_a

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_to_mesh_vertices_mapping'
    character(len=1024), parameter :: test_name_local = 'graph_to_mesh_vertices'
    character(len=1024)            :: test_name
    logical, dimension(mesh%nV)    :: mask_a_tot

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call gather_to_all( mask_a, mask_a_tot)

    call test_graph_to_mesh_vertices_mapping_logical_2D( test_name, mesh, graph, mask_a_tot)
    call test_graph_to_mesh_vertices_mapping_logical_3D( test_name, mesh, graph, mask_a_tot)
    call test_graph_to_mesh_vertices_mapping_int_2D    ( test_name, mesh, graph, mask_a_tot)
    call test_graph_to_mesh_vertices_mapping_int_3D    ( test_name, mesh, graph, mask_a_tot)
    call test_graph_to_mesh_vertices_mapping_dp_2D     ( test_name, mesh, graph, mask_a_tot)
    call test_graph_to_mesh_vertices_mapping_dp_3D     ( test_name, mesh, graph, mask_a_tot)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping

  subroutine test_graph_to_mesh_vertices_mapping_logical_2D( test_name_parent, mesh, graph, mask_a_tot)

    ! In/output variables:
    character(len=*),            intent(in) :: test_name_parent
    type(type_mesh),             intent(in) :: mesh
    type(type_graph),            intent(in) :: graph
    logical, dimension(mesh%nV), intent(in) :: mask_a_tot

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_graph_to_mesh_vertices_mapping_logical_2D'
    character(len=1024), parameter             :: test_name_local = 'logical_2D'
    character(len=1024)                        :: test_name
    logical, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    type(MPI_WIN)                              :: wd_graph_nih_ex
    logical, dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_mesh_nih_mapped
    integer                                    :: vi, ni, d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%vi1 :mesh%vi2 ), source = .false.)
    allocate( d_mesh_loc_mapped( mesh%vi1 :mesh%vi2 ), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_V%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_V%n_nih)
    d_mesh_nih_ex    ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      if (mask_a_tot( vi)) then
        d = ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2)))
        d_mesh_loc_ex( vi) = modulo( d,2) == 0
        d_mesh_nih_ex( vi) = modulo( d,2) == 0
      end if
    end do
    do ni = graph%ni1, graph%ni2
      d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)))
      d_graph_nih_ex( ni) = modulo( d,2) == 0
    end do

    ! ! Map data from the mesh to the graph
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)
    call unit_test( test_eqv( d_mesh_loc_ex, d_mesh_loc_mapped), trim( test_name) // '/hybrid_dist')

    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)
    call unit_test( test_eqv( d_mesh_nih_ex, d_mesh_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping_logical_2D

  subroutine test_graph_to_mesh_vertices_mapping_logical_3D( test_name_parent, mesh, graph, mask_a_tot)

    ! In/output variables:
    character(len=*),            intent(in) :: test_name_parent
    type(type_mesh),             intent(in) :: mesh
    type(type_graph),            intent(in) :: graph
    logical, dimension(mesh%nV), intent(in) :: mask_a_tot

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_graph_to_mesh_vertices_mapping_logical_3D'
    character(len=1024), parameter               :: test_name_local = 'logical_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    type(MPI_WIN)                                :: wd_graph_nih_ex
    logical, dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_mesh_nih_mapped
    integer                                      :: vi, k, ni, d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%vi1 :mesh%vi2, 1:nz), source = .false.)
    allocate( d_mesh_loc_mapped( mesh%vi1 :mesh%vi2, 1:nz), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_V%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_V%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      if (mask_a_tot( vi)) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))) + k
          d_mesh_loc_ex( vi,k) = modulo( d,2) == 0
          d_mesh_nih_ex( vi,k) = modulo( d,2) == 0
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      do k = 1, nz
        d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))) + k
        d_graph_nih_ex( ni,k) = modulo( d,2) == 0
      end do
    end do

    ! ! Map data from the mesh to the graph
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)
    call unit_test( test_eqv( d_mesh_loc_ex, d_mesh_loc_mapped), trim( test_name) // '/hybrid_dist')

    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)
    call unit_test( test_eqv( d_mesh_nih_ex, d_mesh_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping_logical_3D

  subroutine test_graph_to_mesh_vertices_mapping_int_2D( test_name_parent, mesh, graph, mask_a_tot)

    ! In/output variables:
    character(len=*),            intent(in) :: test_name_parent
    type(type_mesh),             intent(in) :: mesh
    type(type_graph),            intent(in) :: graph
    logical, dimension(mesh%nV), intent(in) :: mask_a_tot

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_graph_to_mesh_vertices_mapping_int_2D'
    character(len=1024), parameter             :: test_name_local = 'int_2D'
    character(len=1024)                        :: test_name
    integer, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    type(MPI_WIN)                              :: wd_graph_nih_ex
    integer, dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_mesh_nih_mapped
    integer                                    :: vi, ni, d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%vi1 :mesh%vi2 ), source = 0)
    allocate( d_mesh_loc_mapped( mesh%vi1 :mesh%vi2 ), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_V%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_V%n_nih)
    d_mesh_nih_ex    ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      if (mask_a_tot( vi)) then
        d = ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2)))
        d_mesh_loc_ex( vi) = d
        d_mesh_nih_ex( vi) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)))
      d_graph_nih_ex( ni) = d
    end do

    ! ! Map data from the mesh to the graph
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)
    call unit_test( test_eq( d_mesh_loc_ex, d_mesh_loc_mapped), trim( test_name) // '/hybrid_dist')

    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)
    call unit_test( test_eq( d_mesh_nih_ex, d_mesh_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping_int_2D

  subroutine test_graph_to_mesh_vertices_mapping_int_3D( test_name_parent, mesh, graph, mask_a_tot)

    ! In/output variables:
    character(len=*),            intent(in) :: test_name_parent
    type(type_mesh),             intent(in) :: mesh
    type(type_graph),            intent(in) :: graph
    logical, dimension(mesh%nV), intent(in) :: mask_a_tot

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_graph_to_mesh_vertices_mapping_int_3D'
    character(len=1024), parameter               :: test_name_local = 'int_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    type(MPI_WIN)                                :: wd_graph_nih_ex
    integer, dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_mesh_nih_mapped
    integer                                      :: vi, k, ni, d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%vi1 :mesh%vi2, 1:nz), source = 0)
    allocate( d_mesh_loc_mapped( mesh%vi1 :mesh%vi2, 1:nz), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_V%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_V%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      if (mask_a_tot( vi)) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))) + k
          d_mesh_loc_ex( vi,k) = d
          d_mesh_nih_ex( vi,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      do k = 1, nz
        d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))) + k
        d_graph_nih_ex( ni,k) = d
      end do
    end do

    ! ! Map data from the mesh to the graph
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)
    call unit_test( test_eq( d_mesh_loc_ex, d_mesh_loc_mapped), trim( test_name) // '/hybrid_dist')

    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)
    call unit_test( test_eq( d_mesh_nih_ex, d_mesh_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping_int_3D

  subroutine test_graph_to_mesh_vertices_mapping_dp_2D( test_name_parent, mesh, graph, mask_a_tot)

    ! In/output variables:
    character(len=*),            intent(in) :: test_name_parent
    type(type_mesh),             intent(in) :: mesh
    type(type_graph),            intent(in) :: graph
    logical, dimension(mesh%nV), intent(in) :: mask_a_tot

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_graph_to_mesh_vertices_mapping_dp_2D'
    character(len=1024), parameter              :: test_name_local = 'dp_2D'
    character(len=1024)                         :: test_name
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    type(MPI_WIN)                               :: wd_graph_nih_ex
    real(dp), dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    type(MPI_WIN)                               :: wd_mesh_nih_ex, wd_mesh_nih_mapped
    integer                                     :: vi, ni
    real(dp)                                    :: d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%vi1 :mesh%vi2 ), source = 0._dp)
    allocate( d_mesh_loc_mapped( mesh%vi1 :mesh%vi2 ), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_V%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_V%n_nih)
    d_mesh_nih_ex    ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      if (mask_a_tot( vi)) then
        d = 1000._dp * hypot( mesh%V( vi,1), mesh%V( vi,2))
        d_mesh_loc_ex( vi) = d
        d_mesh_nih_ex( vi) = d
      end if
    end do
    do ni = graph%ni1, graph%ni2
      d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))
      d_graph_nih_ex( ni) = d
    end do

    ! ! Map data from the mesh to the graph
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)
    call unit_test( test_eq( d_mesh_loc_ex, d_mesh_loc_mapped), trim( test_name) // '/hybrid_dist')

    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)
    call unit_test( test_eq( d_mesh_nih_ex, d_mesh_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping_dp_2D

  subroutine test_graph_to_mesh_vertices_mapping_dp_3D( test_name_parent, mesh, graph, mask_a_tot)

    ! In/output variables:
    character(len=*),            intent(in) :: test_name_parent
    type(type_mesh),             intent(in) :: mesh
    type(type_graph),            intent(in) :: graph
    logical, dimension(mesh%nV), intent(in) :: mask_a_tot

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_graph_to_mesh_vertices_mapping_dp_3D'
    character(len=1024), parameter                :: test_name_local = 'dp_3D'
    character(len=1024)                           :: test_name
    integer                                       :: nz
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    type(MPI_WIN)                                 :: wd_graph_nih_ex
    real(dp), dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    type(MPI_WIN)                                 :: wd_mesh_nih_ex, wd_mesh_nih_mapped
    integer                                       :: vi, k, ni
    real(dp)                                      :: d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%vi1 :mesh%vi2, 1:nz), source = 0._dp)
    allocate( d_mesh_loc_mapped( mesh%vi1 :mesh%vi2, 1:nz), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_V%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_V%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_V%i1_nih: mesh%pai_V%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    do vi = mesh%vi1, mesh%vi2
      if (mask_a_tot( vi)) then
        do k = 1, nz
          d = 1000._dp + hypot( mesh%V( vi,1), mesh%V( vi,2)) + real( k,dp)
          d_mesh_loc_ex( vi,k) = d
          d_mesh_nih_ex( vi,k) = d
        end do
      end if
    end do
    do ni = graph%ni1, graph%ni2
      do k = 1, nz
        d = 1000._dp + hypot( graph%V( ni,1), graph%V( ni,2)) + real( k,dp)
        d_graph_nih_ex( ni,k) = d
      end do
    end do

    ! ! Map data from the mesh to the graph
    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)
    call unit_test( test_eq( d_mesh_loc_ex, d_mesh_loc_mapped), trim( test_name) // '/hybrid_dist')

    call map_graph_to_mesh_vertices( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)
    call unit_test( test_eq( d_mesh_nih_ex, d_mesh_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_vertices_mapping_dp_3D

  ! Mesh triangles to graph
  ! =======================

  subroutine test_mesh_triangles_to_graph_mapping( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_mesh_triangles_to_graph_mapping'
    character(len=1024), parameter :: test_name_local = 'mesh_triangles_to_graph'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_mesh_triangles_to_graph_mapping_logical_2D( test_name, mesh, graph)
    call test_mesh_triangles_to_graph_mapping_logical_3D( test_name, mesh, graph)
    call test_mesh_triangles_to_graph_mapping_int_2D    ( test_name, mesh, graph)
    call test_mesh_triangles_to_graph_mapping_int_3D    ( test_name, mesh, graph)
    call test_mesh_triangles_to_graph_mapping_dp_2D     ( test_name, mesh, graph)
    call test_mesh_triangles_to_graph_mapping_dp_3D     ( test_name, mesh, graph)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping

  subroutine test_mesh_triangles_to_graph_mapping_logical_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_mesh_triangles_to_graph_mapping_logical_2D'
    character(len=1024), parameter             :: test_name_local = 'logical_2D'
    character(len=1024)                        :: test_name
    logical, dimension(:), allocatable         :: d_mesh_loc_ex
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    logical, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    logical, dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                    :: ti, ni, d, mi
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex( mesh%ti1 :mesh%ti2), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_Tri%n_nih)
    d_mesh_nih_ex( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      if (graph%mi2ni( ti) > 0) then
        d = ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2)))
        d_mesh_loc_ex( ti) = modulo( d,2) == 0
        d_mesh_nih_ex( ti) = modulo( d,2) == 0
      else
        d_mesh_loc_ex( ti) = .false.
        d_mesh_nih_ex( ti) = .false.
      end if
    end do
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)))
        d_graph_nih_ex( ni) = modulo( d,2) == 0
      else
        d_graph_nih_ex( ni) = .false.
      end if
    end do

    ! Map (distributed) data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        test_result = test_result .and. d_graph_nih_mapped( ni) .eqv. d_graph_nih_ex( ni)
      end if
    end do
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map (hybrid distributed/shared) data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        test_result = test_result .and. d_graph_nih_mapped( ni) .eqv. d_graph_nih_ex( ni)
      end if
    end do
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping_logical_2D

  subroutine test_mesh_triangles_to_graph_mapping_logical_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_mesh_triangles_to_graph_mapping_logical_3D'
    character(len=1024), parameter               :: test_name_local = 'logical_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    logical, dimension(:,:), allocatable         :: d_mesh_loc_ex
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                    :: ti, ni, d, k, mi
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%ti1 :mesh%ti2, nz), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_Tri%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih,1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih,1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih,1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      do k = 1, nz
        d = ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))) + k
        d_mesh_loc_ex( ti,k) = modulo( d,2) == 0
        d_mesh_nih_ex( ti,k) = modulo( d,2) == 0
      end do
    end do
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))) + k
          d_graph_nih_ex( ni,k) = modulo( d,2) == 0
        end do
      else
        d_graph_nih_ex( ni,:) = .false.
      end if
    end do

    ! Map (distributed) data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_mapped( ni,k) .eqv. d_graph_nih_ex( ni,k)
        end do
      end if
    end do
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map (hybrid distributed/shared) data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_mapped( ni,k) .eqv. d_graph_nih_ex( ni,k)
        end do
      end if
    end do
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping_logical_3D

  subroutine test_mesh_triangles_to_graph_mapping_int_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_mesh_triangles_to_graph_mapping_int_2D'
    character(len=1024), parameter             :: test_name_local = 'int_2D'
    character(len=1024)                        :: test_name
    integer, dimension(:), allocatable         :: d_mesh_loc_ex
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    integer, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    integer, dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                    :: ti, ni, d, mi
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex( mesh%ti1 :mesh%ti2), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_Tri%n_nih)
    d_mesh_nih_ex( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      if (graph%mi2ni( ti) > 0) then
        d = ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2)))
        d_mesh_loc_ex( ti) = d
        d_mesh_nih_ex( ti) = d
      else
        d_mesh_loc_ex( ti) = 0
        d_mesh_nih_ex( ti) = 0
      end if
    end do
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)))
        d_graph_nih_ex( ni) = d
      else
        d_graph_nih_ex( ni) = 0
      end if
    end do

    ! Map (distributed) data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        test_result = test_result .and. d_graph_nih_mapped( ni) == d_graph_nih_ex( ni)
      end if
    end do
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map (hybrid distributed/shared) data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        test_result = test_result .and. d_graph_nih_mapped( ni) == d_graph_nih_ex( ni)
      end if
    end do
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping_int_2D

  subroutine test_mesh_triangles_to_graph_mapping_int_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_mesh_triangles_to_graph_mapping_int_3D'
    character(len=1024), parameter               :: test_name_local = 'int_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    integer, dimension(:,:), allocatable         :: d_mesh_loc_ex
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                    :: ti, ni, d, k, mi
    logical                                    :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%ti1 :mesh%ti2, nz), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_Tri%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih,1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih,1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih,1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      do k = 1, nz
        d = ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))) + k
        d_mesh_loc_ex( ti,k) = d
        d_mesh_nih_ex( ti,k) = d
      end do
    end do
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        do k = 1, nz
          d = ceiling( 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))) + k
          d_graph_nih_ex( ni,k) = d
        end do
      else
        d_graph_nih_ex( ni,:) = 0
      end if
    end do

    ! Map (distributed) data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_mapped( ni,k) == d_graph_nih_ex( ni,k)
        end do
      end if
    end do
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map (hybrid distributed/shared) data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_mapped( ni,k) == d_graph_nih_ex( ni,k)
        end do
      end if
    end do
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping_int_3D

  subroutine test_mesh_triangles_to_graph_mapping_dp_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_mesh_triangles_to_graph_mapping_dp_2D'
    character(len=1024), parameter              :: test_name_local = 'dp_2D'
    character(len=1024)                         :: test_name
    real(dp), dimension(:), allocatable         :: d_mesh_loc_ex
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_ex      => null()
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                               :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                     :: ti, ni, mi
    real(dp)                                    :: d
    logical                                     :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex( mesh%ti1 :mesh%ti2), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_Tri%n_nih)
    d_mesh_nih_ex( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      if (graph%mi2ni( ti) > 0) then
        d = 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))
        d_mesh_loc_ex( ti) = d
        d_mesh_nih_ex( ti) = d
      else
        d_mesh_loc_ex( ti) = 0._dp
        d_mesh_nih_ex( ti) = 0._dp
      end if
    end do
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2))
        d_graph_nih_ex( ni) = d
      else
        d_graph_nih_ex( ni) = 0._dp
      end if
    end do

    ! Map (distributed) data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        test_result = test_result .and. d_graph_nih_mapped( ni) == d_graph_nih_ex( ni)
      end if
    end do
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map (hybrid distributed/shared) data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        test_result = test_result .and. d_graph_nih_mapped( ni) == d_graph_nih_ex( ni)
      end if
    end do
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping_dp_2D

  subroutine test_mesh_triangles_to_graph_mapping_dp_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_mesh_triangles_to_graph_mapping_dp_3D'
    character(len=1024), parameter                :: test_name_local = 'dp_3D'
    character(len=1024)                           :: test_name
    integer                                       :: nz
    real(dp), dimension(:,:), allocatable         :: d_mesh_loc_ex
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_ex      => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_mapped => null()
    type(MPI_WIN)                                 :: wd_mesh_nih_ex, wd_graph_nih_ex, wd_graph_nih_mapped
    integer                                       :: ti, ni, k, mi
    real(dp)                                      :: d
    logical                                       :: test_result

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex( mesh%ti1 :mesh%ti2, nz), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex, wd_mesh_nih_ex, mesh%pai_Tri%n_nih, nz)
    d_mesh_nih_ex( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih,1:nz) => d_mesh_nih_ex

    call allocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    , graph%pai%n_nih, nz)
    call allocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped, graph%pai%n_nih, nz)
    d_graph_nih_ex    ( graph%pai%i1_nih: graph%pai%i2_nih,1:nz) => d_graph_nih_ex
    d_graph_nih_mapped( graph%pai%i1_nih: graph%pai%i2_nih,1:nz) => d_graph_nih_mapped

    ! Fill in exact solutions
    do ti = mesh%ti1, mesh%ti2
      do k = 1, nz
        d = 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2)) + k
        d_mesh_loc_ex( ti,k) = d
        d_mesh_nih_ex( ti,k) = d
      end do
    end do
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        do k = 1, nz
          d = 1000._dp * hypot( graph%V( ni,1), graph%V( ni,2)) + k
          d_graph_nih_ex( ni,k) = d
        end do
      else
        d_graph_nih_ex( ni,:) = 0._dp
      end if
    end do

    ! Map (distributed) data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_loc_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_mapped( ni,k) == d_graph_nih_ex( ni,k)
        end do
      end if
    end do
    call unit_test( test_result, trim( test_name) // '/dist_hybrid')

    ! Map (hybrid distributed/shared) data from the mesh to the graph
    call map_mesh_triangles_to_graph( mesh, d_mesh_nih_ex, graph, d_graph_nih_mapped)

    test_result = .true.
    do ni = graph%ni1, graph%ni2
      if (.not. graph%is_ghost( ni)) then
        do k = 1, nz
          test_result = test_result .and. d_graph_nih_mapped( ni,k) == d_graph_nih_ex( ni,k)
        end do
      end if
    end do
    call unit_test( test_result, trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_mesh_nih_ex     , wd_mesh_nih_ex     )
    call deallocate_dist_shared( d_graph_nih_ex    , wd_graph_nih_ex    )
    call deallocate_dist_shared( d_graph_nih_mapped, wd_graph_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_mesh_triangles_to_graph_mapping_dp_3D

  ! Graph to mesh triangles
  ! =======================

  subroutine test_graph_to_mesh_triangles_mapping( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_graph_to_mesh_triangles_mapping'
    character(len=1024), parameter :: test_name_local = 'graph_to_mesh_triangles'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_graph_to_mesh_triangles_mapping_logical_2D( test_name, mesh, graph)
    call test_graph_to_mesh_triangles_mapping_logical_3D( test_name, mesh, graph)
    call test_graph_to_mesh_triangles_mapping_int_2D    ( test_name, mesh, graph)
    call test_graph_to_mesh_triangles_mapping_int_3D    ( test_name, mesh, graph)
    call test_graph_to_mesh_triangles_mapping_dp_2D     ( test_name, mesh, graph)
    call test_graph_to_mesh_triangles_mapping_dp_3D     ( test_name, mesh, graph)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping

  subroutine test_graph_to_mesh_triangles_mapping_logical_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_graph_to_mesh_triangles_mapping_logical_2D'
    character(len=1024), parameter             :: test_name_local = 'logical_2D'
    character(len=1024)                        :: test_name
    logical, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    type(MPI_WIN)                              :: wd_graph_nih_ex
    logical, dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    logical, dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_mesh_nih_mapped
    integer                                    :: ti, ni, d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2 ), source = .false.)
    allocate( d_mesh_loc_mapped( mesh%ti1 :mesh%ti2 ), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_Tri%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_Tri%n_nih)
    d_mesh_nih_ex    ( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2 ) = .false.
    d_mesh_loc_mapped( mesh%ti1 :mesh%ti2 ) = .false.
    d_mesh_nih_ex    ( mesh%ti1 :mesh%ti2 ) = .false.
    d_mesh_nih_mapped( mesh%ti1 :mesh%ti2 ) = .false.
    d_graph_nih_ex   ( graph%ni1:graph%ni2) = .false.

    do ni = 1, graph%n
      if (.not. graph%is_ghost( ni)) then

        ti = graph%ni2mi( ni)
        d = ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2)))

        ! Graph data
        if (ni >= graph%ni1 .and. ni <= graph%ni2) then
          d_graph_nih_ex( ni) = modulo( d,2) == 0
        end if

        ! Mesh data
        if (ti >= mesh%ti1 .and. ti <= mesh%ti2) then
          d_mesh_loc_ex( ti) = modulo( d,2) == 0
          d_mesh_nih_ex( ti) = modulo( d,2) == 0
        end if

      end if
    end do

    ! ! Map data from the mesh to the graph
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)
    call unit_test( test_eqv( d_mesh_loc_ex, d_mesh_loc_mapped), trim( test_name) // '/hybrid_dist')

    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)
    call unit_test( test_eqv( d_mesh_nih_ex, d_mesh_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping_logical_2D

  subroutine test_graph_to_mesh_triangles_mapping_logical_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_graph_to_mesh_triangles_mapping_logical_3D'
    character(len=1024), parameter               :: test_name_local = 'logical_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    logical, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    type(MPI_WIN)                                :: wd_graph_nih_ex
    logical, dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    logical, dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_mesh_nih_mapped
    integer                                      :: ti, ni, k, d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2, 1:nz), source = .false.)
    allocate( d_mesh_loc_mapped( mesh%ti1 :mesh%ti2, 1:nz), source = .false.)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_Tri%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_Tri%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2 ,:) = .false.
    d_mesh_loc_mapped( mesh%ti1 :mesh%ti2 ,:) = .false.
    d_mesh_nih_ex    ( mesh%ti1 :mesh%ti2 ,:) = .false.
    d_mesh_nih_mapped( mesh%ti1 :mesh%ti2 ,:) = .false.
    d_graph_nih_ex   ( graph%ni1:graph%ni2,:) = .false.

    do ni = 1, graph%n
      if (.not. graph%is_ghost( ni)) then
        do k = 1, nz

          ti = graph%ni2mi( ni)
          d = ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))) + k

          ! Graph data
          if (ni >= graph%ni1 .and. ni <= graph%ni2) then
            d_graph_nih_ex( ni,k) = modulo( d,2) == 0
          end if

          ! Mesh data
          if (ti >= mesh%ti1 .and. ti <= mesh%ti2) then
            d_mesh_loc_ex( ti,k) = modulo( d,2) == 0
            d_mesh_nih_ex( ti,k) = modulo( d,2) == 0
          end if

        end do
      end if
    end do

    ! ! Map data from the mesh to the graph
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)
    call unit_test( test_eqv( d_mesh_loc_ex, d_mesh_loc_mapped), trim( test_name) // '/hybrid_dist')

    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)
    call unit_test( test_eqv( d_mesh_nih_ex, d_mesh_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping_logical_3D

  subroutine test_graph_to_mesh_triangles_mapping_int_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter             :: routine_name = 'test_graph_to_mesh_triangles_mapping_int_2D'
    character(len=1024), parameter             :: test_name_local = 'int_2D'
    character(len=1024)                        :: test_name
    integer, dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    type(MPI_WIN)                              :: wd_graph_nih_ex
    integer, dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    integer, dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    type(MPI_WIN)                              :: wd_mesh_nih_ex, wd_mesh_nih_mapped
    integer                                    :: ti, ni, d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2 ), source = 0)
    allocate( d_mesh_loc_mapped( mesh%ti1 :mesh%ti2 ), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_Tri%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_Tri%n_nih)
    d_mesh_nih_ex    ( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2 ) = 0
    d_mesh_loc_mapped( mesh%ti1 :mesh%ti2 ) = 0
    d_mesh_nih_ex    ( mesh%ti1 :mesh%ti2 ) = 0
    d_mesh_nih_mapped( mesh%ti1 :mesh%ti2 ) = 0
    d_graph_nih_ex   ( graph%ni1:graph%ni2) = 0

    do ni = 1, graph%n
      if (.not. graph%is_ghost( ni)) then

        ti = graph%ni2mi( ni)
        d = ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2)))

        ! Graph data
        if (ni >= graph%ni1 .and. ni <= graph%ni2) then
          d_graph_nih_ex( ni) = d
        end if

        ! Mesh data
        if (ti >= mesh%ti1 .and. ti <= mesh%ti2) then
          d_mesh_loc_ex( ti) = d
          d_mesh_nih_ex( ti) = d
        end if

      end if
    end do

    ! ! Map data from the mesh to the graph
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)
    call unit_test( test_eq( d_mesh_loc_ex, d_mesh_loc_mapped), trim( test_name) // '/hybrid_dist')

    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)
    call unit_test( test_eq( d_mesh_nih_ex, d_mesh_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping_int_2D

  subroutine test_graph_to_mesh_triangles_mapping_int_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter               :: routine_name = 'test_graph_to_mesh_triangles_mapping_int_3D'
    character(len=1024), parameter               :: test_name_local = 'int_3D'
    character(len=1024)                          :: test_name
    integer                                      :: nz
    integer, dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    type(MPI_WIN)                                :: wd_graph_nih_ex
    integer, dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    integer, dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    type(MPI_WIN)                                :: wd_mesh_nih_ex, wd_mesh_nih_mapped
    integer                                      :: ti, ni, k, d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2, 1:nz), source = 0)
    allocate( d_mesh_loc_mapped( mesh%ti1 :mesh%ti2, 1:nz), source = 0)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_Tri%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_Tri%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2 ,:) = 0
    d_mesh_loc_mapped( mesh%ti1 :mesh%ti2 ,:) = 0
    d_mesh_nih_ex    ( mesh%ti1 :mesh%ti2 ,:) = 0
    d_mesh_nih_mapped( mesh%ti1 :mesh%ti2 ,:) = 0
    d_graph_nih_ex   ( graph%ni1:graph%ni2,:) = 0

    do ni = 1, graph%n
      if (.not. graph%is_ghost( ni)) then
        do k = 1, nz

          ti = graph%ni2mi( ni)
          d = ceiling( 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))) + k

          ! Graph data
          if (ni >= graph%ni1 .and. ni <= graph%ni2) then
            d_graph_nih_ex( ni,k) = d
          end if

          ! Mesh data
          if (ti >= mesh%ti1 .and. ti <= mesh%ti2) then
            d_mesh_loc_ex( ti,k) = d
            d_mesh_nih_ex( ti,k) = d
          end if

        end do
      end if
    end do

    ! ! Map data from the mesh to the graph
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)
    call unit_test( test_eq( d_mesh_loc_ex, d_mesh_loc_mapped), trim( test_name) // '/hybrid_dist')

    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)
    call unit_test( test_eq( d_mesh_nih_ex, d_mesh_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping_int_3D

  subroutine test_graph_to_mesh_triangles_mapping_dp_2D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter              :: routine_name = 'test_graph_to_mesh_triangles_mapping_dp_2D'
    character(len=1024), parameter              :: test_name_local = 'dp_2D'
    character(len=1024)                         :: test_name
    real(dp), dimension(:), contiguous, pointer :: d_graph_nih_ex     => null()
    type(MPI_WIN)                               :: wd_graph_nih_ex
    real(dp), dimension(:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_ex     => null()
    real(dp), dimension(:), contiguous, pointer :: d_mesh_nih_mapped => null()
    type(MPI_WIN)                               :: wd_mesh_nih_ex, wd_mesh_nih_mapped
    integer                                     :: ti, ni, d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    allocate( d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2 ), source = 0._dp)
    allocate( d_mesh_loc_mapped( mesh%ti1 :mesh%ti2 ), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_Tri%n_nih)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_Tri%n_nih)
    d_mesh_nih_ex    ( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih) => d_graph_nih_ex

    ! Fill in exact solutions
    d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2 ) = 0._dp
    d_mesh_loc_mapped( mesh%ti1 :mesh%ti2 ) = 0._dp
    d_mesh_nih_ex    ( mesh%ti1 :mesh%ti2 ) = 0._dp
    d_mesh_nih_mapped( mesh%ti1 :mesh%ti2 ) = 0._dp
    d_graph_nih_ex   ( graph%ni1:graph%ni2) = 0._dp

    do ni = 1, graph%n
      if (.not. graph%is_ghost( ni)) then

        ti = graph%ni2mi( ni)
        d = 1000._dp * hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2))

        ! Graph data
        if (ni >= graph%ni1 .and. ni <= graph%ni2) then
          d_graph_nih_ex( ni) = d
        end if

        ! Mesh data
        if (ti >= mesh%ti1 .and. ti <= mesh%ti2) then
          d_mesh_loc_ex( ti) = d
          d_mesh_nih_ex( ti) = d
        end if

      end if
    end do

    ! ! Map data from the mesh to the graph
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)
    call unit_test( test_eq( d_mesh_loc_ex, d_mesh_loc_mapped), trim( test_name) // '/hybrid_dist')

    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)
    call unit_test( test_eq( d_mesh_nih_ex, d_mesh_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping_dp_2D

  subroutine test_graph_to_mesh_triangles_mapping_dp_3D( test_name_parent, mesh, graph)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent
    type(type_mesh),  intent(in) :: mesh
    type(type_graph), intent(in) :: graph

    ! Local variables:
    character(len=1024), parameter                :: routine_name = 'test_graph_to_mesh_triangles_mapping_dp_3D'
    character(len=1024), parameter                :: test_name_local = 'dp_3D'
    character(len=1024)                           :: test_name
    integer                                       :: nz
    real(dp), dimension(:,:), contiguous, pointer :: d_graph_nih_ex     => null()
    type(MPI_WIN)                                 :: wd_graph_nih_ex
    real(dp), dimension(:,:), allocatable         :: d_mesh_loc_ex, d_mesh_loc_mapped
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_ex     => null()
    real(dp), dimension(:,:), contiguous, pointer :: d_mesh_nih_mapped => null()
    type(MPI_WIN)                                 :: wd_mesh_nih_ex, wd_mesh_nih_mapped
    integer                                       :: ti, ni, k, d

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    ! Allocate memory
    nz = 12
    allocate( d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2, 1:nz), source = 0._dp)
    allocate( d_mesh_loc_mapped( mesh%ti1 :mesh%ti2, 1:nz), source = 0._dp)

    ! Allocate hybrid distributed/shared memory
    call allocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    , mesh%pai_Tri%n_nih, nz)
    call allocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped, mesh%pai_Tri%n_nih, nz)
    d_mesh_nih_ex    ( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_ex
    d_mesh_nih_mapped( mesh%pai_Tri%i1_nih: mesh%pai_Tri%i2_nih, 1:nz) => d_mesh_nih_mapped

    call allocate_dist_shared( d_graph_nih_ex, wd_graph_nih_ex, graph%pai%n_nih, nz)
    d_graph_nih_ex( graph%pai%i1_nih: graph%pai%i2_nih, 1:nz) => d_graph_nih_ex

    ! Fill in exact solutions
    d_mesh_loc_ex    ( mesh%ti1 :mesh%ti2 ,:) = 0._dp
    d_mesh_loc_mapped( mesh%ti1 :mesh%ti2 ,:) = 0._dp
    d_mesh_nih_ex    ( mesh%ti1 :mesh%ti2 ,:) = 0._dp
    d_mesh_nih_mapped( mesh%ti1 :mesh%ti2 ,:) = 0._dp
    d_graph_nih_ex   ( graph%ni1:graph%ni2,:) = 0._dp

    do ni = 1, graph%n
      if (.not. graph%is_ghost( ni)) then
        do k = 1, nz

          ti = graph%ni2mi( ni)
          d = hypot( mesh%TriGC( ti,1), mesh%TriGC( ti,2)) + real( k,dp)

          ! Graph data
          if (ni >= graph%ni1 .and. ni <= graph%ni2) then
            d_graph_nih_ex( ni,k) = d
          end if

          ! Mesh data
          if (ti >= mesh%ti1 .and. ti <= mesh%ti2) then
            d_mesh_loc_ex( ti,k) = d
            d_mesh_nih_ex( ti,k) = d
          end if

        end do
      end if
    end do

    ! ! Map data from the mesh to the graph
    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_loc_mapped)
    call unit_test( test_eq( d_mesh_loc_ex, d_mesh_loc_mapped), trim( test_name) // '/hybrid_dist')

    call map_graph_to_mesh_triangles( graph, d_graph_nih_ex, mesh, d_mesh_nih_mapped)
    call unit_test( test_eq( d_mesh_nih_ex, d_mesh_nih_mapped), trim( test_name) // '/hybrid_hybrid')

    ! Clean up after yourself
    call deallocate_dist_shared( d_graph_nih_ex   , wd_graph_nih_ex   )
    call deallocate_dist_shared( d_mesh_nih_ex    , wd_mesh_nih_ex    )
    call deallocate_dist_shared( d_mesh_nih_mapped, wd_mesh_nih_mapped)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graph_to_mesh_triangles_mapping_dp_3D

end module ut_mesh_graphs_mapping