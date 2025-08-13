module ut_mesh_graphs

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mpi_basic, only: par
  use mesh_types, only: type_mesh
  use graph_types, only: type_graph
  use parameters, only: pi
  use mesh_memory, only: allocate_mesh_primary, crop_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use create_graphs_from_masked_mesh, only: create_graph_from_masked_mesh_b

  use netcdf_io_main

  implicit none

  private

  public :: test_graphs

contains

  subroutine test_graphs( test_name_parent)

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_graphs'
    character(len=1024), parameter     :: test_name_local = 'graphs'
    character(len=1024)                :: test_name
    real(dp), parameter                :: xmin = -1._dp
    real(dp), parameter                :: xmax =  1._dp
    real(dp), parameter                :: ymin = -1._dp
    real(dp), parameter                :: ymax =  1._dp
    real(dp)                           :: alpha_min
    real(dp)                           :: res_max
    type(type_mesh)                    :: mesh
    logical, dimension(:), allocatable :: mask_a
    integer                            :: vi
    type(type_graph)                   :: graph

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call allocate_mesh_primary( mesh, 'test_mesh', 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 23.2_dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call crop_mesh_primary( mesh)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)

    ! Create a graph from the mesh triangles
    allocate( mask_a( mesh%vi1:mesh%vi2), source = .false.)
    do vi = mesh%vi1, mesh%vi2
      if (hypot( mesh%V( vi,1), mesh%V( vi,2)) < mesh%xmax * 0.5_dp) mask_a( vi) = .true.
    end do
    call create_graph_from_masked_mesh_b( mesh, mask_a, graph)

    call unit_test( .true., trim( test_name))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_graphs

end module ut_mesh_graphs