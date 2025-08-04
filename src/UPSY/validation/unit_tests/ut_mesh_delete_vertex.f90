module ut_mesh_delete_vertex

  use tests_main
  use assertions_basic
  use ut_basic
  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash, warning
  use mpi_basic, only: par
  use mesh_types, only: type_mesh
  use mesh_memory, only: allocate_mesh_primary, crop_mesh_primary
  use mesh_dummy_meshes, only: initialise_dummy_mesh_5
  use mesh_refinement_basic, only: refine_mesh_uniform
  use parameters, only: pi
  use mesh_secondary, only: calc_all_secondary_mesh_data
  use delete_vertices, only: delete_vertex

  use netcdf_io_main

  implicit none

  private

  public :: test_delete_vertex

contains

  subroutine test_delete_vertex( test_name_parent)
    ! Test the trace_line_grid_a subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'test_delete_vertex'
    character(len=1024), parameter :: test_name_local = 'delete_vertex'
    character(len=1024)            :: test_name

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    call test_delete_vertex_nV9( test_name)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_delete_vertex

  subroutine test_delete_vertex_nV9( test_name_parent)
    ! Test the trace_line_grid_a subroutine

    ! In/output variables:
    character(len=*), intent(in) :: test_name_parent

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'test_delete_vertex_nV9'
    character(len=1024), parameter     :: test_name_local = 'nV9'
    character(len=1024)                :: test_name
    logical                            :: verified
    real(dp), parameter                :: xmin = 0._dp
    real(dp), parameter                :: xmax = 1._dp
    real(dp), parameter                :: ymin = 0._dp
    real(dp), parameter                :: ymax = 1._dp
    real(dp)                           :: alpha_min
    real(dp)                           :: res_max
    type(type_mesh)                    :: mesh
    integer                            :: vi, ci, vj, vi_kill
    logical                            :: is_border_vertex, has_border_neighbours
    integer                            :: nV_before, nTri_before
    integer, dimension(:), allocatable :: vi_new2vi_old, vi_old2vi_new, ti_new2ti_old, ti_old2ti_new

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Add test name to list
    test_name = trim( test_name_parent) // '/' // trim( test_name_local)

    verified = .true.

    call allocate_mesh_primary( mesh, 'test_mesh', 100, 200)
    call initialise_dummy_mesh_5( mesh, xmin, xmax, ymin, ymax)

    ! Refine the test mesh
    alpha_min = 25._dp * pi / 180._dp
    res_max = pi / 20._dp
    call refine_mesh_uniform( mesh, res_max, alpha_min)
    call crop_mesh_primary( mesh)
    call calc_all_secondary_mesh_data( mesh, 0._dp, -90._dp, 71._dp)

    ! Find a nice vertex to kill
    ! (must not be a border vertex, nor be adjacent to any border vertices)
    vi_kill = 0
    do vi = 1, mesh%nV
      is_border_vertex = mesh%VBI( vi) > 0
      has_border_neighbours = .false.
      do ci = 1, mesh%nC( vi)
        vj = mesh%C( vi,ci)
        if (mesh%VBI( vj) > 0) then
          has_border_neighbours = .true.
          exit
        end if
      end do
      if (.not. is_border_vertex .and. .not. has_border_neighbours) then
        vi_kill = vi
        exit
      end if
    end do

    nV_before   = mesh%nV
    nTri_before = mesh%nTri

    call delete_vertex( mesh, vi_kill, vi_new2vi_old, vi_old2vi_new, ti_new2ti_old, ti_old2ti_new)

    verified = verified .and. test_mesh_is_self_consistent( mesh)
    verified = verified .and. mesh%nV   == nV_before   - 1
    verified = verified .and. mesh%nTri == ntri_before - 2

    call unit_test( verified, trim(test_name))

    call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine test_delete_vertex_nV9

end module ut_mesh_delete_vertex