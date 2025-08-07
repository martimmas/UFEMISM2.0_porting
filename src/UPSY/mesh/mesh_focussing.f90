module mesh_focussing

  use precisions, only: dp
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use polyline_types, only: type_polyline
  use remapping_types, only: type_single_row_mapping_matrices
  use line_tracing_Voronoi, only: trace_line_Vor
  use delete_vertices, only: delete_vertex
  use split_triangles, only: split_triangle
  use mesh_memory, only: extend_mesh_primary, crop_mesh_primary

  implicit none

  private

  public :: delete_vertices_along_polyline, focus_mesh_on_polyline

contains

  subroutine focus_mesh_on_polyline( mesh, ll)
    !< Focus a mesh on a polyline

    ! In/output variables:
    type(type_mesh),     intent(inout) :: mesh
    type(type_polyline), intent(in   ) :: ll

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'focus_mesh_on_polyline'

    ! Add routine to path
    call init_routine( routine_name)

    call delete_vertices_along_polyline( mesh, ll)
    call add_polyline_vertices_to_mesh( mesh, ll)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine focus_mesh_on_polyline

  subroutine delete_vertices_along_polyline( mesh, ll)
    !< Delete all vertices whose Voronoi cells are crossed by the polyline

    ! In/output variables:
    type(type_mesh),     intent(inout) :: mesh
    type(type_polyline), intent(in   ) :: ll

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'delete_vertices_along_polyline'
    integer, dimension(:), allocatable :: vi_list
    integer                            :: i, vi_kill, ii
    integer, dimension(:), allocatable :: vi_new2vi_old, vi_old2vi_new, ti_new2ti_old, ti_old2ti_new

    ! Add routine to path
    call init_routine( routine_name)

    call list_vertices_crossed_by_polyline( mesh, ll, vi_list)

    do i = 1, size( vi_list)

      ! Delete vertex
      vi_kill = vi_list( i)
      call delete_vertex( mesh, vi_kill, vi_new2vi_old, vi_old2vi_new, ti_new2ti_old, ti_old2ti_new)

      ! Update indices of remaining to-be-deleted vertices
      do ii = i+1, size( vi_list)
        vi_list( ii) = vi_old2vi_new( vi_list( ii))
      end do

      ! Clean up after yourself (not sure if memory is automatically deallocated when arrays are allocated in another procedure...)
      deallocate( vi_new2vi_old)
      deallocate( vi_old2vi_new)
      deallocate( ti_new2ti_old)
      deallocate( ti_old2ti_new)

    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine delete_vertices_along_polyline

  subroutine list_vertices_crossed_by_polyline( mesh, ll, vi_list)
    !< List all vertices whose Voronoi cells are crossed by the polyline

    ! In/output variables:
    type(type_mesh),                    intent(in   ) :: mesh
    type(type_polyline),                intent(in   ) :: ll
    integer, dimension(:), allocatable, intent(  out) :: vi_list

    ! Local variables:
    character(len=1024), parameter         :: routine_name = 'delete_vertices_along_polyline'
    integer                                :: vi_hint
    integer                                :: i_end, i1, i2
    real(dp), dimension(2)                 :: p,q
    type(type_single_row_mapping_matrices) :: single_row

    ! Add routine to path
    call init_routine( routine_name)

    if (ll%is_closed) then
      ! Include the line segment from the last to the first vertex
      i_end = ll%n
    else
      ! Don't
      i_end = ll%n-1
    end if

    ! Allocate memory for single row results
    single_row%n_max = mesh%nV
    single_row%n     = 0
    allocate( single_row%index_left( single_row%n_max))
    allocate( single_row%LI_xdy(     single_row%n_max))
    allocate( single_row%LI_mxydx(   single_row%n_max))
    allocate( single_row%LI_xydy(    single_row%n_max))

    vi_hint = 1
    do i1 = 1, i_end
      i2 = i1 + 1
      if (i2 == ll%n+1) i2 = 1

      p = ll%p( i1,:)
      q = ll%p( i2,:)

      call trace_line_Vor( mesh, p, q, single_row, .true., vi_hint)

    end do

    vi_list = single_row%index_left( 1: single_row%n)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine list_vertices_crossed_by_polyline

  subroutine add_polyline_vertices_to_mesh( mesh, ll)
    !< Add the vertices of the polyline to the mesh

    ! In/output variables:
    type(type_mesh),     intent(inout) :: mesh
    type(type_polyline), intent(in   ) :: ll

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_polyline_vertices_to_mesh'
    integer                        :: i, ti_in
    real(dp), dimension(2)         :: p

    ! Add routine to path
    call init_routine( routine_name)

    ! Make sure there's enough memory for the new vertices
    call extend_mesh_primary( mesh, mesh%nV + ll%n, mesh%nTri + ll%n*3)

    ti_in = 1
    do i = 1, ll%n
      p = ll%p( i,:)
      call split_triangle( mesh, ti_in, p)
    end do

    call crop_mesh_primary( mesh)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine add_polyline_vertices_to_mesh

end module mesh_focussing
