module switch_vertices_triangles

  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use switch_array_elements, only: switch_rows, switch_entries_by_value

  implicit none

  private

  public :: switch_vertices

contains

  subroutine switch_vertices( mesh, vi, vj)
    !< Switch vertices vi and vj

    ! Input variables
    type(type_mesh), intent(inout) :: mesh
    integer,         intent(in   ) :: vi, vj

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'switch_vertices'

    ! Add routine to path
    call init_routine( routine_name)

    ! Switch rows in V, nC, C, niTri, iTri, VBI
    call switch_rows( mesh%V    , vi, vj)
    call switch_rows( mesh%nC   , vi, vj)
    call switch_rows( mesh%C    , vi, vj)
    call switch_rows( mesh%niTri, vi, vj)
    call switch_rows( mesh%iTri , vi, vj)
    call switch_rows( mesh%VBI  , vi, vj)

    ! Switch entries in C, Tri
    call switch_entries_by_value( mesh%C  , vi, vj)
    call switch_entries_by_value( mesh%Tri, vi, vj)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine switch_vertices

end module switch_vertices_triangles
