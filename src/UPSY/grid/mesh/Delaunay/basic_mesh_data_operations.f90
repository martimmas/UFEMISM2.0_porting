module basic_mesh_data_operations

  use precisions, only: dp
  use assertions_basic
  use mesh_types, only: type_mesh

  implicit none

  private

  public :: replace_vj_in_C_vi_with_vks
  public :: remove_vj_in_C_vi
  public :: remove_ti_in_iTri_vi
  public :: replace_ti_in_iTri_vi_with_tj
  public :: add_tjs_in_iTri_vi_after_ti
  public :: replace_vi_in_Tri_ti_with_vj
  public :: replace_tj_in_TriC_ti_with_tk

contains

  subroutine replace_vj_in_C_vi_with_vks( mesh, vi, vj, vks)
    ! Replace vj in C(vi,:) with vk(:)

    ! In/output variables:
    type(type_mesh),       intent(inout) :: mesh
    integer,               intent(in   ) :: vi, vj
    integer, dimension(:), intent(in   ) :: vks

    ! Local variables:
    logical :: found_it
    integer :: ci, i

    found_it = .false.
    do ci = 1, mesh%nC( vi)
      if (mesh%C( vi,ci) == vj) then
        found_it = .true.
        mesh%C( vi,1:mesh%nC( vi) - 1 + size( vks,1)) = [mesh%C( vi,1:ci-1), vks, mesh%C( vi,ci+1:mesh%nC( vi))]
        mesh%nC( vi) = mesh%nC( vi) - 1 + size( vks,1)
        exit
      end if
    end do
#if (DO_ASSERTIONS)
    call assert( found_it, 'couldnt find vj in C(vi,:)')
#endif

  end subroutine replace_vj_in_C_vi_with_vks

  subroutine remove_vj_in_C_vi( mesh, vi, vj)
    ! Remove vj in C(vi,:)

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh
    integer,         intent(in   ) :: vi, vj

    ! Local variables:
    logical :: found_it
    integer :: ci

    found_it = .false.
    do ci = 1, mesh%nC( vi)
      if (mesh%C( vi,ci) == vj) then
        found_it = .true.
        mesh%C( vi,1:mesh%nC( vi)) = [mesh%C( vi,1:ci-1), mesh%C( vi,ci+1:mesh%nC( vi)), 0]
        mesh%nC( vi) = mesh%nC( vi) - 1
        exit
      end if
    end do
#if (DO_ASSERTIONS)
    call assert( found_it, 'couldnt find vj in C(vi,:)')
#endif

  end subroutine remove_vj_in_C_vi

  subroutine remove_ti_in_iTri_vi( mesh, vi, ti)
    ! Remove ti in iTri( vi,:)

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh
    integer,         intent(in   ) :: vi, ti

    ! Local variables:
    logical :: found_it
    integer :: iti

    found_it = .false.
    do iti = 1, mesh%niTri( vi)
      if (mesh%iTri( vi,iti) == ti) then
        found_it = .true.
        mesh%iTri( vi,1:mesh%niTri( vi)) = [mesh%iTri( vi,1:iti-1), mesh%iTri( vi,iti+1:mesh%niTri( vi)), 0]
        mesh%niTri( vi) = mesh%niTri( vi) - 1
        exit
      end if
    end do
#if (DO_ASSERTIONS)
    call assert( found_it, 'couldnt find ti in iTri(vi,:)')
#endif

  end subroutine remove_ti_in_iTri_vi

  subroutine replace_ti_in_iTri_vi_with_tj( mesh, vi, ti, tj)
    ! Replace ti in iTri(vi,:) with tj

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh
    integer,         intent(in   ) :: vi, ti, tj

    ! Local variables:
    logical :: found_it
    integer :: iti

    found_it = .false.
    do iti = 1, mesh%niTri( vi)
      if (mesh%iTri( vi,iti) == ti) then
        found_it = .true.
        mesh%iTri( vi,iti) = tj
        exit
      end if
    end do
#if (DO_ASSERTIONS)
    call assert( found_it, 'couldnt find ti in iTri(vi,:)')
#endif

  end subroutine replace_ti_in_iTri_vi_with_tj

  subroutine add_tjs_in_iTri_vi_after_ti( mesh, vi, ti, tjs)
    ! Add tjs in iTri( vi,:) after tj

    ! In/output variables:
    type(type_mesh),       intent(inout) :: mesh
    integer,               intent(in   ) :: vi, ti
    integer, dimension(:), intent(in   ) :: tjs

    ! Local variables:
    logical :: found_it
    integer :: iti, i

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
#if (DO_ASSERTIONS)
    call assert( found_it, 'couldnt find ti in iTri(vi,:)')
#endif

  end subroutine add_tjs_in_iTri_vi_after_ti

  subroutine replace_vi_in_Tri_ti_with_vj( mesh, ti, vi, vj)
    ! Replace vi in Tri(ti,:) with tj

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh
    integer,         intent(in   ) :: ti, vi, vj

    ! Local variables:
    logical :: found_it
    integer :: n

    found_it = .false.
    do n = 1, 3
      if (mesh%Tri( ti,n) == vi) then
        found_it = .true.
        mesh%Tri( ti,n) = vj
        exit
      end if
    end do
#if (DO_ASSERTIONS)
    call assert( found_it, 'couldnt find vi in Tri(ti,:)')
#endif

  end subroutine replace_vi_in_Tri_ti_with_vj

  subroutine replace_tj_in_TriC_ti_with_tk( mesh, ti, tj, tk)
    ! Replace tj in TriC(ti,:) with tk

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh
    integer,         intent(in   ) :: ti, tj, tk

    ! Local variables:
    logical :: found_it
    integer :: n

    found_it = .false.
    do n = 1, 3
      if (mesh%TriC( ti,n) == tj) then
        found_it = .true.
        mesh%TriC( ti,n) = tk
        exit
      end if
    end do
#if (DO_ASSERTIONS)
    call assert( found_it, 'couldnt find tj in TriC(ti,:)')
#endif

  end subroutine replace_tj_in_TriC_ti_with_tk

end module basic_mesh_data_operations
