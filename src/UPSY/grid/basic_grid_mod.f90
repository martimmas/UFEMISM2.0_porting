module basic_grid_mod

  use precisions, only: dp

  implicit none

  private

  public :: atype_grid

  type, abstract :: atype_grid

    ! Metadata
    character(len=1024) :: name

  end type atype_grid

end module basic_grid_mod