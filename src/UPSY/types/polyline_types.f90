module polyline_types

  use precisions, only: dp

  implicit none

  private

  public :: type_polyline

  type type_polyline
    logical                               :: is_closed           !<     Whether or not the polyline is a closed loop
    integer                               :: n                   !<     The number of vertices spanning this polyline
    real(dp), dimension(:,:), allocatable :: p                   !< [m] x/y-coordinates of the vertices spanning this polyline
  end type type_polyline

contains

end module polyline_types
