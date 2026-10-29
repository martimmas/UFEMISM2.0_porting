module lat_grid_basic

  use precisions, only: dp
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use basic_grid_mod, only: atype_grid

  implicit none

  private

  public :: type_lat_grid, deallocate_lat_grid

  type, extends(atype_grid) :: type_lat_grid
    ! A lat-only "grid"

    ! Basic properties
    integer                                 :: nlat         !                 Number of grid cells in the latitude direction
    integer                                 :: n            !                 Total number of grid cells (= nx * ny)
    real(dp)                                :: dlat         ! [degrees]       Resolution in the latitude direction
    real(dp), dimension(:    ), allocatable :: lat          ! [degrees north] Latitude  of each grid point
    real(dp)                                :: latmin       ! [degrees north]
    real(dp)                                :: latmax       ! [degrees north]

    ! Remapping data
    real(dp)                                :: tol_dist     ! [m]       Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    integer,  dimension(:,:  ), allocatable :: ij2n, n2ij   !           Conversion tables for grid-form vs. vector-form data

    ! Parallelisation
    integer                                 :: n1,n2,n_loc  ! Matrix rows owned by each process

  end type type_lat_grid

contains

  subroutine deallocate_lat_grid( grd)
    ! deallocate memory for a lat-only grid object

    ! In/output variables:
    type(type_lat_grid), intent(inout) :: grd

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_lat_grid'

    ! Add routine to path
    call init_routine( routine_name)

    if (allocated( grd%lat )) deallocate( grd%lat )
    if (allocated( grd%ij2n)) deallocate( grd%ij2n)
    if (allocated( grd%n2ij)) deallocate( grd%n2ij)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_lat_grid

end module lat_grid_basic
