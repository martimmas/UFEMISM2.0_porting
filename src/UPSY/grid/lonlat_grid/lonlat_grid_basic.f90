module lonlat_grid_basic

  ! Functions for working with simple lon/lat-grids
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: warning, crash, happy, init_routine, finalise_routine, colour_string
  use parameters
  use reallocate_mod, only: reallocate
  use interpolation, only: linint_points
  use projections, only: inverse_oblique_sg_projection
  use mpi_distributed_memory, only: partition_list, distribute_from_primary, gather_to_primary

  implicit none

  private

  public :: type_lonlat_grid, setup_lonlat_grid, deallocate_lonlat_grid, &
    calc_lonlat_field_to_vector_form_translation_tables, &
    distribute_lonlat_gridded_data_from_primary_dp_2D, &
    distribute_lonlat_gridded_data_from_primary_dp_3D, &
    gather_lonlat_gridded_data_to_primary_dp_2D, &
    gather_lonlat_gridded_data_to_primary_dp_3D

  type type_lonlat_grid
    ! A regular lon/lat-grid

    ! Basic properties
    character(len=256)                      :: name         !           A nice name tag
    integer                                 :: nlon         !           Number of grid cells in the longitude direction
    integer                                 :: nlat         !           Number of grid cells in the latitude direction
    integer                                 :: n            !           Total number of grid cells (= nx * ny)
    real(dp)                                :: dlon         ! [degrees] Resolution in the longitude direction
    real(dp)                                :: dlat         ! [degrees] Resolution in the latitude direction
    real(dp), dimension(:    ), allocatable :: lon          ! [degrees east ] Longitude of each grid point
    real(dp), dimension(:    ), allocatable :: lat          ! [degrees north] Latitude  of each grid point
    real(dp)                                :: lonmin       ! [degrees east ] Lon/lat range covered by the grid
    real(dp)                                :: lonmax       ! [degrees east ]
    real(dp)                                :: latmin       ! [degrees north]
    real(dp)                                :: latmax       ! [degrees north]

    ! Remapping data
    real(dp)                                :: tol_dist     ! [m]       Horizontal distance tolerance; points closer together than this are assumed to be identical (typically set to a billionth of linear domain size)
    integer,  dimension(:,:  ), allocatable :: ij2n, n2ij   !           Conversion tables for grid-form vs. vector-form data

    ! Parallelisation
    integer                                 :: n1,n2,n_loc  ! Matrix rows owned by each process

  end type type_lonlat_grid

contains

! == Basic square grid functionality

  subroutine setup_lonlat_grid( name, nlon, nlat, grid)
    ! Set up a simple lon/lat-grid that covers the entire globe

    ! In/output variables:
    character(len=*),       intent(in   ) :: name
    integer,                intent(in   ) :: nlon, nlat
    type(type_lonlat_grid), intent(  out) :: grid

    ! Local variables:
    character(len=*), parameter :: routine_name = 'setup_lonlat_grid'
    integer                     :: i,j

    ! Add routine to path
    call init_routine( routine_name)

    ! Name
    grid%name = name

    ! size
    grid%nlon = nlon
    grid%nlat = nlat

    ! allocate memory
    allocate( grid%lon( grid%nlon))
    allocate( grid%lat( grid%nlat))

    ! Fill in lon/lat values
    grid%dlon = 360._dp / real( grid%nlon, dp)
    do i = 1, grid%nlon
      grid%lon( i) = real( i, dp) * grid%dlon
    end do

    grid%dlat = 180._dp / real( grid%nlat-1, dp)
    do j = 1, grid%nlat
      grid%lat( j) = -180._dp + real( j-1, dp) * grid%dlat
    end do

    ! Metadata
    grid%lonmin = grid%lon( 1)
    grid%lonmax = grid%lon( grid%nlon)
    grid%latmin = grid%lat( 1)
    grid%latmax = grid%lat( grid%nlat)

    ! Secondary data
    call calc_lonlat_field_to_vector_form_translation_tables( grid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_lonlat_grid

  subroutine deallocate_lonlat_grid( grid)
    ! deallocate memory for a grid object

    ! In/output variables:
    type(type_lonlat_grid), intent(inout) :: grid

    ! Local variables:
    character(len=*), parameter :: routine_name = 'deallocate_lonlat_grid'

    ! Add routine to path
    call init_routine( routine_name)

    if (allocated( grid%lon )) deallocate( grid%lon )
    if (allocated( grid%lat )) deallocate( grid%lat )
    if (allocated( grid%ij2n)) deallocate( grid%ij2n)
    if (allocated( grid%n2ij)) deallocate( grid%n2ij)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine deallocate_lonlat_grid

  subroutine check_if_lonlat_grids_are_identical( grid1, grid2, isso)
    ! Check if two grids are identical

    ! In/output variables:
    type(type_lonlat_grid), intent(inout) :: grid1, grid2
    logical,                intent(  out) :: isso

    ! Local variables:
    character(len=*), parameter :: routine_name = 'check_if_lonlat_grids_are_identical'
    real(dp), parameter         :: tol = 1E-9_dp
    integer                     :: i,j

    ! Add routine to path
    call init_routine( routine_name)

    isso = .true.

    ! size
    if (grid1%nlon /= grid2%nlon .OR. grid1%nlat /= grid2%nlat) then
      isso = .false.
      return
    end if

    ! Coordinates
    do i = 1, grid1%nlon
      if ( 1._dp - max( 0.001_dp, abs( grid1%lon( i))) / max( 0.001_dp, abs( grid2%lon( i))) > tol) then
        isso = .false.
        return
      end if
    end do
    do j = 1, grid1%nlat
      if ( 1._dp - max( 0.001_dp, abs( grid1%lat( j))) / max( 0.001_dp, abs( grid2%lat( j))) > tol) then
        isso = .false.
        return
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine check_if_lonlat_grids_are_identical

  subroutine calc_lonlat_field_to_vector_form_translation_tables( grid)
    ! Calculate grid-cell-to-matrix-row translation tables

    ! In/output variables
    type(type_lonlat_grid), intent(inout) :: grid

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_lonlat_field_to_vector_form_translation_tables'
    integer                     :: i,j,n

    ! Add routine to path
    call init_routine( routine_name)

    ! Total number of grid cells
    grid%n = grid%nlon * grid%nlat

    ! allocate memory
    if (allocated( grid%ij2n)) deallocate( grid%ij2n)
    if (allocated( grid%n2ij)) deallocate( grid%n2ij)
    allocate( grid%ij2n( grid%nlon, grid%nlat), source = 0)
    allocate( grid%n2ij( grid%n   , 2        ), source = 0)

    ! Fill in tables
    n = 0
    do i = 1, grid%nlon
    do j = 1, grid%nlat
      n = n + 1
      grid%ij2n( i,j) = n
      grid%n2ij( n,:) = [i,j]
    end do
    end do

    ! Parallelisation domains
    call partition_list( grid%n, par%i, par%n, grid%n1, grid%n2)
    grid%n_loc = grid%n2 + 1 - grid%n1

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_lonlat_field_to_vector_form_translation_tables

! == subroutines for manipulating gridded data in distributed memory

  subroutine distribute_lonlat_gridded_data_from_primary_dp_2D( grid, d_grid, d_grid_vec_partial)
    ! Distribute a 2-D gridded data field from the primary.
    ! Input from primary: total data field in field form
    ! Output to all: partial data in vector form

    ! In/output variables:
    type(type_lonlat_grid),             intent(inout) :: grid
    real(dp), dimension(:,:), optional, intent(inout) :: d_grid
    real(dp), dimension(: ),            intent(  out) :: d_grid_vec_partial

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'distribute_lonlat_gridded_data_from_primary_dp_2D'
    integer                             :: n,i,j
    real(dp), dimension(:), allocatable :: d_grid_vec_total

    ! Add routine to path
    call init_routine( routine_name)

    ! Convert gridded data to vector form
    if (par%primary) then
      if (.not. present(d_grid)) call crash('d_grid should be present on the primary process')

      ! allocate memory
      allocate( d_grid_vec_total( grid%n), source = 0._dp)

      ! Convert to vector form
      do n = 1, grid%n
        i = grid%n2ij( n,1)
        j = grid%n2ij( n,2)
        d_grid_vec_total( n) = d_grid( i,j)
      end do

    end if

    ! Distribute vector-form data to the processes
    call distribute_from_primary( d_grid_vec_partial, d_grid_vec_total)

    ! Clean up after yourself
    if (par%primary) deallocate( d_grid_vec_total)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine distribute_lonlat_gridded_data_from_primary_dp_2D

  subroutine distribute_lonlat_gridded_data_from_primary_dp_3D( grid, d_grid, d_grid_vec_partial)
    ! Distribute a 3-D gridded data field from the primary.
    ! Input from primary: total data field in field form
    ! Output to all: partial data in vector form

    ! In/poutput variables:
    type(type_lonlat_grid),     intent(inout) :: grid
    real(dp), dimension(:,:,:), intent(inout) :: d_grid
    real(dp), dimension(:,:),   intent(  out) :: d_grid_vec_partial

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'distribute_lonlat_gridded_data_from_primary_dp_3D'
    integer                               :: k
    real(dp), dimension(:,:), allocatable :: d_grid_2D
    real(dp), dimension(:  ), allocatable :: d_grid_vec_partial_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (par%primary .and. size( d_grid,3) /= size( d_grid_vec_partial,2)) call crash('vector sizes dont match!')

    ! allocate memory
    if (par%primary) allocate( d_grid_2D( size( d_grid,1), size( d_grid,2)), source = 0._dp)
    allocate( d_grid_vec_partial_2D( size( d_grid_vec_partial,1)), source = 0._dp)

    ! Treat each layer as a separate 2-D field
    do k = 1, size( d_grid_vec_partial,2)
      if (par%primary) d_grid_2D = d_grid( :,:,k)
      call distribute_lonlat_gridded_data_from_primary_dp_2D( grid, d_grid_2D, d_grid_vec_partial_2D)
      d_grid_vec_partial( :,k) = d_grid_vec_partial_2D
    end do

    ! Clean up after yourself
    if (par%primary) deallocate( d_grid_2D)
    deallocate( d_grid_vec_partial_2D)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine distribute_lonlat_gridded_data_from_primary_dp_3D

  subroutine gather_lonlat_gridded_data_to_primary_dp_2D( grid, d_grid_vec_partial, d_grid)
    ! Gather a 2-D gridded data field to the primary.
    ! Input from all: partial data in vector form
    ! Output to primary: total data field in field form

    ! In/output variables:
    type(type_lonlat_grid),   intent(inout) :: grid
    real(dp), dimension(:),   intent(inout) :: d_grid_vec_partial
    real(dp), dimension(:,:), intent(  out) :: d_grid

    ! Local variables:
    character(len=*), parameter         :: routine_name = 'gather_lonlat_gridded_data_to_primary_dp_2D'
    integer                             :: n,i,j
    real(dp), dimension(:), allocatable :: d_grid_vec_total

    ! Add routine to path
    call init_routine( routine_name)

    ! allocate memory
    if (par%primary) allocate( d_grid_vec_total( grid%n), source = 0._dp)

    ! Gather data
    call gather_to_primary( d_grid_vec_partial, d_grid_vec_total)

    ! Convert to grid form
    if (par%primary) then
      do n = 1, grid%n
        i = grid%n2ij( n,1)
        j = grid%n2ij( n,2)
        d_grid( i,j) = d_grid_vec_total( n)
      end do
    end if

    ! Clean up after yourself
    if (par%primary) deallocate( d_grid_vec_total)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine gather_lonlat_gridded_data_to_primary_dp_2D

  subroutine gather_lonlat_gridded_data_to_primary_dp_3D( grid, d_grid_vec_partial, d_grid)
    ! Gather a 3-D gridded data field to the primary.
    ! Input from all: partial data in vector form
    ! Output to primary: total data field in field form

    ! In/output variables:
    type(type_lonlat_grid),     intent(inout) :: grid
    real(dp), dimension(:,:),   intent(inout) :: d_grid_vec_partial
    real(dp), dimension(:,:,:), intent(  out) :: d_grid

    ! Local variables:
    character(*), parameter               :: routine_name = 'gather_lonlat_gridded_data_to_primary_dp_3D'
    integer                               :: k
    real(dp), dimension(:,:), allocatable :: d_grid_2D
    real(dp), dimension(:  ), allocatable :: d_grid_vec_partial_2D

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (par%primary .and. size( d_grid,3) /= size( d_grid_vec_partial,2)) call crash('vector sizes dont match!')

    ! allocate memory
    if (par%primary) allocate( d_grid_2D( grid%nlon, grid%nlat), source = 0._dp)
    allocate( d_grid_vec_partial_2D( grid%n_loc), source = 0._dp)

    ! Treat each layer as a separate 2-D field
    do k = 1, size( d_grid_vec_partial,2)
      d_grid_vec_partial_2D = d_grid_vec_partial( :,k)
      call gather_lonlat_gridded_data_to_primary_dp_2D( grid, d_grid_vec_partial_2D, d_grid_2D)
      if (par%primary) d_grid( :,:,k) = d_grid_2D
    end do

    ! Clean up after yourself
    if (par%primary) deallocate( d_grid_2D)
    deallocate( d_grid_vec_partial_2D)

    ! Add routine to path
    call finalise_routine( routine_name)

  end subroutine gather_lonlat_gridded_data_to_primary_dp_3D

end module lonlat_grid_basic
