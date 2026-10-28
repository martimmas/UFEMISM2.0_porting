module grid_contour

  ! Functions for working with simple square x/y-grids

  use precisions, only: dp
  use xy_grid_basic, only: type_xy_grid
  use call_stack_and_comp_time_tracking, only: crash, init_routine, finalise_routine
  use reallocate_mod, only: reallocate
  use interpolation, only: linint_points

  implicit none

  private

  public :: calc_grid_contour_as_line, calc_grid_mask_as_polygons

contains

  subroutine calc_grid_contour_as_line( grid, d, f, line, mask)
    !< Calculate a contour line at level f for data d on a square grid.
    !< Generate the contour line in UFEMISM line-segment format (i.e. unordered
    !< individual line segments).

    ! In/output variables
    type(type_xy_grid),                    intent(in   ) :: grid
    real(dp), dimension(:,:),              intent(in   ) :: d
    real(dp),                              intent(in   ) :: f
    real(dp), dimension(:,:), allocatable, intent(  out) :: line
    logical,  dimension(:,:), optional,    intent(in   ) :: mask

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'calc_grid_contour_as_line'
    real(dp), parameter                   :: tol = 1E-5_dp
    real(dp), dimension(:,:), allocatable :: d_scaled
    integer                               :: i,j
    logical,  dimension(:,:), allocatable :: mask_loc
    integer                               :: n_max, n
    real(dp)                              :: d_sw, d_nw, d_se, d_ne
    real(dp)                              :: xw, xe, xs, xn, yw, ye, ys, yn
    logical                               :: do_cross_w, do_cross_e, do_cross_s, do_cross_n

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( d,1) /= grid%nx .or. size( d,2) /= grid%ny) call crash('d is not nx-by-ny!')

    ! Trivial case: if all values of d are greater or smaller than f,
    ! the contour line is empty
    if (minval( d) >= f .or. maxval( d) <= f) then
      allocate( line( 0,0))
      call finalise_routine( routine_name)
      return
    end if

    ! Shift d so the contour lies at d_scaled = 0
    allocate( d_scaled( grid%nx, grid%ny))
    d_scaled = d - f

    ! Set the mask to optionally skip certain grid cells
    allocate( mask_loc( grid%nx, grid%ny))
    if (present( mask)) then
      mask_loc = mask
    else
      mask_loc =  .true.
    end if

    n_max = 1000
    allocate( line( n_max, 4))

    n = 0
    do i = 1, grid%nx-1
    do j = 1, grid%ny-1

      ! Skip this grid cell if told so
      if ((.not. mask_loc( i  ,j  )) .and. &
          (.not. mask_loc( i  ,j+1)) .and. &
          (.not. mask_loc( i+1,j  )) .and. &
          (.not. mask_loc( i+1,j+1))) cycle

      ! Extend allocated memory IF needed
      if (n > n_max - 10) then
        n_max = n + 1000
        call reallocate( line, n_max, 4)
      end if

      ! The four corners of the b-grid cell
      d_sw = d_scaled( i  ,j  )
      d_nw = d_scaled( i  ,j+1)
      d_se = d_scaled( i+1,j  )
      d_ne = d_scaled( i+1,j+1)

      xw = grid%x( i  )
      xe = grid%x( i+1)
      ys = grid%y( j  )
      yn = grid%y( j+1)

      ! If all four corners are above/below the level, no line here
      if ((d_sw >= 0._dp .and. d_nw >= 0._dp .and. d_se >= 0._dp .and. d_ne >= 0._dp) .or. &
          (d_sw <= 0._dp .and. d_nw <= 0._dp .and. d_se <= 0._dp .and. d_ne <= 0._dp)) then
        cycle
      end if

      ! Add tolerances to keep line lengths finite
      if (d_sw >= 0._dp) then
        d_sw = max( d_sw, tol)
      else
        d_sw = min( d_sw, tol)
      end if
      if (d_nw >= 0._dp) then
        d_nw = max( d_nw, tol)
      else
        d_nw = min( d_nw, tol)
      end if
      if (d_se >= 0._dp) then
        d_se = max( d_se, tol)
      else
        d_se = min( d_se, tol)
      end if
      if (d_ne >= 0._dp) then
        d_ne = max( d_ne, tol)
      else
        d_ne = min( d_ne, tol)
      end if

      ! Find boundary crossings

      do_cross_w = .false.
      do_cross_e = .false.
      do_cross_s = .false.
      do_cross_n = .false.

      if (d_sw * d_nw < 0._dp) then
        ! The contour crosses the western boundary of this b-grid cell
        do_cross_w = .true.
        yw = linint_points( ys, yn, d_sw, d_nw, 0._dp)
      end if

      if (d_se * d_ne < 0._dp) then
        ! The contour crosses the eastern boundary of this b-grid cell
        do_cross_e = .true.
        ye = linint_points( ys, yn, d_se, d_ne, 0._dp)
      end if

      if (d_nw * d_ne < 0._dp) then
        ! The contour crosses the northern boundary of this b-grid cell
        do_cross_n = .true.
        xn = linint_points( xw, xe, d_nw, d_ne, 0._dp)
      end if

      if (d_sw * d_se < 0._dp) then
        ! The contour crosses the southern boundary of this b-grid cell
        do_cross_s = .true.
        xs = linint_points( xw, xe, d_sw, d_se, 0._dp)
      end if

      ! Add line segments
      n = n + 1
      if     (do_cross_w) then
        if     (do_cross_e) then
          ! From west to east
          line( n,:) = [xw,yw,xe,ye]
        elseif (do_cross_s) then
          ! From west to south
          line( n,:) = [xw,yw,xs,ys]
        elseif (do_cross_n) then
          ! From west to north
          line( n,:) = [xw,yw,xn,yn]
        else
          call crash('found only a crossing at the western boundary!')
        end if
      elseif (do_cross_e) then
        IF     (do_cross_s) then
          ! From east to south
          line( n,:) = [xe,ye,xs,ys]
        elseif (do_cross_n) then
          ! From east to north
          line( n,:) = [xe,ye,xn,yn]
        else
          call crash('found only a crossing at the eastern boundary!')
        end if
      elseif (do_cross_s) then
        IF     (do_cross_n) then
          ! From south to north
          line( n,:) = [xs,ys,xn,yn]
        else
          call crash('found only a crossing at the southern boundary!')
        end if
      elseif (do_cross_n) then
          call crash('found only a crossing at the northern boundary!')
      else
        call crash('whaa!')
      end if

    end do
    end do

    ! Crop memory
    call reallocate( line, n, 4)

    ! Clean up after yourself
    deallocate( d_scaled)
    deallocate( mask_loc)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grid_contour_as_line

  subroutine calc_grid_mask_as_polygons( grid, mask, poly_mult)
    !< Calculate a set of polygon enveloping all TRUE-valued mask cells

    ! In/output variables
    type(type_xy_grid),                    intent(in   ) :: grid
    logical,  dimension(:,:),              intent(in   ) :: mask
    real(dp), dimension(:,:), allocatable, intent(  out) :: poly_mult

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'calc_grid_mask_as_polygons'
    logical,  dimension(:,:), allocatable :: mask_loc
    integer                               :: i,j
    real(dp), dimension(:,:), allocatable :: poly
    integer                               :: n_poly, n_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( mask,1) /= grid%nx .or. size( mask,2) /= grid%ny) call crash('incorrect data dimensions!')

    ! Trivial case for no TRUE values at all
    if (.not. any( mask)) then
      allocate( poly_mult( 0,0))
      call finalise_routine( routine_name)
      return
    end if

    ! Make a local copy of the logical mask
    allocate( mask_loc( grid%nx, grid%ny))
    mask_loc = mask

    ! Initialise poly_mult and poly
    allocate( poly_mult( grid%nx*grid%ny,2))
    allocate( poly(      grid%nx*grid%ny,2))
    n_tot = 0

    ! Calculate polygons for all TRUE regions of the mask
    do i = 1, grid%nx
    do j = 1, grid%ny

      if (mask_loc( i,j)) then
        ! Found a seed for a TRUE region

        ! Calculate a polygon enveloping this TRUE region, and
        ! remove the region from the mask
        call calc_grid_mask_as_polygon( grid, mask_loc, i, j, poly, n_poly)

        ! Add this polygon to poly_mult
        poly_mult( n_tot+1,1) = real( n_poly,dp)
        poly_mult( n_tot+1,2) = 0._dp
        poly_mult( n_tot+2:n_tot+1+n_poly,:) = poly( 1:n_poly,:)
        n_tot = n_tot + 1 + n_poly

      end if ! if (mask_loc( i,j)) then

    end do
    end do

    ! Crop memory
    call reallocate( poly_mult, n_tot, 2)

    ! Clean up after yourself
    deallocate( mask_loc)
    deallocate( poly)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grid_mask_as_polygons

  subroutine calc_grid_mask_as_polygon( grid, mask, i0, j0, poly, n_poly)
    !< Calculate a polygon enveloping the set of TRUE-valued mask cells around i0,j0,
    !< and remove that set of grid cells from the mask

    ! In/output variables
    type(type_xy_grid),       intent(in   ) :: grid
    logical,  dimension(:,:), intent(inout) :: mask
    integer,                  intent(in   ) :: i0,j0
    real(dp), dimension(:,:), intent(  out) :: poly
    integer,                  intent(  out) :: n_poly

    ! Local variables:
    character(len=*), parameter           :: routine_name = 'calc_grid_mask_as_polygon'
    logical,  dimension(:,:), allocatable :: mask_ext
    integer,  dimension(:,:), allocatable :: map
    integer,  dimension(:,:), allocatable :: stack
    integer                               :: stackN
    integer                               :: i,j,ii,jj,i_sw,j_sw
    character(len=256)                    :: dir, dir_prev
    integer                               :: it

    ! Add routine to path
    call init_routine( routine_name)

    ! Safety
    if (size( mask,1) /= grid%nx .or. size( mask,2) /= grid%ny) call crash('incorrect data dimensions!')
    if (.not. mask( i0,j0)) call crash('seed at i0,j0 is not TRUE!')

    ! Pad a row of FALSEs around the mask so the set of TRUEs around i0,j0 is really an island
    allocate( mask_ext( 0:grid%nx+1, 0:grid%ny+1), source = .false.)
    mask_ext( 1:grid%nx, 1:grid%ny) = mask

    ! Use a flood-fill algorithm to find the map of same-valued grid cells around mask cell i0,j0
    allocate( map(   0:grid%nx+1, 0:grid%ny+1) , source = 0)
    allocate( stack( (grid%nx+2)*(grid%ny+2),2))

    map( i0,j0) = 1
    stackN = 1
    stack( 1,:) = [i0,j0]

    i_sw = 0
    j_sw = 0

    do while (stackN > 0)

      ! Take the last element from the stack
      i = stack( stackN,1)
      j = stack( stackN,2)
      stackN = stackN - 1

      ! Mark it as mapped
      map( i,j) = 2

      ! Remove it from the input mask
      mask( i,j) = .false.

      ! Add all non-mapped, non-stacked, TRUE-valued neighbours to the stack
      DO ii = max( 0, i-1), min( grid%nx+1, i+1)
      DO jj = max( 0, j-1), min( grid%ny+1, j+1)
        if (ii /= i .and. jj /= j) cycle ! Don't include diagonal neighbours
        if (map( ii,jj) == 0 .and. mask_ext( ii,jj)) then
          ! Add this neighbour to the stack
          stackN = stackN + 1
          stack( stackN,:) = [ii,jj]
          ! Mark this neighbour on the map as stacked
          map( ii,jj) = 1
        end if
      end do
      end do

      ! If it is a southwest corner, save it as a starting point for the outline tracer
      if (.not. mask_ext( i-1,j) .and. .not. mask_ext( i,j-1)) then
        i_sw = i
        j_sw = j
      end if

    end do ! do while (stackN > 0)
    ! Safety
    if (i_sw == 0 .or. j_sw == 0) call crash('couldnt find starting SW corner!')

    ! Start at the southwest corner we found earlier
    ii = i_sw
    jj = j_sw
    n_poly = 2
    poly( 1,:) = [grid%x( ii) - grid%dx/2._dp, grid%y( jj) - grid%dx/2._dp]
    poly( 2,:) = [grid%x( ii) + grid%dx/2._dp, grid%y( jj) - grid%dx/2._dp]
    dir = 'east'

    ! Trace the outline of the mapped grid cells.
    it = 0
    do while (.true.)

      ! Safety
      it = it + 1
      if (it > grid%nx*grid%ny) call crash('outline tracer got stuck!')

      ! cycle direction
      dir_prev = dir

      ! Check which way we go next
      if     (dir_prev == 'east') then
        ! We were going east, so [ii,jj] is to the north of us, and we can
        ! continue north, east, or south

        if     (map( ii+1,jj) == 0) then
          ! Go north
          dir = 'north'
        elseif (map( ii+1,jj-1) == 0) then
          ! Go east
          dir = 'east'
          ii = ii+1
        elseif (map( ii,jj-1) == 0) then
          ! Go south
          dir = 'south'
          ii = ii+1
          jj = jj-1
        else
          call crash('outline tracer got stuck while going east!')
        end if

      elseif (dir_prev == 'north') then
        ! We were going north, so [ii,jj] is to the west of us, and we can
        ! continue west, north, or east

        if     (map( ii,jj+1) == 0) then
          ! Go west
          dir = 'west'
        elseif (map( ii+1,jj+1) == 0) then
          ! Go north
          dir = 'north'
          jj = jj+1
        elseif (map( ii+1,jj) == 0) then
          ! Go east
          dir = 'east'
          ii = ii+1
          jj = jj+1
        else
          call crash('outline tracer got stuck while going north!')
        end if

      elseif (dir_prev == 'west') then
        ! We were going west, so [ii,jj] is to the south of us, and we can
        ! continue south, west, or north

        if     (map( ii-1,jj) == 0) then
          ! Go south
          dir = 'south'
        elseif (map( ii-1,jj+1) == 0) then
          ! Go west
          dir = 'west'
          ii = ii-1
        elseif (map( ii,jj+1) == 0) then
          ! Go north
          dir = 'north'
          ii = ii-1
          jj = jj+1
        else
          call crash('outline tracer got stuck while going north!')
        end if

      elseif (dir_prev == 'south') then
        ! We were going south, so [ii,jj] is to the east of us, and we can
        ! continue east, south, or west

        if     (map( ii,jj-1) == 0) then
          ! Go east
          dir = 'east'
        elseif (map( ii-1,jj-1) == 0) then
          ! Go south
          dir = 'south'
          jj = jj-1
        elseif (map( ii-1,jj) == 0) then
          ! Go west
          dir = 'west'
          ii = ii-1
          jj = jj-1
        else
          call crash('outline tracer got stuck while going north!')
        end if

      else
        call crash('unknown dir_prev "' // trim( dir_prev) // '"!')
      end if

      ! Add new vertex to the polygon
      n_poly = n_poly+1
      if     (dir == 'east') then
        poly( n_poly,:) = poly( n_poly-1,:) + [grid%dx, 0._dp]
      elseif (dir == 'north') then
        poly( n_poly,:) = poly( n_poly-1,:) + [0._dp, grid%dx]
      elseif (dir == 'west') then
        poly( n_poly,:) = poly( n_poly-1,:) + [-grid%dx, 0._dp]
      elseif (dir == 'south') then
        poly( n_poly,:) = poly( n_poly-1,:) + [0._dp, -grid%dx]
      else
        call crash('unknown dir "' // trim( dir) // '"!')
      end if

      ! If we've reached the starting point again, we're done
      if (norm2( poly( n_poly,:) - poly( 1,:)) < grid%dx / 10._dp) then
        ! Don't double count
        n_poly = n_poly-1
        exit
      end if

    end do ! do while (.true.)

    ! Clean up after yourself
    deallocate( mask_ext)
    deallocate( map)
    deallocate( stack)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_grid_mask_as_polygon

end module grid_contour
