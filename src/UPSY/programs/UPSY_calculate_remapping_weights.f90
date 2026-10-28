program UPSY_calculate_remapping_weights

  ! A program for calculating remapping weights between two
  ! meshes/grids. Call syntax:
  !
  ! mpiexec -n [n_cores] UPSY_calculate_remapping_weights [filename_grid_src] [Arakawa_grid_src] [filename_grid_dst] [Arakawa_grid_dst] [filename_output]
  !
  ! n_cores          : number of CPU cores you want to run on
  ! filename_grid_src: path/to/source/grid_or_mesh.nc
  ! Arakawa_grid_src : Arakawa grid of source data (a, b, c, cx, cy)
  ! filename_grid_dst: path/to/destination/grid_or_mesh.nc
  ! Arakawa_grid_dst : Arakawa grid of destination data (a, b, c, cx, cy)
  ! filename_output  : path/to/output/file/containing/remapping_weights.nc

#include <petsc/finclude/petscsys.h>

  use basic_program_info, only: program_name
  use mpi_basic, only: par, sync
  use petsc, only: PetscErrorF, PetscFinalize
  use mpi_basic, only: initialise_parallelisation
  use parameters, only: initialise_constants
  use call_stack_and_comp_time_tracking, only: initialise_control_and_resource_tracker
  use crash_mod, only: crash
  use UPSY_main, only: UPSY
  use netcdf_io_main, only: inquire_xy_grid, inquire_lonlat_grid, inquire_mesh, &
    open_existing_netcdf_file_for_reading, setup_xy_grid_from_file, setup_lonlat_grid_from_file, &
    setup_mesh_from_file, close_netcdf_file
  use grid_types, only: type_grid, type_grid_lonlat
  use mesh_types, only: type_mesh
  use remapping_types, only: type_map

  implicit none

  integer                             :: ierr
  character(len=1024)                 :: filename_grid_src, Arakawa_grid_src
  character(len=1024)                 :: filename_grid_dst, Arakawa_grid_dst
  character(len=1024)                 :: filename_output
  logical                             :: file_exists
  logical                             :: has_xy_grid_src, has_lonlat_grid_src, has_mesh_src
  logical                             :: has_xy_grid_dst, has_lonlat_grid_dst, has_mesh_dst
  type(type_grid),        allocatable :: xy_grid_src, xy_grid_dst
  type(type_grid_lonlat), allocatable :: lonlat_grid_src, lonlat_grid_dst
  type(type_mesh),        allocatable :: mesh_src, mesh_dst
  type(type_map)                      :: map

  program_name = 'UPSY_calculate_remapping_weights'

  ! Initialise MPI parallelisation and PETSc
  call initialise_parallelisation

  ! Initialise constants (pi, NaN, ...)
  call initialise_constants

  ! Initialise the control and resource tracker
  call initialise_control_and_resource_tracker

  ! Get the input arguments
  call get_input_arguments( filename_grid_src, Arakawa_grid_src, filename_grid_dst, Arakawa_grid_dst, filename_output)

  ! Find out what kind of grids we're remapping between
  call read_grid_or_mesh( filename_grid_src, has_xy_grid_src, has_lonlat_grid_src, has_mesh_src, xy_grid_src, lonlat_grid_src, mesh_src)
  call read_grid_or_mesh( filename_grid_dst, has_xy_grid_dst, has_lonlat_grid_dst, has_mesh_dst, xy_grid_dst, lonlat_grid_dst, mesh_dst)

  ! Calculate the remapping operator
  call calc_remapping_operator( &
    has_xy_grid_src, has_lonlat_grid_src, has_mesh_src, xy_grid_src, lonlat_grid_src, mesh_src, &
    has_xy_grid_dst, has_lonlat_grid_dst, has_mesh_dst, xy_grid_dst, lonlat_grid_dst, mesh_dst, map)

  ! Write it to a NetCDF file
  call write_map_to_file( map, filename_output)

  ! Finalise PETSc and MPI parallelisation
  PetscCall( PetscFinalize( ierr))

contains

  subroutine get_input_arguments( filename_grid_src, Arakawa_grid_src, filename_grid_dst, Arakawa_grid_dst, filename_output)

    character(len=*), intent(out) :: filename_grid_src, Arakawa_grid_src, filename_grid_dst, Arakawa_grid_dst, filename_output

    call safety_check_input_arguments()

    call getarg( 1, filename_grid_src)
    call getarg( 2, Arakawa_grid_src )
    call getarg( 3, filename_grid_dst)
    call getarg( 4, Arakawa_grid_dst )
    call getarg( 5, filename_output  )

    if (par%primary) then
      write(0,*) UPSY%stru%colour_string( trim( program_name) // ': calculating remapping weights.', 'green')
      write(0,*) '  Source      mesh grid   : ', UPSY%stru%colour_string( trim( filename_grid_src), 'light blue')
      write(0,*) '  Source      Arakawa grid: ', UPSY%stru%colour_string( trim( Arakawa_grid_src ), 'light blue')
      write(0,*) '  Destination mesh grid   : ', UPSY%stru%colour_string( trim( filename_grid_dst), 'light blue')
      write(0,*) '  Destination Arakawa grid: ', UPSY%stru%colour_string( trim( Arakawa_grid_dst ), 'light blue')
      write(0,*) '  Output file             : ', UPSY%stru%colour_string( trim( filename_output  ), 'light blue')
    end if

    ! Check if the input files actually exist
    inquire( exist = file_exists, file = trim( filename_grid_src))
    if (.not. file_exists) call crash('file ' // trim( filename_grid_src) // ' not found')
    inquire( exist = file_exists, file = trim( filename_grid_dst))
    if (.not. file_exists) call crash('file ' // trim( filename_grid_dst) // ' not found')
    inquire( exist = file_exists, file = trim( filename_output))
    if (      file_exists) call crash('file ' // trim( filename_output) // ' already exists')

  end subroutine get_input_arguments

  subroutine safety_check_input_arguments()
    if (iargc() /= 5) then
      if (par%primary) then
        write(0,*) UPSY%stru%colour_string( trim( program_name) // 'Error: invalid call syntax for ' // trim( program_name) // '. Should be:', 'yellow')
        write(0,*) UPSY%stru%colour_string( '  mpiexec -n [n_cores] UPSY_calculate_remapping_weights [filename_grid_src] [Arakawa_grid_src] [filename_grid_dst] [Arakawa_grid_dst] [filename_output]', 'yellow')
        write(0,*) UPSY%stru%colour_string( '    n_cores          : number of CPU cores you want to run on', 'yellow')
        write(0,*) UPSY%stru%colour_string( '    filename_grid_src: path/to/source/grid_or_mesh.nc', 'yellow')
        write(0,*) UPSY%stru%colour_string( '    Arakawa_grid_src : Arakawa grid of source data (a, b, c, cx, cy)', 'yellow')
        write(0,*) UPSY%stru%colour_string( '    filename_grid_dst: path/to/destination/grid_or_mesh.nc', 'yellow')
        write(0,*) UPSY%stru%colour_string( '    Arakawa_grid_dst : Arakawa grid of destination data (a, b, c, cx, cy)', 'yellow')
        write(0,*) UPSY%stru%colour_string( '    filename_output  : path/to/output/file/containing/remapping_weights.nc', 'yellow')
      end if
      call sync
      error stop
    end if
  end subroutine safety_check_input_arguments

  subroutine read_grid_or_mesh( filename, has_xy_grid, has_lonlat_grid, has_mesh, xy_grid, lonlat_grid, mesh)

    ! In/output variables:
    character(len=*),                    intent(in   ) :: filename
    logical,                             intent(  out) :: has_xy_grid
    logical,                             intent(  out) :: has_lonlat_grid
    logical,                             intent(  out) :: has_mesh
    type(type_grid),        allocatable, intent(inout) :: xy_grid
    type(type_grid_lonlat), allocatable, intent(inout) :: lonlat_grid
    type(type_mesh),        allocatable, intent(inout) :: mesh

    ! Local variables:
    integer :: ncid

    ! Inquire what kind of grid the file is defined on
    call inquire_xy_grid(     filename, has_xy_grid    )
    call inquire_lonlat_grid( filename, has_lonlat_grid)
    call inquire_mesh(        filename, has_mesh       )

    ! Safety
    if (has_xy_grid     .and. has_lonlat_grid) call crash('file ' // trim( filename) // ' contains both an x/y-grid and a lon/lat-grid')
    if (has_xy_grid     .and. has_mesh       ) call crash('file ' // trim( filename) // ' contains both an x/y-grid and a mesh')
    if (has_lonlat_grid .and. has_mesh       ) call crash('file ' // trim( filename) // ' contains both a lon/lat-grid and a mesh')

    ! Read the grid from the file
    call open_existing_netcdf_file_for_reading( filename, ncid)
    if (has_xy_grid) then
      allocate( xy_grid)
      call setup_xy_grid_from_file( filename, ncid, xy_grid)
    elseif (has_lonlat_grid) then
      allocate( lonlat_grid)
      call setup_lonlat_grid_from_file( filename, ncid, lonlat_grid)
    elseif (has_mesh) then
      allocate( mesh)
      call setup_mesh_from_file( filename, ncid, mesh)
    else
      call crash('could not find a valid grid in file ' // trim( filename))
    end if
    call close_netcdf_file( ncid)

  end subroutine read_grid_or_mesh

  subroutine calc_remapping_operator( &
    has_xy_grid_src, has_lonlat_grid_src, has_mesh_src, xy_grid_src, lonlat_grid_src, mesh_src, &
    has_xy_grid_dst, has_lonlat_grid_dst, has_mesh_dst, xy_grid_dst, lonlat_grid_dst, mesh_dst, map)

    ! In/output variables
    logical,                             intent(in   ) :: has_xy_grid_src
    logical,                             intent(in   ) :: has_lonlat_grid_src
    logical,                             intent(in   ) :: has_mesh_src
    type(type_grid),        allocatable, intent(in   ) :: xy_grid_src
    type(type_grid_lonlat), allocatable, intent(in   ) :: lonlat_grid_src
    type(type_mesh),        allocatable, intent(in   ) :: mesh_src
    logical,                             intent(in   ) :: has_xy_grid_dst
    logical,                             intent(in   ) :: has_lonlat_grid_dst
    logical,                             intent(in   ) :: has_mesh_dst
    type(type_grid),        allocatable, intent(in   ) :: xy_grid_dst
    type(type_grid_lonlat), allocatable, intent(in   ) :: lonlat_grid_dst
    type(type_mesh),        allocatable, intent(in   ) :: mesh_dst
    type(type_map),                      intent(  out) :: map

    ! Local variables

    call crash('fixme!')

  end subroutine calc_remapping_operator

  subroutine write_map_to_file( map, filename_output)

    ! In/output variables
    type(type_map),   intent(in) :: map
    character(len=*), intent(in) :: filename_output

    ! Local variables

    call crash('fixme!')

  end subroutine write_map_to_file

end program UPSY_calculate_remapping_weights