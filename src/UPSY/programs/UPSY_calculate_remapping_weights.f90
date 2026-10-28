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

  use ct_basic, only: create_component_tests_output_folder
  use ct_create_test_meshes, only: create_all_test_meshes_and_grids_Antarctica

  implicit none

  integer             :: ierr
  character(len=1024) :: filename_grid_src, Arakawa_grid_src
  character(len=1024) :: filename_grid_dst, Arakawa_grid_dst
  character(len=1024) :: filename_output

  program_name = 'UPSY_calculate_remapping_weights'

  ! Initialise MPI parallelisation and PETSc
  call initialise_parallelisation

  ! Initialise constants (pi, NaN, ...)
  call initialise_constants

  ! Initialise the control and resource tracker
  call initialise_control_and_resource_tracker

  ! Get the input arguments
  if (iargc() == 5) then
    call getarg( 1, filename_grid_src)
    call getarg( 2, Arakawa_grid_src )
    call getarg( 3, filename_grid_dst)
    call getarg( 4, Arakawa_grid_dst )
    call getarg( 5, filename_output  )
  else
    if (par%primary) write(0,*) UPSY%stru%colour_string( trim( program_name) // 'Error: invalid call syntax for ' // trim( program_name) // '. Should be:', 'yellow')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '  mpiexec -n [n_cores] UPSY_calculate_remapping_weights [filename_grid_src] [Arakawa_grid_src] [filename_grid_dst] [Arakawa_grid_dst] [filename_output]', 'yellow')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '    n_cores          : number of CPU cores you want to run on', 'yellow')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '    filename_grid_src: path/to/source/grid_or_mesh.nc', 'yellow')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '    Arakawa_grid_src : Arakawa grid of source data (a, b, c, cx, cy)', 'yellow')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '    filename_grid_dst: path/to/destination/grid_or_mesh.nc', 'yellow')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '    Arakawa_grid_dst : Arakawa grid of destination data (a, b, c, cx, cy)', 'yellow')
    if (par%primary) write(0,*) UPSY%stru%colour_string( '    filename_output  : path/to/output/file/containing/remapping_weights.nc', 'yellow')
    call sync
    error stop
  end if

  ! Finalise PETSc and MPI parallelisation
  PetscCall( PetscFinalize( ierr))

end program UPSY_calculate_remapping_weights