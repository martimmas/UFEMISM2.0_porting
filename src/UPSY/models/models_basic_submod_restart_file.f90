submodule( models_basic) models_basic_submod_restart_file

use netcdf_io_main
use fields_main, only: type_third_dimension

contains

  subroutine write_to_restart_file( model, output_dir)

    ! In/output variables:
    class(atype_model), intent(in) :: model
    character(len=*),   intent(in) :: output_dir

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_restart_file'
    character(len=1024)            :: filename
    integer                        :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    filename = trim(  output_dir) // '/restart_file_' // trim( model%name()) // '.nc'

    call create_new_netcdf_file_for_writing( filename, ncid)

    ! Set up grid/mesh
    select type (grid => model%grid_val)
    class default
      call crash('invalid model%grid type - programming error')
    class is (type_grid)
      call setup_xy_grid_in_netcdf_file( filename, ncid, grid)
    class is (type_mesh)
      call setup_mesh_in_netcdf_file( filename, ncid, grid)
    end select

    call model%flds_reg%write_to_netcdf( filename, ncid)

    call close_netcdf_file( ncid)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_restart_file

end submodule models_basic_submod_restart_file
