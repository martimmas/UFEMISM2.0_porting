submodule (fields_basic) fields_basic_submod_read_from_netcdf

  use netcdf_io_main
  use mpi_distributed_memory_grid, only: distribute_gridded_data_from_primary

contains

  subroutine read_from_netcdf( field, filename, ncid)
    ! NOTE: assumes the NetCDF file already exists, is open, and
    !       already contains the grid/mesh and third dimension

    ! In/output variables:
    class(atype_field), intent(inout) :: field
    character(len=*),   intent(in   ) :: filename
    integer,            intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'read_from_netcdf'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f => field)
    class default
      call crash('invalid field type')
    class is (atype_field_2D)
    class is (atype_field_3D)
      if (.not. f%third_dimension_val%is_dim_and_var_in_netcdf( filename, ncid)) then
        call crash('third dimension of field "' // trim( field%name()) // &
          '" not found in file "' // trim( filename) // '"')
      end if
    end select

    select type (g => field%grid())
    class default
      call crash('invalid field%grid type')
    class is (type_grid)
      call read_from_netcdf_grid( field, g, filename, ncid)
    class is (type_mesh)
      call read_from_netcdf_mesh( field, g, filename, ncid)
    end select

    call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf

  ! ===== Grid-based fields
  ! =======================

  subroutine read_from_netcdf_grid( field, grid, filename, ncid)

    ! In/output variables:
    class(atype_field), intent(inout) :: field
    type(type_grid),    intent(in   ) :: grid
    character(len=*),   intent(in   ) :: filename
    integer,            intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'read_from_netcdf_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f => field)
    class default
      call crash('invalid field type')
    class is (type_field_logical_2D)
      call read_from_netcdf_grid_logical_2D( f, grid, filename, ncid)
    end select

    call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_grid

  subroutine read_from_netcdf_grid_logical_2D( field, grid, filename, ncid)

    ! In/output variables:
    type(type_field_logical_2D), intent(inout) :: field
    type(type_grid),             intent(in   ) :: grid
    character(len=*),            intent(in   ) :: filename
    integer,                     intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter       :: routine_name = 'read_from_netcdf_grid_logical_2D'
    logical, dimension(:,:), allocatable :: d_grid_tot
    logical, dimension(:  ), allocatable :: d_grid_vec_partial

    ! Add routine to call stack
    call init_routine( routine_name)

    ! call distribute_gridded_data_from_primary( grid, d_grid_tot, d_grid_vec_partial)

    call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_grid_logical_2D

  ! ===== Mesh-based fields
  ! =======================

  subroutine read_from_netcdf_mesh( field, mesh, filename, ncid)

    ! In/output variables:
    class(atype_field), intent(inout) :: field
    type(type_mesh),    intent(in   ) :: mesh
    character(len=*),   intent(in   ) :: filename
    integer,            intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'read_from_netcdf_mesh'

    ! Add routine to call stack
    call init_routine( routine_name)

    call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf_mesh

end submodule fields_basic_submod_read_from_netcdf