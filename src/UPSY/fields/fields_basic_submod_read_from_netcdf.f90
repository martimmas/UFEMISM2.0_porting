submodule (fields_basic) fields_basic_submod_read_from_netcdf

  use netcdf_io_main

contains

  subroutine read_from_netcdf( field, filename, ncid)
    ! NOTE: assumes the NetCDF file already exists, is open, and
    !       already contains the grid/mesh and third dimension

    ! In/output variables:
    class(atype_field), intent(inout) :: field
    character(len=*),   intent(in   ) :: filename
    integer,            intent(in   ) :: ncid

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'read_from_netcdf'

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

    call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine read_from_netcdf

end submodule fields_basic_submod_read_from_netcdf