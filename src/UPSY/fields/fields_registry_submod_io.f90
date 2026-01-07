submodule( fields_registry) fields_registry_submod_io

  use fields_dimensions, only: type_third_dimension
  use fields_basic, only: atype_field_2D, atype_field_3D

contains

  subroutine write_to_netcdf( flds_reg, filename, ncid)
    ! NOTE: assumes the NetCDF file already exists, is open, and
    !       already contains the grid/mesh

    ! In/output variables:
    class(type_fields_registry), intent(in) :: flds_reg
    character(len=*),            intent(in) :: filename
    integer,                     intent(in) :: ncid

    ! Local variables:
    character(len=1024), parameter                        :: routine_name = 'write_to_netcdf'
    type(type_third_dimension), dimension(:), allocatable :: dims
    integer                                               :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    do i = 1, flds_reg%n
      call flds_reg%items(i)%p%write_to_netcdf( filename, ncid)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine write_to_netcdf

end submodule fields_registry_submod_io