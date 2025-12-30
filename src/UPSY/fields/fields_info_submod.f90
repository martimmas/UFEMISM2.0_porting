submodule(fields_basic) fields_info_submod

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string, crash
  use Arakawa_grid_mod, only: type_Arakawa_grid

contains

  subroutine print_field_grid_info_logical_2D(self)
    class(type_field_grid_logical_2D), intent(in) :: self
    call print_field_info_general( self%name, '2D', 'logical', 'grid', self%parent%name, &
      self%Arakawa_grid, self%long_name, self%units)
  end subroutine print_field_grid_info_logical_2D

  subroutine print_field_grid_info_logical_3D(self)
    class(type_field_grid_logical_3D), intent(in) :: self
    call print_field_info_general( self%name, '3D', 'logical', 'grid', self%parent%name, &
      self%Arakawa_grid, self%long_name, self%units)
  end subroutine print_field_grid_info_logical_3D

  subroutine print_field_grid_info_int_2D(self)
    class(type_field_grid_int_2D), intent(in) :: self
    call print_field_info_general( self%name, '2D', 'integer', 'grid', self%parent%name, &
      self%Arakawa_grid, self%long_name, self%units)
  end subroutine print_field_grid_info_int_2D

  subroutine print_field_grid_info_int_3D(self)
    class(type_field_grid_int_3D), intent(in) :: self
    call print_field_info_general( self%name, '3D', 'integer', 'grid', self%parent%name, &
      self%Arakawa_grid, self%long_name, self%units)
  end subroutine print_field_grid_info_int_3D

  subroutine print_field_grid_info_dp_2D(self)
    class(type_field_grid_dp_2D), intent(in) :: self
    call print_field_info_general( self%name, '2D', 'double precision', 'grid', self%parent%name, &
      self%Arakawa_grid, self%long_name, self%units)
  end subroutine print_field_grid_info_dp_2D

  subroutine print_field_grid_info_dp_3D(self)
    class(type_field_grid_dp_3D), intent(in) :: self
    call print_field_info_general( self%name, '3D', 'double precision', 'grid', self%parent%name, &
      self%Arakawa_grid, self%long_name, self%units)
  end subroutine print_field_grid_info_dp_3D

  subroutine print_field_mesh_info_logical_2D(self)
    class(type_field_mesh_logical_2D), intent(in) :: self
    call print_field_info_general( self%name, '2D', 'logical', 'mesh', self%parent%name, &
      self%Arakawa_grid, self%long_name, self%units)
  end subroutine print_field_mesh_info_logical_2D

  subroutine print_field_mesh_info_logical_3D(self)
    class(type_field_mesh_logical_3D), intent(in) :: self
    call print_field_info_general( self%name, '3D', 'logical', 'mesh', self%parent%name, &
      self%Arakawa_grid, self%long_name, self%units)
  end subroutine print_field_mesh_info_logical_3D

  subroutine print_field_mesh_info_int_2D(self)
    class(type_field_mesh_int_2D), intent(in) :: self
    call print_field_info_general( self%name, '2D', 'integer', 'mesh', self%parent%name, &
      self%Arakawa_grid, self%long_name, self%units)
  end subroutine print_field_mesh_info_int_2D

  subroutine print_field_mesh_info_int_3D(self)
    class(type_field_mesh_int_3D), intent(in) :: self
    call print_field_info_general( self%name, '3D', 'integer', 'mesh', self%parent%name, &
      self%Arakawa_grid, self%long_name, self%units)
  end subroutine print_field_mesh_info_int_3D

  subroutine print_field_mesh_info_dp_2D(self)
    class(type_field_mesh_dp_2D), intent(in) :: self
    call print_field_info_general( self%name, '2D', 'double precision', 'mesh', self%parent%name, &
      self%Arakawa_grid, self%long_name, self%units)
  end subroutine print_field_mesh_info_dp_2D

  subroutine print_field_mesh_info_dp_3D(self)
    class(type_field_mesh_dp_3D), intent(in) :: self
    call print_field_info_general( self%name, '3D', 'double precision', 'mesh', self%parent%name, &
      self%Arakawa_grid, self%long_name, self%units)
  end subroutine print_field_mesh_info_dp_3D

  subroutine print_field_info_general( field_name, field_dimension, field_type, &
    field_parent_type, field_parent_name, field_Arakawa_grid, field_long_name, field_units)

    ! In/output variables:
    character(len=*),        intent(in) :: field_name
    character(len=*),        intent(in) :: field_dimension
    character(len=*),        intent(in) :: field_type
    character(len=*),        intent(in) :: field_parent_type
    character(len=*),        intent(in) :: field_parent_name
    type(type_Arakawa_grid), intent(in) :: field_Arakawa_grid
    character(len=*),        intent(in) :: field_long_name
    character(len=*),        intent(in) :: field_units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'print_field_info_general'

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) then
      write(0,*) '    Field: ', colour_string( trim( field_name),'light blue')
      write(0,*) '      Long name   : ', trim( field_long_name)
      write(0,*) '      Dimension   : ', trim( field_dimension)
      write(0,*) '      Type        : ', trim( field_type)
      write(0,*) '      Parent type : ', trim( field_parent_type)
      write(0,*) '      Parent name : ', trim( field_parent_name)
      write(0,*) '      Arakawa grid: ', trim( field_Arakawa_grid%str())
      write(0,*) '      Units       : ', trim( field_units)
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine print_field_info_general

end submodule fields_info_submod