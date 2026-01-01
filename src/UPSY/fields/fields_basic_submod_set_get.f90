submodule (fields_basic) fields_basic_submod_set_get

contains

  subroutine set_metadata( field, name, long_name, units)

    ! In/output variables:
    class(atype_field), intent(inout) :: field
    character(len=*),   intent(in   ) :: name
    character(len=*),   intent(in   ) :: long_name
    character(len=*),   intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_metadata'

    ! Add routine to call stack
    call init_routine( routine_name)

    field%name_val      = name
    field%long_name_val = long_name
    field%units_val     = units

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_metadata

  function get_name( field) result( name)
    class(atype_field), intent(in) :: field
    character(:), allocatable      :: name
    name = field%name_val
  end function get_name

  function get_long_name( field) result( long_name)
    class(atype_field), intent(in) :: field
    character(:), allocatable      :: long_name
    long_name = field%long_name_val
  end function get_long_name

  function get_units( field) result( units)
    class(atype_field), intent(in) :: field
    character(:), allocatable      :: units
    units = field%units_val
  end function get_units

  subroutine set_parent_grid( field, grid)

    ! In/output variables:
    class(atype_field_grid), intent(inout) :: field
    type(type_grid), target, intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_parent_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    field%parent => grid

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_parent_grid

  subroutine set_parent_mesh( field, mesh)

    ! In/output variables:
    class(atype_field_mesh), intent(inout) :: field
    type(type_mesh), target, intent(in   ) :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_parent_mesh'

    ! Add routine to call stack
    call init_routine( routine_name)

    field%parent => mesh

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_parent_mesh

  subroutine set_parent_Arakawa_grid( field, field_Arakawa_grid)

    ! In/output variables:
    class(atype_field),      intent(inout) :: field
    type(type_Arakawa_grid), intent(in   ) :: field_Arakawa_grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_parent_Arakawa_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    field%Arakawa_grid = field_Arakawa_grid

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_parent_Arakawa_grid

  subroutine set_third_dimension( field, field_third_dimension)

    ! In/output variables:
    class(atype_field),         intent(inout) :: field
    type(type_third_dimension), intent(in   ) :: field_third_dimension

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_third_dimension'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (p => field)
    class default
      call crash('invalid field type')
    class is (type_field_grid_logical_3D)
      p%third_dimension = field_third_dimension
    class is (type_field_grid_int_3D)
      p%third_dimension = field_third_dimension
    class is (type_field_grid_dp_3D)
      p%third_dimension = field_third_dimension
    class is (type_field_mesh_logical_3D)
      p%third_dimension = field_third_dimension
    class is (type_field_mesh_int_3D)
      p%third_dimension = field_third_dimension
    class is (type_field_mesh_dp_3D)
      p%third_dimension = field_third_dimension
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_third_dimension

end submodule fields_basic_submod_set_get