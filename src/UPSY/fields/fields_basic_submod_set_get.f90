submodule (fields_basic) fields_basic_submod_set_get

contains

  subroutine set_field_metadata( field, name, long_name, units)

    ! In/output variables:
    class(atype_field), intent(inout) :: field
    character(len=*),   intent(in   ) :: name
    character(len=*),   intent(in   ) :: long_name
    character(len=*),   intent(in   ) :: units

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_field_metadata'

    ! Add routine to call stack
    call init_routine( routine_name)

    field%name      = name
    field%long_name = long_name
    field%units     = units

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_field_metadata

  subroutine set_field_parent_grid( field, grid, field_Arakawa_grid)

    ! In/output variables:
    class(atype_field_grid), intent(inout) :: field
    type(type_grid), target, intent(in   ) :: grid
    type(type_Arakawa_grid), intent(in   ) :: field_Arakawa_grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_field_parent_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    field%parent       => grid
    field%Arakawa_grid = field_Arakawa_grid

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_field_parent_grid

  subroutine set_field_parent_mesh( field, mesh, field_Arakawa_grid)

    ! In/output variables:
    class(atype_field_mesh), intent(inout) :: field
    type(type_mesh), target, intent(in   ) :: mesh
    type(type_Arakawa_grid), intent(in   ) :: field_Arakawa_grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_field_parent_mesh'

    ! Add routine to call stack
    call init_routine( routine_name)

    field%parent       => mesh
    field%Arakawa_grid = field_Arakawa_grid

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_field_parent_mesh

end submodule fields_basic_submod_set_get