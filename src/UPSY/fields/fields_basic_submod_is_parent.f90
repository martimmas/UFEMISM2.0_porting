submodule (fields_basic) fields_basic_submod_is_parent

contains

  function is_parent( field, parent) result( res)
    !< Test if a certain object is identical to the parent object of this field
    !< Can be an x/y-grid, a mesh, an Arakawa grid, or a third dimension

    ! In/output variables:
    class(atype_field), intent(in) :: field
    class(*),           intent(in) :: parent
    logical                        :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_parent'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => parent)
    class is (type_grid)
      res = is_parent_grid( field, p)
    class is (type_mesh)
      res = is_parent_mesh( field, p)
    class is (type_Arakawa_grid)
      res = is_parent_Arakawa_grid( field, p)
    class is (type_third_dimension)
      res = is_parent_third_dimension( field, p)
    class default
      res = .false.
      call crash('invalid parent class')
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_parent

  function is_parent_grid( field, grid) result( res)
    !< Test if a certain x/y-grid is the parent x/y-grid of this field

    ! In/output variables:
    class(atype_field), intent(in) :: field
    type(type_grid),    intent(in) :: grid
    logical                        :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_parent_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => field)
    class is (atype_field_grid)
      res = p%parent%name == grid%name
    class default
      res = .false.
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_parent_grid

  function is_parent_mesh( field, mesh) result( res)
    !< Test if a certain mesh is the parent mesh of this field

    ! In/output variables:
    class(atype_field), intent(in) :: field
    type(type_mesh),    intent(in) :: mesh
    logical                        :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_parent_mesh'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => field)
    class is (atype_field_mesh)
      res = p%parent%name == mesh%name
    class default
      res = .false.
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_parent_mesh

  function is_parent_Arakawa_grid( field, field_Arakawa_grid) result( res)
    !< Test if a certain Arakawa grid is the parent Arakawa grid of this field

    ! In/output variables:
    class(atype_field),      intent(in) :: field
    type(type_Arakawa_grid), intent(in) :: field_Arakawa_grid
    logical                             :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_parent_Arakawa_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => field)
    class is (atype_field_grid)
      res = p%Arakawa_grid() == field_Arakawa_grid
    class is (atype_field_mesh)
      res = p%Arakawa_grid() == field_Arakawa_grid
    class default
      res = .false.
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_parent_Arakawa_grid

  function is_parent_third_dimension( field, field_third_dimension) result( res)
    !< Test if a certain third dimension is the parent third dimension of this field

    ! In/output variables:
    class(atype_field),         intent(in) :: field
    type(type_third_dimension), intent(in) :: field_third_dimension
    logical                                :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_third_dimension'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type(p => field)
    class default
      res = .false.
      call crash('class/type not implemented')
    class is (type_field_grid_logical_2D)
      res = .false.
    class is (type_field_grid_logical_3D)
      res = p%third_dimension() == field_third_dimension
    class is (type_field_grid_int_2D)
      res = .false.
    class is (type_field_grid_int_3D)
      res = p%third_dimension() == field_third_dimension
    class is (type_field_grid_dp_2D)
      res = .false.
    class is (type_field_grid_dp_3D)
      res = p%third_dimension() == field_third_dimension
    class is (type_field_mesh_logical_2D)
      res = .false.
    class is (type_field_mesh_logical_3D)
      res = p%third_dimension() == field_third_dimension
    class is (type_field_mesh_int_2D)
      res = .false.
    class is (type_field_mesh_int_3D)
      res = p%third_dimension() == field_third_dimension
    class is (type_field_mesh_dp_2D)
      res = .false.
    class is (type_field_mesh_dp_3D)
      res = p%third_dimension() == field_third_dimension
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_parent_third_dimension

end submodule fields_basic_submod_is_parent