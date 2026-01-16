submodule (fields_basic) fields_basic_submod_remap

contains

  subroutine remap( self, mesh_new)

    ! In/output variables:
    class(atype_field),      intent(inout) :: self
    type(type_mesh), target, intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_field_remap'
    type(type_mesh), pointer       :: mesh_old

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (g => self%grid())
    class default
      call crash('remapping only defined for mesh-based fields')
    class is (type_mesh)
      mesh_old => g
    end select

    select type (f => self)
    class default
      call crash('invalid field class')
    class is (type_field_dp_2D)
      call remap_field_dp_2D( f, mesh_old, mesh_new)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap

  subroutine remap_field_dp_2D( field, mesh_old, mesh_new)

    ! In/output variables:
    type(type_field_dp_2D),  intent(inout) :: field
    type(type_mesh), target, intent(in   ) :: mesh_old, mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'remap_field_dp_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_field_dp_2D

end submodule fields_basic_submod_remap