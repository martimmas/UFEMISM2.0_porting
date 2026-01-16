submodule (fields_basic) fields_basic_submod_reallocate

contains

  subroutine reallocate( self, mesh_new)

    ! In/output variables:
    class(atype_field),      intent(inout) :: self
    type(type_mesh), target, intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_field_reallocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (f => self)
    class default
      call crash('invalid field class')
    class is (type_field_dp_2D)
      call reallocate_field_dp_2D( f, mesh_new)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate

  subroutine reallocate_field_dp_2D( field, mesh_new)

    ! In/output variables:
    type(type_field_dp_2D),  intent(inout) :: field
    type(type_mesh), target, intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'reallocate_field_dp_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call crash('whoopsiedaisy')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine reallocate_field_dp_2D

end submodule fields_basic_submod_reallocate