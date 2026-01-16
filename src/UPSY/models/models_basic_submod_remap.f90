submodule( models_basic) models_basic_submod_remap

contains

  subroutine remap( self, mesh_new)

    ! In/output variables:
    class(atype_model), intent(inout) :: self
    type(type_mesh),    intent(in   ) :: mesh_new

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_model_remap'
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (g => self%grid())
    class default
      call crash('remapping only defined for mesh-based models')
    class is (type_mesh)
    end select

    call self%flds_reg%remap( mesh_new)
    call self%set_grid( mesh_new)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap

end submodule models_basic_submod_remap
