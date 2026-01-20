submodule( fields_registry) fields_registry_submod_remap

contains

  subroutine remap_dp_2D( self, mesh_new, field_name, d_nih)

    ! In/output variables:
    class(type_fields_registry),                 intent(inout) :: self
    character(len=*),                            intent(in   ) :: field_name
    type(type_mesh),                             intent(in   ) :: mesh_new
    real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'type_fields_registry_remap_dp_2D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%items( self%find( field_name))%p%remap( mesh_new, d_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_dp_2D

  subroutine remap_dp_3D( self, mesh_new, field_name, d_nih)

    ! In/output variables:
    class(type_fields_registry),                   intent(inout) :: self
    character(len=*),                              intent(in   ) :: field_name
    type(type_mesh),                               intent(in   ) :: mesh_new
    real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'type_fields_registry_remap_dp_3D'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%items( self%find( field_name))%p%remap( mesh_new, d_nih)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine remap_dp_3D

end submodule fields_registry_submod_remap