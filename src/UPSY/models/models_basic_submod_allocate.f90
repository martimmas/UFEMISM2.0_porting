submodule(models_basic) models_basic_submod_allocate

contains

  subroutine allocate( self, context)
    !< Allocate an instance of a model

    ! In/output variables:
    class(atype_model),                          intent(inout) :: self
    class(atype_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_model_allocate'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Set model metadata and mesh
    call self%set_name       ( context%name)
    call self%set_region_name( context%region_name)
    self%mesh => context%mesh

    ! Call the model-specific allocation routine
    call self%allocate_model( context)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate

end submodule models_basic_submod_allocate