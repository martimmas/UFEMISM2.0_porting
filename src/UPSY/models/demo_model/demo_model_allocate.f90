submodule(demo_model) demo_model_allocate

contains

  subroutine allocate_model_abs( self, context)

    ! In/output variables:
    class(atype_demo_model),                     intent(inout) :: self
    class(atype_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'atype_demo_model_allocate_model_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Downcast context class
    select type (ct => context)
    class default
      call crash('invalid context class; should be atype_demo_model_context_allocate')
    class is (atype_demo_model_context_allocate)
      call allocate_model( self, ct)
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_model_abs

  subroutine allocate_model( self, context)

    ! In/output variables:
    class(atype_demo_model),                          intent(inout) :: self
    class(atype_demo_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_model'

    ! Add routine to call stack
    call init_routine( routine_name)

    call allocate_model_common( self, context)  !< Part common to all models of type demo_model
    call self%allocate_demo_model( context)     !< Part specific to a model class inheriting from type_demo_model

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_model

  subroutine allocate_model_common( self, context)

    ! In/output variables:
    class(atype_demo_model),                          intent(inout) :: self
    class(atype_demo_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_model_common'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%create_field( self%H, self%wH, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'H', &
      long_name = 'ice thickness', &
      units     = 'm', &
      remap_method = '2nd_order_conservative')

    call self%create_field( self%u_3D, self%wu_3D, &
      self%mesh, Arakawa_grid%b(), third_dimension%ice_zeta( context%nz, 'regular'), &
      name      = 'u_3D', &
      long_name = 'depth-dependent horizontal ice velocity in x-direction', &
      units     = 'm yr^-1', &
      remap_method = '2nd_order_conservative')

    call self%create_field( self%v_3D, self%wv_3D, &
      self%mesh, Arakawa_grid%b(), third_dimension%ice_zeta( context%nz, 'regular'), &
      name      = 'v_3D', &
      long_name = 'depth-dependent horizontal ice velocity in y-direction', &
      units     = 'm yr^-1', &
      remap_method = '2nd_order_conservative')

    call self%create_field( self%mask_ice, self%wmask_ice, &
      self%mesh, Arakawa_grid%a(), &
      name      = 'mask_ice', &
      long_name = 'ice mask', &
      units     = '-', &
      remap_method = 'reallocate')

    call self%create_field( self%T2m, self%wT2m, &
      self%mesh, Arakawa_grid%a(), third_dimension%month(), &
      name      = 'T2m', &
      long_name = 'Monthly 2-m air temperature', &
      units     = 'K', &
      remap_method = '2nd_order_conservative')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_model_common

end submodule demo_model_allocate