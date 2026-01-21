submodule( models_basic) models_basic_submod_set_get

contains

  ! Metadata

  subroutine set_name( self, name)
    class(atype_model), intent(inout) :: self
    character(len=*),   intent(in   ) :: name
    self%name_val = name
  end subroutine set_name

  function get_name( self) result( name)
    class(atype_model), intent(in) :: self
    character(:), allocatable      :: name
    name = self%name_val
  end function get_name

  function is_name( self, name) result( res)
    class(atype_model), intent(in) :: self
    character(len=*),   intent(in) :: name
    logical                        :: res
    res = self%name_val == name
  end function is_name

  ! Grid

  subroutine set_grid( self, grid)

    ! In/output variables:
    class(atype_model), intent(inout) :: self
    class(*), target,   intent(in   ) :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'set_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (p => grid)
    class default
      call crash('invalid grid class')
    class is (type_grid)
      self%grid_val => grid
    class is (type_mesh)
      self%grid_val => grid
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine set_grid

  function get_grid( self) result( grid)

    ! In/output variables:
    class(atype_model), intent(in) :: self
    class(*), pointer              :: grid

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'get_grid'

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (p => self%grid_val)
    class default
      call crash('invalid grid class')
    class is (type_grid)
      grid => self%grid_val
    class is (type_mesh)
      grid => self%grid_val
    end select

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function get_grid

  function is_grid( self, grid) result( res)

    ! In/output variables:
    class(atype_model), intent(in) :: self
    class(*),           intent(in) :: grid
    logical                        :: res

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'is_grid'
    character(len=1024)            :: name1, name2

    ! Add routine to call stack
    call init_routine( routine_name)

    select type (p => self%grid_val)
    class default
      call crash('invalid model%grid class')
    class is (type_grid)
      name1 = p%name
    class is (type_mesh)
      name1 = p%name
    end select

    select type (p => grid)
    class default
      call crash('invalid grid class')
    class is (type_grid)
      name2 = p%name
    class is (type_mesh)
      name2 = p%name
    end select

    res = name1 == name2

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function is_grid

end submodule models_basic_submod_set_get
