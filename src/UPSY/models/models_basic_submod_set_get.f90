submodule( models_basic) models_basic_submod_set_get

contains

  ! Metadata

  subroutine set_name( model, name)
    class(atype_model), intent(inout) :: model
    character(len=*),   intent(in   ) :: name
    model%name_val = name
  end subroutine set_name

  function get_name( model) result( name)
    class(atype_model), intent(in) :: model
    character(:), allocatable      :: name
    name = model%name_val
  end function get_name

  function is_name( model, name) result( res)
    class(atype_model), intent(in) :: model
    character(len=*),   intent(in) :: name
    logical                        :: res
    res = model%name_val == name
  end function is_name

end submodule models_basic_submod_set_get
