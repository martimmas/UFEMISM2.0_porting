submodule( fields_registry) fields_registry_submod_basic

contains

  subroutine add_field_to_registry( flds_reg, field)

    ! In/output variables:
    class(type_fields_registry), intent(inout) :: flds_reg
    class(atype_field),          intent(in   ) :: field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_field_to_registry'
    logical                        :: is_in_use
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    ! If the registry is not allocated yet, allocate it.
    if (.not. allocated( flds_reg%items)) then
      flds_reg%n     = 1
      flds_reg%n_max = 1
      allocate( flds_reg%items(1))
      allocate( flds_reg%items(1)%p, source = field)
      call finalise_routine( routine_name)
      return
    end if

    ! If the registry is full, extend it
    if (flds_reg%n == flds_reg%n_max) then
      call flds_reg%extend
    end if

    ! Check that this name is not already in use
    is_in_use = .false.
    do i = 1, flds_reg%n
      is_in_use = is_in_use .or. flds_reg%items(i)%p%name() == field%name()
    end do
    if (is_in_use) call crash('a field of name "' // trim( field%name()) // &
      '" already exists in this registry')

    ! Add field to registry
    flds_reg%n = flds_reg%n+1
    allocate( flds_reg%items( flds_reg%n)%p, source = field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine add_field_to_registry

  subroutine extend_field_registry( flds_reg)

    ! In/output variables:
    class(type_fields_registry), intent(inout) :: flds_reg

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'extend_field_registry'
    type(type_field_box), allocatable :: tmp(:)
    integer                           :: i

    ! Add routine to call stack
    call init_routine( routine_name)


    allocate( tmp( flds_reg%n))
    do i = 1, flds_reg%n
      allocate( tmp(i)%p, source = flds_reg%items(i)%p)
    end do

    flds_reg%n_max = flds_reg%n_max * 2
    deallocate( flds_reg%items)
    allocate( flds_reg%items( flds_reg%n_max))

    do i = 1, flds_reg%n
      allocate( flds_reg%items(i)%p, source = tmp(i)%p)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine extend_field_registry

  subroutine print_info( flds_reg)

    ! In/output variables:
    class(type_fields_registry), intent(in) :: flds_reg

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'print_info'
    integer                        :: i

    if (.not. allocated( flds_reg%items)) return

    ! Add routine to call stack
    call init_routine( routine_name)

    do i = 1, flds_reg%n
       call flds_reg%items(i)%p%print_info
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine print_info

  function find_field_by_name( flds_reg, name) result(i)

    ! In/output variables:
    class(type_fields_registry), intent(in) :: flds_reg
    character(len=*),            intent(in) :: name
    integer                                 :: i

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'find_field_by_name'
    logical                        :: found_it
    integer                        :: ii

    ! Add routine to call stack
    call init_routine( routine_name)

    if (.not. allocated(flds_reg%items)) then
      call crash('field collection not allocated')
    end if

    i = 0
    found_it = .false.
    do ii = 1, size( flds_reg%items)
      if (allocated( flds_reg%items(ii)%p)) then
        if (flds_reg%items(ii)%p%name() == name) then
          found_it = .true.
          i = ii
          exit
        end if
      end if
    end do

    if (.not. found_it) then
      call crash('could not find field "' // trim( name) // '" in field collection')
    end if

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end function find_field_by_name

  function test_fields_registry_equality( flds_reg1, flds_reg2) result( res)

    ! In/output variables:
    class(type_fields_registry), intent(in) :: flds_reg1, flds_reg2
    logical                                 :: res

    ! Local variables:
    integer :: i

    res = flds_reg1%n == flds_reg2%n
    if (.not. res) return
    do i = 1, flds_reg1%n
      res = res .and. flds_reg1%items(i)%p == flds_reg2%items(i)%p
    end do

  end function test_fields_registry_equality

end submodule fields_registry_submod_basic