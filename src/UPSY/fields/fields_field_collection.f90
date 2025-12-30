module fields_field_collection

  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use fields_basic, only: atype_field
  use mpi_f08, only: MPI_WIN, MPI_GATHER, MPI_INTEGER, MPI_COMM_WORLD

  implicit none

  private

  public :: &
    type_field_collection, add_initialised_field_to_collection, find_field_by_name

  ! Wrapper so we can have a mixed-type array
  ! =========================================

  type :: type_field_box
     class(atype_field), allocatable :: p
  end type type_field_box

  ! Field collection
  ! ================

  type :: type_field_collection
     type(type_field_box), allocatable :: items(:)
     integer                           :: n     = 0
     integer                           :: n_max = 0
   contains
     procedure :: add_initialised_field_to_collection
     procedure :: print_fields_info
     procedure :: find_field_by_name
  end type type_field_collection

contains

  subroutine add_initialised_field_to_collection( bof, field)

    ! In/output variables:
    class(type_field_collection), intent(inout) :: bof
    class(atype_field),           intent(in   ) :: field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'add_initialised_field_to_collection'
    logical                        :: is_in_use
    integer                        :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    ! If the collection is not alocated yet, allocate it.
    if (.not. allocated( bof%items)) then
      bof%n     = 1
      bof%n_max = 1
      allocate( bof%items(1))
      allocate( bof%items(1)%p, source = field)
      call finalise_routine( routine_name)
      return
    end if

    ! If the collection is full, extend it
    if (bof%n == bof%n_max) then
      bof%n_max = bof%n_max * 2
      call extend_field_collection( bof%items, bof%n, bof%n_max)
    end if

    ! Check that this name is not already in use
    is_in_use = .false.
    do i = 1, bof%n
      is_in_use = is_in_use .or. bof%items(i)%p%name == field%name
    end do
    if (is_in_use) call crash('a field of name "' // trim( field%name) // '" already exists')

    ! Add field to collection
    bof%n = bof%n+1
    allocate( bof%items( bof%n)%p, source = field)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine add_initialised_field_to_collection

  subroutine extend_field_collection( items, n_keep, n_max_new)

    ! In/output variables:
    type(type_field_box), allocatable, intent(inout) :: items(:)
    integer,                           intent(in   ) :: n_keep, n_max_new

    ! Local variables:
    character(len=1024), parameter    :: routine_name = 'extend_field_collection'
    type(type_field_box), allocatable :: tmp(:)
    integer                           :: i

    ! Add routine to call stack
    call init_routine( routine_name)

    allocate( tmp( n_keep))
    do i = 1, n_keep
      allocate( tmp(i)%p, source = items(i)%p)
    end do

    deallocate( items)
    allocate( items( n_max_new))

    do i = 1, n_keep
      allocate( items(i)%p, source = tmp(i)%p)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine extend_field_collection

  subroutine print_fields_info( bof)

    ! In/output variables:
    class(type_field_collection), intent(in) :: bof

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'print_fields_info'
    integer                        :: i

    if (.not. allocated( bof%items)) return

    ! Add routine to call stack
    call init_routine( routine_name)

    do i = 1, bof%n
       call bof%items(i)%p%print_field_info
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine print_fields_info

  function find_field_by_name( bof, name) result(i)

    ! In/output variables:
    class(type_field_collection), intent(in) :: bof
    character(len=*),             intent(in) :: name
    integer                                  :: i

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'find_field_by_name'
    logical                        :: found_it
    integer                        :: ii

    ! Add routine to call stack
    call init_routine( routine_name)

    if (.not. allocated(bof%items)) then
      call crash('field collection not allocated')
    end if

    i = 0
    found_it = .false.
    do ii = 1, size( bof%items)
      if (allocated( bof%items(ii)%p)) then
        if (bof%items(ii)%p%name == name) then
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

end module fields_field_collection