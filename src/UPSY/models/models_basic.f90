module models_basic

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine
  use fields_main, only: type_fields_registry

  implicit none

  private

  public :: atype_model

  ! Abstract basic model type
  ! =========================

  type, abstract :: atype_model

    ! Metadata
    character(len=1024), private :: name_val

    ! Fields registry
    type(type_fields_registry), public :: flds_reg

    contains

    ! ===== Set/get functions

    ! Metadata
    procedure, public :: set_name
    procedure, public :: name => get_name
    procedure, public :: is_name

  end type atype_model

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  interface

    ! ===== Set/get functions

    ! Metadata

    module subroutine set_name( model, name)
      class(atype_model), intent(inout) :: model
      character(len=*),   intent(in   ) :: name
    end subroutine set_name

    module function get_name( model) result( name)
      class(atype_model), intent(in) :: model
      character(:), allocatable      :: name
    end function get_name

    module function is_name( model, name) result( res)
      class(atype_model), intent(in) :: model
      character(len=*),   intent(in) :: name
      logical                        :: res
    end function is_name

  end interface

end module models_basic
