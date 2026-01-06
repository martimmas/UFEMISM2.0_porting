module models_basic

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use fields_main, only: type_fields_registry

  implicit none

  private

  public :: atype_model

  ! Abstract basic model type
  ! =========================

  type, abstract :: atype_model

    ! Metadata
    character(len=1024), private :: name_val

    ! Grid
    class(*), pointer, private :: grid_val

    ! Fields registry
    type(type_fields_registry), public :: flds_reg

    contains

    ! ===== Set/get functions

    ! Metadata
    procedure, public :: set_name
    procedure, public :: name => get_name
    procedure, public :: is_name

    ! Grid
    procedure, public :: set_grid
    procedure, public :: grid => get_grid
    procedure, public :: is_grid

    ! ===== NetCDF output

    procedure, public :: write_to_restart_file

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

    ! Grid

    module subroutine set_grid( model, grid)
      class(atype_model), intent(inout) :: model
      class(*), target,   intent(in   ) :: grid
    end subroutine set_grid

    module function get_grid( model) result( grid)
      class(atype_model), intent(in) :: model
      class(*), pointer              :: grid
    end function get_grid

    module function is_grid( model, grid) result( res)
      class(atype_model), intent(in) :: model
      class(*),           intent(in) :: grid
      logical                        :: res
    end function is_grid

    ! ===== NetCDF output

    module subroutine write_to_restart_file( model, output_dir)
      class(atype_model), intent(in) :: model
      character(len=*),   intent(in) :: output_dir
    end subroutine write_to_restart_file

  end interface

end module models_basic
