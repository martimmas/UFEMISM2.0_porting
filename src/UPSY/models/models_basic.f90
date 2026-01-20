module models_basic

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, crash
  use grid_types, only: type_grid
  use mesh_types, only: type_mesh
  use fields_main, only: type_fields_registry
  use Arakawa_grid_mod, only: type_Arakawa_grid
  use fields_dimensions, only: type_third_dimension
  use mpi_f08, only: MPI_WIN

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
      type(type_fields_registry), private :: flds_reg

    contains

      generic,   public  :: create_field => &
        create_field_logical_2D, create_field_int_2D, create_field_dp_2D, &
        create_field_logical_3D, create_field_int_3D, create_field_dp_3D
      procedure, private :: create_field_logical_2D
      procedure, private :: create_field_int_2D
      procedure, private :: create_field_dp_2D
      procedure, private :: create_field_logical_3D
      procedure, private :: create_field_int_3D
      procedure, private :: create_field_dp_3D

      generic,   public  :: reallocate_field => &
        reallocate_field_logical_2D, &
        reallocate_field_logical_3D, &
        reallocate_field_int_2D, &
        reallocate_field_int_3D, &
        reallocate_field_dp_2D, &
        reallocate_field_dp_3D
      procedure, private :: reallocate_field_logical_2D
      procedure, private :: reallocate_field_logical_3D
      procedure, private :: reallocate_field_int_2D
      procedure, private :: reallocate_field_int_3D
      procedure, private :: reallocate_field_dp_2D
      procedure, private :: reallocate_field_dp_3D

      procedure, public :: write_to_restart_file
      procedure, public :: read_from_restart_file

      ! ===== Basics

      generic,   public  :: operator(==) => eq
      procedure, private :: eq => test_model_equality

      ! Metadata
      procedure, public :: set_name
      procedure, public :: name => get_name
      procedure, public :: is_name

      ! Grid
      procedure, public :: set_grid
      procedure, public :: grid => get_grid
      procedure, public :: is_grid

  end type atype_model

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  ! create_field
  interface

    module subroutine create_field_logical_2D( self, d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units, remap_method)
      class(atype_model),                         intent(inout) :: self
      logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                              intent(inout) :: w
      class(*), target,                           intent(in   ) :: field_grid
      type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
      character(len=*),                 optional, intent(in   ) :: name
      character(len=*),                 optional, intent(in   ) :: long_name
      character(len=*),                 optional, intent(in   ) :: units
      character(len=*),                 optional, intent(in   ) :: remap_method
    end subroutine create_field_logical_2D

    module subroutine create_field_logical_3D( self, d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)
      class(atype_model),                           intent(inout) :: self
      logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                                intent(inout) :: w
      class(*), target,                             intent(in   ) :: field_grid
      type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
      type(type_third_dimension),                   intent(in   ) :: field_third_dimension
      character(len=*),                   optional, intent(in   ) :: name
      character(len=*),                   optional, intent(in   ) :: long_name
      character(len=*),                   optional, intent(in   ) :: units
      character(len=*),                   optional, intent(in   ) :: remap_method
    end subroutine create_field_logical_3D

    module subroutine create_field_int_2D( self, d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units, remap_method)
      class(atype_model),                         intent(inout) :: self
      integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                              intent(inout) :: w
      class(*), target,                           intent(in   ) :: field_grid
      type(type_Arakawa_grid),                    intent(in   ) :: field_Arakawa_grid
      character(len=*),                 optional, intent(in   ) :: name
      character(len=*),                 optional, intent(in   ) :: long_name
      character(len=*),                 optional, intent(in   ) :: units
      character(len=*),                 optional, intent(in   ) :: remap_method
    end subroutine create_field_int_2D

    module subroutine create_field_int_3D( self, d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)
      class(atype_model),                           intent(inout) :: self
      integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                                intent(inout) :: w
      class(*), target,                             intent(in   ) :: field_grid
      type(type_Arakawa_grid),                      intent(in   ) :: field_Arakawa_grid
      type(type_third_dimension),                   intent(in   ) :: field_third_dimension
      character(len=*),                   optional, intent(in   ) :: name
      character(len=*),                   optional, intent(in   ) :: long_name
      character(len=*),                   optional, intent(in   ) :: units
      character(len=*),                   optional, intent(in   ) :: remap_method
    end subroutine create_field_int_3D

    module subroutine create_field_dp_2D( self, d_nih, w, field_grid, &
      field_Arakawa_grid, name, long_name, units, remap_method)
      class(atype_model),                          intent(inout) :: self
      real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                               intent(inout) :: w
      class(*), target,                            intent(in   ) :: field_grid
      type(type_Arakawa_grid),                     intent(in   ) :: field_Arakawa_grid
      character(len=*),                  optional, intent(in   ) :: name
      character(len=*),                  optional, intent(in   ) :: long_name
      character(len=*),                  optional, intent(in   ) :: units
      character(len=*),                  optional, intent(in   ) :: remap_method
    end subroutine create_field_dp_2D

    module subroutine create_field_dp_3D( self, d_nih, w, field_grid, &
      field_Arakawa_grid, field_third_dimension, name, long_name, units, remap_method)
      class(atype_model),                            intent(inout) :: self
      real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
      type(MPI_WIN),                                 intent(inout) :: w
      class(*), target,                              intent(in   ) :: field_grid
      type(type_Arakawa_grid),                       intent(in   ) :: field_Arakawa_grid
      type(type_third_dimension),                    intent(in   ) :: field_third_dimension
      character(len=*),                    optional, intent(in   ) :: name
      character(len=*),                    optional, intent(in   ) :: long_name
      character(len=*),                    optional, intent(in   ) :: units
      character(len=*),                    optional, intent(in   ) :: remap_method
    end subroutine create_field_dp_3D

  end interface

  ! reallocate
  interface

    module subroutine reallocate_field_logical_2D( self, mesh_new, field_name, d_nih)
      class(atype_model),                         intent(inout) :: self
      character(len=*),                           intent(in   ) :: field_name
      type(type_mesh),                            intent(in   ) :: mesh_new
      logical, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_logical_2D

    module subroutine reallocate_field_logical_3D( self, mesh_new, field_name, d_nih)
      class(atype_model),                           intent(inout) :: self
      character(len=*),                             intent(in   ) :: field_name
      type(type_mesh),                              intent(in   ) :: mesh_new
      logical, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_logical_3D

    module subroutine reallocate_field_int_2D( self, mesh_new, field_name, d_nih)
      class(atype_model),                         intent(inout) :: self
      character(len=*),                           intent(in   ) :: field_name
      type(type_mesh),                            intent(in   ) :: mesh_new
      integer, dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_int_2D

    module subroutine reallocate_field_int_3D( self, mesh_new, field_name, d_nih)
      class(atype_model),                           intent(inout) :: self
      character(len=*),                             intent(in   ) :: field_name
      type(type_mesh),                              intent(in   ) :: mesh_new
      integer, dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_int_3D

    module subroutine reallocate_field_dp_2D( self, mesh_new, field_name, d_nih)
      class(atype_model),                          intent(inout) :: self
      character(len=*),                            intent(in   ) :: field_name
      type(type_mesh),                             intent(in   ) :: mesh_new
      real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_dp_2D

    module subroutine reallocate_field_dp_3D( self, mesh_new, field_name, d_nih)
      class(atype_model),                            intent(inout) :: self
      character(len=*),                              intent(in   ) :: field_name
      type(type_mesh),                               intent(in   ) :: mesh_new
      real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih
    end subroutine reallocate_field_dp_3D

  end interface

  ! NetCDF i/o
  interface

    module subroutine write_to_restart_file( self, output_dir, filename)
      class(atype_model),                  intent(in   ) :: self
      character(len=*),                    intent(in   ) :: output_dir
      character(:), allocatable, optional, intent(  out) :: filename
    end subroutine write_to_restart_file

    module subroutine read_from_restart_file( self, filename)
      class(atype_model), intent(inout) :: self
      character(len=*),   intent(in   ) :: filename
    end subroutine read_from_restart_file

  end interface

  ! Basics
  interface

    module function test_model_equality( model1, model2) result( res)
      class(atype_model), intent(in) :: model1, model2
      logical                        :: res
    end function test_model_equality

    ! ===== Set/get functions

    ! Metadata

    module subroutine set_name( self, name)
      class(atype_model), intent(inout) :: self
      character(len=*),   intent(in   ) :: name
    end subroutine set_name

    module function get_name( self) result( name)
      class(atype_model), intent(in) :: self
      character(:), allocatable      :: name
    end function get_name

    module function is_name( self, name) result( res)
      class(atype_model), intent(in) :: self
      character(len=*),   intent(in) :: name
      logical                        :: res
    end function is_name

    ! Grid

    module subroutine set_grid( self, grid)
      class(atype_model), intent(inout) :: self
      class(*), target,   intent(in   ) :: grid
    end subroutine set_grid

    module function get_grid( self) result( grid)
      class(atype_model), intent(in) :: self
      class(*), pointer              :: grid
    end function get_grid

    module function is_grid( self, grid) result( res)
      class(atype_model), intent(in) :: self
      class(*),           intent(in) :: grid
      logical                        :: res
    end function is_grid

  end interface

end module models_basic
