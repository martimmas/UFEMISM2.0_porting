module demo_model_a

  use precisions, only: dp
  use parameters, only: pi, NaN
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use demo_model, only: atype_demo_model, type_demo_model_context_allocate, &
    type_demo_model_context_initialise
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: type_demo_model_a

  type, extends(atype_demo_model) :: type_demo_model_a

      ! Some additional ice-model-esque data fields specific to demo_model_a
      real(dp), dimension(:  ), contiguous, pointer :: till_friction_angle => null()
      type(MPI_WIN) :: wtill_friction_angle

    contains

      procedure, public :: allocate_demo_model   => allocate_demo_model_a_abs
      procedure, public :: initialise_demo_model => initialise_demo_model_a_abs

  end type type_demo_model_a

contains

  subroutine allocate_demo_model_a_abs( self, context)

    ! In/output variables:
    class(type_demo_model_a),                       intent(inout) :: self
    type(type_demo_model_context_allocate), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_demo_model_a_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call allocate_demo_model_a( self%mesh, self)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_demo_model_a_abs

  subroutine allocate_demo_model_a( mesh, self)

    ! In/output variables:
    type(type_mesh),         intent(in   ) :: mesh
    type(type_demo_model_a), intent(inout) :: self

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'allocate_demo_model_a'

    ! Add routine to call stack
    call init_routine( routine_name)

    call self%create_field( self%till_friction_angle, self%wtill_friction_angle, &
      mesh, Arakawa_grid%a(), &
      name      = 'till_friction_angle', &
      long_name = 'till friction angle', &
      units     = 'degrees', &
      remap_method = 'reallocate')

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine allocate_demo_model_a

  subroutine initialise_demo_model_a_abs( self, context)

    ! In/output variables:
    class(type_demo_model_a),                         intent(inout) :: self
    type(type_demo_model_context_initialise), target, intent(in   ) :: context

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_demo_model_a_abs'

    ! Add routine to call stack
    call init_routine( routine_name)

    ! Retrieve input variables from context object
    call initialise_demo_model_a( self%mesh, self, context%till_friction_angle_uniform)

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_demo_model_a_abs

  subroutine initialise_demo_model_a( mesh, self, till_friction_angle_uniform)

    ! In/output variables:
    type(type_mesh),         intent(in   ) :: mesh
    type(type_demo_model_a), intent(inout) :: self
    real(dp),                intent(in   ) :: till_friction_angle_uniform

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_demo_model_a'

    ! Add routine to call stack
    call init_routine( routine_name)

    self%till_friction_angle( mesh%vi1: mesh%vi2) = till_friction_angle_uniform

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine initialise_demo_model_a

end module demo_model_a
