module demo_model

  use precisions, only: dp
  use parameters, only: pi, NaN
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, warning, crash
  use mesh_types, only: type_mesh
  use Arakawa_grid_mod, only: Arakawa_grid
  use fields_main, only: third_dimension
  use models_basic, only: atype_model, atype_model_context_allocate
  use mpi_f08, only: MPI_WIN

  implicit none

  private

  public :: atype_demo_model, type_demo_model_context_allocate

  type, abstract, extends(atype_model) :: atype_demo_model

      ! Some ice-model-esque data fields
      real(dp), dimension(:  ), contiguous, pointer :: H        => null()
      real(dp), dimension(:,:), contiguous, pointer :: u_3D     => null()
      real(dp), dimension(:,:), contiguous, pointer :: v_3D     => null()
      logical,  dimension(:  ), contiguous, pointer :: mask_ice => null()
      real(dp), dimension(:,:), contiguous, pointer :: T2m      => null()
      type(MPI_WIN) :: wH, wu_3D, wv_3D, wmask_ice, wT2m

    contains

      procedure, public :: allocate_model => allocate_model_abs

      procedure(allocate_demo_model_ifc), deferred :: allocate_demo_model

      procedure, nopass, public :: ct_allocate

  end type atype_demo_model

  ! Context classes for allocate/initialise/run/remap
  ! =================================================

  type, extends(atype_model_context_allocate) :: type_demo_model_context_allocate
    integer  :: nz
  end type type_demo_model_context_allocate

  ! Abstract interfaces for deferred procedures
  ! ===========================================

  abstract interface

    subroutine allocate_demo_model_ifc( self, context)
      import atype_demo_model, type_demo_model_context_allocate
      class(atype_demo_model),                        intent(inout) :: self
      type(type_demo_model_context_allocate), target, intent(in   ) :: context
    end subroutine allocate_demo_model_ifc

  end interface

  ! Interfaces to type-bound procedures defined in submodules
  ! =========================================================

  interface

    module subroutine allocate_model_abs( self, context)
      class(atype_demo_model),                     intent(inout) :: self
      class(atype_model_context_allocate), target, intent(in   ) :: context
    end subroutine allocate_model_abs

    module function ct_allocate( name, region_name, mesh, nz) result( context)
      character(len=*),           intent(in) :: name
      character(len=*),           intent(in) :: region_name
      type(type_mesh), target,    intent(in) :: mesh
      integer,                    intent(in) :: nz
      type(type_demo_model_context_allocate) :: context
    end function ct_allocate

  end interface

  ! subroutine initialise_demo_model( model, mesh, nz)

  !   ! In/output variables:
  !   class(type_demo_model),  intent(  out) :: model
  !   type(type_mesh), target, intent(in   ) :: mesh
  !   integer,                 intent(in   ) :: nz

  !   ! Local variables:
  !   character(len=1024), parameter :: routine_name = 'initialise_demo_model'
  !   integer                        :: vi,ti,k,m
  !   real(dp)                       :: x,y,cx,cy

  !   ! Add routine to call stack
  !   call init_routine( routine_name)

  !   ! Initialise fields with simple analytical functions

  !   cx = mesh%xmax - mesh%xmin
  !   cy = mesh%ymax - mesh%ymin

  !   do vi = mesh%vi1, mesh%vi2

  !     x = mesh%V( vi,1)
  !     y = mesh%V( vi,2)

  !     model%H( vi) = max( 0._dp, cos( x * pi / cx) * cos( y * pi / cy) - 0.2_dp)
  !     model%mask_ice( vi) = model%H( vi) > 0._dp

  !     do m = 1, 12
  !       model%T2m( vi,m) = hypot( x,y) + sin( real(m,dp) * 2._dp * pi / 12._dp)
  !     end do

  !   end do

  !   do ti = mesh%ti1, mesh%ti2

  !     x = mesh%Trigc( ti,1)
  !     y = mesh%Trigc( ti,2)

  !     do k = 1, nz
  !       model%u_3D( ti,k) = max( 0._dp, cos( x * pi / cx) * cos( y * pi / cy) - 0.2_dp) + real( k,dp)
  !       model%v_3D( ti,k) =             sin( x * pi / cx) * sin( y * pi / cy)           + real( k,dp)
  !     end do

  !   end do

  !   ! Remove routine from call stack
  !   call finalise_routine( routine_name)

  ! end subroutine initialise_demo_model

end module demo_model
