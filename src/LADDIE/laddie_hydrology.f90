MODULE laddie_hydrology

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE control_resources_and_error_messaging                  , ONLY: crash, init_routine, finalise_routine, colour_string
  USE model_configuration                                    , ONLY: C
  USE parameters
  USE mesh_types                                             , ONLY: type_mesh
  USE laddie_model_types                                     , ONLY: type_laddie_model, type_laddie_timestep
  USE laddie_forcing_types                                   , ONLY: type_laddie_forcing
  USE reallocate_mod                                         , ONLY: reallocate_bounds
  USE mpi_distributed_memory                                 , ONLY: gather_to_all
  USE mesh_disc_apply_operators                              , ONLY: ddx_a_b_2D, ddy_a_b_2D, map_a_b_2D, map_b_a_2D
  USE laddie_utilities                                       , ONLY: compute_ambient_TS, map_H_a_b, map_H_a_c
  USE laddie_physics                                         , ONLY: compute_buoyancy
  USE CSR_matrix_vector_multiplication                       , ONLY: multiply_CSR_matrix_with_vector_1D
  USE mesh_halo_exchange                                     , ONLY: exchange_halos
  USE checksum_mod                                           , ONLY: checksum
  USE transects_main                                         
  USE string_module 
  USE transect_types
  USE mesh_utilities, ONLY:find_containing_vertex

  IMPLICIT NONE

  public 

CONTAINS

  subroutine initialise_transects_SGD( mesh, forcing)

    ! In/output variables
    type(type_mesh),                 intent(in   ) :: mesh
    type(type_laddie_forcing),       intent(inout) :: forcing

    ! Local variables:
    character(len=1024), parameter                 :: routine_name = 'initialise_transects_SGD'
    character(len=1024)                            :: transects_str
    character(len=1024), dimension(:), allocatable :: transect_strs
    integer                                        :: it

    ! Add routine to path
    call init_routine( routine_name)

    transects_str = C%transects_SGD

    ! If empty, finalise routine
    if (transects_str == '') then
      allocate( forcing%transects( 0))
      call finalise_routine( routine_name)
      return
    end if

    ! Separate by double vertical bars
    call separate_strings_by_double_vertical_bars( transects_str, transect_strs)

    ! Allocate the right size
    allocate( forcing%transects( size(transect_strs)))

    ! Initialise each transact
    do it = 1, size(forcing%transects)
      call initialise_transect_SGD( mesh, forcing%transects(it), transect_strs(it))
    end do

    print*, 'nr of channels=', size(forcing%transects)
    print*, 'name=', forcing%transects%name
    print*, 'SF=',   forcing%transects%flux_strength
    print*, 'nV=', forcing%transects%nV
    print*, 'nV_loc=', forcing%transects%nV_loc
    print*, 'vi1=', forcing%transects%vi1
    print*, 'vi2=', forcing%transects%vi2
    ! print*, 'index_point', forcing%transects%index_point(1)
    ! print*, 'V=', forcing%transects%V

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_transects_SGD


  subroutine initialise_transect_SGD( mesh, transect, transect_str)

    ! In/output variables
    type(type_mesh),     intent(in   ) :: mesh
    type(type_transect), intent(  out) :: transect
    character(len=*),    intent(in   ) :: transect_str

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'initialise_transect_SGD'
    character(len=1024)                   :: source, name, filename
    real(dp), dimension(:,:), allocatable :: waypoints
    real(dp), dimension(2) :: p
    real(dp)                              :: flux_strength
    integer                               :: index_point, vi

    ! Add routine to path
    call init_routine( routine_name)

    call parse_transect_str_SGD( transect_str, source, name, filename, flux_strength)

    transect%name            = name
    transect%flux_strength   = flux_strength

    if (par%primary) write(0,*) '  Initialising output transect ', &
      colour_string( trim( transect%name),'light blue'), '...'

    select case (source)
    case default
      call crash('invalid transect source "' // trim( source) // '"')
    case ('hardcoded')
      ! do nothing
    case ('read_from_file')
      call initialise_transect_waypoints_from_file( filename, waypoints)
    end select

    ! call calc_vertices_from_waypoints(transect, waypoints)
    
    ! First make sure that transect is connected by calculating transect vertices
    call calc_transect_vertices_from_waypoints( transect, waypoints, 1000._dp)

    index_point = 1
    ! Find vertex in mesh that contains vertex in transect
    do vi = transect%vi1, transect%vi2
       p = transect%V(vi,:)
       print*, 'p_tran=', p
       call find_containing_vertex( mesh, p, index_point)
       print*, 'p_mesh=', mesh%V(index_point, :)
       
       ! Save index_points of mesh to transect type
       ! transect%index_point(vi) = index_point
       
    end do

    ! p(1) = waypoints(1, 1)
    ! p(2) = waypoints(1, 2)

    ! index_point = 1
    ! print*, 'p=', p
    ! call find_containing_vertex( mesh, p, index_point)

    ! print*, index_point
    ! print*, mesh%V(index_point, :)

    ! transect%nz = mesh%nz
    ! allocate( transect%zeta( mesh%nz))
    ! transect%zeta = mesh%zeta

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_transect_SGD


  subroutine parse_transect_str_SGD( transect_str, source, name, filename, flux_strength)

    ! The transect should be specified by a name (either one of the hard-coded options,
    ! or a direction to an external file), and a resolution, e.g.:
    !
    !   positiveyaxis,dx=5e3
    !
    ! This indicates using the hard-coded option "positiveyaxis" with a resolution of 5 km
    !
    !   file:transect_MISMIP+_crossshelf.cfg,dx=2e3
    !

    ! In/output variables:
    character(len=*), intent(in   ) :: transect_str
    character(len=*), intent(  out) :: source, name, filename
    real(dp),         intent(  out) :: flux_strength

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'parse_transect_str_SGD'
    integer                        :: i

    ! Add routine to path
    call init_routine( routine_name)

    i = index( transect_str,',SF=')
    if (i==0) call crash('invalid transect string "' // trim(transect_str) // '" - could not find resolution specification!')

    name = transect_str( 1:i-1)
    read( transect_str( i+4:len_trim(transect_str)),*) flux_strength

    ! Separate source and name (and filename)
    i = index(name,'file:')
    if (i==0) then
      source = 'hardcoded'
    elseif (i==1) then
      source = 'read_from_file'
      filename = name(6:len_trim(name))
      i = index( filename,'/',back=.true.)
      name = filename(i+1:len_trim(filename)-4)
    else
      call crash('invalid transect string "' // trim(transect_str))
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine parse_transect_str_SGD

  ! SUBROUTINE compute_SGD_at_transects( mesh, laddie, forcing)

  ! ! Compute subglacial discharge (SGD) from intersect of transects with grounding line

  !   ! In- and output variables

  !   TYPE(type_mesh),                        INTENT(IN   ) :: mesh
  !   TYPE(type_laddie_model),                INTENT(INOUT) :: laddie
  !   TYPE(type_laddie_forcing),              INTENT(IN   ) :: forcing

  !   ! Local variables:
  !   CHARACTER(LEN=256), PARAMETER                         :: routine_name = 'compute_SGD_at_transects'
  !   INTEGER                                               :: vi, ierr
  !   REAL(dp)                                              :: total_area 

  !   ! Add routine to path
  !   CALL init_routine( routine_name)

  !   ! Initialise transects:
  !   CALL initialise_transects_SGD(mesh, forcing)

  !   ! Loop over transect
  !     ! Initialise waypoints from file 
  !       ! 
  !   ! DO transect in compute_SGD_at_transects
  !   ! CALL 


  !   ! If laddie-stand alone: only apply at initialisation
  !   ! If laddie not stand alone: apply at every ice dt


  !   ! Find intersect of transect_i with GL --> Should be one vertex
  !   ! Read in 
  !   ! Loop over all transects 

  !   ! Read in all transects 
  !   ! Need to get: forcing%mask_SGD( vi) / 
  !   ! or laddie%SGD( vi) = C%laddie_SGD_flux / total_area


  !   ! Finalise routine path
  !   CALL finalise_routine( routine_name)

  ! END SUBROUTINE compute_SGD_at_transects



END MODULE laddie_hydrology

