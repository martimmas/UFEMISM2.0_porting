module laddie_forcing_main

  ! Set up forcing of the standalone LADDIE model

! ===== Preamble =====
! ====================

  use parameters
  use precisions, only: dp
  use control_resources_and_error_messaging, only: crash, init_routine, finalise_routine, colour_string, warning
  use model_configuration, only: C
  use reference_geometry_types, only: type_reference_geometry
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use ocean_model_types, only: type_ocean_model
  use laddie_forcing_types, only: type_laddie_forcing
  use netcdf_io_main
  use reference_geometries_main, only: initialise_reference_geometry_raw_from_file, remap_reference_geometry_to_mesh, &
    reallocate_reference_geometry_on_mesh
  use laddie_utilities , only: allocate_laddie_forcing
  use mesh_creation_main, only: create_mesh_from_gridded_geometry, create_mesh_from_meshed_geometry, write_mesh_success
  use ice_model_memory, only: allocate_ice_model
  use thermodynamics_main, only: initialise_ice_temperature_uniform
  use masks_mod, only: determine_masks, calc_mask_ROI, calc_mask_noice, calc_mask_SGD
  use conservation_of_mass_main, only: apply_ice_thickness_BC_explicit, apply_mask_noice_direct
  use ice_geometry_basics, only: ice_surface_elevation, thickness_above_floatation, Hi_from_Hb_Hs_and_SL
  use subgrid_ice_margin, only: calc_effective_thickness
  use subgrid_grounded_fractions_main, only: calc_grounded_fractions
  use ice_shelf_base_slopes_onesided, only: calc_ice_shelf_base_slopes_onesided
  use ocean_main, only: initialise_ocean_model
  use projections, only: inverse_oblique_sg_projection

  implicit none

contains

! ===== Main routine =====
! ========================

  subroutine initialise_forcing( mesh, forcing)
    ! Set up all forcing

    ! In/output variables
    type(type_mesh)          ,  intent(out) :: mesh
    type(type_laddie_forcing),  intent(out) :: forcing

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'initialise_forcing'
    type(type_reference_geometry)  :: refgeo
    type(type_ice_model)           :: ice
    type(type_ocean_model)         :: ocean
    integer                        :: vi
    character(len=256), parameter  :: mesh_name = 'mesh'

    ! Add routine to path
    call init_routine( routine_name)

    ! Raw geometry
    ! ============

    ! Clean up memory if necessary
    if (allocated( refgeo%Hi_grid_raw)) deallocate( refgeo%Hi_grid_raw)
    if (allocated( refgeo%Hb_grid_raw)) deallocate( refgeo%Hb_grid_raw)
    if (allocated( refgeo%Hs_grid_raw)) deallocate( refgeo%Hs_grid_raw)
    if (allocated( refgeo%SL_grid_raw)) deallocate( refgeo%SL_grid_raw)

    call initialise_reference_geometry_raw_from_file( 'ANT', 'refgeo', refgeo, C%filename_refgeo_PD_ANT, C%timeframe_refgeo_PD_ANT)

    ! Setup mesh
    ! ==========

    ! Determine if the initial geometry is provided gridded or meshed
    if (allocated( refgeo%grid_raw%x)) then 
      ! Gridded

      ! Safety
      if (allocated( refgeo%mesh_raw%V)) call crash('found both grid and mesh in refgeo!')

      ! Create mesh from gridded initial geometry data
      call create_mesh_from_gridded_geometry( 'ANT', mesh_name, &
        refgeo%grid_raw, &
        refgeo%Hi_grid_raw, &
        refgeo%Hb_grid_raw, &
        refgeo%Hs_grid_raw, &
        refgeo%SL_grid_raw, &
        C%xmin_ANT, C%xmax_ANT, C%ymin_ANT, C%ymax_ANT, C%lambda_M_ANT, C%phi_M_ANT, C%beta_stereo_ANT, &
        mesh)

    elseif (allocated( refgeo%mesh_raw%V)) then 
      ! Meshed

      ! Safety
      if (allocated( refgeo%grid_raw%x)) call crash('found both grid and mesh in refgeo!')

      ! Create mesh from meshed initial geometry data
      call create_mesh_from_meshed_geometry( 'ANT', mesh_name, &
        refgeo%mesh_raw, &
        refgeo%Hi_mesh_raw, &
        refgeo%Hb_mesh_raw, &
        refgeo%Hs_mesh_raw, &
        refgeo%SL_mesh_raw, &
        C%xmin_ANT, C%xmax_ANT, C%ymin_ANT, C%ymax_ANT, C%lambda_M_ANT, C%phi_M_ANT, C%beta_stereo_ANT, &
        mesh)

    else 
      call crash('no grid or mesh is found in refgeo!')
    end if

    ! Write the mesh creation success message to the terminal
    call write_mesh_success( mesh)

    ! Geometry on mesh
    ! ================

    call reallocate_reference_geometry_on_mesh( mesh, refgeo)
    call remap_reference_geometry_to_mesh( mesh, refgeo)

    ! Use ice model
    ! =============

    call allocate_ice_model( mesh, ice)

    ! Basic geometry
    do vi = mesh%vi1, mesh%vi2
      ice%Hb( vi) = refgeo%Hb( vi)
      ice%Hs( vi) = refgeo%Hs ( vi)
      ice%Hi( vi) = Hi_from_Hb_Hs_and_SL( ice%Hb( vi), ice%Hs( vi), ice%SL( vi)) 
    end do

    ! Calculate the no-ice mask
    call calc_mask_noice( mesh, ice) 

    ! Apply no-ice mask
    call apply_mask_noice_direct( mesh, ice%mask_noice, ice%Hi)

    ! Apply boundary conditions at the domain border
    call apply_ice_thickness_BC_explicit( mesh, ice%mask_noice, ice%Hb, ice%SL, ice%Hi)

    do vi = mesh%vi1, mesh%vi2

      ! Derived geometry
      ice%Hs ( vi) = ice_surface_elevation( ice%Hi( vi), ice%Hb( vi), ice%SL( vi)) 
      ice%Hib( vi) = ice%Hs( vi) - ice%Hi( vi)
      ice%TAF( vi) = thickness_above_floatation( ice%Hi( vi), ice%Hb( vi), ice%SL( vi)) 

    end do

    ! Initialise masks
    ! ================

    ! call it twice so also the "prev" versions are set
    call determine_masks( mesh, ice)
    call determine_masks( mesh, ice)

    ! Compute mask_ROI only at initialisation, (NOTE: This works only for one single ROI right now)
    call calc_mask_ROI( mesh, ice, 'ANT')

    ! Compute mask_SGD only at initialisation
    call calc_mask_SGD( mesh, ice)

    ! Effective ice thickness
    ! =======================

    ! Compute effective thickness at calving fronts
    call calc_effective_thickness( mesh, ice, ice%Hi, ice%Hi_eff, ice%fraction_margin)

    ! Calculate ice shelf draft gradients
    call calc_ice_shelf_base_slopes_onesided( mesh, ice)

    ! Ice temperature
    ! ===============

    ! Only uniform for now
    call initialise_ice_temperature_uniform( mesh, ice, 'ANT')


    ! Ocean forcing
    ! =============

    call initialise_ocean_model( mesh, ice, ocean, 'ANT', C%start_time_of_run)

    ! Allocate forcing
    ! ================

    call allocate_laddie_forcing( mesh, forcing)

    ! Write all to laddie forcing
    ! ===========================

    forcing%Hi                ( mesh%vi1:mesh%vi2  ) = ice%Hi                ( mesh%vi1:mesh%vi2  )
    forcing%Hib               ( mesh%vi1:mesh%vi2  ) = ice%Hib               ( mesh%vi1:mesh%vi2  )
    forcing%TAF               ( mesh%vi1:mesh%vi2  ) = ice%TAF               ( mesh%vi1:mesh%vi2  )
    forcing%dHib_dx_b         ( mesh%ti1:mesh%ti2  ) = ice%dHib_dx_b         ( mesh%ti1:mesh%ti2  )
    forcing%dHib_dy_b         ( mesh%ti1:mesh%ti2  ) = ice%dHib_dy_b         ( mesh%ti1:mesh%ti2  )
    forcing%mask_icefree_land ( mesh%vi1:mesh%vi2  ) = ice%mask_icefree_land ( mesh%vi1:mesh%vi2  )
    forcing%mask_icefree_ocean( mesh%vi1:mesh%vi2  ) = ice%mask_icefree_ocean( mesh%vi1:mesh%vi2  )
    forcing%mask_grounded_ice ( mesh%vi1:mesh%vi2  ) = ice%mask_grounded_ice ( mesh%vi1:mesh%vi2  )
    forcing%mask_floating_ice ( mesh%vi1:mesh%vi2  ) = ice%mask_floating_ice ( mesh%vi1:mesh%vi2  )

    forcing%mask_gl_fl        ( mesh%vi1:mesh%vi2  ) = ice%mask_gl_fl        ( mesh%vi1:mesh%vi2  )
    forcing%mask_SGD          ( mesh%vi1:mesh%vi2  ) = ice%mask_SGD          ( mesh%vi1:mesh%vi2  )
    forcing%mask              ( mesh%vi1:mesh%vi2  ) = ice%mask              ( mesh%vi1:mesh%vi2  )
        
    forcing%Ti                ( mesh%vi1:mesh%vi2,:) = ice%Ti                ( mesh%vi1:mesh%vi2,:) - 273.15 ! [degC]
    forcing%T_ocean           ( mesh%vi1:mesh%vi2,:) = ocean%T               ( mesh%vi1:mesh%vi2,:)
    forcing%S_ocean           ( mesh%vi1:mesh%vi2,:) = ocean%S               ( mesh%vi1:mesh%vi2,:)

    ! Determine coriolis parameter
    ! ============================

    call calculate_coriolis_parameter( mesh, forcing, C%lambda_M_ANT, C%phi_M_ANT, C%beta_stereo_ANT)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine initialise_forcing

  subroutine calculate_coriolis_parameter( mesh, forcing, lambda_M, phi_M, beta_stereo)
    ! Set up all forcing

    ! In/output variables
    type(type_mesh)          , intent(in   ) :: mesh
    type(type_laddie_forcing), intent(inout) :: forcing
    real(dp)                 , intent(in   ) :: lambda_M, phi_M, beta_stereo

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'calculate_coriolis_parameter'
    real(dp)                       :: lon_b, lat_b
    integer                        :: ti

    ! Add routine to path
    call init_routine( routine_name)

    select case (C%choice_laddie_coriolis)
      case default
        call crash('unknown choice_laddie_coriolis "' // trim( C%choice_laddie_coriolis) // '"!')
      case ('uniform')
        forcing%f_coriolis( mesh%ti1:mesh%ti2) = C%uniform_laddie_coriolis_parameter
      case ('realistic')
        do ti = mesh%ti1, mesh%ti2
          call inverse_oblique_sg_projection( mesh%Tricc( ti,1), mesh%Tricc( ti,2), lambda_M, phi_M, beta_stereo, lon_b, lat_b)
          forcing%f_coriolis( ti) = 2.0_dp * earth_rotation_rate * sin(2.0_dp * pi * lat_b / 360.0_dp)
        end do

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calculate_coriolis_parameter

end module laddie_forcing_main
