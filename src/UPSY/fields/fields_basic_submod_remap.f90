submodule (fields_basic) fields_basic_submod_remap

contains

    subroutine remap_dp_2D( self, mesh_new, d_nih)

      ! In/output variables:
      class(atype_field),                          intent(inout) :: self
      type(type_mesh),                             intent(in   ) :: mesh_new
      real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'atype_field_remap_dp_2D'

      ! Add routine to call stack
      call init_routine( routine_name)

      ! Downcast field and grid
      select type (f => self)
      class default
        call crash('invalid field type')
      class is (type_field_dp_2D)

        select type (g => self%grid())
        class default
          call crash('invalid field%grid type')
        class is (type_mesh)
          call remap_dp_2D_applied( f, g, mesh_new, d_nih)
        end select

      end select

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine remap_dp_2D

    subroutine remap_dp_2D_applied( field, mesh_old, mesh_new, d_nih)

      ! In/output variables:
      type(type_field_dp_2D),                      intent(inout) :: field
      type(type_mesh),                             intent(in   ) :: mesh_old
      type(type_mesh),                             intent(in   ) :: mesh_new
      real(dp), dimension(:), contiguous, pointer, intent(inout) :: d_nih

      ! Local variables:
      character(len=1024), parameter      :: routine_name = 'remap_dp_2D_applied'
      real(dp), dimension(:), allocatable :: d_dist

      ! Add routine to call stack
      call init_routine( routine_name)

      ! Copy actual array to distributed memory
      allocate( d_dist( mesh_old%vi1: mesh_old%vi2), &
        source = d_nih( mesh_old%vi1: mesh_old%vi2))
      call hybrid_to_dist( mesh_old%pai_V, d_nih, d_dist)

      ! Reallocate hybrid memory and update bounds and pointers
      call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_V%n_nih)
      d_nih      ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih) => d_nih
      field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih) => d_nih

      ! Remap the distributed array (remapping code doesn't yet work on hybrid memory)
      call map_from_mesh_to_mesh_with_reallocation_2D( mesh_old, mesh_new, &
        C%output_dir, d_dist, field%remap_method())

      ! Copy the result back to hybrid memory
      call dist_to_hybrid( mesh_new%pai_V, d_dist, d_nih)

      ! Update field grid and parallel array info
      call field%set_grid( mesh_new)
      call field%set_pai ( mesh_new%pai_V)

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine remap_dp_2D_applied

    subroutine remap_dp_3D( self, mesh_new, d_nih)

      ! In/output variables:
      class(atype_field),                            intent(inout) :: self
      type(type_mesh),                               intent(in   ) :: mesh_new
      real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

      ! Local variables:
      character(len=1024), parameter :: routine_name = 'atype_field_remap_dp_3D'

      ! Add routine to call stack
      call init_routine( routine_name)

      ! Downcast field and grid
      select type (f => self)
      class default
        call crash('invalid field type')
      class is (type_field_dp_3D)

        select type (g => self%grid())
        class default
          call crash('invalid field%grid type')
        class is (type_mesh)
          call remap_dp_3D_applied( f, g, mesh_new, d_nih)
        end select

      end select

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine remap_dp_3D

    subroutine remap_dp_3D_applied( field, mesh_old, mesh_new, d_nih)

      ! In/output variables:
      type(type_field_dp_3D),                        intent(inout) :: field
      type(type_mesh),                               intent(in   ) :: mesh_old
      type(type_mesh),                               intent(in   ) :: mesh_new
      real(dp), dimension(:,:), contiguous, pointer, intent(inout) :: d_nih

      ! Local variables:
      character(len=1024), parameter        :: routine_name = 'remap_dp_3D_applied'
      integer                               :: nz
      real(dp), dimension(:,:), allocatable :: d_dist

      ! Add routine to call stack
      call init_routine( routine_name)

      nz = size( d_nih,2)

      ! Copy actual array to distributed memory
      allocate( d_dist( mesh_old%vi1: mesh_old%vi2, 1:nz), &
        source = d_nih( mesh_old%vi1: mesh_old%vi2, 1:nz))
      call hybrid_to_dist( mesh_old%pai_V, nz, d_nih, d_dist)

      ! Reallocate hybrid memory and update bounds and pointers
      call reallocate_dist_shared( d_nih, field%w, mesh_new%pai_V%n_nih, nz)
      d_nih      ( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:nz) => d_nih
      field%d_nih( mesh_new%pai_V%i1_nih: mesh_new%pai_V%i2_nih, 1:nz) => d_nih

      ! Remap the distributed array (remapping code doesn't yet work on hybrid memory)
      call map_from_mesh_to_mesh_with_reallocation_3D( mesh_old, mesh_new, &
        C%output_dir, d_dist, field%remap_method())

      ! Copy the result back to hybrid memory
      call dist_to_hybrid( mesh_new%pai_V, nz, d_dist, d_nih)

      ! Update field grid and parallel array info
      call field%set_grid( mesh_new)
      call field%set_pai ( mesh_new%pai_V)

      ! Remove routine from call stack
      call finalise_routine( routine_name)

    end subroutine remap_dp_3D_applied

end submodule fields_basic_submod_remap