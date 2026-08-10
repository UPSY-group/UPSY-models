module apply_maps

  ! Apply pre-created mapping operators to data fields to remap
  ! data fields between different grids/meshes.

  use petscksp
  use precisions, only: dp
  use mpi_basic, only: par, sync
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
  use mesh_types, only: type_mesh
  use grid_types, only: type_grid, type_grid_lonlat
  use remapping_types, only: type_map
  use CSR_matrix_mod, only: type_CSR_matrix_dp
  use petsc_basic, only: multiply_PETSc_matrix_with_vector_1D, multiply_PETSc_matrix_with_vector_2D, &
    mat_petsc2CSR
  use mesh_utilities, only: set_border_vertices_to_interior_mean_dp_2D, set_border_vertices_to_interior_mean_dp_3D, &
    set_border_triangles_to_interior_mean_dp_2D, set_border_triangles_to_interior_mean_dp_3D
  use mpi_distributed_memory, only: gather_to_all
  use mpi_distributed_memory_grid, only: gather_gridded_data_to_primary, distribute_gridded_data_from_primary
  use CSR_matrix_vector_multiplication, only: multiply_CSR_matrix_with_vector_1D_wrapper, &
    multiply_CSR_matrix_with_vector_2D_wrapper
  use dist_to_hybrid_mod, only: dist_to_hybrid, hybrid_to_dist

  implicit none

  private

  public :: Atlas, clear_all_maps_involving_this_mesh, &
    apply_map_xy_grid_to_mesh_2D, apply_map_xy_grid_to_mesh_3D, &
    apply_map_xy_grid_to_mesh_triangles_2D, apply_map_xy_grid_to_mesh_triangles_3D, &
    apply_map_lonlat_grid_to_mesh_2D, apply_map_lonlat_grid_to_mesh_3D, &
    apply_map_mesh_vertices_to_xy_grid_2D, apply_map_mesh_vertices_to_xy_grid_3D, apply_map_mesh_vertices_to_xy_grid_2D_minval, &
    apply_map_mesh_triangles_to_xy_grid_2D, apply_map_mesh_triangles_to_xy_grid_3D, &
    apply_map_mesh_to_mesh_2D, apply_map_mesh_to_mesh_3D, &
    apply_map_mesh_tri_to_mesh_tri_2D, apply_map_mesh_tri_to_mesh_tri_3D

  ! The Atlas: the complete collection of all mapping objects.
  type(type_map), dimension(1000) :: Atlas

contains

  !> Clear all mapping objects involving a mesh by this name from the Atlas
  !> (used after a mesh update, as the mapping objects for the old mesh are useless anyway)
  subroutine clear_all_maps_involving_this_mesh( mesh)

    ! In/output variables
    type(type_mesh), intent(in)    :: mesh

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'map_from_xy_grid_to_mesh_2D'
    integer                        :: mi, perr

    ! Add routine to path
    call init_routine( routine_name)

    do mi = 1, size( Atlas,1)
      if (Atlas( mi)%is_in_use .and. &
        (Atlas( mi)%name_src == mesh%name .or. &
         Atlas( mi)%name_dst == mesh%name .or. &
         Atlas( mi)%name_src == (trim( mesh%name) // '_triangles') .or. &
         Atlas( mi)%name_dst == (trim( mesh%name) // '_triangles'))) then
        ! This map involves the current mesh
        Atlas( mi)%is_in_use = .false.
        Atlas( mi)%name_src  = ''
        Atlas( mi)%name_dst  = ''
        call MatDestroy( Atlas( mi)%M, perr)
      end if
    end do

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine clear_all_maps_involving_this_mesh

  ! ===== x/y-grid to mesh vertices =====
  ! =====================================

  !> Map a 2-D data field from an x/y-grid to a mesh.
  subroutine apply_map_xy_grid_to_mesh_2D( grid, mesh, map, d_grid_vec_partial, d_mesh)

    ! In/output variables
    type(type_grid),        intent(in   ) :: grid
    type(type_mesh),        intent(in   ) :: mesh
    type(type_map),         intent(in   ) :: map
    real(dp), dimension(:), intent(in   ) :: d_grid_vec_partial
    real(dp), dimension(:), intent(  out) :: d_mesh

    ! Local variables:
    character(len=*), parameter            :: routine_name = 'apply_map_xy_grid_to_mesh_2D'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: d_mesh_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_1D( map%M, d_grid_vec_partial, d_mesh_loc)

    ! Support d_mesh as both distributed and hybrid distributed/shared memory
    if (size( d_mesh,1) == mesh%pai_V%n_loc) then
      d_mesh = d_mesh_loc
    elseif (size( d_mesh,1) == mesh%pai_V%n_nih) then
      call dist_to_hybrid( mesh%pai_V, d_mesh_loc, d_mesh)
    else
      call crash('invalid size for d_mesh')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_xy_grid_to_mesh_2D

  !> Map a 3-D data field from an x/y-grid to a mesh.
  subroutine apply_map_xy_grid_to_mesh_3D( grid, mesh, map, d_grid_vec_partial, d_mesh)

    ! In/output variables
    type(type_grid),          intent(in   ) :: grid
    type(type_mesh),          intent(in   ) :: mesh
    type(type_map),           intent(in   ) :: map
    real(dp), dimension(:,:), intent(in   ) :: d_grid_vec_partial
    real(dp), dimension(:,:), intent(  out) :: d_mesh

    ! Local variables:
    character(len=*), parameter                             :: routine_name = 'apply_map_xy_grid_to_mesh_3D'
    real(dp), dimension(mesh%vi1:mesh%vi2,1:size(d_mesh,2)) :: d_mesh_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_2D( map%M, d_grid_vec_partial, d_mesh_loc)

    ! Support d_mesh as both distributed and hybrid distributed/shared memory
    if (size( d_mesh,1) == mesh%pai_V%n_loc) then
      d_mesh = d_mesh_loc
    elseif (size( d_mesh,1) == mesh%pai_V%n_nih) then
      call dist_to_hybrid( mesh%pai_V, size( d_mesh,2), d_mesh_loc, d_mesh)
    else
      call crash('invalid size for d_mesh')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_xy_grid_to_mesh_3D

  ! ===== x/y-grid to mesh triangles =====
  ! ======================================

  ! Map a 2-D data field from an x/y-grid to a mesh triangles.
  subroutine apply_map_xy_grid_to_mesh_triangles_2D( grid, mesh, map, d_grid_vec_partial, d_mesh)

    ! In/output variables
    type(type_grid),                intent(in)  :: grid
    type(type_mesh),                intent(in)  :: mesh
    type(type_map),                 intent(in)  :: map
    real(dp), dimension(:),         intent(in)  :: d_grid_vec_partial
    real(dp), dimension(:), target, intent(out) :: d_mesh

    ! Local variables:
    character(len=*), parameter            :: routine_name = 'apply_map_xy_grid_to_mesh_triangles_2D'
    real(dp), dimension(mesh%ti1:mesh%ti2) :: d_mesh_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_1D( map%M, d_grid_vec_partial, d_mesh_loc)

    ! Support d_mesh as both distributed and hybrid distributed/shared memory
    if (size( d_mesh,1) == mesh%pai_Tri%n_loc) then
      d_mesh = d_mesh_loc
    elseif (size( d_mesh,1) == mesh%pai_Tri%n_nih) then
      call dist_to_hybrid( mesh%pai_Tri, d_mesh_loc, d_mesh)
    else
      call crash('invalid size for d_mesh')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_xy_grid_to_mesh_triangles_2D

  !> Map a 3-D data field from an x/y-grid to a mesh triangles.
  subroutine apply_map_xy_grid_to_mesh_triangles_3D( grid, mesh, map, d_grid_vec_partial, d_mesh)

    ! In/output variables
    type(type_grid),                  intent(in)  :: grid
    type(type_mesh),                  intent(in)  :: mesh
    type(type_map),                   intent(in)  :: map
    real(dp), dimension(:,:),         intent(in)  :: d_grid_vec_partial
    real(dp), dimension(:,:), target, intent(out) :: d_mesh

    ! Local variables:
    character(len=*), parameter                             :: routine_name = 'apply_map_xy_grid_to_mesh_triangles_3D'
    real(dp), dimension(mesh%ti1:mesh%ti2,1:size(d_mesh,2)) :: d_mesh_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_2D( map%M, d_grid_vec_partial, d_mesh_loc)

    ! Support d_mesh as both distributed and hybrid distributed/shared memory
    if (size( d_mesh,1) == mesh%pai_Tri%n_loc) then
      d_mesh = d_mesh_loc
    elseif (size( d_mesh,1) == mesh%pai_Tri%n_nih) then
      call dist_to_hybrid( mesh%pai_Tri, size( d_mesh,2), d_mesh_loc, d_mesh)
    else
      call crash('invalid size for d_mesh')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_xy_grid_to_mesh_triangles_3D

  ! ===== lon/lat-grid to mesh =====
  ! ================================

  !> Map a 2-D data field from a lon/lat-grid to a mesh.
  subroutine apply_map_lonlat_grid_to_mesh_2D( grid, mesh, map, d_grid_vec_partial, d_mesh)

    ! In/output variables
    type(type_grid_lonlat),         intent(in)  :: grid
    type(type_mesh),                intent(in)  :: mesh
    type(type_map),                 intent(in)  :: map
    real(dp), dimension(:),         intent(in)  :: d_grid_vec_partial
    real(dp), dimension(:), target, intent(out) :: d_mesh

    ! Local variables:
    character(len=*), parameter            :: routine_name = 'apply_map_lonlat_grid_to_mesh_2D'
    real(dp), dimension(mesh%vi1:mesh%vi2) :: d_mesh_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_1D( map%M, d_grid_vec_partial, d_mesh_loc)

    ! Support d_mesh as both distributed and hybrid distributed/shared memory
    if (size( d_mesh,1) == mesh%pai_V%n_loc) then
      d_mesh = d_mesh_loc
    elseif (size( d_mesh,1) == mesh%pai_V%n_nih) then
      call dist_to_hybrid( mesh%pai_V, d_mesh_loc, d_mesh)
    else
      call crash('invalid size for d_mesh')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_lonlat_grid_to_mesh_2D

  !> Map a 3-D data field from a lon/lat-grid to a mesh.
  subroutine apply_map_lonlat_grid_to_mesh_3D( grid, mesh, map, d_grid_vec_partial, d_mesh)

    ! In/output variables
    type(type_grid_lonlat),           intent(in)  :: grid
    type(type_mesh),                  intent(in)  :: mesh
    type(type_map),                   intent(in)  :: map
    real(dp), dimension(:,:),         intent(in)  :: d_grid_vec_partial
    real(dp), dimension(:,:), target, intent(out) :: d_mesh

    ! Local variables:
    character(len=*), parameter                             :: routine_name = 'apply_map_lonlat_grid_to_mesh_3D'
    real(dp), dimension(mesh%vi1:mesh%vi2,1:size(d_mesh,2)) :: d_mesh_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_2D( map%M, d_grid_vec_partial, d_mesh_loc)

    ! Support d_mesh as both distributed and hybrid distributed/shared memory
    if (size( d_mesh,1) == mesh%pai_V%n_loc) then
      d_mesh = d_mesh_loc
    elseif (size( d_mesh,1) == mesh%pai_V%n_nih) then
      call dist_to_hybrid( mesh%pai_V, size( d_mesh,2), d_mesh_loc, d_mesh)
    else
      call crash('invalid size for d_mesh')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_lonlat_grid_to_mesh_3D

  ! ===== mesh to x/y-grid =====
  ! ============================

  subroutine apply_map_mesh_vertices_to_xy_grid_2D( mesh, grid, map, d_mesh_partial, d_grid_vec_partial)
    !< Map a 2-D data field from the vertices of a mesh to an x/y-grid.

    ! In/output variables
    type(type_mesh),        intent(in   )  :: mesh
    type(type_grid),        intent(in   )  :: grid
    type(type_map),         intent(in   )  :: map
    real(dp), dimension(:), intent(in   )  :: d_mesh_partial
    real(dp), dimension(:), intent(  out) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'apply_map_mesh_vertices_to_xy_grid_2D'
    type(type_CSR_matrix_dp)         :: M_CSR
    real(dp), dimension(:,:  ), allocatable :: d_grid

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call mat_petsc2CSR( map%M, M_CSR)
    call multiply_CSR_matrix_with_vector_1D_wrapper( M_CSR, &
      mesh%pai_V, d_mesh_partial, grid%pai, d_grid_vec_partial, &
      buffer_xx_nih = mesh%buffer1_d_a_nih)

    ! == Because the remapping operators are sometimes inaccurate at the
    !     domain boundary, set values in the outermost row of grid cells
    !    equal to those in the second-outermost row

    if (par%primary) then
      ! allocate memory for complete gridded data
      allocate( d_grid( grid%nx, grid%ny))
      ! Gather complete gridded data
      call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)
      ! Set values in the outermost row of grid cells
      ! equal to those in the second-outermost row
      d_grid( 1      ,:) = d_grid( 2        ,:)
      d_grid( grid%nx,:) = d_grid( grid%nx-1,:)
      d_grid( :,1      ) = d_grid( :,2        )
      d_grid( :,grid%ny) = d_grid( :,grid%ny-1)
      ! Distribute complete gridded data back over the processes
      call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid)
      ! Clean up after yourself
      deallocate( d_grid)
    else ! if (par%primary) then
      ! allocate zero memory for complete gridded data (only the primary needs this)
      allocate( d_grid( 0,0))
      ! Gather complete gridded data
      call gather_gridded_data_to_primary( grid, d_grid_vec_partial)
      ! Distribute complete gridded data back over the processes
      call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid)
      ! Clean up after yourself
      deallocate( d_grid)
    end if ! if (par%primary) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_vertices_to_xy_grid_2D

  subroutine apply_map_mesh_vertices_to_xy_grid_3D( mesh, grid, map, d_mesh_partial, d_grid_vec_partial)
    !< Map a 3-D data field from the vertices of a mesh to an x/y-grid.

    ! In/output variables
    type(type_mesh),          intent(in   ) :: mesh
    type(type_grid),          intent(in   ) :: grid
    type(type_map),           intent(in   ) :: map
    real(dp), dimension(:,:), intent(in   ) :: d_mesh_partial
    real(dp), dimension(:,:), intent(  out) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'apply_map_mesh_vertices_to_xy_grid_3D'
    type(type_CSR_matrix_dp)         :: M_CSR
    real(dp), dimension(:,:,:), allocatable :: d_grid

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call mat_petsc2CSR( map%M, M_CSR)
    call multiply_CSR_matrix_with_vector_2D_wrapper( M_CSR, &
      mesh%pai_V, d_mesh_partial, grid%pai, d_grid_vec_partial, &
      buffer_xx_nih = mesh%buffer1_d_ak_nih)

    ! == Because the remapping operators are sometimes inaccurate at the
    !     domain boundary, set values in the outermost row of grid cells
    !    equal to those in the second-outermost row

    if (par%primary) then
      ! allocate memory for complete gridded data
      allocate( d_grid( grid%nx, grid%ny, size( d_mesh_partial,2)))
      ! Gather complete gridded data
      call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)
      ! Set values in the outermost row of grid cells
      ! equal to those in the second-outermost row
      d_grid( 1      ,:,:) = d_grid( 2        ,:,:)
      d_grid( grid%nx,:,:) = d_grid( grid%nx-1,:,:)
      d_grid( :,1      ,:) = d_grid( :,2        ,:)
      d_grid( :,grid%ny,:) = d_grid( :,grid%ny-1,:)
      ! Distribute complete gridded data back over the processes
      call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid)
      ! Clean up after yourself
      deallocate( d_grid)
    else ! if (par%primary) then
      ! allocate zero memory for complete gridded data (only the primary needs this)
      allocate( d_grid( 0,0,0))
      ! Gather complete gridded data
      call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)
      ! Distribute complete gridded data back over the processes
      call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid)
      ! Clean up after yourself
      deallocate( d_grid)
    end if ! if (par%primary) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_vertices_to_xy_grid_3D

  subroutine apply_map_mesh_vertices_to_xy_grid_2D_minval( mesh, grid, map, d_mesh_partial, d_grid_vec_partial)
    !< Map a 2-D data field from the vertices of a mesh to an x/y-grid.

    ! For each grid cell, get the minimum value of all overlapping mesh vertices

    ! In/output variables
    type(type_mesh),                        intent(in)  :: mesh
    type(type_grid),                        intent(in)  :: grid
    type(type_map),                         intent(in)  :: map
    real(dp), dimension(mesh%vi1:mesh%vi2), intent(in)  :: d_mesh_partial
    real(dp), dimension(grid%n1 :grid%n2 ), intent(out) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter  :: routine_name = 'apply_map_mesh_vertices_to_xy_grid_2D_minval'
    real(dp), dimension(mesh%nV)    :: d_mesh_tot
    type(type_CSR_matrix_dp) :: M_CSR
    integer                         :: n,k1,k2,k,col,vi
    real(dp)                        :: d_min

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather global mesh data
    call gather_to_all( d_mesh_partial, d_mesh_tot)

    ! Convert mapping matrix from PETSc format to UFEMISM CSR format
    call mat_petsc2CSR( map%M, M_CSR)

    ! Map data
    do n = grid%n1, grid%n2

      d_min = huge( d_min)

      ! Loop over all mesh vertices that this grid cell overlaps with
      k1 = M_CSR%ptr( n)
      k2 = M_CSR%ptr( n+1)-1
      do k = k1, k2
        col = M_CSR%ind( k)
        ! This matrix row corresponds to this mesh vertex
        vi = mesh%n2vi( col)
        ! Update minimum value
        d_min = min( d_min, d_mesh_tot( vi))
      end do

      ! Fill into array
      d_grid_vec_partial( n) = d_min

    end do ! do n = grid%n1, grid%n2

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_vertices_to_xy_grid_2D_minval

  subroutine apply_map_mesh_triangles_to_xy_grid_2D( mesh, grid, map, d_mesh_partial, d_grid_vec_partial)
    !< Map a 2-D data field from the triangles of a mesh to an x/y-grid.

    ! In/output variables
    type(type_mesh),        intent(in   ) :: mesh
    type(type_grid),        intent(in   ) :: grid
    type(type_map),         intent(in   ) :: map
    real(dp), dimension(:), intent(in   ) :: d_mesh_partial
    real(dp), dimension(:), intent(  out) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'apply_map_mesh_triangles_to_xy_grid_2D'
    type(type_CSR_matrix_dp)       :: M_CSR
    real(dp), dimension(:,:), allocatable :: d_grid

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call mat_petsc2CSR( map%M, M_CSR)
    call multiply_CSR_matrix_with_vector_1D_wrapper( M_CSR, &
      mesh%pai_Tri, d_mesh_partial, grid%pai, d_grid_vec_partial, &
      buffer_xx_nih = mesh%buffer1_d_b_nih)

    ! == Because the remapping operators are sometimes inaccurate at the
    !     domain boundary, set values in the outermost row of grid cells
    !    equal to those in the second-outermost row

    if (par%primary) then
      ! allocate memory for complete gridded data
      allocate( d_grid( grid%nx, grid%ny))
      ! Gather complete gridded data
      call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)
      ! Set values in the outermost row of grid cells
      ! equal to those in the second-outermost row
      d_grid( 1      ,:) = d_grid( 2        ,:)
      d_grid( grid%nx,:) = d_grid( grid%nx-1,:)
      d_grid( :,1      ) = d_grid( :,2        )
      d_grid( :,grid%ny) = d_grid( :,grid%ny-1)
      ! Distribute complete gridded data back over the processes
      call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid)
      ! Clean up after yourself
      deallocate( d_grid)
    else ! if (par%primary) then
      ! allocate zero memory for complete gridded data (only the primary needs this)
      allocate( d_grid( 0,0))
      ! Gather complete gridded data
      call gather_gridded_data_to_primary( grid, d_grid_vec_partial)
      ! Distribute complete gridded data back over the processes
      call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid)
      ! Clean up after yourself
      deallocate( d_grid)
    end if ! if (par%primary) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_triangles_to_xy_grid_2D

  subroutine apply_map_mesh_triangles_to_xy_grid_3D( mesh, grid, map, d_mesh_partial, d_grid_vec_partial)
    !< Map a 3-D data field from the triangles of a mesh to an x/y-grid.

    ! In/output variables
    type(type_mesh),          intent(in   ) :: mesh
    type(type_grid),          intent(in   ) :: grid
    type(type_map),           intent(in   ) :: map
    real(dp), dimension(:,:), intent(in   ) :: d_mesh_partial
    real(dp), dimension(:,:), intent(  out) :: d_grid_vec_partial

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'apply_map_mesh_triangles_to_xy_grid_3D'
    type(type_CSR_matrix_dp)         :: M_CSR
    real(dp), dimension(:,:,:), allocatable :: d_grid

    ! Add routine to path
    call init_routine( routine_name)

    ! Perform the mapping operation as a matrix multiplication
    call mat_petsc2CSR( map%M, M_CSR)
    call multiply_CSR_matrix_with_vector_2D_wrapper( M_CSR, &
      mesh%pai_Tri, d_mesh_partial, grid%pai, d_grid_vec_partial, &
      buffer_xx_nih = mesh%buffer1_d_bk_nih)

    ! == Because the remapping operators are sometimes inaccurate at the
    !     domain boundary, set values in the outermost row of grid cells
    !    equal to those in the second-outermost row

    if (par%primary) then
      ! allocate memory for complete gridded data
      allocate( d_grid( grid%nx, grid%ny, size( d_mesh_partial,2)))
      ! Gather complete gridded data
      call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)
      ! Set values in the outermost row of grid cells
      ! equal to those in the second-outermost row
      d_grid( 1      ,:,:) = d_grid( 2        ,:,:)
      d_grid( grid%nx,:,:) = d_grid( grid%nx-1,:,:)
      d_grid( :,1      ,:) = d_grid( :,2        ,:)
      d_grid( :,grid%ny,:) = d_grid( :,grid%ny-1,:)
      ! Distribute complete gridded data back over the processes
      call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid)
      ! Clean up after yourself
      deallocate( d_grid)
    else ! if (par%primary) then
      ! allocate zero memory for complete gridded data (only the primary needs this)
      allocate( d_grid( 0,0,0))
      ! Gather complete gridded data
      call gather_gridded_data_to_primary( grid, d_grid_vec_partial, d_grid)
      ! Distribute complete gridded data back over the processes
      call distribute_gridded_data_from_primary( grid, d_grid_vec_partial, d_grid)
      ! Clean up after yourself
      deallocate( d_grid)
    end if ! if (par%primary) then

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_triangles_to_xy_grid_3D

  ! ===== mesh to mesh =====
  ! ========================

  subroutine apply_map_mesh_to_mesh_2D( mesh_src, mesh_dst, map, d_src, d_dst)

    ! In/output variables
    type(type_mesh),        intent(in)  :: mesh_src
    type(type_mesh),        intent(in)  :: mesh_dst
    type(type_map),         intent(in)  :: map
    real(dp), dimension(:), intent(in)  :: d_src
    real(dp), dimension(:), intent(out) :: d_dst

    ! Local variables:
    character(len=*), parameter                    :: routine_name = 'apply_map_mesh_to_mesh_2D'
    real(dp), dimension(mesh_src%vi1:mesh_src%vi2) :: d_src_loc
    real(dp), dimension(mesh_dst%vi1:mesh_dst%vi2) :: d_dst_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Support d_src as both distributed and hybrid distributed/shared memory
    if (size( d_src,1) == mesh_src%pai_V%n_loc) then
      d_src_loc = d_src
    elseif (size( d_src,1) == mesh_src%pai_V%n_nih) then
      call hybrid_to_dist( mesh_src%pai_V, d_src, d_src_loc)
    else
      call crash('invalid size for d_src')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_1D( map%M, d_src_loc, d_dst_loc)

    ! Set values of border vertices to mean of interior neighbours
    ! Used to fix problems with conservative remapping on the border
    call set_border_vertices_to_interior_mean_dp_2D( mesh_dst, d_dst_loc)

    ! Support d_dst as both distributed and hybrid distributed/shared memory
    if (size( d_dst,1) == mesh_dst%pai_V%n_loc) then
      d_dst = d_dst_loc
    elseif (size( d_dst,1) == mesh_dst%pai_V%n_nih) then
      call dist_to_hybrid( mesh_dst%pai_V, d_dst_loc, d_dst)
    else
      call crash('invalid size for d_dst')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_to_mesh_2D

  subroutine apply_map_mesh_to_mesh_3D( mesh_src, mesh_dst, map, d_src, d_dst)

    ! In/output variables
    type(type_mesh),          intent(in)  :: mesh_src
    type(type_mesh),          intent(in)  :: mesh_dst
    type(type_map),           intent(in)  :: map
    real(dp), dimension(:,:), intent(in)  :: d_src
    real(dp), dimension(:,:), intent(out) :: d_dst

    ! Local variables:
    character(len=*), parameter                                      :: routine_name = 'apply_map_mesh_to_mesh_3D'
    real(dp), dimension(mesh_src%vi1:mesh_src%vi2, 1:size( d_src,2)) :: d_src_loc
    real(dp), dimension(mesh_dst%vi1:mesh_dst%vi2, 1:size( d_dst,2)) :: d_dst_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Support d_src as both distributed and hybrid distributed/shared memory
    if (size( d_src,1) == mesh_src%pai_V%n_loc) then
      d_src_loc = d_src
    elseif (size( d_src,1) == mesh_src%pai_V%n_nih) then
      call hybrid_to_dist( mesh_src%pai_V, size( d_src,2), d_src, d_src_loc)
    else
      call crash('invalid size for d_src')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_2D( map%M, d_src_loc, d_dst_loc)

    ! Set values of border vertices to mean of interior neighbours
    ! Used to fix problems with conservative remapping on the border
    call set_border_vertices_to_interior_mean_dp_3D( mesh_dst, d_dst_loc)

    ! Support d_dst as both distributed and hybrid distributed/shared memory
    if (size( d_dst,1) == mesh_dst%pai_V%n_loc) then
      d_dst = d_dst_loc
    elseif (size( d_dst,1) == mesh_dst%pai_V%n_nih) then
      call dist_to_hybrid( mesh_dst%pai_V, size( d_dst,2), d_dst_loc, d_dst)
    else
      call crash('invalid size for d_dst')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_to_mesh_3D

  subroutine apply_map_mesh_tri_to_mesh_tri_2D( mesh_src, mesh_dst, map, d_src, d_dst)

    ! In/output variables
    type(type_mesh),        intent(in)  :: mesh_src
    type(type_mesh),        intent(in)  :: mesh_dst
    type(type_map),         intent(in)  :: map
    real(dp), dimension(:), intent(in)  :: d_src
    real(dp), dimension(:), intent(out) :: d_dst

    ! Local variables:
    character(len=*), parameter                    :: routine_name = 'apply_map_mesh_tri_to_mesh_tri_2D'
    real(dp), dimension(mesh_src%ti1:mesh_src%ti2) :: d_src_loc
    real(dp), dimension(mesh_dst%ti1:mesh_dst%ti2) :: d_dst_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Support d_src as both distributed and hybrid distributed/shared memory
    if (size( d_src,1) == mesh_src%pai_Tri%n_loc) then
      d_src_loc = d_src
    elseif (size( d_src,1) == mesh_src%pai_Tri%n_nih) then
      call hybrid_to_dist( mesh_src%pai_Tri, d_src, d_src_loc)
    else
      call crash('invalid size for d_src')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_1D( map%M, d_src_loc, d_dst_loc)

    ! Set values of border triangles to mean of interior neighbours
    ! Used to fix problems with conservative remapping on the border
    call set_border_triangles_to_interior_mean_dp_2D( mesh_dst, d_dst_loc)

    ! Support d_dst as both distributed and hybrid distributed/shared memory
    if (size( d_dst,1) == mesh_dst%pai_Tri%n_loc) then
      d_dst = d_dst_loc
    elseif (size( d_dst,1) == mesh_dst%pai_Tri%n_nih) then
      call dist_to_hybrid( mesh_dst%pai_Tri, d_dst_loc, d_dst)
    else
      call crash('invalid size for d_dst')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_tri_to_mesh_tri_2D

  subroutine apply_map_mesh_tri_to_mesh_tri_3D( mesh_src, mesh_dst, map, d_src, d_dst)

    ! In/output variables
    type(type_mesh),          intent(in)  :: mesh_src
    type(type_mesh),          intent(in)  :: mesh_dst
    type(type_map),           intent(in)  :: map
    real(dp), dimension(:,:), intent(in)  :: d_src
    real(dp), dimension(:,:), intent(out) :: d_dst

    ! Local variables:
    character(len=*), parameter                                      :: routine_name = 'apply_map_mesh_tri_to_mesh_tri_3D'
    real(dp), dimension(mesh_src%ti1:mesh_src%ti2, 1:size( d_src,2)) :: d_src_loc
    real(dp), dimension(mesh_dst%ti1:mesh_dst%ti2, 1:size( d_dst,2)) :: d_dst_loc

    ! Add routine to path
    call init_routine( routine_name)

    ! Support d_src as both distributed and hybrid distributed/shared memory
    if (size( d_src,1) == mesh_src%pai_Tri%n_loc) then
      d_src_loc = d_src
    elseif (size( d_src,1) == mesh_src%pai_Tri%n_nih) then
      call hybrid_to_dist( mesh_src%pai_Tri, size( d_src,2), d_src, d_src_loc)
    else
      call crash('invalid size for d_src')
    end if

    ! Perform the mapping operation as a matrix multiplication
    call multiply_PETSc_matrix_with_vector_2D( map%M, d_src_loc, d_dst_loc)

    ! Set values of border vertices to mean of interior neighbours
    ! Used to fix problems with conservative remapping on the border
    call set_border_triangles_to_interior_mean_dp_3D( mesh_dst, d_dst_loc)

    ! Support d_dst as both distributed and hybrid distributed/shared memory
    if (size( d_dst,1) == mesh_dst%pai_Tri%n_loc) then
      d_dst = d_dst_loc
    elseif (size( d_dst,1) == mesh_dst%pai_Tri%n_nih) then
      call dist_to_hybrid( mesh_dst%pai_Tri, size( d_dst,2), d_dst_loc, d_dst)
    else
      call crash('invalid size for d_dst')
    end if

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine apply_map_mesh_tri_to_mesh_tri_3D

end module apply_maps
