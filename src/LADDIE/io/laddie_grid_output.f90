module laddie_grid_output

  use precisions, only: dp
  use mpi_basic, only: par
  use control_resources_and_error_messaging, only: init_routine, finalise_routine, colour_string, warning, crash
  use model_configuration, only: C
  use laddie_model_types, only: type_laddie_model
  use laddie_forcing_types, only: type_laddie_forcing
  use mesh_types, only: type_mesh
  use grid_types, only: type_grid
  use netcdf_io_main
  use remapping_main, only: map_from_mesh_vertices_to_xy_grid_2D, &
    map_from_mesh_vertices_to_xy_grid_3D, map_from_mesh_vertices_to_xy_grid_2D_minval, &
    map_from_mesh_triangles_to_xy_grid_2D, map_from_mesh_triangles_to_xy_grid_3D

  implicit none 

  private

  public :: create_laddie_output_file_grid, write_to_laddie_output_file_grid !, &
            !create_laddie_output_file_grid_ROI, write_to_laddie_output_file_grid_ROI

contains

  subroutine write_to_laddie_output_file_grid( mesh, laddie, forcing, time)
    !< Write to the main regional output NetCDF file - grid version

    ! In/output variables:
    type(type_mesh),           intent(in   ) :: mesh
    type(type_laddie_model),   intent(in   ) :: laddie
    type(type_laddie_forcing), intent(in   ) :: forcing
    real(dp),                  intent(in   ) :: time

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_laddie_output_file_mesh'
    integer                        :: ncid 

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then 
      call finalise_routine( routine_name)
      return
    end if

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to grid output file "' // colour_string( trim( laddie%output_grid_filename), 'light blue') // '"...' 

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( laddie%output_grid_filename, ncid)

    ! write the time to the file
    call write_time_to_file( laddie%output_grid_filename, ncid, time)

    ! write the default data fields to the file
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, 'H_lad')
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, 'U_lad')
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, 'V_lad')
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, 'T_lad')
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, 'S_lad')
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, 'melt')

    ! write all user-defined data fields to the file
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_01)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_02)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_03)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_04)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_05)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_06)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_07)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_08)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_09)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_10)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_11)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_12)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_13)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_14)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_15)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_16)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_17)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_18)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_19)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_20)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_21)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_22)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_23)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_24)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_25)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_26)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_27)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_28)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_29)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_30)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_31)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_32)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_33)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_34)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_35)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_36)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_37)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_38)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_39)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_40)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_41)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_42)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_43)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_44)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_45)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_46)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_47)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_48)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_49)
    call write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, C%choice_output_field_50)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_laddie_output_file_grid

  subroutine write_to_laddie_output_file_grid_field( mesh, laddie, forcing, ncid, choice_output_field)
    !< Write a single field to the main regional output NetCDF file - grid version

    ! In/output variables:
    type(type_mesh),           intent(in   ) :: mesh
    type(type_laddie_model),   intent(in   ) :: laddie
    type(type_laddie_forcing), intent(in   ) :: forcing
    integer,                   intent(in   ) :: ncid
    character(len=*),          intent(in   ) :: choice_output_field

    ! Local variables:
    character(len=1024), parameter        :: routine_name = 'write_to_laddie_output_file_grid_field'
    real(dp), dimension(:),   allocatable :: d_mesh_vec_partial_2D
    real(dp), dimension(:),   allocatable :: d_grid_vec_partial_2D
    real(dp), dimension(:,:), allocatable :: d_grid_vec_partial_2D_monthly
    real(dp), dimension(:,:), allocatable :: d_grid_vec_partial_3D
    real(dp), dimension(:,:), allocatable :: d_grid_vec_partial_3D_ocean
    real(dp), dimension(:),   allocatable :: mask_int
    type(type_grid)                       :: grid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    grid = laddie%output_grid

    ! allocate memory
    allocate( d_mesh_vec_partial_2D( mesh%vi1:mesh%vi2))
    allocate( d_grid_vec_partial_2D(         grid%n_loc                ))
    allocate( d_grid_vec_partial_2D_monthly( grid%n_loc, 12            ))
    allocate( d_grid_vec_partial_3D(         grid%n_loc, mesh%nz))
    allocate( d_grid_vec_partial_3D_ocean(   grid%n_loc, C%nz_ocean    ))
    allocate( mask_int( mesh%vi1:mesh%vi2), source = 0._dp)

    ! Add the specified data field to the file
    select case (choice_output_field)
      case default
        call crash('unknown choice_output_field "' // trim( choice_output_field) // '"')
      case ('none')
        ! Do nothing

    ! ===== Mesh properties =====
    ! ===========================

      case ('resolution')
        d_mesh_vec_partial_2D = mesh%R( mesh%vi1:mesh%vi2)
        call map_from_mesh_vertices_to_xy_grid_2D_minval( mesh, grid, C%output_dir, d_mesh_vec_partial_2D, d_grid_vec_partial_2D)
        call write_to_field_multopt_grid_dp_2D_notime( grid, laddie%output_grid_filename, ncid, 'resolution', d_grid_vec_partial_2D)

    ! ===== Forcing =====
    ! ===================

      ! Ice geometry
      case ('Hi')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, forcing%Hi, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D_notime( grid, laddie%output_grid_filename, ncid, 'Hi', d_grid_vec_partial_2D)
      case ('Hs')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, forcing%Hs, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D_notime( grid, laddie%output_grid_filename, ncid, 'Hs', d_grid_vec_partial_2D)
      case ('Hb')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, forcing%Hb, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D_notime( grid, laddie%output_grid_filename, ncid, 'Hb', d_grid_vec_partial_2D)
      case ('Hib')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, forcing%Hib, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D_notime( grid, laddie%output_grid_filename, ncid, 'Hib', d_grid_vec_partial_2D)
      case ('TAF')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, forcing%TAF, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D_notime( grid, laddie%output_grid_filename, ncid, 'TAF', d_grid_vec_partial_2D)

      ! Ice temperature
      case ('Ti')
        call map_from_mesh_vertices_to_xy_grid_3D( mesh, grid, C%output_dir, forcing%Ti, d_grid_vec_partial_3D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_3D_notime( grid, laddie%output_grid_filename, ncid, 'Ti', d_grid_vec_partial_3D)

      ! Main ocean variables
      case ('T_ocean')
        call map_from_mesh_vertices_to_xy_grid_3D( mesh, grid, C%output_dir, forcing%T_ocean, d_grid_vec_partial_3D_ocean, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_3D_ocean( grid, laddie%output_grid_filename, ncid, 'T_ocean', d_grid_vec_partial_3D_ocean)
      case ('S_ocean')
        call map_from_mesh_vertices_to_xy_grid_3D( mesh, grid, C%output_dir, forcing%S_ocean, d_grid_vec_partial_3D_ocean, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_3D_ocean( grid, laddie%output_grid_filename, ncid, 'S_ocean', d_grid_vec_partial_3D_ocean)

      case ('f_coriolis')
        call map_from_mesh_triangles_to_xy_grid_2D( mesh, grid, C%output_dir, forcing%f_coriolis, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D_notime( grid, laddie%output_grid_filename, ncid, 'f_coriolis', d_grid_vec_partial_2D)
    ! ===== Masks =====
    ! =================

      ! NOTE: logical/integer fields cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      case ('mask_icefree_land')
      case ('mask_icefree_ocean')
      case ('mask_grounded_ice')
      case ('mask_floating_ice')
      case ('mask_gl_fl')
      case ('mask_SGD')
      case ('mask')


      case ('grounding_line')
        ! Do nothing; only written to mesh files
      case ('ice_margin')
        ! Do nothing; only written to mesh files
      case ('calving_front')
        ! Do nothing; only written to mesh files
      case ('coastline')
        ! Do nothing; only written to mesh files
      case ('grounded_ice_contour')
        ! Do nothing; only written to mesh files

    ! == LADDIE ==
    ! ============

      ! Main laddie variables
      case ('H_lad')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%now%H, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'H_lad', d_grid_vec_partial_2D)
      case ('U_lad')
        call map_from_mesh_triangles_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%now%U, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'U_lad', d_grid_vec_partial_2D)
      case ('V_lad')
        call map_from_mesh_triangles_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%now%V, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'V_lad', d_grid_vec_partial_2D)
      case ('T_lad')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%now%T, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'T_lad', d_grid_vec_partial_2D)
      case ('S_lad')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%now%S, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'S_lad', d_grid_vec_partial_2D)

      ! Useful laddie fields
      case ('drho_amb')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%drho_amb, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'drho_amb', d_grid_vec_partial_2D)
      case ('drho_base')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%drho_base, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'drho_base', d_grid_vec_partial_2D)
      case ('entr')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%entr, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'entr', d_grid_vec_partial_2D)
      case ('entr_dmin')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%entr_dmin, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'entr_dmin', d_grid_vec_partial_2D)
      case ('SGD')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%SGD, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'SGD', d_grid_vec_partial_2D)
      case ('melt')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%melt, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'melt', d_grid_vec_partial_2D)
      case ('divQH')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%divQH, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'divQH', d_grid_vec_partial_2D)
      case ('divQT')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%divQT, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'divQT', d_grid_vec_partial_2D)
      case ('divQS')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%divQS, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'divQS', d_grid_vec_partial_2D)
      case ('diffT')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%diffT, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'diffT', d_grid_vec_partial_2D)
      case ('diffS')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%diffS, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'diffS', d_grid_vec_partial_2D)
      case ('T_base')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%T_base, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'T_base', d_grid_vec_partial_2D)
      case ('T_amb')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%T_amb, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'T_amb', d_grid_vec_partial_2D)
      case ('T_freeze')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%T_freeze, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'T_freeze', d_grid_vec_partial_2D)
      case ('u_star')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%u_star, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'u_star', d_grid_vec_partial_2D)
      case ('gamma_T')
        call map_from_mesh_vertices_to_xy_grid_2D( mesh, grid, C%output_dir, laddie%gamma_T, d_grid_vec_partial_2D, d_mesh_is_hybrid = .true.)
        call write_to_field_multopt_grid_dp_2D( grid, laddie%output_grid_filename, ncid, 'gamma_T', d_grid_vec_partial_2D)

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_laddie_output_file_grid_field

  subroutine create_laddie_output_file_grid( mesh, laddie, forcing)
    !< Create the main regional output NetCDF file - grid version

    ! In/output variables:
    type(type_mesh),           intent(in   ) :: mesh
    type(type_laddie_model),   intent(inout) :: laddie
    type(type_laddie_forcing), intent(in   ) :: forcing

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_laddie_output_file_grid'
    integer                        :: ncid
    type(type_grid)                       :: grid

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    grid = laddie%output_grid

    ! Set the filename
    laddie%output_grid_filename = trim( C%output_dir) // 'laddie_output_grid.nc'

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating grid output file "' // colour_string( trim( laddie%output_grid_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( laddie%output_grid_filename, ncid)

    ! Set up the grid in the file
    call setup_xy_grid_in_netcdf_file( laddie%output_grid_filename, ncid, laddie%output_grid)

    ! Add time, zeta, and month dimensions+variables to the file
    call add_time_dimension_to_file(  laddie%output_grid_filename, ncid)
    call add_zeta_dimension_to_file(  laddie%output_grid_filename, ncid, mesh%zeta)
    call add_month_dimension_to_file( laddie%output_grid_filename, ncid)
    call add_depth_dimension_to_file( laddie%output_grid_filename, ncid, C%z_ocean)

    ! Add the default data fields to the file
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, 'H_lad')
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, 'U_lad')
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, 'V_lad')
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, 'T_lad')
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, 'S_lad')
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, 'melt')

    ! Add all user-defined data fields to the file
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_01)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_02)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_03)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_04)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_05)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_06)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_07)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_08)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_09)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_10)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_11)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_12)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_13)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_14)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_15)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_16)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_17)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_18)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_19)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_20)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_21)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_22)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_23)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_24)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_25)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_26)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_27)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_28)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_29)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_30)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_31)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_32)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_33)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_34)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_35)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_36)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_37)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_38)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_39)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_40)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_41)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_42)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_43)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_44)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_45)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_46)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_47)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_48)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_49)
    call create_laddie_output_file_grid_field( laddie%output_grid_filename, ncid, C%choice_output_field_50)

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_laddie_output_file_grid

  subroutine create_laddie_output_file_grid_field( filename, ncid, choice_output_field)
    !< Create a single field in the main regional output NetCDF file - grid version

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: choice_output_field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_laddie_output_file_grid_field'

    ! Add routine to path
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if

    ! Add the specified data field to the file
    select case (choice_output_field)
      case default
        call crash('unknown choice_output_field "' // trim( choice_output_field) // '"')
      case ('none')
        ! Do nothing

    ! ===== Mesh properties =====
    ! ===========================

      case ('resolution')
        call add_field_grid_dp_2D_notime( filename, ncid, 'resolution', long_name = 'Mesh resolution (distance to nearest neighbour)', units = 'm')

    ! ===== Forcing =====
    ! ===================

      ! Ice geometry
      case ('Hi')
        call add_field_grid_dp_2D_notime( filename, ncid, 'Hi', long_name = 'Ice thickness', units = 'm')
      case ('Hs')
        call add_field_grid_dp_2D_notime( filename, ncid, 'Hs', long_name = 'Ice surface elevation', units = 'm')
      case ('Hb')
        call add_field_grid_dp_2D_notime( filename, ncid, 'Hb', long_name = 'Bedrock elevation', units = 'm')
      case ('Hib')
        call add_field_grid_dp_2D_notime( filename, ncid, 'Hib', long_name = 'Ice base elevation', units = 'm w.r.t. sea level')
      case ('TAF')
        call add_field_grid_dp_2D_notime( filename, ncid, 'TAF', long_name = 'Ice thickness above floatation', units = 'm w.r.t. sea level')
      case ('Ti')
        call add_field_grid_dp_3D_notime( filename, ncid, 'Ti', long_name = 'Englacial temperature', units = 'K')
      case ('T_ocean')
        call add_field_grid_dp_3D_ocean_notime( filename, ncid, 'T_ocean', long_name = 'Ocean temperature', units = 'deg C')
      case ('S_ocean')
        call add_field_grid_dp_3D_ocean_notime( filename, ncid, 'S_ocean', long_name = 'Ocean salinity', units = 'psu')

      case ('f_coriolis')
        call add_field_grid_dp_2D_notime( filename, ncid, 'f_coriolis', long_name = 'Coriolis parameter', units = 's^-1')

    ! ===== Masks =====
    ! =================

      ! NOTE: logical/integer fields cannot be meaningfully remapped;
      !       if you want these as gridded data, you will have to compute
      !       them yourself in post-processing

      case ('mask_icefree_land')
      case ('mask_icefree_ocean')
      case ('mask_grounded_ice')
      case ('mask_floating_ice')
      case ('mask_gl_fl')
      case ('mask_SGD')
      case ('mask')


      case ('grounding_line')
        ! Do nothing; only written to mesh files
      case ('ice_margin')
        ! Do nothing; only written to mesh files
      case ('calving_front')
        ! Do nothing; only written to mesh files
      case ('coastline')
        ! Do nothing; only written to mesh files
      case ('grounded_ice_contour')
        ! Do nothing; only written to mesh files

    ! == LADDIE ==
    ! ============

      ! Main laddie variables
      case ('H_lad')
        call add_field_grid_dp_2D( filename, ncid, 'H_lad', long_name = 'Laddie layer thickness', units = 'm')
      case ('U_lad')
        call add_field_grid_dp_2D( filename, ncid, 'U_lad', long_name = 'Laddie U velocity', units = 'm s^-1')
      case ('V_lad')
        call add_field_grid_dp_2D( filename, ncid, 'V_lad', long_name = 'Laddie V velocity', units = 'm s^-1')
      case ('T_lad')
        call add_field_grid_dp_2D( filename, ncid, 'T_lad', long_name = 'Laddie temperature', units = 'deg C')
      case ('S_lad')
        call add_field_grid_dp_2D( filename, ncid, 'S_lad', long_name = 'Laddie salinity', units = 'PSU')

      ! Useful laddie fields
      case ('drho_amb')
        call add_field_grid_dp_2D( filename, ncid, 'drho_amb', long_name = 'Depth integrated buoyancy', units = 'kg m^-2')
      case ('drho_base')
        call add_field_grid_dp_2D( filename, ncid, 'drho_base', long_name = 'Depth integrated buoyancy', units = 'kg m^-2')
      case ('entr')
        call add_field_grid_dp_2D( filename, ncid, 'entr', long_name = 'Entrainment rate', units = 'm s^-1')
      case ('entr_dmin')
        call add_field_grid_dp_2D( filename, ncid, 'entr_dmin', long_name = 'Entrainment rate for Dmin', units = 'm s^-1')
      case ('SGD')
        call add_field_grid_dp_2D( filename, ncid, 'SGD', long_name = 'Subglacial discharge rate', units = 'm s^-1')
      case ('melt')
        call add_field_grid_dp_2D( filename, ncid, 'melt', long_name = 'melt rate', units = 'm s^-1')
      case ('divQH')
        call add_field_grid_dp_2D( filename, ncid, 'divQH', long_name = 'Thickness divergence', units = 'm s^-1')
      case ('divQT')
        call add_field_grid_dp_2D( filename, ncid, 'divQT', long_name = 'Heat divergence', units = 'degC m s^-1')
      case ('divQS')
        call add_field_grid_dp_2D( filename, ncid, 'divQS', long_name = 'Salt divergence', units = 'PSU m s^-1')
      case ('diffT')
        call add_field_grid_dp_2D( filename, ncid, 'diffT', long_name = 'Heat diffusion', units = 'degC m s^-1')
      case ('diffS')
        call add_field_grid_dp_2D( filename, ncid, 'diffS', long_name = 'Salt diffusion', units = 'PSU m s^-1')
      case ('T_base')
        call add_field_grid_dp_2D( filename, ncid, 'T_base', long_name = 'Temperature at ice/ocean interface', units = 'deg C')
      case ('T_amb')
        call add_field_grid_dp_2D( filename, ncid, 'T_amb', long_name = 'Temperature at interface with ambient ocean', units = 'deg C')
      case ('T_freeze')
        call add_field_grid_dp_2D( filename, ncid, 'T_freeze', long_name = 'Freezing temperature', units = 'deg C')
      case ('u_star')
        call add_field_grid_dp_2D( filename, ncid, 'u_star', long_name = 'Friction velocity', units = 'm s^-1')
      case ('gamma_T')
        call add_field_grid_dp_2D( filename, ncid, 'gamma_T', long_name = 'Heat exchange coefficient', units = 'm s^-1')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_laddie_output_file_grid_field


end module laddie_grid_output
