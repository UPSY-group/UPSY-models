module laddie_mesh_output

  use precisions, only: dp
  use mpi_basic, only: par, sync
  use UPSY_main, only: UPSY
  use parameters
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, warning, crash
  use model_configuration, only: C
  use mesh_types, only: type_mesh
  use laddie_model_types, only: type_laddie_model
  use laddie_forcing_types, only: type_laddie_forcing
  use netcdf_io_main
  use reallocate_mod
  use netcdf, only: NF90_DOUBLE
  use mesh_contour, only: calc_mesh_contour
  use mpi_f08, only: MPI_WIN
  use mpi_distributed_shared_memory, only: allocate_dist_shared, deallocate_dist_shared
  use mesh_disc_apply_operators , only: map_a_b_2D, ddx_a_a_2D, ddy_a_a_2D

  implicit none

  private

  public :: create_laddie_mesh_output_file, write_to_laddie_mesh_output_file

contains

  subroutine write_to_laddie_mesh_output_file( mesh, laddie, forcing, time, is_standalone)

    ! In/output variables
    type(type_mesh),           intent(in   ) :: mesh
    type(type_laddie_model),   intent(inout) :: laddie
    type(type_laddie_forcing), intent(in   ) :: forcing
    real(dp),                  intent(in   ) :: time
    logical,                   intent(in   ) :: is_standalone

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'write_to_laddie_mesh_output_file'
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! If the mesh has been updated, create a new output file
    if (.not. laddie%mesh_output_file_matches_current_mesh) then
      call create_laddie_mesh_output_file( mesh, laddie, is_standalone)
    end if

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Writing to mesh output file "' // UPSY%stru%colour_string( trim( laddie%output_mesh_filename), 'light blue') // '"...'

    ! Open the NetCDF file
    call open_existing_netcdf_file_for_writing( laddie%output_mesh_filename, ncid)

    ! write the time to the file
    call write_time_to_file( laddie%output_mesh_filename, ncid, time)

    ! write the default data fields to the file
    call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, 'H_lad')
    call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, 'U_lad')
    call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, 'V_lad')
    call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, 'T_lad')
    call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, 'S_lad')
    call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, 'melt')

    if (is_standalone) then
      ! write all user-defined data fields to the file
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_01)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_02)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_03)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_04)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_05)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_06)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_07)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_08)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_09)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_10)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_11)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_12)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_13)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_14)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_15)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_16)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_17)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_18)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_19)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_20)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_21)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_22)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_23)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_24)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_25)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_26)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_27)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_28)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_29)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_30)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_31)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_32)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_33)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_34)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_35)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_36)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_37)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_38)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_39)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_40)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_41)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_42)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_43)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_44)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_45)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_46)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_47)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_48)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_49)
      call write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, C%choice_output_field_50)
    end if

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_laddie_mesh_output_file

  subroutine write_to_laddie_output_file_mesh_field( mesh, laddie, forcing, ncid, choice_output_field)
    !< Write a specific field to the main regional output NetCDF file - mesh version

    ! In/output variables:
    type(type_mesh),           intent(in   ) :: mesh
    type(type_laddie_model),   intent(in   ) :: laddie
    type(type_laddie_forcing), intent(in   ) :: forcing
    integer,                   intent(in   ) :: ncid
    character(len=*),          intent(in   ) :: choice_output_field

    ! Local variables:
    character(len=1024), parameter     :: routine_name = 'write_to_laddie_output_file_mesh_field'
    integer, dimension(:  ), pointer     :: mask_int => null()
    type(MPI_WIN)                        :: wmask_int
    real(dp), dimension(:  ), pointer    :: d_partial_a => null()
    real(dp), dimension(:  ), pointer    :: d_partial_b => null()
    type(MPI_WIN)                        :: wd_partial_a
    type(MPI_WIN)                        :: wd_partial_b
    call init_routine( routine_name)

    ! if no NetCDF output should be created, do nothing
    if (.not. C%do_create_netcdf_output) then
      call finalise_routine( routine_name)
      return
    end if


    ! Allocate mask_int
    call allocate_dist_shared( mask_int, wmask_int, mesh%pai_V%n_nih)
    mask_int( mesh%pai_V%i1_nih : mesh%pai_V%i2_nih) => mask_int

    call allocate_dist_shared( d_partial_a, wd_partial_a, mesh%pai_V%n_nih)
    d_partial_a( mesh%pai_V%i1_nih : mesh%pai_V%i2_nih) => d_partial_a

    call allocate_dist_shared( d_partial_b, wd_partial_b, mesh%pai_Tri%n_nih)
    d_partial_b( mesh%pai_Tri%i1_nih : mesh%pai_Tri%i2_nih) => d_partial_b

    ! Add the specified data field to the file
    select case (choice_output_field)
      case default
        call crash('unknown choice_output_field "' // trim( choice_output_field) // '"')
      case ('none')
        ! Do nothing

    ! ===== Mesh properties =====
    ! ===========================

      case ('resolution')
        ! Do nothing - this is already part of the regular mesh data; only write this to the square grid output

    ! ===== Forcing =====
    ! ===================

      ! Ice geometry
      case ('Hi')
        call write_to_field_multopt_mesh_dp_2D_notime( mesh, laddie%output_mesh_filename, ncid, 'Hi', forcing%Hi)
      case ('Hs')
        call write_to_field_multopt_mesh_dp_2D_notime( mesh, laddie%output_mesh_filename, ncid, 'Hs', forcing%Hs)
      case ('Hb')
        call write_to_field_multopt_mesh_dp_2D_notime( mesh, laddie%output_mesh_filename, ncid, 'Hb', forcing%Hb)
      case ('Hib')
        call write_to_field_multopt_mesh_dp_2D_notime( mesh, laddie%output_mesh_filename, ncid, 'Hib', forcing%Hib)
      case ('TAF')
        call write_to_field_multopt_mesh_dp_2D_notime( mesh, laddie%output_mesh_filename, ncid, 'TAF', forcing%TAF)
      case ('grounding_line')
        call write_grounding_line_to_file( laddie%output_mesh_filename, ncid, mesh, forcing)

      ! Geometry on triangles
      case ('Hs_b')
        call map_a_b_2D( mesh, forcing%Hs, d_partial_b, d_a_is_hybrid = .true., d_b_is_hybrid = .true.)
        call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, laddie%output_mesh_filename, ncid, 'Hs_b', d_partial_b, d_is_hybrid = .true.)
      case ('Hb_b')
        call map_a_b_2D( mesh, forcing%Hb, d_partial_b, d_a_is_hybrid = .true., d_b_is_hybrid = .true.)
        call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, laddie%output_mesh_filename, ncid, 'Hb_b', d_partial_b, d_is_hybrid = .true.)
      case ('Hib_b')
        call map_a_b_2D( mesh, forcing%Hib, d_partial_b, d_a_is_hybrid = .true., d_b_is_hybrid = .true.)
        call write_to_field_multopt_mesh_dp_2D_b_notime( mesh, laddie%output_mesh_filename, ncid, 'Hib_b', d_partial_b, d_is_hybrid = .true.)

      ! Geometry gradients
      case ('dHs_dx')
        call ddx_a_a_2D( mesh, forcing%Hs, d_partial_a, d_a_is_hybrid = .true., ddx_a_is_hybrid = .true.)
        call write_to_field_multopt_mesh_dp_2D_notime( mesh, laddie%output_mesh_filename, ncid, 'dHs_dx', d_partial_a, d_is_hybrid = .true.)
      case ('dHs_dy')
        call ddy_a_a_2D( mesh, forcing%Hs, d_partial_a, d_a_is_hybrid = .true., ddy_a_is_hybrid = .true.)
        call write_to_field_multopt_mesh_dp_2D_notime( mesh, laddie%output_mesh_filename, ncid, 'dHs_dy', d_partial_a, d_is_hybrid = .true.)
      case ('dHb_dx')
        call ddx_a_a_2D( mesh, forcing%Hb, d_partial_a, d_a_is_hybrid = .true., ddx_a_is_hybrid = .true.)
        call write_to_field_multopt_mesh_dp_2D_notime( mesh, laddie%output_mesh_filename, ncid, 'dHb_dx', d_partial_a, d_is_hybrid = .true.)
      case ('dHb_dy')
        call ddy_a_a_2D( mesh, forcing%Hb, d_partial_a, d_a_is_hybrid = .true., ddy_a_is_hybrid = .true.)
        call write_to_field_multopt_mesh_dp_2D_notime( mesh, laddie%output_mesh_filename, ncid, 'dHb_dy', d_partial_a, d_is_hybrid = .true.)
      case ('dHib_dx')
        call ddx_a_a_2D( mesh, forcing%Hib, d_partial_a, d_a_is_hybrid = .true., ddx_a_is_hybrid = .true.)
        call write_to_field_multopt_mesh_dp_2D_notime( mesh, laddie%output_mesh_filename, ncid, 'dHib_dx', d_partial_a, d_is_hybrid = .true.)
      case ('dHib_dy')
        call ddy_a_a_2D( mesh, forcing%Hib, d_partial_a, d_a_is_hybrid = .true., ddy_a_is_hybrid = .true.)
        call write_to_field_multopt_mesh_dp_2D_notime( mesh, laddie%output_mesh_filename, ncid, 'dHib_dy', d_partial_a, d_is_hybrid = .true.)

      ! Ice temperature
      case ('Ti')
        call write_to_field_multopt_mesh_dp_3D_notime( mesh, laddie%output_mesh_filename, ncid, 'Ti', forcing%Ti)

      ! Main ocean variables
      case ('T_ocean')
        call write_to_field_multopt_mesh_dp_3D_ocean_notime( mesh, laddie%output_mesh_filename, ncid, 'T_ocean', forcing%T_ocean)
      case ('S_ocean')
        call write_to_field_multopt_mesh_dp_3D_ocean_notime( mesh, laddie%output_mesh_filename, ncid, 'S_ocean', forcing%S_ocean)

      case ('f_coriolis')
        call write_to_field_multopt_mesh_dp_2D_notime( mesh, laddie%output_mesh_filename, ncid, 'f_coriolis', forcing%f_coriolis)

    ! ===== Masks =====
    ! =================

      case ('mask_SGD')
        where (forcing%mask_SGD)
          mask_int = 1
        elsewhere
          mask_int = 0
        end where
        call write_to_field_multopt_mesh_int_2D( mesh, laddie%output_mesh_filename, ncid, 'mask_SGD', mask_int)
      case ('mask')
        call write_to_field_multopt_mesh_int_2D( mesh, laddie%output_mesh_filename, ncid, 'mask', forcing%mask)

    ! == LADDIE ==
    ! ============

      ! Main laddie variables
      case ('H_lad')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'H_lad', laddie%now%H)
      case ('U_lad')
        call write_to_field_multopt_mesh_dp_2D_b( mesh, laddie%output_mesh_filename, ncid, 'U_lad', laddie%now%U)
      case ('V_lad')
        call write_to_field_multopt_mesh_dp_2D_b( mesh, laddie%output_mesh_filename, ncid, 'V_lad', laddie%now%V)
      case ('T_lad')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'T_lad', laddie%now%T)
      case ('S_lad')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'S_lad', laddie%now%S)

      case ('H_lad_b')
        call write_to_field_multopt_mesh_dp_2D_b( mesh, laddie%output_mesh_filename, ncid, 'H_lad_b', laddie%now%H_b, d_is_hybrid = .true.)

      ! Useful laddie fields
      case ('drho_amb')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'drho_amb', laddie%drho_amb)
      case ('drho_base')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'drho_base', laddie%drho_base)
      case ('entr')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'entr', laddie%entr)
      case ('entr_dmin')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'entr_dmin', laddie%entr_dmin)
      case ('SGD')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'SGD', laddie%SGD)
      case ('melt')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'melt', laddie%melt)
      case ('divQH')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'divQH', laddie%divQH)
      case ('divQT')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'divQT', laddie%divQT)
      case ('divQS')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'divQS', laddie%divQS)
      case ('diffT')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'diffT', laddie%diffT)
      case ('diffS')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'diffS', laddie%diffS)
      case ('viscU')
        call write_to_field_multopt_mesh_dp_2D_b( mesh, laddie%output_mesh_filename, ncid, 'viscU', laddie%viscU)
      case ('viscV')
        call write_to_field_multopt_mesh_dp_2D_b( mesh, laddie%output_mesh_filename, ncid, 'viscV', laddie%viscV)
      case ('T_base')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'T_base', laddie%T_base)
      case ('T_amb')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'T_amb', laddie%T_amb)
      case ('T_freeze')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'T_freeze', laddie%T_freeze)
      case ('u_star')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'u_star', laddie%u_star)
      case ('gamma_T')
        call write_to_field_multopt_mesh_dp_2D( mesh, laddie%output_mesh_filename, ncid, 'gamma_T', laddie%gamma_T)
      case ('divQU')
        call write_to_field_multopt_mesh_dp_2D_b( mesh, laddie%output_mesh_filename, ncid, 'divQU', laddie%divQU)
      case ('divQV')
        call write_to_field_multopt_mesh_dp_2D_b( mesh, laddie%output_mesh_filename, ncid, 'divQV', laddie%divQV)
      case ('HU_lad')
        call write_to_field_multopt_mesh_dp_2D_b( mesh, laddie%output_mesh_filename, ncid, 'HU_lad', laddie%now%H_b*laddie%now%U)
      case ('HV_lad')
        call write_to_field_multopt_mesh_dp_2D_b( mesh, laddie%output_mesh_filename, ncid, 'HV_lad', laddie%now%H_b*laddie%now%V)

    end select

    ! Deallocate mask_int
    call deallocate_dist_shared( mask_int, wmask_int)
    call deallocate_dist_shared( d_partial_a, wd_partial_a)
    call deallocate_dist_shared( d_partial_b, wd_partial_b)


    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_to_laddie_output_file_mesh_field

  subroutine create_laddie_mesh_output_file( mesh, laddie, is_standalone)

    ! In/output variables
    type(type_mesh),         intent(in   ) :: mesh
    type(type_laddie_model), intent(inout) :: laddie
    logical,                 intent(in   ) :: is_standalone

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_laddie_mesh_output_file'
    character(len=1024)            :: filename_base
    integer                        :: ncid

    ! Add routine to path
    call init_routine( routine_name)

    ! Set filename
    filename_base = trim( C%output_dir) // 'laddie_output'
    call generate_filename_XXXXXdotnc( filename_base, laddie%output_mesh_filename)

    ! Print to terminal
    if (par%primary) write(0,'(A)') '   Creating laddie output file "' // UPSY%stru%colour_string( trim( laddie%output_mesh_filename), 'light blue') // '"...'

    ! Create the NetCDF file
    call create_new_netcdf_file_for_writing( laddie%output_mesh_filename, ncid)

    ! Set up the mesh in the file
    call setup_mesh_in_netcdf_file( laddie%output_mesh_filename, ncid, mesh)

    ! Add time dimension+variable to the file
    call add_time_dimension_to_file(  laddie%output_mesh_filename, ncid)

    ! Add the default data fields to the file
    call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, 'H_lad')
    call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, 'U_lad')
    call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, 'V_lad')
    call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, 'T_lad')
    call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, 'S_lad')
    call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, 'melt')

    if (is_standalone) then
      ! Add all user-defined data fields to the file
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_01)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_02)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_03)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_04)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_05)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_06)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_07)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_08)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_09)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_10)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_11)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_12)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_13)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_14)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_15)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_16)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_17)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_18)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_19)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_20)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_21)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_22)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_23)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_24)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_25)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_26)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_27)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_28)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_29)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_30)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_31)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_32)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_33)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_34)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_35)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_36)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_37)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_38)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_39)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_40)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_41)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_42)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_43)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_44)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_45)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_46)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_47)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_48)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_49)
      call create_laddie_output_file_mesh_field( laddie%output_mesh_filename, ncid, C%choice_output_field_50)
    end if

    ! Confirm that the current output file match the current model mesh
    ! (set to false whenever a new mesh is created,
    ! and set to true whenever a new output file is created)
    laddie%mesh_output_file_matches_current_mesh = .true.

    ! Close the file
    call close_netcdf_file( ncid)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_laddie_mesh_output_file

  subroutine create_laddie_output_file_mesh_field( filename, ncid, choice_output_field)
    !< Create a single field in the main regional output NetCDF file - mesh version

    ! In/output variables:
    character(len=*), intent(in   ) :: filename
    integer,          intent(in   ) :: ncid
    character(len=*), intent(in   ) :: choice_output_field

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'create_laddie_output_file_mesh_field'
    integer                        :: int_dummy, id_dim_ei, id_dim_two, id_dim_time
    integer                        :: id_var_grounding_line

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
        ! Do nothing - this is already part of the regular mesh data; only write this to the square grid output

    ! ===== Forcing =====
    ! ===================

      ! Ice geometry
      case ('Hi')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'Hi', long_name = 'Ice thickness', units = 'm')
      case ('Hs')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'Hs', long_name = 'Ice surface elevation', units = 'm')
      case ('Hb')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'Hb', long_name = 'Bedrock elevation', units = 'm w.r.t. sea level')
      case ('Hib')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'Hib', long_name = 'Ice base elevation', units = 'm w.r.t. sea level')
      case ('TAF')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'TAF', long_name = 'Ice thickness above floatation', units = 'm w.r.t. sea level')
      case ('grounding_line')
        call inquire_dim( filename, ncid, 'ei', int_dummy, id_dim_ei)
        call inquire_dim( filename, ncid, 'two', int_dummy, id_dim_two)
        call inquire_dim( filename, ncid, 'time', int_dummy, id_dim_time)
        call create_variable( filename, ncid, 'grounding_line', NF90_DOUBLE, (/ id_dim_ei, id_dim_two/), id_var_grounding_line)
        call add_attribute_char( filename, ncid, id_var_grounding_line, 'long_name', 'Grounding-line coordinates')
        call add_attribute_char( filename, ncid, id_var_grounding_line, 'units', 'm')
        call add_attribute_char( filename, ncid, id_var_grounding_line, 'format', 'Matlab/Python contour format')

      case ('Hs_b')
        call add_field_mesh_dp_2D_b_notime( filename, ncid, 'Hs_b', long_name = 'Ice surface elevation on triangles', units = 'm')
      case ('Hb_b')
        call add_field_mesh_dp_2D_b_notime( filename, ncid, 'Hb_b', long_name = 'Bedrock elevation on triangles', units = 'm')
      case ('Hib_b')
        call add_field_mesh_dp_2D_b_notime( filename, ncid, 'Hib_b', long_name = 'Ice base elevation on triangles', units = 'm')

      case ('dHs_dx')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'dHs_dx', long_name = 'Ice surface elevation gradient x-direction', units = 'm')
      case ('dHs_dy')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'dHs_dy', long_name = 'Ice surface elevation gradient y-direction', units = 'm')
      case ('dHb_dx')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'dHb_dx', long_name = 'Bedrock elevation gradient x-direction', units = 'm')
      case ('dHb_dy')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'dHb_dy', long_name = 'Bedrock elevation gradient y-direction', units = 'm')
      case ('dHib_dx')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'dHib_dx', long_name = 'Ice base elevation gradient x-direction', units = 'm')
      case ('dHib_dy')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'dHib_dy', long_name = 'Ice base elevation gradient y-direction', units = 'm')

      ! Ice temperature
      case ('Ti')
        call add_field_mesh_dp_3D_notime( filename, ncid, 'Ti', long_name = 'Englacial temperature', units = 'deg C')

      ! Main ocean variables
      case ('T_ocean')
        call add_field_mesh_dp_3D_ocean_notime( filename, ncid, 'T_ocean', long_name = 'Ocean temperature', units = 'deg C')
      case ('S_ocean')
        call add_field_mesh_dp_3D_ocean_notime( filename, ncid, 'S_ocean', long_name = 'Ocean salinity', units = 'PSU')

      case ('f_coriolis')
        call add_field_mesh_dp_2D_notime( filename, ncid, 'f_coriolis', long_name = 'Coriolis parameter', units = 's^-1')

      case ('mask_SGD')
        call add_field_mesh_int_2d( filename, ncid, 'mask_SGD', long_name = 'Mask indicating potential subglacial discharge cells')
      case ('mask')
        call add_field_mesh_int_2d( filename, ncid, 'mask', long_name = 'General mask')

    ! == LADDIE ==
    ! ============

      ! Main laddie variables
      case ('H_lad')
        call add_field_mesh_dp_2D( filename, ncid, 'H_lad', long_name = 'Laddie layer thickness', units = 'm')
      case ('U_lad')
        call add_field_mesh_dp_2D_b( filename, ncid, 'U_lad', long_name = 'Laddie U velocity', units = 'm s^-1')
      case ('V_lad')
        call add_field_mesh_dp_2D_b( filename, ncid, 'V_lad', long_name = 'Laddie V velocity', units = 'm s^-1')
      case ('T_lad')
        call add_field_mesh_dp_2D( filename, ncid, 'T_lad', long_name = 'Laddie temperature', units = 'deg C')
      case ('S_lad')
        call add_field_mesh_dp_2D( filename, ncid, 'S_lad', long_name = 'Laddie salinity', units = 'PSU')

      case ('H_lad_b')
        call add_field_mesh_dp_2D_b( filename, ncid, 'H_lad_b', long_name = 'Laddie layer thickness on triangles', units = 'm')

      ! Useful laddie fields
      case ('drho_amb')
        call add_field_mesh_dp_2D( filename, ncid, 'drho_amb', long_name = 'Depth integrated buoyancy', units = 'kg m^-2')
      case ('drho_base')
        call add_field_mesh_dp_2D( filename, ncid, 'drho_base', long_name = 'Depth integrated buoyancy', units = 'kg m^-2')
      case ('entr')
        call add_field_mesh_dp_2D( filename, ncid, 'entr', long_name = 'Entrainment rate', units = 'm s^-1')
      case ('entr_dmin')
        call add_field_mesh_dp_2D( filename, ncid, 'entr_dmin', long_name = 'Entrainment rate for Dmin', units = 'm s^-1')
      case ('SGD')
        call add_field_mesh_dp_2D( filename, ncid, 'SGD', long_name = 'Subglacial discharge rate', units = 'm s^-1')
      case ('melt')
        call add_field_mesh_dp_2D( filename, ncid, 'melt', long_name = 'Melt rate', units = 'm s^-1')
      case ('divQH')
        call add_field_mesh_dp_2D( filename, ncid, 'divQH', long_name = 'Thickness divergence', units = 'm s^-1')
      case ('divQT')
        call add_field_mesh_dp_2D( filename, ncid, 'divQT', long_name = 'Heat divergence', units = 'degC m s^-1')
      case ('divQS')
        call add_field_mesh_dp_2D( filename, ncid, 'divQS', long_name = 'Salt divergence', units = 'PSU m s^-1')
      case ('diffT')
        call add_field_mesh_dp_2D( filename, ncid, 'diffT', long_name = 'Heat diffusion', units = 'degC m s^-1')
      case ('diffS')
        call add_field_mesh_dp_2D( filename, ncid, 'diffS', long_name = 'Salt diffusion', units = 'PSU m s^-1')
      case ('viscU')
        call add_field_mesh_dp_2D_b( filename, ncid, 'viscU', long_name = 'Laddie U viscosity', units = 'm^2 s^-2')
      case ('viscV')
        call add_field_mesh_dp_2D_b( filename, ncid, 'viscV', long_name = 'Laddie V viscosity', units = 'm^2 s^-2')
      case ('T_base')
        call add_field_mesh_dp_2D( filename, ncid, 'T_base', long_name = 'Temperature at ice/ocean interface', units = 'deg C')
      case ('T_amb')
        call add_field_mesh_dp_2D( filename, ncid, 'T_amb', long_name = 'Temperature at interface with ambient ocean', units = 'deg C')
      case ('T_freeze')
        call add_field_mesh_dp_2D( filename, ncid, 'T_freeze', long_name = 'Feezing temperature', units = 'deg C')
      case ('u_star')
        call add_field_mesh_dp_2D( filename, ncid, 'u_star', long_name = 'Friction velocity', units = 'm s^-1')
      case ('gamma_T')
        call add_field_mesh_dp_2D( filename, ncid, 'gamma_T', long_name = 'Heat exchange coefficient', units = 'm s^-1')
      case ('divQU')
        call add_field_mesh_dp_2D_b( filename, ncid, 'divQU', long_name = 'Laddie U divergence', units = 'm^2 s^-2')
      case ('divQV')
        call add_field_mesh_dp_2D_b( filename, ncid, 'divQV', long_name = 'Laddie V divergence', units = 'm^2 s^-2')
      case ('HU_lad')
        call add_field_mesh_dp_2D_b( filename, ncid, 'HU_lad', long_name = 'Laddie HU ', units = 'm^2 s^-1')
      case ('HV_lad')
        call add_field_mesh_dp_2D_b( filename, ncid, 'HV_lad', long_name = 'Laddie HV ', units = 'm^2 s^-1')

    end select

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine create_laddie_output_file_mesh_field

  subroutine write_grounding_line_to_file( filename, ncid, mesh, forcing)

    ! In/output variables:
    character(len=*),          intent(in   ) :: filename
    integer,                   intent(in   ) :: ncid
    type(type_mesh),           intent(in   ) :: mesh
    type(type_laddie_forcing), intent(in   ) :: forcing

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_grounding_line_to_file'
    real(dp), dimension(mesh%vi1:mesh%vi2)  :: TAF_for_GL
    integer                                 :: vi
    real(dp), dimension(:,:  ), allocatable :: CC

    ! Add routine to path
    call init_routine( routine_name)

    ! Replace thickness above floatation with NaN in ice-free vertices so GL wont be found there
    do vi = mesh%vi1, mesh%vi2
      if (forcing%Hi( vi) > 0.1_dp) then
        TAF_for_GL( vi) = forcing%TAF( vi)
      else
        TAF_for_GL( vi) = NaN
      end if
    end do

    ! Calculate grounding line contour
    if (par%primary) allocate( CC( mesh%nE,2))
    call calc_mesh_contour( mesh, TAF_for_GL, 0._dp, CC)

    ! Write to NetCDF
    call write_contour_to_file_notime( filename, ncid, mesh, CC, 'grounding_line')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_grounding_line_to_file

  subroutine write_contour_to_file_notime( filename, ncid, mesh, CC, var_name)

    ! In/output variables:
    character(len=*),         intent(in   ) :: filename
    integer,                  intent(in   ) :: ncid
    type(type_mesh),          intent(in   ) :: mesh
    real(dp), dimension(:,:), intent(in   ) :: CC
    character(len=*),         intent(in   ) :: var_name

    ! Local variables:
    character(len=1024), parameter          :: routine_name = 'write_contour_to_file_notime'
    real(dp), dimension(:,:,:), allocatable :: CC_with_time
    integer                                 :: id_dim_time, ti, id_var

    ! Add routine to path
    call init_routine( routine_name)

    ! Write to NetCDF
    call inquire_var( filename, ncid, var_name, id_var)
    call write_var_primary( filename, ncid, id_var, CC)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine write_contour_to_file_notime

end module laddie_mesh_output
