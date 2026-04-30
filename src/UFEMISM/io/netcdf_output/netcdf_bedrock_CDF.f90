module netcdf_bedrock_CDF

  use precisions, only: dp
  use mpi_basic, only: par
  use mpi_distributed_memory, only: gather_to_primary
  use model_configuration, only: C
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use ice_model_types, only: type_ice_model
  use netcdf_basic
  use netcdf, only: NF90_DOUBLE

  implicit none

  private

  public :: setup_bedrock_CDF_in_netcdf_file

contains

  subroutine setup_bedrock_CDF_in_netcdf_file( mesh, filename, ncid, ice)
    !< Set up a bedrock CDF in an existing NetCDF file

    ! In/output variables:
    type(type_mesh),      intent(in   ) :: mesh
    character(len=*),     intent(in   ) :: filename
    integer,              intent(in   ) :: ncid
    type(type_ice_model), intent(in   ) :: ice

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'setup_bedrock_CDF_in_netcdf_file'
    integer                        :: id_dim_vi, id_dim_ti, id_dim_bin
    integer                        :: id_var_cdf, id_var_cdf_b
    real(dp), dimension(:,:), allocatable :: bedrock_cdf_tot, bedrock_cdf_b_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Gather bedrock cdfs to primary
    if (par%primary) then
      allocate( bedrock_cdf_tot  ( mesh%nV  , C%subgrid_bedrock_cdf_nbins))
      allocate( bedrock_cdf_b_tot( mesh%nTri, C%subgrid_bedrock_cdf_nbins))
    end if
    call gather_to_primary( ice%bedrock_cdf  , d_tot = bedrock_cdf_tot)
    call gather_to_primary( ice%bedrock_cdf_b, d_tot = bedrock_cdf_b_tot)

    ! Create CDF bin dimension
    call create_dimension( filename, ncid, 'bin', C%subgrid_bedrock_cdf_nbins, id_dim_bin)

    ! Inquire mesh dimensions
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nV,   id_dim_vi)
    call inquire_dim_multopt( filename, ncid, field_name_options_dim_nTri, id_dim_ti)

    ! Vertex data
    call create_variable( filename, ncid, 'bedrock_cdf', NF90_DOUBLE, (/ id_dim_vi, id_dim_bin /), id_var_cdf)
    call add_attribute_char( filename, ncid, id_var_cdf, 'long_name', 'Bedrock CDF of vertices')
    call add_attribute_char( filename, ncid, id_var_cdf, 'units'    , '%'                 )
    call write_var_primary( filename, ncid, id_var_cdf, bedrock_cdf_tot)

    ! Triangle data
    call create_variable( filename, ncid, 'bedrock_cdf_b', NF90_DOUBLE, (/ id_dim_ti, id_dim_bin /), id_var_cdf_b)
    call add_attribute_char( filename, ncid, id_var_cdf_b, 'long_name', 'Bedrock CDF of triangles')
    call add_attribute_char( filename, ncid, id_var_cdf_b, 'units', '%')
    call write_var_primary( filename, ncid, id_var_cdf_b, bedrock_cdf_b_tot)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine setup_bedrock_CDF_in_netcdf_file

end module netcdf_bedrock_CDF
