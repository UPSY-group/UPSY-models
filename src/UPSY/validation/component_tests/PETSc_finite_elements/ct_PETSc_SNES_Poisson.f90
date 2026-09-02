module ct_PETSc_SNES_Poisson

#include <petsc/finclude/petscsys.h>
  use precisions, only: dp
  use iso_c_binding, only: c_bool, c_char, c_double, c_funloc, c_funptr, c_int, c_intptr_t, c_loc, c_null_char, &
    c_null_funptr, c_null_ptr, c_ptr, c_f_pointer
  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use crash_mod, only: warning
  use mpi_basic, only: par
  use mesh_types, only: type_mesh
  use netcdf_io_main, only: open_existing_netcdf_file_for_reading, setup_mesh_from_file, &
    close_netcdf_file
  use petsc, only: PetscErrorF, PETSC_COMM_SELF, PETSC_COMM_WORLD, PETSC_TRUE, PETSC_FALSE, &
    PETSC_NULL_DMLABEL, PETSC_NULL_VEC, DM_BC_ESSENTIAL, tDM, tVec, tSNES, tKSP, tPC, &
    tPetscObject, tPetscFE, tPetscDS, tDMLabel, DMDestroy, &
    VecDestroy, SNESDestroy, PetscFEDestroy, PetscFECreateDefault, PetscObjectSetName, DMSetField, &
    DMCreateDS, DMGetDS, PetscDSSetConstants, DMCreateLabel, DMGetLabel, DMPlexMarkBoundaryFaces, &
    SNESCreate, SNESSetDM, SNESSetType, SNESSetTolerances, SNESGetKSP, SNESNEWTONLS, &
    KSPSetType, KSPGetPC, KSPPREONLY, PCSetType, PCLU, DMCreateGlobalVector, VecSet, SNESSolve
  use petsc_dmplex, only: mesh_to_dmplex
  use string_module, only: strrep

  implicit none

  private

  public :: ct_solve_Poisson_eq_with_PETSc_SNES

  real(dp), dimension(4) :: poisson_domain
  integer(c_int), parameter :: dm_bc_essential_value = 1_c_int
  real(dp), parameter :: poisson_snes_relative_tolerance = 1.0e-10_dp
  real(dp), parameter :: poisson_snes_absolute_tolerance = 1.0e-10_dp
  real(dp), parameter :: poisson_snes_step_tolerance     = 1.0e-12_dp
  integer,  parameter :: poisson_snes_max_iterations     = 20
  integer,  parameter :: poisson_snes_max_function_evaluations = 1000
  character(kind=c_char), dimension(9), parameter :: boundary_label_name = [character(kind=c_char) :: &
    'b', 'o', 'u', 'n', 'd', 'a', 'r', 'y', c_null_char]

  interface
    subroutine DMPlexSetSNESLocalFEM( dm, has_boundary, ctx, ierr)
      import :: c_bool, c_intptr_t, tDM
      type(tDM),              intent(inout) :: dm
      logical(kind=c_bool),   intent(in)    :: has_boundary
      integer(c_intptr_t),    intent(in)    :: ctx
      integer,                intent(out)   :: ierr
    end subroutine DMPlexSetSNESLocalFEM

    integer(c_int) function petsc_ds_set_residual( ds, field, f0, f1) bind(C, name='PetscDSSetResidual')
      import :: c_funptr, c_int, c_intptr_t
      integer(c_intptr_t), value :: ds
      integer(c_int),      value :: field
      type(c_funptr),      value :: f0, f1
    end function petsc_ds_set_residual

    integer(c_int) function petsc_ds_set_jacobian( ds, field_test, field_trial, g0, g1, g2, g3) &
      bind(C, name='PetscDSSetJacobian')
      import :: c_funptr, c_int, c_intptr_t
      integer(c_intptr_t), value :: ds
      integer(c_int),      value :: field_test, field_trial
      type(c_funptr),      value :: g0, g1, g2, g3
    end function petsc_ds_set_jacobian

    integer(c_int) function petsc_ds_set_exact_solution( ds, field, function, ctx) &
      bind(C, name='PetscDSSetExactSolution')
      import :: c_funptr, c_int, c_intptr_t, c_ptr
      integer(c_intptr_t), value :: ds
      integer(c_int),      value :: field
      type(c_funptr),      value :: function
      type(c_ptr),          value :: ctx
    end function petsc_ds_set_exact_solution

    integer(c_int) function dm_add_boundary( dm, boundary_type, name, label, nvalues, values, field, ncomponents, &
      components, boundary_function, boundary_time_derivative, ctx, boundary_index) bind(C, name='DMAddBoundary')
      import :: c_char, c_funptr, c_int, c_intptr_t, c_ptr
      integer(c_intptr_t),                value :: dm, label
      integer(c_int),                     value :: boundary_type, nvalues, field, ncomponents
      character(kind=c_char), dimension(*), intent(in) :: name
      type(c_ptr),                        value :: values, components, ctx
      type(c_funptr),                     value :: boundary_function, boundary_time_derivative
      integer(c_int),                   intent(out) :: boundary_index
    end function dm_add_boundary
  end interface

contains

  subroutine ct_solve_Poisson_eq_with_PETSc_SNES( foldername_output, test_mesh_filenames, test_grid_filenames)

    ! In/output variables:
    character(len=*),               intent(in) :: foldername_output
    character(len=*), dimension(:), intent(in) :: test_mesh_filenames
    character(len=*), dimension(:), intent(in) :: test_grid_filenames

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'ct_solve_Poisson_eq_with_PETSc_SNES'
    integer                       :: i_mesh
    type(type_mesh), allocatable  :: mesh
    character(len=:), allocatable :: filename_mesh
    integer                       :: ncid

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '   Running PETSc SNES test: solve a simple Poisson equation'
    if (par%primary) write(0,*) ''

    ! Meshes read from files
    do i_mesh = 1, size( test_mesh_filenames)
      filename_mesh = test_mesh_filenames( i_mesh)
      allocate( mesh)
      call open_existing_netcdf_file_for_reading( filename_mesh, ncid)
      call setup_mesh_from_file( filename_mesh, ncid, mesh)
      call close_netcdf_file( ncid)
      call ct_solve_Poisson_eq_with_PETSc_SNES_on_mesh( foldername_output, mesh)
      deallocate( mesh)
    end do

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ct_solve_Poisson_eq_with_PETSc_SNES

  subroutine ct_solve_Poisson_eq_with_PETSc_SNES_on_mesh( foldername_output, mesh)

    ! In/output variables:
    character(len=*), intent(in) :: foldername_output
    type(type_mesh),  intent(in) :: mesh

    ! Local variables:
    character(len=*), parameter   :: routine_name = 'ct_solve_Poisson_eq_with_PETSc_SNES_on_mesh'
    type(tDM)                     :: dm
    type(tVec)                    :: solution
    type(tSNES)                   :: snes
    type(tPetscFE)                :: fe
    type(tPetscObject)            :: fe_object
    type(tPetscDS)                :: ds
    type(tDMLabel)                :: boundary_label
    integer                       :: ierr
    integer(c_int)                :: boundary_index
    integer(c_int), target, dimension(1) :: boundary_ids, unused_components
    integer(c_intptr_t)           :: no_context
    character(len=:), allocatable :: mesh_name_cleaned
    character(len=:), allocatable :: filename

    ! Add routine to call stack
    call init_routine( routine_name)

    if (par%primary) write(0,*) '    Using PETSc SNES to solve the Poisson equation on mesh ', trim( mesh%name), '...'

    call mesh_to_dmplex( mesh, dm)

    ! The analytical solution is defined on the mesh bounding box. For a rectangular
    ! mesh its value is zero on every boundary face.
    poisson_domain( 1) = minval( mesh%V(:,1))
    poisson_domain( 2) = minval( mesh%V(:,2))
    poisson_domain( 3) = maxval( mesh%V(:,1)) - poisson_domain( 1)
    poisson_domain( 4) = maxval( mesh%V(:,2)) - poisson_domain( 2)
    no_context = 0_c_intptr_t

    ! Attach one scalar P1 field, u, to the DMPlex.
    PetscCall( PetscFECreateDefault( PETSC_COMM_SELF, 2, 1, PETSC_TRUE, 'poisson_', -1, fe, ierr))
    PetscCall( PetscObjectSetName( fe, 'u', ierr))
    PetscObjectSpecificCast( fe_object, fe)
    PetscCall( DMSetField( dm, 0, PETSC_NULL_DMLABEL, fe_object, ierr))
    PetscCall( DMCreateDS( dm, ierr))
    PetscCall( DMGetDS( dm, ds, ierr))
    PetscCall( PetscDSSetConstants( ds, 4, poisson_domain, ierr))

    ! Define Poisson weak form
    !     Strong form: -grad^2 u = f within the domain, u = 0 on the border
    !            with: f = pi^2 (1 / w_x^2 + 1 / w_y^2) sin( pi x / w_x) sin( pi y / w_y)
    !            where [w_x,w_y] is the size of the rectangular domain.
    ierr = petsc_ds_set_residual( ds%v, 0_c_int, c_funloc( poisson_f0), c_funloc( poisson_f1))
    CHKERRQ( ierr)
    ierr = petsc_ds_set_jacobian( ds%v, 0_c_int, 0_c_int, c_null_funptr, c_null_funptr, c_null_funptr, &
      c_funloc( poisson_g3))
    CHKERRQ( ierr)
    ierr = petsc_ds_set_exact_solution( ds%v, 0_c_int, c_funloc( poisson_exact_solution), c_null_ptr)
    CHKERRQ( ierr)

    ! Mark all exterior faces and impose the manufactured solution on them. The
    ! prescribed value is identically zero for the rectangular-domain problem.
    PetscCall( DMCreateLabel( dm, 'boundary', ierr))
    PetscCall( DMGetLabel( dm, 'boundary', boundary_label, ierr))
    PetscCall( DMPlexMarkBoundaryFaces( dm, 1, boundary_label, ierr))
    boundary_ids = [1_c_int]
    unused_components = 0_c_int
    ierr = dm_add_boundary( dm%v, dm_bc_essential_value, boundary_label_name, boundary_label%v, 1_c_int, &
      c_loc( boundary_ids), 0_c_int, 0_c_int, c_loc( unused_components), c_funloc( poisson_dirichlet_boundary), &
      c_null_funptr, c_null_ptr, boundary_index)
    CHKERRQ( ierr)

    ! Let DMPlex assemble the residual and Jacobian from the PetscDS weak form.
    PetscCall( SNESCreate( PETSC_COMM_WORLD, snes, ierr))
    PetscCall( SNESSetDM( snes, dm, ierr))
    PetscCall( DMPlexSetSNESLocalFEM( dm, PETSC_FALSE, no_context, ierr))
    call configure_PETSc_SNES_for_Poisson( snes)

    PetscCall( DMCreateGlobalVector( dm, solution, ierr))
    PetscCall( VecSet( solution, 0._dp, ierr))
    PetscCall( SNESSolve( snes, PETSC_NULL_VEC, solution, ierr))

    ! ========================================================================
    ! NEXT STEPS
    ! ========================================================================
    ! 1. Run this component test on one small mesh and inspect PETSc's solver
    !    behaviour. Solver choices and tolerances are set in
    !    configure_PETSc_SNES_for_Poisson below; move those named parameters to
    !    the UPSY configuration type when this becomes model functionality.
    !
    ! 2. Add an L2-error calculation against poisson_exact_solution. PETSc's
    !    Fortran module in the current environment does not expose
    !    DMComputeL2Diff, but the function is exported by libpetsc. Follow the
    !    direct bind(C) pattern used above for the PetscDS registrations, pass
    !    dm%v and solution%v, and use poisson_exact_solution via c_funloc.
    !
    ! 3. Check convergence under mesh refinement. The default PetscFE is a P1
    !    simplex element, so the L2 error should decrease at approximately
    !    second order for this smooth manufactured solution.
    !
    ! 4. The manufactured solution and the zero Dirichlet condition assume an
    !    axis-aligned rectangular boundary. Replace them before using meshes
    !    with another outer boundary, or prescribe the exact value on that
    !    boundary instead of zero.
    !
    ! 5. This is currently a serial DMPlex experiment. Verify the behaviour
    !    with more than one MPI rank before treating it as a parallel solve;
    !    mesh_to_dmplex constructs the full mesh locally and may need explicit
    !    DMPlex distribution and ownership handling.
    !
    ! 6. The direct C bindings in this module work around PETSc routines whose
    !    legacy Fortran wrappers are absent from the installed PETSc 3.25.5
    !    package. Prefer generated PETSc Fortran bindings if a later PETSc
    !    installation supplies them.


    ! Clean up after yourself
    PetscCall( VecDestroy( solution, ierr))
    PetscCall( SNESDestroy( snes, ierr))
    PetscCall( PetscFEDestroy( fe, ierr))
    PetscCall( DMDestroy( dm, ierr))

    ! Remove routine from call stack
    call finalise_routine( routine_name)

  end subroutine ct_solve_Poisson_eq_with_PETSc_SNES_on_mesh

  subroutine configure_PETSc_SNES_for_Poisson( snes)

    type(tSNES), intent(inout) :: snes

    type(tKSP) :: ksp
    type(tPC)  :: preconditioner
    integer    :: ierr

    PetscCall( SNESSetType( snes, SNESNEWTONLS, ierr))
    call SNESSetTolerances( snes, poisson_snes_absolute_tolerance, poisson_snes_relative_tolerance, &
      poisson_snes_step_tolerance, poisson_snes_max_iterations, poisson_snes_max_function_evaluations, ierr)
    CHKERRQ( ierr)

    PetscCall( SNESGetKSP( snes, ksp, ierr))
    PetscCall( KSPSetType( ksp, KSPPREONLY, ierr))
    PetscCall( KSPGetPC( ksp, preconditioner, ierr))
    PetscCall( PCSetType( preconditioner, PCLU, ierr))

  end subroutine configure_PETSc_SNES_for_Poisson

  integer(c_int) function poisson_exact_solution( dim, time, x, ncomp, u, ctx) bind(C)

    integer(c_int),        value :: dim, ncomp
    real(c_double),        value :: time
    type(c_ptr),           value :: x, u, ctx
    real(c_double), pointer      :: x_values(:), u_values(:)

    call c_f_pointer( x, x_values, [int( dim)])
    call c_f_pointer( u, u_values, [int( ncomp)])
    u_values( 1) = sin( acos( -1._c_double) * (x_values( 1) - poisson_domain( 1)) / poisson_domain( 3)) * &
                    sin( acos( -1._c_double) * (x_values( 2) - poisson_domain( 2)) / poisson_domain( 4))
    poisson_exact_solution = 0_c_int

  end function poisson_exact_solution

  integer(c_int) function poisson_dirichlet_boundary( dim, time, x, ncomp, u, ctx) bind(C)

    integer(c_int),        value :: dim, ncomp
    real(c_double),        value :: time
    type(c_ptr),           value :: x, u, ctx
    real(c_double), pointer      :: u_values(:)

    call c_f_pointer( u, u_values, [int( ncomp)])
    u_values( 1) = 0._c_double
    poisson_dirichlet_boundary = 0_c_int

  end function poisson_dirichlet_boundary

  subroutine poisson_f0( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, x, nconstants, constants, f0) bind(C)

    integer(c_int), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants, f0
    real(c_double), value :: time
    real(c_double), pointer :: x_values(:), constants_values(:), f0_values(:)
    real(c_double)          :: pi, solution

    call c_f_pointer( x, x_values, [int( dim)])
    call c_f_pointer( constants, constants_values, [int( nconstants)])
    call c_f_pointer( f0, f0_values, [1])
    pi = acos( -1._c_double)
    solution = sin( pi * (x_values( 1) - constants_values( 1)) / constants_values( 3)) * &
               sin( pi * (x_values( 2) - constants_values( 2)) / constants_values( 4))

    ! PETSc assembles f0 * v + f1 . grad(v), hence f0 is minus the RHS.
    f0_values( 1) = -pi**2 * (1._c_double / constants_values( 3)**2 + &
      1._c_double / constants_values( 4)**2) * solution

  end subroutine poisson_f0

  subroutine poisson_f1( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, x, nconstants, constants, f1) bind(C)

    integer(c_int), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants, f1
    real(c_double), value :: time
    real(c_double), pointer :: u_x_values(:), f1_values(:)

    call c_f_pointer( u_x, u_x_values, [int( dim)])
    call c_f_pointer( f1, f1_values, [int( dim)])
    f1_values( 1:dim) = u_x_values( 1:dim)

  end subroutine poisson_f1

  subroutine poisson_g3( dim, nf, nfaux, uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, &
    time, u_tshift, x, nconstants, constants, g3) bind(C)

    integer(c_int), value :: dim, nf, nfaux, nconstants
    type(c_ptr),    value :: uoff, uoff_x, u, u_t, u_x, aoff, aoff_x, a, a_t, a_x, x, constants, g3
    real(c_double), value :: time, u_tshift
    real(c_double), pointer :: g3_values(:)
    integer                 :: d

    call c_f_pointer( g3, g3_values, [int( dim * dim)])
    g3_values( 1:dim*dim) = 0._c_double
    do d = 1, dim
      g3_values( (d - 1) * dim + d) = 1._c_double
    end do

  end subroutine poisson_g3

end module ct_PETSc_SNES_Poisson