MODULE mesh_memory

  ! Routines for allocating, deallocating, extending and cropping the memory for the mesh data.

! ===== Preamble =====
! ====================

  USE precisions                                             , ONLY: dp
  USE mpi_basic                                              , ONLY: par, sync
  USE call_stack_and_comp_time_tracking                  , ONLY: warning, crash, happy, init_routine, finalise_routine, colour_string
  USE mesh_types                                             , ONLY: type_mesh
  USE reallocate_mod                                         , ONLY: reallocate
  USE CSR_matrix_basics                            , ONLY: deallocate_matrix_CSR_dist
  use mpi_distributed_shared_memory, only: deallocate_dist_shared

  IMPLICIT NONE

CONTAINS

! ===== Subroutines =====
! =======================

  SUBROUTINE allocate_mesh_primary( mesh, name, nV_mem, nTri_mem)
    ! Allocate memory for primary mesh data (everything thats's needed for mesh creation & refinement)

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    CHARACTER(LEN=*),                INTENT(IN)        :: name
    INTEGER,                         INTENT(IN)        :: nV_mem, nTri_mem

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'allocate_mesh_primary'

    ! Add routine to path
    CALL init_routine( routine_name)

    mesh%name     = trim(name)
    mesh%nV_mem   = nV_mem
    mesh%nTri_mem = nTri_mem

    ! Safety: check to make sure that no memory is allocated for this mesh yet
    IF (ALLOCATED( mesh%V                   ) .OR. &
        ALLOCATED( mesh%nC                  ) .OR. &
        ALLOCATED( mesh%C                   ) .OR. &
        ALLOCATED( mesh%niTri               ) .OR. &
        ALLOCATED( mesh%iTri                ) .OR. &
        ALLOCATED( mesh%VBI                 ) .OR. &
        ALLOCATED( mesh%Tri                 ) .OR. &
        ALLOCATED( mesh%Tricc               ) .OR. &
        ALLOCATED( mesh%TriC                ) .OR. &
        ALLOCATED( mesh%check_Delaunay_map  ) .OR. &
        ALLOCATED( mesh%check_Delaunay_stack) .OR. &
        ALLOCATED( mesh%refinement_map      ) .OR. &
        ALLOCATED( mesh%refinement_stack    ) .OR. &
        ALLOCATED( mesh%Tri_li              )) THEN
      CALL crash('memory is already allocated for mesh "' // TRIM( name) // '"!')
    END IF

    ! Allocate memory
    ! ===============

    ! Vertex data
    ALLOCATE( mesh%V                   (nV_mem,   2          ), source = 0._dp  )
    ALLOCATE( mesh%nC                  (nV_mem               ), source = 0      )
    ALLOCATE( mesh%C                   (nV_mem,   mesh%nC_mem), source = 0      )
    ALLOCATE( mesh%niTri               (nV_mem               ), source = 0      )
    ALLOCATE( mesh%iTri                (nV_mem,   mesh%nC_mem), source = 0      )
    ALLOCATE( mesh%VBI                 (nV_mem               ), source = 0      )

    ! Triangle data
    ALLOCATE( mesh%Tri                 (nTri_mem, 3          ), source = 0      )
    ALLOCATE( mesh%Tricc               (nTri_mem, 2          ), source = 0._dp  )
    ALLOCATE( mesh%TriC                (nTri_mem, 3          ), source = 0      )

    ! Mesh generation/refinement data
    allocate( mesh%check_Delaunay_map  (nTri_mem, 3          ), source = .false.)
    allocate( mesh%check_Delaunay_stack(nTri_mem, 4          ), source = 0      )
    mesh%check_Delaunay_stackN = 0
    ALLOCATE( mesh%refinement_map      (nTri_mem             ), source = 0      )
    ALLOCATE( mesh%refinement_stack    (nTri_mem             ), source = 0      )
    mesh%refinement_stackN = 0
    ALLOCATE( mesh%Tri_li              (nTri_mem, 2          ), source = 0      )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE allocate_mesh_primary

  SUBROUTINE extend_mesh_primary( mesh, nV_mem_new, nTri_mem_new)
    ! For when we didn't allocate enough. Field by field, copy the data to a temporary array,
    ! deallocate the old field, allocate a new (bigger) one, and copy the data back.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh
    INTEGER,                         INTENT(IN)        :: nV_mem_new, nTri_mem_new

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'extend_mesh_primary'

    ! Add routine to path
    CALL init_routine( routine_name)

    mesh%nV_mem   = nV_mem_new
    mesh%nTri_mem = nTri_mem_new

    ! Vertex data
    CALL reallocate( mesh%V                   , nV_mem_new  , 2          )
    CALL reallocate( mesh%nC                  , nV_mem_new               )
    CALL reallocate( mesh%C                   , nV_mem_new  , mesh%nC_mem)
    CALL reallocate( mesh%niTri               , nV_mem_new               )
    CALL reallocate( mesh%iTri                , nV_mem_new  , mesh%nC_mem)
    CALL reallocate( mesh%VBI                 , nV_mem_new               )

    ! Triangle data
    CALL reallocate( mesh%Tri                 , nTri_mem_new, 3          )
    CALL reallocate( mesh%Tricc               , nTri_mem_new, 2          )
    CALL reallocate( mesh%TriC                , nTri_mem_new, 3          )

    ! Mesh generation/refinement data
    CALL reallocate( mesh%check_Delaunay_map  , nTri_mem_new, 3          )
    CALL reallocate( mesh%check_Delaunay_stack, nTri_mem_new, 4          )
    CALL reallocate( mesh%refinement_map      , nTri_mem_new             )
    CALL reallocate( mesh%refinement_stack    , nTri_mem_new             )
    CALL reallocate( mesh%Tri_li              , nTri_mem_new, 2          )

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE extend_mesh_primary

  SUBROUTINE crop_mesh_primary( mesh)
    ! For when we allocated too much. Field by field, copy the data to a temporary array,
    ! deallocate the old field, allocate a new (smaller) one, and copy the data back.

    IMPLICIT NONE

    ! In/output variables:
    TYPE(type_mesh),                 INTENT(INOUT)     :: mesh

    ! Local variables:
    CHARACTER(LEN=256), PARAMETER                      :: routine_name = 'crop_mesh_primary'

    ! Add routine to path
    CALL init_routine( routine_name)

    IF (mesh%nV_mem > mesh%nV) THEN

      mesh%nV_mem   = mesh%nV

      ! Vertex data
      CALL reallocate( mesh%V               , mesh%nV,   2          )
      CALL reallocate( mesh%nC              , mesh%nV               )
      CALL reallocate( mesh%C               , mesh%nV,   mesh%nC_mem)
      CALL reallocate( mesh%niTri           , mesh%nV               )
      CALL reallocate( mesh%iTri            , mesh%nV,   mesh%nC_mem)
      CALL reallocate( mesh%VBI             , mesh%nV               )

    END IF ! IF (mesh%nV_mem > mesh%nV) THEN

    IF (mesh%nTri_mem > mesh%nTri) THEN

      mesh%nTri_mem = mesh%nTri

      ! Triangle data
      CALL reallocate( mesh%Tri                 , mesh%nTri, 3          )
      CALL reallocate( mesh%Tricc               , mesh%nTri, 2          )
      CALL reallocate( mesh%TriC                , mesh%nTri, 3          )

      ! Mesh generation/refinement data
      CALL reallocate( mesh%check_Delaunay_map  , mesh%nTri, 3          )
      CALL reallocate( mesh%check_Delaunay_stack, mesh%nTri, 4          )
      CALL reallocate( mesh%refinement_map      , mesh%nTri             )
      CALL reallocate( mesh%refinement_stack    , mesh%nTri             )
      CALL reallocate( mesh%Tri_li              , mesh%nTri, 2          )

    END IF ! IF (mesh%nTri_mem > mesh%nTri) THEN

    ! Finalise routine path
    CALL finalise_routine( routine_name)

  END SUBROUTINE crop_mesh_primary

  subroutine duplicate_mesh_primary( mesh, mesh2)
    !< Create a copy of a mesh

    ! In/output variables:
    type(type_mesh), intent(in   ) :: mesh
    type(type_mesh), intent(  out) :: mesh2

    ! Local variables:
    character(len=1024), parameter :: routine_name = 'duplicate_mesh_primary'

    ! Add routine to path
    call init_routine( routine_name)

    call allocate_mesh_primary( mesh2, trim( mesh%name) // '_copy', mesh%nV, mesh%nTri)

    ! Metadata
    mesh2%lambda_M    = mesh%lambda_M
    mesh2%phi_M       = mesh%phi_M
    mesh2%beta_stereo = mesh%beta_stereo
    mesh2%xmin        = mesh%xmin
    mesh2%xmax        = mesh%xmax
    mesh2%ymin        = mesh%ymin
    mesh2%ymax        = mesh%ymax
    mesh2%tol_dist    = mesh%tol_dist
    mesh2%nV          = mesh%nV
    mesh2%nTri        = mesh%nTri

    ! Vertex data
    mesh2%vi_SW  = mesh%vi_SW
    mesh2%vi_SE  = mesh%vi_SE
    mesh2%vi_NW  = mesh%vi_NW
    mesh2%vi_NE  = mesh%vi_NE
    mesh2%V      = mesh%V    ( 1:mesh%nV,:)
    mesh2%nC     = mesh%nC   ( 1:mesh%nV  )
    mesh2%C      = mesh%C    ( 1:mesh%nV,:)
    mesh2%niTri  = mesh%niTri( 1:mesh%nV  )
    mesh2%iTri   = mesh%iTri ( 1:mesh%nV,:)
    mesh2%VBI    = mesh%VBI  ( 1:mesh%nV  )

    ! Triangle data
    mesh2%Tri    = mesh%Tri  ( 1:mesh%nTri,:)
    mesh2%Tricc  = mesh%Tricc( 1:mesh%nTri,:)
    mesh2%TriC   = mesh%TriC ( 1:mesh%nTri,:)

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine duplicate_mesh_primary

END MODULE mesh_memory
