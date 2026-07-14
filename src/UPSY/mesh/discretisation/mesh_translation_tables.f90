module mesh_translation_tables

  ! Calculate translation tables relating grid points to matrix rows.

  use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine
  use mesh_types, only: type_mesh
  use allocate_dist_shared_mod, only: allocate_dist_shared
  use mpi_basic, only: par, sync

  implicit none

contains

  subroutine calc_field_to_vector_form_translation_tables( mesh)
    ! Calculate grid-cell-to-matrix-row translation tables

    ! In/output variables:
    type(type_mesh), intent(inout) :: mesh

    ! Local variables:
    character(len=*), parameter :: routine_name = 'calc_field_to_vector_form_translation_tables'
    integer                     :: nz,vi,ti,ei,k,ks,uv,n

    ! Add routine to path
    call init_routine( routine_name)

    nz = mesh%nz

    ! Grid sizes
    mesh%nna     = mesh%nV
    mesh%nnauv   = mesh%nV            * 2
    mesh%nnak    = mesh%nV   *  nz
    mesh%nnakuv  = mesh%nV   *  nz    * 2
    mesh%nnaks   = mesh%nV   * (nz-1)
    mesh%nnaksuv = mesh%nV   * (nz-1) * 2

    mesh%nnb     = mesh%nTri
    mesh%nnbuv   = mesh%nTri          * 2
    mesh%nnbk    = mesh%nTri *  nz
    mesh%nnbkuv  = mesh%nTri *  nz    * 2
    mesh%nnbks   = mesh%nTri * (nz-1)
    mesh%nnbksuv = mesh%nTri * (nz-1) * 2

    mesh%nnc     = mesh%nE
    mesh%nncuv   = mesh%nE            * 2
    mesh%nnck    = mesh%nE   *  nz
    mesh%nnckuv  = mesh%nE   *  nz    * 2
    mesh%nncks   = mesh%nE   * (nz-1)
    mesh%nncksuv = mesh%nE   * (nz-1) * 2

    ! Allocate node-shared memory
    call allocate_dist_shared( mesh%n2vi    , mesh%wn2vi    , [1, mesh%nna    ]        )
    call allocate_dist_shared( mesh%n2viuv  , mesh%wn2viuv  , [1, mesh%nnauv  ], [1, 2])
    call allocate_dist_shared( mesh%n2vik   , mesh%wn2vik   , [1, mesh%nnak   ], [1, 2])
    call allocate_dist_shared( mesh%n2vikuv , mesh%wn2vikuv , [1, mesh%nnakuv ], [1, 3])
    call allocate_dist_shared( mesh%n2viks  , mesh%wn2viks  , [1, mesh%nnaks  ], [1, 2])
    call allocate_dist_shared( mesh%n2viksuv, mesh%wn2viksuv, [1, mesh%nnaksuv], [1, 3])

    call allocate_dist_shared( mesh%n2ti    , mesh%wn2ti    , [1, mesh%nnb    ]        )
    call allocate_dist_shared( mesh%n2tiuv  , mesh%wn2tiuv  , [1, mesh%nnbuv  ], [1, 2])
    call allocate_dist_shared( mesh%n2tik   , mesh%wn2tik   , [1, mesh%nnbk   ], [1, 2])
    call allocate_dist_shared( mesh%n2tikuv , mesh%wn2tikuv , [1, mesh%nnbkuv ], [1, 3])
    call allocate_dist_shared( mesh%n2tiks  , mesh%wn2tiks  , [1, mesh%nnbks  ], [1, 2])
    call allocate_dist_shared( mesh%n2tiksuv, mesh%wn2tiksuv, [1, mesh%nnbksuv], [1, 3])

    call allocate_dist_shared( mesh%n2ei    , mesh%wn2ei    , [1, mesh%nnc    ]        )
    call allocate_dist_shared( mesh%n2eiuv  , mesh%wn2eiuv  , [1, mesh%nncuv  ], [1, 2])
    call allocate_dist_shared( mesh%n2eik   , mesh%wn2eik   , [1, mesh%nnck   ], [1, 2])
    call allocate_dist_shared( mesh%n2eikuv , mesh%wn2eikuv , [1, mesh%nnckuv ], [1, 3])
    call allocate_dist_shared( mesh%n2eiks  , mesh%wn2eiks  , [1, mesh%nncks  ], [1, 2])
    call allocate_dist_shared( mesh%n2eiksuv, mesh%wn2eiksuv, [1, mesh%nncksuv], [1, 3])

    call allocate_dist_shared( mesh%vi2n    , mesh%wvi2n    , [1, mesh%nV  ]                 )
    call allocate_dist_shared( mesh%viuv2n  , mesh%wviuv2n  , [1, mesh%nV  ],          [1, 2])
    call allocate_dist_shared( mesh%vik2n   , mesh%wvik2n   , [1, mesh%nV  ], [1, nz]        )
    call allocate_dist_shared( mesh%vikuv2n , mesh%wvikuv2n , [1, mesh%nV  ], [1, nz], [1, 2])
    call allocate_dist_shared( mesh%viks2n  , mesh%wviks2n  , [1, mesh%nV  ], [1, nz]        )
    call allocate_dist_shared( mesh%viksuv2n, mesh%wviksuv2n, [1, mesh%nV  ], [1, nz], [1, 2])

    call allocate_dist_shared( mesh%ti2n    , mesh%wti2n    , [1, mesh%nTri]                 )
    call allocate_dist_shared( mesh%tiuv2n  , mesh%wtiuv2n  , [1, mesh%nTri],          [1, 2])
    call allocate_dist_shared( mesh%tik2n   , mesh%wtik2n   , [1, mesh%nTri], [1, nz]        )
    call allocate_dist_shared( mesh%tikuv2n , mesh%wtikuv2n , [1, mesh%nTri], [1, nz], [1, 2])
    call allocate_dist_shared( mesh%tiks2n  , mesh%wtiks2n  , [1, mesh%nTri], [1, nz]        )
    call allocate_dist_shared( mesh%tiksuv2n, mesh%wtiksuv2n, [1, mesh%nTri], [1, nz], [1, 2])

    call allocate_dist_shared( mesh%ei2n    , mesh%wei2n    , [1, mesh%nE  ]                 )
    call allocate_dist_shared( mesh%eiuv2n  , mesh%weiuv2n  , [1, mesh%nE  ],          [1, 2])
    call allocate_dist_shared( mesh%eik2n   , mesh%weik2n   , [1, mesh%nE  ], [1, nz]        )
    call allocate_dist_shared( mesh%eikuv2n , mesh%weikuv2n , [1, mesh%nE  ], [1, nz], [1, 2])
    call allocate_dist_shared( mesh%eiks2n  , mesh%weiks2n  , [1, mesh%nE  ], [1, nz]        )
    call allocate_dist_shared( mesh%eiksuv2n, mesh%weiksuv2n, [1, mesh%nE  ], [1, nz], [1, 2])

    if (par%node_primary) then

    ! == a-grid (vertices)

      ! == 2-D

        ! == scalar

      n = 0
      do vi = 1, mesh%nV
        n = n+1
        mesh%vi2n( vi) = n
        mesh%n2vi( n ) = vi
      end do

        ! == vector

      n = 0
      do vi = 1, mesh%nV
        do uv = 1, 2
          n = n+1
          mesh%viuv2n( vi,uv) = n
          mesh%n2viuv( n,1) = vi
          mesh%n2viuv( n,2) = uv
        end do
      end do

      ! == 3-D regular

        ! == scalar

      n = 0
      do vi = 1, mesh%nV
        do k = 1, nz
          n = n+1
          mesh%vik2n( vi,k) = n
          mesh%n2vik( n,1) = vi
          mesh%n2vik( n,2) = k
        end do
      end do

        ! == vector

      n = 0
      do vi = 1, mesh%nV
        do k = 1, nz
          do uv = 1, 2
            n = n+1
            mesh%vikuv2n( vi,k,uv) = n
            mesh%n2vikuv( n,1) = vi
            mesh%n2vikuv( n,2) = k
            mesh%n2vikuv( n,3) = uv
          end do
        end do
      end do

      ! == 3-D staggered

        ! == scalar

      n = 0
      do vi = 1, mesh%nV
        do ks = 1, nz-1
          n = n+1
          mesh%viks2n( vi,ks) = n
          mesh%n2viks( n,1) = vi
          mesh%n2viks( n,2) = ks
        end do
      end do

        ! == vector

      n = 0
      do vi = 1, mesh%nV
        do ks = 1, nz-1
          do uv = 1, 2
            n = n+1
            mesh%viksuv2n( vi,ks,uv) = n
            mesh%n2viksuv( n,1) = vi
            mesh%n2viksuv( n,2) = ks
            mesh%n2viksuv( n,3) = uv
          end do
        end do
      end do

    ! == b-grid (triangles)

      ! == 2-D

        ! == scalar

      n = 0
      do ti = 1, mesh%nTri
        n = n+1
        mesh%ti2n( ti) = n
        mesh%n2ti( n ) = ti
      end do

        ! == vector

      n = 0
      do ti = 1, mesh%nTri
        do uv = 1, 2
          n = n+1
          mesh%tiuv2n( ti,uv) = n
          mesh%n2tiuv( n,1) = ti
          mesh%n2tiuv( n,2) = uv
        end do
      end do

      ! == 3-D regular

        ! == scalar

      n = 0
      do ti = 1, mesh%nTri
        do k = 1, nz
          n = n+1
          mesh%tik2n( ti,k) = n
          mesh%n2tik( n,1) = ti
          mesh%n2tik( n,2) = k
        end do
      end do

        ! == vector

      n = 0
      do ti = 1, mesh%nTri
        do k = 1, nz
          do uv = 1, 2
            n = n+1
            mesh%tikuv2n( ti,k,uv) = n
            mesh%n2tikuv( n,1) = ti
            mesh%n2tikuv( n,2) = k
            mesh%n2tikuv( n,3) = uv
          end do
        end do
      end do

      ! == 3-D staggered

        ! == scalar

      n = 0
      do ti = 1, mesh%nTri
        do ks = 1, nz-1
          n = n+1
          mesh%tiks2n( ti,ks) = n
          mesh%n2tiks( n,1) = ti
          mesh%n2tiks( n,2) = ks
        end do
      end do

        ! == vector

      n = 0
      do ti = 1, mesh%nTri
        do ks = 1, nz-1
          do uv = 1, 2
            n = n+1
            mesh%tiksuv2n( ti,ks,uv) = n
            mesh%n2tiksuv( n,1) = ti
            mesh%n2tiksuv( n,2) = ks
            mesh%n2tiksuv( n,3) = uv
          end do
        end do
      end do

    ! == c-grid (edges)

      ! == 2-D

        ! == scalar

      n = 0
      do ei = 1, mesh%nE
        n = n+1
        mesh%ei2n( ei) = n
        mesh%n2ei( n ) = ei
      end do

        ! == vector

      n = 0
      do ei = 1, mesh%nE
        do uv = 1, 2
          n = n+1
          mesh%eiuv2n( ei,uv) = n
          mesh%n2eiuv( n,1) = ei
          mesh%n2eiuv( n,2) = uv
        end do
      end do

      ! == 3-D regular

        ! == scalar

      n = 0
      do ei = 1, mesh%nE
        do k = 1, nz
          n = n+1
          mesh%eik2n( ei,k) = n
          mesh%n2eik( n,1) = ei
          mesh%n2eik( n,2) = k
        end do
      end do

        ! == vector

      n = 0
      do ei = 1, mesh%nE
        do k = 1, nz
          do uv = 1, 2
            n = n+1
            mesh%eikuv2n( ei,k,uv) = n
            mesh%n2eikuv( n,1) = ei
            mesh%n2eikuv( n,2) = k
            mesh%n2eikuv( n,3) = uv
          end do
        end do
      end do

      ! == 3-D staggered

        ! == scalar

      n = 0
      do ei = 1, mesh%nE
        do ks = 1, nz-1
          n = n+1
          mesh%eiks2n( ei,ks) = n
          mesh%n2eiks( n,1) = ei
          mesh%n2eiks( n,2) = ks
        end do
      end do

        ! == vector

      n = 0
      do ei = 1, mesh%nE
        do ks = 1, nz-1
          do uv = 1, 2
            n = n+1
            mesh%eiksuv2n( ei,ks,uv) = n
            mesh%n2eiksuv( n,1) = ei
            mesh%n2eiksuv( n,2) = ks
            mesh%n2eiksuv( n,3) = uv
          end do
        end do
      end do

    end if
    call sync

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_field_to_vector_form_translation_tables

end module mesh_translation_tables
