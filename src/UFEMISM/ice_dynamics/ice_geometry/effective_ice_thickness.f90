submodule( ice_geometry_model_basic) effective_ice_thickness

contains

  subroutine calc_effective_thickness( self)
    !< Determine the ice-filled fraction and effective ice thickness of floating margin pixels

    ! In- and output variables
    class(type_ice_geometry_model), intent(inout) :: self

    ! Local variables:
    character(len=*), parameter                      :: routine_name = 'calc_effective_thickness'
    integer                                          :: vi, ci, vc
    real(dp)                                         :: Hi_neighbour_max
    real(dp), dimension(self%mesh%nV)                :: Hi_tot, Hb_tot, SL_tot
    logical,  dimension(self%mesh%vi1:self%mesh%vi2) :: mask_margin, mask_floating
    logical,  dimension(self%mesh%nV)                :: mask_margin_tot, mask_floating_tot

    ! Add routine to path
    call init_routine( routine_name)

    ! Collect Hi from all processes
    call gather_to_all( self%Hi, Hi_tot)
    call gather_to_all( self%Hb, Hb_tot)
    call gather_to_all( self%SL, SL_tot)

    ! == Margin mask
    ! ==============

    ! Initialise
    mask_margin = .false.

    do vi = self%mesh%vi1, self%mesh%vi2
      do ci = 1, self%mesh%nC( vi)
        vc = self%mesh%C( vi,ci)
        if (Hi_tot( vi) > 0._dp .and. Hi_tot( vc) == 0._dp) then
          mask_margin( vi) = .true.
        end if
      end do
    end do

    ! Gather mask values from all processes
    call gather_to_all( mask_margin, mask_margin_tot)

    ! == Floating mask
    ! ================

    ! Initialise
    mask_floating = .false.

    do vi = self%mesh%vi1, self%mesh%vi2
      if (is_floating( self%Hi( vi), self%Hb( vi), self%SL( vi))) then
        mask_floating( vi) = .true.
      end if
    end do

    ! Gather mask values from all processes
    call gather_to_all( mask_floating, mask_floating_tot)

    ! == default values
    ! =================

    ! Initialise values
    do vi = self%mesh%vi1, self%mesh%vi2
      ! Grounded regions: ice-free land and grounded ice
      ! NOTE: important to let ice-free land have a non-zero
      ! margin fraction to let it get SMB in the ice thickness equation
      if (.not. mask_floating( vi)) then
        self%fraction_margin( vi) = 1._dp
        self%Hi_eff( vi) = Hi_tot( vi)
      ! Old and new floating regions
      elseif (Hi_tot( vi) > 0._dp) then
        self%fraction_margin( vi) = 1._dp
        self%Hi_eff( vi) = Hi_tot( vi)
      ! New ice-free ocean regions
      else
        self%fraction_margin( vi) = 0._dp
        self%Hi_eff( vi) = 0._dp
      end if
    end do

    ! === Compute ===
    ! ===============

    do vi = self%mesh%vi1, self%mesh%vi2

      ! Only check margin vertices
      if (.not. mask_margin_tot( vi)) then
        ! Simply use initialised values
        cycle
      end if

      ! === Max neighbour thickness ===
      ! ===============================

      ! Find the max ice thickness among non-margin neighbours
      Hi_neighbour_max = 0._dp

      do ci = 1, self%mesh%nC( vi)
        vc = self%mesh%C( vi,ci)

        ! Ignore margin neighbours
        if (mask_margin_tot( vc)) then
          cycle
        end if

        ! Floating margins check for neighbours
        if (mask_floating( vi)) then
          Hi_neighbour_max = max( Hi_neighbour_max, Hi_tot( vc))
        end if

      end do

      ! === Effective ice thickness ===
      ! ===============================

      ! Only apply if the thickest non-margin neighbour is thicker than
      ! this vertex. Otherwise, simply use initialised values.
      if (Hi_neighbour_max > Hi_tot( vi)) then
        ! Calculate sub-grid ice-filled fraction
        self%Hi_eff( vi) = Hi_neighbour_max
        self%fraction_margin( vi) = Hi_tot( vi) / Hi_neighbour_max
      end if

    end do

    call checksum( self%mesh%pai_V, self%Hi_eff         , 'Hi_eff')
    call checksum( self%mesh%pai_V, self%fraction_margin, 'fraction_margin')

    ! Finalise routine path
    call finalise_routine( routine_name)

  end subroutine calc_effective_thickness

end submodule effective_ice_thickness