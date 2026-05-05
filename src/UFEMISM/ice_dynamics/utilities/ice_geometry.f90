module ice_geometry_calculations
 
    use precisions, only: dp
    use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
    use mesh_types, only: type_mesh
    use subgrid_ice_margin, only: calc_effective_thickness
    use subgrid_grounded_fractions_main, only: calc_grounded_fractions
    use masks_mod, only: determine_masks
    use ice_geometry_basics, only: ice_surface_elevation, thickness_above_floatation, height_of_water_column_at_ice_front
    use mpi_distributed_memory, only: gather_to_all
    use model_configuration, only: C

    implicit none

    private 

    public :: type_ice_geometry
 
    type type_ice_geometry
        ! Defines the derives type "ice geometry" that includes all the needed variables and the procedures that will be applied 
        
        ! Variables 
        type(type_mesh), pointer                    :: mesh 
        real(dp), dimension(:), allocatable, public :: Hi                    ! [m]       Ice thickness
        real(dp), dimension(:), allocatable, public :: Hb                    ! [m]       Bedrock elevation 
        real(dp), dimension(:), allocatable, public :: SL                    ! [m]       Water surface elevation 

        real(dp), dimension(:), allocatable, public :: Hs                    ! [m]       Ice surface elevation
        real(dp), dimension(:), allocatable, public :: Hib                   ! [m]       Ice base elevation
        real(dp), dimension(:), allocatable, public :: TAF                   ! [m]       Thickness above floatation
        real(dp), dimension(:), allocatable, public :: Ho                    ! [m]       Height of water column at ice front
        real(dp), dimension(:), allocatable, public :: Hi_eff                ! [m]       Effective ice thickness
        real(dp), dimension(:), allocatable, public :: dHb
        real(dp), dimension(:), allocatable, public :: fraction_margin       ! [0-1]     Sub-grid ice-filled fraction
        real(dp), dimension(:), allocatable, public :: fraction_gr
        real(dp), dimension(:), allocatable, public :: fraction_gr_b
        real(dp), dimension(:,:), allocatable, public :: bedrock_cdf
        real(dp), dimension(:,:), allocatable, public :: bedrock_cdf_b

        logical,  dimension(:), allocatable, public :: mask_icefree_land       ! T: ice-free land , F: otherwise
        logical,  dimension(:), allocatable, public :: mask_icefree_ocean      ! T: ice-free ocean, F: otherwise
        logical,  dimension(:), allocatable, public :: mask_grounded_ice       ! T: grounded ice  , F: otherwise
        logical,  dimension(:), allocatable, public :: mask_floating_ice       ! T: floating ice  , F: otherwise
        logical,  dimension(:), allocatable, public :: mask_margin             ! T: ice next to ice-free, F: otherwise
        logical,  dimension(:), allocatable, public :: mask_gl_gr              ! T: grounded ice next to floating ice, F: otherwise
        logical,  dimension(:), allocatable, public :: mask_gl_fl              ! T: floating ice next to grounded ice, F: otherwise
        logical,  dimension(:), allocatable, public :: mask_cf_gr              ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
        logical,  dimension(:), allocatable, public :: mask_cf_fl              ! T: floating ice next to ice-free water (sea or lake), F: otherwise
        logical,  dimension(:), allocatable, public :: mask_coastline          ! T: ice-free land next to ice-free ocean, F: otherwise
        integer,  dimension(:), allocatable, public :: mask        

        contains

        ! Procedures 
        procedure, public :: allocate_ice_geometry, calc_ice_geometry, calc_ice_geometry_primary_fields
    
    end type type_ice_geometry  
 
    contains
 
    subroutine calc_ice_geometry( self, mesh, Hi, SL, Hb)
        ! Calculates the ice geometry fields (Hs, Hib, TAF, Ho) from the primary input fields (Hi, SL, Hb), secondary fields (fraction margin ect) and calculates the masks.         

        !In/out variables
        class(type_ice_geometry),                        intent(inout) :: self
        type(type_mesh), target,                         intent(in   ) :: mesh 
        real(dp),dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: Hi
        real(dp),dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: SL
        real(dp),dimension(self%mesh%vi1:self%mesh%vi2), intent(in   ) :: Hb

        !Local variables 
        character(len=1024), parameter :: routine_name = 'calc_ice_geometry'

        call init_routine( routine_name)

        !Correspond the "self" to it's input so that the other subroutines in the ice geometry type have access to these inputs
        self%mesh => mesh
        self%Hi   = Hi 
        self%Hb   = Hb
        self%SL   = SL

        ! Apply no ice mask on Hi then calculate basic geometry
        call self%calc_ice_geometry_primary_fields()

        ! Update the masks
       call determine_masks(self%mesh, self%Hi, self%Hb, self%SL, self%mask, self%mask_icefree_land, self%mask_icefree_ocean, self%mask_grounded_ice, self%mask_floating_ice, self%mask_margin, self%mask_gl_fl, self%mask_gl_gr, self%mask_cf_gr, self%mask_cf_fl, self%mask_coastline)
        
        ! Grounded fraction 
       call calc_grounded_fractions( self%mesh, self%Hi, self%Hb, self%SL, self%dHb, self%fraction_gr, self%fraction_gr_b, self%mask_floating_ice, self%bedrock_cdf,self%bedrock_cdf_b)

        ! Fraction margin and effective thickness 
        call calc_effective_thickness(self%mesh, self%Hi, self%Hb, self%SL, self%Hi_eff, self%fraction_margin)

        call finalise_routine( routine_name)

    end subroutine calc_ice_geometry

  
    subroutine  calc_ice_geometry_primary_fields(self)
        ! From Hi (after making sure no ice mask is applied), Hb and SL, calculate Hs, Hib, TAF and Ho

        class(type_ice_geometry),intent(inout) :: self
        integer :: vi

        do vi = self%mesh%vi1, self%mesh%vi2
            self%Hs ( vi) = ice_surface_elevation( self%Hi( vi), self%Hb( vi), self%SL( vi))
            self%Hib( vi) = self%Hs( vi) - self%Hi( vi)
            self%TAF( vi) = thickness_above_floatation( self%Hi( vi), self%Hb( vi), self%SL( vi))
            self%Ho ( vi) = height_of_water_column_at_ice_front( self%Hi( vi), self%Hb( vi), self%SL( vi))
        end do

    end subroutine  calc_ice_geometry_primary_fields


    subroutine allocate_ice_geometry (self, mesh)
        
        class(type_ice_geometry), intent(inout) :: self
        type(type_mesh), target,  intent(in   ) :: mesh

        character(len=1024), parameter :: routine_name = 'allocate_ice_geometry'

        call init_routine( routine_name)

        self%mesh => mesh

        ! Check if config value is valid
        ! if (C%subgrid_bedrock_cdf_nbins <= 0) then
            ! Use a default value or print error
            ! print *, "WARNING: C%subgrid_bedrock_cdf_nbins not set, using default 11", C%subgrid_bedrock_cdf_nbins
          
        ! end if

        ! Primary input fields 
        if (.not. allocated(self%Hi)) allocate(self%Hi(self%mesh%vi1:self%mesh%vi2)) 
        if (.not. allocated(self%SL)) allocate(self%SL(self%mesh%vi1:self%mesh%vi2)) 
        if (.not. allocated(self%Hb)) allocate(self%Hb(self%mesh%vi1:self%mesh%vi2)) 

        ! Primary output fields
        if (.not. allocated(self%Hs)) allocate(self%Hs(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%Hib)) allocate(self%Hib(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%TAF)) allocate(self%TAF(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%Ho)) allocate(self%Ho(self%mesh%vi1:self%mesh%vi2))

        ! Secondary fields
        if (.not. allocated(self%Hi_eff)) allocate(self%Hi_eff(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%dHb)) allocate(self%dHb(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%fraction_margin)) allocate(self%fraction_margin(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%fraction_gr)) allocate(self%fraction_gr(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%fraction_gr_b)) allocate(self%fraction_gr_b(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%bedrock_cdf)) allocate(self%bedrock_cdf( mesh%vi1:mesh%vi2, C%subgrid_bedrock_cdf_nbins))
        if (.not. allocated(self%bedrock_cdf_b)) allocate(self%bedrock_cdf_b( mesh%ti1:mesh%ti2, C%subgrid_bedrock_cdf_nbins))
        
        ! Masks
        if (.not. allocated(self%mask)) allocate(self%mask(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%mask_icefree_land)) allocate(self%mask_icefree_land(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%mask_icefree_ocean)) allocate(self%mask_icefree_ocean(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%mask_grounded_ice)) allocate(self%mask_grounded_ice(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%mask_floating_ice)) allocate(self%mask_floating_ice(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%mask_margin)) allocate(self%mask_margin(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%mask_gl_gr)) allocate(self%mask_gl_gr(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%mask_gl_fl)) allocate(self%mask_gl_fl(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%mask_cf_gr)) allocate(self%mask_cf_gr(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%mask_cf_fl)) allocate(self%mask_cf_fl(self%mesh%vi1:self%mesh%vi2))
        if (.not. allocated(self%mask_coastline)) allocate(self%mask_coastline(self%mesh%vi1:self%mesh%vi2))

        ! Initialize arrays to zero/false
        self%Hi = 0.0_dp
        self%SL = 0.0_dp
        self%Hb = 0.0_dp
        self%Hs = 0.0_dp
        self%Hib = 0.0_dp
        self%TAF = 0.0_dp
        self%Ho = 0.0_dp
        self%Hi_eff = 0.0_dp
        self%dHb = 0.0_dp
        self%fraction_margin = 0.0_dp
        self%fraction_gr = 0.0_dp
        self%fraction_gr_b = 0.0_dp

        self%mask = 0
        self%mask_icefree_land = .false.
        self%mask_icefree_ocean = .false.
        self%mask_grounded_ice = .false.
        self%mask_floating_ice = .false.
        self%mask_margin = .false.
        self%mask_gl_gr = .false.
        self%mask_gl_fl = .false.
        self%mask_cf_gr = .false.
        self%mask_cf_fl = .false.
        self%mask_coastline = .false.

        call finalise_routine( routine_name)

    end subroutine allocate_ice_geometry
    
end module ice_geometry_calculations
 
 
 