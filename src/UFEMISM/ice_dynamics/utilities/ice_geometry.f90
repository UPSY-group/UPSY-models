module ice_geometry_calculations_mod
 
    use precisions, only: dp
    use call_stack_and_comp_time_tracking, only: init_routine, finalise_routine, crash
    use mesh_types, only: type_mesh
    use subgrid_ice_margin, only: calc_effective_thickness
    use masks_mod, only: determine_masks
    use ice_geometry_basics, only: ice_surface_elevation, thickness_above_floatation, height_of_water_column_at_ice_front
    use mpi_distributed_memory, only: gather_to_all

    implicit none

    private 

    public :: type_ice_geometry
 
    type type_ice_geometry
        ! Defines the derives type "ice geometry" that includes all the needed variables and the procedures that will be applied 
        
        ! Variables 
        type(type_mesh), pointer                     :: mesh 
        real(dp), dimension(:), allocatable, private :: Hi                    ! [m]       Ice thickness
        real(dp), dimension(:), allocatable, private :: Hb                    ! [m]       Bedrock elevation 
        real(dp), dimension(:), allocatable, private :: SL                    ! [m]       Water surface elevation 

        real(dp), dimension(:), allocatable, private :: Hs                    ! [m]       Ice surface elevation
        real(dp), dimension(:), allocatable, private :: Hib                   ! [m]       Ice base elevation
        real(dp), dimension(:), allocatable, private :: TAF                   ! [m]       Thickness above floatation
        real(dp), dimension(:), allocatable, private :: Ho                    ! [m]       Height of water column at ice front
        real(dp), dimension(:), allocatable, private :: Hi_eff                ! [m]       Effective ice thickness
        real(dp), dimension(:), allocatable, private :: dHb
        real(dp), dimension(:), allocatable, private :: fraction_margin       ! [0-1]     Sub-grid ice-filled fraction
        real(dp), dimension(:), allocatable, private :: fraction_gr
        real(dp), dimension(:), allocatable, private :: fraction_gr_b
        ! real(dp), dimension(:,:), allocatable, private :: bedrock_cdf
        ! real(dp), dimension(:,:), allocatable, private :: bedrock_cdf_b

        logical,  dimension(:), allocatable, private :: mask_icefree_land       ! T: ice-free land , F: otherwise
        logical,  dimension(:), allocatable, private :: mask_icefree_ocean      ! T: ice-free ocean, F: otherwise
        logical,  dimension(:), allocatable, private :: mask_grounded_ice       ! T: grounded ice  , F: otherwise
        logical,  dimension(:), allocatable, private :: mask_floating_ice       ! T: floating ice  , F: otherwise
        logical,  dimension(:), allocatable, private :: mask_margin             ! T: ice next to ice-free, F: otherwise
        logical,  dimension(:), allocatable, private :: mask_gl_gr              ! T: grounded ice next to floating ice, F: otherwise
        logical,  dimension(:), allocatable, private :: mask_gl_fl              ! T: floating ice next to grounded ice, F: otherwise
        logical,  dimension(:), allocatable, private :: mask_cf_gr              ! T: grounded ice next to ice-free water (sea or lake), F: otherwise
        logical,  dimension(:), allocatable, private :: mask_cf_fl              ! T: floating ice next to ice-free water (sea or lake), F: otherwise
        logical,  dimension(:), allocatable, private :: mask_coastline          ! T: ice-free land next to ice-free ocean, F: otherwise
        integer,  dimension(:), allocatable, private :: mask        

        contains

        ! Procedures 
        procedure, public :: allocate_ice_geometry, calc_ice_geometry, calc_ice_geometry_primary_fields, get_ice_geometry, get_mask
    
    end type type_ice_geometry  
 
    contains
 
    subroutine calc_ice_geometry( self, mesh, Hi, SL, Hb)
        ! Calculates the ice geometry fields (Hs, Hib, TAF, Ho) from the primary input fields (Hi, SL, Hb), secondary fiels ( fraction margin ect ) and calculates the masks.         

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
       ! call determine_masks(self%mesh, self%Hi, self%Hb, self%SL, self%mask, self%mask_icefree_land, self%mask_icefree_ocean, self%mask_grounded_ice, self%mask_floating_ice, self%mask_margin, self%mask_gl_fl, self%mask_gl_gr, slef%mask_cf_gr, self%mask_cf_fl, self%mask_coastline)

        ! Grounded fraction 
       ! call calc_grounded_fractions( mesh, self%Hi, self%Hb, self%SL, self%dHb, self%fraction_gr, self%fraction_gr_b, self%mask_floating_ice, self%bedrock_cdf,self%bedrock_cdf_b)

        ! Fraction margin and effectif thickness 
        call calc_effective_thickness(mesh, self%Hi, self%Hb, self%SL, self%Hi_eff, self%fraction_margin)

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

        ! Primary input fields 
        allocate(self%Hi (self%mesh%vi1:self%mesh%vi2)) 
        allocate(self%SL (self%mesh%vi1:self%mesh%vi2)) 
        allocate(self%Hb (self%mesh%vi1:self%mesh%vi2)) 

        ! Primary output fields
        allocate(self%Hs(self%mesh%vi1:self%mesh%vi2))
        allocate(self%Hib(self%mesh%vi1:self%mesh%vi2))
        allocate(self%TAF(self%mesh%vi1:self%mesh%vi2))
        allocate(self%Ho(self%mesh%vi1:self%mesh%vi2))

        ! Secondary fields
        allocate(self%Hi_eff(self%mesh%vi1:self%mesh%vi2))
        allocate(self%dHb(self%mesh%vi1:self%mesh%vi2))
        allocate(self%fraction_margin(self%mesh%vi1:self%mesh%vi2))
        allocate(self%fraction_gr(self%mesh%vi1:self%mesh%vi2))
        allocate(self%fraction_gr_b(self%mesh%vi1:self%mesh%vi2))

        ! Masks
        allocate(self%mask(self%mesh%vi1:self%mesh%vi2))
        allocate(self%mask_icefree_land(self%mesh%vi1:self%mesh%vi2))
        allocate(self%mask_icefree_ocean(self%mesh%vi1:self%mesh%vi2))
        allocate(self%mask_grounded_ice(self%mesh%vi1:self%mesh%vi2))
        allocate(self%mask_floating_ice(self%mesh%vi1:self%mesh%vi2))
        allocate(self%mask_margin(self%mesh%vi1:self%mesh%vi2))
        allocate(self%mask_gl_gr(self%mesh%vi1:self%mesh%vi2))
        allocate(self%mask_gl_fl(self%mesh%vi1:self%mesh%vi2))
        allocate(self%mask_cf_gr(self%mesh%vi1:self%mesh%vi2))
        allocate(self%mask_cf_fl(self%mesh%vi1:self%mesh%vi2))
        allocate(self%mask_coastline(self%mesh%vi1:self%mesh%vi2))

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


    subroutine get_ice_geometry(self, name, array_out)
        ! Get a field from the ice geometry type

        class(type_ice_geometry),            intent(in ) :: self
        character(len=1024),                 intent(in ) :: name
        real(dp), dimension(:), allocatable, intent(out) :: array_out

        character(len=1024), parameter :: routine_name = 'get_ice_geometry'
        
        call init_routine(routine_name)
        
        ! Select field based on name
        select case(trim(name))
            ! Fields 
            case('Hs')
                array_out = self%Hs
            case('Hib')
                array_out = self%Hib
            case('TAF')
                array_out = self%TAF
            case('Ho')
                array_out = self%Ho
            case('Hi_eff')
                array_out = self%Hi_eff
            case('dHb')
                array_out = self%dHb
            case('fraction_margin')
                array_out = self%fraction_margin
            case('fraction_gr')
                array_out = self%fraction_gr
            case('fraction_gr_b')
                array_out = self%fraction_gr_b
                        
            case default
                call crash("Unknown field name: " // trim(name))
        end select

        call finalise_routine(routine_name)

    end subroutine get_ice_geometry


    subroutine get_mask(self, mask_name, mask_out)
        ! Get a mask from the ice geometry object

        class(type_ice_geometry),           intent(in ) :: self
        character(len=1024),                intent(in ) :: mask_name
        logical, dimension(:), allocatable, intent(out) :: mask_out
        
        character(len=1024), parameter :: routine_name = 'get_mask'
        
        call init_routine(routine_name)
        
        ! Select mask based on name
        select case(trim(mask_name))
            case('mask')
                mask_out = (self%mask == 1)  ! Convert integer to logical
            case('mask_icefree_land')
                mask_out = self%mask_icefree_land
            case('mask_icefree_ocean')
                mask_out = self%mask_icefree_ocean
            case('mask_grounded_ice')
                mask_out = self%mask_grounded_ice
            case('mask_floating_ice')
                mask_out = self%mask_floating_ice
            case('mask_margin')
                mask_out = self%mask_margin
            case('mask_gl_gr')
                mask_out = self%mask_gl_gr
            case('mask_gl_fl')
                mask_out = self%mask_gl_fl
            case('mask_cf_gr')
                mask_out = self%mask_cf_gr
            case('mask_cf_fl')
                mask_out = self%mask_cf_fl
            case('mask_coastline')
                mask_out = self%mask_coastline
            case default
                call crash("Unknown mask name: " // trim(mask_name))
        end select
        
        call finalise_routine(routine_name)
        
    end subroutine get_mask


end module ice_geometry_calculations_mod
 
 
 


