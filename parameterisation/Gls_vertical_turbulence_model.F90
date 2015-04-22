!    Copyright (C) 2006 Imperial College London and others.
!    
!    Please see the AUTHORS file in the main source directory for a full list
!    of copyright holders.
!
!    Prof. C Pain
!    Applied Modelling and Computation Group
!    Department of Earth Science and Engineering
!    Imperial College London
!
!    amcgsoftware@imperial.ac.uk
!    
!    This library is free software; you can redistribute it and/or
!    modify it under the terms of the GNU Lesser General Public
!    License as published by the Free Software Foundation,
!    version 2.1 of the License.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied arranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module gls
  use quadrature
  use elements
  use field_derivatives
  use fields
  use sparse_matrices_fields
  use state_module
  use spud
  use global_parameters, only:   OPTION_PATH_LEN
  use equation_of_state
  use state_fields_module
  use boundary_conditions
  use Coordinates
  use FLDebug
  use fefields
  use vertical_extrapolation_module

  implicit none

  private

  ! These variables are the parameters requried by GLS. 
  ! They are all private to prevent tampering
  ! and save'd so that we don't need to call init every time some GLS-y
  ! happens. These are *all* private
  real, save               :: gls_n, gls_m, gls_p
  real, save               :: sigma_psi, sigma_k, kappa
  integer, save            :: nNodes
  type(scalar_field), save :: tke_old, ll
  type(scalar_field), save :: local_tke ! our local copy of TKE. We amend this
  !to add the values of the Dirichlet BC onto the field for calculating the
  !diagnostic quantities and for output. See Warner et al 2005.
  type(scalar_field), save :: MM2, NN2, eps, Fwall, S_H, S_M
  type(scalar_field), save :: K_H, K_M, density, P, B
  real, save               :: eps_min = 1e-10, psi_min, k_min
  real, save               :: cm0, cde
  !  the a_i's for the ASM
  real, save               :: a1,a2,a3,a4,a5
  real, save               :: at1,at2,at3,at4,at5
  real, save               :: cc1
  real, save               :: ct1,ctt
  real, save               :: cc2,cc3,cc4,cc5,cc6
  real, save               :: ct2,ct3,ct4,ct5
  real, save               :: cPsi1,cPsi2,cPsi3,cPsi3_plus,cPsi3_minus
  real, save               :: relaxation


  ! these are the fields and variables for the surface values
  type(scalar_field), save             :: top_surface_values, bottom_surface_values ! these are used to populate the bcs
  type(scalar_field), save             :: top_surface_tke_values, bottom_surface_tke_values ! for the Psi BC
  type(scalar_field), save             :: top_surface_km_values, bottom_surface_km_values ! for the Psi BC
  integer, save                        :: NNodes_sur, NNodes_bot
  integer, dimension(:), pointer, save :: bottom_surface_nodes, top_surface_nodes
  integer, dimension(:), pointer, save :: top_surface_element_list, bottom_surface_element_list
  logical, save                        :: calculate_bcs, calc_fwall
  character(len=FIELD_NAME_LEN), save  :: gls_wall_option, gls_stability_option, gls_option

  ! Switch for on sphere simulations to rotate the required tensors
  logical, save :: on_sphere
  
  ! The following are the public subroutines
  public :: gls_init, gls_cleanup, gls_tke, gls_diffusivity, gls_psi, gls_adapt_mesh, gls_check_options

  ! General plan is:
  !  - Init in main/Fluids.F90
  !  - Populate_State and BoundaryConditionsFromOptions also contain some set up
  !  routines such as looking for the GLS fields and setting up the automatic
  !  boundary conditions
  !  - If solve is about to do TKE, call gls_tke (which calculates NN, MM and set source/absorption for solve)
  !  - If solve is about to do Psi, call gls_psi (which fixes TKE surfaces, set source/absorption for solve)
  !  - After Psi solve, call gls_diffusivity, which sets the diffusivity and viscosity, via the lengthscale
  !  - If we adapt, call gls_adapt, which deallocates and re-allocates the
  !  fields on the new mesh. This also sets up the diagnostic fields again,
  !  which aren't interpolated
  !  - When done, clean-up

contains

!----------
! gls_init does the following:
!    - check we have the right fields (if not abort)
!    - initialise GLS parameters based on options
!    - allocate space for optional fields, which are module level variables (to save passing them around)
!----------
subroutine gls_init(state)

    type(state_type), intent(inout) :: state

    real                           :: N,rad,rcm,cmsf
    integer                        :: stat
    type(scalar_field), pointer    :: psi, tke

    psi => extract_scalar_field(state, "GLSGenericSecondQuantity")
    tke => extract_scalar_field(state, "GLSTurbulentKineticEnergy")

    ! Allocate the temporary, module-level variables
    call gls_allocate_temps(state)

    ! Check if we're on the sphere
    on_sphere = have_option('/geometry/spherical_earth/')
   
    ! populate some useful variables from options
    calculate_bcs = have_option("/material_phase[0]/subgridscale_parameterisations/GLS/calculate_boundaries/")
    calc_fwall = .false.
    call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                    &scalar_field::GLSTurbulentKineticEnergy/prognostic/minimum_value", k_min, stat)

    ! these lot are global
    call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/relax_diffusivity", relaxation, default=0.0)
    call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/stability_function", gls_stability_option)
    call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/option", gls_option)    
    ! there are lots of alternative formulae for this wall function, so let the
    ! user choose!
    if (have_option("/material_phase[0]/subgridscale_parameterisations/GLS/wall_function")) then
        call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/wall_function",gls_wall_option)
    else
        gls_wall_option = "none"
    end if

    ! Check the model used - we have four choices - then set the parameters appropriately
    select case (gls_option)
    case ("k-kl")  ! Mellor-Yamada 2.5
        gls_p = 0.0
        gls_m = 1.0
        gls_n = 1.0
        sigma_k = 2.44   ! turbinent Schmidt number
        sigma_psi = 2.44  ! turbulent Schmidt number
        cPsi1 = 0.9
        cPsi2 = 0.5
        cPsi3_plus = 1.0
        ! c3 depends on which stability function has been choosen
        select case (trim(gls_stability_option))
        case ("KanthaClayson-94")
            cPsi3_minus = 2.53
        case ("Canuto-01-A")
            cPsi3_minus = 2.681
        case ("Canuto-01-B")
            FLExit("GLS - Stability function combination not supported")
        case ("GibsonLaunder-78")
            FLExit("GLS - Stability function combination not supported")
        end select
        psi_min = 1.e-8
        calc_fwall = .true.
    case ("k-epsilon")
        gls_p = 3.0
        gls_m = 1.5
        gls_n = -1.0
        sigma_k = 1.3  ! turbulent Schmidt number
        sigma_psi = 1.0  ! turbulent Schmidt number
        cPsi1 = 1.44
        cPsi2 = 1.92
        cPsi3_plus = 1.0 
        select case (trim(gls_stability_option))
        case ("KanthaClayson-94")
            cPsi3_minus = -0.41
        case ("Canuto-01-A")
            cPsi3_minus = -0.63
        case ("Canuto-01-B")
            cPsi3_minus = -0.57
        case ("GibsonLaunder-78")
            cPsi3_minus = -0.3700
        end select
        psi_min = 1.e-12
    case ("k-omega")
        gls_p = -1.0
        gls_m = 0.5
        gls_n = -1.0
        sigma_k = 2.0  ! turbulent Schmidt number
        sigma_psi = 2.0  ! turbulent Schmidt number
        cPsi1 = 0.555
        cPsi2 = 0.833
        cPsi3_plus = 1.0 
        select case (trim(gls_stability_option))
        case ("KanthaClayson-94")
            cPsi3_minus = -0.58
        case ("Canuto-01-A")
            cPsi3_minus = -0.64
        case ("Canuto-01-B")
            cPsi3_minus = -0.61 
        case ("GibsonLaunder-78")
            cPsi3_minus = -0.4920
        end select
        psi_min = 1.e-12
    case ("gen")
        gls_p = 2.0
        gls_m = 1.0
        gls_n = -0.67
        sigma_k = 0.8  ! turbulent Schmidt number
        sigma_psi = 1.07  ! turbulent Schmidt number
        cPsi1 = 1.0
        cPsi2 = 1.22
        cPsi3_plus = 1.0 
        select case (trim(gls_stability_option))
        case ("KanthaClayson-94")
            cPsi3_minus = 0.1
        case ("Canuto-01-A")
            cPsi3_minus = 0.05
        case ("Canuto-01-B")
            cPsi3_minus = 0.08
        case ("GibsonLaunder-78")
            cPsi3_minus = 0.1704
        end select
        psi_min = 1.e-12
    case default
        FLAbort("Unknown gls_option")           
    end select

    select case (trim(gls_stability_option))
    case ("KanthaClayson-94")
         ! parameters for Kantha and Clayson (2004)
         cc1 =  6.0000
         cc2 =  0.3200
         cc3 =  0.0000
         cc4 =  0.0000
         cc5 =  0.0000
         cc6 =  0.0000
         ct1 =  3.7280
         ct2 =  0.7000
         ct3 =  0.7000
         ct4 =  0.0000
         ct5 =  0.2000
         ctt =  0.6102
    case("Canuto-01-B")
         cc1 =  5.0000
         cc2 =  0.6983
         cc3 =  1.9664
         cc4 =  1.0940
         cc5 =  0.0000
         cc6 =  0.4950
         ct1 =  5.6000
         ct2 =  0.6000
         ct3 =  1.0000
         ct4 =  0.0000
         ct5 =  0.3333
         ctt =  0.4770
    case("Canuto-01-A")
         cc1 =  5.0000
         cc2 =  0.8000
         cc3 =  1.9680
         cc4 =  1.1360
         cc5 =  0.0000
         cc6 =  0.4000
         ct1 =  5.9500
         ct2 =  0.6000
         ct3 =  1.0000
         ct4 =  0.0000
         ct5 =  0.3333
         ctt =  0.720     
     case("GibsonLaunder-78")
         cc1 =  3.6000
         cc2 =  0.8000
         cc3 =  1.2000
         cc4 =  1.2000
         cc5 =  0.0000
         cc6 =  0.5000
         ct1 =  3.0000
         ct2 =  0.3333
         ct3 =  0.3333
         ct4 =  0.0000
         ct5 =  0.3333
         ctt =  0.8000
    case default
        FLAbort("Unknown gls_stability_function") 
    end select

    ! compute the a_i's for the Algebraic Stress Model
    a1 =  2./3. - cc2/2.
    a2 =  1.    - cc3/2.
    a3 =  1.    - cc4/2.
    a4 =          cc5/2.
    a5 =  1./2. - cc6/2.

    at1 =          1. - ct2
    at2 =          1. - ct3
    at3 = 2. *   ( 1. - ct4)
    at4 = 2. *   ( 1. - ct5)
    at5 = 2.*ctt*( 1. - ct5)
    
    ! compute cm0
    N   =  cc1/2.
    cm0 =  ( (a2**2. - 3.*a3**2. + 3.*a1*N)/(3.* N**2.) )**0.25
    cmsf  =   a1/N/cm0**3
    rad=sigma_psi*(cPsi2-cPsi1)/(gls_n**2.)
    kappa = 0.41
    if (gls_option .ne. "k-kl") then    
        kappa=cm0*sqrt(rad)
    end if
    rcm  = cm0/cmsf
    cde = cm0**3.

    ewrite(1,*) "GLS Parameters"
    ewrite(1,*) "--------------------------------------------"
    ewrite(1,*) "cm0: ",cm0
    ewrite(1,*) "kappa: ",kappa
    ewrite(1,*) "p: ",gls_p
    ewrite(1,*) "m: ",gls_m
    ewrite(1,*) "n: ",gls_n
    ewrite(1,*) "sigma_k: ",sigma_k
    ewrite(1,*) "sigma_psi: ",sigma_psi
    ewrite(1,*) "Calculating BCs: ", calculate_bcs
    ewrite(1,*) "Using wall function: ", gls_wall_option
    ewrite(1,*) "Smoothing NN2: ", have_option('/material_phase[0]/subgridscale_parameterisations/GLS/smooth_buoyancy/')
    ewrite(1,*) "Smoothing MM2: ", have_option('/material_phase[0]/subgridscale_parameterisations/GLS/smooth_shear/')
    ewrite(1,*) "--------------------------------------------"

    ! intilise surface - only if we need to though
    if (calculate_bcs) then
        call gls_init_surfaces(state)
    end if

    call gls_init_diagnostics(state)
    
    call gls_calc_wall_function(state)

    ! we're all done!
end subroutine gls_init

!----------
! Update TKE
!----------
subroutine gls_tke(state)

    type(state_type), intent(inout)  :: state 
    
    type(scalar_field), pointer      :: source, absorption, scalarField
    type(tensor_field), pointer      :: tke_diff, background_diff
    type(scalar_field), pointer      :: tke
    real                             :: prod, buoyan, diss
    integer                          :: i, stat, ele
    character(len=FIELD_NAME_LEN)    :: bc_type
    type(scalar_field), pointer      :: scalar_surface
    type(vector_field), pointer      :: positions
    type(scalar_field), pointer      :: lumped_mass
    type(scalar_field)               :: inverse_lumped_mass


    ! Temporary tensor to hold  rotated values if on the sphere (note: must be a 3x3 mat)
    real, dimension(3,3) :: K_M_sphere_node

    ewrite(1,*) "In gls_tke"

    ! Get N^2 and M^2 -> NN2 and MM2
    call gls_buoyancy(state)

    do i=1,NNodes_sur
        call set(NN2,top_surface_nodes(i),0.0)
    end do

    ! calculate stability function
    call gls_stability_function(state)

    source => extract_scalar_field(state, "GLSTurbulentKineticEnergySource")
    absorption  => extract_scalar_field(state, "GLSTurbulentKineticEnergyAbsorption")
    tke => extract_scalar_field(state, "GLSTurbulentKineticEnergy")
    tke_diff => extract_tensor_field(state, "GLSTurbulentKineticEnergyDiffusivity")
    positions => extract_vector_field(state, "Coordinate")

    ! Create a local_tke in which we can mess about with the surface values
    ! for creating some of the diagnostic values later
    call set(tke,local_tke)
    ! Assembly loop
    call zero(P)
    call zero(B)
    call allocate(inverse_lumped_mass, P%mesh, "InverseLumpedMass")
    lumped_mass => get_lumped_mass(state, P%mesh)
    call invert(lumped_mass, inverse_lumped_mass)
    ! create production terms, P, B
    do ele=1, ele_count(P)
        call assemble_tke_prodcution_terms(ele, P, B, mesh_dim(P))
    end do
    call scale(P,inverse_lumped_mass)
    call scale(B,inverse_lumped_mass)
    call deallocate(inverse_lumped_mass)

    call zero(source)
    call zero(absorption)
    do ele = 1, ele_count(tke)
        call assemble_kk_src_abs(ele,tke, mesh_dim(tke))
    end do
    call allocate(inverse_lumped_mass, tke%mesh, "InverseLumpedMass")
    lumped_mass => get_lumped_mass(state, tke%mesh)
    call invert(lumped_mass, inverse_lumped_mass)
    ! source and absorption terms are set, apart from the / by lumped mass
    call scale(source,inverse_lumped_mass)
    call scale(absorption,inverse_lumped_mass)
    call deallocate(inverse_lumped_mass)

    ! set diffusivity for tke
    call zero(tke_diff)
    background_diff => extract_tensor_field(state, "GLSBackgroundDiffusivity")
    if (on_sphere) then
      do i=1,nNodes
        K_M_sphere_node=align_with_radial(node_val(positions,i),node_val(K_M,i))
        K_M_sphere_node=K_M_sphere_node*1./sigma_k
        call set(tke_diff,i,K_M_sphere_node)
      end do
    else
      call set(tke_diff,tke_diff%dim(1),tke_diff%dim(2),K_M,scale=1./sigma_k)
    end if
    call addto(tke_diff,background_diff) 

    ! boundary conditions
    if (calculate_bcs) then
        ewrite(1,*) "Calculating BCs"
        call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/calculate_boundaries/", bc_type)
        call gls_tke_bc(state,bc_type)
        ! above puts the BC boundary values in top_surface_values and bottom_surface_values module level variables
        ! map these onto the actual BCs in tke
        scalar_surface => extract_surface_field(tke, 'tke_bottom_boundary', "value")
        call remap_field(bottom_surface_values, scalar_surface)
        scalar_surface => extract_surface_field(tke, 'tke_top_boundary', "value")
        call remap_field(top_surface_values, scalar_surface)
    end if
    
    ! finally, we need a copy of the old TKE for Psi, so grab it before we solve
    call set(tke_old,tke)

    ! that's the TKE set up ready for the solve which is the next thing to happen (see Fluids.F90)
    ewrite_minmax(source)
    ewrite_minmax(absorption)
    ! set source and absorption terms in optional output fields
    scalarField => extract_scalar_field(state, "GLSSource1", stat)
    if(stat == 0) then
       call set(scalarField,source)  
    end if
    scalarField => extract_scalar_field(state, "GLSAbsorption1", stat)
    if(stat == 0) then
       call set(scalarField,absorption)  
    end if

    contains

    subroutine assemble_tke_prodcution_terms(ele, P, B, dim)

        integer, intent(in)               :: ele, dim
        type(scalar_field), intent(inout) :: P, B

        real, dimension(ele_loc(P,ele),ele_ngi(P,ele),dim)  :: dshape_P
        real, dimension(ele_ngi(P,ele))                     :: detwei
        real, dimension(ele_loc(P,ele))                     :: rhs_addto_vel, rhs_addto_buoy
        type(element_type), pointer                         :: shape_p
        integer, pointer, dimension(:)                      :: nodes_p
    
        nodes_p => ele_nodes(p, ele)
        shape_p => ele_shape(p, ele)
        call transform_to_physical( positions, ele, shape_p, dshape=dshape_p, detwei=detwei )

        ! Shear production term:
        rhs_addto_vel = shape_rhs(shape_p, detwei*ele_val_at_quad(K_M,ele)*ele_val_at_quad(MM2,ele))
        ! Buoyancy production term:
        rhs_addto_buoy = shape_rhs(shape_p, -detwei*ele_val_at_quad(K_H,ele)*ele_val_at_quad(NN2,ele))
 
        call addto(P, nodes_p, rhs_addto_vel)
        call addto(B, nodes_p, rhs_addto_buoy)
    
    end subroutine assemble_tke_prodcution_terms

    subroutine assemble_kk_src_abs(ele, kk, dim)

        integer, intent(in)            :: ele, dim
        type(scalar_field), intent(in) :: kk

        real, dimension(ele_loc(kk,ele),ele_ngi(kk,ele),dim) :: dshape_kk
        real, dimension(ele_ngi(kk,ele))                     :: detwei
        real, dimension(ele_loc(kk,ele))                     :: rhs_addto_disip, rhs_addto_src
        type(element_type), pointer                          :: shape_kk
        integer, pointer, dimension(:)                       :: nodes_kk
    
        nodes_kk => ele_nodes(kk, ele)
        shape_kk => ele_shape(kk, ele)
        call transform_to_physical( positions, ele, shape_kk, dshape=dshape_kk, detwei=detwei )

        ! ROMS and GOTM hide the absorption term in the source if the
        ! total is > 0. Kinda hard to do that over an element (*what* should
        ! be > 0?). So we don't pull this trick. Doesn't seem to make a
        ! difference to anything.
        rhs_addto_src = shape_rhs(shape_kk, detwei * (&
                                 (ele_val_at_quad(P,ele))) &
                                 )
        rhs_addto_disip = shape_rhs(shape_kk, detwei * ( &
                                 (ele_val_at_quad(eps,ele) - &
                                 ele_val_at_quad(B,ele)) / &
                                 ele_val_at_quad(tke,ele)) &
                                 )
           
        call addto(source, nodes_kk, rhs_addto_src)
        call addto(absorption, nodes_kk, rhs_addto_disip)


    end subroutine assemble_kk_src_abs

end subroutine gls_tke


!----------
! Calculate the second quantity
!----------
subroutine gls_psi(state)

    type(state_type), intent(inout)  :: state
    
    type(scalar_field), pointer      :: source, absorption, tke,  psi, scalarField
    type(tensor_field), pointer      :: psi_diff, background_diff
    real                             :: prod, buoyan,diss,PsiOverTke
    character(len=FIELD_NAME_LEN)    :: bc_type
    integer                          :: i, stat, ele
    type(scalar_field), pointer      :: scalar_surface
    type(vector_field), pointer      :: positions
    type(scalar_field), pointer      :: lumped_mass
    type(scalar_field)               :: inverse_lumped_mass, vel_prod, buoy_prod
    ! variables for the ocean parameterisation
    type(csr_matrix)                   :: face_normal_gravity
    integer, dimension(:), allocatable :: ordered_elements
    logical, dimension(:), allocatable :: node_list
    logical                            :: got_surface
    real                               :: lengthscale, percentage, tke_surface
    type(scalar_field), pointer        :: distanceToTop, distanceToBottom
    type(vector_field), pointer        :: vertical_normal

    ! Temporary tensor to hold  rotated values (note: must be a 3x3 mat)
    real, dimension(3,3)             :: psi_sphere_node
    real                             :: src, absn

    ewrite(1,*) "In gls_psi"

    source => extract_scalar_field(state, "GLSGenericSecondQuantitySource")
    absorption  => extract_scalar_field(state, "GLSGenericSecondQuantityAbsorption")
    tke => extract_scalar_field(state, "GLSTurbulentKineticEnergy")
    psi => extract_scalar_field(state, "GLSGenericSecondQuantity")
    psi_diff => extract_tensor_field(state, "GLSGenericSecondQuantityDiffusivity")
    positions => extract_vector_field(state, "Coordinate")
    vertical_normal => extract_vector_field(state, "GravityDirection")
    call allocate(vel_prod, psi%mesh, "_vel_prod_psi")
    call allocate(buoy_prod, psi%mesh, "_buoy_prod_psi")
    ewrite(2,*) "In gls_psi: setting up"

    ! store the tke in an internal field and then
    ! add the dirichlet conditions to the upper and lower surfaces. This
    ! helps stabilise the diffusivity (e.g. rapid heating cooling of the surface
    ! can destabilise the run)
    call set(local_tke,tke)

    ! clip at k_min
    do i=1,nNodes
        call set(tke,i, max(node_val(tke,i),k_min))
        call set(tke_old,i, max(node_val(tke_old,i),k_min))
        call set(local_tke,i, max(node_val(local_tke,i),k_min))
    end do

    ! This is the extra term meant to add in internal wave breaking and the
    ! like. Based on a similar term in NEMO
    if (have_option("/material_phase[0]/subgridscale_parameterisations/GLS/ocean_parameterisation")) then

        distanceToTop => extract_scalar_field(state, "DistanceToTop")
        distanceToBottom => extract_scalar_field(state, "DistanceToBottom")

        call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/ocean_parameterisation/lengthscale",lengthscale)
        call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/ocean_parameterisation/percentage",percentage)
        ewrite(2,*) "Computing extra ocean parameterisation"
        allocate(node_list(NNodes))
        node_list = .false.
        ! create gravity face normal
        call compute_face_normal_gravity(face_normal_gravity, positions, vertical_normal)
        allocate(ordered_elements(size(face_normal_gravity,1)))

        ! create an element ordering from the mesh that moves vertically downwards
        call vertical_element_ordering(ordered_elements, face_normal_gravity)
        ! I assume the above fails gracefully if the mesh isn't suitable, hence
        ! no options checks are carried out

        got_surface = .false.
        do i=1, size(ordered_elements)
            if (.not. got_surface) then
                ! First time around grab, the surface TKE
                ! for this column
                tke_surface = maxval(ele_val(tke,i))
                got_surface = .true.
            end if
            if (got_surface) then
                call ocean_tke(i,tke,distanceToTop,lengthscale,percentage,tke_surface, node_list)
            end if
            if (minval(ele_val(distanceToBottom,i)) < 1e-6) then
                got_surface = .false.
            end if
        end do

        deallocate(ordered_elements)
        call deallocate(face_normal_gravity)
        deallocate(node_list)

    end if

    call allocate(inverse_lumped_mass, psi%mesh, "InverseLumpedMass")
    lumped_mass => get_lumped_mass(state, psi%mesh)
    call invert(lumped_mass, inverse_lumped_mass)
    ! Set Psi from previous timestep
    do i=1,nNodes
        call set(psi,i, cm0**gls_p * node_val(tke_old,i)**gls_m * node_val(ll,i)**gls_n)
        call set(psi,i,max(node_val(psi,i),psi_min))
    end do

    ewrite(2,*) "In gls_psi: computing RHS"
    call zero(vel_prod)
    call zero(buoy_prod)
    call zero(source)
    call zero(absorption)    
    do ele = 1, ele_count(psi)
        call assemble_production_terms_psi(ele, vel_prod, buoy_prod, psi, mesh_dim(psi))
    end do
    call scale(vel_prod,inverse_lumped_mass)
    call scale(buoy_prod,inverse_lumped_mass)

    do ele = 1, ele_count(psi)
        call assemble_psi_src_abs(ele, psi, tke_old, mesh_dim(psi))
    end do
    call scale(source,inverse_lumped_mass)
    call scale(absorption,inverse_lumped_mass)
    call deallocate(inverse_lumped_mass)
    
    ewrite(2,*) "In gls_psi: setting diffusivity"
    ! Set diffusivity for Psi
    call zero(psi_diff)
    background_diff => extract_tensor_field(state, "GLSBackgroundDiffusivity")
    if (on_sphere) then
      do i=1,nNodes
        psi_sphere_node=align_with_radial(node_val(positions,i),node_val(K_M,i))
        psi_sphere_node=psi_sphere_node*1./sigma_psi
        call set(psi_diff,i,psi_sphere_node)
      end do
    else
      call set(psi_diff,psi_diff%dim(1),psi_diff%dim(2),K_M,scale=1./sigma_psi)
    end if
    call addto(psi_diff,background_diff) 

    ewrite(2,*) "In gls_psi: setting BCs"
    ! boundary conditions
    if (calculate_bcs) then
        call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/calculate_boundaries/", bc_type)
        call gls_psi_bc(state,bc_type)
        ! above puts the BC boundary values in top_surface_values and bottom_surface_values module level variables
        ! map these onto the actual BCs in Psi
        scalar_surface => extract_surface_field(psi, 'psi_bottom_boundary', "value")
        call remap_field(bottom_surface_values, scalar_surface)
        scalar_surface => extract_surface_field(psi, 'psi_top_boundary', "value")
        call remap_field(top_surface_values, scalar_surface)
    end if


    ewrite(2,*) "In gls_psi: tearing down"
    ! Psi is now ready for solving (see Fluids.F90)
    ewrite_minmax(source)
    ewrite_minmax(absorption)
    ! set source and absorption terms in optional output fields
    scalarField => extract_scalar_field(state, "GLSSource2", stat)
    if(stat == 0) then
       call set(scalarField,source)  
    end if
    scalarField => extract_scalar_field(state, "GLSAbsorption2", stat)
    if(stat == 0) then
       call set(scalarField,absorption)  
    end if

    call deallocate(vel_prod)
    call deallocate(buoy_prod)


    contains 
    
    subroutine ocean_tke(ele, tke, distanceToTop, lengthscale, percentage, tke_surface, nodes_done)
        type(scalar_field),pointer, intent(in)  :: distanceToTop
        type(scalar_field),pointer, intent(out) :: tke
        real, intent(in)                        :: lengthscale, percentage, tke_surface
        integer, intent(in)                     :: ele
        logical, dimension(:), intent(inout)    :: nodes_done

        integer, dimension(:), pointer  :: element_nodes
        integer                         :: i, node
        real                            :: current_TKE, depth

        element_nodes => ele_nodes(tke, ele)


        ! smooth out TKE according to length scale
        do i = 1, size(element_nodes)
            node = element_nodes(i)
            depth = node_val(distanceToTop,node)
            current_TKE = node_val(tke,node)
            if (nodes_done(node)) then
                cycle
            end if
            current_TKE = current_TKE + &
                           & percentage*TKE_surface * EXP( -depth / lengthscale )
            call set(tke,node,current_TKE)  

            nodes_done(node) = .true.
        end do

    end subroutine ocean_tke

    subroutine reconstruct_psi(ele, psi, dim)

        integer, intent(in)               :: ele, dim
        type(scalar_field), intent(inout) :: psi

        real, dimension(ele_loc(psi,ele),ele_ngi(psi,ele),dim) :: dshape_psi
        real, dimension(ele_ngi(psi,ele))                      :: detwei
        real, dimension(ele_loc(psi,ele))                      :: rhs_addto
        type(element_type), pointer                            :: shape_psi
        integer, pointer, dimension(:)                         :: nodes_psi
    
        nodes_psi => ele_nodes(psi, ele)
        shape_psi => ele_shape(psi, ele)
        call transform_to_physical( positions, ele, shape_psi, dshape=dshape_psi, detwei=detwei )

        rhs_addto = shape_rhs(shape_psi, detwei* &
                                         (cm0**gls_p) * &
                                         ele_val_at_quad(tke_old,ele)**gls_m * &
                                         ele_val_at_quad(ll,ele)**gls_n)

        call addto(psi, nodes_psi, rhs_addto)

    end subroutine reconstruct_psi

    subroutine assemble_production_terms_psi(ele, vel_prod, buoy_prod, psi, dim)


        integer, intent(in)               :: ele, dim
        type(scalar_field), intent(inout) :: psi, vel_prod, buoy_prod

        real, dimension(ele_loc(psi,ele),ele_ngi(psi,ele),dim) :: dshape_psi
        real, dimension(ele_ngi(psi,ele))                      :: detwei
        real, dimension(ele_loc(psi,ele))                      :: rhs_addto_vel, rhs_addto_buoy
        real, dimension(ele_ngi(psi,ele))                      :: cPsi3
        type(element_type), pointer                            :: shape_psi
        integer, pointer, dimension(:)                         :: nodes_psi
        
        nodes_psi => ele_nodes(psi, ele)
        shape_psi => ele_shape(psi, ele)
        call transform_to_physical( positions, ele, shape_psi, dshape=dshape_psi, detwei=detwei )

        ! Buoyancy production term:
        ! First we need to work out if cPsi3 is for stable or unstable
        ! stratification
        where(ele_val_at_quad(B,ele) .gt. 0.0)
            cPsi3 = cPsi3_plus  ! unstable strat
        elsewhere
            cPsi3 = cPsi3_minus ! stable strat
        end where
        rhs_addto_buoy = shape_rhs(shape_psi, detwei*(cPsi3*ele_val_at_quad(B,ele)*&
                                   (ele_val_at_quad(psi, ele)/ele_val_at_quad(local_tke,ele))))

        ! shear production term:
        rhs_addto_vel = shape_rhs(shape_psi, detwei*(cPsi1*ele_val_at_quad(P,ele)*&
                                    (ele_val_at_quad(psi, ele)/ele_val_at_quad(local_tke,ele))))
 
        call addto(vel_prod, nodes_psi, rhs_addto_vel)
        call addto(buoy_prod, nodes_psi, rhs_addto_buoy)


    end subroutine assemble_production_terms_psi

    subroutine assemble_psi_src_abs(ele, psi, tke, dim)

        integer, intent(in)            :: ele, dim
        type(scalar_field), intent(in) :: psi, tke

        real, dimension(ele_loc(psi,ele),ele_ngi(psi,ele),dim) :: dshape_psi
        real, dimension(ele_ngi(psi,ele))                      :: detwei
        real, dimension(ele_loc(psi,ele))                      :: rhs_addto_src, rhs_addto_disip
        type(element_type), pointer                            :: shape_psi
        integer, pointer, dimension(:)                         :: nodes_psi
    
        nodes_psi => ele_nodes(psi, ele)
        shape_psi => ele_shape(psi, ele)
        call transform_to_physical( positions, ele, shape_psi, dshape=dshape_psi, detwei=detwei )

        where (ele_val_at_quad(vel_prod,ele) + ele_val_at_quad(buoy_prod,ele) .gt. 0)
            rhs_addto_src = shape_rhs(shape_psi, detwei* ( &
                                (ele_val_at_quad(vel_prod,ele)) + &
                                (ele_val_at_quad(buoy_prod,ele)) &
                                 ) & !detwei
                                 ) ! shape_rhs
            rhs_addto_disip = shape_rhs(shape_psi, detwei* (&
                                  (cPsi2*ele_val_at_quad(eps,ele) * &
                                  (ele_val_at_quad(Fwall,ele)*ele_val_at_quad(psi, ele)/ele_val_at_quad(tke,ele))) / &
                                  ele_val_at_quad(psi,ele) &
                                  ) & ! detwei
                                  ) !shape_rhs
        elsewhere
            rhs_addto_src = shape_rhs(shape_psi, detwei * ( &
                                (ele_val_at_quad(vel_prod,ele)) &
                                 ) & !detwei
                                 ) !shape_rhs
            rhs_addto_disip = shape_rhs(shape_psi, detwei * (&
                                  ((cPsi2*ele_val_at_quad(eps,ele) * &
                                  (ele_val_at_quad(Fwall,ele)*ele_val_at_quad(psi, ele)/ele_val_at_quad(tke,ele))) - &! disipation term
                                  (ele_val_at_quad(buoy_prod,ele)))/ & ! buoyancy term
                                  ele_val_at_quad(psi,ele)&
                                  ) & !detwei 
                                  ) !shape_rhs
        end where
        call addto(source, nodes_psi, rhs_addto_src)
        call addto(absorption, nodes_psi, rhs_addto_disip)

    end subroutine assemble_psi_src_abs

end subroutine gls_psi


!----------
! gls_diffusivity fixes the top/bottom boundaries of Psi
! then calulates the lengthscale, and then uses those to calculate the 
! diffusivity and viscosity
! These are placed in the GLS fields ready for other tracer fields to use
! Viscosity is placed in the velocity viscosity
!----------
subroutine gls_diffusivity(state)

    type(state_type), intent(inout)  :: state
    
    type(scalar_field), pointer      :: tke_state, psi
    type(tensor_field), pointer      :: eddy_diff_KH,eddy_visc_KM,viscosity,background_diff,background_visc
    real                             :: exp1, exp2, exp3, x
    integer                          :: i, stat
    real                             :: psi_limit, tke_cur, limit, epslim
    real, parameter                  :: galp = 0.748331 ! sqrt(0.56)
    type(vector_field), pointer      :: positions, velocity
    type(scalar_field)               :: remaped_K_M, tke
    type(tensor_field)               :: remaped_background_visc


    ! Temporary tensors to hold  rotated values (note: must be a 3x3 mat)
    real, dimension(3,3) :: eddy_diff_KH_sphere_node, eddy_visc_KM_sphere_node, viscosity_sphere_node

    ewrite(1,*) "In gls_diffusivity"

    tke_state => extract_scalar_field(state, "GLSTurbulentKineticEnergy")
    psi => extract_scalar_field(state, "GLSGenericSecondQuantity")
    eddy_visc_KM  => extract_tensor_field(state, "GLSEddyViscosityKM",stat)
    eddy_diff_KH => extract_tensor_field(state, "GLSEddyDiffusivityKH",stat)
    viscosity => extract_tensor_field(state, "Viscosity",stat)
    positions => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "Velocity")

    call allocate(tke, tke_state%mesh, name="MyLocalTKE")
    !if (gls_n > 0) then
        ! set the TKE to use below to the unaltered TKE
        ! with no changes to the upper/lower surfaces
        ! Applies to k-kl model only
    !    call set (tke, local_tke)
    !else
        ! Use the altered TKE to get the surface diffusivity correct
        call set (tke, tke_state)
    !end if

    exp1 = 3.0 + gls_p/gls_n
    exp2 = 1.5 + gls_m/gls_n
    exp3 =       - 1.0/gls_n

    if (gls_n > 0) then
        do i=1,nNodes
            tke_cur = node_val(tke,i)
            psi_limit = (sqrt(0.56) * tke_cur**(exp2) * (1./sqrt(max(node_val(NN2,i)+1e-10,0.))) &
                        & * cm0**(gls_p / gls_n))**(-gls_n)
            call set(psi,i,max(psi_min,min(node_val(psi,i),psi_limit)))
        end do
    end if
   
    do i=1,nNodes

        tke_cur = node_val(tke,i)

        ! recover dissipation rate from k and psi
        call set(eps,i, cm0**exp1 * tke_cur**exp2 * node_val(psi,i)**exp3)

        ! limit dissipation rate under stable stratification,
        ! see Galperin et al. (1988)
        if (node_val(NN2,i) > 0) then
            epslim = (cde*tke_cur*sqrt(node_val(NN2,i)))/galp
        else
            epslim = eps_min
        end if
        call set(eps,i, max(node_val(eps,i),max(eps_min,epslim)))

        ! compute dissipative scale
        call set(ll,i,cde*sqrt(tke_cur**3.)/node_val(eps,i))
        !if (gls_n > 0) then
        !    if (node_val(NN2,i) > 0) then
        !        limit = sqrt(0.56 * tke_cur / node_val(NN2,i))
        !        call set(ll,i,min(limit,node_val(ll,i)))
        !    end if
        !end if

    end do

    ! calc fwall
    ewrite(2,*) "Calculating the wall function for GLS"
    call gls_calc_wall_function(state)

    ! calculate diffusivities for next step and for use in other fields
    do i=1,nNodes
        x = sqrt(node_val(tke,i))*node_val(ll,i)
        ! momentum
        call set(K_M,i, relaxation*node_val(K_M,i) + (1-relaxation)*node_val(S_M,i)*x)
        ! tracer
        call set(K_H,i, relaxation*node_val(K_H,i) + (1-relaxation)*node_val(S_H,i)*x)
    end do

    ! put KM onto surface fields for Psi_bc
    if (calculate_bcs) then
        call remap_field_to_surface(K_M, top_surface_km_values, top_surface_element_list)
        call remap_field_to_surface(K_M, bottom_surface_km_values, bottom_surface_element_list)
    end if

    !set the eddy_diffusivity and viscosity tensors for use by other fields
    call zero(eddy_diff_KH) ! zero it first as we're using an addto below
    call zero(eddy_visc_KM)

    if (on_sphere) then
      do i=1,nNodes
        eddy_diff_KH_sphere_node=align_with_radial(node_val(positions,i),node_val(K_H,i))
        eddy_visc_KM_sphere_node=align_with_radial(node_val(positions,i),node_val(K_M,i))
        call set(eddy_diff_KH,i,eddy_diff_KH_sphere_node)
        call set(eddy_visc_KM,i,eddy_visc_KM_sphere_node)
      end do
    else
      call set(eddy_diff_KH,eddy_diff_KH%dim(1),eddy_diff_KH%dim(2),K_H)
      call set(eddy_visc_KM,eddy_visc_KM%dim(1),eddy_visc_KM%dim(2),K_M)
    end if

    background_diff => extract_tensor_field(state, "GLSBackgroundDiffusivity")
    call addto(eddy_diff_KH,background_diff)
    background_visc => extract_tensor_field(state, "GLSBackgroundViscosity")
    call addto(eddy_visc_KM,background_visc)

    ewrite_minmax(K_H)
    ewrite_minmax(K_M)
    ewrite_minmax(S_H)
    ewrite_minmax(S_M)
    ewrite_minmax(ll)
    ewrite_minmax(eps)
    ewrite_minmax(tke)
    ewrite_minmax(psi) 

    ! Set viscosity
    call allocate(remaped_K_M,velocity%mesh,name="remaped_Km")
    call allocate(remaped_background_visc,velocity%mesh,name="remaped_viscosity")
    if (K_M%mesh%continuity /= viscosity%mesh%continuity) then
        ! remap
        call remap_field(K_M,remaped_K_M)
        call remap_field(background_visc,remaped_background_visc)
    else
        ! copy
        call set(remaped_K_M,K_M)
        call set(remaped_background_visc,background_visc)
    end if
    call zero(viscosity)
    if (on_sphere) then
      do i=1,nNodes
        viscosity_sphere_node=align_with_radial(node_val(positions,i),node_val(remaped_K_M,i))
        call set(viscosity,i,viscosity_sphere_node)
      end do
    else
      call set(viscosity,viscosity%dim(1),viscosity%dim(2),remaped_K_M)
    end if
    call addto(viscosity,remaped_background_visc)
    
    ! Set output on optional fields - if the field exists, stick something in it
    ! We only need to do this to those fields that we haven't pulled from state, but
    ! allocated ourselves
    call gls_output_fields(state)

    call deallocate(remaped_background_visc)
    call deallocate(remaped_K_M) 
    call deallocate(tke)

end subroutine gls_diffusivity

!----------
! gls_cleanup does...have a guess...go on.
!----------
subroutine gls_cleanup()

    ewrite(1,*) "Cleaning up GLS variables"
    ! deallocate all our variables
    if (calculate_bcs) then
        ewrite(1,*) "Cleaning up GLS surface variables"
        call deallocate(bottom_surface_values)
        call deallocate(bottom_surface_tke_values) 
        call deallocate(bottom_surface_km_values)
        call deallocate(top_surface_values)
        call deallocate(top_surface_tke_values)
        call deallocate(top_surface_km_values)
    end if
    call deallocate(NN2)
    call deallocate(MM2)  
    call deallocate(B)
    call deallocate(P)  
    call deallocate(S_H)    
    call deallocate(S_M)    
    call deallocate(K_H)    
    call deallocate(K_M)        
    call deallocate(eps) 
    call deallocate(Fwall) 
    call deallocate(density)
    call deallocate(tke_old)
    call deallocate(local_tke)
    call deallocate(ll)
    ewrite(1,*) "Finished gls_cleanup"

end subroutine gls_cleanup

!---------
! Needs to be called after an adapt to reset the fields
! and arrays within the module
! Note that clean_up has already been called in the pre-adapt hook
!----------
subroutine gls_adapt_mesh(state)

    type(state_type), intent(inout) :: state

    ewrite(1,*) "In gls_adapt_mesh"
    call gls_allocate_temps(state) ! reallocate everything
    if (calculate_bcs) then
        call gls_init_surfaces(state) ! re-do the boundaries
    end if
    call gls_init_diagnostics(state)

end subroutine gls_adapt_mesh


subroutine gls_check_options
    
    character(len=FIELD_NAME_LEN) :: buffer
    integer                       :: stat
    real                          :: min_tke, relax, nbcs
    integer                       :: dimension

    ! Don't do GLS if it's not included in the model!
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/")) return

    ! one dimensional problems not supported
    call get_option("/geometry/dimension/", dimension) 
    if (dimension .eq. 1 .and. have_option("/material_phase[0]/subgridscale_parameterisations/GLS/")) then
        FLExit("GLS modelling is only supported for dimension > 1")
    end if

    call get_option("/problem_type", buffer)
    if (.not. (buffer .eq. "oceans" .or. buffer .eq. "large_scale_ocean_options")) then
        FLExit("GLS modelling is only supported for problem type oceans or large_scale_oceans.")
    end if

    if (.not.have_option("/physical_parameters/gravity")) then
        ewrite(-1, *) "GLS modelling requires gravity" 
        FLExit("(otherwise buoyancy is a bit meaningless)")
    end if

    ! checking for required fields
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSTurbulentKineticEnergy")) then
        FLExit("You need GLSTurbulentKineticEnergy field for GLS")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSGenericSecondQuantity")) then
        FLExit("You need GLSGenericSecondQuantity field for GLS")
    end if

    ! check that the diffusivity is on for the two turbulent fields and is
    ! diagnostic
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSTurbulentKineticEnergy/prognostic/&
                          &tensor_field::Diffusivity")) then
        FLExit("You need GLSTurbulentKineticEnergy Diffusivity field for GLS")
    end if    
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSTurbulentKineticEnergy/prognostic/&
                          &tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
        FLExit("You need GLSTurbulentKineticEnergy Diffusivity field set to diagnostic/internal")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSGenericSecondQuantity/prognostic/&
                          &tensor_field::Diffusivity")) then
        FLExit("You need GLSGenericSecondQuantity Diffusivity field for GLS")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSGenericSecondQuantity/prognostic/&
                          &tensor_field::Diffusivity/diagnostic/algorithm::Internal")) then
        FLExit("You need GLSGenericSecondQuantity Diffusivity field set to diagnostic/internal")
    end if


    ! source and absorption terms...
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSTurbulentKineticEnergy/prognostic/&
                          &scalar_field::Source")) then
        FLExit("You need GLSTurbulentKineticEnergy Source field for GLS")
    end if 
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSTurbulentKineticEnergy/prognostic/&
                          &scalar_field::Source/diagnostic/algorithm::Internal")) then
        FLExit("You need GLSTurbulentKineticEnergy Source field set to diagnostic/internal")
    end if 

    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSGenericSecondQuantity/prognostic/&
                          &scalar_field::Source")) then
        FLExit("You need GLSGenericSecondQuantity Source field for GLS")
    end if      
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSGenericSecondQuantity/prognostic/&
                          &scalar_field::Source/diagnostic/algorithm::Internal")) then
        FLExit("You need GLSGenericSecondQuantity Source field set to diagnostic/internal")
    end if   

    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSTurbulentKineticEnergy/prognostic/&
                          &scalar_field::Absorption")) then
        FLExit("You need GLSTurbulentKineticEnergy Absorption field for GLS")
    end if    
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSTurbulentKineticEnergy/prognostic/&
                          &scalar_field::Absorption/diagnostic/algorithm::Internal")) then
        FLExit("You need GLSTurbulentKineticEnergy Source field set to diagnostic/internal")
    end if 
    
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSGenericSecondQuantity/prognostic/&
                          &scalar_field::Absorption")) then
        FLExit("You need GLSGenericSecondQuantity Absorption field for GLS")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &scalar_field::GLSGenericSecondQuantity/prognostic/&
                          &scalar_field::Absorption/diagnostic/algorithm::Internal")) then
        FLExit("You need GLSGenericSecondQuantity Source field set to diagnostic/internal")
    end if 


    ! background diffusivities are also needed
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &tensor_field::GLSBackgroundDiffusivity/prescribed")) then
        FLExit("You need GLSBackgroundDiffusivity tensor field for GLS")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &tensor_field::GLSBackgroundViscosity/prescribed")) then
        FLExit("You need GLSBackgroundViscosity tensor field for GLS")
    end if

    ! check for some purturbation density and velocity
    if (.not.have_option("/material_phase[0]/scalar_field::PerturbationDensity")) then
        FLExit("You need PerturbationDensity field for GLS")
    end if
    if (.not.have_option("/material_phase[0]/vector_field::Velocity")) then
        FLExit("You need Velocity field for GLS")
    end if

    ! these two fields allow the new diffusivities/viscosities to be used in
    ! other calculations - we need them!
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &tensor_field::GLSEddyViscosityKM")) then
        FLExit("You need GLSEddyViscosityKM field for GLS")
    end if
    if (.not.have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                          &tensor_field::GLSEddyViscosityKM")) then
        FLExit("You need GLSEddyViscosityKM field for GLS")
    end if

    
    ! check there's a viscosity somewhere
    if (.not.have_option("/material_phase[0]/vector_field::Velocity/prognostic/&
                          &tensor_field::Viscosity/")) then
        FLExit("Need viscosity switched on under the Velcotiy field for GLS.") 
    end if
    ! check that the user has switch Velocity/viscosity to diagnostic
    if (.not.have_option("/material_phase[0]/vector_field::Velocity/prognostic/&
                          &tensor_field::Viscosity/diagnostic/")) then
        FLExit("You need to switch the viscosity field under Velocity to diagnostic/internal")
    end if

  
    ! check a minimum value of TKE has been set
    call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                    &scalar_field::GLSTurbulentKineticEnergy/prognostic/minimum_value", min_tke, stat)
    if (stat/=0) then
        FLExit("You need to set a minimum TKE value - recommend a value of around 1e-6")
    end if

    ! check if priorities have been set - if so warn the user this might screw
    ! things up
    if (have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                     &scalar_field::GLSTurbulentKineticEnergy/prognostic/priority")) then
        ewrite(-1,*)("WARNING: Priorities for the GLS fields are set internally. Setting them in the FLML might mess things up")
    end if
    if (have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                     &scalar_field::GLSGenericSecondQuantity/prognostic/priority")) then
        ewrite(-1,*)("WARNING: Priorities for the GLS fields are set internally. Setting them in the FLML might mess things up")
    end if

    ! check the relax option is valid
    if (have_option("/material_phase[0]/subgridscale_parameterisations/GLS/relax_diffusivity")) then
        call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/relax_diffusivity", relax)
        if (relax < 0 .or. relax >= 1.0) then
            FLExit("The GLS diffusivity relaxation value should be greater than or equal to zero, but less than 1.0")
        end if
        if (.not. have_option("/material_phase[0]/subgridscale_parameterisations/GLS/scalar_field::GLSVerticalViscosity/")) then
            FLExit("You will need to switch on the GLSVerticalViscosity field when using relaxation")
        end if
        if (.not. have_option("/material_phase[0]/subgridscale_parameterisations/GLS/scalar_field::GLSVerticalDiffusivity/")) then
            FLExit("You will need to switch on the GLSVerticalDiffusivity field when using relaxation")
        end if
    end if
   
    ! Check that the we don't have auto boundaries and user-defined boundaries
    if (have_option("/material_phase[0]/subgridscale_parameterisations/GLS/calculate_boundaries")) then
        nbcs=option_count(trim("/material_phase[0]/subgridscale_parameterisations/GLS/scalar_field::GLSTurbulentKineticEnergy/prognostic/boundary_conditions"))
        if (nbcs > 0) then
            FLExit("You have automatic boundary conditions on, but some boundary conditions on the GLS TKE field. Not allowed")
        end if
        nbcs=option_count(trim("/material_phase[0]/subgridscale_parameterisations/GLS/scalar_field::GLSGenericSecondQuantity/prognostic/boundary_conditions"))
        if (nbcs > 0) then
            FLExit("You have automatic boundary conditions on, but some boundary conditions on the GLS Psi field. Not allowed")
        end if
    end if

    ! If the user has selected k-kl we need the ocean surface and bottom fields
    ! on in ocean_boundaries
    call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/option", buffer)
    if (trim(buffer) .eq. "k-kl") then
        if (.not. have_option("/geometry/ocean_boundaries")) then 
            FLExit("If you use the k-kl option under GLS, you need to switch on ocean_boundaries under /geometry/ocean_boundaries")
       end if
    end if

  end subroutine gls_check_options

!------------------------------------------------------------------!
!------------------------------------------------------------------!
!                                                                  !
!                       Private subroutines                        !
!                                                                  !
!------------------------------------------------------------------!
!------------------------------------------------------------------!

!---------
! Initilise the surface meshes used for the BCS
! Called at startup and after an adapt
!----------
subroutine gls_init_surfaces(state)
    type(state_type), intent(in)     :: state  

    type(scalar_field), pointer      :: tke
    type(vector_field), pointer      :: position
    type(mesh_type), pointer         :: ocean_mesh
    type(mesh_type)                  :: meshy
  
    ewrite(1,*) "Initialising the GLS surfaces required for BCs"

    ! grab hold of some essential field
    tke => extract_scalar_field(state, "GLSTurbulentKineticEnergy")
    position => extract_vector_field(state, "Coordinate")

    ! create a surface mesh to place values onto. This is for the top surface
    call get_boundary_condition(tke, 'tke_top_boundary', surface_mesh=ocean_mesh, &
        surface_element_list=top_surface_element_list)
    NNodes_sur = node_count(ocean_mesh) 
    call allocate(top_surface_values, ocean_mesh, name="top_surface")
    call allocate(top_surface_tke_values,ocean_mesh, name="surface_tke")
    call allocate(top_surface_km_values,ocean_mesh, name="surface_km")
    ! Creating a surface mesh gives a mapping between to global node number
    call create_surface_mesh(meshy, top_surface_nodes, tke%mesh, &
                                     top_surface_element_list, 'OceanTop')
    call deallocate(meshy)

    ! bottom
    call get_boundary_condition(tke, 'tke_bottom_boundary', surface_mesh=ocean_mesh, &
        surface_element_list=bottom_surface_element_list)
    NNodes_bot = node_count(ocean_mesh) 
    call allocate(bottom_surface_values, ocean_mesh, name="bottom_surface")
    call allocate(bottom_surface_tke_values,ocean_mesh, name="bottom_tke")
    call allocate(bottom_surface_km_values,ocean_mesh, name="bottom_km")
    call create_surface_mesh(meshy, bottom_surface_nodes, tke%mesh, &
                                     bottom_surface_element_list, 'OceanBottom')
    call deallocate(meshy)

end subroutine gls_init_surfaces

!----------------------
! Initialise the diagnostic fields, such as diffusivity, length
! scale, etc. This is called during initialisation and after an
! adapt
!----------------------
subroutine gls_init_diagnostics(state)
    type(state_type), intent(inout)     :: state

    type(scalar_field), pointer :: tke

    tke => extract_scalar_field(state, "GLSTurbulentKineticEnergy")

    ! put tke onto surface fields if we need to
    if (calculate_bcs) then
        call remap_field_to_surface(tke, top_surface_tke_values, top_surface_element_list)
        call remap_field_to_surface(tke, bottom_surface_tke_values, bottom_surface_element_list)
    end if

    call set(tke_old,tke)
    call set(FWall,1.0)

    ! bit complicated here - we need to repopulate the fields internal to this
    ! module, post adapt or at initialisation. We need the diffusivity for the first iteration to
    ! calculate the TKE src/abs terms, but for diffusivity, we need stability
    ! functions, for those we need epsilon, which is calculated in the
    ! diffusivity subroutine, but first we need the buoyancy freq.
    ! So, working backwards...
    call gls_buoyancy(state) ! buoyancy for epsilon calculation
    call gls_diffusivity(state) ! gets us epsilon, but K_H and K_M are wrong
    call gls_stability_function(state) ! requires espilon, but sets S_H and S_M
    call gls_diffusivity(state) ! sets K_H, K_M to correct values
    ! and this one sets up the diagnostic fields for output
    call gls_output_fields(state)

end subroutine gls_init_diagnostics

!----------
! Calculate the buoyancy frequency and shear velocities
!----------
subroutine gls_buoyancy(state)

    type(state_type), intent(inout)       :: state

    type(scalar_field), pointer           :: pert_rho
    type(vector_field), pointer           :: positions, gravity
    type(vector_field), pointer           :: velocity
    type(scalar_field)                    :: NU, NV, MM2_av, NN2_av, inverse_lumpedmass
    type(scalar_field), pointer           :: lumpedmass
    real                                  :: g
    logical                               :: on_sphere, smooth_buoyancy, smooth_shear
    integer                               :: ele, i, dim
    type(csr_matrix), pointer             :: mass

    ! grab variables required from state - already checked in init, so no need to check here
    positions => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "Velocity")
    pert_rho => extract_scalar_field(state, "PerturbationDensity")
    gravity => extract_vector_field(state, "GravityDirection")

    ! now allocate our temp fields
    call allocate(NU, velocity%mesh, "NU")    
    call allocate(NV, velocity%mesh, "NV")
    call set(NU, extract_scalar_field(velocity, 1))
    call set(NV, extract_scalar_field(velocity, 2))

    call get_option("/physical_parameters/gravity/magnitude", g)
    on_sphere = have_option('/geometry/spherical_earth/')
    smooth_buoyancy = have_option('/material_phase[0]/subgridscale_parameterisations/GLS/smooth_buoyancy/')
    smooth_shear = have_option('/material_phase[0]/subgridscale_parameterisations/GLS/smooth_shear/')
    dim = mesh_dim(NN2)
    
    call zero(NN2)
    call zero(MM2)
    element_loop: do ele=1, element_count(velocity)
        call assemble_elements(ele,MM2,NN2,velocity,pert_rho,NU,NV,on_sphere,dim)
    end do element_loop
  
    ! Solve
    lumpedmass => get_lumped_mass(state, NN2%mesh)
    NN2%val = NN2%val / lumpedmass%val
    lumpedmass => get_lumped_mass(state, MM2%mesh)
    MM2%val = MM2%val / lumpedmass%val

    if (smooth_shear) then
        call allocate(MM2_av, MM2%mesh, "MM2_averaged")
        call allocate(inverse_lumpedmass, MM2%mesh, "InverseLumpedMass")
        mass => get_mass_matrix(state, MM2%mesh)
        lumpedmass => get_lumped_mass(state, MM2%mesh)
        call invert(lumpedmass, inverse_lumpedmass)
        call mult( MM2_av, mass, MM2) 
        call scale(MM2_av, inverse_lumpedmass) ! so the averaging operator is [inv(ML)*M*]
        call set(MM2, MM2_av)
        call deallocate(inverse_lumpedmass)
        call deallocate(MM2_av)
    end if

    if (smooth_buoyancy) then
        call allocate(NN2_av, NN2%mesh, "NN2_averaged")
        call allocate(inverse_lumpedmass, NN2%mesh, "InverseLumpedMass")
        mass => get_mass_matrix(state, NN2%mesh)
        lumpedmass => get_lumped_mass(state, NN2%mesh)
        call invert(lumpedmass, inverse_lumpedmass)
        call mult( NN2_av, mass, NN2) 
        call scale(NN2_av, inverse_lumpedmass) ! so the averaging operator is [inv(ML)*M*]
        call set(NN2, NN2_av)
        call deallocate(NN2_av)
        call deallocate(inverse_lumpedmass)
    end if

    call deallocate(NU)
    call deallocate(NV)

    contains 
        subroutine assemble_elements(ele,MM2,NN2,velocity,rho,NU,NV,on_sphere,dim)

            type(vector_field), intent(in), pointer    :: velocity
            type(scalar_field), intent(in)             :: rho
            type(scalar_field), intent(inout)          :: NN2, MM2
            type(scalar_field), intent(in)             :: NU, NV
            logical, intent(in)                        :: on_sphere
            integer, intent(in)                        :: ele, dim

            type(element_type), pointer                :: NN2_shape, MM2_shape
            real, dimension(ele_ngi(velocity,ele))     :: detwei, shear, drho_dz
            real, dimension(dim, ele_ngi(velocity,ele)) :: grad_theta_gi, du_dz
            real, dimension(dim,ele_ngi(velocity,ele)) :: grav_at_quads
            type(element_type), pointer                :: theta_shape, velocity_shape
            integer, dimension(:), pointer             :: element_nodes
            real, dimension(ele_loc(velocity,ele),ele_ngi(velocity,ele),dim)  :: dn_t
            real, dimension(ele_loc(rho,ele),ele_ngi(rho,ele),dim)            :: dtheta_t
            real, dimension(ele_loc(velocity, ele),ele_ngi(velocity, ele),dim):: du_t
    
            NN2_shape => ele_shape(NN2, ele)
            MM2_shape => ele_shape(MM2, ele)
            velocity_shape => ele_shape(velocity, ele)
            theta_shape => ele_shape(rho, ele)

            call transform_to_physical(positions, ele, NN2_shape, &
              & dshape = dn_t, detwei = detwei)
            if(NN2_shape == velocity_shape) then
              du_t = dn_t
            else
              call transform_to_physical(positions, ele, velocity_shape, dshape = du_t)
            end if
            if(theta_shape == velocity_shape) then
              dtheta_t = dn_t
            else
              call transform_to_physical(positions, ele, theta_shape, dshape = dtheta_t)
            end if
          
            if (on_sphere) then
                grav_at_quads=radial_inward_normal_at_quad_ele(positions, ele)
            else
                grav_at_quads=ele_val_at_quad(gravity, ele)
            end if
            grad_theta_gi=ele_grad_at_quad(rho, ele, dtheta_t)
            do i=1,ele_ngi(velocity,ele)
                drho_dz(i)=dot_product(grad_theta_gi(:,i),grav_at_quads(:,i)) ! Divide this by rho_0 for non-Boussinesq?
            end do
            grad_theta_gi=ele_grad_at_quad(NU, ele, dtheta_t)
            do i=1,ele_ngi(velocity,ele)
                du_dz(1,i)=dot_product(grad_theta_gi(:,i),grav_at_quads(:,i)) ! Divide this by rho_0 for non-Boussinesq?
            end do        
            grad_theta_gi=ele_grad_at_quad(NV, ele, dtheta_t)
            do i=1,ele_ngi(velocity,ele)
                du_dz(2,i)=dot_product(grad_theta_gi(:,i),grav_at_quads(:,i)) ! Divide this by rho_0 for non-Boussinesq?
            end do
            shear = 0.0
            do i = 1, dim - 1
              shear = shear + du_dz(i,:) ** 2
            end do
              
            element_nodes => ele_nodes(NN2, ele)
            
            call addto(NN2, element_nodes, &
              ! already in the right direction due to multipling by grav_at_quads
              & shape_rhs(NN2_shape, detwei * g * drho_dz) &
              & )

            call addto(MM2, element_nodes, &
              & shape_rhs(MM2_shape,detwei * shear) &
              & )

        end subroutine assemble_elements
end subroutine gls_buoyancy

!----------
! Stability function based on Caunto et al 2001
!----------
subroutine gls_stability_function(state)

    type(state_type), intent(in)     :: state   
 
    integer                          :: i
    real                             :: N,Nt,an,anMin,anMinNum,anMinDen
    real, parameter                  :: anLimitFact = 0.5
    real                             :: d0,d1,d2,d3,d4,d5
    real                             :: n0,n1,n2,nt0,nt1,nt2
    real                             :: dCm,nCm,nCmp,cm3_inv
    real                             :: tmp0,tmp1,tmp2,tau2,as
    type(scalar_field), pointer      :: KK

    ewrite(1,*) "Calculating GLS stability functions"

    ! grab stuff from state
    KK => extract_scalar_field(state, 'GLSTurbulentKineticEnergy')
    
    ! This is written out verbatim as in GOTM v4.3.1 (also GNU GPL)
    N    =   0.5*cc1
    Nt   =   ct1
    d0   =   36.* N**3. * Nt**2.
    d1   =   84.*a5*at3 * N**2. * Nt  + 36.*at5 * N**3. * Nt
    d2   =   9.*(at2**2.-at1**2.) * N**3. - 12.*(a2**2.-3.*a3**2.) * N * Nt**2.
    d3   =   12.*a5*at3*(a2*at1-3.*a3*at2) * N + 12.*a5*at3*(a3**2.-a2**2.) * Nt       &
           + 12.*at5*(3.*a3**2.-a2**2.) * N * Nt
    d4   =   48.*a5**2.*at3**2. * N + 36.*a5*at3*at5 * N**2.
    d5   =   3.*(a2**2.-3.*a3**2.)*(at1**2.-at2**2.) * N
    n0   =   36.*a1 * N**2. * Nt**2.
    n1   = - 12.*a5*at3*(at1+at2) * N**2. + 8.*a5*at3*(6.*a1-a2-3.*a3) * N * Nt        &
           + 36.*a1*at5 * N**2. * Nt
    n2   =   9.*a1*(at2**2.-at1**2.) * N**2.
    nt0  =   12.*at3 * N**3. * Nt
    nt1  =   12.*a5*at3**2.  * N**2.
    nt2  =   9.*a1*at3*(at1-at2) * N**2. + (  6.*a1*(a2-3.*a3)                         &
           - 4.*(a2**2.-3.*a3**2.) )*at3 * N * Nt
    cm3_inv = 1./cm0**3


    ! mininum value of "an" to insure that "as" > 0 in equilibrium
    anMinNum  = -(d1 + nt0) + sqrt((d1+nt0)**2. - 4.*d0*(d4+nt1))
    anMinDen  = 2.*(d4+nt1)
    anMin     = anMinNum / anMinDen

    if (abs(n2-d5) .lt. 1e-7) then
        ! (special treatment to  avoid a singularity)
        do i=1,nNodes
            tau2   = node_val(KK,i)*node_val(KK,i) / ( node_val(eps,i)*node_val(eps,i) )
            an = tau2 * node_val(NN2,i)
            ! clip an at minimum value
            an = max(an,anLimitFact*anMin)
            ! compute the equilibrium value of as
            tmp0  = -d0 - (d1 + nt0)*an - (d4 + nt1)*an*an
            tmp1  = -d2 + n0 +  (n1-d3-nt2)*an

            as = -tmp0 / tmp1

            ! compute stability function
            dCm  = d0  +  d1*an +  d2*as + d3*an*as + d4*an*an + d5*as*as
            nCm  = n0  +  n1*an +  n2*as
            nCmp = nt0 + nt1*an + nt2*as
            call set(S_M,i, cm3_inv*nCm /dCm)
            call set(S_H,i, cm3_inv*nCmp/dCm)

        end do

     else
        do i=1,nNodes

            tau2 = node_val(KK,i)*node_val(KK,i) / ( node_val(eps,i)*node_val(eps,i) )
            an = tau2 * node_val(NN2,i)
            ! clip an at minimum value
            an = max(an,anLimitFact*anMin)

            ! compute the equilibrium value of as
            tmp0  = -d0 - (d1 + nt0)*an - (d4 + nt1)*an*an
            tmp1  = -d2 + n0 + (n1-d3-nt2)*an
            tmp2  =  n2-d5
            as =  (-tmp1 + sqrt(tmp1*tmp1-4.*tmp0*tmp2) ) / (2.*tmp2)

            ! compute stability function
            dCm  = d0  +  d1*an +  d2*as + d3*an*as + d4*an*an + d5*as*as
            nCm  = n0  +  n1*an +  n2*as
            nCmp = nt0 + nt1*an + nt2*as
            call set(S_M,i, cm3_inv*nCm /dCm)
            call set(S_H,i, cm3_inv*nCmp/dCm)
        end do

    endif

end subroutine gls_stability_function


!----------
! gls_tke_bc calculates the boundary conditions on the TKE (tke) field
! Boundary can be either Dirichlet or Neumann.
!----------
subroutine gls_tke_bc(state, bc_type)

    type(state_type), intent(in)     :: state  
    character(len=*), intent(in)     :: bc_type

    type(vector_field), pointer      :: positions 
    real                             :: gravity_magnitude
    integer                          :: i
    real, allocatable, dimension(:)  :: z0s, z0b, u_taus_squared, u_taub_squared
 
    ! grab hold of some essential field
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    positions => extract_vector_field(state, "Coordinate") 

    ! Top boundary condition
    select case(bc_type)
    case("neumann")
        ! Top TKE flux BC
        call set(top_surface_values,0.0)
        call set(bottom_surface_values,0.0)
    case("dirichlet") 
        allocate(z0s(NNodes_sur))
        allocate(z0b(NNodes_bot))
        allocate(u_taus_squared(NNodes_sur))
        allocate(u_taub_squared(NNodes_bot))
        call gls_friction(state,z0s,z0b,gravity_magnitude,u_taus_squared,u_taub_squared)

        ! Top TKE value set
        do i=1,NNodes_sur
            call set(top_surface_values,i,max(u_taus_squared(i)/(cm0**2),k_min))
        end do 
        do i=1,NNodes_bot
            call set(bottom_surface_values,i,max(u_taub_squared(i)/(cm0**2),k_min))
        end do   
        deallocate(z0s)
        deallocate(z0b)
        deallocate(u_taus_squared)
        deallocate(u_taub_squared)
    case default
        FLAbort('Unknown BC for TKE')
    end select 

end subroutine gls_tke_bc

!----------
! gls_psi_bc calculates the boundary conditions on the Psi (psi) field
! Boundary can be either Dirichlet or Neumann.
!----------
subroutine gls_psi_bc(state, bc_type)

    type(state_type), intent(in)     :: state  
    character(len=*), intent(in)     :: bc_type

    type(vector_field), pointer      :: positions 
    real                             :: gravity_magnitude
    integer                          :: i
    real, allocatable, dimension(:)  :: z0s, z0b, u_taus_squared, u_taub_squared
    type(scalar_field), pointer      :: tke, psi
    real                             :: value


    ewrite(2,*) "In gls_psi_bc: setting up"
    ! grab hold of some essential fields
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    positions => extract_vector_field(state, "Coordinate") 
    tke => extract_scalar_field(state, "GLSTurbulentKineticEnergy")    
    psi => extract_scalar_field(state, "GLSGenericSecondQuantity")

    allocate(z0s(NNodes_sur))
    allocate(z0b(NNodes_bot))
    allocate(u_taus_squared(NNodes_sur))
    allocate(u_taub_squared(NNodes_bot))

    ewrite(2,*) "In gls_psi_bc: friction"
    ! get friction
    call gls_friction(state,z0s,z0b,gravity_magnitude,u_taus_squared,u_taub_squared)
    
    ! put tke onto surface fields
    call remap_field_to_surface(tke, top_surface_tke_values, top_surface_element_list)
    call remap_field_to_surface(tke, bottom_surface_tke_values, bottom_surface_element_list)

    ewrite(2,*) "In gls_psi_bc: setting values"
    select case(bc_type)
    case("neumann")
        do i=1,NNodes_sur
            ! GOTM Boundary
            value = -(gls_n*(cm0**(gls_p+1.))*(kappa**(gls_n+1.)))/sigma_psi     &
                     *node_val(top_surface_tke_values,i)**(gls_m+0.5)*(z0s(i))**gls_n  
            ! Warner 2005 - left here for posterity and debugging
            !value = -gls_n*(cm0**(gls_p))*(node_val(top_surface_tke_values,i)**gls_m)* &
            !         (kappa**gls_n)*(z0s(i)**(gls_n-1))*((node_val(top_surface_km_values,i)/sigma_psi))
            call set(top_surface_values,i,value)
        end do
        do i=1,NNodes_bot
            if (u_taub_squared(i) < 1e-16) then
                value = 0.0
            else
                ! GOTM Boundary
                value = - gls_n*cm0**(gls_p+1.)*(kappa**(gls_n+1.)/sigma_psi)      &
                           *node_val(bottom_surface_tke_values,i)**(gls_m+0.5)*(z0b(i))**gls_n
                ! Warner 2005 - as above
                !value = gls_n*cm0**(gls_p)*node_val(bottom_surface_tke_values,i)**(gls_m)* &
                !         kappa**gls_n*(z0b(i)**(gls_n-1))*(node_val(bottom_surface_km_values,i)/sigma_psi)
            end if
            call set(bottom_surface_values,i,value)
        end do
    case("dirichlet")
        do i=1,NNodes_sur
            value = max(cm0**(gls_p-2.*gls_m)*kappa**gls_n*u_taus_squared(i)**gls_m * &
                    (z0s(i))**gls_n,psi_min)
            call set(top_surface_values,i,value)
        end do
        do i=1,NNodes_bot
            value = max(cm0**(gls_p-2.*gls_m)*kappa**gls_n*u_taub_squared(i)**gls_m * &
                    (z0b(i))**gls_n,psi_min)
            call set(bottom_surface_values,i,value)
        end do    
    case default
        FLAbort('Unknown boundary type for Psi')
    end select
  
    
    deallocate(z0s)
    deallocate(z0b)
    deallocate(u_taus_squared)
    deallocate(u_taub_squared)

end subroutine gls_psi_bc

!----------
! gls_frction works out the depth of the friction layer
! either due to bottom topography roughness or the shear stress
! on the surface
!---------
subroutine gls_friction(state,z0s,z0b,gravity_magnitude,u_taus_squared,u_taub_squared)

    type(state_type), intent(in)         :: state
    real, intent(in)                     :: gravity_magnitude
    real, dimension(:), intent(inout)    :: z0s,z0b,u_taus_squared,u_taub_squared

    integer                              :: nobcs
    integer                              :: i,ii, MaxIter
    real                                 :: rr
    real                                 :: charnock_val=18500.
    character(len=OPTION_PATH_LEN)       :: bctype
    type(vector_field), pointer          :: wind_surface_field, positions, velocity
    type(scalar_field), pointer          :: tke
    type(vector_field)                   :: bottom_velocity, surface_forcing, cont_vel, surface_pos
    type(mesh_type)                      :: ocean_mesh
    real                                 :: u_taub, z0s_min
    real, dimension(1)                   :: temp_vector_1D ! Obviously, not really a vector, but lets keep the names consistant
    real, dimension(2)                   :: temp_vector_2D
    real, dimension(3)                   :: temp_vector_3D
    logical                              :: surface_allocated
    integer                              :: stat

    MaxIter = 10
    z0s_min = 0.003
    surface_allocated = .false.
  
    ! get meshes
    velocity => extract_vector_field(state, "Velocity")
    positions => extract_vector_field(state, "Coordinate")
    tke => extract_scalar_field(state, "GLSTurbulentKineticEnergy")
    wind_surface_field => null()
   
    ! grab stresses from velocity field - Surface
    nobcs = get_boundary_condition_count(velocity)
    do i=1, nobcs
        call get_boundary_condition(velocity, i, type=bctype)
        if (bctype=='wind_forcing') then
            wind_surface_field => extract_surface_field(velocity, i, "WindSurfaceField")
            call create_surface_mesh(ocean_mesh, top_surface_nodes, tke%mesh, &
                                     top_surface_element_list, 'OceanTop')
            call allocate(surface_forcing, wind_surface_field%dim, ocean_mesh, name="surface_velocity")
            surface_pos = get_coordinates_remapped_to_surface(positions, ocean_mesh, top_surface_element_list) 
            call deallocate(ocean_mesh)

            if (tke%mesh%continuity == velocity%mesh%continuity) then
                call set(surface_forcing,wind_surface_field)
            else
                ! remap onto same mesh as TKE
                call project_field(wind_surface_field, surface_forcing, surface_pos)
            end if

            surface_allocated = .true.
            call deallocate(surface_pos)
            exit

        end if
    end do

    ! sort out bottom surface velocity
    call create_surface_mesh(ocean_mesh, bottom_surface_nodes, tke%mesh, &
                             bottom_surface_element_list, 'OceanBottom')
    call allocate(bottom_velocity, velocity%dim, ocean_mesh, name="bottom_velocity")
    call allocate(cont_vel, velocity%dim, tke%mesh, name="ContVel")
    call deallocate(ocean_mesh)
    ! Do we need to project or copy?
    if (velocity%mesh%continuity == tke%mesh%continuity) then
        call set(cont_vel,velocity)
    else            
        ! remap onto same mesh as TKE
        call project_field(velocity, cont_vel, positions)
    end if

    call remap_field_to_surface(cont_vel, bottom_velocity, &
                                bottom_surface_element_list)
    call deallocate(cont_vel)
    ! we now have a bottom velocity surface and a top surface
    ! with the wind stress on (note easier to to zero the output array
    ! below than set wind_forcing to zero and work through all the calcs

            
    ! work out the friction in either 3 or 2 dimensions.
    if (positions%dim .eq. 3) then

        if (surface_allocated) then
            do i=1,NNodes_sur
                temp_vector_2D = node_val(surface_forcing,i)
                ! big hack! Assumes that the wind stress forcing has ALREADY been divded by ocean density
                ! Note that u_taus = sqrt(wind_stress/rho0)
                ! we assume here that the wind stress in diamond is already
                ! wind_stress/rho0, hence here:
                ! u_taus = sqrt(wind_stress)
                u_taus_squared(i) = max(1e-12,sqrt(((temp_vector_2D(1))**2+(temp_vector_2D(2))**2)))
                !  use the Charnock formula to compute the surface roughness
                z0s(i)=charnock_val*u_taus_squared(i)/gravity_magnitude
                if (z0s(i).lt.z0s_min) z0s(i)=z0s_min
            end do
        else
            z0s = z0s_min
            u_taus_squared = 0.0
        end if

        do i=1,NNodes_bot
            temp_vector_3D = node_val(bottom_velocity,i)
            u_taub = sqrt(temp_vector_3D(1)**2+temp_vector_3D(2)**2+temp_vector_3D(3)**2)
            if (u_taub <= 1e-12) then
                z0b(i) = z0s_min
            else


            !  iterate bottom roughness length MaxIter times
            do ii=1,MaxIter
                    z0b(i)=(1e-7/max(1e-6,u_taub)+0.03*0.1)

                ! compute the factor r
                rr=kappa/log(z0b(i))

                ! compute the friction velocity at the bottom
                u_taub = rr*sqrt(temp_vector_3D(1)**2+temp_vector_3D(2)**2+temp_vector_3D(3)**2)

            end do
            end if

            u_taub_squared(i) = u_taub**2
        end do

    else if (positions%dim .eq. 2) then
        if (surface_allocated) then
            do i=1,NNodes_sur
                temp_vector_1D = node_val(surface_forcing,i)
                u_taus_squared(i) = max(1e-12,abs(temp_vector_1D(1)))
                !  use the Charnock formula to compute the surface roughness
                z0s(i)=charnock_val*u_taus_squared(i)/gravity_magnitude
                if (z0s(i).lt.z0s_min) z0s(i)=z0s_min

            end do
        else
            z0s = z0s_min
            u_taus_squared = 0.0
        end if

        do i=1,NNodes_bot
            temp_vector_2D = node_val(bottom_velocity,i)
            u_taub = sqrt(temp_vector_2D(1)**2+temp_vector_2D(2)**2)

            !  iterate bottom roughness length MaxIter times
            do ii=1,MaxIter
                z0b(i)=(1e-7/(max(1e-6,u_taub)+0.03*0.1))
                rr=kappa/log(z0b(i))

                ! compute the friction velocity at the bottom
                u_taub = rr*sqrt((temp_vector_2D(1)**2+temp_vector_2D(2)**2))

            end do

            u_taub_squared(i) = u_taub**2
        end do

    else
        FLAbort("Unsupported dimension in GLS friction")
    end if


    call deallocate(bottom_velocity)
    if (surface_allocated) then
        call deallocate(surface_forcing)
    end if
        
    return

end subroutine gls_friction


!---------
! Output the optional fields if they exist in state
!---------
subroutine gls_output_fields(state)

    type(state_type), intent(in)     :: state

    type(scalar_field), pointer      :: scalarField
    type(tensor_field), pointer      :: tensorField
    integer                          :: stat

    scalarField => extract_scalar_field(state, "GLSLengthScale", stat)
    if(stat == 0) then
        call set(scalarField,ll) 
    end if

    scalarField => extract_scalar_field(state,"GLSTurbulentKineticEnergyOriginal", stat)
    if(stat == 0) then
        call set(scalarField,local_tke) 
    end if  

    scalarField => extract_scalar_field(state, "GLSBuoyancyFrequency", stat)
    if(stat == 0) then
        call set(scalarField,NN2) 
    end if
    
    scalarField => extract_scalar_field(state, "GLSVelocityShear", stat)
    if(stat == 0) then
        call set(scalarField,MM2) 
    end if
    
    scalarField => extract_scalar_field(state, "GLSShearProduction", stat)
    if(stat == 0) then
        call set(scalarField,P) 
    end if
    
    scalarField => extract_scalar_field(state, "GLSBuoyancyProduction", stat)
    if(stat == 0) then
        call set(scalarField,B) 
    end if
    
    scalarField => extract_scalar_field(state, "GLSDissipationEpsilon", stat)
    if(stat == 0) then
        call set(scalarField,eps) 
    end if

    scalarField => extract_scalar_field(state, "GLSStabilityFunctionSH", stat)
    if(stat == 0) then
        call set(scalarField,S_H) 
    end if    

    scalarField => extract_scalar_field(state, "GLSStabilityFunctionSM", stat)
    if(stat == 0) then
        call set(scalarField,S_M) 
    end if           
    
    scalarField => extract_scalar_field(state, "GLSWallFunction", stat)
    if(stat == 0) then
        call set(scalarField,Fwall) 
    end if    

    scalarField => extract_scalar_field(state, "GLSVerticalViscosity", stat)
    if(stat == 0) then
        ! add vertical background
        tensorField => extract_tensor_field(state, "GLSBackgroundDiffusivity")
       call set(scalarField,K_M)  
       call addto(scalarField, extract_scalar_field(tensorField, tensorField%dim(1), tensorField%dim(2)))
    end if     
      
    scalarField => extract_scalar_field(state, "GLSVerticalDiffusivity", stat)
    if(stat == 0) then
        ! add vertical background
        tensorField => extract_tensor_field(state, "GLSBackgroundDiffusivity")
        call set(scalarField,K_H)
        call addto(scalarField, extract_scalar_field(tensorField,tensorField%dim(1), tensorField%dim(2)))
    end if  


end subroutine gls_output_fields


!---------
! Calculate the wall function as set by the user
! Each wall function has been designed with a 
! particular problem in mind, so best to have a choice here
!---------
subroutine gls_calc_wall_function(state)

    type(state_type), intent(in)     :: state

    type(scalar_field), pointer    :: distanceToBottom, distanceToTop, tke
    real                           :: LLL, distTop, distBot
    type(scalar_field)             :: top, bottom
    real, parameter                :: E2 = 1.33, E4 = 0.25
    integer                        :: i, stat

 
    ! FWall is initialised in gls_init to 1, so no need to do anything
    if (gls_wall_option .eq. "none") return

    tke => extract_scalar_field(state,"GLSTurbulentKineticEnergy")
    distanceToTop => extract_scalar_field(state, "DistanceToTop")
    distanceToBottom => extract_scalar_field(state, "DistanceToBottom") 
    call allocate(top,tke%mesh,"TopOnTKEMesh")
    call allocate(bottom,tke%mesh,"BottomOnTKEMesh")
    call remap_field(distanceToTop,top,stat)
    call remap_field(distanceToBottom,bottom,stat)
    select case (gls_wall_option)
    case ("MellorYamda")
        do i=1,nNodes
            distTop = max(1.0,node_val(top,i))
            distBot = max(1.0,node_val(bottom,i))
            LLL = (distBot +  distTop) / (distTop * distBot)
            call set( Fwall, i, 1.0 + E2*( ((node_val(ll,i)/kappa)*( LLL ))**2 ))       
        end do
    case ("Burchard98")
        do i=1,nNodes
            distTop = max(1.0,node_val(top,i))
            distBot = max(1.0,node_val(bottom,i))
            LLL = 1.0 / min(distTop,distBot)
            call set( Fwall, i, 1.0 + E2*( ((node_val(ll,i)/kappa)*( LLL ))**2 ))       
        end do
    case ("Burchard01")
        do i=1,nNodes
            distTop = max(1.0,node_val(top,i))
            distBot = max(1.0,node_val(bottom,i))
            LLL = 1.0 / distTop
            call set( Fwall, i, 1.0 + E2*( ((node_val(ll,i)/kappa)*( LLL ))**2 ))       
        end do
    case ("Blumberg")   
        do i=1,nNodes
            distTop = max(0.1,node_val(top,i))
            distBot = max(0.1,node_val(bottom,i))
            LLL = E2 * (node_val(ll,i) / (kappa *  distBot)) ** 2
            LLL = LLL + E4 * (node_val(ll,i) / (kappa *  distTop)) ** 2
            call set( Fwall, i, 1.0 + LLL)      
        end do
    case default
        FLAbort("Unknown wall function") 
    end select
    call deallocate(top)
    call deallocate(bottom)


end subroutine gls_calc_wall_function

subroutine gls_allocate_temps(state)

    type(state_type), intent(inout) :: state
    type(scalar_field), pointer    :: tkeField

    tkeField => extract_scalar_field(state,"GLSTurbulentKineticEnergy")

    ! allocate some space for the fields we need for calculations, but are optional in the model
    ! we're going to allocate these on the velocity mesh as we need one of these...
    call allocate(ll,         tkeField%mesh, "LengthScale")    
    call allocate(NN2,        tkeField%mesh, "BuoyancyFrequency")
    call allocate(MM2,        tkeField%mesh, "VelocityShear")  
    call allocate(B,          tkeField%mesh, "BuoyancyFrequency")
    call allocate(P,          tkeField%mesh, "ShearProduction")  
    call allocate(S_H,        tkeField%mesh, "StabilityH")    
    call allocate(S_M,        tkeField%mesh, "StabilityM")    
    call allocate(K_H,        tkeField%mesh, "EddyDiff")    
    call allocate(K_M,        tkeField%mesh, "EddyVisc")        
    call allocate(eps,        tkeField%mesh, "GLS_TKE_Dissipation") 
    call allocate(Fwall,      tkeField%mesh, "GLS_WallFunction") 
    call allocate(density,    tkeField%mesh, "Density")
    call allocate(tke_old,    tkeField%mesh, "Old_TKE")
    call allocate(local_tke,  tkeField%mesh, "Local_TKE")

    call set(ll,0.)
    call set(NN2,0.)
    call set(MM2,0.)
    call set(B,0.)
    call set(P,0.)
    call set(S_H,0.)
    call set(S_M,0.)
    call set(K_H,0.)
    call set(K_M,0.)
    call set(eps,0.)
    call set(density,0.)
    call set(tke_old,0.)
    call set(local_tke,tkeField)

    nNodes = node_count(tkeField)

end subroutine gls_allocate_temps

!---------
! Align the diff/visc tensors with gravity when on the sphere
!---------
function align_with_radial(position, scalar) result(rotated_tensor)
    ! Function to align viscosities/diffusivities in the radial direction when on
    ! the sphere
    real, dimension(:), intent(in) :: position
    real, intent(in) :: scalar
    real, dimension(size(position),size(position)) :: rotated_tensor
    real :: rad, phi, theta

    assert(size(position)==3)

    rad=sqrt(sum(position(:)**2))
    phi=atan2(position(2),position(1))
    theta=acos(position(3)/rad)

    rotated_tensor(1,1)=scalar*sin(theta)**2*cos(phi)**2
    rotated_tensor(1,2)=scalar*sin(theta)**2*sin(phi)*cos(phi)
    rotated_tensor(1,3)=scalar*sin(theta)*cos(theta)*cos(phi)
    rotated_tensor(2,1)=rotated_tensor(1,2)
    rotated_tensor(2,2)=scalar*sin(theta)**2*sin(phi)**2
    rotated_tensor(2,3)=scalar*sin(theta)*cos(theta)*sin(phi)
    rotated_tensor(3,1)=rotated_tensor(1,3)
    rotated_tensor(3,2)=rotated_tensor(2,3)
    rotated_tensor(3,3)=scalar*cos(theta)**2

end function align_with_radial

end module gls

