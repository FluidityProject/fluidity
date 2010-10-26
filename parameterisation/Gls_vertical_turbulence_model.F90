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
!    C.Pain@Imperial.ac.uk
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

  implicit none

  private

  ! These variables are the parameters requried by GLS. 
  ! They are all private to prevent tampering
  ! and save'd so that we don't need to call init every time some GLS-y
  ! happens.
  real, save               :: gls_n, gls_m, gls_p
  real, save               :: sigma_psi, sigma_k, kappa
  integer, save            :: nNodes
  type(scalar_field), save :: tke_old, ll
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
  type(scalar_field), save             :: top_surface_KK_values, bottom_surface_KK_values ! for the Psi BC
  integer, save                        :: NNodes_sur, NNodes_bot
  integer, dimension(:), pointer, save :: bottom_surface_nodes, top_surface_nodes
  integer, dimension(:), pointer, save :: top_surface_element_list, bottom_surface_element_list
  logical, save                        :: calculate_bcs, fix_surface_values, calc_fwall
  real, allocatable, dimension(:)      :: dzb, dzs

  ! Switch for on sphere simulations to rotate the required tensors
  logical :: on_sphere
  
  ! The following are the public subroutines
  public :: gls_init, gls_cleanup, gls_tke, gls_diffusivity, gls_psi, gls_adapt_mesh, gls_check_options

  ! General plan is:
  !  - Init in populate_state
  !  - If solve is about to do TKE, call gls_tke (which calculates NN, MM and set source/absorption for solve)
  !  - If solve is about to do Psi, call gls_psi (which fixes TKE surfaces, set source/absorption for solve)
  !  - After Psi solve, call gls_diffusivity, which sets the diffusivity and viscosity, via the lengthscale
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

    real                           :: N, gen_l, gen_alpha, gen_d, rad,rcm,cmsf
    integer                        :: stat
    character(len=FIELD_NAME_LEN)  :: gls_option, gls_stability_function
    type(scalar_field), pointer    :: s_cur
    integer                        :: tke_pri, psi_pri

    ! Allocate the temporary, module-level variables
    call gls_allocate_temps(state)

    ! Check if we're on the sphere
    on_sphere = have_option('/geometry/spherical_earth/')
   
    ! populate some useful variables
    kappa = 0.41
    calculate_bcs = have_option("/material_phase[0]/subgridscale_parameterisations/GLS/calculate_boundaries/")
    fix_surface_values = have_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                                      &calculate_boundaries/fix_surface_values")
    calc_fwall = .false.
    call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/&
                    &scalar_field::GLSTurbulentKineticEnergy/prognostic/minimum_value", k_min, stat)

    call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/relax_diffusivity", relaxation, default=0.0)

    
    call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/option", gls_option)    
    ! Check the model used - we have four choices - then set the parameters appropriately
    select case (gls_option)
    case ("k-kl")  ! Mellor-Yamada 2.5
        gls_p = 0.0
        gls_m = 1.0
        gls_n = 1.0
        sigma_k = 2.44   ! turbulent Schmidt number
        sigma_psi = 2.44  ! turbulent Schmidt number
        cPsi1 = 0.9
        cPsi2 = 0.5
        cPsi3_plus = 1.0 
        cPsi3_minus = 2.53
        psi_min = 1.e-8
        calc_fwall = .false.
        call set( Fwall, 1.0 )  
    case ("k-epsilon")
        gls_p = 3.0
        gls_m = 1.5
        gls_n = -1.0
        sigma_k = 1.3  ! turbulent Schmidt number
        sigma_psi = 1.0  ! turbulent Schmidt number
        cPsi1 = 1.44
        cPsi2 = 1.92
        cPsi3_plus = 1.0 
        cPsi3_minus = -0.52
        psi_min = 1.e-12
        call set( Fwall, 1.0 )  
    case ("k-omega")
        gls_p = -1.0
        gls_m = 0.5
        gls_n = -1.0
        sigma_k = 2.0  ! turbulent Schmidt number
        sigma_psi = 2.0  ! turbulent Schmidt number
        cPsi1 = 0.555
        cPsi2 = 0.833
        cPsi3_plus = 1.0 
        cPsi3_minus = -0.58
        psi_min = 1.e-12
        call set( Fwall, 1.0 ) 
    case ("gen")
        gls_p = 2.0
        gls_m = 1.0
        gls_n = -0.67
        sigma_k = 0.8  ! turbulent Schmidt number
        sigma_psi = 1.07  ! turbulent Schmidt number
        cPsi1 = 1.0
        cPsi2 = 1.22
        cPsi3_plus = 1.0 
        cPsi3_minus = 0.1
        psi_min = 1.e-12
        call set( Fwall, 1.0 ) 
    case default
        FLAbort("Unknown gls_option")           
    end select

    call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/stability_function", gls_stability_function)
    select case (trim(gls_stability_function))
    ! initilise the constant based on several well-known schemes
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
    kappa=cm0*sqrt(rad)
    if (gls_option .eq. "k-kl") then
        kappa = 0.41
    end if
    rcm  = cm0/cmsf
    cde = cm0**3
    gen_d = -2.*gls_n/(2.*gls_m + gls_n - 2.*cPsi2)
    gen_alpha  = -4.*gls_n*sqrt(sigma_k) / &
                 ( (1.+4.*gls_m)*sqrt(sigma_k) &
                 - sqrt(sigma_k + 24.*sigma_psi*cPsi2 ) )
    gen_l      = cm0*sqrt(rcm)* &
                 sqrt( ( (1.+4.*gls_m+8.*gls_m**2)*sigma_k &
                 + 12.*sigma_psi*cPsi2 &
                 - (1.+4.*gls_m) &
                 *sqrt(sigma_k*(sigma_k+24.*sigma_psi*cPsi2)) ) &
                 /(12.*gls_n**2.) )

    ewrite(1,*) "GLS Parameters"
    ewrite(1,*) "--------------------------------------------"
    ewrite(1,*) "gen_d", gen_d
    ewrite(1,*) "gen_alpha", gen_alpha
    ewrite(1,*) "gen_l", gen_l
    ewrite(1,*) "cm0: ",cm0
    ewrite(1,*) "kappa: ",kappa
    ewrite(1,*) "p: ",gls_p
    ewrite(1,*) "m: ",gls_m
    ewrite(1,*) "n: ",gls_n
    ewrite(1,*) "sigma_k: ",sigma_k
    ewrite(1,*) "sigma_psi: ",sigma_psi
    ewrite(1,*) "Fixing surface values: ", fix_surface_values
    ewrite(1,*) "Calculating BCs: ", calculate_bcs
    ewrite(1,*) "--------------------------------------------"
    
    ! initilise 2 GLS fields with minimum values
    s_cur => extract_scalar_field(state, "GLSTurbulentKineticEnergy")
    call set(s_cur,k_min)
    s_cur => extract_scalar_field(state, "GLSGenericSecondQuantity")
    call set(s_cur,psi_min)
    ! init other fields
    call set(NN2,0.0)
    call set(MM2,0.0)
    call set(ll,cde*k_min**1.5/psi_min)
    call set(K_H,1e-6)
    call set(K_M,1e-6)
    call set(eps,eps_min)
    if (calc_fwall) then
        call gls_calc_wall_function(state)
    end if

    ! intilise surface - only if we need to though
    if (calculate_bcs) then
        call gls_init_surfaces(state)
    end if


    ! we're all done!
end subroutine gls_init

!----------
! Update TKE
!----------
subroutine gls_tke(state)

    type(state_type), intent(inout)  :: state 
    
    type(scalar_field), pointer      :: source, absorption, kk, scalarField
    type(tensor_field), pointer      :: kk_diff, background_diff
    real                             :: prod, buoyan, diss
    integer                          :: i, stat
    character(len=FIELD_NAME_LEN)    :: bc_type
    type(scalar_field), pointer      :: scalar_surface
    type(vector_field), pointer      :: positions

    ! Temporary tensor to hold  rotated values if on the sphere (note: must be a 3x3 mat)
    real, dimension(3,3) :: K_M_sphere_node

    ewrite(1,*) "In gls_tke"

    call gls_buoyancy(state)
    ! calculate stability function
    call gls_stability_function(state)

    source => extract_scalar_field(state, "GLSTurbulentKineticEnergySource")
    absorption  => extract_scalar_field(state, "GLSTurbulentKineticEnergyAbsorption")
    kk => extract_scalar_field(state, "GLSTurbulentKineticEnergy")
    kk_diff => extract_tensor_field(state, "GLSTurbulentKineticEnergyDiffusivity")
    positions => extract_vector_field(state, "Coordinate")

    do i=1,nNodes
        prod = node_val(K_M,i)*node_val(MM2,i) 
        buoyan = -1.*node_val(K_H,i)*node_val(NN2,i)
        diss = node_val(eps,i)
        if (prod+buoyan.gt.0) then
            call set(source,i, prod + buoyan)
            call set(absorption,i, diss / node_val(KK,i))
        else
            call set(source,i, prod)
            call set(absorption,i, (diss-buoyan)/node_val(KK,i))
        end if
        call set(P, i, prod)
        call set(B, i, buoyan)
    enddo

    ! set diffusivity for KK
    call zero(kk_diff)
    background_diff => extract_tensor_field(state, "GLSBackgroundDiffusivity")
    if (on_sphere) then
      do i=1,nNodes
        K_M_sphere_node=align_with_radial(node_val(positions,i),node_val(K_M,i))
        K_M_sphere_node=K_M_sphere_node*1./sigma_k
        call set(kk_diff,i,K_M_sphere_node)
      end do
    else
      call set(kk_diff,kk_diff%dim,kk_diff%dim,K_M,scale=1./sigma_k)
    end if
    call addto(KK_diff,background_diff) 

    ! boundary conditions
    if (calculate_bcs) then
        ewrite(1,*) "Calculating BCs"
        call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/calculate_boundaries/", bc_type)
        call gls_tke_bc(state,bc_type)
        ! above puts the BC boundary values in top_surface_values and bottom_surface_values module level variables
        ! map these onto the actual BCs in kk
        scalar_surface => extract_surface_field(KK, 'tke_bottom_boundary', "value")
        call remap_field(bottom_surface_values, scalar_surface)
        scalar_surface => extract_surface_field(KK, 'tke_top_boundary', "value")
        call remap_field(top_surface_values, scalar_surface)
    end if
    
    ! finally, we need a copy of the old TKE for Psi, so grab it before we solve
    call set(tke_old,KK)

    ! that's the TKE set up ready for the solve which is the next thing to happen (see Fluids.F90)
    ewrite_minmax(source%val(:))
    ewrite_minmax(absorption%val(:))
    ! set source and absorption terms in optional output fields
    scalarField => extract_scalar_field(state, "GLSSource1", stat)
    if(stat == 0) then
       call set(scalarField,source)  
    end if
    scalarField => extract_scalar_field(state, "GLSAbsorption1", stat)
    if(stat == 0) then
       call set(scalarField,absorption)  
    end if

end subroutine gls_tke


!----------
! Calculate the second quantity
!----------
subroutine gls_psi(state)

    type(state_type), intent(inout)  :: state
    
    type(scalar_field), pointer      :: source, absorption, kk,  psi, scalarField
    type(tensor_field), pointer      :: psi_diff, background_diff
    real                             :: prod, buoyan,diss,PsiOverTke
    character(len=FIELD_NAME_LEN)    :: bc_type
    integer                          :: i, stat
    type(scalar_field), pointer      :: scalar_surface
    type(vector_field), pointer      :: positions

    ! Temporary tensor to hold  rotated values (note: must be a 3x3 mat)
    real, dimension(3,3) :: psi_sphere_node

    ewrite(1,*) "In gls_psi"

    source => extract_scalar_field(state, "GLSGenericSecondQuantitySource")
    absorption  => extract_scalar_field(state, "GLSGenericSecondQuantityAbsorption")
    kk => extract_scalar_field(state, "GLSTurbulentKineticEnergy")
    psi => extract_scalar_field(state, "GLSGenericSecondQuantity")
    psi_diff => extract_tensor_field(state, "GLSGenericSecondQuantityDiffusivity")
    positions => extract_vector_field(state, "Coordinate")

    ! clip at k_min
    do i=1,nNodes
        call set(KK,i, max(node_val(KK,i),k_min))
        call set(psi,i,max(node_val(psi,i),psi_min))
    end do

    ! re-construct psi at "old" timestep
    do i=1,nNodes
        call set(psi,i, cm0**gls_p * node_val(tke_old,i)**gls_m * node_val(ll,i)**gls_n)
    end do

    ! calc fwall if kkl
    if (calc_fwall) then
        ewrite(1,*) "Calculating the wall function for GLS"
        call gls_calc_wall_function(state)
    else
        ! this is only needed such that Fwall is not NAN after an adapt
        call set( Fwall, 1.0 )
    end if
    ! compute RHS
    do i=1,nNodes

        ! compute production terms in psi-equation
        if (node_val(B,i).gt.0) then ! note that we have already set B when setting up the RHS for TKE
            cPsi3=cPsi3_plus  ! unstable strat
        else
            cPsi3=cPsi3_minus ! stable strat
        end if

        ! compute production terms in psi-equation
        PsiOverTke  = node_val(psi,i)/node_val(tke_old,i)
        prod        = cPsi1*PsiOverTke*node_val(P,i)
        buoyan      = cPsi3*PsiOverTke*node_val(B,i)
        diss        = cPsi2*PsiOverTke*node_val(eps,i)*node_val(Fwall,i)
        if (prod+buoyan.gt.0) then
            call set(source,i, prod+buoyan)
            call set(absorption,i, diss/node_val(psi,i))
        else
            call set(source,i, prod)
            call set(absorption,i, (diss-buoyan)/node_val(psi,i))
        end if
    end do

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
      call set(psi_diff,psi_diff%dim,psi_diff%dim,K_M,scale=1./sigma_psi)
    end if
    call addto(psi_diff,background_diff) 

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


    ! Psi is now ready for solving (see Fluids.F90)
    ewrite_minmax(source%val(:))
    ewrite_minmax(absorption%val(:))
    ! set source and absorption terms in optional output fields
    scalarField => extract_scalar_field(state, "GLSSource2", stat)
    if(stat == 0) then
       call set(scalarField,source)  
    end if
    scalarField => extract_scalar_field(state, "GLSAbsorption2", stat)
    if(stat == 0) then
       call set(scalarField,absorption)  
    end if

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
    
    type(scalar_field), pointer      :: KK_state, psi_state
    type(scalar_field)               :: KK, psi, kk_copy
    type(tensor_field), pointer      :: eddy_diff_KH,eddy_visc_KM,viscosity,background_diff,background_visc
    real                             :: exp1, exp2, exp3, x
    integer                          :: i, stat
    real                             :: epslim, tke
    real, parameter                  :: galp = 0.748331 ! sqrt(0.56)
    logical                          :: limit_length = .true.
    type(vector_field), pointer      :: positions

    ! Temporary tensors to hold  rotated values (note: must be a 3x3 mat)
    real, dimension(3,3) :: eddy_diff_KH_sphere_node, eddy_visc_KM_sphere_node, viscosity_sphere_node

    ewrite(1,*) "In gls_diffusivity"

    KK_state => extract_scalar_field(state, "GLSTurbulentKineticEnergy")
    psi_state => extract_scalar_field(state, "GLSGenericSecondQuantity")
    eddy_visc_KM  => extract_tensor_field(state, "GLSEddyViscosityKM",stat)
    eddy_diff_KH => extract_tensor_field(state, "GLSEddyDiffusivityKH",stat)
    viscosity => extract_tensor_field(state, "Viscosity",stat)
    positions => extract_vector_field(state, "Coordinate")

    call allocate(KK,KK_state%mesh,"TKE")
    call allocate(psi,psi_state%mesh,"PSI")
    call allocate(kk_copy,kk_state%mesh,"TKE_COPY")
    call set(kk,kk_state)
    call set(psi,psi_state)
    
    if (fix_surface_values) then
    
        ! call the bc code, but specify we want dirichlet
        call gls_tke_bc(state, 'dirichlet')
     
        ! copy the values onto the mesh using the global node id
        do i=1,NNodes_sur
            call set(kk,top_surface_nodes(i),node_val(top_surface_values,i))
        end do   

        call set(kk_copy,kk_state)
        call set(kk_state,kk)

        ! call the bc code, but specify we want dirichlet
        call gls_psi_bc(state, 'dirichlet')

        ! copy the values onto the mesh using the global node id
        do i=1,NNodes_sur
            call set(psi,top_surface_nodes(i),node_val(top_surface_values,i))
        end do   

        call set(kk_state,kk_copy)

    end if


    exp1 = 3.0 + gls_p/gls_n
    exp2 = 1.5 + gls_m/gls_n
    exp3 =       - 1.0/gls_n
   
    do i=1,nNodes

        tke = node_val(KK,i)

        ! recover dissipation rate from k and psi
        call set(eps,i, cm0**exp1 * tke**exp2 * node_val(psi,i)**exp3)
        
        ! clip at eps_min
        call set(eps,i, max(node_val(eps,i),eps_min))

        ! compute dissipative scale
        call set(ll,i,cde*sqrt(tke**3.)/node_val(eps,i))
    end do

    ! limit dissipation rate under stable stratification,
    ! see Galperin et al. (1988)
    if (limit_length) then
        do i=1,nNodes
            if (node_val(NN2,i) .gt. 0) then
                ! compute limit (the sqrt(2) comes from the fact that we've
                ! defined k = tke, but the limit uses q where tke = 1/2*(q^2)
                epslim = (cde*node_val(KK,i)*sqrt(node_val(NN2,i))) / (sqrt(2.)*galp)

                ! clip at new limit
                call set(eps,i, max(node_val(eps,i),epslim))

                ! re-compute dissipative scale
                call set(ll,i, cde*sqrt(node_val(kk,i)**3)/node_val(eps,i))
            end if

        end do
    endif
    

    ! calculate diffusivities for next step and for use in other fields
    do i=1,nNodes
        x = sqrt(node_val(KK,i))*node_val(ll,i)
        ! momentum
        call set(K_M,i, min(7.0,relaxation*node_val(K_M,i) + (1-relaxation)*node_val(S_M,i)*x))
        ! tracer
        call set(K_H,i, min(7.0,relaxation*node_val(K_H,i) + (1-relaxation)*node_val(S_H,i)*x))
    end do


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
      call set(eddy_diff_KH,eddy_diff_KH%dim,eddy_diff_KH%dim,K_H)
      call set(eddy_visc_KM,eddy_visc_KM%dim,eddy_visc_KM%dim,K_M)
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
    ewrite_minmax(kk)
    ewrite_minmax(psi) 

    ! Set viscosity
    call zero(viscosity)
    if (on_sphere) then
      do i=1,nNodes
        viscosity_sphere_node=align_with_radial(node_val(positions,i),node_val(K_M,i))
        call set(viscosity,i,viscosity_sphere_node)
      end do
    else
      call set(viscosity,viscosity%dim,viscosity%dim,K_M)
    end if
    call addto(viscosity,background_visc)
    
    ! Set output on optional fields - if the field exists, stick something in it
    ! We only need to do this to those fields that we haven't pulled from state, but
    ! allocated ourselves
    call gls_output_fields(state)

    call deallocate(kk)
    call deallocate(psi)
    call deallocate(kk_copy)
    
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
        call deallocate(bottom_surface_kk_values)    
        call deallocate(top_surface_values)
        call deallocate(top_surface_kk_values)
        deallocate(dzb)
        deallocate(dzs)
    end if
    call deallocate(ll)    
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
    ewrite(1,*) "Finished gls_cleanup"

end subroutine gls_cleanup

!---------
! Needs to be called after an adapt to reset the fields
! and arrays within the module
!----------
subroutine gls_adapt_mesh(state)

    type(state_type), intent(inout) :: state

    ewrite(1,*) "In gls_adapt_mesh"
    call gls_allocate_temps(state) ! reallocate everything
    if (calculate_bcs) then
        call gls_init_surfaces(state) ! re-do the boundaries
    end if

    ! bit complicated here - we need to repopulate the fields internal to this
    ! module, post adapt. We need the diffusivity for the first iteration to
    ! calculate the TKE src/abs terms, but for diffusivity, we need stability
    ! functions, for those we need epsilon, which is calculated in the
    ! diffusivity subroutine, but first we need the buoyancy freq.
    ! So, working backwards...
    call gls_buoyancy(state) ! buoyancy for epsilon calculation
    call gls_diffusivity(state) ! gets us epsilon, but K_H and K_M are wrong
    call gls_stability_function(state) ! requires espilon, but sets S_H and S_M
    call gls_diffusivity(state) ! sets K_H, K_M to correct values
    ! calling these two, sets the source/abs correctly
    call gls_tke(state)
    call gls_psi(state)
    ! and this one sets up the diagnostic fields again
    call gls_output_fields(state)

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
    if (buffer/="oceans") then
        FLExit("GLS modelling is only supported for problem type oceans.")
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


    ! If the user has selected kkl we need the ocean surface and bottom fields
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

    type(scalar_field), pointer      :: kk
    type(vector_field), pointer      :: position
    type(mesh_type), pointer         :: ocean_mesh
  
    ewrite(1,*) "Initialising the GLS surfaces required for BCs"

    ! grab hold of some essential field
    kk => extract_scalar_field(state, "GLSTurbulentKineticEnergy")
    position => extract_vector_field(state, "Coordinate")

    ! create a surface mesh to place values onto. This is for the top surface
    call get_boundary_condition(kk, 'tke_top_boundary', surface_mesh=ocean_mesh, &
        surface_element_list=top_surface_element_list)
    NNodes_sur = node_count(ocean_mesh) 
    call allocate(top_surface_values, ocean_mesh, name="top_surface")
    call allocate(top_surface_kk_values,ocean_mesh, name="surface_tke")
    
    ! bottom
    call get_boundary_condition(kk, 'tke_bottom_boundary', surface_mesh=ocean_mesh, &
        surface_element_list=bottom_surface_element_list)
    NNodes_bot = node_count(ocean_mesh) 
    call allocate(bottom_surface_values, ocean_mesh, name="bottom_surface")
    call allocate(bottom_surface_kk_values,ocean_mesh, name="bottom_tke")

    allocate(dzb(NNodes_bot))
    allocate(dzs(NNodes_sur))

    call gls_calculate_dz(state)

end subroutine gls_init_surfaces

!----------
! Calculate the buoyancy frequency and shear velocities
!----------
subroutine gls_buoyancy(state)

    type(state_type), intent(inout)       :: state
    type(scalar_field), pointer           :: pert_rho
    type(vector_field), pointer           :: positions, gravity
    type(vector_field), pointer           :: velocity
    type(scalar_field)                    :: pert_rho_averaged,inverse_lumpedmass,NU_averaged,NV_averaged
    type(scalar_field), pointer           :: lumpedmass
    type(csr_matrix), pointer             :: mass
    real                                  :: g
    logical                               :: on_sphere
    integer                               :: ele, i, dim
    integer, dimension(:), pointer        :: element_nodes
    type(element_type), pointer           :: NN2_shape, MM2_shape
    real, allocatable, dimension(:)       :: detwei, shear, drho_dz
    real, allocatable, dimension(:,:)     :: grad_theta_gi, grav_at_quads, du_dz
    real, allocatable, dimension(:,:,:)   :: dn_t
    real, allocatable, dimension(:,:,:)   :: dtheta_t
    real, allocatable, dimension(:,:,:)   :: du_t
    type(element_type), pointer           :: theta_shape, velocity_shape


    ! grab variables required from state - already checked in init, so no need to check here
    positions => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "Velocity")
    pert_rho => extract_scalar_field(state, "PerturbationDensity")
    gravity => extract_vector_field(state, "GravityDirection")

    ! now allocate our temp fields
    call allocate(pert_rho_averaged, velocity%mesh, "pert_rho_averaged")   
    call allocate(NU_averaged, velocity%mesh, "NU_averaged")    
    call allocate(NV_averaged, velocity%mesh, "NV_averaged")


    ! Small smoothing filter to iron out any big wiggles that might 
    ! be spurious unstable strat
    call allocate(inverse_lumpedmass, velocity%mesh, "InverseLumpedMass")
    mass => get_mass_matrix(state, velocity%mesh)
    lumpedmass => get_lumped_mass(state, velocity%mesh)
    call invert(lumpedmass, inverse_lumpedmass)
    call mult( pert_rho_averaged, mass, pert_rho)
    call scale(pert_rho_averaged, inverse_lumpedmass) ! so the averaging operator is [inv(ML)*M*]
    call mult( NU_averaged, mass, extract_scalar_field(velocity, 1) )
    call scale(NU_averaged, inverse_lumpedmass) ! so the averaging operator is [inv(ML)*M*]
    call mult( NV_averaged, mass, extract_scalar_field(velocity, 2) )
    call scale(NV_averaged, inverse_lumpedmass) ! so the averaging operator is [inv(ML)*M*]
    call get_option("/physical_parameters/gravity/magnitude", g)
    on_sphere = have_option('/geometry/spherical_earth/')
    dim = mesh_dim(NN2)
    NN2_shape => ele_shape(NN2, ele)
    MM2_shape => ele_shape(MM2, ele)
    velocity_shape => ele_shape(velocity, ele)
    theta_shape => ele_shape(pert_rho_averaged, ele)
    
    call zero(NN2)
    call zero(MM2)
    element_loop: do ele=1, element_count(velocity)

        allocate(grad_theta_gi(ele_ngi(velocity, ele), dim))
        allocate(grav_at_quads(dim,ele_ngi(velocity, ele)))
        allocate(dn_t(ele_loc(velocity, ele), ele_ngi(velocity, ele), dim))
        allocate(dtheta_t(ele_loc(pert_rho_averaged, ele), ele_ngi(pert_rho_averaged, ele), dim))
        allocate(du_t(ele_loc(velocity, ele), ele_ngi(velocity, ele), dim))
        allocate(detwei(ele_ngi(velocity, ele)))
        allocate(shear(ele_ngi(velocity, ele)))
        allocate(drho_dz(ele_ngi(velocity, ele)))
        allocate(du_dz(ele_ngi(velocity, ele),dim))

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
            grav_at_quads=sphere_inward_normal_at_quad_ele(positions, ele)
        else
            grav_at_quads=ele_val_at_quad(gravity, ele)
        end if
        grad_theta_gi=ele_grad_at_quad(pert_rho_averaged, ele, dtheta_t)
        do i=1,ele_ngi(velocity,ele)
            drho_dz(i)=dot_product(grad_theta_gi(i,:),grav_at_quads(:,i)) ! Divide this by rho_0 for non-Boussinesq?
        end do
        grad_theta_gi=ele_grad_at_quad(NU_averaged, ele, dtheta_t)
        do i=1,ele_ngi(velocity,ele)
            du_dz(i,1)=dot_product(grad_theta_gi(i,:),grav_at_quads(:,i)) ! Divide this by rho_0 for non-Boussinesq?
        end do        
        grad_theta_gi=ele_grad_at_quad(NV_averaged, ele, dtheta_t)
        do i=1,ele_ngi(velocity,ele)
            du_dz(i,2)=dot_product(grad_theta_gi(i,:),grav_at_quads(:,i)) ! Divide this by rho_0 for non-Boussinesq?
        end do
        shear = 0.0
        do i = 1, dim - 1
          shear = shear + du_dz(:,i) ** 2
        end do
          
        element_nodes => ele_nodes(NN2, ele)
        
        !call addto(NN2, element_nodes, &
        !  & shape_rhs(NN2_shape, -detwei * g * grad_theta_gi(:,dim)*grav_at_quads(:,:)) &
        !  & )
        call addto(NN2, element_nodes, &
          ! already in the right direction due to multipling by grav_at_quads
          & shape_rhs(NN2_shape, detwei * g * drho_dz) &
          & )


        call addto(MM2, element_nodes, &
          & shape_rhs(MM2_shape,detwei * shear) &
          & )

        deallocate(grad_theta_gi)
        deallocate(grav_at_quads)
        deallocate(dn_t)
        deallocate(dtheta_t)
        deallocate(du_t)
        deallocate(detwei)
        deallocate(shear)
        deallocate(drho_dz)
        deallocate(du_dz)

    end do element_loop
  
    ! Solve
    NN2%val = NN2%val / lumpedmass%val
    MM2%val = MM2%val / lumpedmass%val

 
    call deallocate(pert_rho_averaged)
    call deallocate(inverse_lumpedmass)
    call deallocate(NU_averaged)    
    call deallocate(NV_averaged)

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
    
    ! This is copied verbatim from GOTM v4.3.1 AUTHORS
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

            tau2   = node_val(KK,i)*node_val(KK,i) / ( node_val(eps,i)*node_val(eps,i) )
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
! gls_tke_bc calculates the boundary conditions on the TKE (kk) field
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

    allocate(z0s(NNodes_sur))
    allocate(z0b(NNodes_bot))
    allocate(u_taus_squared(NNodes_sur))
    allocate(u_taub_squared(NNodes_bot))

    ! Top boundary condition
    select case(bc_type)
    case("neumann")
        ! Top TKE flux BC
        do i=1,NNodes_sur
            call set(top_surface_values,i,0.0)
        end do 
        do i=1,NNodes_bot
            call set(bottom_surface_values,i,0.0)
        end do
    case("dirichlet") 
        call gls_friction(state,z0s,z0b,gravity_magnitude,u_taus_squared,u_taub_squared) 
        ! Top TKE value set
        do i=1,NNodes_sur
            call set(top_surface_values,i,u_taus_squared(i)/cm0**2)
        end do 
        do i=1,NNodes_bot
            call set(bottom_surface_values,i,u_taub_squared(i)/cm0**2)
        end do
    case default
        FLAbort('Unknown BC for TKE')
    end select 
    
    deallocate(z0s)
    deallocate(z0b)
    deallocate(u_taus_squared)
    deallocate(u_taub_squared)

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
    type(scalar_field), pointer      :: kk, kkOld
    real                             :: value
 
    ! grab hold of some essential field
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
    positions => extract_vector_field(state, "Coordinate") 
    kk => extract_scalar_field(state, "GLSTurbulentKineticEnergy")    
    kkOld => extract_scalar_field(state, "OldGLSTurbulentKineticEnergy")

    allocate(z0s(NNodes_sur))
    allocate(z0b(NNodes_bot))
    allocate(u_taus_squared(NNodes_sur))
    allocate(u_taub_squared(NNodes_bot))

    ! get friction
    call gls_friction(state,z0s,z0b,gravity_magnitude,u_taus_squared,u_taub_squared)
    
    ! put KK onto surface fields
    call remap_field_to_surface(kk, top_surface_kk_values, top_surface_element_list)
    call remap_field_to_surface(kk, bottom_surface_kk_values, bottom_surface_element_list)

    select case(bc_type)
    case("neumann")
        if (fix_surface_values) then
            do i=1,NNodes_sur
                value = (-gls_n*cm0**(gls_p+1.)*kappa**(gls_n+1.))/sigma_psi      &
                           *(u_taus_squared(i)/(cm0**2))**(gls_m+0.5)*(z0s(i)+z0s(i))**gls_n
                call set(top_surface_values,i,value)
            end do
            do i=1,NNodes_bot
                value = - gls_n*cm0**(gls_p+1.)*kappa**(gls_n+1.)/sigma_psi      &
                           *(u_taub_squared(i)/(cm0**2))**(gls_m+0.5)*(z0b(i)+z0b(i))**gls_n
                call set(bottom_surface_values,i,value)
            end do
        else
            do i=1,NNodes_sur
                value = (-gls_n*cm0**(gls_p+1.)*kappa**(gls_n+1.))/sigma_psi      &
                           *node_val(top_surface_kk_values,i)**(gls_m+0.5)*(z0s(i)+z0s(i))**gls_n
                call set(top_surface_values,i,value)
            end do
            do i=1,NNodes_bot
                value = - gls_n*cm0**(gls_p+1.)*kappa**(gls_n+1.)/sigma_psi      &
                           *node_val(bottom_surface_kk_values,i)**(gls_m+0.5)*(z0b(i)+z0b(i))**gls_n
                call set(bottom_surface_values,i,value)
            end do
        end if
    case("dirichlet")
        if (fix_surface_values) then
            do i=1,NNodes_sur
                value = cm0**gls_p*kappa**gls_n*(u_taus_squared(i)/(cm0**2))**gls_m * &
                    (z0s(i)+z0s(i))**gls_n
                call set(top_surface_values,i,value)
            end do
            do i=1,NNodes_bot
                value = cm0**gls_p*kappa**gls_n*(u_taub_squared(i)/(cm0**2))**gls_m * &
                    (z0b(i)+z0b(i))**gls_n
                call set(bottom_surface_values,i,value)
            end do            
        else
            do i=1,NNodes_sur
                value = cm0**gls_p*kappa**gls_n*(node_val(top_surface_kk_values,i))**gls_m * &
                    (z0s(i)+z0s(i))**gls_n
                call set(top_surface_values,i,value)
            end do
            do i=1,NNodes_bot
                value = cm0**gls_p*kappa**gls_n*(node_val(bottom_surface_kk_values,i))**gls_m * &
                    (z0b(i)+z0b(i))**gls_n
                call set(bottom_surface_values,i,value)
            end do
        end if
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
    real                                 :: charnock_val=1400.
    character(len=OPTION_PATH_LEN)       :: bctype
    type(vector_field), pointer          :: wind_surface_field, positions, velocity
    type(vector_field)                   :: bottom_velocity
    type(mesh_type)                      :: ocean_mesh
    real                                 :: u_taub, z0s_min
    real, dimension(1)                   :: temp_vector_1D ! Obviously, not really a vector, but lets keep the names consistant
    real, dimension(2)                   :: temp_vector_2D
    real, dimension(3)                   :: temp_vector_3D

    MaxIter = 10
    z0s_min = 0.003
   
    ! get meshes
    velocity => extract_vector_field(state, "Velocity")
    positions => extract_vector_field(state, "Coordinate")
    wind_surface_field => null()
   
    ! grab stresses from velocity field - Surface
    nobcs = get_boundary_condition_count(velocity)
    do i=1, nobcs
        call get_boundary_condition(velocity, i, type=bctype)
        if (bctype=='wind_forcing') then
            wind_surface_field => extract_surface_field(velocity, i, "WindSurfaceField")
        end if
    end do
            
    if (positions%dim .eq. 3) then

        if (associated(wind_surface_field)) then
            do i=1,NNodes_sur
                temp_vector_2D = node_val(wind_surface_field,i)
                ! big hack! Assumes that the wind stress forcing has ALREADY been divded by ocean density
                ! u_taus = sqrt(wind_stress/rho0)
                ! we assume here that the wind stress in diamond is already
                ! wind_stress/rho0
                u_taus_squared(i) = sqrt(((temp_vector_2D(1))**2+(temp_vector_2D(2))**2))
                !  use the Charnock formula to compute the surface roughness
                z0s(i)=charnock_val*u_taus_squared(i)/gravity_magnitude
                if (z0s(i).lt.z0s_min) z0s(i)=z0s_min
            end do
        else
            z0s = 0.0
        end if

        ! grab values of velocity from bottom surface
        call create_surface_mesh(ocean_mesh, bottom_surface_nodes, velocity%mesh, bottom_surface_element_list, 'OceanBottom')
        call allocate(bottom_velocity, velocity%dim, ocean_mesh, name="bottom_velocity")
        call remap_field_to_surface(velocity, bottom_velocity, &
                                    bottom_surface_element_list)

        do i=1,NNodes_bot
            temp_vector_3D = node_val(bottom_velocity,i)
            u_taub = sqrt(temp_vector_3D(1)**2+temp_vector_3D(2)**2+temp_vector_3D(3)**2)


            !  iterate bottom roughness length MaxItz0b times
            do ii=1,MaxIter
                z0b(i)=1e-7/max(1e-6,u_taub)+0.03*0.1

                ! compute the factor r (version 1, with log-law)
                ! Note that we do this in the closest 1m to surface - irrespective of grid size
                rr=kappa/(log((z0b(i)+dzb(i))/z0b(i)))

                ! compute the friction velocity at the bottom
                u_taub = rr*sqrt(temp_vector_3D(1)**2+temp_vector_3D(2)**2+temp_vector_3D(3)**2)

            end do

            u_taub_squared(i) = u_taub**2
        end do

    else if (positions%dim .eq. 2) then
        if (associated(wind_surface_field)) then
            do i=1,NNodes_sur
                temp_vector_1D = node_val(wind_surface_field,i)
                ! big hack! Assumes that the wind stress forcing has ALREADY been divded by ocean density
                ! u_taus = sqrt(wind_stress/rho0)
                ! we assume here that the wind stress in diamond is already
                ! wind_stress/rho0
                u_taus_squared(i) = temp_vector_1D(1)
                !  use the Charnock formula to compute the surface roughness
                z0s(i)=charnock_val*u_taus_squared(i)/gravity_magnitude
                if (z0s(i).lt.z0s_min) z0s(i)=z0s_min

            end do
        else
            z0s = 0.0
        end if

        ! grab values of velocity from bottom surface
        call create_surface_mesh(ocean_mesh, bottom_surface_nodes, velocity%mesh, bottom_surface_element_list, 'OceanBottom')
        call allocate(bottom_velocity, velocity%dim, ocean_mesh, name="bottom_velocity")
        call remap_field_to_surface(velocity, bottom_velocity, &
                                    bottom_surface_element_list)

        do i=1,NNodes_bot
            temp_vector_2D = node_val(bottom_velocity,i)
            u_taub = sqrt(temp_vector_2D(1)**2+temp_vector_2D(2)**2)


            !  iterate bottom roughness length MaxItz0b times
            do ii=1,MaxIter
                z0b(i)=1e-7/max(1e-6,u_taub)+0.03*0.1

                ! compute the factor r (version 1, with log-law)
                ! Note that we do this in the closest 1m to surface - irrespective of grid size
                rr=kappa/(log((z0b(i)+dzb(i))/z0b(i)))

                ! compute the friction velocity at the bottom
                u_taub = rr*sqrt((temp_vector_2D(1)**2+temp_vector_2D(2)**2))

            end do

            u_taub_squared(i) = u_taub**2
        end do

    else
        FLAbort("Unsupported dimension in GLS friction")
    end if


    call deallocate(bottom_velocity)
    call deallocate(ocean_mesh)
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
        call set(scalarField,FWall) 
    end if           

    scalarField => extract_scalar_field(state, "GLSVerticalViscosity", stat)
    if(stat == 0) then
        ! add vertical background
        tensorField => extract_tensor_field(state, "GLSBackgroundDiffusivity")
       call set(scalarField,K_M)  
       call addto(scalarField, extract_scalar_field(tensorField, tensorField%dim, tensorField%dim))
    end if     
      
    scalarField => extract_scalar_field(state, "GLSVerticalDiffusivity", stat)
    if(stat == 0) then
        ! add vertical background
        tensorField => extract_tensor_field(state, "GLSBackgroundDiffusivity")
        call set(scalarField,K_H)
        call addto(scalarField, extract_scalar_field(tensorField,tensorField%dim, tensorField%dim))
    end if  


end subroutine gls_output_fields

!----------
! Calculates an approximation to "dz" in finite difference code
!----------
subroutine gls_calculate_dz(state)
    
    type(state_type), intent(in)     :: state

    type(vector_field), pointer      :: positions, velocity
    integer                          :: i, j, nele, ele, sele
    real                             :: node_dz, h
    type(patch_type)                 :: current_patch
    type(mesh_type)                  :: sur_mesh

    positions => extract_vector_field(state, "Coordinate")
    velocity => extract_vector_field(state, "Velocity")

    ! for each node in the surface, get the elements and find the dz
    ! average out for this (i.e. add and divide by nodes per element)
    ! Top surface
    call create_surface_mesh(sur_mesh, top_surface_nodes, velocity%mesh, top_surface_element_list, 'Top') 

    do i=1,NNodes_sur
        current_patch = get_patch_ele(sur_mesh,i)
        nele = current_patch%count
        node_dz = 100000.
        do j=1,nele
            sele = top_surface_element_list(current_patch%elements(j))
            ele = face_ele(positions, sele)
            h = gls_get_normal_element_size(ele, sele, positions, velocity)
            if (h < node_dz) then
                node_dz = h
            end if
        end do
        dzs(i) = node_dz
        deallocate(current_patch%elements)
    end do

    call deallocate(sur_mesh)
 
    ! Bottom surface
    call create_surface_mesh(sur_mesh, bottom_surface_nodes, velocity%mesh, bottom_surface_element_list, 'Bottom') 

    do i=1,NNodes_bot
        current_patch = get_patch_ele(sur_mesh,i)
        nele = current_patch%count
        node_dz = 100000.
        do j=1,nele
            sele = bottom_surface_element_list(current_patch%elements(j))
            ele = face_ele(positions, sele)
            h = gls_get_normal_element_size(ele, sele, positions, velocity)
            if (h < node_dz) then
                node_dz = h
            end if
        end do
        dzb(i) = node_dz
        deallocate(current_patch%elements)
    end do

    call deallocate(sur_mesh)
    
end subroutine gls_calculate_dz


subroutine gls_calc_wall_function(state)

    type(state_type), intent(in)     :: state

    type(scalar_field), pointer    :: distanceToBottom, distanceToTop
    real                           :: E2, LLL
    integer                        :: i
    
    E2 = 1.33

    distanceToTop => extract_scalar_field(state, "DistanceToTop")
    distanceToBottom => extract_scalar_field(state, "DistanceToBottom")  
    do i=1,nNodes ! there are lots of alternative formulae for this wall function
        !LLL = (node_val(distanceToTop,i) + node_val(distanceToBottom,i)) / &
        !      (node_val(distanceToTop,i) * node_val(distanceToBottom,i)) 
        LLL = 1.0 / node_val(distanceToBottom,i)
        !if( (node_val(distanceToBottom,i).lt.1.0) .or.  (node_val(distanceToTop,i).lt.1.0) ) then
        if( (node_val(distanceToBottom,i).lt.1.0)) then
            call set( Fwall, i, 1.0 + E2 ) ! hanert-ish       
        else
            call set( Fwall, i, 1.0 + E2*( ((node_val(ll,i)/kappa)*( LLL ))**2 ))       
        end if
    end do


end subroutine gls_calc_wall_function

subroutine gls_allocate_temps(state)

    type(state_type), intent(inout) :: state
    type(vector_field), pointer    :: vectorField

    vectorField => extract_vector_field(state,"Velocity")

    ! allocate some space for the fields we need for calculations, but are optional in the model
    ! we're going to allocate these on the velocity mesh as we need one of these...
    call allocate(ll,         vectorField%mesh, "LengthScale")    
    call allocate(NN2,        vectorField%mesh, "BuoyancyFrequency")
    call allocate(MM2,        vectorField%mesh, "VelocityShear")  
    call allocate(B,          vectorField%mesh, "BuoyancyFrequency")
    call allocate(P,          vectorField%mesh, "ShearProduction")  
    call allocate(S_H,        vectorField%mesh, "StabilityH")    
    call allocate(S_M,        vectorField%mesh, "StabilityM")    
    call allocate(K_H,        vectorField%mesh, "EddyDiff")    
    call allocate(K_M,        vectorField%mesh, "EddyVisc")        
    call allocate(eps,        vectorField%mesh, "GLS_TKE_Dissipation") 
    call allocate(Fwall,      vectorField%mesh, "GLS_WallFunction") 
    call allocate(density,    vectorField%mesh, "Density")
    call allocate(tke_old,    vectorField%mesh, "Old_TKE")

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
    call set(FWall,0.)
    call set(density,0.)
    call set(tke_old,0.)

    nNodes = node_count(vectorField)

end subroutine gls_allocate_temps

function gls_get_normal_element_size(ele, sele, x, u) &
                result (h)

    integer              :: ele, sele
    type(vector_field)   :: x, u

    real, dimension(1,1)                            :: hb
    integer                                         :: ndim, snloc
    integer, dimension(face_loc(u, sele))           :: u_nodes_bdy
    real, dimension(x%dim,x%dim)                    :: G
    real, dimension(x%dim,1)                        :: n
    real                                            :: h
    real, dimension(x%dim, x%dim, ele_ngi(u, sele)) :: invJ
    real, dimension(x%dim, face_ngi(u, sele))       :: normal_bdy
    real, dimension(face_ngi(u, sele))              :: detwei_bdy
    integer, dimension(:), pointer                  :: ele_nodes_u

    ndim        = x%dim
    snloc       = face_loc(u, sele)
    u_nodes_bdy = face_global_nodes(u, sele)
    ele_nodes_u => ele_nodes(u, ele)

    call compute_inverse_jacobian( ele_val(x, ele), &
            ele_shape(x, ele), invJ )

    call transform_facet_to_physical( x, sele, detwei_f=detwei_bdy, normal=normal_bdy )
    
    ! calculate wall-normal element mesh size
    G = matmul(transpose(invJ(:,:,1)), invJ(:,:,1))
    n(:,1) = normal_bdy(:,1)
    hb = 1. / sqrt( matmul(matmul(transpose(n), G), n) )
    h  = hb(1,1) 
end function gls_get_normal_element_size


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

