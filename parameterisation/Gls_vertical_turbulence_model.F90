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
!
module gls
  use quadrature
  use elements
  use field_derivatives
  use fields
  use sparse_matrices_fields
  use state_module
  use spud
  use allsorts
  use global_parameters, only:   OPTION_PATH_LEN
  use equation_of_state
  use state_fields_module
  use boundary_conditions
  use FLDebug
!
  implicit none
!
  private
!  
  public :: gls_vertical_turbulence_model
!
contains
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine gls_vertical_turbulence_model(state) 
  
  ! The generic length scale turbulence closure model
  ! implementation based upon
  ! Ocean Modelling 8 (2005) 81-113, by J.C. Warner et al,
  ! "Performance of four turbulence closure models implemented using a 
  ! generic lenght scale model"
    type(state_type), intent(inout) :: state

    type(vector_field), pointer :: XX,NU,gravity_direction
    type(tensor_field), pointer :: viscosity,kk_diff,psi_diff
    type(scalar_field), pointer :: kk,psi,pert_rho,topdis,botdis,source1,source2,absorption1,absorption2
    type(scalar_field), pointer :: buoyancy_freq,vel_shear,stab_SH,stab_SM,src1,src2,abs1,abs2,length_scale
    type(scalar_field), pointer :: shear_production,buoyancy_production,GLSGH,GLSFwall,dissipation_epsilon,GLSVertVisc,GLSVertDiff
    type(tensor_field), pointer :: eddy_visc_KM,eddy_diff_KH,background_diff,background_visc
    type(scalar_field) :: ll,MM2,NN2,epsilon,Fwall,c3,S_H,S_M,K_H,K_M,density,PP,BB
    type(scalar_field), dimension(1) :: du_dz,dv_dz,drho_dz
    type(scalar_field) :: bottom_kk, surface_kk
    character(len=OPTION_PATH_LEN) :: gls_option, gls_stability_function, gls_bc_type
    real :: gravity_magnitude
    real :: gls_p,gls_m,gls_n,gls_sigma_k,gls_sigma_psi,k_min,psi_min
    real :: c,c1,c2,c3plus,c3minus,c_mu_zero,cc,rho_0
    integer :: nonods,ii,stat
    logical :: output_buoyancy_freq,output_vel_shear,output_stab_SH,output_stab_SM,output_src1,output_src2,output_abs1,output_abs2
    logical :: output_shear_production,output_buoyancy_production,output_GLSGH,output_GLSFwall,output_GLSVertVisc,output_GLSVertDiff
    logical :: output_length_scale,output_dissipation_epsilon
    logical :: got_background_visc,got_background_diff
    real :: value
    type(scalar_field) :: boundary_condition_surface, boundary_condition_bottom
    type(scalar_field) :: surface_dz, bottom_dz
    type(scalar_field), pointer :: scalar_surface
    type(scalar_field), pointer :: distanceToTop, distanceToBottom
    integer :: i, NNodes_sur, NNodes_bot
    type(mesh_type) :: ocean_mesh, input_mesh
    type(mesh_type), pointer :: surface_mesh
    character(len=FIELD_NAME_LEN) input_mesh_name
    integer, dimension(:), pointer :: surface_element_list, surface_nodes
    real :: kappa 
    type(vector_field), pointer :: positions
    real, dimension(:,:,:), allocatable :: J
    real :: node_dz
    integer :: ele, nele, gi
    type(patch_type) :: current_patch
    type(scalar_field) :: dz
    real, allocatable, dimension(:) :: z0s, z0b, u_taus_squared, u_taub_squared
    real, allocatable, dimension(:,:) :: X_val
    type(element_type), pointer :: X_shape


    ewrite(1,*)'Now in gls_turbulence_model - MDP'

 
    XX => extract_vector_field(state, "Coordinate")
    NU => extract_vector_field(state, "Velocity")

    kk  => extract_scalar_field(state, "GLSTurbulentKineticEnergy",stat)
    if(stat/=0) FLAbort("Need GLSTurbulentKineticEnergy field")

    kk_diff  => extract_tensor_field(state, "GLSTurbulentKineticEnergyDiffusivity",stat)
    if(stat/=0) FLAbort("Need GLSTurbulentKineticEnergyDiffusivity field")

    psi => extract_scalar_field(state, "GLSGenericSecondQuantity",stat)
    if(stat/=0) FLAbort("Need GLSGenericSecondQuantity field")

    psi_diff => extract_tensor_field(state, "GLSGenericSecondQuantityDiffusivity",stat)
    if(stat/=0) FLAbort("Need GLSGenericSecondQuantityDiffusivity field")

    source1  => extract_scalar_field(state, "GLSTurbulentKineticEnergySource",stat)
    if(stat/=0) FLAbort("Need source for GLSTurbulentKineticEnergySource field")

    source2  => extract_scalar_field(state, "GLSGenericSecondQuantitySource",stat)
    if(stat/=0) FLAbort("Need source for GLSGenericSecondQuantitySource field")
    
    absorption1  => extract_scalar_field(state, "GLSTurbulentKineticEnergyAbsorption",stat)
    if(stat/=0) FLAbort("Need source for GLSTurbulentKineticEnergyAbsorption field")

    absorption2  => extract_scalar_field(state, "GLSGenericSecondQuantityAbsorption",stat)
    if(stat/=0) FLAbort("Need source for GLSGenericSecondQuantityAbsorption field")


!    temperature_diffusivity => extract_tensor_field(state, "TemperatureDiffusivity",stat)
!    if(stat/=0) FLAbort("Need diffusivity for Temperature field")

!    salinity_diffusivity => extract_tensor_field(state, "SalinityDiffusivity",stat)
!    if(stat/=0) FLAbort("Need diffusivity for Salinity field")
    
    viscosity => extract_tensor_field(state, "Viscosity",stat)
    if(stat/=0) FLAbort("Need viscosity")
        
    gravity_direction => extract_vector_field(state,"GravityDirection",stat)
    if(stat/=0) FLAbort("Need gravity")
    
    pert_rho => extract_scalar_field(state, "PerturbationDensity",stat)
    if(stat/=0) FLAbort("Need PerturbationDensity")

    eddy_visc_KM  => extract_tensor_field(state, "GLSEddyViscosityKM",stat)
    if(stat/=0) FLAbort("Need GLSEddyViscosityKM")
    
    eddy_diff_KH => extract_tensor_field(state, "GLSEddyDiffusivityKH",stat)
    if(stat/=0) FLAbort("Need GLSEddyDiffusivityKH")
    

    length_scale => extract_scalar_field(state, "GLSLengthScale", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSLengthScale field'
      output_length_scale = .true.
    else
      output_length_scale = .false.
    end if
    buoyancy_freq => extract_scalar_field(state, "GLSBuoyancyFrequency", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSBuoyancyFrequency field'
      output_buoyancy_freq = .true.
    else
      output_buoyancy_freq = .false.
    end if
    vel_shear => extract_scalar_field(state, "GLSVelocityShear", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSVelocityShear field'
      output_vel_shear = .true.
    else
      output_vel_shear = .false.
    end if
    shear_production => extract_scalar_field(state, "GLSShearProduction", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSShearProduction field'
      output_shear_production = .true.
    else
      output_shear_production = .false.
    end if
    buoyancy_production => extract_scalar_field(state, "GLSBuoyancyProduction", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSBuoyancyProduction field'
      output_buoyancy_production = .true.
    else
      output_buoyancy_production = .false.
    end if
    dissipation_epsilon => extract_scalar_field(state, "GLSDissipationEpsilon", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSDissipationEpsilon field'
      output_dissipation_epsilon = .true.
    else
      output_dissipation_epsilon = .false.
    end if


    stab_SH => extract_scalar_field(state, "GLSStabilityFunctionSH", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSStabilityFunctionSH field'
      output_stab_SH = .true.
    else
      output_stab_SH = .false.
    end if    
    stab_SM => extract_scalar_field(state, "GLSStabilityFunctionSM", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSStabilityFunctionSM field'
      output_stab_SM = .true.
    else
      output_stab_SM = .false.
    end if        


    src1 => extract_scalar_field(state, "GLSSource1", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSSource1 field'
      output_src1 = .true.
    else
      output_src1 = .false.
    end if        
    src2 => extract_scalar_field(state, "GLSSource2", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSSource2 field'
      output_src2 = .true.
    else
      output_src2 = .false.
    end if        
    abs1 => extract_scalar_field(state, "GLSAbsorption1", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSAbsorption1 field'
      output_abs1 = .true.
    else
      output_abs1 = .false.
    end if        
    abs2 => extract_scalar_field(state, "GLAbsorption2", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSAbsorption2 field'
      output_abs2 = .true.
    else
      output_abs2 = .false.
    end if        


    GLSGH => extract_scalar_field(state, "GLSGH", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSGH field'
      output_GLSGH = .true.
    else
      output_GLSGH = .false.
    end if      
    GLSFwall => extract_scalar_field(state, "GLSWallFunction", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSWallFunction field'
      output_GLSFwall = .true.
    else
      output_GLSFwall = .false.
    end if           


    background_diff => extract_tensor_field(state, "GLSBackgroundDiffusivity",stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSBackgroundDiffusivity field'
      got_background_diff = .true.
    else
      got_background_diff = .false.
    end if      
    background_visc => extract_tensor_field(state, "GLSBackgroundViscosity",stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSBackgroundViscosity field'
      got_background_visc = .true.
    else
      got_background_visc = .false.
    end if          



    GLSVertVisc => extract_scalar_field(state, "GLSVerticalViscosity", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSVerticalViscosity field'
      output_GLSVertVisc = .true.
    else
      output_GLSVertVisc = .false.
    end if           
    GLSVertDiff => extract_scalar_field(state, "GLSVerticalDiffusivity", stat)
    if(stat == 0) then
      ewrite(3,*)'Found GLSVerticalDiffusivity field'
      output_GLSVertDiff = .true.
    else
      output_GLSVertDiff = .false.
    end if           


        
    nonods = node_count(NU)


! assume for the moment that everything is on the velocity mesh until I figure out the correct thing to do
    call allocate(ll, NU%mesh, "LengthScale")    
    call allocate(NN2, NU%mesh, "BuoyancyFrequency")
    call allocate(MM2, NU%mesh, "VelocityShear")  
    call allocate(BB, NU%mesh, "BuoyancyFrequency")
    call allocate(PP, NU%mesh, "ShearProduction")  
    call allocate(S_H, NU%mesh, "StabilityH")    
    call allocate(S_M, NU%mesh, "StabilityM")    
    call allocate(K_H, NU%mesh, "EddyDiff")    
    call allocate(K_M, NU%mesh, "EddyVisc")        
    call allocate(epsilon, NU%mesh, "GLS_TKE_Dissipation") 
    call allocate(Fwall, NU%mesh, "GLS_WallFunction") 
    call allocate(c3, NU%mesh, "GLS_c3") 
    call allocate(du_dz(1), NU%mesh, "DuDz")
    call allocate(dv_dz(1), NU%mesh, "DvDz")
    call allocate(drho_dz(1), NU%mesh, "DRhoDz")
    call allocate(density,NU%mesh, "Density")

   
    kappa = 0.4
    call calculate_perturbation_density(state, density)
    ewrite_minmax(density%val(:))

    call zero(ll)
    call zero(epsilon)
    
    ! Calculate the bouyancy frequency and the velocity shear in the vertical
    call calc_vel_shear_and_buoyancy_freq(state,MM2,NN2,NU,XX,density,du_dz,dv_dz,drho_dz)
    ewrite_minmax(MM2%val(:))
    ewrite_minmax(NN2%val(:))


    ! which model are we using?  - again, assume first material_phase
    call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/option", gls_option)
    ! which stability function option are we using?
    call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/stability_function", gls_stability_function)
    ! which boundary conditions are set?
    if (have_option("/material_phase[0]/subgridscale_parameterisations/GLS/calculate_boundaries/")) then
        call get_option("/material_phase[0]/subgridscale_parameterisations/GLS/calculate_boundaries/", gls_bc_type)
    end if

    ! fix negative TKS values - for the second time in here for the second equation    
    do ii = 1,node_count(kk)
      call set(kk,ii,max(0.0,node_val(kk,ii)))
    end do
 
    ! calculate some gls parameters and the wall proximity function for k-kl (MY25)
    call gls_parameters_and_wall_proximity(gls_option,gls_stability_function,topdis,botdis,epsilon,kk,psi,ll,NN2,MM2, &
                                gls_p,gls_m,gls_n,gls_sigma_k,gls_sigma_psi,k_min,psi_min,Fwall, &
                                c1,c2,c3plus,c3minus,c_mu_zero,cc,kappa)
   
    ! apply some hard limiting 
    call gls_limiting(ll,kk,psi,epsilon,NN2,k_min,psi_min,gls_n,gls_p,gls_m,c_mu_zero)        

    ! calculate the stability functions
    call gls_stability_functions(gls_stability_function,MM2,NN2,ll,kk,epsilon,S_H,S_M,c_mu_zero,GLSGH,output_GLSGH,XX)

    ! apply some hard limiting again
    call gls_limiting(ll,kk,psi,epsilon,NN2,k_min,psi_min,gls_n,gls_p,gls_m,c_mu_zero)

    !caclulate the buoyancy parameter c3
    call gls_buoyancy_parameter(c3,c3plus,c3minus,NN2)

    ! calculate GLS source terms, viscosities etc
    call gls_source_absor_visc_diff(source1,source2,absorption1,absorption2,epsilon,psi,cc,c1,c2,c3,MM2,NN2,kk,ll,S_M,S_H,K_M,K_H,Fwall, &
                                        gls_sigma_k,gls_sigma_psi,PP,BB, &
                                        viscosity,kk_diff,psi_diff,eddy_visc_KM,eddy_diff_KH)

    !set the eddy_diffusivity and viscosoty tensors for use by other fields
    call zero(eddy_diff_KH) ! zero it first as we're using an addto below
    call zero(eddy_visc_KM)
    
    call set(eddy_diff_KH,eddy_diff_KH%dim,eddy_diff_KH%dim,K_H)     
    call set(eddy_visc_KM,eddy_visc_KM%dim,eddy_visc_KM%dim,K_M) 

    if(got_background_diff) then
      call addto(eddy_diff_KH,background_diff)    
    endif
    if(got_background_visc) then
      call addto(eddy_visc_KM,background_visc)    
    endif

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! set boundary conditions !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (have_option('/material_phase[0]/subgridscale_parameterisations/GLS/calculate_boundaries')) then
      distanceToTop => extract_scalar_field(state, "DistanceToTop")
      distanceToBottom => extract_scalar_field(state, "DistanceToBottom")
      call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)
      positions => extract_vector_field(state, "Coordinate")

      ! create field for dz
      call allocate(dz, positions%mesh, name="mesh_spacing_z")

      ! surface
      call get_option(trim(kk%option_path)//'/prognostic/mesh/name', input_mesh_name)
      input_mesh = extract_mesh(state, input_mesh_name);
      call get_boundary_condition(distanceToTop, name='top', surface_element_list=surface_element_list)
      call create_surface_mesh(ocean_mesh, surface_nodes, input_mesh, surface_element_list, 'OceanSurface')
      NNodes_sur = node_count(ocean_mesh) 
      call allocate(boundary_condition_surface, ocean_mesh, name="boundary_condition_surface")
      call allocate(surface_kk,ocean_mesh, name="surface_tke")
      call allocate(surface_dz, ocean_mesh, name="surface_dz")
      call remap_field_to_surface(kk, surface_kk, surface_element_list)
      call deallocate(ocean_mesh)

      ! bottom
      call get_boundary_condition(distanceToBottom, name='bottom', surface_element_list=surface_element_list)
      call create_surface_mesh(ocean_mesh, surface_nodes, input_mesh, surface_element_list, 'OceanBottom')
      NNodes_bot = node_count(ocean_mesh) 
      call allocate(boundary_condition_bottom, ocean_mesh, name="boundary_condition_bottom")
      call allocate(bottom_kk,ocean_mesh, name="bottom_tke")
      call allocate(bottom_dz, ocean_mesh, name="bottom_dz")
      call remap_field_to_surface(kk, bottom_kk, surface_element_list)
      call deallocate(ocean_mesh)

      allocate(z0s(NNodes_sur))
      allocate(z0b(NNodes_bot))
      allocate(u_taus_squared(NNodes_sur))
      allocate(u_taub_squared(NNodes_bot))

      ! calculate "dz" - put into zb and zs
      ! we only need an estimate, so this method will do for now
      ! I don't think it's the best method though...      
      ! loop over nodes in position mesh
      ! For each element that contains that node, calculate
      ! the Jacobian of positions, which will contain
      ! dx, dy and dz. Then average the estimates for each element
      ! and place on the node.
      ! This method is non-optimal as we calculate the same thing several
      ! times...
      do ii=1,positions%mesh%nodes
        current_patch = get_patch_ele(positions%mesh,ii)
        nele = current_patch%count
        node_dz = 0
        do ele=1,nele
          allocate(X_val(positions%dim,ele_loc(positions,current_patch%elements(ele))))
          X_val=ele_val(positions,current_patch%elements(ele))
          X_shape => ele_shape(positions,current_patch%elements(ele))
          allocate(J(size(X_val,1), size(X_val,1), X_shape%ngi))
          do gi=1,X_shape%ngi
            J(:,:,gi)=matmul(X_val(:,:), X_shape%dn(:, gi, :))
          end do
          ! dx is an abs non-zero item from column 1
          ! dy is an abs non-zero item from column 2
          ! dz is an abs non zero item from column 3
          ! We grab the first non-zero number...
          do i=1,3
            if (J(3,1,i) /= 0) then
              node_dz = node_dz + abs(J(3,1,i))
              exit
            end if
          end do
          deallocate(X_val)
          deallocate(J)
        end do
        node_dz = node_dz / current_patch%count
        node_dz = node_dz / 2
        call set(dz,ii,node_dz)
        deallocate(current_patch%elements)
      end do

      call get_boundary_condition(distanceToTop, name='top', surface_element_list=surface_element_list)
      call remap_field_to_surface(dz, surface_dz, &
                                surface_element_list)
      call get_boundary_condition(distanceToBottom, name='bottom', surface_element_list=surface_element_list)
      call remap_field_to_surface(dz, bottom_dz, &
                                surface_element_list)


      call deallocate(dz)

      ! get frictions, z0s and z0b
      call friction(kappa,state,z0s,z0b,gravity_magnitude,bottom_dz,u_taus_squared,u_taub_squared)

      select case(gls_bc_type)
      case("neumann")
        ! Top TKE flux BC
        do i=1,NNodes_sur
            call set(boundary_condition_surface,i,0.0)
        end do 
    case("dirichlet")  
        ! Top TKE value set
        do i=1,NNodes_sur
            call set(boundary_condition_surface,i,u_taus_squared(i)/c_mu_zero**2)
        end do 
    case default
        FLAbort('Unknown surface BC for TKE')
    end select
      scalar_surface => extract_surface_field(kk, 'tke_top_boundary', "value")
      call remap_field(boundary_condition_surface, scalar_surface)

      select case(gls_bc_type)
      case("neumann")

        select case(gls_option)
        case ("k-kl","gen") ! Mellor-Yamada 2.5 & general
            do i=1,NNodes_sur
                value = - gls_n*c_mu_zero**(gls_p+1.)*0.41**(gls_n+1.)/gls_sigma_psi      &
                       *node_val(surface_kk,i)**(gls_m+0.5)*(node_val(surface_dz,i)+z0s(i))**gls_n
                call set(boundary_condition_surface,i,value)
            end do
        case ("k-epsilon","k-omega")  
            do i=1,NNodes_sur
                value = ((c_mu_zero**4.)*(node_val(surface_kk,i))**2.)/(gls_sigma_psi*(node_val(surface_dz,i)+z0s(i)))
                call set(boundary_condition_surface,i,value)
            end do 
        case default
            FLAbort("Unknown gls_option")           
        end select

    case("dirichlet")
        select case(gls_option)
        case ("k-kl","gen") ! Mellor-Yamada 2.5 & general
            do i=1,NNodes_sur
                value = c_mu_zero**gls_p*kappa**gls_n*(node_val(surface_kk,i))**gls_m * &
                        (node_val(surface_dz,i)+z0s(i))**gls_n
                call set(boundary_condition_surface,i,value)
            end do
        case ("k-epsilon","k-omega")  
            do i=1,NNodes_sur
                value = c_mu_zero**3*(node_val(surface_kk,i))**1.5/(kappa*(node_val(surface_dz,i)+z0s(i)))
                call set(boundary_condition_surface,i,value)
            end do 
        case default
            FLAbort("Unknown gls_option")           
        end select
    
        case default
        FLAbort('Unknown surface BC for Psi')
    end select
            scalar_surface => extract_surface_field(psi, 'psi_top_boundary', "value")
        call remap_field(boundary_condition_surface, scalar_surface)
        call deallocate(boundary_condition_surface)

    select case(gls_bc_type)
    case("neumann")  
        do i=1,NNodes_bot
            call set(boundary_condition_bottom,i,0.0)
        end do 
    case("dirichlet")
      do i=1,NNodes_bot
        call set(boundary_condition_bottom,i,u_taub_squared(i)/c_mu_zero**2)
      end do 
  case default
        FLAbort('Unknown bottom BC for TKE')
    end select
      scalar_surface => extract_surface_field(kk, 'tke_bottom_boundary', "value")
      call remap_field(boundary_condition_bottom, scalar_surface)

    select case(gls_bc_type)
    case("neumann")
      select case(gls_option)
      case ("k-kl","gen") ! Mellor-Yamada 2.5 & general
        do i=1,NNodes_bot
            value = - gls_n*c_mu_zero**(gls_p+1.)*0.41**(gls_n+1.)/gls_sigma_psi      &
                       *node_val(bottom_kk,i)**(gls_m+0.5)*(node_val(bottom_dz,i)+z0b(i))**gls_n
            call set(boundary_condition_bottom,i,value)
        end do
      case ("k-epsilon","k-omega")  
        do i=1,NNodes_bot
            value = ((c_mu_zero**4.)*(node_val(bottom_kk,i))**2.)/(gls_sigma_psi*(node_val(bottom_dz,i)+z0b(i)))
            call set(boundary_condition_bottom,i,value)
        end do 
      case default
        FLAbort("Unknown gls_option")           
      end select
      case("dirichlet")
        select case(gls_option)
        case ("k-kl","gen") ! Mellor-Yamada 2.5 & general
            do i=1,NNodes_sur
                value = c_mu_zero**gls_p*kappa**gls_n*(node_val(bottom_kk,i))**gls_m * &
                        (node_val(bottom_dz,i)+z0s(i))**gls_n
                call set(boundary_condition_bottom,i,value)
            end do
         case ("k-epsilon","k-omega")  
            do i=1,NNodes_bot
                value = c_mu_zero**3*(node_val(surface_kk,i))**1.5/(kappa*(node_val(bottom_dz,i)+z0s(i)))
                call set(boundary_condition_bottom,i,value)
            end do 
        case default
            FLAbort("Unknown gls_option")           
        end select
    
    case default
        FLAbort('Unknown bottom BC for Psi')
    end select
      
      scalar_surface => extract_surface_field(psi, 'psi_bottom_boundary', "value")
      call remap_field(boundary_condition_bottom, scalar_surface)
      call deallocate(boundary_condition_bottom)
      
      deallocate(z0s)
      deallocate(z0b)
      deallocate(u_taus_squared)
      deallocate(u_taub_squared)
      call deallocate(bottom_kk)
      call deallocate(surface_kk)
      call deallocate(surface_dz)
      call deallocate(bottom_dz)

    end if

    ewrite_minmax(source1%val(:))
    ewrite_minmax(source2%val(:))
    ewrite_minmax(absorption1%val(:))
    ewrite_minmax(absorption2%val(:))
    ewrite_minmax(kk%val(:))
    ewrite_minmax(psi%val(:))
    ewrite_minmax(epsilon%val(:))   
    ewrite_minmax(ll%val(:))
    ewrite_minmax(Fwall%val(:))
    ewrite_minmax(S_H%val(:))
    ewrite_minmax(S_M%val(:))
    ewrite_minmax(K_H%val(:))
    ewrite_minmax(K_M%val(:))

    
    
    if(output_length_scale) then
      call set(length_scale,ll)
    end if

    if(output_buoyancy_freq) then
      call set(buoyancy_freq,NN2)
    end if

    if(output_vel_shear) then
      call set(vel_shear,MM2)
    end if

    if(output_stab_SH) then
      call set(stab_SH,S_H)
    end if

    if(output_stab_SM) then
      call set(stab_SM,S_M)
    end if    

    if(output_src1) then
      call set(src1,source1)
    end if  

    if(output_src2) then
      call set(src2,source2)
    end if     

    if(output_abs1) then
      call set(abs1,absorption1)
    end if  

    if(output_abs2) then
      call set(abs2,absorption2)
    end if  
    
    if(output_buoyancy_production) then
      call set(buoyancy_production,BB)
    end if    

    if(output_shear_production) then
      call set(shear_production,PP)
    end if     

    if(output_dissipation_epsilon) then
      call set(dissipation_epsilon,epsilon)
    end if     

    if(output_GLSFwall) then
      call set(GLSFwall,Fwall)
    end if     

    if(output_GLSVertVisc) then
      call set(GLSVertVisc,K_M)
    end if     

    if(output_GLSVertDiff) then
      call set(GLSVertDiff,K_H)
    end if       
!           
    call deallocate(ll)
    call deallocate(MM2)
    call deallocate(NN2)
    call deallocate(BB)
    call deallocate(PP)
    call deallocate(S_H)
    call deallocate(S_M)
    call deallocate(K_H)
    call deallocate(K_M)
    call deallocate(epsilon)
    call deallocate(Fwall)
    call deallocate(c3)   
    call deallocate(du_dz(1))
    call deallocate(dv_dz(1))
    call deallocate(drho_dz(1))
    call deallocate(density)

!         
  end subroutine gls_vertical_turbulence_model
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine calc_vel_shear_and_buoyancy_freq(state,MM2,NN2,NU,XX,pert_rho,du_dz,dv_dz,drho_dz)
    type(state_type), intent(inout) :: state
    type(scalar_field), intent(inout) :: MM2,NN2
    type(scalar_field), dimension(1), intent(inout) :: du_dz,dv_dz,drho_dz
    type(scalar_field), intent(in) :: pert_rho
    type(vector_field), intent(in) :: XX
    type(vector_field), intent(inout) :: NU
    type(scalar_field) :: NU_averaged,NV_averaged,pert_rho_averaged,inverse_lumpedmass
    type(scalar_field), pointer :: lumpedmass
    type(csr_matrix), pointer :: mass
    real :: gravity_magnitude
    integer :: ii
    logical :: average
  

! Should we do some averaging of velocities and density first?

    call allocate(NU_averaged, NU%mesh, "NU_averaged")    
    call allocate(NV_averaged, NU%mesh, "NV_averaged")   
    call allocate(pert_rho_averaged, NU%mesh, "pert_rho_averaged")   

    average = .false.
    if (average) then
      call allocate(inverse_lumpedmass, NU%mesh, "InverseLumpedMass")
      
      mass => get_mass_matrix(state, NU%mesh)
      lumpedmass => get_lumped_mass(state, NU%mesh)
      call invert(lumpedmass, inverse_lumpedmass)


      call mult( pert_rho_averaged, mass, pert_rho )
      call scale(pert_rho_averaged, inverse_lumpedmass) ! so the averaging operator is [inv(ML)*M*]
        
      call mult( NU_averaged, mass, extract_scalar_field(NU, 1) )
      call scale(NU_averaged, inverse_lumpedmass) ! so the averaging operator is [inv(ML)*M*]

      call mult( NV_averaged, mass, extract_scalar_field(NU, 2) )
      call scale(NV_averaged, inverse_lumpedmass) ! so the averaging operator is [inv(ML)*M*]

    else
      call set( NU_averaged, extract_scalar_field(NU, 1) )     
      call set( NV_averaged, extract_scalar_field(NU, 2) )     
      call set( pert_rho_averaged, pert_rho )     
    endif
    
    if(NU%dim==2) then
      call differentiate_field( NU_averaged, XX, (/.false., .true./), du_dz )
!!      call differentiate_field( NV_averaged, XX, (/.false., .true./), dv_dz )
      call differentiate_field( pert_rho_averaged, XX, (/.false., .true./), drho_dz )
    else
      call differentiate_field( NU_averaged, XX, (/.false., .false., .true./), du_dz )
      call differentiate_field( NV_averaged, XX, (/.false., .false., .true./), dv_dz )
      call differentiate_field( pert_rho_averaged, XX, (/.false., .false., .true./), drho_dz )    
    endif
    call get_option("/physical_parameters/gravity/magnitude", gravity_magnitude)

    do ii = 1, node_count(MM2)
      if(NU%dim==2) then
        call set( MM2, ii, (node_val(du_dz(1),ii))**2  ) ! velocity shear (squared)
      else
        call set( MM2, ii, (node_val(du_dz(1),ii))**2 + (node_val(dv_dz(1),ii))**2 ) ! velocity shear (squared)
      end if    
      call set( NN2, ii, max(-gravity_magnitude * node_val(drho_dz(1),ii),0.0)  )    ! bouyancy frequency (squared)
    end do      
    
    call deallocate(NU_averaged)    
    call deallocate(NV_averaged)   
    call deallocate(pert_rho_averaged)   
    if (average) then
      call deallocate(inverse_lumpedmass)
    endif
!
  end subroutine calc_vel_shear_and_buoyancy_freq
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine gls_parameters_and_wall_proximity(gls_option,gls_stability_function,topdis,botdis,epsilon,kk,psi,ll,NN2,MM2, &
                                gls_p,gls_m,gls_n,gls_sigma_k,gls_sigma_psi,k_min,psi_min,Fwall, &
                                c1,c2,c3plus,c3minus,c_mu_zero,cc,kappa)
                                               
    ! Warner et al 2005, Table 1
    character(len=OPTION_PATH_LEN), intent(in) :: gls_option,gls_stability_function
    type(scalar_field), intent(in)  :: topdis,botdis,NN2,MM2
    type(scalar_field), intent(inout) :: kk,psi,epsilon,ll,Fwall
    real, intent(out) :: gls_p,gls_m,gls_n,gls_sigma_k,gls_sigma_psi,k_min,psi_min
    real, intent(out) :: c1,c2,c3plus,c3minus,c_mu_zero,cc
    real, intent(in) :: kappa
    real :: E2, LLL
    integer :: ii, nonods
        
    nonods = node_count(Fwall)       
  
    select case (gls_option)
    case ("k-kl")  ! Mellor-Yamada 2.5
      gls_p = 0.0
      gls_m = 1.0
      gls_n = 1.0
      gls_sigma_k = 2.44   ! turbulent Schmidt number
      gls_sigma_psi = 2.44  ! turbulent Schmidt number
      c1 = 0.9
      c2 = 0.5
      c3plus = 1.0 
      k_min = 5.e-6
      psi_min = 1.e-8
      E2 = 1.33 ! 
      do ii = 1, nonods   ! there are lots of alternative formulae for this wall function
        LLL = 1./max(1.,node_val(topdis,ii)) + 1./max(1.,node_val(botdis,ii))
        LLL = min(1./(0.0001*(node_val(botdis,ii) + node_val(topdis,ii))),LLL)

        if( (node_val(botdis,ii).lt.1.0) .or.  (node_val(topdis,ii).lt.1.0) ) then
          call set( Fwall, ii, 1.0 + E2 ) ! hanert-ish       
        else
          call set( Fwall, ii, 1.0 + E2*( ((node_val(ll,ii)/kappa)*( LLL ))**2 ))       
        end if
      end do

      if((trim(gls_stability_function)=="KanthaClayson-94").or.(trim(gls_stability_function)=="Galperin-88")) then      
        c3minus = 2.53  ! Warner et al 2005, Table 2
        c_mu_zero = 0.5544      
      else if(gls_stability_function=="Canuto-01-A") then
        c3minus = 2.38  ! Warner et al 2005, Table 2
        c_mu_zero = 0.5270   
      else if(gls_stability_function=="Canuto-01-B") then 
        !c3minus =  not_defined  ! Warner et al 2005, Table 2          
        c_mu_zero = 0.5540   
      end if
      cc = 1.0  ! Warner et al 2005, top of page 85       
    case ("k-epsilon")
      gls_p = 3.0
      gls_m = 1.5
      gls_n = -1.0
      gls_sigma_k = 1.0  ! turbulent Schmidt number
      gls_sigma_psi = 1.3  ! turbulent Schmidt number
      c1 = 1.44
      c2 = 1.92
      c3plus = 1.0 
      k_min = 7.6e-6
      psi_min = 1.e-12
      call set( Fwall, 1.0 )  
      if((trim(gls_stability_function)=="KanthaClayson-94").or.(trim(gls_stability_function)=="Galperin-88")) then      
        c3minus = -0.52  ! Warner et al 2005, Table 2
        c_mu_zero = 0.5544      
      else if(trim(gls_stability_function)=="Canuto-01-A") then
        c3minus = -0.63  ! Warner et al 2005, Table 2
        c_mu_zero = 0.5270   
      else if(trim(gls_stability_function)=="Canuto-01-B") then 
        c3minus =  -0.57  ! Warner et al 2005, Table 2          
        c_mu_zero = 0.5540   
      end if 
      cc = 1.0  ! Warner et al 2005, top of page 85      
    case ("k-omega")
      gls_p = -1.0
      gls_m = 0.5
      gls_n = -1.0
      gls_sigma_k = 2.0  ! turbulent Schmidt number
      gls_sigma_psi = 2.0  ! turbulent Schmidt number
      c1 = 0.555
      c2 = 0.833
      c3plus = 1.0 
      k_min = 7.6e-6
      psi_min = 1.e-12
      call set( Fwall, 1.0 ) 
      if((trim(gls_stability_function)=="KanthaClayson-94").or.(trim(gls_stability_function)=="Galperin-88")) then      
        c3minus = -0.58  ! Warner et al 2005, Table 2
        c_mu_zero = 0.5544      
      else if(trim(gls_stability_function)=="Canuto-01-A") then
        c3minus = -0.64  ! Warner et al 2005, Table 2
        c_mu_zero = 0.5270   
      else if(trim(gls_stability_function)=="Canuto-01-B") then 
        !c3minus =  not_defined  ! Warner et al 2005, Table 2          
        c_mu_zero = 0.5540   
      end if
      cc = sqrt(2.)*c_mu_zero**3 ! Warner et al 2005, top of page 85
    case ("gen")
      gls_p = 2.0
      gls_m = 1.0
      gls_n = -0.67
      gls_sigma_k = 0.8  ! turbulent Schmidt number
      gls_sigma_psi = 1.07  ! turbulent Schmidt number
      c1 = 1.0
      c2 = 1.22
      c3plus = 1.0 
      k_min = 7.6e-6
      psi_min = 1.e-12
      call set( Fwall, 1.0 ) 
      if((trim(gls_stability_function)=="KanthaClayson-94").or.(trim(gls_stability_function)=="Galperin-88")) then      
        c3minus = 0.1  ! Warner et al 2005, Table 2
        c_mu_zero = 0.5544      
      else if(trim(gls_stability_function)=="Canuto-01-A") then
        c3minus = 0.05  ! Warner et al 2005, Table 2
        c_mu_zero = 0.5270   
      else if(trim(gls_stability_function)=="Canuto-01-B") then 
        !c3minus =  not_defined  ! Warner et al 2005, Table 2          
        c_mu_zero = 0.5540     
      end if
      cc = sqrt(2.)*c_mu_zero**3  ! Warner et al 2005, top of page 85           
    case default
      FLAbort("Unknown gls_option")           
    end select


    !Do some limiting before we calculate extra fields
    call gls_limiting(ll,kk,psi,epsilon,NN2,k_min,psi_min,gls_n,gls_p,gls_m,c_mu_zero)
    
    do ii = 1,nonods

      !!!!call set(psi,     ii, (c_mu_zero**gls_p)*(node_val(kk,ii)**gls_m)*(node_val(ll,ii)**gls_n)  !Warner et al 2005, Eqn. (14)

      call set(epsilon, ii, (c_mu_zero**(3.+(gls_p/gls_n)))*(node_val(kk,ii)**(1.5 + gls_m/gls_n))*(max(1.0e-10,node_val(psi,ii))**(-1./gls_n)))   !Warner et al 2005, Eqn. (12)

!!!      call set(ll,      ii, (c_mu_zero**3)*(node_val(kk,ii)**1.5)*(1./max(1.e-10,node_val(epsilon,ii))))    !Warner et al 2005, Eqn. (15)
      call set(ll,      ii, (c_mu_zero**(-gls_p/gls_n))*(node_val(kk,ii)**(-gls_m/gls_n))*(max(1.0e-10,node_val(psi,ii))**(1./gls_n)) )     !Warner et al 2005, Table 5
     
    end do 
!
  end subroutine gls_parameters_and_wall_proximity
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine gls_stability_functions(gls_stability_function,MM2,NN2,ll,kk,epsilon,S_H,S_M,c_mu_zero,GLSGH,output_GLSGH,XX)

    character(len=OPTION_PATH_LEN), intent(in) :: gls_stability_function
    type(scalar_field), intent(in) :: MM2,NN2,ll,kk,epsilon
    type(scalar_field), intent(inout) :: S_H,S_M   ! the stability functions
    type(scalar_field), intent(inout) :: GLSGH
    type(vector_field), intent(in) :: XX
    real, intent(in) :: c_mu_zero
    ! for use with KanthaClayson-94 and Galperin-88
    real :: G_h_unlimited,G_h,G_h_limited,AA1,AA2,BB1,BB2,CC2,CC3,G_h_min,G_h0,G_h_crit
    ! for use with Canuto
    real :: s0,s1,s2,s4,s5,s6,b0,b1,b2,b3,b4,b5,cff,f6,G_m_unlimited,G_m
    real :: lam1,lam2,lam3,lam4,lam5,lam6,lam7,lam8
    !
    integer :: nonods,ii
    logical :: output_GLSGH

    select case (trim(gls_stability_function))
    case ("KanthaClayson-94")
! parameters for Kantha and Clayson (2004)
      AA1 = 0.92
      AA2 = 0.74
      BB1 = 16.6
      BB2 = 10.1
      CC2 = 0.7
      CC3 = 0.2
!      G_h_min = -0.56
!      G_h0 = 0.0466 ! Warner et al 2005 cite Kantha and Clayson (1994) G_h0 = 0.028(?), used to "assure positive defiteness of all velocity variances"
!      G_h_crit = 0.02
      G_h_min = -0.28
      G_h0 = 0.0233 ! Warner et al 2005 cite Kantha and Clayson (1994) G_h0 = 0.028(?), used to "assure positive defiteness of all velocity variances"
      G_h_crit = 0.02
    case("Galperin-88")
      AA1 = 0.92
      AA2 = 0.74
      BB1 = 16.6
      BB2 = 10.1
      CC2 = 0.0       !Difference to KC
      CC3 = 0.0       !Difference to KC
      G_h_min = -0.28
      G_h0 = 0.0233 ! Warner et al 2005 cite Kantha and Clayson (1994) G_h0 = 0.028, used to "assure positive defiteness of all velocity variances"
      G_h_crit = 0.02  
    case("Canuto-01-A")
      lam1 = 0.107
      lam2 = 0.0032 
      lam3 = 0.0864 
      lam4 = 0.12 
      lam5 = 11.9 
      lam6 = 0.4 
      lam7 = 0.0 
      lam8 = 0.48 
      G_h_min = -0.28
      G_h0 = 0.0329
      G_h_crit = 0.03
    case("Canuto-01-B")
      lam1 = 0.127    !Difference to CA
      lam2 = 0.00336  !Difference to CA
      lam3 = 0.0906   !Difference to CA 
      lam4 = 0.101    !Difference to CA
      lam5 = 11.2     !Difference to CA 
      lam6 = 0.4 
      lam7 = 0.0 
      lam8 = 0.318    !Difference to CA
      G_h_min = -0.28
      G_h0 = 0.0444    !Difference to CA
      G_h_crit = 0.0414  !Difference to CA    
    case default
      FLAbort("Unknown gls_stability_function")       
    end select
!

    nonods = node_count(S_H)
    do ii = 1,nonods
      if((trim(gls_stability_function)=="KanthaClayson-94").or.(trim(gls_stability_function)=="Galperin-88")) then

        ! This is a buoyancy parameter
        G_h_unlimited = -(node_val(NN2,ii))*(node_val(ll,ii)**2)/(2*max(1.e-10,node_val(kk,ii)))   !Warner et al 2005, Eqn. (32)

        ! Doing some smoothing of G_h_unlimited
        G_h = min(G_h0,G_h_unlimited)
        if(G_h .gt. G_h_crit) G_h = min(G_h, G_h - ((G_h - G_h_crit)**2)/(G_h + G_h0 - 2*G_h_crit))   ! This is quite different to Warner et al !!!

        G_h_limited = min(max(G_h,G_h_min),G_h0)


        if(output_GLSGH) call set(GLSGH,ii,G_h_limited)
        ! KC quasi-equilibrium stability functions, Warner et al 2005 cite Kantha and Clayson (1994)   

        call set(S_H, ii, AA2*(1.-6.*AA1/BB1) / (1.-(3.*AA2*G_h_limited*(6.*AA1 + BB2*(1.-CC3)))) )  !Warner et al 2005, Eqn. (30)
        call set(S_M, ii, (BB1**(-1./3.) + (18.*AA1**2.+9.*AA1*AA2*(1.-CC2))*node_val(S_H,ii)*G_h_limited) / (1.-9.*AA1*AA2*G_h_limited) )  !Warner et al 2005, Eqn. (31)

      else
        s0 = 1.5*lam1*(lam5**2)
        s1 = -lam4*(lam6 + lam7) + 2.*lam4*lam5*(lam1 - lam2/3. - lam3) + 1.5*lam1*lam5*lam8
        s2 = -0.375*lam1*(lam6**2 - lam7**2)
        !Apparently there isn't an s3
        s4 = 2.*lam5
        s5 = 2.*lam4
        s6 = (2./3.)*lam5*(3.*lam3**2 - lam2**2) - 0.5*lam5*lam1*(3.*lam3 - lam2) + 0.75*lam1*(lam6 - lam7)
        b0 = 3.*lam5**2
        b1 = lam5*(7.*lam4 + 3.*lam8)
        b2 = (lam5**2)*(3.*lam3**2 - lam2**2) - 0.75*(lam6**2 - lam7**2)
        b3 = lam4*(4.*lam4 + 3.*lam8)
        b4 = lam4*(lam2*lam6 - 3.*lam3*lam7 - lam5*(lam2**2 - lam3**2)) + lam5*lam8*(3.*lam3**2 - lam2**2)
        b5 = 0.25*(lam2**2 - 3.*lam3**2)*(lam6**2 - lam7**2)

        f6 = 8./(c_mu_zero**6)  !Warner et al 2005, Eqn. (39)

     
        !These are the same as KanthaClayson-94 and Galperin-88
        G_h_unlimited = -(node_val(NN2,ii))*(node_val(ll,ii)**2)/(2.*max(1.e-10,node_val(kk,ii)))   !Warner et al 2005, Eqn. (32)


        G_h = min(G_h0,G_h_unlimited)
        if(G_h .gt. G_h_crit) G_h = min(G_h, G_h - ((G_h - G_h_crit)**2)/(G_h + G_h0 - 2*G_h_crit))

        G_h_limited = min(max(G_h,G_h_min),G_h0)

        if(output_GLSGH) call set(GLSGH,ii,G_h_limited)

        G_m_unlimited = (node_val(MM2,ii))*(node_val(ll,ii)**2)/(2.*max(1.e-10,node_val(kk,ii)))   !Warner et al 2005, Eqn. (40a)
        G_m = min( (b0/f6 - b1*G_h_limited + b3*f6*(G_h_limited**2))/(b2 - b4*f6*G_h_limited), G_m_unlimited )   !Warner et al 2005, Eqn. (40b)

        !
        cff = b0 - b1*f6*G_h_limited + b2*f6*G_m + b3*(f6**2)*(G_h_limited**2) - b4*(f6**2)*G_h_limited*G_m + b5*(f6**2)*(G_m**2) !Warner et al 2005, Eqn. (36)

        call set(S_M, ii, (sqrt(2.)/c_mu_zero**3)*max(0.0,(s0 - s1*f6*G_h_limited + s2*f6*G_m)/cff))   !Warner et al 2005, Eqn. (37)
        call set(S_H, ii, (sqrt(2.)/c_mu_zero**3)*max(0.0,(s4 - s5*f6*G_h_limited + s6*f6*G_m)/cff))   !Warner et al 2005, Eqn. (38)
        
       
      end if
    end do
!
  end subroutine gls_stability_functions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine gls_source_absor_visc_diff(source1,source2,absorption1,absorption2,epsilon,psi,cc,c1,c2,c3,MM2,NN2,kk,ll,S_M,S_H,K_M,K_H,Fwall, &
                                        gls_sigma_k,gls_sigma_psi,PP,BB, &
                                        viscosity,kk_diff,psi_diff,eddy_visc_KM,eddy_diff_KH)
    type(scalar_field), intent(inout) :: source1,source2,absorption1,absorption2,K_M,K_H,PP,BB
    type(scalar_field), intent(in) :: epsilon,psi,c3,MM2,NN2,kk,ll,S_M,S_H,Fwall
    type(tensor_field), intent(inout) :: viscosity,kk_diff,psi_diff,eddy_visc_KM,eddy_diff_KH 
    real, intent(in) :: gls_sigma_k,gls_sigma_psi,cc,c1,c2
    integer :: ii
!
    
    do ii = 1, node_count(K_H)    

! do some under-relaxation,    grab the previous value of the dim,dim entry of eddy viscosity and diffusivity
      call set(K_M, ii, 0.8*node_val(eddy_visc_KM,eddy_visc_KM%dim,eddy_visc_KM%dim,ii) + 0.2*cc*sqrt(2.*node_val(kk,ii))*node_val(ll,ii)*node_val(S_M,ii)  )    
      call set(K_H, ii, 0.8*node_val(eddy_diff_KH,eddy_diff_KH%dim,eddy_diff_KH%dim,ii) + 0.2*cc*sqrt(2.*node_val(kk,ii))*node_val(ll,ii)*node_val(S_H,ii)  )  

    end do
! 
    do ii = 1, node_count(source1)
      call set(PP, ii, node_val(K_M,ii)*node_val(MM2,ii))       !Warner et al 2005, Eqn. (10)
      call set(BB, ii, -1.0*node_val(K_H,ii)*node_val(NN2,ii))  !Warner et al 2005, Eqn. (11)


      if(.false.) then
      call set(source1, ii, ( node_val(PP,ii) + node_val(BB,ii) - node_val(epsilon,ii)) )    
      call set(source2, ii, ( (node_val(psi,ii)/node_val(kk,ii))*(c1*node_val(PP,ii) + node_val(c3,ii)*node_val(BB,ii) - c2*node_val(epsilon,ii)*node_val(Fwall,ii))) )

      else

      if( (node_val(PP,ii) + node_val(BB,ii)) .ge. 0.0 ) then  ! Safe to treat P+B in the explicit source term (cf ROMS code: gls_corstep.F)

        call set(source1,     ii, ( node_val(PP,ii) + node_val(BB,ii) ) )    ! P+B

        call set(absorption1, ii,   node_val(epsilon,ii)/node_val(kk,ii) )   ! epsilon/k

      else ! just add the positive P term and include B term with dissipation term (epsilon) in absorption (i.e. treat implicitly) see ROMS code, GOTM literature and Patanaker 1980

        call set(source1,     ii, ( node_val(PP,ii) ) )    !P

        call set(absorption1, ii, ( node_val(epsilon,ii) - node_val(BB,ii))/node_val(kk,ii) )  !(epsilon -B)/k

      endif

      endif


      if( (c1*node_val(PP,ii) + node_val(c3,ii)*node_val(BB,ii)) .ge. 0.0 ) then  ! Safe to treat P+B in the explicit source term (cf ROMS code: gls_corstep.F)

        call set(source2, ii, ( (node_val(psi,ii)/node_val(kk,ii))*(c1*node_val(PP,ii) + node_val(c3,ii)*node_val(BB,ii)) ) ) !c1*P + c3*B

        call set(absorption2, ii, c2*node_val(epsilon,ii)*node_val(Fwall,ii)/node_val(kk,ii) ) !c2*epsilon*Fwall/k

      else ! just add the positive P term and include B term with dissipation term (epsilon) in absorption (i.e. treat implicitly) see ROMS code, GOTM literature and Patanaker 1980

        call set(source2, ii, ( (node_val(psi,ii)/node_val(kk,ii))*(c1*node_val(PP,ii))) )  !c1*P

        call set(absorption2, ii,  (c2*node_val(epsilon,ii)*node_val(Fwall,ii) -  node_val(c3,ii)*node_val(BB,ii))/node_val(kk,ii) )   ! (c2*epsilon*Fwall - c3*B)/k

      endif


    end do




    ! add in the turbulent visc and diff in the dim,dim entry of tensors only at present

    call set(viscosity,viscosity%dim,viscosity%dim,K_M)  

!    call set(temperature_diffusivity,temperature_diffusivity%dim,temperature_diffusivity%dim,K_H)

    call set(kk_diff,kk_diff%dim,kk_diff%dim,K_M,scale=1./gls_sigma_k)
    
    call set(psi_diff,psi_diff%dim,psi_diff%dim,K_M,scale=1./gls_sigma_psi)

! add in a constant background diff/visc - this needs to be fixed to reflect a background value set in diamond
    do ii = 1, node_count(source1)
      call addto(viscosity,viscosity%dim,viscosity%dim,ii,1.0e-6)  
!      call addto(temperature_diffusivity,temperature_diffusivity%dim,temperature_diffusivity%dim,ii,1.0e-6)
      call addto(kk_diff,kk_diff%dim,kk_diff%dim,ii,1.0e-6)
      call addto(psi_diff,psi_diff%dim,psi_diff%dim,ii,1.0e-6)
    end do 

!
  end subroutine gls_source_absor_visc_diff
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine gls_limiting(ll,kk,psi,epsilon,NN2,k_min,psi_min,gls_n,gls_p,gls_m,c_mu_zero)
    type(scalar_field), intent(in) :: NN2
    type(scalar_field), intent(inout) :: ll,kk,psi,epsilon
    real, intent(in) :: k_min,psi_min,gls_n,gls_p,gls_m,c_mu_zero
    integer :: ii
    real :: limit
      
  ! Limit the length scale ll to represent the limiting effects of stable stratification
  ! Also limit k and psi 
  
    do ii = 1,node_count(ll)
      call set(kk,  ii, max( node_val(kk,ii), k_min ))
      call set(psi, ii, max( node_val(psi,ii), psi_min ))
      call set(ll,       ii, min( node_val(ll,ii),      sqrt(0.56*node_val(kk,ii)/max(1.0e-10,node_val(NN2,ii))) ))  !Warner et al 2005, Eqn. (42), GOTM manual (2.71)
      limit = (0.56**(gls_n/2.))*(c_mu_zero**(gls_p))*(max(1.0e-10,node_val(kk,ii))**(gls_m + gls_n/2.0))*(max(1.0e-10,node_val(NN2,ii))**(-gls_n/2.0))
      if(gls_n>0.0) then
        call set(psi, ii, min( node_val(psi,ii), limit )) !Warner et al 2005, Eqn. (43)/table 5
      else
        call set(psi, ii, max( node_val(psi,ii), limit )) !Warner et al 2005, Eqn. (43)/table 5
      end if
    end do
!  
  end subroutine gls_limiting
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  subroutine gls_buoyancy_parameter(c3,c3plus,c3minus,NN2)
    type(scalar_field), intent(inout) :: c3
    real, intent(in) :: c3plus,c3minus
    type(scalar_field), intent(in) :: NN2
    integer :: ii
        
    do ii = 1,node_count(c3)
      if (node_val(NN2,ii) < 0.0 ) then
        call set(c3, ii, c3plus)   ! unstable stratification
      else
        call set(c3, ii, c3minus)  ! stable stratification
      end if
    end do

  end subroutine gls_buoyancy_parameter

  subroutine friction(kappa,state,z0s,z0b,gravity,bottom_dz,u_taus_squared,u_taub_squared)

   IMPLICIT NONE
   real, intent(in) :: kappa, gravity
   type(scalar_field) :: bottom_dz
   type(state_type), intent(in):: state
   real, intent(out), dimension(:) :: z0s,z0b,u_taus_squared,u_taub_squared
   integer :: NNodes, nobcs
   integer          :: i,ii, MaxIter
   real :: rr
   real :: charnock_val=1400.
   character(len=OPTION_PATH_LEN) bctype
   type(vector_field), pointer:: wind_surface_field, positions, velocity
   type(vector_field) :: bottom_velocity
   type(mesh_type) :: ocean_mesh, input_mesh
   real :: u_taub, z0s_min
   real, dimension(2) :: temp_vector_2D
   real, dimension(3) :: temp_vector_3D
   type(scalar_field), pointer :: distanceToTop, distanceToBottom
   integer, dimension(:), pointer :: surface_element_list, surface_nodes

   MaxIter = 10
   z0s_min = 0.03
   
   ! get meshes
   velocity => extract_vector_field(state, "Velocity")
   positions => extract_vector_field(state, "Coordinate")
   ! explicitly using the VelocityMesh here as the nodes on the BC should line
   ! up with the windSurfaceField
   input_mesh = extract_mesh(state, "VelocityMesh")
   
   ! grab stresses from velocity field - Surface
   distanceToTop => extract_scalar_field(state, "DistanceToTop")
   call get_boundary_condition(distanceToTop, name='top', surface_element_list=surface_element_list)
   call create_surface_mesh(ocean_mesh, surface_nodes, input_mesh, surface_element_list, 'OceanSurface')
   NNodes = node_count(ocean_mesh) 
   nobcs = get_boundary_condition_count(velocity)
   do i=1, nobcs
     call get_boundary_condition(velocity, i, type=bctype, &
            surface_element_list=surface_element_list)
     if (bctype=='wind_forcing') then
         wind_surface_field => extract_surface_field(velocity, i, "WindSurfaceField")
     end if
   end do

   do i=1,NNodes
     temp_vector_2D = node_val(wind_surface_field,i)
     u_taus_squared(i) = (temp_vector_2D(1)**2+temp_vector_2D(2)**2)

     !  use the Charnock formula to compute the surface roughness
     z0s(i)=charnock_val*u_taus_squared(i)/gravity
     if (z0s(i).lt.z0s_min) z0s(i)=z0s_min

   end do
   call deallocate(ocean_mesh)

   ! grab values of velocity from bottom surface
   distanceToBottom => extract_scalar_field(state, "DistanceToBottom")
   call get_boundary_condition(distanceToBottom, name='bottom', surface_element_list=surface_element_list)
   call create_surface_mesh(ocean_mesh, surface_nodes, input_mesh, surface_element_list, 'OceanBottom')
   NNodes = node_count(ocean_mesh) 
   call allocate(bottom_velocity, 3, ocean_mesh, name="bottom_velocity")
   call remap_field_to_surface(velocity, bottom_velocity, &
                                surface_element_list)

   do i=1,NNodes
     temp_vector_3D = node_val(bottom_velocity,i)
     u_taub = sqrt(temp_vector_3D(1)**2+temp_vector_3D(2)**2+temp_vector_3D(3)**2)


     !  iterate bottom roughness length MaxItz0b times
     do ii=1,MaxIter
       z0b(i)=1e-7/max(1e-6,u_taub)+0.03*0.05

       !  compute the factor r (version 1, with log-law)
       ! Note that bottom_dz is already divided by 2
       rr=kappa/(log((z0b(i)+node_val(bottom_dz,i))/z0b(i)))

       !  compute the friction velocity at the bottom
       u_taub = rr*sqrt((temp_vector_2D(1)**2+temp_vector_2D(2)**2))

     end do

     u_taub_squared(i) = u_taub**2
  end do
  call deallocate(bottom_velocity)
  call deallocate(ocean_mesh)


  return
  end subroutine friction


end module gls

