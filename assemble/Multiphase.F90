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
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

   module multiphase_module
      !! This module contains various subroutines and functions for 
      !! multiphase flow simulations
      use fldebug
      use spud
      use global_parameters, only: OPTION_PATH_LEN
      use sparse_tools
      use parallel_fields
      use fetools
      use fields
      use profiler
      use sparse_tools_petsc
      use state_module
      use field_options
      use field_priority_lists

      implicit none

      private
      public :: get_phase_submaterials, get_nonlinear_volume_fraction, &
                calculate_diagnostic_phase_volume_fraction, &
                add_fluid_particle_drag, add_heat_transfer

   contains

      subroutine get_phase_submaterials(state, istate, submaterials, phase_istate, submaterials_indices)
         !!< Sets up an array of the submaterials of a phase.
         !!< NB: This includes the current state itself (i.e. state(istate)).

         type(state_type), dimension(:), target, intent(inout) :: state
         integer, intent(in) :: istate
         type(state_type), dimension(:), pointer :: submaterials
         integer, intent(inout), optional :: phase_istate
         integer, dimension(:), pointer, optional :: submaterials_indices

         !! Local variables
         integer :: i, next, stat, material_count
         type(vector_field), pointer :: u
         character(len=OPTION_PATH_LEN) :: phase_name, target_name
         logical, dimension(:), pointer :: is_submaterial
         
         ewrite(1,*) 'Entering get_phase_submaterials'

         !! Store whether state(i) is a material or not in an array of logicals
         !! to save on computations in the second loop
         allocate(is_submaterial(size(state)))

         !! Get the number of submaterials so we can make submaterials the correct size

         material_count = 1 ! We will include state(istate) as one of the materials

         phase_name = trim(state(istate)%name)

         do i = 1, size(state)
        
            if(i == istate) then
               is_submaterial(i) = .true.
               cycle ! Already counted the current state
            end if
            
            u => extract_vector_field(state(i), "Velocity", stat)
            is_submaterial(i) = .false.

            ! If velocity field exists and is aliased to a phase's velocity field...
            if(stat == 0) then
               if(aliased(u)) then
                  ! ...then find out which phase it is aliased to. If it's the current phase,
                  ! then we have found one more submaterial.

                  ! Save the name of the phase that the current state is aliased to
                  call get_option(trim(state(i)%option_path)//"/vector_field::Velocity/aliased/material_phase_name", target_name)

                  if(target_name == phase_name) then
                     ! Found one more submaterial!
                     material_count = material_count + 1
                     is_submaterial(i) = .true.
                  end if
               end if
            end if

         end do

         ewrite(1,*) 'Number of sub-materials = ', material_count

         !! Allocate submaterials array
         allocate(submaterials(material_count))
         if (present(submaterials_indices)) then
           allocate(submaterials_indices(material_count))
         end if
         
         !! Assign the states to the submaterials array
         next = 1 ! Keep track of where we are in the submaterials array
         do i = 1, size(state)
            if(is_submaterial(i)) then
               submaterials(next) = state(i)
               if (present(submaterials_indices)) then
                 submaterials_indices(next) = i
               end if

               ! Keep track of the phase's index in the new submaterials array
               if(present(phase_istate) .and. (i == istate)) then
                  phase_istate = next
               end if

               next = next + 1
            end if
         end do

         deallocate(is_submaterial)

         ewrite(1,*) 'Exiting get_phase_submaterials'

      end subroutine get_phase_submaterials


      subroutine get_nonlinear_volume_fraction(state, nvfrac)
         !!< Computes the nonlinear approximation to the phase volume fraction
         !!< and stores it in a locally allocated field before assembling the momentum equation.

         type(state_type), intent(inout) :: state
         type(scalar_field), intent(inout) :: nvfrac

         type(scalar_field), pointer :: volumefraction, oldvolumefraction
         type(vector_field), pointer :: velocity
         type(scalar_field) :: remapvfrac

         integer :: stat
         real :: theta

         logical :: cap
         real:: u_cap_val, l_cap_val

         ewrite(1,*) 'Entering get_nonlinear_volume_fraction'


         volumefraction => extract_scalar_field(state, 'PhaseVolumeFraction')

         ! Calculate the non-linear PhaseVolumeFraction
         call remap_field(volumefraction, nvfrac)
              
         velocity => extract_vector_field(state, 'Velocity', stat=stat)
         if(stat==0) then
            call get_option(trim(velocity%option_path)//'/prognostic/temporal_discretisation/theta', &
                              theta, stat)
            if(stat==0) then
               call allocate(remapvfrac, nvfrac%mesh, "RemppedPhaseVolumeFraction")
               
               oldvolumefraction => extract_scalar_field(state, 'OldPhaseVolumeFraction')

               ewrite_minmax(oldvolumefraction)
               ewrite_minmax(volumefraction)

               call remap_field(oldvolumefraction, remapvfrac)
               
               call scale(nvfrac, theta)
               call addto(nvfrac, remapvfrac, (1.0-theta))
               
               call deallocate(remapvfrac)
            end if
         end if

         ! Cap the volume fraction to take care of under or overshoots.
         ! This will have typically occurred during advection.
         cap = (have_option(trim(complete_field_path(volumefraction%option_path))//"/cap_values"))  
         if(cap) then
            call get_option(trim(complete_field_path(volumefraction%option_path))//"/cap_values/upper_cap", &
                           u_cap_val, default=huge(0.0)*epsilon(0.0))
            call get_option(trim(complete_field_path(volumefraction%option_path))//"/cap_values/lower_cap", &
                           l_cap_val, default=-huge(0.0)*epsilon(0.0))
                  
            call bound(nvfrac, l_cap_val, u_cap_val)               
         end if

         ewrite_minmax(nvfrac)

         ewrite(1,*) 'Exiting get_nonlinear_volume_fraction'
      
      end subroutine get_nonlinear_volume_fraction


      subroutine calculate_diagnostic_phase_volume_fraction(state)
         !!< Searches for the state with the diagnostic PhaseVolumeFraction field,
         !!< and then computes it using the formula:
         !!< diagnostic volume fraction = 1.0 - (sum of all other volume fractions)

         type(state_type), dimension(:), intent(inout) :: state
         
         ! Local variables
         type(scalar_field), pointer :: phasevolumefraction
         integer :: i, stat, diagnostic_count
         type(scalar_field) :: sumvolumefractions
         type(scalar_field), pointer :: sfield
         logical :: diagnostic

         ewrite(1,*) 'Entering calculate_diagnostic_phase_volume_fraction'

         diagnostic_count = option_count("/material_phase/scalar_field::PhaseVolumeFraction/diagnostic")
         if(diagnostic_count>1) then
            ewrite(0,*) diagnostic_count, ' diagnostic PhaseVolumeFractions'
            FLExit("Only one diagnostic PhaseVolumeFraction permitted.")
         end if

         if(diagnostic_count==1) then
            ! Find the diagnostic volume fraction
            state_loop: do i = 1, size(state)
               phasevolumefraction=>extract_scalar_field(state(i), 'PhaseVolumeFraction', stat)
               if (stat==0) then
                  diagnostic = (have_option(trim(phasevolumefraction%option_path)//'/diagnostic'))
                  if((.not. aliased(phasevolumefraction)).and. diagnostic) then
                     exit state_loop
                  end if
               end if
            end do state_loop
            
            call allocate(sumvolumefractions, phasevolumefraction%mesh, 'Sum of volume fractions')
            call zero(sumvolumefractions)
            
            do i = 1, size(state)
               sfield=>extract_scalar_field(state(i),'PhaseVolumeFraction',stat)
               diagnostic=(have_option(trim(sfield%option_path)//'/diagnostic'))
               if ( (stat==0).and.(.not. aliased(sfield)).and.(.not.diagnostic)) then
                  call addto(sumvolumefractions, sfield)
               end if
            end do
            
            call set(phasevolumefraction, 1.0)
            call addto(phasevolumefraction, sumvolumefractions, -1.0)
            call deallocate(sumvolumefractions)
         end if

         ewrite(1,*) 'Exiting calculate_diagnostic_phase_volume_fraction'

      end subroutine calculate_diagnostic_phase_volume_fraction
      

      !! Multiphase interaction terms, F_i
      subroutine add_fluid_particle_drag(state, istate, u, x, big_m, mom_rhs)
         !!< This computes the fluid-particle drag force term.
         !!< Note that this assumes only one fluid phase, and one or more particle phases.
         
         type(state_type), dimension(:), intent(inout) :: state     
         integer, intent(in) :: istate
         type(vector_field), intent(in) :: u, x
         type(petsc_csr_matrix), intent(inout) :: big_m
         type(vector_field), intent(inout) :: mom_rhs
         
         ! Local variables              
         integer :: ele
         type(element_type) :: test_function
         type(element_type), pointer :: u_shape
         integer, dimension(:), pointer :: u_nodes
         logical :: dg
             
         type(vector_field), pointer :: velocity_fluid
                 
         logical :: is_particle_phase
         
         real :: dt, theta
         logical, dimension(u%dim, u%dim) :: block_mask ! Control whether the off diagonal entries are used
                 
         integer :: i, dim
         logical :: not_found ! Error flag. Have we found the fluid phase?
         integer :: istate_fluid, istate_particle
   
         ! Types of drag correlation
         integer, parameter :: DRAG_CORRELATION_TYPE_STOKES = 1, DRAG_CORRELATION_TYPE_WEN_YU = 2, DRAG_CORRELATION_TYPE_ERGUN = 3, DRAG_CORRELATION_TYPE_SCHILLER_NAUMANN = 4, DRAG_CORRELATION_TYPE_LAIN_1_1999= 5, DRAG_CORRELATION_TYPE_LAIN_2_2002 = 6         
         
         ewrite(1, *) "Entering add_fluid_particle_drag"
         
         ! Let's check whether we actually have at least one particle phase.
         if( (option_count("/material_phase/multiphase_properties/particle_diameter") == 0) .and. (option_count("/material_phase/multiphase_properties/particle_dia_use_scalar_field") == 0) ) then
            FLExit("Fluid-particle drag enabled but no particle_diameter has been specified for the particle phase(s).")
         end if
               
         ! Get the timestepping options
         call get_option("/timestepping/timestep", dt)
         call get_option(trim(u%option_path)//"/prognostic/temporal_discretisation/theta", theta)
                                          
         ! For the big_m matrix. Controls whether the off diagonal entries are used    
         block_mask = .false.
         do dim = 1, u%dim
            block_mask(dim, dim) = .true.
         end do

         ! Are we using a discontinuous Galerkin discretisation?
         dg = continuity(u) < 0

         ! Is this phase a particle phase?
         is_particle_phase = have_option(trim(state(istate)%option_path)//"/multiphase_properties/particle_diameter") .or. have_option(trim(state(istate)%option_path)//"/multiphase_properties/particle_dia_use_scalar_field")
              
         ! Retrieve the index of the fluid phase in the state array.
         not_found = .true.
         if(is_particle_phase) then    
            do i = 1, size(state)
               if(.not.(have_option(trim(state(i)%option_path)//"/multiphase_properties/particle_diameter") .or. have_option(trim(state(i)%option_path)//"/multiphase_properties/particle_dia_use_scalar_field"))) then

                  velocity_fluid => extract_vector_field(state(i), "Velocity")
                  ! Aliased material_phases will also not have a particle_diameter,
                  ! so here we make sure that we don't count these as the fluid phase
                  if(.not.aliased(velocity_fluid)) then
                     istate_fluid = i
                     if(.not.not_found) then
                        FLExit("Fluid-particle drag does not currently support more than one fluid phase.")
                     end if
                     not_found = .false.
                  end if

               end if
            end do
         else
            istate_fluid = istate
            not_found = .false.
         end if

         if(not_found) then
            FLExit("No fluid phase found for the fluid-particle drag.")
         end if

         ! If we have a fluid-particle pair, then assemble the drag term
         if(is_particle_phase) then
            call assemble_fluid_particle_drag(istate_fluid, istate)
         else
            state_loop: do i = 1, size(state)
               if(i /= istate_fluid) then
                  call assemble_fluid_particle_drag(istate_fluid, i)
               end if
            end do state_loop
         end if
         
         ewrite(1, *) "Exiting add_fluid_particle_drag"
         
         contains

            subroutine assemble_fluid_particle_drag(istate_fluid, istate_particle)

               integer, intent(in) :: istate_fluid, istate_particle

               type(scalar_field), pointer :: vfrac_fluid, vfrac_particle
               type(scalar_field), pointer :: density_fluid, density_particle
               type(scalar_field), pointer :: d_field ! scalar field defining particle diameter 
               type(vector_field), pointer :: velocity_fluid, velocity_particle
               type(vector_field), pointer :: oldu_fluid, oldu_particle
               type(vector_field), pointer :: nu_fluid, nu_particle ! Non-linear approximation to the Velocities
               type(tensor_field), pointer :: viscosity_fluid
               type(scalar_field) :: nvfrac_fluid, nvfrac_particle
               real :: d ! Particle diameter
               real :: d_cap_lower
               character(len=OPTION_PATH_LEN) :: drag_correlation_name
               character(len=OPTION_PATH_LEN) :: d_field_name ! name of scalar field that defines particle diameter (can be the sauter mean dia)
               integer :: drag_correlation
               logical :: have_constant_d ! checks if the particle diameter is a constant or not 

               ! Get the necessary fields to calculate the drag force
               velocity_fluid => extract_vector_field(state(istate_fluid), "Velocity")
               velocity_particle => extract_vector_field(state(istate_particle), "Velocity")
               if(.not.aliased(velocity_particle)) then ! Don't count the aliased material_phases
                  
                  vfrac_fluid => extract_scalar_field(state(istate_fluid), "PhaseVolumeFraction")
                  vfrac_particle => extract_scalar_field(state(istate_particle), "PhaseVolumeFraction")
                  density_fluid => extract_scalar_field(state(istate_fluid), "Density")
                  density_particle => extract_scalar_field(state(istate_particle), "Density")
                  ! Make sure that the molecular viscosity is used. Not the effective viscosity
                  if(have_option(trim(state(istate_fluid)%option_path)//'/subgridscale_parameterisations/k-epsilon')) then
                     viscosity_fluid => extract_tensor_field(state(istate_fluid),"BackgroundViscosity")
                  else
                     viscosity_fluid => extract_tensor_field(state(istate_fluid), "Viscosity")
                  end if
                  if(have_option(trim(state(istate_particle)%option_path)//"/multiphase_properties/particle_diameter")) then
                     call get_option(trim(state(istate_particle)%option_path)//"/multiphase_properties/particle_diameter", d)
                     have_constant_d = .true.
                  else if(have_option(trim(state(istate_particle)%option_path)//"/multiphase_properties/particle_diameter_use_scalar_field")) then
                     call get_option(trim(state(istate_particle)%option_path)//"/multiphase_properties/particle_diameter_use_scalar_field", d_field_name)
                     d_field => extract_scalar_field(state(istate_particle), d_field_name)
                     have_constant_d = .false.
                  end if

                  call get_option(trim(state(istate_particle)%option_path)//"/multiphase_properties/apply_diameter_cap/lower_cap", d_cap_lower, default = 1.0e-12)

                  ! Calculate the non-linear approximation to the PhaseVolumeFractions
                  call allocate(nvfrac_fluid, vfrac_fluid%mesh, "NonlinearPhaseVolumeFraction")
                  call allocate(nvfrac_particle, vfrac_particle%mesh, "NonlinearPhaseVolumeFraction")
                  call zero(nvfrac_fluid)
                  call zero(nvfrac_particle)
                  call get_nonlinear_volume_fraction(state(istate_fluid), nvfrac_fluid)
                  call get_nonlinear_volume_fraction(state(istate_particle), nvfrac_particle)
                  
                  ! Get the non-linear approximation to the Velocities
                  nu_fluid => extract_vector_field(state(istate_fluid), "NonlinearVelocity")
                  nu_particle => extract_vector_field(state(istate_particle), "NonlinearVelocity")
                  oldu_fluid => extract_vector_field(state(istate_fluid), "OldVelocity")
                  oldu_particle => extract_vector_field(state(istate_particle), "OldVelocity")
                  
                  call get_option("/multiphase_interaction/fluid_particle_drag/drag_correlation/name", drag_correlation_name)
                  select case(trim(drag_correlation_name))
                     case("stokes")
                        drag_correlation = DRAG_CORRELATION_TYPE_STOKES
                     case("wen_yu")
                        drag_correlation = DRAG_CORRELATION_TYPE_WEN_YU
                     case("ergun")
                        drag_correlation = DRAG_CORRELATION_TYPE_ERGUN
                     case("schiller_naumann")
                        drag_correlation = DRAG_CORRELATION_TYPE_SCHILLER_NAUMANN
                     case("lain_1_1999")
                        drag_correlation = DRAG_CORRELATION_TYPE_LAIN_1_1999
                     case("lain_2_2002")
                        drag_correlation = DRAG_CORRELATION_TYPE_LAIN_2_2002
                     case default
                        FLAbort("Unknown correlation for fluid-particle drag")
                  end select
                  
                  ! ----- Volume integrals over elements -------------           
                  call profiler_tic(u, "element_loop")
                  element_loop: do ele = 1, element_count(u)

                     if(.not.dg .or. (dg .and. element_owned(u,ele))) then
                        u_nodes => ele_nodes(u, ele)
                        u_shape => ele_shape(u, ele)
                        test_function = u_shape                         

                        call add_fluid_particle_drag_element(ele, test_function, u_shape, &
                                                            x, u, big_m, mom_rhs, &
                                                            nvfrac_fluid, nvfrac_particle, &
                                                            density_fluid, density_particle, &
                                                            nu_fluid, nu_particle, &
                                                            oldu_fluid, oldu_particle, &
                                                            viscosity_fluid, &
                                                            have_constant_d, d, d_field, d_cap_lower, &
                                                            drag_correlation)
                     end if

                  end do element_loop
                  call profiler_toc(u, "element_loop")

                  call deallocate(nvfrac_fluid)
                  call deallocate(nvfrac_particle)
               end if

            end subroutine assemble_fluid_particle_drag
         
            subroutine add_fluid_particle_drag_element(ele, test_function, u_shape, &
                                                      x, u, big_m, mom_rhs, &
                                                      vfrac_fluid, vfrac_particle, &
                                                      density_fluid, density_particle, &
                                                      nu_fluid, nu_particle, &
                                                      oldu_fluid, oldu_particle, &
                                                      viscosity_fluid, &
                                                      have_constant_d, d, d_field, d_cap_lower, &
                                                      drag_correlation)
                                                         
               integer, intent(in) :: ele
               type(element_type), intent(in) :: test_function
               type(element_type), intent(in) :: u_shape
               type(vector_field), intent(in) :: u, x
               type(petsc_csr_matrix), intent(inout) :: big_m
               type(vector_field), intent(inout) :: mom_rhs
                    
               type(scalar_field), intent(in) :: vfrac_fluid, vfrac_particle
               type(scalar_field), intent(in) :: density_fluid, density_particle
               type(scalar_field), intent(in) :: d_field ! Scalar field representing particle diameter
               type(vector_field), intent(in) :: nu_fluid, nu_particle
               type(vector_field), intent(in) :: oldu_fluid, oldu_particle
               type(tensor_field), intent(in) :: viscosity_fluid    
               real, intent(in) :: d ! Constant particle diameter 
               real, intent(in) :: d_cap_lower
               integer, intent(in) :: drag_correlation
               logical, intent(in) :: have_constant_d ! is true if particle diameter is a constant. is false if it is a scalar field (e.g. sauter mean dia)

               ! Local variables
               real, dimension(ele_ngi(u,ele)) :: vfrac_fluid_gi, vfrac_particle_gi
               real, dimension(ele_ngi(u,ele)) :: density_fluid_gi, density_particle_gi
               real, dimension(u%dim, ele_ngi(u,ele)) :: nu_fluid_gi, nu_particle_gi
               real, dimension(u%dim, u%dim, ele_ngi(u,ele)) :: viscosity_fluid_gi
               
               real, dimension(u%dim, ele_loc(u, ele)) :: oldu_val
               
               real, dimension(ele_loc(u, ele), ele_ngi(u, ele), x%dim) :: du_t
               real, dimension(ele_ngi(u, ele)) :: detwei
               real, dimension(u%dim, ele_loc(u,ele)) :: interaction_rhs
               real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: interaction_big_m_mat
               real, dimension(u%dim, ele_loc(u,ele)) :: rhs_addto
               real, dimension(u%dim, u%dim, ele_loc(u,ele), ele_loc(u,ele)) :: big_m_tensor_addto
               
               real, dimension(ele_ngi(u,ele)) :: particle_re_gi ! Particle Reynolds number
               real, dimension(ele_ngi(u,ele)) :: drag_coefficient_gi
               real, dimension(ele_ngi(u,ele)) :: magnitude_gi ! |v_f - v_p|
               real, dimension(ele_ngi(u,ele)) :: d_gi ! particle diameter at the Gauss points
               
               real, dimension(ele_ngi(u,ele)) :: K
               real, dimension(ele_ngi(u,ele)) :: drag_force_big_m
               real, dimension(u%dim, ele_ngi(u,ele)) :: drag_force_rhs ! drag_force = K*(v_f - v_p)
                             
               integer :: dim, gi
               
               ! Compute detwei
               call transform_to_physical(x, ele, u_shape, dshape=du_t, detwei=detwei)
               
               ! Get the values of the necessary fields at the Gauss points
               vfrac_fluid_gi = ele_val_at_quad(vfrac_fluid, ele)
               vfrac_particle_gi = ele_val_at_quad(vfrac_particle, ele)
               density_fluid_gi = ele_val_at_quad(density_fluid, ele)
               density_particle_gi = ele_val_at_quad(density_particle, ele)
               nu_fluid_gi = ele_val_at_quad(nu_fluid, ele)
               nu_particle_gi = ele_val_at_quad(nu_particle, ele)
               viscosity_fluid_gi = ele_val_at_quad(viscosity_fluid, ele)               
         
               ! Compute the magnitude of the relative velocity
               do gi = 1, ele_ngi(u,ele)
                  magnitude_gi(gi) = norm2(nu_fluid_gi(:,gi) - nu_particle_gi(:,gi))
               end do

               ! Compute the particle diameter at the Gauss points
               if(have_constant_d) then
                  d_gi = d
               else 
                  d_gi = ele_val_at_quad(d_field, ele)
               end if

               ! Cap diameter on the lower side to prevent large drag forces
               WHERE (d_gi<d_cap_lower)
                  d_gi = d_cap_lower
               END WHERE
         
               ! Compute the particle Reynolds number
               ! (Assumes isotropic viscosity for now)
               particle_re_gi = (vfrac_fluid_gi*density_fluid_gi*magnitude_gi*d_gi) / viscosity_fluid_gi(1,1,:)
           
               ! Compute the drag coefficient
               select case(drag_correlation)
                  case(DRAG_CORRELATION_TYPE_STOKES)
                     ! Stokes drag correlation
                     do gi = 1, ele_ngi(u,ele)
                        if(particle_re_gi(gi) < 1000) then
                           drag_coefficient_gi(gi) = (24.0/particle_re_gi(gi))
                        else
                           drag_coefficient_gi(gi) = 0.44
                        end if
                     end do
                     
                  case(DRAG_CORRELATION_TYPE_WEN_YU)
                     ! Wen & Yu (1966) drag correlation
                     do gi = 1, ele_ngi(u,ele)
                        if(particle_re_gi(gi) < 1000) then
                           drag_coefficient_gi(gi) = (24.0/particle_re_gi(gi))*(1.0+0.15*particle_re_gi(gi)**0.687)
                        else
                           drag_coefficient_gi(gi) = 0.44
                        end if
                     end do
                     
                  case(DRAG_CORRELATION_TYPE_ERGUN)
                     ! No drag coefficient is needed here.                  

                  case(DRAG_CORRELATION_TYPE_SCHILLER_NAUMANN)
                     ! Schiller & Naumann (1933) drag correlation, also same as the one implemented in Fluent
                     ! Since the particle Reynolds number definition currently contains vfrac_fluid in numerator,
                     ! we need to take that out by dividing as done below. 
                     do gi = 1, ele_ngi(u,ele)
                        if((particle_re_gi(gi)/vfrac_fluid_gi(gi)) < 1000) then
                           drag_coefficient_gi(gi) = (24.0/(particle_re_gi(gi)/vfrac_fluid_gi(gi)))*(1.0+0.15*(particle_re_gi(gi)/vfrac_fluid_gi(gi))**0.687)
                        else
                           drag_coefficient_gi(gi) = 0.44
                        end if
                     end do
                  case(DRAG_CORRELATION_TYPE_LAIN_1_1999)
                     ! Lain 1999 - Rigid Bubble
                     do gi = 1, ele_ngi(u,ele)
                        if((particle_re_gi(gi)/vfrac_fluid_gi(gi)) < 500) then
                           drag_coefficient_gi(gi) = (24.0/(particle_re_gi(gi)/vfrac_fluid_gi(gi)))*(1.0+0.15*(particle_re_gi(gi)/vfrac_fluid_gi(gi))**0.687)
                        else if ((particle_re_gi(gi)/vfrac_fluid_gi(gi)) < 1500) then
                           drag_coefficient_gi(gi) = 9.5E-5 * (particle_re_gi(gi)/vfrac_fluid_gi(gi))**1.397
                        else 
                           drag_coefficient_gi(gi) = 2.61
                        end if
                     end do
                  case(DRAG_CORRELATION_TYPE_LAIN_2_2002)
                     ! Lain 2002
                     do gi = 1, ele_ngi(u,ele)
                        if((particle_re_gi(gi)/vfrac_fluid_gi(gi)) <= 1.5) then
                           drag_coefficient_gi(gi) = 16.0/(particle_re_gi(gi)/vfrac_fluid_gi(gi))
                        else if ((particle_re_gi(gi)/vfrac_fluid_gi(gi)) <= 80.0) then
                           drag_coefficient_gi(gi) = 14.9/(particle_re_gi(gi)/vfrac_fluid_gi(gi))**0.78
                        else if ((particle_re_gi(gi)/vfrac_fluid_gi(gi)) <= 1500.0) then
                           drag_coefficient_gi(gi) = (48.0/(particle_re_gi(gi)/vfrac_fluid_gi(gi)))*(1.0-2.21*(particle_re_gi(gi)/vfrac_fluid_gi(gi))**(-0.5))&
                                                     & + 1.86E-15*(particle_re_gi(gi)/vfrac_fluid_gi(gi))**4.756
                        else
                           drag_coefficient_gi(gi) = 2.61
                        end if
                     end do
               end select
                      
               ! Don't let the drag_coefficient_gi be NaN
               do gi = 1, ele_ngi(u,ele)
                  if(particle_re_gi(gi) == 0) then
                     drag_coefficient_gi(gi) = 1e16
                  end if
               end do
           
               select case(drag_correlation)
                  case(DRAG_CORRELATION_TYPE_STOKES, DRAG_CORRELATION_TYPE_SCHILLER_NAUMANN)
                     K = vfrac_particle_gi*(3.0/4.0)*drag_coefficient_gi*(vfrac_fluid_gi*density_fluid_gi*magnitude_gi)/(d_gi)
                  case(DRAG_CORRELATION_TYPE_WEN_YU)
                     ! Wen & Yu (1966) drag term
                     K = vfrac_particle_gi*(3.0/4.0)*drag_coefficient_gi*(vfrac_fluid_gi*density_fluid_gi*magnitude_gi)/(d_gi*vfrac_fluid_gi**2.7)
                  case(DRAG_CORRELATION_TYPE_ERGUN)
                     K = 150.0*((vfrac_particle_gi**2)*viscosity_fluid_gi(1,1,:))/(vfrac_fluid_gi*(d_gi**2)) + &
                                1.75*(vfrac_particle_gi*density_fluid_gi*magnitude_gi/d_gi)
                  case(DRAG_CORRELATION_TYPE_LAIN_1_1999, DRAG_CORRELATION_TYPE_LAIN_2_2002)
                     K = vfrac_particle_gi*(3.0/4.0)*drag_coefficient_gi*(density_fluid_gi*magnitude_gi)/(d_gi)
               end select               
               
               if(is_particle_phase) then
                  drag_force_big_m = -K
                  do dim = 1, u%dim
                     drag_force_rhs(dim,:) = -K*(-nu_fluid_gi(dim,:))
                  end do
               else
                  drag_force_big_m = -K
                  do dim = 1, u%dim
                     drag_force_rhs(dim,:) = K*(nu_particle_gi(dim,:))
                  end do
               end if
               
               ! Form the element interaction/drag matrix
               interaction_big_m_mat = shape_shape(test_function, u_shape, detwei*drag_force_big_m)
               interaction_rhs = shape_vector_rhs(test_function, drag_force_rhs, detwei)
              
               ! Add contribution  
               big_m_tensor_addto = 0.0            
               rhs_addto = 0.0
               if(is_particle_phase) then
                  oldu_val = ele_val(oldu_particle, ele)
               else
                  oldu_val = ele_val(oldu_fluid, ele)
               end if
               do dim = 1, u%dim
                  big_m_tensor_addto(dim, dim, :, :) = big_m_tensor_addto(dim, dim, :, :) - dt*theta*interaction_big_m_mat
                  rhs_addto(dim,:) = rhs_addto(dim,:) + matmul(interaction_big_m_mat, oldu_val(dim,:)) + interaction_rhs(dim,:)
               end do
               
               ! Add the contribution to mom_rhs
               call addto(mom_rhs, u_nodes, rhs_addto) 
               ! Add to the big_m matrix
               call addto(big_m, u_nodes, u_nodes, big_m_tensor_addto, block_mask=block_mask)

            end subroutine add_fluid_particle_drag_element

      end subroutine add_fluid_particle_drag

      
      !! Multiphase energy interaction term Q_i
      !! to be added to the RHS of the internal energy equation.
      subroutine add_heat_transfer(state, istate, internal_energy, matrix, rhs)
         !!< This computes the inter-phase heat transfer term.
         !!< Only between fluid and particle phase pairs.
         !!< Uses the empirical correlation by Gunn (1978).
         
         type(state_type), dimension(:), intent(inout) :: state     
         integer, intent(in) :: istate
         type(scalar_field), intent(in) :: internal_energy         
         type(csr_matrix), intent(inout) :: matrix
         type(scalar_field), intent(inout) :: rhs
         
         ! Local variables              
         integer :: ele
         type(element_type) :: test_function
         type(element_type), pointer :: internal_energy_shape
         integer, dimension(:), pointer :: internal_energy_nodes
         logical :: dg
             
         type(vector_field), pointer :: x, velocity_fluid
         
         real :: dt, theta
         
         logical :: is_particle_phase
         logical :: not_found ! Error flag. Have we found the fluid phase?
         integer :: i, istate_fluid, istate_particle
         
         
         ewrite(1, *) "Entering add_heat_transfer"
               
         ! Get the timestepping options
         call get_option("/timestepping/timestep", dt)
         call get_option(trim(internal_energy%option_path)//"/prognostic/temporal_discretisation/theta", &
                        theta)
         
         ! Get the coordinate field from state(istate)
         x => extract_vector_field(state(istate), "Coordinate")
         ewrite_minmax(x)
         assert(x%dim == mesh_dim(internal_energy))
         assert(ele_count(x) == ele_count(internal_energy))

         ! Are we using a discontinuous Galerkin discretisation?
         dg = continuity(internal_energy) < 0

         ! Is this phase a particle phase?
         is_particle_phase = have_option(trim(state(istate)%option_path)//"/multiphase_properties/particle_diameter") .or. have_option(trim(state(istate)%option_path)//"/multiphase_properties/particle_dia_use_scalar_field")
              
         ! Retrieve the index of the fluid phase in the state array.
         not_found = .true.
         if(is_particle_phase) then    
            do i = 1, size(state)
               if(.not.(have_option(trim(state(i)%option_path)//"/multiphase_properties/particle_diameter") .or. have_option(trim(state(i)%option_path)//"/multiphase_properties/particle_dia_use_scalar_field"))) then

                  velocity_fluid => extract_vector_field(state(i), "Velocity")
                  ! Aliased material_phases will also not have a particle_diameter,
                  ! so here we make sure that we don't count these as the fluid phase
                  if(.not.aliased(velocity_fluid)) then
                     istate_fluid = i
                     if(.not.not_found) then
                        FLExit("Heat transfer term does not currently support more than one fluid phase.")
                     end if
                     not_found = .false.
                  end if

               end if
            end do
         else
            istate_fluid = istate
            not_found = .false.
         end if

         if(not_found) then
            FLExit("No fluid phase found for the heat transfer term.")
         end if

         ! If we have a fluid-particle pair, then assemble the heat transfer term
         if(is_particle_phase) then
            call assemble_heat_transfer(istate_fluid, istate)
         else
            state_loop: do i = 1, size(state)
               if(i /= istate_fluid) then
                  call assemble_heat_transfer(istate_fluid, i)
               end if
            end do state_loop
         end if
         
         ewrite(1, *) "Exiting add_heat_transfer"
         
         contains

            subroutine assemble_heat_transfer(istate_fluid, istate_particle)

               integer, intent(in) :: istate_fluid, istate_particle

               type(scalar_field), pointer :: vfrac_fluid, vfrac_particle
               type(scalar_field), pointer :: density_fluid, density_particle
               type(scalar_field), pointer :: d_field ! scalar field defining particle diameter
               type(vector_field), pointer :: velocity_fluid, velocity_particle
               type(scalar_field), pointer :: internal_energy_fluid, internal_energy_particle
               type(scalar_field), pointer :: old_internal_energy_fluid, old_internal_energy_particle
               type(vector_field), pointer :: nu_fluid, nu_particle ! Non-linear approximation to the Velocities
               type(tensor_field), pointer :: viscosity_fluid
               type(scalar_field) :: nvfrac_fluid, nvfrac_particle
               real :: d ! Particle diameter
               real :: k ! Effective gas conductivity
               real :: C_fluid, C_particle ! Specific heat of the fluid and particle phases at constant volume
               real :: gamma ! Ratio of specific heats for compressible phase
               integer :: kstat, cstat_fluid, cstat_particle, gstat
               character(len=OPTION_PATH_LEN) :: d_field_name ! name of scalar field that defines particle diameter (can be the sauter mean dia)
               logical :: have_constant_d ! checks if the particle diameter is a constant or not 

               ! Get the necessary fields to calculate the heat transfer term
               velocity_fluid => extract_vector_field(state(istate_fluid), "Velocity")
               velocity_particle => extract_vector_field(state(istate_particle), "Velocity")
               if(.not.aliased(velocity_particle)) then ! Don't count the aliased material_phases
                  
                  vfrac_fluid => extract_scalar_field(state(istate_fluid), "PhaseVolumeFraction")
                  vfrac_particle => extract_scalar_field(state(istate_particle), "PhaseVolumeFraction")
                  density_fluid => extract_scalar_field(state(istate_fluid), "Density")
                  density_particle => extract_scalar_field(state(istate_particle), "Density")
                  viscosity_fluid => extract_tensor_field(state(istate_fluid), "Viscosity")
         
                  if(have_option(trim(state(istate_particle)%option_path)//"/multiphase_properties/particle_diameter")) then
                     call get_option(trim(state(istate_particle)%option_path)//"/multiphase_properties/particle_diameter", d)
                     have_constant_d = .true.
                  else if(have_option(trim(state(istate_particle)%option_path)//"/multiphase_properties/particle_dia_use_scalar_field")) then
                     call get_option(trim(state(istate_particle)%option_path)//"/multiphase_properties/particle_dia_use_scalar_field", d_field_name)
                     d_field => extract_scalar_field(state(istate_particle), d_field_name)
                     have_constant_d = .false.
                  end if
                           
                  call get_option(trim(state(istate_fluid)%option_path)//"/multiphase_properties/effective_conductivity", k, kstat)
                           
                  call get_option(trim(state(istate_fluid)%option_path)//"/multiphase_properties/specific_heat", C_fluid, cstat_fluid)
                  
                  call get_option(trim(state(istate_particle)%option_path)//"/multiphase_properties/specific_heat", C_particle, cstat_particle)
                  
                  call get_option(trim(state(istate_fluid)%option_path)//"/equation_of_state/compressible/stiffened_gas/ratio_specific_heats", gamma, gstat)
                           
                  if(kstat /= 0) then
                     FLExit("For inter-phase heat transfer, an effective_conductivity needs to be specified for the fluid phase.")
                  end if                  
                  if(cstat_fluid /= 0 .or. cstat_particle /= 0) then
                     FLExit("For inter-phase heat transfer, a specific_heat needs to be specified for each phase.")
                  end if
                  if(gstat /= 0) then
                     FLExit("For inter-phase heat transfer, ratio_specific_heats needs to be specified for the compressible phase.")
                  end if

                  ! Calculate the non-linear approximation to the PhaseVolumeFractions
                  call allocate(nvfrac_fluid, vfrac_fluid%mesh, "NonlinearPhaseVolumeFraction")
                  call allocate(nvfrac_particle, vfrac_particle%mesh, "NonlinearPhaseVolumeFraction")
                  call zero(nvfrac_fluid)
                  call zero(nvfrac_particle)
                  call get_nonlinear_volume_fraction(state(istate_fluid), nvfrac_fluid)
                  call get_nonlinear_volume_fraction(state(istate_particle), nvfrac_particle)
                  
                  ! Get the non-linear approximation to the Velocities
                  nu_fluid => extract_vector_field(state(istate_fluid), "NonlinearVelocity")
                  nu_particle => extract_vector_field(state(istate_particle), "NonlinearVelocity")
                  
                  ! Get the current and old internal energy fields
                  internal_energy_fluid => extract_scalar_field(state(istate_fluid), "InternalEnergy")
                  internal_energy_particle => extract_scalar_field(state(istate_particle), "InternalEnergy")
                  old_internal_energy_fluid => extract_scalar_field(state(istate_fluid), "OldInternalEnergy")
                  old_internal_energy_particle => extract_scalar_field(state(istate_particle), "OldInternalEnergy")
      
                  ! ----- Volume integrals over elements -------------           
                  call profiler_tic(internal_energy, "element_loop")
                  element_loop: do ele = 1, element_count(internal_energy)

                     if(.not.dg .or. (dg .and. element_owned(internal_energy,ele))) then
                        internal_energy_nodes => ele_nodes(internal_energy, ele)
                        internal_energy_shape => ele_shape(internal_energy, ele)
                        test_function = internal_energy_shape                         

                        call add_heat_transfer_element(ele, test_function, internal_energy_shape, &
                                                      x, internal_energy, matrix, rhs, &
                                                      nvfrac_fluid, nvfrac_particle, &
                                                      density_fluid, density_particle, &
                                                      nu_fluid, nu_particle, &
                                                      internal_energy_fluid, &
                                                      internal_energy_particle, &
                                                      old_internal_energy_fluid, &
                                                      old_internal_energy_particle, &
                                                      viscosity_fluid, &
                                                      have_constant_d, d, d_field, &
                                                      k, C_fluid, &
                                                      C_particle, gamma)
                     end if

                  end do element_loop
                  call profiler_toc(internal_energy, "element_loop")

                  call deallocate(nvfrac_fluid)
                  call deallocate(nvfrac_particle)
               end if

            end subroutine assemble_heat_transfer
         
            subroutine add_heat_transfer_element(ele, test_function, internal_energy_shape, &
                                                x, internal_energy, matrix, rhs, &
                                                vfrac_fluid, vfrac_particle, &
                                                density_fluid, density_particle, &
                                                nu_fluid, nu_particle, &
                                                internal_energy_fluid, &
                                                internal_energy_particle, &
                                                old_internal_energy_fluid, &
                                                old_internal_energy_particle, &
                                                viscosity_fluid, &
                                                have_constant_d, d, d_field, &
                                                k, C_fluid, &
                                                C_particle, gamma)
                                                         
               integer, intent(in) :: ele
               type(element_type), intent(in) :: test_function
               type(element_type), intent(in) :: internal_energy_shape
               type(vector_field), intent(in) :: x
               type(scalar_field), intent(in) :: internal_energy
               type(csr_matrix), intent(inout) :: matrix
               type(scalar_field), intent(inout) :: rhs
                    
               type(scalar_field), intent(in) :: vfrac_fluid, vfrac_particle
               type(scalar_field), intent(in) :: density_fluid, density_particle
               type(scalar_field), intent(in) :: d_field ! Scalar field representing particle diameter
               type(vector_field), intent(in) :: nu_fluid, nu_particle
               type(scalar_field), intent(in) :: internal_energy_fluid, internal_energy_particle
               type(scalar_field), intent(in) :: old_internal_energy_fluid, old_internal_energy_particle
               type(tensor_field), intent(in) :: viscosity_fluid    
               real, intent(in) :: d ! Constant particle diameter 
               real, intent(in) :: k, C_fluid, C_particle, gamma
               logical, intent(in) :: have_constant_d ! is true if particle diameter is a constant. is false if it is a scalar field (e.g. sauter mean dia)
               
               ! Local variables
               real, dimension(ele_ngi(x,ele)) :: internal_energy_fluid_gi, internal_energy_particle_gi
               real, dimension(ele_ngi(x,ele)) :: vfrac_fluid_gi, vfrac_particle_gi
               real, dimension(ele_ngi(x,ele)) :: density_fluid_gi, density_particle_gi
               real, dimension(x%dim, ele_ngi(x,ele)) :: nu_fluid_gi, nu_particle_gi
               real, dimension(x%dim, x%dim, ele_ngi(x,ele)) :: viscosity_fluid_gi
               
               real, dimension(ele_loc(internal_energy,ele)) :: old_internal_energy_val
               
               real, dimension(ele_ngi(x,ele)) :: detwei
               real, dimension(ele_ngi(x,ele)) :: coefficient_for_matrix, coefficient_for_rhs
               real, dimension(ele_loc(internal_energy,ele)) :: rhs_addto
               real, dimension(ele_loc(internal_energy,ele), ele_loc(internal_energy,ele)) :: matrix_addto
               
               real, dimension(ele_ngi(x,ele)) :: d_gi ! particle diameter at quadrature points
               real, dimension(ele_ngi(x,ele)) :: particle_Re ! Particle Reynolds number
               real, dimension(ele_ngi(x,ele)) :: Pr ! Prandtl number
               real, dimension(ele_ngi(x,ele)) :: particle_Nu ! Particle Nusselt number
               real, dimension(ele_ngi(x,ele)) :: velocity_magnitude ! |v_f - v_p|
               
               real, dimension(ele_ngi(x,ele)) :: Q ! heat transfer term = Q*(T_p - T_f)
               real, dimension(ele_loc(internal_energy,ele), ele_loc(internal_energy,ele)) :: heat_transfer_matrix
               real, dimension(ele_loc(internal_energy,ele)) :: heat_transfer_rhs
               
               integer ::  gi
               
               ! Compute detwei
               call transform_to_physical(x, ele, detwei)
               
               ! Get the values of the necessary fields at the Gauss points
               vfrac_fluid_gi = ele_val_at_quad(vfrac_fluid, ele)
               vfrac_particle_gi = ele_val_at_quad(vfrac_particle, ele)
               density_fluid_gi = ele_val_at_quad(density_fluid, ele)
               density_particle_gi = ele_val_at_quad(density_particle, ele)
               nu_fluid_gi = ele_val_at_quad(nu_fluid, ele)
               nu_particle_gi = ele_val_at_quad(nu_particle, ele)
               viscosity_fluid_gi = ele_val_at_quad(viscosity_fluid, ele)  
               internal_energy_fluid_gi = ele_val_at_quad(internal_energy_fluid, ele)
               internal_energy_particle_gi = ele_val_at_quad(internal_energy_particle, ele)
         
               ! Compute the magnitude of the relative velocity
               do gi = 1, ele_ngi(x,ele)
                  velocity_magnitude(gi) = norm2(nu_fluid_gi(:,gi) - nu_particle_gi(:,gi))
               end do

               ! Compute the particle diameter at quadrature points
               if(have_constant_d) then
                  d_gi = d
               else 
                  d_gi = ele_val_at_quad(d_field, ele)
               end if

               ! Compute the particle Reynolds number
               ! (Assumes isotropic viscosity for now)
               particle_Re = (density_fluid_gi*velocity_magnitude*d_gi) / viscosity_fluid_gi(1,1,:)
           
               ! Compute the Prandtl number
               ! (Assumes isotropic viscosity for now)
               ! Note: C_fluid (at constant volume) multiplied by gamma = C_fluid at constant pressure
               Pr = C_fluid*gamma*viscosity_fluid_gi(1,1,:)/k
               
               particle_Nu = (7.0 - 10.0*vfrac_fluid_gi + 5.0*vfrac_fluid_gi**2)*(1.0 + 0.7*(particle_Re**0.2)*(Pr**(1.0/3.0))) + &
                              (1.33 - 2.4*vfrac_fluid_gi + 1.2*vfrac_fluid_gi**2)*(particle_Re**0.7)*(Pr**(1.0/3.0))

               Q = (6.0*k*vfrac_particle_gi*particle_Nu)/(d_gi**2)
               
               ! Note that the transfer term is defined in terms of temperatures (T_fluid and T_particle)
               ! Let's convert the temperatures to internal energy (E) using E_i = C_i*T_i,
               ! where C is the specific heat of phase i at constant volume.
               if(is_particle_phase) then
                  coefficient_for_matrix = -Q/C_particle
                  coefficient_for_rhs = -Q*(-internal_energy_fluid_gi/C_fluid)
               else
                  coefficient_for_matrix = -Q/C_fluid
                  coefficient_for_rhs = Q*(internal_energy_particle_gi/C_particle)
               end if
               
               ! Form the element heat transfer matrix and RHS
               heat_transfer_matrix = shape_shape(test_function, internal_energy_shape, detwei*coefficient_for_matrix)
               heat_transfer_rhs = shape_rhs(test_function, coefficient_for_rhs*detwei)
              
               ! Add contribution  
               matrix_addto = 0.0            
               rhs_addto = 0.0
               if(is_particle_phase) then
                  old_internal_energy_val = ele_val(old_internal_energy_particle, ele)
               else
                  old_internal_energy_val = ele_val(old_internal_energy_fluid, ele)
               end if

               matrix_addto = matrix_addto - dt*theta*heat_transfer_matrix
               rhs_addto = rhs_addto + matmul(heat_transfer_matrix, old_internal_energy_val) + heat_transfer_rhs
               
               ! Add to the internal energy equation's RHS
               call addto(rhs, internal_energy_nodes, rhs_addto) 
               ! Add to the matrix
               call addto(matrix, internal_energy_nodes, internal_energy_nodes, matrix_addto)

            end subroutine add_heat_transfer_element

      end subroutine add_heat_transfer
      
   end module multiphase_module
