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
      use state_module
      use fields
      use spud
      use global_parameters, only: OPTION_PATH_LEN
      use field_priority_lists
      use field_options
      use fetools
      use sparse_tools_petsc
      use solvers
      use sparsity_patterns
      use sparsity_patterns_meshes
      use profiler

      implicit none

      private
      public :: get_phase_submaterials, get_nonlinear_volume_fraction, &
                calculate_diagnostic_phase_volume_fraction, calculate_bulk_velocity, &
                add_fluid_particle_drag, add_heat_transfer, &
                add_particle_particle_drag

   contains

      subroutine get_phase_submaterials(state, istate, submaterials, phase_istate)
         !!< Sets up an array of the submaterials of a phase.
         !!< NB: This includes the current state itself (i.e. state(istate)).

         type(state_type), dimension(:), target, intent(inout) :: state
         integer, intent(in) :: istate
         type(state_type), dimension(:), pointer :: submaterials
         integer, intent(inout), optional :: phase_istate

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
         
         !! Assign the states to the submaterials array
         next = 1 ! Keep track of where we are in the submaterials array
         do i = 1, size(state)
            if(is_submaterial(i)) then
               submaterials(next) = state(i)

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
      
      
      subroutine calculate_bulk_velocity(states, v_field)
         !!< Calculates the bulk velocity \sum_{i=1}^N {vfrac_i*u_i}

         type(state_type), dimension(:), intent(inout) :: states
         type(vector_field), intent(inout) :: v_field ! Bulk velocity field

         ! Velocities of the continuous and particle phases
         type(vector_field), pointer :: u, x
         type(scalar_field), pointer :: vfrac

         integer :: i, dim, stat
         ! Counters over the elements and Gauss points
         integer :: ele, gi
         ! Transformed quadrature weights.
         real, dimension(ele_ngi(v_field, 1)) :: detwei

         ! Field values at each quadrature point.
         real, dimension(mesh_dim(v_field), ele_ngi(v_field, 1)) :: bulk_velocity_gi
         real, dimension(mesh_dim(v_field), ele_ngi(v_field, 1)) :: u_gi
         real, dimension(ele_ngi(v_field, 1)) :: vfrac_gi

         ! Current element global node numbers.
         integer, dimension(:), pointer :: bulk_velocity_nodes
         ! Current bulk_velocity element shape
         type(element_type), pointer :: bulk_velocity_shape

         ! Local bulk_velocity matrix and rhs for the current element.
         real, dimension(ele_loc(v_field, 1), ele_loc(v_field, 1)) :: mass_mat_addto
         real, dimension(mesh_dim(v_field), ele_loc(v_field, 1)) :: rhs_addto
         
         type(csr_sparsity), pointer :: sparsity
         type(csr_matrix) :: mass_mat
         type(vector_field) :: rhs

         character(len=OPTION_PATH_LEN) :: option_path


         ewrite(1,*) 'Entering calculate_bulk_velocity'

         ! Allocate / extract from state
         sparsity => get_csr_sparsity_firstorder(states, v_field%mesh, v_field%mesh)
         call allocate(mass_mat, sparsity, name = trim(v_field%name) // "MassMatrix")
         call allocate(rhs, v_field%dim, v_field%mesh, trim(v_field%name) // "Rhs")

         ! Zero bulk velocity field
         call zero(v_field)
         call zero(rhs)
         call zero(mass_mat)

         ! Loop through and integrate over each element
         do ele = 1, element_count(v_field)

            bulk_velocity_nodes => ele_nodes(v_field, ele)
            bulk_velocity_shape => ele_shape(v_field, ele)

            ! Get the base Coordinate field, and compute detwei
            x => extract_vector_field(states(1), "Coordinate")
            call transform_to_physical(x, ele, detwei = detwei)

            ! Zero arrays
            bulk_velocity_gi = 0.0
            rhs_addto = 0.0
            mass_mat_addto = 0.0
            
            do i = 1, size(states)
               u => extract_vector_field(states(i), "Velocity", stat)

               ! If there's no velocity then cycle
               if(stat/=0) cycle
               ! If this is an aliased velocity then cycle
               if(aliased(u)) cycle
               ! If the velocity isn't prognostic then cycle
               if(.not.have_option(trim(u%option_path)//"/prognostic")) cycle

               ! In case the v_field (the MaterialVolumeFraction field for now) is solved for explicitly,
               ! pass a Velocity field's solver options to PETSc instead.
               option_path=trim(u%option_path)

               u => extract_vector_field(states(i), "NonlinearVelocity")
               vfrac => extract_scalar_field(states(i), "PhaseVolumeFraction")

               ! Calculate the bulk_velocity at each quadrature point.
               u_gi = ele_val_at_quad(u, ele)
               vfrac_gi = ele_val_at_quad(vfrac, ele)

               ! Compute the bulk velocity
               ! (Assumes velocities are on the same mesh)
               do dim = 1, mesh_dim(v_field)
                  bulk_velocity_gi(dim,:) = bulk_velocity_gi(dim,:) + vfrac_gi*u_gi(dim,:)
               end do
     
            end do
            
            ! Compute the local mass matrix and rhs
            mass_mat_addto = shape_shape(bulk_velocity_shape, bulk_velocity_shape, detwei)
            rhs_addto = shape_vector_rhs(bulk_velocity_shape, bulk_velocity_gi, detwei)
            
            ! Add it to the global mass matrix and rhs
            call addto(mass_mat, bulk_velocity_nodes, bulk_velocity_nodes, mass_mat_addto)
            call addto(rhs, bulk_velocity_nodes, rhs_addto)

         end do

         ! Perform PETSc solve on the system of equations to get the bulk velocity field
         call petsc_solve(v_field, mass_mat, rhs, option_path=option_path)

         ewrite_minmax(v_field)

         call deallocate(mass_mat)
         call deallocate(rhs)

         ewrite(1,*) 'Exiting calculate_bulk_velocity'

      end subroutine calculate_bulk_velocity
      

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
         integer, parameter :: DRAG_CORRELATION_TYPE_STOKES = 1, DRAG_CORRELATION_TYPE_WEN_YU = 2, DRAG_CORRELATION_TYPE_ERGUN = 3
         
         ewrite(1, *) "Entering add_fluid_particle_drag"
         
         ! Let's check whether we actually have at least one particle phase.
         if(option_count("/material_phase/multiphase_properties/particle_diameter") == 0) then
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
         is_particle_phase = have_option(trim(state(istate)%option_path)//"/multiphase_properties/particle_diameter")
              
         ! Retrieve the index of the fluid phase in the state array.
         not_found = .true.
         if(is_particle_phase) then    
            do i = 1, size(state)
               if(.not.have_option(trim(state(i)%option_path)//"/multiphase_properties/particle_diameter")) then

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
               type(vector_field), pointer :: velocity_fluid, velocity_particle
               type(vector_field), pointer :: oldu_fluid, oldu_particle
               type(vector_field), pointer :: nu_fluid, nu_particle ! Non-linear approximation to the Velocities
               type(tensor_field), pointer :: viscosity_fluid
               type(scalar_field) :: nvfrac_fluid, nvfrac_particle
               real :: d ! Particle diameter
               character(len=OPTION_PATH_LEN) :: drag_correlation_name
               integer :: drag_correlation

               ! Get the necessary fields to calculate the drag force
               velocity_fluid => extract_vector_field(state(istate_fluid), "Velocity")
               velocity_particle => extract_vector_field(state(istate_particle), "Velocity")
               if(.not.aliased(velocity_particle)) then ! Don't count the aliased material_phases
                  
                  vfrac_fluid => extract_scalar_field(state(istate_fluid), "PhaseVolumeFraction")
                  vfrac_particle => extract_scalar_field(state(istate_particle), "PhaseVolumeFraction")
                  density_fluid => extract_scalar_field(state(istate_fluid), "Density")
                  density_particle => extract_scalar_field(state(istate_particle), "Density")
                  viscosity_fluid => extract_tensor_field(state(istate_fluid), "Viscosity")
         
                  call get_option(trim(state(istate_particle)%option_path)//"/multiphase_properties/particle_diameter", d)

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
                     case("default")
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
                                                            viscosity_fluid, d, drag_correlation)
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
                                                      viscosity_fluid, d, drag_correlation)
                                                         
               integer, intent(in) :: ele
               type(element_type), intent(in) :: test_function
               type(element_type), intent(in) :: u_shape
               type(vector_field), intent(in) :: u, x
               type(petsc_csr_matrix), intent(inout) :: big_m
               type(vector_field), intent(inout) :: mom_rhs
                    
               type(scalar_field), intent(in) :: vfrac_fluid, vfrac_particle
               type(scalar_field), intent(in) :: density_fluid, density_particle
               type(vector_field), intent(in) :: nu_fluid, nu_particle
               type(vector_field), intent(in) :: oldu_fluid, oldu_particle
               type(tensor_field), intent(in) :: viscosity_fluid    
               real, intent(in) :: d ! Particle diameter 
               integer, intent(in) :: drag_correlation
               
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

               ! Compute the particle Reynolds number
               ! (Assumes isotropic viscosity for now)
               particle_re_gi = (vfrac_fluid_gi*density_fluid_gi*magnitude_gi*d) / viscosity_fluid_gi(1,1,:)
           
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
               end select
                      
               ! Don't let the drag_coefficient_gi be NaN
               do gi = 1, ele_ngi(u,ele)
                  if(particle_re_gi(gi) == 0) then
                     drag_coefficient_gi(gi) = 1e16
                  end if
               end do
           
               select case(drag_correlation)
                  case(DRAG_CORRELATION_TYPE_STOKES)
                     K = vfrac_particle_gi*(3.0/4.0)*drag_coefficient_gi*(vfrac_fluid_gi*density_fluid_gi*magnitude_gi)/(d)
                  case(DRAG_CORRELATION_TYPE_WEN_YU)
                     ! Wen & Yu (1966) drag term
                     K = vfrac_particle_gi*(3.0/4.0)*drag_coefficient_gi*(vfrac_fluid_gi*density_fluid_gi*magnitude_gi)/(d*vfrac_fluid_gi**2.7)
                  case(DRAG_CORRELATION_TYPE_ERGUN)
                     K = 150.0*((vfrac_particle_gi**2)*viscosity_fluid_gi(1,1,:))/(vfrac_fluid_gi*(d**2)) + 1.75*(vfrac_particle_gi*density_fluid_gi*magnitude_gi/d)
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
         is_particle_phase = have_option(trim(state(istate)%option_path)//"/multiphase_properties/particle_diameter")
              
         ! Retrieve the index of the fluid phase in the state array.
         not_found = .true.
         if(is_particle_phase) then    
            do i = 1, size(state)
               if(.not.have_option(trim(state(i)%option_path)//"/multiphase_properties/particle_diameter")) then

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

               ! Get the necessary fields to calculate the heat transfer term
               velocity_fluid => extract_vector_field(state(istate_fluid), "Velocity")
               velocity_particle => extract_vector_field(state(istate_particle), "Velocity")
               if(.not.aliased(velocity_particle)) then ! Don't count the aliased material_phases
                  
                  vfrac_fluid => extract_scalar_field(state(istate_fluid), "PhaseVolumeFraction")
                  vfrac_particle => extract_scalar_field(state(istate_particle), "PhaseVolumeFraction")
                  density_fluid => extract_scalar_field(state(istate_fluid), "Density")
                  density_particle => extract_scalar_field(state(istate_particle), "Density")
                  viscosity_fluid => extract_tensor_field(state(istate_fluid), "Viscosity")
         
                  call get_option(trim(state(istate_particle)%option_path)//"/multiphase_properties/particle_diameter", d)
                           
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
                                                      viscosity_fluid, d, k, C_fluid, &
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
                                                viscosity_fluid, d, k, C_fluid, &
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
               type(vector_field), intent(in) :: nu_fluid, nu_particle
               type(scalar_field), intent(in) :: internal_energy_fluid, internal_energy_particle
               type(scalar_field), intent(in) :: old_internal_energy_fluid, old_internal_energy_particle
               type(tensor_field), intent(in) :: viscosity_fluid    
               real, intent(in) :: d, k, C_fluid, C_particle, gamma
               
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

               ! Compute the particle Reynolds number
               ! (Assumes isotropic viscosity for now)
               particle_Re = (density_fluid_gi*velocity_magnitude*d) / viscosity_fluid_gi(1,1,:)
           
               ! Compute the Prandtl number
               ! (Assumes isotropic viscosity for now)
               ! Note: C_fluid (at constant volume) multiplied by gamma = C_fluid at constant pressure
               Pr = C_fluid*gamma*viscosity_fluid_gi(1,1,:)/k
               
               particle_Nu = (7.0 - 10.0*vfrac_fluid_gi + 5.0*vfrac_fluid_gi**2)*(1.0 + 0.7*(particle_Re**0.2)*(Pr**(1.0/3.0))) + &
                              (1.33 - 2.4*vfrac_fluid_gi + 1.2*vfrac_fluid_gi**2)*(particle_Re**0.7)*(Pr**(1.0/3.0))

               Q = (6.0*k*vfrac_particle_gi*particle_Nu)/(d**2)
               
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
      
      subroutine add_particle_particle_drag(state, istate, u, x, big_m, mom_rhs)
         !!< This computes the particle-particle drag force term, using
         !!< an extension of the expression by Syamlal (1985). 
         !!< See Neri et al. (2003) for the extended expression.
         
         type(state_type), dimension(:), intent(inout) :: state     
         integer, intent(in) :: istate
         type(vector_field), intent(in) :: u, x
         type(petsc_csr_matrix), intent(inout) :: big_m
         type(vector_field), intent(inout) :: mom_rhs
         
         ! Local variables              
         integer :: ele
         type(element_type) :: test_function
         type(element_type), pointer :: u_shape
         integer, dimension(:), pointer :: u_ele
         logical :: dg
         
         real :: dt, theta
         logical, dimension(u%dim, u%dim) :: block_mask ! Control whether the off diagonal entries are used
                 
         integer :: i, dim, ipair
         
         integer :: istate_particle_i, istate_particle_j
         character(len=OPTION_PATH_LEN) :: particle_i_name, particle_j_name
         
         real :: e ! The restitution coefficient (e)
         real :: alpha
         
         logical :: current_phase_is_particle_i
         
         ewrite(1, *) "Entering add_particle_particle_drag"
               
         ! Get the timestepping options
         call get_option("/timestepping/timestep", dt)
         call get_option(trim(u%option_path)//"/prognostic/temporal_discretisation/theta", &
                        theta)
                                          
         ! For the big_m matrix. Controls whether the off diagonal entries are used    
         block_mask = .false.
         do dim = 1, u%dim
            block_mask(dim, dim) = .true.
         end do

         ! Are we using a discontinuous Galerkin discretisation?
         dg = continuity(u) < 0
         
         ! Get the restitution coefficient (e) which is independent of the phases.
         call get_option("/multiphase_interaction/particle_particle_drag/restitution_coefficient", e)
         call get_option("/multiphase_interaction/particle_particle_drag/non_headon_collision_coefficient", alpha)

         ! Loop over all of the particle-particle phase pairs
         pair_loop: do ipair = 1, option_count("/multiphase_interaction/particle_particle_drag/particle_particle_pair")

            call get_option("/multiphase_interaction/particle_particle_drag/particle_particle_pair["//int2str(ipair-1)//"]/particle_phase_i_name", particle_i_name)
            call get_option("/multiphase_interaction/particle_particle_drag/particle_particle_pair["//int2str(ipair-1)//"]/particle_phase_j_name", particle_j_name)
                                 
            ! If the current phase is one of the phases in the pair, then we need to compute the drag term
            if((particle_i_name == trim(state(istate)%name)) .or. (particle_j_name == trim(state(istate)%name))) then

               ! Determine the indices of particle phases i and j
               do i=1,size(state)
                  if(trim(state(i)%name) == particle_i_name) then
                     istate_particle_i = i
                  end if

                  if(trim(state(i)%name) == particle_j_name) then
                     istate_particle_j = i
                  end if
               end do

               ! Is the current phase particle i? If not, we might have to swap some variables 
               ! around in the assembly to compute D_ij and not D_ji
               current_phase_is_particle_i = (trim(state(istate)%name) == particle_i_name)
               
               if(istate_particle_i == istate_particle_j) then
                  cycle ! No point assembling when the particle phases are the same
               else
                  call assemble_particle_particle_drag(istate_particle_i, istate_particle_j)
               end if
            
            end if

         end do pair_loop

         ewrite(1, *) "Exiting add_particle_particle_drag"

         contains

            subroutine assemble_particle_particle_drag(istate_particle_i, istate_particle_j)

               integer, intent(in) :: istate_particle_i, istate_particle_j

               type(scalar_field), pointer :: vfrac_i, vfrac_j
               type(scalar_field), pointer :: density_i, density_j
               type(vector_field), pointer :: velocity_i, velocity_j
               type(vector_field), pointer :: oldu_i, oldu_j
               type(vector_field), pointer :: nu_i, nu_j ! Non-linear approximation to the Velocities
               type(scalar_field) :: nvfrac_i, nvfrac_j
               real :: d_i, d_j ! Particle diameter
               real :: phi_i, phi_j ! Maximum packing volume fractions for the particle phases

               ewrite(1, *) "Entering assemble_particle_particle_drag"

               ! Get the necessary fields to calculate the drag force
               velocity_i => extract_vector_field(state(istate_particle_i), "Velocity")
               velocity_j => extract_vector_field(state(istate_particle_j), "Velocity")
               if(.not.(aliased(velocity_i) .or. aliased(velocity_j))) then ! Don't count the aliased material_phases
                  
                  vfrac_i => extract_scalar_field(state(istate_particle_i), "PhaseVolumeFraction")
                  vfrac_j => extract_scalar_field(state(istate_particle_j), "PhaseVolumeFraction")
                  density_i => extract_scalar_field(state(istate_particle_i), "Density")
                  density_j => extract_scalar_field(state(istate_particle_j), "Density")
         
                  call get_option("/material_phase["//int2str(istate_particle_i-1)//&
                           &"]/multiphase_properties/particle_diameter", d_i)
                  call get_option("/material_phase["//int2str(istate_particle_j-1)//&
                           &"]/multiphase_properties/particle_diameter", d_j)
                           
                  call get_option("/material_phase["//int2str(istate_particle_i-1)//&
                           &"]/multiphase_properties/max_packing_volume_fraction", phi_i)
                  call get_option("/material_phase["//int2str(istate_particle_j-1)//&
                           &"]/multiphase_properties/max_packing_volume_fraction", phi_j)

                  ! Calculate the non-linear approximation to the PhaseVolumeFractions
                  call allocate(nvfrac_i, vfrac_i%mesh, "NonlinearPhaseVolumeFraction")
                  call allocate(nvfrac_j, vfrac_j%mesh, "NonlinearPhaseVolumeFraction")
                  call zero(nvfrac_i)
                  call zero(nvfrac_j)
                  call get_nonlinear_volume_fraction(state(istate_particle_i), nvfrac_i)
                  call get_nonlinear_volume_fraction(state(istate_particle_j), nvfrac_j)
                  
                  ! Get the non-linear approximation to the Velocities
                  nu_i => extract_vector_field(state(istate_particle_i), "NonlinearVelocity")
                  nu_j => extract_vector_field(state(istate_particle_j), "NonlinearVelocity")
                  oldu_i => extract_vector_field(state(istate_particle_i), "OldVelocity")
                  oldu_j => extract_vector_field(state(istate_particle_j), "OldVelocity")

                  ! ----- Volume integrals over elements -------------           
                  call profiler_tic(u, "element_loop")
                  element_loop: do ele = 1, element_count(u)

                     if(.not.dg .or. (dg .and. element_owned(u,ele))) then
                        u_ele=>ele_nodes(u, ele)
                        u_shape => ele_shape(u, ele)
                        test_function = u_shape                         

                        call add_particle_particle_drag_element(ele, test_function, u_shape, &
                                                               x, u, big_m, mom_rhs, &
                                                               nvfrac_i, nvfrac_j, &
                                                               density_i, density_j, &
                                                               nu_i, nu_j, &
                                                               oldu_i, oldu_j, &
                                                               d_i, d_j, &
                                                               phi_i, phi_j, &
                                                               e)
                     end if

                  end do element_loop
                  call profiler_toc(u, "element_loop")

                  call deallocate(nvfrac_i)
                  call deallocate(nvfrac_j)

               end if

               ewrite(1, *) "Exiting assemble_particle_particle_drag"

            end subroutine assemble_particle_particle_drag
            
            subroutine add_particle_particle_drag_element(ele, test_function, u_shape, &
                                                      x, u, big_m, mom_rhs, &
                                                      vfrac_i, vfrac_j, &
                                                      density_i, density_j, &
                                                      nu_i, nu_j, &
                                                      oldu_i, oldu_j, &
                                                      d_i, d_j, &
                                                      phi_i, phi_j, &
                                                      e)
                                                         
               integer, intent(in) :: ele
               type(element_type), intent(in) :: test_function
               type(element_type), intent(in) :: u_shape
               type(vector_field), intent(in) :: u, x
               type(petsc_csr_matrix), intent(inout) :: big_m
               type(vector_field), intent(inout) :: mom_rhs
                    
               type(scalar_field), intent(in) :: vfrac_i, vfrac_j
               type(scalar_field), intent(in) :: density_i, density_j
               type(vector_field), intent(in) :: nu_i, nu_j
               type(vector_field), intent(in) :: oldu_i, oldu_j
               
               real, intent(in) :: d_i, d_j
               real, intent(in) :: phi_i, phi_j
               real, intent(in) :: e
               
               ! Local variables
               real, dimension(ele_ngi(u,ele)) :: vfrac_i_gi, vfrac_j_gi
               real, dimension(ele_ngi(u,ele)) :: density_i_gi, density_j_gi
               real, dimension(u%dim, ele_ngi(u,ele)) :: nu_i_gi, nu_j_gi
               
               real, dimension(u%dim, ele_loc(u, ele)) :: oldu_val
               
               real, dimension(ele_loc(u, ele), ele_ngi(u, ele), x%dim) :: du_t
               real, dimension(ele_ngi(u, ele)) :: detwei
               real, dimension(u%dim, ele_loc(u,ele)) :: interaction_rhs_mat
               real, dimension(ele_loc(u, ele), ele_loc(u, ele)) :: interaction_big_m_mat
               real, dimension(u%dim, ele_loc(u,ele)) :: rhs_addto
               real, dimension(u%dim, u%dim, ele_loc(u,ele), ele_loc(u,ele)) :: big_m_tensor_addto
               
               real, dimension(ele_ngi(u,ele)) :: magnitude ! |u_i - u_j|
               real, dimension(ele_ngi(u,ele)) :: C
               real :: a, temp_phi_i, temp_phi_j
               real, dimension(ele_ngi(u,ele)) :: vfrac_ij, F
               real, dimension(ele_ngi(u,ele)) :: K ! drag force = K*(u_j - u_i)
               
               real, dimension(ele_ngi(u,ele)) :: drag_force_big_m
               real, dimension(u%dim, ele_ngi(u,ele)) :: drag_force_rhs ! drag_force = K*(u_j - u_i)
                             
               integer :: dim, gi
               
               ! Compute detwei
               call transform_to_physical(x, ele, u_shape, dshape=du_t, detwei=detwei)
               
               ! Get the values of the necessary fields at the Gauss points
               vfrac_i_gi = ele_val_at_quad(vfrac_i, ele)
               vfrac_j_gi = ele_val_at_quad(vfrac_j, ele)
               density_i_gi = ele_val_at_quad(density_i, ele)
               density_j_gi = ele_val_at_quad(density_j, ele)
               nu_i_gi = ele_val_at_quad(nu_i, ele)
               nu_j_gi = ele_val_at_quad(nu_j, ele)             
         
               ! Compute the magnitude of the relative (non-linear) velocity |nu_i - nu_j|
               do gi = 1, ele_ngi(u,ele)
                  magnitude(gi) = norm2(nu_i_gi(:,gi) - nu_j_gi(:,gi))
               end do

               ! We need to use the indices i and j in the correct way so that d_i >= d_j
               if(d_i < d_j) then
                  ! Swap the variables phi_i and phi_j, and compute variable 'a' such that d_i >= d_j
                  ! Note that temporary variables are used here instead of doing: temp=phi_j; phi_j=phi_i; phi_i=temp
                  ! because we don't want to have to reset phi_i and phi_j each time in the calling function.
                  temp_phi_i = phi_j
                  temp_phi_j = phi_i

                  a = sqrt(d_i/d_j)
                  C = vfrac_j_gi / (vfrac_j_gi + vfrac_i_gi)
               else
                  ! d_i >= d_j already, so use the variables as they are
                  temp_phi_i = phi_i
                  temp_phi_j = phi_j

                  a = sqrt(d_j/d_i)
                  C = vfrac_i_gi / (vfrac_i_gi + vfrac_j_gi)
               end if

               ! vfrac_ij is the "maximum solids volume fraction of a random closely packed structure" (Syamlal, 1985)
               do gi = 1, ele_ngi(u, ele)
                  if(C(gi) <= temp_phi_i/(temp_phi_i + (1.0 - temp_phi_i)*temp_phi_j)) then
                     vfrac_ij(gi) = ((temp_phi_i - temp_phi_j) + (1.0-a)*(1.0-temp_phi_i)*temp_phi_j)*&
                                    &((temp_phi_i + (1.0-temp_phi_j)*temp_phi_i)/temp_phi_i)*C(gi) + temp_phi_j
                  else
                     vfrac_ij(gi) = (1.0-a)*(temp_phi_i + (1.0-temp_phi_i)*temp_phi_j)*(1.0 - C(gi)) + temp_phi_i
                  end if
               end do
               
               ! F is a function involving the volume fraction of both phases and vfrac_ij defined above
               F = (3.0*vfrac_ij**(1.0/3.0) + (vfrac_i_gi + vfrac_j_gi)**(1.0/3.0)) / & 
                   & (2.0*(vfrac_ij**(1.0/3.0) - (vfrac_i_gi + vfrac_j_gi)**(1.0/3.0)))
               
               ! Drag force = K*(u_i - u_j)
               K = F*alpha*(1.0+e)*density_i_gi*vfrac_i_gi*density_j_gi*vfrac_j_gi*&
                   &( ((d_i + d_j)**2) / (density_i_gi*d_i**3 + density_j_gi*d_j**3) )*magnitude

               drag_force_big_m = -K
               
               if(current_phase_is_particle_i) then
                  do dim = 1, u%dim
                     drag_force_rhs(dim,:) = K*(nu_j_gi(dim,:))
                  end do
               else
                  do dim = 1, u%dim
                     drag_force_rhs(dim,:) = K*(nu_i_gi(dim,:))
                  end do
               end if
               
               ! Form the element interaction/drag matrix
               interaction_big_m_mat = shape_shape(test_function, u_shape, detwei*drag_force_big_m)
               interaction_rhs_mat = shape_vector_rhs(test_function, drag_force_rhs, detwei)
              
               ! Add contribution  
               big_m_tensor_addto = 0.0            
               rhs_addto = 0.0
               if(current_phase_is_particle_i) then
                  oldu_val = ele_val(oldu_i, ele)
               else
                  oldu_val = ele_val(oldu_j, ele)
               end if

               do dim = 1, u%dim
                  big_m_tensor_addto(dim, dim, :, :) = big_m_tensor_addto(dim, dim, :, :) - dt*theta*interaction_big_m_mat
                  rhs_addto(dim,:) = rhs_addto(dim,:) + matmul(interaction_big_m_mat, oldu_val(dim,:)) + interaction_rhs_mat(dim,:)
               end do
               
               ! Add the contribution to mom_rhs
               call addto(mom_rhs, u_ele, rhs_addto) 
               ! Add to the big_m matrix
               call addto(big_m, u_ele, u_ele, big_m_tensor_addto, block_mask=block_mask)
               
            end subroutine add_particle_particle_drag_element

      end subroutine add_particle_particle_drag
      
   end module multiphase_module
