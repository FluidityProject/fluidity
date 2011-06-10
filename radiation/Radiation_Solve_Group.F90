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

module radiation_assemble_solve_group

   !!< This module contains procedures associated with solving the 
   !!< group g particle matrix problem

   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use fields
   use boundary_conditions
   use advection_diffusion_cg 
      
   use radiation_particle_data_type
   use radiation_materials_interpolation
   use radiation_extract_flux_field
   use radiation_energy_group_set_tools
     
   implicit none
   
   private 

   public :: particle_assemble_solve_group

contains

   ! --------------------------------------------------------------------------

   subroutine particle_assemble_solve_group(particle, &
                                            g, &
                                            invoke_eigenvalue_group_solve, &
                                            petsc_iterations_taken_group_g) 
   
      !!< Assemble and solve the group g particle matrix problem

      type(particle_type), intent(inout) :: particle
      integer, intent(in) :: g
      logical, intent(in) :: invoke_eigenvalue_group_solve
      integer, intent(out) :: petsc_iterations_taken_group_g
      
      ! local variables
      character(len=OPTION_PATH_LEN) :: field_name
      real :: dt
      type(scalar_field) :: extra_discretised_source
      type(scalar_field) :: mass_matrix_coeff      
      type(scalar_field), pointer :: t 
            
      ! assemble the diffusivity and absorption 
      ! scalar fields and insert into particle%state
      ! Also form the extra production/scatter/delayed discretised source
      ! and the 'density' term for time dependent problems
      call assemble_coeff_source_group_g(particle, &
                                         extra_discretised_source, &
                                         mass_matrix_coeff, &
                                         g, &
                                         invoke_eigenvalue_group_solve)
      
      ! form the field name to solve for
      field_name = 'ParticleFluxGroup'//int2str(g)//'Moment1'//trim(particle%name)

      eig_or_time: if (invoke_eigenvalue_group_solve) then
      
         dt = 1.0
      
      else eig_or_time
      
         call get_option('/timestepping/timestep',dt)
      
      end if eig_or_time
      
      ! this procedure is in assemble/Advection_Diffusion_CG.F90
      call solve_field_equation_cg(trim(field_name), &
                                   particle%state, &
                                   dt, &
                                   extra_discretised_source = extra_discretised_source, &
                                   iterations_taken = petsc_iterations_taken_group_g)
      
      ! set the direchlet bc nodes to be consistent
        
      t => extract_scalar_field(particle%state, &
                                trim(field_name))

      call set_dirichlet_consistent(t)
      
      ! deallocate the extra_discretised_source
      call deallocate(extra_discretised_source)
      
      ! deallocate the mass_matrix_coeff
      call deallocate(mass_matrix_coeff)
            
   end subroutine particle_assemble_solve_group

   ! --------------------------------------------------------------------------
   
   subroutine assemble_coeff_source_group_g(particle, &
                                            extra_discretised_source, &
                                            mass_matrix_coeff, &
                                            g, &
                                            invoke_eigenvalue_group_solve)
      
      !!< Assemble the coeff and discretised source fields for group g that will then be used 
      !!< somewhere else to assemble the linear system. For a time run the velocity coeff
      !!< are used to form a mass matrix coeff term (this does not happen
      !!< if the option is chosen to exclude the mass terms)

      type(particle_type), intent(inout) :: particle
      type(scalar_field), intent(out) :: extra_discretised_source
      type(scalar_field), intent(out) :: mass_matrix_coeff      
      integer, intent(in) :: g
      logical, intent(in) :: invoke_eigenvalue_group_solve
      
      ! local variables
      logical :: form_mass_matrix_coeff
      integer :: stat
      integer :: vele
      integer :: g_dash
      integer :: g_set
      real :: theta    
      real, dimension(:), allocatable :: detwei_vele
      real, dimension(:), allocatable :: rhs_addto      
      type(scalar_field) :: production_coeff
      type(scalar_field) :: prompt_spectrum_coeff
      type(scalar_field) :: scatter_coeff
      type(scalar_field), pointer :: absorption_coeff
      type(tensor_field), pointer :: diffusivity_coeff
      type(mesh_type), pointer :: material_fn_space
      type(vector_field), pointer :: positions      
      type(scalar_field_pointer), dimension(:), pointer :: particle_flux
      type(scalar_field_pointer), dimension(:), pointer :: particle_flux_old
      type(scalar_field_pointer), dimension(:), pointer :: particle_flux_iter     
      character(len=OPTION_PATH_LEN) :: material_fn_space_name
      character(len=OPTION_PATH_LEN) :: energy_group_set_path
      character(len=OPTION_PATH_LEN) :: field_name
      character(len=OPTION_PATH_LEN) :: mass_terms_path
                 
      ! extract the particle flux fields for all energy groups
      call extract_flux_all_group(particle, & 
                                  particle_flux      = particle_flux, &
                                  particle_flux_old  = particle_flux_old, &
                                  particle_flux_iter = particle_flux_iter)

      ! currently assume all energy groups point to the same coordinate mesh
      
      ! if the above assumption was to be broken then the there would need to be a sweep of the other groups
      ! fluxes (new, old and iter) to see which needs supermeshing to a new field and then the pointer within the 
      ! scalar_field_pointer particle_flux (or _old or _iter) is then changed such that no code change below is needed
      
      ! determine which group set this g belongs to
      call which_group_set_contains_g(g, &
                                      trim(particle%option_path), &
                                      g_set)

      ! set the energy_group_set path
      energy_group_set_path = trim(particle%option_path)//'/energy_discretisation/energy_group_set['//int2str(g_set - 1)//']'
         
      ! get the material fn space name for this group set
      call get_option(trim(energy_group_set_path)//'/angular_discretisation/method/parity/angular_moment_set[0]/mesh/name',material_fn_space_name)
      
      ! extract the material fn_space of this energy group set of this particle type 
      material_fn_space => extract_mesh(particle%state, &
                                        trim(material_fn_space_name))
      
      ! get the positions field for this energy group set
      positions => extract_vector_field(particle%state, &
                                        'Coordinate')
      
      ! deduce if for a time equation the mass terms will be included
      time_run: if (.not. invoke_eigenvalue_group_solve) then
         
         ! set the mass terms option path
         discr: if (have_option(trim(particle_flux(g)%ptr%option_path)//'/prognostic/spatial_discretisation/continuous_galerkin')) then
            
            mass_terms_path = trim(particle_flux(g)%ptr%option_path)//'/prognostic/spatial_discretisation/continuous_galerkin/mass_terms'
         
         else discr
            
            FLAbort('Unknown spatial discretisation for radiation model')
         
         end if discr
         
         exclude_mass: if (have_option(trim(mass_terms_path)//'/exclude_mass_terms')) then

            form_mass_matrix_coeff = .false.

         else exclude_mass

            form_mass_matrix_coeff = .true.

         end if exclude_mass
         
      else time_run
      
         form_mass_matrix_coeff = .false.
         
      end if time_run
      
      ! allocate the input discretised_source using the mesh of the group to solve for
      field_name = trim(particle_flux(g)%ptr%name)//'DiscretisedSource'
      call allocate(extra_discretised_source, particle_flux(g)%ptr%mesh, trim(field_name))
      call zero(extra_discretised_source)
      
      ! allocate the local fields as needed

      field_name = trim(particle_flux(g)%ptr%name)//'Scatter'
      call allocate(scatter_coeff, material_fn_space, trim(field_name))

      field_name = trim(particle_flux(g)%ptr%name)//'Production'
      call allocate(production_coeff, material_fn_space, trim(field_name))

      field_name = trim(particle_flux(g)%ptr%name)//'PromptSpectrum'
      call allocate(prompt_spectrum_coeff, material_fn_space, trim(field_name))
      
      ! extract the assemble fields as needed
      absorption_coeff => extract_scalar_field(particle%state, &
                                               trim(particle_flux(g)%ptr%name) // 'Absorption', &
                                               stat = stat)
                                        
      diffusivity_coeff => extract_tensor_field(particle%state, &
                                                trim(particle_flux(g)%ptr%name) // 'Diffusivity', &
                                                stat = stat)
      
      ! form the velocity coeff field for this energy group if time run
      field_name = trim(particle_flux(g)%ptr%name)//'MassMatrixCoeff'
      form_velocity: if (form_mass_matrix_coeff) then
                 
         call allocate(mass_matrix_coeff, material_fn_space, trim(field_name))
                  
         call form(particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                   particle%particle_radmat, &
                   mass_matrix_coeff, &
                   g, &
                   component = 'velocity')
         
         ! set as 1.0/velocity
         call invert(mass_matrix_coeff)
         
      else form_velocity
        
        ! form a constant 1.0 field
        call allocate(mass_matrix_coeff, material_fn_space, trim(field_name), field_type=FIELD_TYPE_CONSTANT)
        
        call set(mass_matrix_coeff, 1.0)
      
      end if form_velocity
            
      ! form the absorption coeff field for this energy group (which is actually the removal cross section)            
      call form(particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                particle%particle_radmat, &
                absorption_coeff, &
                g, &
                component = 'removal')
            
      ! form the diffusivity tensor coeff field for this energy group
      call form(particle%state, &
                particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                particle%particle_radmat, &
                diffusivity_coeff, &
                g, &
                positions%dim)
            
      ! form the discretised source - scatter and production            
      
      ! form the prompt spectrum coeff field for this energy group
      call form(particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                particle%particle_radmat, &
                prompt_spectrum_coeff, &
                g, &
                component = 'prompt_spectrum')
                  
      group_loop: do g_dash = 1,size(particle_flux)

         ! form the production coeff field for this energy group
         call form(particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                   particle%particle_radmat, &
                   production_coeff, &
                   g_dash, &
                   component = 'production')
         
         not_within_group: if (g /= g_dash) then
 
            ! form the scatter field for this energy group g_dash to g            
            call form(particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                      particle%particle_radmat, &
                      scatter_coeff, &
                      g, &
                      component = 'scatter', &
                      g_dash = g_dash)
                                           
         end if not_within_group

         ! get the time theta value for this energy group g to solve for
         get_theta: if (.not. invoke_eigenvalue_group_solve) then
         
            call get_option(trim(particle_flux(g)%ptr%option_path)//'/prognostic/temporal_discretisation/theta',theta)
         
         end if get_theta
         
         vele_loop: do vele = 1,ele_count(extra_discretised_source)
           
            ! allocate the jacobian transform and gauss weight array for this vele
            allocate(detwei_vele(ele_ngi(extra_discretised_source,vele)))
                        
            ! form the velement jacobian transform and gauss weight
            call transform_to_physical(positions, vele, detwei = detwei_vele)
            
            ! allocate the rhs_addto
            allocate(rhs_addto(ele_loc(extra_discretised_source,vele)))
            
            rhs_addto = 0.0
            
            keff_or_time: if (invoke_eigenvalue_group_solve) then
            
               ! add the scatter - not the within group
               not_within_group_eig: if (g /= g_dash) then
                  
                  rhs_addto = rhs_addto + matmul(shape_shape(ele_shape(particle_flux(g)%ptr,vele), &
                                                             ele_shape(particle_flux(g_dash)%ptr,vele), &
                                                             detwei_vele*ele_val_at_quad(scatter_coeff,vele)), &
                                                 ele_val(particle_flux(g_dash)%ptr,vele))
          
               end if not_within_group_eig
            
               ! add the spectrum production - this is the eigenvector
               rhs_addto = rhs_addto + matmul(shape_shape(ele_shape(particle_flux(g)%ptr,vele), &
                                                          ele_shape(particle_flux(g_dash)%ptr,vele), &
                                                          detwei_vele*ele_val_at_quad(production_coeff,vele) &
                                                                     *ele_val_at_quad(prompt_spectrum_coeff,vele) &
                                                                     *(1.0/particle%keff%keff_new) ), &
                                              ele_val(particle_flux_iter(g_dash)%ptr,vele))
                                          
            else keff_or_time

               ! add the scatter - not the within group
               not_within_group_time: if (g /= g_dash) then

                  rhs_addto = rhs_addto + matmul(shape_shape(ele_shape(particle_flux(g)%ptr,vele), &
                                                             ele_shape(particle_flux(g_dash)%ptr,vele), &
                                                             detwei_vele*ele_val_at_quad(scatter_coeff,vele)), &
                                                 (ele_val(particle_flux(g_dash)%ptr,vele)*theta + &
                                                  ele_val(particle_flux_old(g_dash)%ptr,vele)*(1.0-theta)))

               end if not_within_group_time

               ! add the spectrum production 
               rhs_addto = rhs_addto + matmul(shape_shape(ele_shape(particle_flux(g)%ptr,vele), &
                                                          ele_shape(particle_flux(g_dash)%ptr,vele), &
                                                          detwei_vele*ele_val_at_quad(production_coeff,vele) &
                                                                     *ele_val_at_quad(prompt_spectrum_coeff,vele)), &
                                              (ele_val(particle_flux_iter(g_dash)%ptr,vele)*theta + &
                                               ele_val(particle_flux_old(g_dash)%ptr,vele)*(1.0-theta)))

            end if keff_or_time
            
            ! add contribution to discretised source for this vele for this group g
            call addto(extra_discretised_source, &
                       ele_nodes(extra_discretised_source, vele), &
                       rhs_addto)
            
            deallocate(detwei_vele)
            deallocate(rhs_addto)
            
         end do vele_loop
               
      end do group_loop
      
      ! deallocate the local fields as needed
            
      call deallocate(scatter_coeff)
      call deallocate(production_coeff)
      call deallocate(prompt_spectrum_coeff)
      
      ! deallocate the particle flux pointer fields
      if (associated(particle_flux)) deallocate(particle_flux)
      if (associated(particle_flux_old)) deallocate(particle_flux_old)
      if (associated(particle_flux_iter)) deallocate(particle_flux_iter)
                        
   end subroutine assemble_coeff_source_group_g
   
   ! --------------------------------------------------------------------------

end module radiation_assemble_solve_group
