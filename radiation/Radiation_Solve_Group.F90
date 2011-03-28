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

   use sparse_tools
   use fields

   use boundary_conditions
   use boundary_conditions_from_options
   use field_options   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use petsc_solve_state_module   
   use sparsity_patterns_meshes 
    
   use radiation_particle_data_type
   use radiation_materials_interpolation
   use radiation_extract_flux_field
   use radiation_energy_group_set_tools
     
   implicit none
   
   private 

   public :: particle_assemble_solve_group

  ! Source?
  logical :: have_source
  ! Absorption?
  logical :: have_absorption
  ! Diffusivity?
  logical :: have_diffusivity
  ! Discretised source?
  logical :: have_discretised_source

contains

   ! --------------------------------------------------------------------------

   subroutine particle_assemble_solve_group(particle, &
                                            g, &
                                            invoke_eigenvalue_group_solve) 
   
      !!< Assemble and solve the group g particle matrix problem

      type(particle_type), intent(inout) :: particle
      integer, intent(in) :: g
      logical, intent(in) :: invoke_eigenvalue_group_solve
      
      ! local variables
      character(len=OPTION_PATH_LEN) :: field_name
            
      ! assemble the diffusivity, absorption and already discretised rhs term 
      ! scalar fields and insert into particle%state
      call assemble_coeff_source_group_g(particle, &
                                         g, &
                                         invoke_eigenvalue_group_solve)
      
      ! form the field name to solve for
      field_name = 'ParticleFluxGroup'//int2str(g)//'Moment1'//trim(particle%name)
            
      call assemble_matrix_solve_group_g(field_name, &
                                         particle%state)
            
   end subroutine particle_assemble_solve_group

   ! --------------------------------------------------------------------------
   
   subroutine assemble_coeff_source_group_g(particle, &
                                            g, &
                                            invoke_eigenvalue_group_solve)
      
      !!< Assemble the coeff and discretised source fields for group g that will then be used 
      !!< somewhere else to assemble the linear system

      type(particle_type), intent(inout) :: particle
      integer, intent(in) :: g
      logical, intent(in) :: invoke_eigenvalue_group_solve
      
      ! local variables
      integer :: stat
      integer :: vele
      integer :: inode
      integer :: g_dash
      integer :: g_set      
      real :: data_value
      real, dimension(:), allocatable :: detwei_vele
      real, dimension(:), allocatable :: rhs_addto      
      type(scalar_field) :: production_coeff 
      type(scalar_field) :: prompt_spectrum_coeff 
      type(scalar_field) :: scatter_coeff            
      type(scalar_field), pointer :: absorption_coeff
      type(tensor_field), pointer :: diffusivity_coeff
      type(scalar_field), pointer :: discretised_source
      type(mesh_type), pointer :: material_fn_space
      type(vector_field), pointer :: positions      
      type(scalar_field_pointer), dimension(:), pointer :: particle_flux 
      type(scalar_field_pointer), dimension(:), pointer :: particle_flux_old      
      character(len=OPTION_PATH_LEN) :: material_fn_space_name
      character(len=OPTION_PATH_LEN) :: energy_group_set_path
      character(len=OPTION_PATH_LEN) :: field_name
                 
      ! extract the particle flux fields for all energy groups
      call extract_flux_all_group(particle, & 
                                  particle_flux     = particle_flux, &
                                  particle_flux_old = particle_flux_old)

      ! currently assume all energy groups point to the same coordinate mesh
      
      ! if the above assumption was to be broken then the there would need to be a sweep of the other groups
      ! fluxes (new and old) to see which needs supermeshing to a new field and then the pointer within the 
      ! scalar_field_pointer particle_flux (or _old) is then changed such that no code change below is needed
      
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

      ! allocate the local fields as needed

      field_name = trim(particle_flux(g)%ptr%name)//'Scatter'
      call allocate(scatter_coeff, material_fn_space, trim(field_name))
      call zero(scatter_coeff)

      field_name = trim(particle_flux(g)%ptr%name)//'Production'
      call allocate(production_coeff, material_fn_space, trim(field_name))
      call zero(production_coeff)

      field_name = trim(particle_flux(g)%ptr%name)//'PromptSpectrum'
      call allocate(prompt_spectrum_coeff, material_fn_space, trim(field_name))
      call zero(prompt_spectrum_coeff)
      
      ! extract the assemble fields as needed
      absorption_coeff => extract_scalar_field(particle%state, &
                                               trim(particle_flux(g)%ptr%name) // 'Absorption', &
                                               stat = stat)
                                        
      diffusivity_coeff => extract_tensor_field(particle%state, &
                                                trim(particle_flux(g)%ptr%name) // 'Diffusivity', &
                                                stat = stat)

      discretised_source => extract_scalar_field(particle%state, &
                                                 trim(particle_flux(g)%ptr%name) // 'DiscretisedSource', &
                                                 stat = stat)
      
      ! zero the extracted fields
      call zero(absorption_coeff)
      call zero(diffusivity_coeff)
      call zero(discretised_source)
            
      ! form the absorption coeff field for this energy group (which is actually the removal cross section)            
      call form(material_fn_space, &
                particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                particle%particle_radmat, &
                absorption_coeff, &
                g, &
                component = 'removal')
            
      ! form the diffusivity tensor coeff field for this energy group (only fill in diagonals)            
      node_loop_d: do inode = 1,node_count(material_fn_space)                                                                      
            
         ! get the inode diffusionx data for this energy group
         call form(particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                   particle%particle_radmat, &
                   inode, &
                   g, &
                   data_value, &
                   component = 'diffusionx')
               
         ! set the interpolated value into the scalar field
         call set(diffusivity_coeff, &
                  1, &
                  1, &
                  inode, &
                  data_value)

         if (positions%dim > 1) then
         
            ! get the inode diffusiony data for this energy group
            call form(particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                      particle%particle_radmat, &
                      inode, &
                      g, &
                      data_value, &
                      component = 'diffusiony')
               
            ! set the interpolated value into the scalar field
            call set(diffusivity_coeff, &
                     2, &
                     2, &
                     inode, &
                     data_value)

         end if 
         
         if (positions%dim > 2) then

            ! get the inode diffusionz data for this energy group
            call form(particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                      particle%particle_radmat, &
                      inode, &
                      g, &
                      data_value, &
                      component = 'diffusionz')
               
            ! set the interpolated value into the scalar field
            call set(diffusivity_coeff, &
                     3, &
                     3, &
                     inode, &
                     data_value)

         end if 
               
      end do node_loop_d
            
      ! form the discretised source - scatter and production            
      
      ! form the prompt spectrum coeff field for this energy group
      call form(material_fn_space, &
                particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                particle%particle_radmat, &
                prompt_spectrum_coeff, &
                g, &
                component = 'prompt_spectrum')
                  
      group_loop: do g_dash = 1,size(particle_flux)

         ! form the production coeff field for this energy group
         call form(material_fn_space, &
                   particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                   particle%particle_radmat, &
                   production_coeff, &
                   g_dash, &
                   component = 'production')
         
         not_within_group: if (g /= g_dash) then
 
            ! form the scatter field for this energy group g_dash to g            
            call form(material_fn_space, &
                      particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                      particle%particle_radmat, &
                      scatter_coeff, &
                      g, &
                      component = 'scatter', &
                      g_dash = g_dash)
                                           
         end if not_within_group
         
         vele_loop: do vele = 1,ele_count(discretised_source)
           
            ! allocate the jacobian transform and gauss weight array for this vele
            allocate(detwei_vele(ele_ngi(discretised_source,vele)))
                        
            ! form the velement jacobian transform and gauss weight
            call transform_to_physical(positions, vele, detwei = detwei_vele)
            
            ! allocate the rhs_addto
            allocate(rhs_addto(ele_loc(discretised_source,vele)))
            
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
                                                          ele_shape(particle_flux_old(g_dash)%ptr,vele), &
                                                          detwei_vele*ele_val_at_quad(production_coeff,vele) &
                                                                     *ele_val_at_quad(prompt_spectrum_coeff,vele) &
                                                                     *(1.0/particle%keff%keff_new) ), &
                                              ele_val(particle_flux_old(g_dash)%ptr,vele))
                                          
            else keff_or_time
  
               ! add the scatter - not the within group
               not_within_group_time: if (g /= g_dash) then
                  
                  ! fill in ...
                  
               end if not_within_group_time

               ! add the spectrum production 

               ! fill in ...

            end if keff_or_time
            
            ! add contribution to discretised source for this vele for this group g
            call addto(discretised_source, &
                       ele_nodes(discretised_source, vele), &
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
      call deallocate_flux_all_group(particle_flux     = particle_flux, &
                                     particle_flux_old = particle_flux_old)
                        
   end subroutine assemble_coeff_source_group_g
   
   ! --------------------------------------------------------------------------
   
   subroutine assemble_matrix_solve_group_g(field_name, &
                                            state)
      
      !!< Assemble the matrix system for this group g and solve
      !!< This is only set up for even parity spherical harmonic P1 (ie diffusion) 
            
      character(len=*), intent(in) :: field_name
      type(state_type), intent(inout) :: state
      
      ! local variables
      integer :: stat
      integer :: vele,sele
      type(csr_matrix) :: matrix
      type(scalar_field) :: rhs
      type(csr_sparsity), pointer :: sparsity
      integer, dimension(:), allocatable :: bc_types
      type(scalar_field) :: bc_value
      type(tensor_field), pointer :: diffusivity
      type(scalar_field), pointer :: absorption
      type(scalar_field), pointer :: source      
      type(scalar_field), pointer :: discretised_source
      type(vector_field), pointer :: positions      
      type(scalar_field), pointer :: t 
        
      t => extract_scalar_field(state, field_name)
      
      if(t%mesh%continuity /= 0) then
         FLExit('Radiation model requires a continuous solution mesh')
      end if
      
      ! determine the sparsity pattern of matrix assuming first order connections 
      sparsity => get_csr_sparsity_firstorder(state, &
                                              t%mesh, &
                                              t%mesh)
      
      ! allocate the matrix and rhs vector used for solving
      call allocate(matrix, &
                    sparsity)
                    
      call allocate(rhs, &
                    t%mesh)
      
      call zero(matrix)
      call zero(rhs)

      ! get the positions field for this energy group set
      positions => extract_vector_field(state, &
                                        'Coordinate')
      
      ! extract the absorption, diffusivity, source and discretised source fields
      absorption => extract_scalar_field(state, &
                                         trim(t%name) // 'Absorption', &
                                         stat = stat)

      have_absorption = stat == 0         
                                        
      diffusivity => extract_tensor_field(state, &
                                          trim(t%name) // 'Diffusivity', &
                                          stat = stat)

      have_diffusivity = stat == 0         

      source => extract_scalar_field(state, &
                                     trim(t%name)//'Source', &
                                     stat=stat)
             
      have_source = stat == 0         

      discretised_source => extract_scalar_field(state, &
                                                 trim(t%name) // 'DiscretisedSource', &
                                                 stat = stat)

      have_discretised_source = stat == 0         
  
      ! volume integrations 
      volume_element_loop: do vele = 1, ele_count(t)      
                  
         call assemble_particle_group_g_vele(vele, &
                                             t, &
                                             matrix, &
                                             rhs, &
                                             positions, &
                                             source, &
                                             absorption, &
                                             diffusivity)

      end do volume_element_loop
      
      ! include the already discretised_source into rhs
      if (have_discretised_source) call addto(rhs, discretised_source)
      
      ! surface integrations for albedo and source BC's for diffusion theory
      
      allocate(bc_types(surface_element_count(t)))

      call get_entire_boundary_condition(t, &
                                         (/"albedo", &
                                           "source"/), &
                                         bc_value, &
                                         bc_types)

      surface_element_loop: do sele = 1,surface_element_count(t)
         
         diffusion_bc: if (bc_types(sele) == 1) then
         
            call assemble_particle_group_g_sele_albedo(sele, &
                                                       positions, &
                                                       matrix, & 
                                                       bc_value, &
                                                       t)
         
         else if (bc_types(sele) == 2) then 
            
            call assemble_particle_group_g_sele_source(sele, &
                                                       positions, &
                                                       rhs, & 
                                                       bc_value, &
                                                       t)
                        
         end if diffusion_bc
         
      end do surface_element_loop

      call deallocate(bc_value)

      if (allocated(bc_types)) deallocate(bc_types)
            
      ! apply direchlet BC
      call apply_dirichlet_conditions(matrix, &
                                      rhs, &
                                      t)
            
      ! solve the linear system       
      call petsc_solve(t, &
                       matrix, &
                       rhs, &
                       state, &
                       option_path = trim(t%option_path) )
      
      ! set the direchlet bc nodes to be consistent
      call set_dirichlet_consistent(t)
   
      ! deallocate the matrix and rhs vector used for solving
      call deallocate(matrix)
      call deallocate(rhs)
              
   end subroutine assemble_matrix_solve_group_g

   ! --------------------------------------------------------------------------

   subroutine assemble_particle_group_g_vele(vele, &
                                             t, &
                                             matrix, &
                                             rhs, &
                                             positions, &
                                             source, &
                                             absorption, &
                                             diffusivity)

      !!< Assemble the volume element vele contribution to the global 
      !!< matrix and rhs for this particle type for group g
            
      integer, intent(in) :: vele
      type(scalar_field), intent(in) :: t 
      type(csr_matrix), intent(inout) :: matrix
      type(scalar_field), intent(inout) :: rhs
      type(vector_field), intent(in) :: positions
      type(scalar_field), intent(in) :: source      
      type(scalar_field), intent(in) :: absorption
      type(tensor_field), intent(in) :: diffusivity
      
      ! local variables 
      real, dimension(ele_ngi(t, vele)) :: detwei
      real, dimension(ele_loc(t, vele), ele_ngi(t, vele), mesh_dim(t)) :: dshape
      real, dimension(ele_loc(t, vele), ele_loc(t, vele)) :: matrix_addto
      real, dimension(ele_loc(t, vele)) :: rhs_addto      
      
      ! initialise the local vele matrix and vector that are added to the global arrays
      matrix_addto = 0.0
      rhs_addto    = 0.0
           
      ! form the velement jacobian transform and gauss weight
      call transform_to_physical(positions, &
                                 vele, &
                                 ele_shape(t, vele), &
                                 dshape = dshape, &
                                 detwei = detwei)
            
      ! add the diffusion term            
      have_diffusivity_if: if (have_diffusivity) then
         
         matrix_addto = matrix_addto + dshape_tensor_dshape(dshape, &
                                                            ele_val_at_quad(diffusivity, vele), &
                                                            dshape, &
                                                            detwei)

      end if have_diffusivity_if
               
      ! add the absortion term 
      have_absorption_if: if (have_absorption) then
         
         matrix_addto = matrix_addto + shape_shape(ele_shape(t,vele), &
                                                   ele_shape(t,vele), &
                                                   detwei*ele_val_at_quad(absorption,vele))
         
      end if have_absorption_if 

      ! add the source term 
      have_source_if: if (have_source) then
         
             rhs_addto = rhs_addto + shape_rhs(ele_shape(t,vele), &
                                               detwei*ele_val_at_quad(source,vele))
         
      end if have_source_if 

      ! insert the local vele matrix and rhs into the global arrays 
      call addto(matrix, &
                 ele_nodes(t, vele), &
                 ele_nodes(t, vele), &
                 matrix_addto)

      call addto(rhs, &
                 ele_nodes(t, vele), &
                 rhs_addto)
              
   end subroutine assemble_particle_group_g_vele

   ! --------------------------------------------------------------------------

   subroutine assemble_particle_group_g_sele_albedo(face, &
                                                    positions, &
                                                    matrix, & 
                                                    bc_value, &
                                                    t)
      
      !!< Assemble the diffusion theory albedo boundary condition that involves adding 
      !!< a surface mass matrix * coefficient to the main matrix
      
      ! this is actually a robin BC for angular diffusion
      
      integer, intent(in) :: face
      type(vector_field), intent(in) :: positions
      type(csr_matrix), intent(inout) :: matrix
      type(scalar_field), intent(in) :: bc_value
      type(scalar_field), pointer :: t 
      
      ! local variables
      real, dimension(face_ngi(t, face)) :: detwei    
      real, dimension(face_loc(t, face), face_loc(t, face)) :: matrix_addto
      
      assert(face_ngi(positions, face) == face_ngi(t, face))
            
      call transform_facet_to_physical(positions, &
                                       face, &
                                       detwei_f = detwei)  
      
      ! the 0.5 is a necessary part of this BC formulation, the bc_value is the albedo coeff
      ! where a value of 1.0 is a perfect reflective and a value of 0.0 is a perfect vacuum
      matrix_addto = shape_shape(face_shape(t, face), &
                                 face_shape(t, face), &
                                 detwei*0.5*( (1.0 - ele_val_at_quad(bc_value,face) ) / &
                                              (1.0 + ele_val_at_quad(bc_value,face) ) ) )
      
      call addto(matrix, &
                 face_global_nodes(t,face), &
                 face_global_nodes(t,face), &
                 matrix_addto)
                
   end subroutine assemble_particle_group_g_sele_albedo

   ! --------------------------------------------------------------------------

   subroutine assemble_particle_group_g_sele_source(face, &
                                                    positions, &
                                                    rhs, & 
                                                    bc_value, &
                                                    t)
      
      !!< Include the particle source boundary condition into the rhs
      
      ! this is actually a neumann BC for angular diffusion
      
      integer, intent(in) :: face
      type(vector_field), intent(in) :: positions
      type(scalar_field), intent(inout) :: rhs
      type(scalar_field), intent(in) :: bc_value
      type(scalar_field), pointer :: t 
      
      ! local variables
      real, dimension(face_ngi(t, face)) :: detwei    
      real, dimension(ele_loc(t, face)) :: rhs_addto      
      
      assert(face_ngi(positions, face) == face_ngi(t, face))
            
      call transform_facet_to_physical(positions, &
                                       face, &
                                       detwei_f = detwei)  
      
      rhs_addto = shape_rhs(face_shape(t, face), &
                            detwei*ele_val_at_quad(bc_value,face) )
      
      call addto(rhs, &
                 face_global_nodes(t,face), &
                 rhs_addto)
                
   end subroutine assemble_particle_group_g_sele_source

   ! --------------------------------------------------------------------------

end module radiation_assemble_solve_group
