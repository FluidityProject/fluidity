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

module radiation_calculate_current_module

   !!< This module contains procedures associated with
   !!< calculating the radiation particle current

   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use fields
   use vector_tools      
   
   use radiation_particle_data_type
   use radiation_materials_interpolation
   use radiation_extract_flux_field
   use radiation_energy_group_set_tools
     
   implicit none
   
   private 

   public :: radiation_calculate_current

contains

   ! --------------------------------------------------------------------------

   subroutine radiation_calculate_current(particle) 
   
      !!< Calculate the radiation particle current that is
      !!< assumed on a discontinuous mesh
      
      ! this currently assumes the angular discretisation is diffusion theory
      
      type(particle_type), intent(inout) :: particle
      
      ! local variables
      integer :: stat
      integer :: idim1, idim2
      integer :: g
      integer :: g_set
      integer :: vele
      integer :: number_of_energy_groups
      real, dimension(:), allocatable :: detwei_vele
      real, dimension(:,:,:), allocatable :: dshape_vele
      real, dimension(:), allocatable :: rhs_vele
      real, dimension(:,:,:), allocatable :: rhs_matrix_vele
      real, dimension(:,:), allocatable :: mass_matrix_vele
      type(tensor_field), pointer :: diffusivity_coeff
      type(vector_field), pointer :: positions      
      type(scalar_field), pointer :: particle_flux
      type(vector_field), pointer :: particle_current
      character(len=OPTION_PATH_LEN) :: current_path

      ! get the number of energy groups via summing the number within each energy group set      
      call find_total_number_energy_groups(trim(particle%option_path), &
                                           number_of_energy_groups)      
      
      group_loop: do g = 1,number_of_energy_groups

         ! determine which group set this g belongs to
         call which_group_set_contains_g(g, &
                                         trim(particle%option_path), &
                                         g_set)
         
         ! form the current path for this energy group set
         current_path = trim(particle%option_path)//'/energy_discretisation/energy_group_set['//int2str(g_set - 1)//']/angular_discretisation/method/parity/vector_field::ParticleCurrent'

         ! if current path not exist cycle this group
         if (.not. have_option(trim(current_path))) cycle group_loop

         ! extract the particle flux field for this energy group
         call extract_flux_group_g(particle, & 
                                   g, &
                                   particle_flux = particle_flux)
         
         ! extract the particle current field for this energy group
         call extract_current_group_g(particle, & 
                                      g, &
                                      particle_current = particle_current)

         ! get the positions field for this energy group set
         ! - assume all the same positions mesh for now
         positions => extract_vector_field(particle%state, &
                                           'Coordinate')
         
         ! assert that the particle current mesh continuity is discontinuous
         assert(continuity(particle_current)<0)
         
         ! extract the diffusivity field 
         diffusivity_coeff => extract_tensor_field(particle%state, &
                                                   trim(particle_flux%name) // 'Diffusivity', &
                                                   stat = stat)

         ! zero the diffusivity field
         call zero(diffusivity_coeff)

         ! form the diffusivity tensor coeff field for this energy group
         call form(particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                   particle%particle_radmat, &
                   diffusivity_coeff, &
                   g, &
                   positions%dim)
         
         vele_loop: do vele = 1,ele_count(positions)

            ! allocate the jacobian transform and gauss weight array for this vele
            allocate(detwei_vele(ele_ngi(positions,vele)))
            
            ! allocate the shape derivative for the particle flux
            allocate(dshape_vele(ele_loc(particle_flux,vele), &
                                 ele_ngi(particle_flux,vele), &
                                 mesh_dim(positions))) 
            
            ! form the velement jacobian transform and gauss weight
            call transform_to_physical(positions, &
                                       vele, &
                                       ele_shape(particle_flux,vele), &
                                       dshape = dshape_vele, &
                                       detwei = detwei_vele)
                              
            ! allocate the mass_matrix_vele
            allocate(mass_matrix_vele(ele_loc(particle_current,vele), &
                                      ele_loc(particle_current,vele)))

            ! form the mass matrix values for this vele      
            mass_matrix_vele = shape_shape(ele_shape(particle_current,vele), &
                                           ele_shape(particle_current,vele), &
                                           detwei_vele)

            ! form the mass inverse into the same variable
            call invert(mass_matrix_vele)
            
            ! allocate the rhs_vele
            allocate(rhs_vele(ele_loc(particle_current,vele)))

            ! allocate the rhs_matrix_vele
            allocate(rhs_matrix_vele(mesh_dim(positions), &
                                     ele_loc(particle_current,vele), &
                                     ele_loc(particle_flux,vele)))
            
            ! find the current in each coordinate direction
            dim_loop1: do idim1 = 1,mesh_dim(positions)
            
               dim_loop2: do idim2 = 1,mesh_dim(positions)
                        
                  ! form the rhs matrix values for this vele      
                  rhs_matrix_vele(idim2:idim2,:,:) = - shape_dshape(ele_shape(particle_current,vele), &
                                                                    dshape_vele(:,:,idim2:idim2), &      
                                                                    detwei_vele*ele_val_at_quad(diffusivity_coeff, &
                                                                                                idim1, &
                                                                                                idim2, &
                                                                                                vele))

               end do dim_loop2
               
               ! form the rhs_vele
               rhs_vele = 0.0
               
               dim_loop3: do idim2 = 1,mesh_dim(positions)
               
                  rhs_vele = rhs_vele + matmul(rhs_matrix_vele(idim2,:,:), &
                                               ele_val(particle_flux,vele))
                  
               end do dim_loop3
                        
               call set(particle_current, &
                        idim1, &
                        ele_nodes(particle_current, vele), &
                        matmul(mass_matrix_vele, rhs_vele))

            end do dim_loop1
            
            deallocate(detwei_vele)
            deallocate(dshape_vele)
            deallocate(rhs_vele)
            deallocate(rhs_matrix_vele)
            deallocate(mass_matrix_vele)
            
         end do vele_loop
               
      end do group_loop
      
   end subroutine radiation_calculate_current

   ! --------------------------------------------------------------------------

end module radiation_calculate_current_module 
