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

module radiation_materials_interpolation_form_radmat

   !!< Procedures associated with the forming of the interpolated radiation material set for 
   !!< a particular material mesh inode from the already determined interpolation instructions.

   use futils   
   use global_parameters, only : OPTION_PATH_LEN
   use fields
  
   use radiation_materials  
   use radiation_materials_interpolation_data_types

   implicit none
   
   private 

   public :: form
   
   interface form
      module procedure form_radiation_material_scalar_field, &
                       form_radmat_inode_component_group_g
   end interface form
      
contains

   ! --------------------------------------------------------------------------

   subroutine form_radiation_material_scalar_field(material_fn_space, &
                                                   energy_group_set_ii, &
                                                   particle_radmat, &
                                                   material_scalar_field, &                                                   
                                                   g, &
                                                   component, &
                                                   g_dash) 
      
      !!< Form a radiation material scalar field for the component for group g (and perhaps associated with g_dash)
      
      type(mesh_type), intent(in) :: material_fn_space
      type(energy_group_set_ii_type), intent(in) :: energy_group_set_ii  
      type(particle_radmat_type), intent(in) :: particle_radmat  
      type(scalar_field), intent(inout) :: material_scalar_field
      integer, intent(in) :: g
      character(len=*), intent(in) :: component
      integer, intent(in), optional :: g_dash
      
      ! local variables
      integer :: inode
      real :: data_value
      
      node_loop_a: do inode = 1,node_count(material_fn_space)                                                                      
            
         ! get the inode material data for this energy group at this spatial dof
         call form(energy_group_set_ii, &
                   particle_radmat, &
                   inode, &
                   g, &
                   data_value, &
                   component, &
                   g_dash = g_dash)
               
         ! set the interpolated value into the scalar field
         call set(material_scalar_field, &
                  inode, &
                  data_value)
                        
      end do node_loop_a      
      
   end subroutine form_radiation_material_scalar_field
   
   ! --------------------------------------------------------------------------

   subroutine form_radmat_inode_component_group_g(energy_group_set_ii, &
                                                  particle_radmat, &
                                                  inode, &
                                                  g, &
                                                  value, &
                                                  component, &
                                                  g_dash)
      
      !!< Form the interpolated material data of a particular component for this inode for group g 
      !!< For scatter it is the data from g_dash to g
   
      type(energy_group_set_ii_type), intent(in) :: energy_group_set_ii  
      type(particle_radmat_type), intent(in) :: particle_radmat 
      integer, intent(in) :: inode 
      integer, intent(in) :: g
      real, intent(inout) :: value 
      character(len=*), intent(in) :: component  
      integer, intent(in), optional :: g_dash   
      
      ! local variables 
      integer :: dmat
      integer :: pmat
      
      ! initialise values
      value = 0.0
      
      ! use the region id ii
      region_id: if (allocated(energy_group_set_ii%region_id_ii)) then
         
         dmat = energy_group_set_ii%region_id_ii(inode)%dataset_ii%dataset_radmat_number 
         
         pmat = energy_group_set_ii%region_id_ii(inode)%dataset_ii%physical_radmat_ii%physical_radmat_number
                  
         call form_radmat_inode_component_group_g_base(energy_group_set_ii%region_id_ii(inode)%dataset_ii%physical_radmat_ii, &
                                                       particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat), &
                                                       g, &
                                                       value, &
                                                       component, &
                                                       g_dash = g_dash)

      end if region_id       
   
   end subroutine form_radmat_inode_component_group_g
   
   ! --------------------------------------------------------------------------
   
   subroutine form_radmat_inode_component_group_g_base(physical_radmat_ii, &
                                                       physical_radmat, &
                                                       g, &
                                                       value, &
                                                       component, &
                                                       g_dash)
      
      !!< Add an interpolated value to the component term for group g
      !!< For scatter it is the data from g_dash to g
      
      type(physical_radmat_ii_type), intent(in) :: physical_radmat_ii
      type(physical_radmat_type), intent(in) :: physical_radmat
      integer, intent(in) :: g
      real, intent(inout) :: value
      character(len=*), intent(in) :: component     
      integer, intent(in), optional :: g_dash   

      ! hard wire in 1d interpolation :(
      one_dim: if (size(physical_radmat_ii%fraction) == 1) then

         one_radmats: if (size(physical_radmat%radmats) == 1) then
                           
            if (trim(component) == 'production') then
            
               value = value + physical_radmat%radmats(1)%production(g)

            else if (trim(component) == 'fission') then
            
               value = value + physical_radmat%radmats(1)%fission(g)

            else if (trim(component) == 'power') then
            
               value = value + physical_radmat%radmats(1)%power(g)

            else if (trim(component) == 'velocity') then
            
               value = value + physical_radmat%radmats(1)%velocity(g)

            else if (trim(component) == 'removal') then
            
               value = value + physical_radmat%radmats(1)%removal(g,1)

            else if (trim(component) == 'diffusionx') then
            
               value = value + physical_radmat%radmats(1)%diffusion(g,1)

            else if (trim(component) == 'diffusiony') then
            
               value = value + physical_radmat%radmats(1)%diffusion(g,2)

            else if (trim(component) == 'diffusionz') then
            
               value = value + physical_radmat%radmats(1)%diffusion(g,3)

            else if (trim(component) == 'prompt_spectrum') then
            
               value = value + physical_radmat%radmats(1)%prompt_spectrum(g)

            else if (trim(component) == 'scatter') then
            
               value = value + physical_radmat%radmats(1)%scatter(g_dash,g,1)
            
            end if 
               
         else one_radmats

            if (trim(component) == 'production') then
                              
               value = value + &
                       (1.0 - physical_radmat_ii%fraction(1))* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1))%production(g) + &
                       physical_radmat_ii%fraction(1)* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1) + 1)%production(g)

            else if (trim(component) == 'fission') then
                              
               value = value + &
                       (1.0 - physical_radmat_ii%fraction(1))* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1))%fission(g) + &
                       physical_radmat_ii%fraction(1)* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1) + 1)%fission(g)

            else if (trim(component) == 'power') then
                              
               value = value + &
                       (1.0 - physical_radmat_ii%fraction(1))* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1))%power(g) + &
                       physical_radmat_ii%fraction(1)* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1) + 1)%power(g)

            else if (trim(component) == 'velocity') then
                              
               value = value + &
                       (1.0 - physical_radmat_ii%fraction(1))* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1))%velocity(g) + &
                       physical_radmat_ii%fraction(1)* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1) + 1)%velocity(g)

            else if (trim(component) == 'removal') then
                              
               value = value + &
                       (1.0 - physical_radmat_ii%fraction(1))* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1))%removal(g,1) + &
                       physical_radmat_ii%fraction(1)* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1) + 1)%removal(g,1)

            else if (trim(component) == 'diffusionx') then
                              
               value = value + &
                       (1.0 - physical_radmat_ii%fraction(1))* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1))%diffusion(g,1) + &
                       physical_radmat_ii%fraction(1)* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1) + 1)%diffusion(g,1)

            else if (trim(component) == 'diffusiony') then
                              
               value = value + &
                       (1.0 - physical_radmat_ii%fraction(1))* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1))%diffusion(g,2) + &
                       physical_radmat_ii%fraction(1)* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1) + 1)%diffusion(g,2)


            else if (trim(component) == 'diffusionz') then
                              
               value = value + &
                       (1.0 - physical_radmat_ii%fraction(1))* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1))%diffusion(g,3) + &
                       physical_radmat_ii%fraction(1)* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1) + 1)%diffusion(g,3)

            else if (trim(component) == 'prompt_spectrum') then
                              
               value = value + &
                       (1.0 - physical_radmat_ii%fraction(1))* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1))%prompt_spectrum(g) + &
                       physical_radmat_ii%fraction(1)* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1) + 1)%prompt_spectrum(g)

            else if (trim(component) == 'scatter') then
                              
               value = value + &
                       (1.0 - physical_radmat_ii%fraction(1))* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1))%scatter(g_dash,g,1) + &
                       physical_radmat_ii%fraction(1)* &
                       physical_radmat%radmats(physical_radmat_ii%radmat_base_coordinate(1) + 1)%scatter(g_dash,g,1)

            end if 
            
         end if one_radmats       
      
      end if one_dim              
                                                 
   end subroutine form_radmat_inode_component_group_g_base
   
   ! --------------------------------------------------------------------------

end module radiation_materials_interpolation_form_radmat
