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
   use spud
   use fields
   use state_module
     
   use radiation_materials  
   use radiation_materials_interpolation_data_types

   implicit none
   
   private 

   public :: form
   
   interface form
      module procedure form_radiation_material_diffusivity_field, &
                       form_radiation_material_scalar_field, &
                       form_radmat_inode_component_group_g
   end interface form
      
contains

   ! --------------------------------------------------------------------------
   
   subroutine form_radiation_material_diffusivity_field(state, &
                                                        energy_group_set_ii, &
                                                        particle_radmat, &
                                                        diffusivity_field, &                                                   
                                                        g, &
                                                        geom_dim)

      !!< Form a radiation material diffusivity field for group g

      type(state_type), intent(in) :: state
      type(energy_group_set_ii_type), intent(in) :: energy_group_set_ii  
      type(particle_radmat_type), intent(in) :: particle_radmat  
      type(tensor_field), intent(inout) :: diffusivity_field
      integer, intent(in) :: g
      integer, intent(in) :: geom_dim
      
      ! local variables
      integer :: status
      integer :: inode,r
      integer :: number_of_rotations
      real :: cos_angle_node_val, sin_angle_node_val
      real, dimension(:,:), allocatable :: diff_node_values
      real, dimension(:,:), allocatable :: rotation_matrix
      character(len=OPTION_PATH_LEN) :: rotation_path
      character(len=OPTION_PATH_LEN) :: field_name
      character, dimension(:), allocatable :: rotation_dimension
      type(scalar_field_pointer), dimension(:), pointer :: angle
      
      number_of_rotations = 0
      
      ! determine if the diffusivity tensor is to be rotated
      not_1d: if (geom_dim > 1) then
      
         rotation_path = trim(diffusivity_field%option_path)//'/diagnostic/rotation'
      
         number_of_rotations = option_count(trim(rotation_path))
      
         rotations_extract: if (number_of_rotations > 0) then
         
            allocate(rotation_matrix(geom_dim,geom_dim))

            allocate(rotation_dimension(geom_dim))
         
            allocate(angle(number_of_rotations))
         
            rotation_extract_loop: do r = 1,number_of_rotations
            
               allocate(angle(r)%ptr)

               field_name = trim(diffusivity_field%name)//'Rotation'//int2str(r)

               angle(r)%ptr => extract_scalar_field(state, &
                                                    trim(field_name), &
                                                    stat=status)

               if (status /= 0) FLAbort('Failed to extract angle rotation field from state')
            
               ! find the rotation dimension (or axis) for this angle field
               call get_option(trim(diffusivity_field%option_path)//'/diagnostic/rotation['//int2str(r-1)//']/axis/name', &
                               rotation_dimension(r))
            
               only_Z_2d: if ((trim(rotation_dimension(r)) /= 'Z') .and. (geom_dim == 2)) then
               
                  FLAbort('For 2 dimension geometry only a Z axis rotation is appropriate')
               
               end if only_Z_2d
            
            end do rotation_extract_loop
                        
         end if rotations_extract
      
      end if not_1d
            
      allocate(diff_node_values(geom_dim,geom_dim))
      
      node_loop: do inode = 1,node_count(diffusivity_field)                                                                      
         
         diff_node_values = 0.0
                  
         ! get the inode diffusionx data for this energy group
         call form(energy_group_set_ii, &
                   particle_radmat, &
                   inode, &
                   g, &
                   diff_node_values(1,1), &
                   component = 'diffusionx')

         second: if (geom_dim > 1) then
         
            ! get the inode diffusiony data for this energy group
            call form(energy_group_set_ii, &
                      particle_radmat, &
                      inode, &
                      g, &
                      diff_node_values(2,2), &
                      component = 'diffusiony')

         end if second
         
         third: if (geom_dim > 2) then

            ! get the inode diffusionz data for this energy group
            call form(energy_group_set_ii, &
                      particle_radmat, &
                      inode, &
                      g, &
                      diff_node_values(3,3), &
                      component = 'diffusionz')

         end if third
         
         ! rotate the diffusivity tensor for this node if required
         
         rotations_apply: if (number_of_rotations > 0) then
         
            rotation_apply_loop: do r = 1,number_of_rotations
               
               cos_angle_node_val = cos(node_val(angle(r)%ptr,inode))
               sin_angle_node_val = sin(node_val(angle(r)%ptr,inode))
               
               dim: if (geom_dim == 2) then
                  
                  rotation_matrix(1,1) =  cos_angle_node_val
                  rotation_matrix(1,2) = -sin_angle_node_val
                  rotation_matrix(2,1) =  sin_angle_node_val
                  rotation_matrix(2,2) =  cos_angle_node_val
               
               else dim

                  axis: if (trim(rotation_dimension(r)) == 'X') then

                     rotation_matrix(1,1) =  1.0
                     rotation_matrix(1,2) =  0.0
                     rotation_matrix(1,3) =  0.0
                     rotation_matrix(2,1) =  0.0
                     rotation_matrix(2,2) =  cos_angle_node_val
                     rotation_matrix(2,3) = -sin_angle_node_val
                     rotation_matrix(3,1) =  0.0
                     rotation_matrix(3,2) =  sin_angle_node_val
                     rotation_matrix(3,3) =  cos_angle_node_val

                  else if (trim(rotation_dimension(r)) == 'Y') then axis

                     rotation_matrix(1,1) =  cos_angle_node_val
                     rotation_matrix(1,2) =  0.0
                     rotation_matrix(1,3) =  sin_angle_node_val
                     rotation_matrix(2,1) =  0.0
                     rotation_matrix(2,2) =  1.0
                     rotation_matrix(2,3) =  0.0
                     rotation_matrix(3,1) = -sin_angle_node_val
                     rotation_matrix(3,2) =  0.0
                     rotation_matrix(3,3) =  cos_angle_node_val

                  else if (trim(rotation_dimension(r)) == 'Z') then axis
                 
                     rotation_matrix(1,1) =  cos_angle_node_val
                     rotation_matrix(1,2) = -sin_angle_node_val
                     rotation_matrix(1,3) =  0.0
                     rotation_matrix(2,1) =  sin_angle_node_val
                     rotation_matrix(2,2) =  cos_angle_node_val
                     rotation_matrix(2,3) =  0.0
                     rotation_matrix(3,1) =  0.0
                     rotation_matrix(3,2) =  0.0
                     rotation_matrix(3,3) =  1.0
                 
                  end if axis
               
               end if dim
               
               diff_node_values = matmul(rotation_matrix,matmul(diff_node_values,transpose(rotation_matrix)))
               
            end do rotation_apply_loop
            
         end if rotations_apply
               
         ! set the interpolated node values into the tensor field
         call set(diffusivity_field, &
                  inode, &
                  diff_node_values)
               
      end do node_loop
      
      if (allocated(diff_node_values))   deallocate(diff_node_values)
      if (allocated(rotation_matrix))    deallocate(rotation_matrix)
      if (allocated(rotation_dimension)) deallocate(rotation_dimension)
      
      if (associated(angle)) deallocate(angle)
   
   end subroutine form_radiation_material_diffusivity_field
   
   ! --------------------------------------------------------------------------

   subroutine form_radiation_material_scalar_field(energy_group_set_ii, &
                                                   particle_radmat, &
                                                   material_scalar_field, &                                                   
                                                   g, &
                                                   component, &
                                                   g_dash) 
      
      !!< Form a radiation material scalar field for the component for group g 
      !!< (and perhaps associated with g_dash)
      
      type(energy_group_set_ii_type), intent(in) :: energy_group_set_ii  
      type(particle_radmat_type), intent(in) :: particle_radmat  
      type(scalar_field), intent(inout) :: material_scalar_field
      integer, intent(in) :: g
      character(len=*), intent(in) :: component
      integer, intent(in), optional :: g_dash
      
      ! local variables
      integer :: inode
      real :: data_value
      
      node_loop: do inode = 1,node_count(material_scalar_field)

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
                        
      end do node_loop      
      
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
      real, intent(out) :: value 
      character(len=*), intent(in) :: component  
      integer, intent(in), optional :: g_dash   
      
      ! local variables 
      integer :: dmat
      integer :: pmat
      
      ! initialise value
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
