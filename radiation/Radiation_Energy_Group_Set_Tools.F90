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

module radiation_energy_group_set_tools

   !!< This module contains procedures associated with tools applied to the energy group sets
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   
   implicit none
   
   private 

   public :: which_group_set_contains_g, &
             first_g_within_group_set, &
             find_total_number_energy_groups

contains

   ! --------------------------------------------------------------------------

   subroutine which_group_set_contains_g(g, &
                                         particle_option_path, &
                                         g_set) 
   
      !!< Determine which group set this g belongs to

      integer, intent(in) :: g
      character(len=*), intent(in) :: particle_option_path
      integer, intent(out) :: g_set
      
      ! local variables
      integer :: number_of_energy_group_set
      integer :: end_g
      integer :: g_set_index
      integer :: number_of_energy_groups_g_set
      character(len=OPTION_PATH_LEN) :: energy_group_set_path

      ! deduce the number of energy group sets
      number_of_energy_group_set = option_count(trim(particle_option_path)//'/energy_discretisation/energy_group_set')
      
      end_g = 1
      energy_group_set_loop: do g_set_index = 1,number_of_energy_group_set
            
         ! set the energy_group_set path
         energy_group_set_path = trim(particle_option_path)//'/energy_discretisation/energy_group_set['//int2str(g_set_index - 1)//']'
                                             
         ! get the number_energy_groups within this set
         call get_option(trim(energy_group_set_path)//'/number_of_energy_groups',number_of_energy_groups_g_set)         
                  
         end_g = end_g + number_of_energy_groups_g_set - 1 
         
         if (g <= end_g) then
         
            g_set = g_set_index
            
            exit energy_group_set_loop
         
         end if 
         
      end do energy_group_set_loop
                 
   end subroutine which_group_set_contains_g

   ! --------------------------------------------------------------------------

   subroutine first_g_within_group_set(g_set, &
                                       particle_option_path, &
                                       first_g_in_g_set)
      
      !!< Determine the first energy group in this group set

      integer, intent(in) :: g_set
      character(len=*), intent(in) :: particle_option_path
      integer, intent(out) :: first_g_in_g_set
      
      ! local variables
      integer :: number_of_energy_group_set
      integer :: g_set_index
      integer :: number_of_energy_groups_g_set
      character(len=OPTION_PATH_LEN) :: energy_group_set_path

      ! deduce the number of energy group sets
      number_of_energy_group_set = option_count(trim(particle_option_path)//'/energy_discretisation/energy_group_set')
      
      first_g_in_g_set = 1
      energy_group_set_loop: do g_set_index = 1,number_of_energy_group_set

         if (g_set_index == g_set) exit energy_group_set_loop
            
         ! set the energy_group_set path
         energy_group_set_path = trim(particle_option_path)//'/energy_discretisation/energy_group_set['//int2str(g_set_index - 1)//']'
                                             
         ! get the number_energy_groups within this set
         call get_option(trim(energy_group_set_path)//'/number_of_energy_groups',number_of_energy_groups_g_set)         
                  
         first_g_in_g_set = first_g_in_g_set + number_of_energy_groups_g_set
         
      end do energy_group_set_loop
      
   end subroutine first_g_within_group_set
   
   ! --------------------------------------------------------------------------

   subroutine find_total_number_energy_groups(particle_option_path, &
                                              total_number_energy_groups)
      
      !!< Determine the total number of energy groups for this particle type
      !!< via counting how many in each group set 
      
      character(len=*), intent(in) :: particle_option_path
      integer, intent(out) :: total_number_energy_groups
      
      ! local variables
      integer :: g_set
      integer :: number_of_energy_group_set
      integer :: number_of_energy_groups_g_set
      character(len=OPTION_PATH_LEN) :: energy_group_set_path
      
      ! deduce the number of energy group sets
      number_of_energy_group_set = option_count(trim(particle_option_path)//'/energy_discretisation/energy_group_set')
      
      total_number_energy_groups = 0
      
      energy_group_set_loop: do g_set = 1,number_of_energy_group_set
            
         ! set the energy_group_set path
         energy_group_set_path = trim(particle_option_path)//'/energy_discretisation/energy_group_set['//int2str(g_set - 1)//']'
            
         ! get the number_energy_groups within this set
         call get_option(trim(energy_group_set_path)//'/number_of_energy_groups',number_of_energy_groups_g_set)         
         
         total_number_energy_groups = total_number_energy_groups +  &
                                      number_of_energy_groups_g_set
         
      end do energy_group_set_loop      

   end subroutine find_total_number_energy_groups
   
   ! --------------------------------------------------------------------------

end module radiation_energy_group_set_tools
