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

module radiation_copy_flux_values

   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use fields

   use radiation_particle_data_type
   use radiation_extract_flux_field
   use radiation_energy_group_set_tools

   implicit none
   
   private 

   public :: copy_to_old_values_particle_flux

contains

   ! --------------------------------------------------------------------------

   subroutine copy_to_old_values_particle_flux(particle) 
      
      !!< Copy the latest flux to the old values 

      type(particle_type), intent(inout) :: particle
      
      ! local variables
      integer :: g
      integer :: number_of_energy_groups
      type(scalar_field), pointer :: particle_flux 
      type(scalar_field), pointer :: particle_flux_old
      
      ewrite(1,*) 'Copy latest radiation particle flux values to old for type ',trim(particle%name)
      
      ! find the number of energy groups
      call find_total_number_energy_groups(trim(particle%option_path), &
                                           number_of_energy_groups)
                  
      ! set the old flux
      group_loop: do g = 1,number_of_energy_groups
         
         call extract_flux_group_g(particle, &
                                   g, &  
                                   particle_flux = particle_flux, &
                                   particle_flux_old = particle_flux_old)
            
         call set(particle_flux_old, particle_flux)
                  
      end do group_loop
    
   end subroutine copy_to_old_values_particle_flux 

   ! --------------------------------------------------------------------------

end module radiation_copy_flux_values
