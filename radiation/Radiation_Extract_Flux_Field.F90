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

module radiation_extract_flux_field

   !!< This module contains procedures associated with extracting the particle flux fields from state
   
   !!< It is currently hard wired for only one angular moment per group (ie diffusion theory) :(
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use fields
   
   use radiation_energy_group_set_tools
   
   implicit none
   
   private 

   public :: extract_flux_all_group, &
             deallocate_flux_all_group, &
             extract_flux_group_g, &
             extract_flux_group_from_group_set 

contains

   ! --------------------------------------------------------------------------

   subroutine extract_flux_all_group(state, &
                                     particle_name, &
                                     number_of_energy_groups, &
                                     particle_flux, &
                                     particle_flux_old)

      !!< Extract the pointer array to the particle flux fields for a particular type

      type(state_type), intent(in) :: state      
      character(len=*), intent(in) :: particle_name   
      integer, intent(in) :: number_of_energy_groups         
      type(scalar_field_pointer), dimension(:), pointer, optional :: particle_flux 
      type(scalar_field_pointer), dimension(:), pointer, optional :: particle_flux_old
      
      ! local variables
      integer :: g
      
      check_present: if (present(particle_flux)) then
      
         check_associated: if (.not. associated(particle_flux)) then
      
            allocate(particle_flux(number_of_energy_groups))
         
            group_loop: do g = 1,number_of_energy_groups
         
               allocate(particle_flux(g)%ptr)
         
            end do group_loop
                  
         else check_associated
         
            FLAbort('Cannot extract flux all group as already associated')
      
         end if check_associated
      
      end if check_present
      
      check_present_old: if (present(particle_flux_old)) then
      
         check_associated_old: if (.not. associated(particle_flux_old)) then
      
            allocate(particle_flux_old(number_of_energy_groups))
         
            group_loop_old: do g = 1,number_of_energy_groups
         
               allocate(particle_flux_old(g)%ptr)
         
            end do group_loop_old
                  
         else check_associated_old
         
            FLAbort('Cannot extract flux old all group as already associated')
      
         end if check_associated_old
      
      end if check_present_old
      
      group_loop_extract: do g = 1,number_of_energy_groups
         
         call extract_flux_group_g(state, &
                                   trim(particle_name), & 
                                   g, &
                                   particle_flux = particle_flux(g)%ptr, &
                                   particle_flux_old = particle_flux_old(g)%ptr)
         
      end do group_loop_extract
        
   end subroutine extract_flux_all_group
   
   ! --------------------------------------------------------------------------

   subroutine deallocate_flux_all_group(particle_flux, &
                                        particle_flux_old)
      
      !!< Deallocate the pointer array to the particle flux fields
      
      type(scalar_field_pointer), dimension(:), pointer, optional :: particle_flux 
      type(scalar_field_pointer), dimension(:), pointer, optional :: particle_flux_old
      
      check_present: if (present(particle_flux)) then
      
         if (associated(particle_flux)) deallocate(particle_flux)
            
      end if check_present

      check_present_old: if (present(particle_flux_old)) then
      
         if (associated(particle_flux_old)) deallocate(particle_flux_old)
            
      end if check_present_old

   end subroutine deallocate_flux_all_group

   ! --------------------------------------------------------------------------

   subroutine extract_flux_group_g(state, &
                                   particle_name, & 
                                   g, &
                                   particle_flux, &
                                   particle_flux_old) 
   
      !!< Extract the particle latest and old flux fields for group g from state

      type(state_type), intent(in) :: state
      character(len=*), intent(in) :: particle_name
      integer, intent(in) :: g
      type(scalar_field), pointer, optional :: particle_flux 
      type(scalar_field), pointer, optional :: particle_flux_old
      
      ! local variables
      integer :: status
      character(len=OPTION_PATH_LEN) :: field_name
      character(len=OPTION_PATH_LEN) :: field_name_old      

      extract_latest: if (present(particle_flux)) then
            
         ! form the field name using the particle_name and group number g for moment 1
         field_name = 'ParticleFluxGroup'//int2str(g)//'Moment1'//trim(particle_name)

         ! extract the field
         particle_flux => extract_scalar_field(state, &
                                         trim(field_name), &
                                         stat=status)
         
         if (status /= 0) FLAbort('Failed to extract particle flux field from state')
         
      end if extract_latest
       
      extract_old: if (present(particle_flux_old)) then
       
         ! form the old field name using the particle_name and group g number for moment 1
         field_name_old = 'OldParticleFluxGroup'//int2str(g)//'Moment1'//trim(particle_name)

         ! extract the field
         particle_flux_old => extract_scalar_field(state, &
                                             trim(field_name_old), &
                                             stat=status)

         if (status /= 0) FLAbort('Failed to extract old particle flux field from state')      
      
      end if extract_old
            
   end subroutine extract_flux_group_g

   ! --------------------------------------------------------------------------

   subroutine extract_flux_group_from_group_set(state, &
                                                particle_option_path, & 
                                                g, &
                                                g_set, &
                                                particle_flux, &
                                                particle_flux_old) 
   
      !!< Extract the particle latest and old flux fields for group g
      !!< within group set g_set from state. g is the within set number not the global number

      type(state_type), intent(in) :: state
      character(len=*), intent(in) :: particle_option_path
      integer, intent(in) :: g
      integer, intent(in) :: g_set
      type(scalar_field), pointer, optional :: particle_flux 
      type(scalar_field), pointer, optional :: particle_flux_old
      
      ! local variables
      integer :: status
      integer :: first_g_in_g_set
      integer :: global_g
      character(len=OPTION_PATH_LEN) :: field_name
      character(len=OPTION_PATH_LEN) :: field_name_old      
      character(len=OPTION_PATH_LEN) :: particle_name      

      ! determine the first energy group in this group set
      call first_g_within_group_set(g_set, &
                                    trim(particle_option_path), &
                                    first_g_in_g_set)
       
      ! determine the global group from the first group number in g_set 
      global_g = first_g_in_g_set + g - 1
   
      ! determine the particle name
      call get_option(trim(particle_option_path)//'/name',particle_name)
      
      extract_latest: if (present(particle_flux)) then
            
         ! form the field name using the particle_name and group number global_g for moment 1
         field_name = 'ParticleFluxGroup'//int2str(global_g)//'Moment1'//trim(particle_name)

         ! extract the field
         particle_flux => extract_scalar_field(state, &
                                               trim(field_name), &
                                               stat=status)

         if (status /= 0) FLAbort('Failed to extract particle flux field from state')      

      end if extract_latest
       
      extract_old: if (present(particle_flux_old)) then
       
         ! form the old field name using the particle_name and group global_g number for moment 1
         field_name_old = 'OldParticleFluxGroup'//int2str(global_g)//'Moment1'//trim(particle_name)

         ! extract the field
         particle_flux_old => extract_scalar_field(state, &
                                                   trim(field_name_old), &
                                                   stat=status)

         if (status /= 0) FLAbort('Failed to extract old particle flux field from state')      
      
      end if extract_old
            
   end subroutine extract_flux_group_from_group_set

   ! --------------------------------------------------------------------------

end module radiation_extract_flux_field
