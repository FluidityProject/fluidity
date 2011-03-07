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

   !!< This module contains procedures associated with extracting the neutral particle flux fields from state
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use fields
   
   implicit none
   
   private 

   public :: extract_flux_all_group, &
             deallocate_flux_all_group, &
             extract_flux_group_g 

contains

   ! --------------------------------------------------------------------------

   subroutine extract_flux_all_group(state, &
                                     np_radmat_name, &
                                     number_of_energy_groups, &
                                     np_flux, &
                                     np_flux_old)

      !!< Extract the pointer array to the neutral particle flux fields for a particular np object
      !!< If a general version of this routine doesnt already exist then it needs creating
      !!<  - this would take a field name list and return them in sequential order

      type(state_type), intent(in) :: state      
      character(len=*), intent(in) :: np_radmat_name   
      integer, intent(in) :: number_of_energy_groups         
      type(scalar_field_pointer), dimension(:), pointer, optional :: np_flux 
      type(scalar_field_pointer), dimension(:), pointer, optional :: np_flux_old
      
      ! local variables
      integer :: g
      
      check_present: if (present(np_flux)) then
      
         check_associated: if (.not. associated(np_flux)) then
      
            allocate(np_flux(number_of_energy_groups))
         
            group_loop: do g = 1,number_of_energy_groups
         
               allocate(np_flux(g)%ptr)
         
            end do group_loop
                  
         else check_associated
         
            FLAbort('Cannot extract flux all group as already associated')
      
         end if check_associated
      
      end if check_present
      
      check_present_old: if (present(np_flux_old)) then
      
         check_associated_old: if (.not. associated(np_flux_old)) then
      
            allocate(np_flux_old(number_of_energy_groups))
         
            group_loop_old: do g = 1,number_of_energy_groups
         
               allocate(np_flux_old(g)%ptr)
         
            end do group_loop_old
                  
         else check_associated_old
         
            FLAbort('Cannot extract flux old all group as already associated')
      
         end if check_associated_old
      
      end if check_present_old
      
      group_loop_extract: do g = 1,number_of_energy_groups
         
         call extract_flux_group_g(state, &
                                   trim(np_radmat_name), & 
                                   g, &
                                   np_flux = np_flux(g)%ptr, &
                                   np_flux_old = np_flux_old(g)%ptr)
         
      end do group_loop_extract
        
   end subroutine extract_flux_all_group
   
   ! --------------------------------------------------------------------------

   subroutine deallocate_flux_all_group(np_flux, &
                                        np_flux_old)
      
      !!< Deallocate the pointer array to the neutral particle flux fields
      !!< If a general version of this routine doesnt already exist then it needs creating
      
      type(scalar_field_pointer), dimension(:), pointer, optional :: np_flux 
      type(scalar_field_pointer), dimension(:), pointer, optional :: np_flux_old
      
      check_present: if (present(np_flux)) then
      
         if (associated(np_flux)) deallocate(np_flux)
            
      end if check_present

      check_present_old: if (present(np_flux_old)) then
      
         if (associated(np_flux_old)) deallocate(np_flux_old)
            
      end if check_present_old

   end subroutine deallocate_flux_all_group

   ! --------------------------------------------------------------------------

   subroutine extract_flux_group_g(state, &
                                   np_radmat_name, & 
                                   g, &
                                   np_flux, &
                                   np_flux_old) 
   
      !!< Extract the neutral particle latest and old flux fields for group g from state

      type(state_type), intent(in) :: state
      character(len=*), intent(in) :: np_radmat_name
      integer, intent(in) :: g
      type(scalar_field), pointer, optional :: np_flux 
      type(scalar_field), pointer, optional :: np_flux_old
      
      ! local variables
      integer :: status
      character(len=OPTION_PATH_LEN) :: field_name
      character(len=OPTION_PATH_LEN) :: field_name_old      

      extract_latest: if (present(np_flux)) then
            
         ! form the field name using the np_radmat_name and group number g 
         field_name = 'NeutralParticleFluxGroup'//int2str(g)//trim(np_radmat_name)

         ! extract the field
         np_flux => extract_scalar_field(state, &
                                         trim(field_name), &
                                         stat=status)

      end if extract_latest
       
      extract_old: if (present(np_flux_old)) then
       
         ! form the old field name using the np_radmat_name and group g number
         field_name_old = 'OldNeutralParticleFluxGroup'//int2str(g)//trim(np_radmat_name)

         ! extract the field
         np_flux_old => extract_scalar_field(state, &
                                             trim(field_name_old), &
                                             stat=status)
      
      end if extract_old
            
   end subroutine extract_flux_group_g

   ! --------------------------------------------------------------------------

end module radiation_extract_flux_field
