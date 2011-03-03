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

module radiation_check_flux_convergence

   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use fields

   use radiation_extract_flux_field

   implicit none
   
   private 

   public ::check_np_flux_convergence 

contains

   ! --------------------------------------------------------------------------

   subroutine check_np_flux_convergence(state, &
                                        np_radmat_name, &
                                        number_of_energy_groups, &
                                        flux_tolerance, &
                                        max_change_flux, &
                                        flux_converged)  
      
      !!< Check the neutral particle flux convergence

      type(state_type), intent(in) :: state
      character(len=*), intent(in) :: np_radmat_name
      integer, intent(in) :: number_of_energy_groups  
      real, intent(in) :: flux_tolerance
      real, intent(out) :: max_change_flux
      logical, intent(out) :: flux_converged          
      
      ! local variable
      integer :: g
      integer :: n
      real :: difference
      type(scalar_field), pointer :: np_flux 
      type(scalar_field), pointer :: np_flux_old
      
      max_change_flux = 0.0
      flux_converged = .false.
      group_loop: do g = 1,number_of_energy_groups
         
         call extract_flux_group_g(state, &
                                   trim(np_radmat_name), & 
                                   g, &  
                                   np_flux = np_flux, &
                                   np_flux_old = np_flux_old)
         
         node_loop: do n = 1,node_count(np_flux)

            difference = abs(node_val(np_flux,n) - node_val(np_flux_old,n))
               
            find_max_change_flux: if ( node_val(np_flux_old,n) > 0) then
               
               max_change_flux = max( max_change_flux, difference/abs(node_val(np_flux_old,n)) )
                
            end if find_max_change_flux
               
            check_flux: if (difference < abs(node_val(np_flux_old,n)*flux_tolerance)) then
               
               flux_converged = .true.
               
            else check_flux
               
               flux_converged = .false.
                  
               ! exit the entire flux loop over groups and nodes if found a failed flux convereged
               exit group_loop
                                    
            end if check_flux
               
          end do node_loop
                  
      end do group_loop
      
      ewrite(1,*) 'max_change_flux,flux_converged: ',max_change_flux,flux_converged  
      
   end subroutine check_np_flux_convergence  

   ! --------------------------------------------------------------------------

end module radiation_check_flux_convergence
