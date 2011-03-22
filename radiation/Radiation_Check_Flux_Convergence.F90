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
   use parallel_tools

   use radiation_extract_flux_field

   implicit none
   
   private 

   public ::check_particle_flux_convergence 

contains

   ! --------------------------------------------------------------------------

   subroutine check_particle_flux_convergence(state, &
                                              particle_name, &
                                              number_of_energy_groups, &
                                              flux_tolerance, &
                                              max_change_flux, &
                                              flux_converged)  
      
      !!< Check the particle flux relative convergence for all energy groups 

      type(state_type), intent(in) :: state
      character(len=*), intent(in) :: particle_name
      integer, intent(in) :: number_of_energy_groups  
      real, intent(in) :: flux_tolerance
      real, intent(out) :: max_change_flux
      logical, intent(out) :: flux_converged          
      
      ! local variable
      integer :: g
      integer :: n
      logical :: found_non_zero_value
      real :: difference
      type(scalar_field), pointer :: particle_flux 
      type(scalar_field), pointer :: particle_flux_old
      
      max_change_flux      = 0.0
      flux_converged       = .false.
      found_non_zero_value = .false.
      
      group_loop: do g = 1,number_of_energy_groups
         
         call extract_flux_group_g(state, &
                                   trim(particle_name), & 
                                   g, &  
                                   particle_flux = particle_flux, &
                                   particle_flux_old = particle_flux_old)
         
         node_loop: do n = 1,node_count(particle_flux)
            
            ! if the flux old is zero (say for a zero BC) then cycle as no need to find relative error
            if (abs(node_val(particle_flux_old,n)) <= epsilon(0.0)) cycle node_loop
            
            ! on first find of non zero value set a flag
            if (.not. found_non_zero_value) found_non_zero_value = .true.
            
            ! the absolute difference
            difference = abs(node_val(particle_flux,n) - node_val(particle_flux_old,n))
            
            ! the maximum change over all groups and space                  
            max_change_flux = max( max_change_flux, difference/abs(node_val(particle_flux_old,n)) )
                               
            check_flux: if (difference < abs(node_val(particle_flux_old,n)*flux_tolerance)) then
               
               flux_converged = .true.
               
            else check_flux
               
               flux_converged = .false.
                  
               exit node_loop
                                    
            end if check_flux
               
          end do node_loop
          
          ! determine if all processes are converged
          call alland(flux_converged)
          
          ! if not converged at any found node exit - no need to check all off them
          if (.not. flux_converged) exit group_loop
                  
      end do group_loop
      
      no_non_zero: if (.not. found_non_zero_value) then
         
         flux_converged = .true.
         
         ewrite(1,*) 'WARNING no non zero old particle flux values found in check converegence'
      
      end if no_non_zero
      
      ewrite(1,*) 'max_change_flux,flux_converged: ',max_change_flux,flux_converged  
      
   end subroutine check_particle_flux_convergence  

   ! --------------------------------------------------------------------------

end module radiation_check_flux_convergence
