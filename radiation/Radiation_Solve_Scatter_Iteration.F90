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

module radiation_solve_scatter_iteration

   !!< This module contains procedures associated with scatter iteration solver for radiation problems
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  

   use radiation_particle_data_type
   use radiation_assemble_solve_group
   use radiation_energy_group_set_tools

   implicit none
   
   private 

   public :: scatter_iteration
   
   type scatter_iteration_options_type
      integer :: max_scatter_iteration
      integer :: highest_upscatter_group
      real :: flux_tolerance_absolute
      logical :: terminate_if_not_converged_scatter
      character(len=OPTION_PATH_LEN) :: whole_domain_rebalance_scatter
      character(len=OPTION_PATH_LEN) :: energy_solve_direction         
   end type scatter_iteration_options_type
   
contains

   ! --------------------------------------------------------------------------

   subroutine scatter_iteration(particle, &
                                count_group_solves, &
                                petsc_iterations_taken_scatter_iter_all, &
                                invoke_eigenvalue_scatter_solve)
      
      !!< Perform iterations around the group loop to solve scatter
      !!< which could involve a dual sweep in both directions being
      !!< high energy to low and low energy to high. A special case 
      !!< of one energy group exists where this is not applicable.
      
      type(particle_type), intent(inout) :: particle
      integer, intent(inout) :: count_group_solves
      logical, intent(in) :: invoke_eigenvalue_scatter_solve
      integer, intent(out) :: petsc_iterations_taken_scatter_iter_all
      
      ! local variables
      integer :: sweeps
      integer :: isweep
      integer :: iscatter
      integer :: start_group
      integer :: end_group
      integer :: group_increment
      integer :: number_of_energy_groups
      integer :: petsc_iterations_taken_scatter_iter
      logical :: scatter_iteration_converged
      type(scatter_iteration_options_type) :: scatter_iteration_options
      character(len=OPTION_PATH_LEN) :: scatter_group_iteration_option_path
       
      ! set the scatter iteration option path
      scatter_path_if: if (invoke_eigenvalue_scatter_solve) then
      
         scatter_group_iteration_option_path = &
         trim(particle%option_path)//'/equation/power_iteration/scatter_group_iteration'
      
      else scatter_path_if
      
         scatter_group_iteration_option_path = &
         trim(particle%option_path)//'/equation/energy_group_iteration/scatter_group_iteration'
      
      end if scatter_path_if
      
      call get_scatter_iteration_options(scatter_iteration_options, &
                                         trim(scatter_group_iteration_option_path))
      
      ! find the number of energy groups
      call find_total_number_energy_groups(trim(particle%option_path), &
                                           number_of_energy_groups)
           
      ! first iteration sweep down groups from 1 through all groups with increments of 1
      start_group     = 1
      end_group       = number_of_energy_groups
      group_increment = 1
      sweeps          = 1
      
      ! if there is only one energy group there is only a need for one scatter iteration
      one_group: if (number_of_energy_groups == 1) then
         
         scatter_iteration_options%max_scatter_iteration = 1
      
      end if one_group
      
      petsc_iterations_taken_scatter_iter_all = 0
      
      scatter_iteration_loop: do iscatter = 1,scatter_iteration_options%max_scatter_iteration
         
         ewrite(1,*) 'Scatter iteration: ',iscatter
         
         petsc_iterations_taken_scatter_iter = 0
         
         choose_groups: if (iscatter > 1) then
         
            call choose_group_first_sweep_bounds(iscatter, &
                                                 number_of_energy_groups, &
                                                 scatter_iteration_options%highest_upscatter_group, &
                                                 sweeps, &
                                                 group_increment, &
                                                 start_group, &
                                                 end_group, &
                                                 trim(scatter_iteration_options%energy_solve_direction))
         
         end if choose_groups
         
         group_sweeps_loop: do isweep = 1,sweeps
            
            to_many_sweeps: if (isweep > 2) then
               
               FLAbort('Greater than 2 energy direction sweeps')
               
            end if to_many_sweeps
            
            second_sweep: if (isweep == 2) then
            
               call choose_group_second_sweep_bounds(number_of_energy_groups, &
                                                     scatter_iteration_options%highest_upscatter_group, &
                                                     group_increment, &
                                                     start_group, &
                                                     end_group, &
                                                     trim(scatter_iteration_options%energy_solve_direction))
            
            end if second_sweep
            
            ! loop the necessary energy groups
            call energy_group_loop(particle, &
                                   start_group, &
                                   end_group, &
                                   group_increment, &
                                   count_group_solves, &
                                   invoke_eigenvalue_scatter_solve, &
                                   petsc_iterations_taken_scatter_iter)
         
         end do group_sweeps_loop
                  
         ! rebalance accelerate the scatter iteration if options chosen
         rebalance: if (trim(scatter_iteration_options%whole_domain_rebalance_scatter) /= 'none') then
         
            call scatter_iteration_whole_domain_rebalance()
         
         end if rebalance
         
         call check_scatter_iteration_convergence(scatter_iteration_converged)
         
         petsc_iterations_taken_scatter_iter_all = &
         petsc_iterations_taken_scatter_iter_all + &
         petsc_iterations_taken_scatter_iter
         
         if (scatter_iteration_converged) exit scatter_iteration_loop
         
      end do scatter_iteration_loop

      need_to_terminate: if ((iscatter == scatter_iteration_options%max_scatter_iteration) .and. &
                              scatter_iteration_options%terminate_if_not_converged_scatter .and. &
                              (.not. scatter_iteration_converged)) then
         
         FLExit('Terminating as scatter iteration finished and solution NOT converged')
         
      end if need_to_terminate
            
   end subroutine scatter_iteration 

   ! --------------------------------------------------------------------------

   subroutine get_scatter_iteration_options(scatter_iteration_options, &
                                            scatter_group_iteration_option_path)
      
      !!< Get the scatter iteration options
      
      type(scatter_iteration_options_type), intent(out) :: scatter_iteration_options
      character(len=*), intent(in) :: scatter_group_iteration_option_path
                
      call get_option(trim(scatter_group_iteration_option_path)//'/maximum', &
                      scatter_iteration_options%max_scatter_iteration, &
                      default = 1)

      call get_option(trim(scatter_group_iteration_option_path)//'/flux_tolerance_absolute', &
                      scatter_iteration_options%flux_tolerance_absolute, &
                      default = 1.0)
                      
      call get_option(trim(scatter_group_iteration_option_path)//'/energy_solve_direction/name', &
                      scatter_iteration_options%energy_solve_direction, &
                      default = 'HighToLow')
         
      scatter_iteration_options%terminate_if_not_converged_scatter = &
      have_option(trim(scatter_group_iteration_option_path)//'/terminate_if_not_converged')
         
      call get_option(trim(scatter_group_iteration_option_path)//'/highest_upscatter_group', &
                      scatter_iteration_options%highest_upscatter_group, &
                      default = 1)
                      
      call get_option(trim(scatter_group_iteration_option_path)//'/whole_domain_group_rebalance/name', &
                      scatter_iteration_options%whole_domain_rebalance_scatter, &
                      default = 'none')
       
   end subroutine get_scatter_iteration_options 

   ! --------------------------------------------------------------------------
   
   subroutine choose_group_first_sweep_bounds(iscatter, &
                                              number_of_energy_groups, &
                                              highest_upscatter_group, &
                                              sweeps, &
                                              group_increment, &
                                              start_group, &
                                              end_group, &
                                              energy_solve_direction)
      
      !!< Choose the group first sweep bounds options for this scatter iteration iscatter
      
      integer, intent(in) :: iscatter
      integer, intent(in) :: number_of_energy_groups
      integer, intent(in) :: highest_upscatter_group
      integer, intent(out) :: sweeps
      integer, intent(out) :: group_increment
      integer, intent(out) :: start_group
      integer, intent(out) :: end_group
      character(len=*), intent(in) :: energy_solve_direction
            
      group_direction: if (trim(energy_solve_direction) == 'HighToLow') then
               
         sweeps = 1
               
         start_group = highest_upscatter_group
         
         end_group = number_of_energy_groups
         
         group_increment = +1      
      
      else if (trim(energy_solve_direction) == 'LowToHigh') then group_direction
               
         sweeps = 1
               
         if (iscatter == 2) then
               
            start_group = number_of_energy_groups - 1
               
         else 
                  
            start_group = number_of_energy_groups
                  
         end if
               
         end_group = highest_upscatter_group
               
         group_increment = -1
              
      else if(trim(energy_solve_direction) == 'HighToLowToHigh') then group_direction
               
         sweeps = 2
                              
         if (iscatter == 2) then
               
            start_group = highest_upscatter_group
               
         else 
               
            start_group = highest_upscatter_group + 1
               
         end if 
               
         end_group = number_of_energy_groups
               
         group_increment = +1
               
      else if(trim(energy_solve_direction) == 'LowToHighToLow') then group_direction
               
         sweeps = 2
               
         start_group = number_of_energy_groups - 1
               
         end_group = highest_upscatter_group
               
         group_increment = -1               
            
      else group_direction
               
         FLAbort('Unknown radiation energy dof solve direction')
               
      end if group_direction
      
   end subroutine choose_group_first_sweep_bounds 
   
   ! --------------------------------------------------------------------------
   
   subroutine choose_group_second_sweep_bounds(number_of_energy_groups, &
                                               highest_upscatter_group, &
                                               group_increment, &
                                               start_group, &
                                               end_group, &
                                               energy_solve_direction)
      
      !!< Choose the group second sweep bounds options for this scatter iteration iscatter
      
      integer, intent(in) :: number_of_energy_groups
      integer, intent(in) :: highest_upscatter_group
      integer, intent(out) :: group_increment
      integer, intent(out) :: start_group
      integer, intent(out) :: end_group
      character(len=*), intent(in) :: energy_solve_direction
               
      dual_group_direction: if (trim(energy_solve_direction) == 'HighToLowToHigh') then
               
         start_group = number_of_energy_groups - 1
               
         end_group = highest_upscatter_group
               
         group_increment = -1
               
      else dual_group_direction
               
         start_group = highest_upscatter_group + 1
               
         end_group = number_of_energy_groups
               
         group_increment = +1
            
      end if dual_group_direction               
      
   end subroutine choose_group_second_sweep_bounds 
   
   ! --------------------------------------------------------------------------
   
   subroutine energy_group_loop(particle, &
                                start_group, &
                                end_group, &
                                group_increment, &
                                count_group_solves, &
                                invoke_eigenvalue_group_solve, &
                                petsc_iterations_taken_scatter_iter)
      
      !!< Sweep the energy groups from the start_group to the end_group
      !!< in steps of group_increment (which could be negative)
      !!< solving the within group particle balance
      
      type(particle_type), intent(inout) :: particle
      integer, intent(in) :: start_group
      integer, intent(in) :: end_group
      integer, intent(in) :: group_increment
      integer, intent(inout) :: count_group_solves
      logical, intent(in) :: invoke_eigenvalue_group_solve
      integer, intent(inout) :: petsc_iterations_taken_scatter_iter
      
      ! local variables
      integer :: g
      integer :: petsc_iterations_taken_group_g
      
      ! first sweep of energy dof from start to end
      group_loop: do g = start_group,end_group,group_increment
         
         ewrite(1,*) 'Assemble and solve energy group ',g
         
         count_group_solves = count_group_solves + 1
         
         ! Assemble and solve the group g particle balance
         call particle_assemble_solve_group(particle, &
                                            g, &
                                            invoke_eigenvalue_group_solve, &
                                            petsc_iterations_taken_group_g)
         
         petsc_iterations_taken_scatter_iter = &
         petsc_iterations_taken_scatter_iter + &
         petsc_iterations_taken_group_g
         
      end do group_loop
      
   end subroutine energy_group_loop

   ! --------------------------------------------------------------------------
   
   subroutine scatter_iteration_whole_domain_rebalance()  
      
      !!< Rebalance accelerate the scatter iteration over the whole domain
      
      
   end subroutine scatter_iteration_whole_domain_rebalance

   ! --------------------------------------------------------------------------
   
   subroutine check_scatter_iteration_convergence(scatter_iteration_converged)  
      
      !!< Check the convergence of the scatter iterations
      
      logical, intent(out) :: scatter_iteration_converged
      
      scatter_iteration_converged = .false.
                  
   end subroutine check_scatter_iteration_convergence

   ! --------------------------------------------------------------------------

end module radiation_solve_scatter_iteration
