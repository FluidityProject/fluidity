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

module radiation_solve_energy_group_iteration

   !!< This module contains procedures associated with the energy group iteration solver for radiation time problems
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use fields
   
   use radiation_particle_data_type
   use radiation_materials_interpolation
   use radiation_solve_scatter_iteration
   use radiation_copy_flux_values
   use radiation_check_flux_convergence
   
   implicit none
   
   private 

   public :: time_energy_group_iteration
   
   type energy_group_iteration_options_type
      integer :: max_energy_group_iteration
      real :: flux_tolerance_absolute
      logical :: terminate_if_not_converged_energy_group
      character(len=OPTION_PATH_LEN) :: whole_domain_rebalance_energy_group
   end type energy_group_iteration_options_type
   
contains

   ! --------------------------------------------------------------------------

   subroutine time_energy_group_iteration(particle) 
   
      !!< Solve a radiation time step problem via a simple energy group iteration
      !!< with embedded scatter group iterations and acceleration schemes

      type(particle_type), intent(inout) :: particle
      
      ! local variables
      integer :: iter
      integer :: count_group_solves
      integer :: total_count_group_solves
      integer :: petsc_iterations_taken_scatter_iter_all
      integer :: petsc_iterations_taken_energy_group_iter_all
      logical :: energy_group_iteration_converged
      type(energy_group_iteration_options_type) :: energy_group_iteration_options
      
      ewrite(1,*) 'In time_energy_group_iteration'

      ! initialise the total number of group solves performed counter
      total_count_group_solves = 0
      
      call get_energy_group_iteration_options(particle, &
                                              energy_group_iteration_options)
      
      petsc_iterations_taken_energy_group_iter_all = 0
      
      ! copy the latest solution values to the Old fields
      call copy_to_old_values_particle_flux(particle)
                  
      energy_group_iteration: do iter = 1,energy_group_iteration_options%max_energy_group_iteration
         
         ewrite(1,*) 'Energy Group iteration: ',iter
         
         ! copy the latest solution values to the Iterated fields
         call copy_to_iter_values_particle_flux(particle)
         
         ! initialise the number of group solves performed for this energy_group iteration counter
         count_group_solves = 0

         call scatter_iteration(particle, &
                                count_group_solves, &
                                petsc_iterations_taken_scatter_iter_all, &
                                invoke_eigenvalue_scatter_solve = .false.)

         ! check the flux convergence
         call check_particle_flux_convergence(particle, &
                                              energy_group_iteration_options%flux_tolerance_absolute, &
                                              energy_group_iteration_converged)
         
         ewrite(1,*) 'Number of group solves for this energy_group iteration: ',count_group_solves
         
         ewrite(1,*) 'Total number of petsc iterations taken for this energy_group iteration only: ', &
                     &petsc_iterations_taken_scatter_iter_all
         
         total_count_group_solves = total_count_group_solves + count_group_solves
         
         petsc_iterations_taken_energy_group_iter_all = &
         petsc_iterations_taken_energy_group_iter_all + &
         petsc_iterations_taken_scatter_iter_all

         if (energy_group_iteration_converged) exit energy_group_iteration

      end do energy_group_iteration
      
      ewrite(1,*) 'Total number of group solves for energy_group iteration method: ', &
                  &total_count_group_solves
      
      ewrite(1,*) 'Total number of petsc iterations taken for energy_group iteration method: ',&
                  &petsc_iterations_taken_energy_group_iter_all

      need_to_terminate: if ((iter == energy_group_iteration_options%max_energy_group_iteration) .and.  &
                             energy_group_iteration_options%terminate_if_not_converged_energy_group .and. &
                             (.not. energy_group_iteration_converged)) then
         
         FLExit('Terminating as energy_group iteration finished and solution NOT converged')
         
      end if need_to_terminate
            
   end subroutine time_energy_group_iteration

   ! --------------------------------------------------------------------------
   
   subroutine get_energy_group_iteration_options(particle, &
                                                 energy_group_iteration_options)
      
      !!< Get the energy_group iteration options
      
      type(particle_type), intent(in) :: particle
      type(energy_group_iteration_options_type), intent(out) :: energy_group_iteration_options
           
      ! local variable
      character(len=OPTION_PATH_LEN) :: energy_group_iteration_option_path
                        
      ! get the maximum number of energy_group iterations for this np 
      energy_group_iteration_option_path = trim(particle%option_path)//'/equation/energy_group_iteration'
      
      call get_option(trim(energy_group_iteration_option_path)//'/maximum', &
                      energy_group_iteration_options%max_energy_group_iteration)

      ! get the flux energy_group iteration absolute tolerance
      call get_option(trim(energy_group_iteration_option_path)//'/flux_tolerance_absolute', &
                      energy_group_iteration_options%flux_tolerance_absolute)

      energy_group_iteration_options%terminate_if_not_converged_energy_group = &
      have_option(trim(energy_group_iteration_option_path)//'/terminate_if_not_converged')
            
   end subroutine get_energy_group_iteration_options

   ! --------------------------------------------------------------------------

end module radiation_solve_energy_group_iteration
