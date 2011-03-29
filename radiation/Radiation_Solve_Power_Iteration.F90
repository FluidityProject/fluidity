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

module radiation_solve_power_iteration

   !!< This module contains procedures associated with power iteration solver for radiation eigenvalue problems
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use fields
   
   use diagnostic_integrate_fields

   use radiation_particle_data_type
   use radiation_materials_interpolation
   use radiation_solve_scatter_iteration
   use radiation_extract_flux_field
   use radiation_copy_flux_values
   use radiation_check_flux_convergence
   use radiation_energy_group_set_tools
   
   implicit none
   
   private 

   public :: eigenvalue_power_iteration
   
   type power_iteration_options_type
      integer :: max_power_iteration
      real :: keff_tolerance_relative
      real :: flux_tolerance_absolute
      logical :: terminate_if_not_converged_power            
   end type power_iteration_options_type
   
contains

   ! --------------------------------------------------------------------------

   subroutine eigenvalue_power_iteration(particle) 
   
      !!< Solve a radiation eigenvalue problem via a power iteration

      type(particle_type), intent(inout) :: particle
      
      ! local variables
      integer :: ipower
      integer :: number_all_convereged_pass
      integer :: count_group_solves
      integer :: total_count_group_solves
      logical :: power_iteration_converged
      type(power_iteration_options_type) :: power_iteration_options
      type(scalar_field), pointer :: keff_field
      
      ewrite(1,*) 'In eigenvalue_power_iteration'
      
      ! keff type already intialised to 1.0
      
      ! intialise the number of all convergence passes in power iterations
      number_all_convereged_pass = 0
      
      ! initialise the total number of groups solves performed counter
      total_count_group_solves = 0
      
      call get_power_iteration_options(particle, &
                                       power_iteration_options)                                       
                  
      power_iteration: do ipower = 1,power_iteration_options%max_power_iteration
         
         ewrite(1,*) 'Power iteration: ',ipower
         
         ! initialise the number of group solves performed for this power iteration counter
         count_group_solves = 0
         
         call copy_to_old_values_eig(particle)
                  
         call scatter_iteration(particle, &
                                count_group_solves, &
                                invoke_eigenvalue_scatter_solve = .true.)
                           
         call calculate_eigenvalue(particle)

         call check_power_iteration_convergence(power_iteration_converged, &
                                                number_all_convereged_pass, &
                                                particle, &
                                                power_iteration_options)
         
         ewrite(1,*) 'Number of group solves for this power iteration: ',count_group_solves
         
         total_count_group_solves = total_count_group_solves + count_group_solves
                           
         if (power_iteration_converged) exit power_iteration
                  
      end do power_iteration
      
      ewrite(1,*) 'Total number of group solves for power iteration method: ',total_count_group_solves
      
      ! set the Keff into the constant scalar field for output
      keff_field => extract_scalar_field(particle%state, &
                                         'ParticleKeff'//trim(particle%name))    
       
      call set(keff_field, &
               particle%keff%keff_new) 
       
      need_to_terminate: if ((ipower == power_iteration_options%max_power_iteration) .and.  &
                             power_iteration_options%terminate_if_not_converged_power .and. &
                             (.not. power_iteration_converged)) then
         
         FLExit('Terminating as power iteration finished and solution NOT converged')
         
      end if need_to_terminate
            
   end subroutine eigenvalue_power_iteration

   ! --------------------------------------------------------------------------
   
   subroutine get_power_iteration_options(particle, &
                                          power_iteration_options)
      
      !!< Get the power iteration options
      
      type(particle_type), intent(in) :: particle
      type(power_iteration_options_type), intent(out) :: power_iteration_options
           
      ! local variable
      character(len=OPTION_PATH_LEN) :: power_iteration_option_path
                        
      ! get the maximum number of power iterations for this np 
      power_iteration_option_path = trim(particle%option_path)//'/equation/power_iteration'
      
      call get_option(trim(power_iteration_option_path)//'/maximum', &
                      power_iteration_options%max_power_iteration)

      ! get the keff power iteration relative tolerance
      call get_option(trim(power_iteration_option_path)//'/keff_tolerance_relative', &
                      power_iteration_options%keff_tolerance_relative)

      ! get the flux power iteration absolute tolerance
      call get_option(trim(power_iteration_option_path)//'/flux_tolerance_absolute', &
                      power_iteration_options%flux_tolerance_absolute)

      power_iteration_options%terminate_if_not_converged_power = &
      have_option(trim(power_iteration_option_path)//'/terminate_if_not_converged')         
            
   end subroutine get_power_iteration_options

   ! --------------------------------------------------------------------------

   subroutine copy_to_old_values_eig(particle) 
      
      !!< Copy the latest flux and eigenvalue to the old values 
      
      type(particle_type), intent(inout) :: particle      
                  
      ! set the eigenvalue old
      particle%keff%keff_old = particle%keff%keff_new
      
      call copy_to_old_values_particle_flux(particle)     
            
   end subroutine copy_to_old_values_eig 

   ! --------------------------------------------------------------------------

   subroutine calculate_eigenvalue(particle) 
      
      !!< Calculate the latest eigenvalue as:
      !!< keff_new = keff_old*(integral_dv(production_source*production_source)/integral_dv(production_source*production_source_old))
      
      type(particle_type), intent(inout) :: particle
      
      ! local variables
      integer :: status
      integer :: g
      integer :: g_set
      integer :: g_global
      integer :: number_of_energy_groups
      integer :: number_of_energy_group_set
      real :: k_top
      real :: k_bottom
      real :: k_top_group
      real :: k_bottom_group
      type(vector_field), pointer :: positions 
      type(scalar_field), pointer :: particle_flux 
      type(scalar_field), pointer :: particle_flux_old
      type(scalar_field), target :: production_coeff      
      type(mesh_type), pointer :: material_fn_space
      type(scalar_field_pointer), dimension(:), pointer :: scalar_fields 
      character(len=OPTION_PATH_LEN) :: positions_mesh_name
      character(len=OPTION_PATH_LEN) :: material_fn_space_name
      character(len=OPTION_PATH_LEN) :: energy_group_set_path
      
      ewrite(1,*) 'Calculate Keff'
      
      ! intialise variables that are summed up
      k_top    = 0.0
      k_bottom = 0.0
      
      ! deduce the number of energy group sets
      number_of_energy_group_set = option_count(trim(particle%option_path)//'/energy_discretisation/energy_group_set')
      
      ! initialise the global group counter
      g_global = 0

      ! allocate the fields to integrate for k_top and k_bottom
      allocate(scalar_fields(4))

      allocate(scalar_fields(1)%ptr)
      allocate(scalar_fields(2)%ptr)
      allocate(scalar_fields(3)%ptr)
      allocate(scalar_fields(4)%ptr)
      
      energy_group_set_loop: do g_set = 1,number_of_energy_group_set

         ! set the positions mesh name for this group set - currently assume all group sets have same positions
         positions_mesh_name = 'Coordinate'
         
         ! extract the positions 
         positions => extract_vector_field(particle%state, &
                                           trim(positions_mesh_name), &
                                           stat=status)  

         ! set the energy_group_set path
         energy_group_set_path = trim(particle%option_path)//'/energy_discretisation/energy_group_set['//int2str(g_set - 1)//']'
         
         ! get the material fn space name for this group set
         call get_option(trim(energy_group_set_path)//'/angular_discretisation/method/parity/angular_moment_set[0]/mesh/name',material_fn_space_name)
      
         ! extract the material fn_space of this energy group set of this particle type 
         material_fn_space => extract_mesh(particle%state, &
                                           trim(material_fn_space_name))

         ! get the number_energy_groups within this set
         call get_option(trim(energy_group_set_path)//'/number_of_energy_groups',number_of_energy_groups)         
         
         ! allocate the the production field for this group set
         call allocate(production_coeff, &
                       material_fn_space, & 
                       'ParticleProduction')
                           
         ! integrate each energy group within this group set
         group_loop: do g = 1,number_of_energy_groups
            
            g_global = g_global + 1
            
            ! extract the group particle flux
            call extract_flux_group_g(particle, &
                                      g_global, &
                                      particle_flux = particle_flux, &
                                      particle_flux_old = particle_flux_old)
            
            call zero(production_coeff)
            
            ! form the production coeff field for this energy group
            call form(material_fn_space, &
                      particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                      particle%particle_radmat, &
                      production_coeff, &
                      g, &
                      component = 'production')
            
            ! set the fields to integrate for k_top
            scalar_fields(1)%ptr => production_coeff
            scalar_fields(2)%ptr => particle_flux
            scalar_fields(3)%ptr => production_coeff
            scalar_fields(4)%ptr => particle_flux

            call integrate(scalar_fields, &
                           positions, &
                           k_top_group)
            
            k_top = k_top + k_top_group

            ! set the fields to integrate for k_bottom
            scalar_fields(1)%ptr => production_coeff
            scalar_fields(2)%ptr => particle_flux
            scalar_fields(3)%ptr => production_coeff
            scalar_fields(4)%ptr => particle_flux_old

            call integrate(scalar_fields, &
                           positions, &
                           k_bottom_group)

            k_bottom = k_bottom + k_bottom_group
            
         end do group_loop
         
         ! deallocate this group set production coeff field
         call deallocate(production_coeff)
      
      end do energy_group_set_loop

      if (associated(scalar_fields)) deallocate(scalar_fields)
      
      ! form the latest keff estimate
      particle%keff%keff_new = particle%keff%keff_old*k_top/k_bottom

      ewrite(1,*) 'Keff: ',particle%keff%keff_new     
                
   end subroutine calculate_eigenvalue  

   ! --------------------------------------------------------------------------

   subroutine check_power_iteration_convergence(power_iteration_converged, &
                                                number_all_convereged_pass, &
                                                particle, &
                                                power_iteration_options) 
   
      !!< Check the convergence of the power iterations via the keff and flux
      !!< Always insist on atleast 2 successive convergence passes   
      
      logical, intent(out) :: power_iteration_converged
      integer, intent(inout) :: number_all_convereged_pass 
      type(particle_type), intent(in) :: particle
      type(power_iteration_options_type), intent(in) :: power_iteration_options
      
      ! local variable
      real :: rel_difference_keff
      real :: abs_difference_keff
      logical :: keff_converged
      logical :: flux_converged
            
      power_iteration_converged = .true.
      
      ! find the keff absolute difference
      abs_difference_keff = abs(particle%keff%keff_new - particle%keff%keff_old)
      
      ! find the keff relative difference
      find_change_keff: if (particle%keff%keff_old > 0.0) then
         
         rel_difference_keff = abs_difference_keff / abs(particle%keff%keff_old)
         
      else find_change_keff
         
         FLAbort('ERROR the radiation eigenvalue Keff old is not > 0.0')
                  
      end if find_change_keff
      
      ! check the keff convergence via relative change
      check_keff: if (rel_difference_keff < power_iteration_options%keff_tolerance_relative) then
         
         keff_converged = .true.
      
      else check_keff
         
         keff_converged = .false.
      
      end if check_keff 
      
      ewrite(1,*) 'rel_difference_keff,abs_difference_keff,keff_converged: ', &
                   rel_difference_keff,abs_difference_keff,keff_converged     
      
      ! check the flux convergence
      call check_particle_flux_convergence(particle, &
                                           power_iteration_options%flux_tolerance_absolute, &
                                           flux_converged)
            
      ! decide if power iteration converged
      power_iteration_converged = .true.
      
      if (.not. keff_converged) power_iteration_converged = .false.     

      if (.not. flux_converged) power_iteration_converged = .false.     
      
      ! count the number of all converged passes to insist on two consecutive passes
      if (power_iteration_converged) number_all_convereged_pass = number_all_convereged_pass + 1
      
      if (.not. power_iteration_converged) number_all_convereged_pass = 0
      
      if (number_all_convereged_pass < 2) power_iteration_converged = .false.    

      ewrite(1,*) 'power_iteration_converged: ',power_iteration_converged  
          
   end subroutine check_power_iteration_convergence  

   ! --------------------------------------------------------------------------

end module radiation_solve_power_iteration
