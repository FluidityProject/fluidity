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
      character(len=OPTION_PATH_LEN) :: flux_acceleration
      real :: flux_relaxation_coeff
      character(len=OPTION_PATH_LEN) :: keff_acceleration
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
      integer :: petsc_iterations_taken_scatter_iter_all
      integer :: petsc_iterations_taken_power_iter_all
      logical :: power_iteration_converged
      type(power_iteration_options_type) :: power_iteration_options
      type(scalar_field), pointer :: keff_field
      
      ewrite(1,*) 'In eigenvalue_power_iteration'
      
      ! keff type already intialised to 1.0
      
      ! intialise the number of all convergence passes in power iterations
      number_all_convereged_pass = 0
      
      ! initialise the total number of group solves performed counter
      total_count_group_solves = 0
      
      call get_power_iteration_options(particle, &
                                       power_iteration_options)
      
      petsc_iterations_taken_power_iter_all = 0
                  
      power_iteration: do ipower = 1,power_iteration_options%max_power_iteration
         
         ewrite(1,*) 'Power iteration: ',ipower
         
         ! initialise the number of group solves performed for this power iteration counter
         count_group_solves = 0
         
         call copy_to_iter_values_eig(particle)
                  
         call scatter_iteration(particle, &
                                count_group_solves, &
                                petsc_iterations_taken_scatter_iter_all, &
                                invoke_eigenvalue_scatter_solve = .true.)

         call power_iteration_flux_acceleration(particle, &
                                                power_iteration_options)
         
         call calculate_eigenvalue(particle)
         
         call power_iteration_keff_acceleration(particle, &
                                                power_iteration_options, &
                                                ipower)
         
         call check_power_iteration_convergence(power_iteration_converged, &
                                                number_all_convereged_pass, &
                                                particle, &
                                                power_iteration_options)
         
         ewrite(1,*) 'Number of group solves for this power iteration: ',count_group_solves
         
         ewrite(1,*) 'Total number of petsc iterations taken for this power iteration only: ', &
                     petsc_iterations_taken_scatter_iter_all
         
         total_count_group_solves = total_count_group_solves + count_group_solves
         
         petsc_iterations_taken_power_iter_all = &
         petsc_iterations_taken_power_iter_all + &
         petsc_iterations_taken_scatter_iter_all

         if (power_iteration_converged) exit power_iteration

      end do power_iteration
      
      ewrite(1,*) 'Total number of group solves for power iteration method: ', &
                  total_count_group_solves
      
      ewrite(1,*) 'Total number of petsc iterations taken for power iteration method: ',&
                  &petsc_iterations_taken_power_iter_all
      
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
      
      ! get the flux acceleration method name if any
      call get_option(trim(power_iteration_option_path)//'/flux_acceleration/name', &
                      power_iteration_options%flux_acceleration, &
                      default = 'none')
      
      ! if Relaxation method get the coefficient
      relax_coeff_flux: if (trim(power_iteration_options%flux_acceleration) == 'Relaxation') then
         
         call get_option(trim(power_iteration_option_path)//'/flux_acceleration', &
                         power_iteration_options%flux_relaxation_coeff, &
                         default = 1.0)
      
      else relax_coeff_flux
         
         power_iteration_options%flux_relaxation_coeff = 1.0
         
      end if relax_coeff_flux
      
      ! get the keff acceleration method name if any
      call get_option(trim(power_iteration_option_path)//'/keff_acceleration/name', &
                      power_iteration_options%keff_acceleration, &
                      default = 'none')

   end subroutine get_power_iteration_options

   ! --------------------------------------------------------------------------

   subroutine copy_to_iter_values_eig(particle) 
      
      !!< Copy the latest flux and eigenvalue to the iter values 
      
      type(particle_type), intent(inout) :: particle
                  
      ! set the eigenvalue iter
      particle%keff%keff_iter = particle%keff%keff_new
      
      ! set the iter flux
      call copy_to_iter_values_particle_flux(particle)
      
      ! set the old flux, this is needed for the generic assemble-solve
      call copy_to_old_values_particle_flux(particle)
            
   end subroutine copy_to_iter_values_eig 

   ! --------------------------------------------------------------------------
      
   subroutine power_iteration_flux_acceleration(particle, &
                                                power_iteration_options)
      
      !!< Power iteration flux acceleration schemes
      
      type(particle_type), intent(inout) :: particle
      type(power_iteration_options_type), intent(in) :: power_iteration_options
      
      ! local variables
      real :: relax
      integer :: g
      integer :: number_of_energy_groups
      type(scalar_field), pointer :: particle_flux 
      type(scalar_field), pointer :: particle_flux_iter
      type(scalar_field) :: particle_flux_tmp

      acceleration_name: if (trim(power_iteration_options%flux_acceleration) == 'none') then
         
         return
      
      else if (trim(power_iteration_options%flux_acceleration) == 'Relaxation') then acceleration_name
      
         ! simple relaxation acceleration
         relax = power_iteration_options%flux_relaxation_coeff

         ! find the number of energy groups
         call find_total_number_energy_groups(trim(particle%option_path), &
                                              number_of_energy_groups)

         ! set the acceleted flux
         group_loop: do g = 1,number_of_energy_groups
         
            call extract_flux_group_g(particle, &
                                      g, &
                                      particle_flux = particle_flux, &
                                      particle_flux_iter = particle_flux_iter)
            
            call allocate(particle_flux_tmp, &
                          particle_flux_iter%mesh)
                                      
            call set(particle_flux_tmp, &
                     particle_flux_iter)

            call scale(particle_flux_tmp, &
                       1.0 - relax)
            
            call scale(particle_flux, &
                       relax)
            
            call addto(particle_flux, &
                       particle_flux_tmp)
            
            call deallocate(particle_flux_tmp)
            
         end do group_loop
               
      end if acceleration_name
      
   end subroutine power_iteration_flux_acceleration   

   ! --------------------------------------------------------------------------

   subroutine calculate_eigenvalue(particle) 
      
      !!< Calculate the latest eigenvalue as:
      !!< keff_new = keff_iter*(integral_dv(production_source*production_source)/integral_dv(production_source*production_source_iter))
      
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
      type(vector_field), pointer :: positions 
      type(scalar_field), pointer :: particle_flux 
      type(scalar_field), pointer :: particle_flux_iter
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
                                      particle_flux_iter = particle_flux_iter)
            
            call zero(production_coeff)
            
            ! form the production coeff field for this energy group
            call form(particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                      particle%particle_radmat, &
                      production_coeff, &
                      g, &
                      component = 'production')
            
            ! set the fields to integrate for k_top
            scalar_fields(1)%ptr => production_coeff
            scalar_fields(2)%ptr => particle_flux
            scalar_fields(3)%ptr => production_coeff
            scalar_fields(4)%ptr => particle_flux
            
            k_top = k_top + fields_integral(scalar_fields, &
                                            positions)

            ! set the fields to integrate for k_bottom
            scalar_fields(1)%ptr => production_coeff
            scalar_fields(2)%ptr => particle_flux
            scalar_fields(3)%ptr => production_coeff
            scalar_fields(4)%ptr => particle_flux_iter
            
            k_bottom = k_bottom + fields_integral(scalar_fields, &
                                                  positions)
            
         end do group_loop
         
         ! deallocate this group set production coeff field
         call deallocate(production_coeff)
      
      end do energy_group_set_loop

      if (associated(scalar_fields)) deallocate(scalar_fields)
      
      ! form the latest keff estimate
      particle%keff%keff_new = particle%keff%keff_iter*k_top/k_bottom

      ewrite(1,*) 'Keff: ',particle%keff%keff_new

   end subroutine calculate_eigenvalue  

   ! --------------------------------------------------------------------------
      
   subroutine power_iteration_keff_acceleration(particle, &
                                                power_iteration_options, &
                                                ipower)
      
      !!< Power iteration keff acceleration schemes
      
      type(particle_type), intent(inout) :: particle
      type(power_iteration_options_type), intent(in) :: power_iteration_options
      integer, intent(in) :: ipower
      
      ! local variables
      real, save :: keff1, keff2, keff3
      real :: delta_k
      
      acceleration_name: if (trim(power_iteration_options%keff_acceleration) == 'none') then
         
         return

      else if (trim(power_iteration_options%keff_acceleration) == 'AitkensDeltaSquared') then acceleration_name
      
         ! initialise
         if (ipower == 1) then
            
            keff1 = 0.0
            keff2 = 0.0
            keff3 = 0.0
         
         end if 
         
         ! accelerate
         if (ipower > 3) then
         
            keff1 = particle%keff%keff_new
                        
            delta_k = ( (keff1 - keff2)**2.0 ) / &
                      ( keff3 - 2.0*keff2 + keff1 )
            
            particle%keff%keff_new = particle%keff%keff_new - &
                                     delta_k
            
         end if 
         
         ! back save eigenvalues
         keff3 = keff2
         keff2 = keff1
         keff1 = particle%keff%keff_new
         
         ewrite(1,*) 'AitkensDeltaSquared accelerated Keff: ',particle%keff%keff_new
      
      end if acceleration_name
      
   end subroutine power_iteration_keff_acceleration   

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
      abs_difference_keff = abs(particle%keff%keff_new - particle%keff%keff_iter)
      
      ! find the keff relative difference
      find_change_keff: if (particle%keff%keff_iter > 0.0) then
         
         rel_difference_keff = abs_difference_keff / abs(particle%keff%keff_iter)
         
      else find_change_keff
         
         FLAbort('ERROR the radiation eigenvalue Keff iter is not > 0.0')
                  
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
