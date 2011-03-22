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

   ! keep in this order, please:
   use quadrature
   use elements
   use sparse_tools
   use fields
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use parallel_tools

   use radiation_particle
   use radiation_materials_interpolation
   use radiation_solve_scatter_iteration
   use radiation_extract_flux_field
   use radiation_copy_flux_values
   use radiation_check_flux_convergence
   
   implicit none
   
   private 

   public :: eigenvalue_power_iteration
   
   type power_iteration_options_type
      integer :: max_power_iteration
      real :: keff_tolerance
      real :: flux_tolerance
      logical :: terminate_if_not_converged_power            
   end type power_iteration_options_type
   
contains

   ! --------------------------------------------------------------------------

   subroutine eigenvalue_power_iteration(state, &
                                         particle) 
   
      !!< Solve a radiation eigenvalue problem via a power iteration

      type(state_type), intent(inout) :: state
      type(particle_type), intent(inout) :: particle
      
      ! local variables
      integer :: ipower
      integer :: number_all_convereged_pass
      integer :: number_of_energy_groups
      logical :: power_iteration_converged
      type(power_iteration_options_type) :: power_iteration_options
      type(scalar_field), pointer :: keff_field
      
      ewrite(1,*) 'In eigenvalue_power_iteration'
      
      ! keff type already intialised to 1.0
      
      ! intialise the number of all convergence passes in power iterations
      number_all_convereged_pass = 0
      
      call get_power_iteration_options(particle, &
                                       number_of_energy_groups, &
                                       power_iteration_options)                                       
                  
      power_iteration: do ipower = 1,power_iteration_options%max_power_iteration
         
         ewrite(1,*) 'Power iteration: ',ipower
         
         call copy_to_old_values_eig(state, &
                                     particle, &
                                     number_of_energy_groups)
                  
         call scatter_iteration(particle, &
                                state, &
                                number_of_energy_groups, & 
                                invoke_eigenvalue_scatter_solve = .true.)
                           
         call calculate_eigenvalue(state, &
                                   particle)

         call check_power_iteration_convergence(power_iteration_converged, &
                                                number_all_convereged_pass, &
                                                state, &
                                                particle, &
                                                number_of_energy_groups, &
                                                power_iteration_options)
                           
         if (power_iteration_converged) exit power_iteration
                  
      end do power_iteration
      
      ! set the Keff into the constant scalar field for output
      keff_field => extract_scalar_field(state,'ParticleKeff'//trim(particle%name))    
       
      call set(keff_field,particle%keff%keff_new) 
       
      need_to_terminate: if ((ipower == power_iteration_options%max_power_iteration) .and.  &
                             power_iteration_options%terminate_if_not_converged_power .and. &
                             (.not. power_iteration_converged)) then
         
         FLExit('Terminating as power iteration finished and solution NOT converged')
         
      end if need_to_terminate
            
   end subroutine eigenvalue_power_iteration

   ! --------------------------------------------------------------------------
   
   subroutine get_power_iteration_options(particle, &
                                          number_of_energy_groups, &
                                          power_iteration_options)
      
      !!< Get the power iteration options
      
      type(particle_type), intent(in) :: particle
      integer, intent(out) :: number_of_energy_groups
      type(power_iteration_options_type), intent(out) :: power_iteration_options
           
      ! local variable
      integer :: g_set
      integer :: number_of_energy_group_set
      integer :: number_of_energy_groups_g_set
      character(len=OPTION_PATH_LEN) :: energy_group_set_path
      character(len=OPTION_PATH_LEN) :: power_iteration_option_path
      
      ! get the number of energy groups via summing the number within each energy group set

      ! deduce the number of energy group sets
      number_of_energy_group_set = option_count(trim(particle%option_path)//'/energy_discretisation/energy_group_set')
      
      number_of_energy_groups = 0
      
      energy_group_set_loop: do g_set = 1,number_of_energy_group_set
            
         ! set the energy_group_set path
         energy_group_set_path = trim(particle%option_path)//'/energy_discretisation/energy_group_set['//int2str(g_set - 1)//']'
            
         ! get the number_energy_groups within this set
         call get_option(trim(energy_group_set_path)//'/number_of_energy_groups',number_of_energy_groups_g_set)         
         
         number_of_energy_groups = number_of_energy_groups + number_of_energy_groups_g_set
         
      end do energy_group_set_loop
                  
      ! get the maximum number of power iterations for this np 
      power_iteration_option_path = trim(particle%option_path)//'/equation/power_iteration'
      call get_option(trim(power_iteration_option_path)//'/maximum',power_iteration_options%max_power_iteration)

      ! get the keff power iteration tolerance
      call get_option(trim(power_iteration_option_path)//'/keff_tolerance',power_iteration_options%keff_tolerance)

      ! get the flux power iteration tolerance
      call get_option(trim(power_iteration_option_path)//'/flux_tolerance',power_iteration_options%flux_tolerance)

      power_iteration_options%terminate_if_not_converged_power = have_option(trim(power_iteration_option_path)//'/terminate_if_not_converged')         
            
   end subroutine get_power_iteration_options

   ! --------------------------------------------------------------------------

   subroutine copy_to_old_values_eig(state, &
                                     particle, &
                                     number_of_energy_groups) 
      
      !!< Copy the latest flux and eigenvalue to the old values 
      
      type(state_type), intent(inout) :: state
      type(particle_type), intent(inout) :: particle      
      integer, intent(in) :: number_of_energy_groups
                  
      ! set the eigenvalue old
      particle%keff%keff_old = particle%keff%keff_new
      
      call copy_to_old_values_particle_flux(state, &
                                            trim(particle%name), &
                                            number_of_energy_groups)     
            
   end subroutine copy_to_old_values_eig 

   ! --------------------------------------------------------------------------

   subroutine calculate_eigenvalue(state,&
                                   particle) 
      
      !!< Calculate the latest eigenvalue as:
      !!< keff_new = keff_old*(integral_dv(production_source*production_source)/integral_dv(production_source*production_source_old))
      
      type(state_type), intent(in) :: state      
      type(particle_type), intent(inout) :: particle
      
      ! local variables
      integer :: status
      integer :: g
      integer :: g_set
      integer :: g_global
      integer :: vele
      integer :: number_of_energy_groups
      integer :: number_of_energy_group_set
      real :: k_top
      real :: k_bottom
      real, dimension(:,:), allocatable :: mass_matrix_vele
      real, dimension(:), allocatable :: detwei_vele
      real, dimension(:), allocatable :: production_val_quad_vele      
      type(vector_field), pointer :: positions 
      type(scalar_field), pointer :: particle_flux 
      type(scalar_field), pointer :: particle_flux_old
      type(scalar_field) :: production_coeff      
      type(mesh_type), pointer :: material_fn_space
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
      
      energy_group_set_loop: do g_set = 1,number_of_energy_group_set

         ! set the positions mesh name for this group set - currently assume all group sets have same positions
         positions_mesh_name = 'Coordinate'
         
         ! extract the positions 
         positions => extract_vector_field(state, trim(positions_mesh_name), stat=status)  

         ! set the energy_group_set path
         energy_group_set_path = trim(particle%option_path)//'/energy_discretisation/energy_group_set['//int2str(g_set - 1)//']'
         
         ! get the material fn space name for this group set
         call get_option(trim(energy_group_set_path)//'/angular_discretisation/method/parity/angular_moment_set[0]/mesh/name',material_fn_space_name)
      
         ! extract the material fn_space of this energy group set of this particle type 
         material_fn_space => extract_mesh(state, trim(material_fn_space_name))

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
            call extract_flux_group_g(state, &
                                      trim(trim(particle%name)), & 
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
            
            ! loop the volume elements to perform the integration
            velement_loop: do vele = 1,ele_count(particle_flux)
        
               ! allocate the jacobian transform and gauss weight array for this vele
               allocate(detwei_vele(ele_ngi(particle_flux,vele)))
         
               ! allocate the local mass matrix for this vele
               allocate(mass_matrix_vele(ele_loc(particle_flux,vele),ele_loc(particle_flux,vele)))
                 
               ! form the velement jacobian transform and gauss weight
               call transform_to_physical(positions, vele, detwei = detwei_vele)
               
               allocate(production_val_quad_vele(ele_ngi(production_coeff, vele)))
               
               production_val_quad_vele = 0.0
               
               ! get the production field values at the quadrature points
               production_val_quad_vele = ele_val_at_quad(production_coeff, vele)
               
               ! form the square of each of these values
               production_val_quad_vele = production_val_quad_vele**2
                  
               ! form the mass matrix
               mass_matrix_vele = shape_shape(ele_shape(particle_flux,vele),ele_shape(particle_flux,vele),detwei_vele*production_val_quad_vele)
                                    
               k_top = k_top + dot_product(ele_val(particle_flux,vele),matmul(mass_matrix_vele,ele_val(particle_flux,vele)))

               k_bottom = k_bottom + dot_product(ele_val(particle_flux,vele),matmul(mass_matrix_vele,ele_val(particle_flux_old,vele)))
         
               deallocate(detwei_vele)
               deallocate(mass_matrix_vele)
               deallocate(production_val_quad_vele)
                                          
            end do velement_loop
            
            ! sum the value of the processes
            call allsum(k_top)
            call allsum(k_bottom)

         end do group_loop
         
         ! deallocate this group set production coeff field
         call deallocate(production_coeff)
      
      end do energy_group_set_loop
      
      ! form the latest keff estimate
      particle%keff%keff_new = particle%keff%keff_old*k_top/k_bottom

      ewrite(1,*) 'Keff: ',particle%keff%keff_new     
                
   end subroutine calculate_eigenvalue  

   ! --------------------------------------------------------------------------

   subroutine check_power_iteration_convergence(power_iteration_converged, &
                                                number_all_convereged_pass, &
                                                state, &
                                                particle, &
                                                number_of_energy_groups, &
                                                power_iteration_options) 
   
      !!< Check the convergence of the power iterations via the keff and flux
      !!< Always insist on atleast 2 successive convergence passes   
      
      logical, intent(out) :: power_iteration_converged
      integer, intent(inout) :: number_all_convereged_pass 
      type(state_type), intent(in) :: state
      type(particle_type), intent(in) :: particle
      integer, intent(in) :: number_of_energy_groups      
      type(power_iteration_options_type), intent(in) :: power_iteration_options
      
      ! local variable
      real :: change_keff
      real :: max_change_flux
      logical :: keff_converged
      logical :: flux_converged
            
      power_iteration_converged = .true.
      
      ! find the keff change
      find_change_keff: if (particle%keff%keff_old > 0.0) then
         
         change_keff = abs(particle%keff%keff_new - particle%keff%keff_old)/ &
                       abs(particle%keff%keff_old)
         
      else find_change_keff
         
         FLAbort('ERROR the radiation eigenvalue Keff old is not > 0.0')
                  
      end if find_change_keff
      
      ! check the keff convergence
      check_keff: if (abs(particle%keff%keff_new - particle%keff%keff_old) <  &
                      abs(particle%keff%keff_old*power_iteration_options%keff_tolerance)) then
         
         keff_converged = .true.
      
      else check_keff
         
         keff_converged = .false.
      
      end if check_keff 
      
      ewrite(1,*) 'change_keff,keff_converged: ',change_keff,keff_converged     
      
      ! check the flux convergence
      call check_particle_flux_convergence(state, &
                                           trim(particle%name), &
                                           number_of_energy_groups, &
                                           power_iteration_options%flux_tolerance, &
                                           max_change_flux, &
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
