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

   use radiation_materials
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
                                         np_radmat, &
                                         np_radmat_ii, &
                                         keff) 
   
      !!< Solve a radiation eigenvalue problem via a power iteration

      type(state_type), intent(inout) :: state
      type(np_radmat_type), intent(in) :: np_radmat
      type(np_radmat_ii_type), intent(in) :: np_radmat_ii
      real, intent(inout) :: keff
      
      ! local variables
      integer :: ipower
      integer :: number_all_convereged_pass
      integer :: number_of_energy_groups
      real :: keff_old
      logical :: power_iteration_converged
      character(len=OPTION_PATH_LEN) :: np_radmat_name
      type(power_iteration_options_type) :: power_iteration_options
      type(scalar_field), pointer :: keff_field
      
      ewrite(1,*) 'In eigenvalue_power_iteration'
        
      ! initialise the keff
      keff = 1.0
      
      ! intialise the number of all convergence passes in power iterations
      number_all_convereged_pass = 0
      
      call get_power_iteration_options(trim(np_radmat%option_path), &
                                       np_radmat_name, &
                                       number_of_energy_groups, &
                                       power_iteration_options)                                       
                  
      power_iteration: do ipower = 1,power_iteration_options%max_power_iteration
         
         ewrite(1,*) 'Power iteration: ',ipower
         
         call copy_to_old_values_eig(keff, &
                                     keff_old, &
                                     state, &
                                     trim(np_radmat_name), &
                                     number_of_energy_groups)
                  
         call scatter_iteration(np_radmat, &
                                np_radmat_ii, &
                                state, &
                                trim(np_radmat_name), &
                                number_of_energy_groups, & 
                                keff = keff)
                           
         call calculate_eigenvalue(keff, &
                                   keff_old, &
                                   state, &
                                   np_radmat, &
                                   np_radmat_ii, &
                                   trim(np_radmat_name), &
                                   number_of_energy_groups)

         call check_power_iteration_convergence(power_iteration_converged, &
                                                keff, &
                                                keff_old, &
                                                number_all_convereged_pass, &
                                                state, &
                                                trim(np_radmat_name), &
                                                number_of_energy_groups, &
                                                power_iteration_options)
                  
         ewrite(1,*) 'keff: ',keff
         
         if (power_iteration_converged) exit power_iteration
                  
      end do power_iteration
      
      ! set the Keff into the constant scalar field for output
      keff_field => extract_scalar_field(state,'NeutralParticleKeff'//trim(np_radmat_name))    
       
      call set(keff_field,keff) 
       
      need_to_terminate: if ((ipower == power_iteration_options%max_power_iteration) .and.  &
                             power_iteration_options%terminate_if_not_converged_power .and. &
                             (.not. power_iteration_converged)) then
         
         FLExit('Terminating as power iteration finished and solution NOT converged')
         
      end if need_to_terminate
            
   end subroutine eigenvalue_power_iteration

   ! --------------------------------------------------------------------------
   
   subroutine get_power_iteration_options(np_radmat_option_path, &
                                          np_radmat_name, &
                                          number_of_energy_groups, &
                                          power_iteration_options)
      
      !!< Get the power iteration options
      
      character(len=*), intent(in) :: np_radmat_option_path
      character(len=*), intent(out) :: np_radmat_name
      integer, intent(out) :: number_of_energy_groups
      type(power_iteration_options_type), intent(out) :: power_iteration_options
           
      ! local variable
      character(len=OPTION_PATH_LEN) :: power_iteration_option_path
      
      ! get the np_radmat name
      call get_option(trim(np_radmat_option_path)//'/name',np_radmat_name)
      
      ! get the number of energy groups
      call get_option(trim(np_radmat_option_path)//'/method/number_of_energy_groups',number_of_energy_groups)
                  
      ! get the maximum number of power iterations for this np 
      power_iteration_option_path = trim(np_radmat_option_path)//'/eigenvalue_run/power_iteration'
      call get_option(trim(power_iteration_option_path)//'/maximum',power_iteration_options%max_power_iteration)

      ! get the keff power iteration tolerance
      call get_option(trim(power_iteration_option_path)//'/keff_tolerance',power_iteration_options%keff_tolerance)

      ! get the flux power iteration tolerance
      call get_option(trim(power_iteration_option_path)//'/flux_tolerance',power_iteration_options%flux_tolerance)

      power_iteration_options%terminate_if_not_converged_power = have_option(trim(power_iteration_option_path)//'/terminate_if_not_converged')         
            
   end subroutine get_power_iteration_options

   ! --------------------------------------------------------------------------

   subroutine copy_to_old_values_eig(keff, &
                                     keff_old, &
                                     state, &
                                     np_radmat_name, &
                                     number_of_energy_groups) 
      
      !!< Copy the latest flux and eigenvalue to the old values 
      
      real, intent(in) :: keff
      real, intent(out) :: keff_old
      type(state_type), intent(inout) :: state
      character(len=*), intent(in) :: np_radmat_name
      integer, intent(in) :: number_of_energy_groups
                  
      ! set the eigenvalue old
      keff_old = keff
      
      call copy_to_old_values_np_flux(state, &
                                      np_radmat_name, &
                                      number_of_energy_groups)     
            
   end subroutine copy_to_old_values_eig 

   ! --------------------------------------------------------------------------

   subroutine calculate_eigenvalue(keff, &
                                   keff_old, &
                                   state,&
                                   np_radmat, &
                                   np_radmat_ii, &
                                   np_radmat_name, &
                                   number_of_energy_groups) 
      
      !!< Calculate the latest eigenvalue as:
      !!< keff_new = keff_old*(integral_dv(production_source*production_source)/integral_dv(production_source*production_source_old))
      
      real, intent(out) :: keff
      real, intent(in) :: keff_old
      type(state_type), intent(in) :: state      
      type(np_radmat_type), intent(in) :: np_radmat
      type(np_radmat_ii_type), intent(in) :: np_radmat_ii
      character(len=*), intent(in) :: np_radmat_name      
      integer, intent(in) :: number_of_energy_groups      
      
      ! local variables
      integer :: status
      integer :: g
      integer :: vele
      real :: k_top
      real :: k_bottom
      real, dimension(:,:), allocatable :: mass_matrix_vele
      real, dimension(:), allocatable :: detwei_vele
      type(vector_field), pointer :: positions 
      type(radmat_type) :: radmat_vele 
      type(scalar_field_pointer), dimension(:), pointer :: np_flux 
      type(scalar_field_pointer), dimension(:), pointer :: np_flux_old

      ! intialise variables that are summed up
      k_top    = 0.0
      k_bottom = 0.0
      
      ! extract the positions - assuming each group has the same topological mesh now
      positions => extract_vector_field(state, "Coordinate", stat=status)  
      
      ! extract the neutral particle flux fields for all energy groups
      call extract_flux_all_group(state, &
                                  trim(np_radmat_name), & 
                                  number_of_energy_groups, &
                                  np_flux = np_flux, &
                                  np_flux_old = np_flux_old)
      
      ! integrate each energy group
      group_loop: do g = 1,number_of_energy_groups

         ! loop the volume elements to perform the integration
         velement_loop: do vele = 1,ele_count(np_flux(g)%ptr)
                                                                      
            ! get the vele production cross section for all energy groups
            call form(radmat_vele, &
                      np_radmat_ii, &
                      np_radmat, &
                      vele, &
                      form_production = .true.)
        
            ! allocate the jacobian transform and gauss weight array for this vele
            allocate(detwei_vele(ele_ngi(np_flux(g)%ptr,vele)))
         
            ! allocate the local mass matrix for this vele
            allocate(mass_matrix_vele(ele_loc(np_flux(g)%ptr,vele),ele_loc(np_flux(g)%ptr,vele)))
                 
            ! form the velement jacobian transform and gauss weight
            call transform_to_physical(positions, vele, detwei = detwei_vele)
                  
            ! form the mass matrix
            mass_matrix_vele = shape_shape(ele_shape(np_flux(g)%ptr,vele),ele_shape(np_flux(g)%ptr,vele),detwei_vele)
                                    
            k_top = k_top + dot_product(ele_val(np_flux(g)%ptr,vele),matmul(mass_matrix_vele,ele_val(np_flux(g)%ptr,vele)))*radmat_vele%production(g)**2

            k_bottom = k_bottom + dot_product(ele_val(np_flux(g)%ptr,vele),matmul(mass_matrix_vele,ele_val(np_flux_old(g)%ptr,vele)))*radmat_vele%production(g)**2
         
            deallocate(detwei_vele)
            deallocate(mass_matrix_vele)
                                 
            call destroy(radmat_vele)
         
         end do velement_loop

      end do group_loop
      
      ! form the latest keff estimate
      keff = keff_old*k_top/k_bottom
      
      ! deallocate the np flux pointer fields
      call deallocate_flux_all_group(np_flux = np_flux, &
                                     np_flux_old = np_flux_old)
            
   end subroutine calculate_eigenvalue  

   ! --------------------------------------------------------------------------

   subroutine check_power_iteration_convergence(power_iteration_converged, &
                                                keff, &
                                                keff_old, &
                                                number_all_convereged_pass, &
                                                state, &
                                                np_radmat_name, &
                                                number_of_energy_groups, &
                                                power_iteration_options) 
   
      !!< Check the convergence of the power iterations via the keff and flux
      !!< Always insist on atleast 2 successive convergence passes   
      
      logical, intent(out) :: power_iteration_converged
      real, intent(in) :: keff
      real, intent(in) :: keff_old
      integer, intent(inout) :: number_all_convereged_pass 
      type(state_type), intent(in) :: state
      character(len=*), intent(in) :: np_radmat_name
      integer, intent(in) :: number_of_energy_groups      
      type(power_iteration_options_type), intent(in) :: power_iteration_options
      
      ! local variable
      real :: change_keff
      real :: max_change_flux
      logical :: keff_converged
      logical :: flux_converged
            
      power_iteration_converged = .true.
      
      ! find the keff change
      find_change_keff: if (keff_old > 0) then
         
         change_keff = abs(keff - keff_old)/abs(keff_old)
         
      else find_change_keff
         
         FLAbort('ERROR the radiation eigenvalue Keff is not > 0')
                  
      end if find_change_keff
      
      ! check the keff convergence
      check_keff: if (abs(keff - keff_old) < abs(keff_old*power_iteration_options%keff_tolerance)) then
         
         keff_converged = .true.
      
      else check_keff
         
         keff_converged = .false.
      
      end if check_keff 
      
      ewrite(1,*) 'change_keff,keff_converged: ',change_keff,keff_converged     
      
      ! check the flux convergence
      call check_np_flux_convergence(state, &
                                     trim(np_radmat_name), &
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
