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

module radiation_materials_interpolation_form_radmat_vele

   !!< Procedures associated with the forming of the interpolated radiation material set for 
   !!< a particular volume element from the already determined interpolation instructions.

   use futils   
   use global_parameters, only : OPTION_PATH_LEN
  
   use radiation_materials  
   use radiation_materials_interpolation_data_types

   implicit none
   
   private 

   public :: form
   
   interface form
      module procedure form_radmat_vele
   end interface form
      
contains

   ! --------------------------------------------------------------------------

   subroutine form_radmat_vele(radmat_vele, &
                               particle_radmat_ii, &
                               particle_radmat, &
                               vele, &
                               form_all, &
                               form_total, &
                               form_absorption , &
                               form_scatter , &
                               form_removal , &
                               form_transport, &
                               form_diffusion, &
                               form_fission, &
                               form_production, &
                               form_power, &
                               form_energy_released_per_fission, &
                               form_particle_released_per_fission , &
                               form_prompt_spectrum, &
                               form_velocity , &
                               form_beta) 
   
      !!< Form the interpolated radmat for this vele using the already formed intructions
      !!< Optional arguments are used to manipulate which components of radmat_vele are formed
      !!< The formed components are allocated, zeroed then set via interpolation of the databases

      type(radmat_type), intent(out) :: radmat_vele 
      type(particle_radmat_ii_type), intent(in) :: particle_radmat_ii  
      type(particle_radmat_type), intent(in) :: particle_radmat 
      integer, intent(in) :: vele 
      logical, intent(in), optional :: form_all
      logical, intent(in), optional :: form_total
      logical, intent(in), optional :: form_absorption
      logical, intent(in), optional :: form_scatter
      logical, intent(in), optional :: form_removal
      logical, intent(in), optional :: form_transport
      logical, intent(in), optional :: form_diffusion
      logical, intent(in), optional :: form_fission
      logical, intent(in), optional :: form_production
      logical, intent(in), optional :: form_power
      logical, intent(in), optional :: form_energy_released_per_fission
      logical, intent(in), optional :: form_particle_released_per_fission
      logical, intent(in), optional :: form_prompt_spectrum
      logical, intent(in), optional :: form_velocity
      logical, intent(in), optional :: form_beta
      
      ! local variables
      integer :: dmat
      integer :: pmat
      logical :: allocate_all
      logical :: allocate_total
      logical :: allocate_absorption
      logical :: allocate_scatter
      logical :: allocate_removal
      logical :: allocate_transport
      logical :: allocate_diffusion
      logical :: allocate_fission
      logical :: allocate_production
      logical :: allocate_power
      logical :: allocate_energy_released_per_fission
      logical :: allocate_particle_released_per_fission
      logical :: allocate_prompt_spectrum
      logical :: allocate_velocity
      logical :: allocate_beta
       
      ! allocate the radmat_vele - only the components as needed
      
      allocate_all                           = .false.
      allocate_total                         = .false.
      allocate_absorption                    = .false.
      allocate_scatter                       = .false.
      allocate_removal                       = .false.
      allocate_transport                     = .false.
      allocate_diffusion                     = .false.
      allocate_fission                       = .false.
      allocate_production                    = .false.
      allocate_power                         = .false.
      allocate_energy_released_per_fission   = .false.
      allocate_particle_released_per_fission = .false.
      allocate_prompt_spectrum               = .false.
      allocate_velocity                      = .false.
      allocate_beta                          = .false. 
      
      if (present(form_all))                           allocate_all                           = form_all
      if (present(form_total))                         allocate_total                         = form_total
      if (present(form_absorption))                    allocate_absorption                    = form_absorption
      if (present(form_scatter))                       allocate_scatter                       = form_scatter
      if (present(form_removal))                       allocate_removal                       = form_removal
      if (present(form_transport))                     allocate_transport                     = form_transport
      if (present(form_diffusion))                     allocate_diffusion                     = form_diffusion
      if (present(form_fission))                       allocate_fission                       = form_fission
      if (present(form_production))                    allocate_production                    = form_production
      if (present(form_power))                         allocate_power                         = form_power
      if (present(form_energy_released_per_fission))   allocate_energy_released_per_fission   = form_energy_released_per_fission
      if (present(form_particle_released_per_fission)) allocate_particle_released_per_fission = form_particle_released_per_fission
      if (present(form_prompt_spectrum))               allocate_prompt_spectrum               = form_prompt_spectrum
      if (present(form_velocity))                      allocate_velocity                      = form_velocity
      if (present(form_beta))                          allocate_beta                          = form_beta
            
      call allocate(radmat_vele, &
                    particle_radmat%particle_radmat_size%number_of_energy_groups, &
                    particle_radmat%max_number_of_scatter_moments, &
                    particle_radmat%particle_radmat_size%number_of_delayed_groups, &
                    allocate_all                           = allocate_all, &
                    allocate_total                         = allocate_total, &
                    allocate_absorption                    = allocate_absorption, &
                    allocate_scatter                       = allocate_scatter, &
                    allocate_removal                       = allocate_removal, &
                    allocate_transport                     = allocate_transport, &
                    allocate_diffusion                     = allocate_diffusion, &
                    allocate_fission                       = allocate_fission, &
                    allocate_production                    = allocate_production, &
                    allocate_power                         = allocate_power, &
                    allocate_energy_released_per_fission   = allocate_energy_released_per_fission, &
                    allocate_particle_released_per_fission = allocate_particle_released_per_fission, &
                    allocate_prompt_spectrum               = allocate_prompt_spectrum, &
                    allocate_velocity                      = allocate_velocity, &
                    allocate_beta                          = allocate_beta)

      ! zero the radmat_vele - only the components allocated      
       call zero(radmat_vele)
       
      ! now ADD into the radmat_vele the necessary interpolated material data
       
      ! use the region id ii
      region_id: if (allocated(particle_radmat_ii%region_id_vele_ii)) then
         
         dmat = particle_radmat_ii%region_id_vele_ii(vele)%dataset_vele_ii%dataset_radmat_number 
         
         pmat = particle_radmat_ii%region_id_vele_ii(vele)%dataset_vele_ii%physical_radmat_vele_ii%physical_radmat_number
         
         call form_radmat_vele_region_id_ii(radmat_vele, &
                                            particle_radmat_ii%region_id_vele_ii(vele)%dataset_vele_ii%physical_radmat_vele_ii, &
                                            particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat))  
               
      end if region_id             
      
   end subroutine form_radmat_vele

   ! --------------------------------------------------------------------------

   subroutine form_radmat_vele_region_id_ii(radmat_vele, &
                                            physical_radmat_vele_ii, &
                                            physical_radmat)
      
      !!< Use the physical_radmat_vele_ii (that acts on the associated physical material dummy variable) 
      !!< for this vele to ADD into the radmat_vele the necessary material data 
      !!< A special case exists if there is only one radmat within the physical_radmat
      
      type(radmat_type), intent(inout) :: radmat_vele
      type(physical_radmat_vele_ii_type), intent(in) :: physical_radmat_vele_ii
      type(physical_radmat_type), intent(in) :: physical_radmat
      
      ! hard wire in 1d interpolation :(
      one_dim: if (size(physical_radmat_vele_ii%fraction) == 1) then
         
         ! only find the radmat components that are allocated (if not needed then not allocated)
         
         ! this code below needs formulating better if possible as it is copy+paste a lot :(
         
         total: if (allocated(radmat_vele%total)) then
            
            if (size(physical_radmat%radmats) == 1) then
               
               radmat_vele%total = physical_radmat%radmats(1)%total
               
            else 
                              
               radmat_vele%total = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%total + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%total
            
            end if 
               
         end if total
         
         absorption: if (allocated(radmat_vele%absorption)) then
            
            if (size(physical_radmat%radmats) == 1) then
               
               radmat_vele%absorption = physical_radmat%radmats(1)%absorption
               
            else 
                              
               radmat_vele%absorption = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%absorption + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%absorption
            
            end if 
               
         end if absorption
         
         scatter: if (allocated(radmat_vele%scatter)) then
            
            if (size(physical_radmat%radmats) == 1) then
              
               radmat_vele%scatter = physical_radmat%radmats(1)%scatter
               
            else 
                              
               radmat_vele%scatter = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%scatter + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%scatter
            
            end if 
               
         end if scatter
         
         removal: if (allocated(radmat_vele%removal)) then
            
            if (size(physical_radmat%radmats) == 1) then
               
               radmat_vele%removal = physical_radmat%radmats(1)%removal
               
            else 
                              
               radmat_vele%removal = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%removal + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%removal
            
            end if 
               
         end if removal
         
         transport: if (allocated(radmat_vele%transport)) then
            
            if (size(physical_radmat%radmats) == 1) then
               
               radmat_vele%transport = physical_radmat%radmats(1)%transport
               
            else 
                              
               radmat_vele%transport = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%transport + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%transport
            
            end if 
               
         end if transport
         
         diffusion: if (allocated(radmat_vele%diffusion)) then
            
            if (size(physical_radmat%radmats) == 1) then
               
               radmat_vele%diffusion = physical_radmat%radmats(1)%diffusion
               
            else 
                              
               radmat_vele%diffusion = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%diffusion + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%diffusion
            
            end if 
               
         end if diffusion
         
         fission: if (allocated(radmat_vele%fission)) then
            
            if (size(physical_radmat%radmats) == 1) then
               
               radmat_vele%fission = physical_radmat%radmats(1)%fission
               
            else 
                              
               radmat_vele%fission = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%fission + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%fission
            
            end if 
               
         end if fission
         
         production: if (allocated(radmat_vele%production)) then
            
            if (size(physical_radmat%radmats) == 1) then
               
               radmat_vele%production = physical_radmat%radmats(1)%production
               
            else 
                              
               radmat_vele%production = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%production + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%production
            
            end if 
               
         end if production
         
         power: if (allocated(radmat_vele%power)) then
            
            if (size(physical_radmat%radmats) == 1) then
               
               radmat_vele%power = physical_radmat%radmats(1)%power
               
            else 
                              
               radmat_vele%power = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%power + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%power
            
            end if 
               
         end if power
         
         energy_released_per_fission: if (allocated(radmat_vele%energy_released_per_fission)) then
            
            if (size(physical_radmat%radmats) == 1) then
               
               radmat_vele%energy_released_per_fission = physical_radmat%radmats(1)%energy_released_per_fission
               
            else 
                              
               radmat_vele%energy_released_per_fission = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%energy_released_per_fission + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%energy_released_per_fission
            
            end if 
               
         end if energy_released_per_fission
         
         particle_released_per_fission: if (allocated(radmat_vele%particle_released_per_fission)) then
            
            if (size(physical_radmat%radmats) == 1) then
               
               radmat_vele%particle_released_per_fission = physical_radmat%radmats(1)%particle_released_per_fission
               
            else 
                              
               radmat_vele%particle_released_per_fission = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%particle_released_per_fission + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%particle_released_per_fission
            
            end if 
               
         end if particle_released_per_fission
         
         prompt_spectrum: if (allocated(radmat_vele%prompt_spectrum)) then
            
            if (size(physical_radmat%radmats) == 1) then
               
               radmat_vele%prompt_spectrum = physical_radmat%radmats(1)%prompt_spectrum
               
            else 
                              
               radmat_vele%prompt_spectrum = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%prompt_spectrum + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%prompt_spectrum
            
            end if 
               
         end if prompt_spectrum
         
         velocity: if (allocated(radmat_vele%velocity)) then
            
            if (size(physical_radmat%radmats) == 1) then
               
               radmat_vele%velocity = physical_radmat%radmats(1)%velocity
               
            else 
                              
               radmat_vele%velocity = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%velocity + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%velocity
            
            end if 
               
         end if velocity
         
         beta: if (allocated(radmat_vele%beta)) then
            
            if (size(physical_radmat%radmats) == 1) then
               
               radmat_vele%beta = physical_radmat%radmats(1)%beta
               
            else 
                              
               radmat_vele%beta = (1.0 - physical_radmat_vele_ii%fraction(1))* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1))%beta + &
                                   physical_radmat_vele_ii%fraction(1)* &
                                   physical_radmat%radmats(physical_radmat_vele_ii%radmat_base_coordinate(1) + 1)%beta
            
            end if 
               
         end if beta
      
      else one_dim
      
         FLExit('Only 1 dimensional radiation material interpolation supported currently')
      
      end if one_dim
      
   end subroutine form_radmat_vele_region_id_ii

   ! --------------------------------------------------------------------------

end module radiation_materials_interpolation_form_radmat_vele
