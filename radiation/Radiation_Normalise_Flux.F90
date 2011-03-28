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

module radiation_normalise_flux

   !!< This module contains procedures associated with normalising the radiation particle flux 
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use fields
   
   use radiation_particle_data_type
   use radiation_extract_flux_field
   use radiation_reaction_rate
   use radiation_energy_group_set_tools

   implicit none
   
   private 

   public :: normalise_particle_flux
   
   type norm_options_type
      ! the norm value
      real :: value
      ! the norm name
      character(len=OPTION_PATH_LEN) :: name
      ! the domain symmetry factor
      integer :: domain_symmetry_factor
   end type norm_options_type

contains

   ! --------------------------------------------------------------------------

   subroutine normalise_particle_flux(particle) 
   
      !!< Normalise the particle flux

      type(particle_type), intent(in) :: particle
      
      ! local variables
      integer :: g
      integer :: number_of_energy_groups
      real :: total_reaction_rate
      real :: norm_factor
      type(scalar_field), pointer :: particle_flux 
      type(norm_options_type) :: norm_options
      character(len=OPTION_PATH_LEN) :: reaction_rate_name
      
      ewrite(1,*) 'Normalise particle flux'
      
      ! get the norm options
      call get_norm_options(particle%option_path, &
                            norm_options)
      
      ! form the reaction_rate_name
      reaction_name_if: if (trim(norm_options%name) == 'TotalFlux') then
         
         reaction_rate_name = 'flux'
      
      else if (trim(norm_options%name) == 'TotalProduction') then reaction_name_if

         reaction_rate_name = 'production'
      
      else if (trim(norm_options%name) == 'TotalFission') then reaction_name_if

         reaction_rate_name = 'fission'
      
      else if (trim(norm_options%name) == 'TotalPower') then reaction_name_if

         reaction_rate_name = 'power'
            
      else reaction_name_if
         
         FLAbort('Unknown particle normalisation name')
         
      end if reaction_name_if
      
      ! find the necessary total reaction rate
      call calculate_reaction_rate(trim(reaction_rate_name), &
                                   particle, &
                                   total_reaction_rate, &
                                   domain_symmetry_factor = norm_options%domain_symmetry_factor)

      ! form the normalisation factor
      norm_factor = norm_options%value / total_reaction_rate

      ! get the number of energy groups via summing the number within each energy group set      
      call find_total_number_energy_groups(trim(particle%option_path), &
                                           number_of_energy_groups)      

      ! normalise the flux for each group
      group_norm_loop: do g = 1,number_of_energy_groups     

         ! extract the group particle flux
         call extract_flux_group_g(particle, &
                                   g, &
                                   particle_flux = particle_flux)
         
         call scale(particle_flux, &
                    norm_factor)

      end do group_norm_loop

      ewrite(1,*) 'Finished normalise particle flux'
            
   end subroutine normalise_particle_flux

   ! --------------------------------------------------------------------------
   
   subroutine get_norm_options(particle_option_path, &
                               norm_options)
   
      !!< Get the normalisation options 
      
      character(len=*), intent(in) :: particle_option_path
      type(norm_options_type), intent(out) :: norm_options
      
      ! local variables
      character(len=OPTION_PATH_LEN) :: flux_normalisation_path
      
      ! set the flux normalisation path
      flux_normalisation_path = trim(particle_option_path)//'/equation/flux_normalisation'
      
      ! get the norm name
      call get_option(trim(flux_normalisation_path)//'/name', &
                      norm_options%name)
      
      ! get the norm value
      call get_option(trim(flux_normalisation_path)//'/value', &
                      norm_options%value)
      
      ! get the norm domain symmetry factor
      call get_option(trim(flux_normalisation_path)//'/domain_symmetry_factor', &
                      norm_options%domain_symmetry_factor, &
                      default = 1)
   
   end subroutine get_norm_options   
   
   ! --------------------------------------------------------------------------

end module radiation_normalise_flux
