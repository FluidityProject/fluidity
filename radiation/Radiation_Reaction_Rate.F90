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

module radiation_reaction_rate

   !!< This module contains procedures associated with calculating particle reaction rates 
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use fields
   
   use diagnostic_integrate_fields
   
   use radiation_particle
   use radiation_materials_interpolation
   use radiation_extract_flux_field

   implicit none
   
   private 

   public :: calculate_reaction_rate

contains

   ! --------------------------------------------------------------------------

   subroutine calculate_reaction_rate(reaction_rate_name, &
                                      state, &
                                      particle, &
                                      reaction_rate, &
                                      domain_symmetry_factor, &
                                      region_ids) 
   
      !!< Calculate a particle reaction rate 
      !!< If region_ids is present then only the corresponding positions field regions are integrated
      !!< else the whole domain is included
      
      character(len=*),                  intent(in)  :: reaction_rate_name
      type(state_type),                  intent(in)  :: state
      type(particle_type),               intent(in)  :: particle
      real,                              intent(out) :: reaction_rate
      integer,                           intent(in), optional  :: domain_symmetry_factor
      integer,             dimension(:), intent(in), optional :: region_ids
      
      ! local variables
      integer :: status
      integer :: g
      integer :: g_set
      integer :: g_global
      integer :: number_of_energy_groups
      integer :: number_of_energy_group_set
      real :: reaction_rate_group
      type(vector_field), pointer :: positions 
      type(scalar_field), pointer :: particle_flux 
      type(scalar_field), target :: reaction_rate_coeff      
      type(mesh_type), pointer :: material_fn_space
      type(scalar_field_pointer), dimension(:), pointer :: scalar_fields 
      character(len=OPTION_PATH_LEN) :: positions_mesh_name
      character(len=OPTION_PATH_LEN) :: material_fn_space_name
      character(len=OPTION_PATH_LEN) :: energy_group_set_path
      
      ewrite(1,*) 'Calculate particle reaction rate ',trim(reaction_rate_name), &
                  ' for particle ',trim(particle%name)
     
      ! intialise the reaction_rate
      reaction_rate = 0.0
      
      ! deduce the number of energy group sets
      number_of_energy_group_set = option_count(trim(particle%option_path)//'/energy_discretisation/energy_group_set')
      
      ! initialise the global group counter
      g_global = 0

      ! allocate the fields to integrate for
      allocate(scalar_fields(2))

      allocate(scalar_fields(1)%ptr)
      allocate(scalar_fields(2)%ptr)
      
      energy_group_set_loop: do g_set = 1,number_of_energy_group_set

         ! set the positions mesh name for this group set - currently assume all group sets have same positions
         positions_mesh_name = 'Coordinate'
         
         ! extract the positions 
         positions => extract_vector_field(state, trim(positions_mesh_name), stat=status)  
         
         if (status /= 0) FLAbort('Could not extract positions mesh for radiation flux normalisation')
         
         ! set the energy_group_set path
         energy_group_set_path = trim(particle%option_path)//'/energy_discretisation/energy_group_set['//int2str(g_set - 1)//']'
         
         ! get the material fn space name for this group set
         call get_option(trim(energy_group_set_path)//'/angular_discretisation/method/parity/angular_moment_set[0]/mesh/name',material_fn_space_name)
      
         ! extract the material fn_space of this energy group set of this particle type 
         material_fn_space => extract_mesh(state, trim(material_fn_space_name))

         ! get the number_energy_groups within this set
         call get_option(trim(energy_group_set_path)//'/number_of_energy_groups',number_of_energy_groups)
         
         ! allocate the reaction rate field for this group set
         call allocate(reaction_rate_coeff, &
                       material_fn_space, & 
                       'ParticleNormReactionRateCoeff')
                           
         ! integrate each energy group within this group set
         group_loop: do g = 1,number_of_energy_groups
            
            g_global = g_global + 1
            
            ! extract the group particle flux
            call extract_flux_group_g(state, &
                                      trim(trim(particle%name)), & 
                                      g_global, &
                                      particle_flux = particle_flux)
            
            call zero(reaction_rate_coeff)
            
            ! form the necessary reaction rate coeff field for this energy group            
            form_reaction_rate_coeff: if (trim(reaction_rate_name) == 'flux') then
            
              ! no reaction rate coeff so set to 1.0
              call set(reaction_rate_coeff,1.0)
            
            else form_reaction_rate_coeff
               
               ! the reaction rate names are the same as the material ii names         
               call form(material_fn_space, &
                         particle%particle_radmat_ii%energy_group_set_ii(g_set), &
                         particle%particle_radmat, &
                         reaction_rate_coeff, &
                         g, &
                         trim(reaction_rate_name))
            
            end if form_reaction_rate_coeff 
            
            ! integrate the material_coeff*solution over the necessary domain defined by positions mesh
            
            ! set the fields to be integrated as a product
            scalar_fields(1)%ptr => reaction_rate_coeff
            scalar_fields(2)%ptr => particle_flux

            call integrate(scalar_fields, &
                           positions, &
                           reaction_rate_group, &
                           region_ids = region_ids)


            reaction_rate = reaction_rate + reaction_rate_group
          
         end do group_loop
         
         ! deallocate this group set reaction_rate_coeff field
         call deallocate(reaction_rate_coeff)
      
      end do energy_group_set_loop

      if (associated(scalar_fields)) deallocate(scalar_fields)
      
      ! take account of any implicit domain symmetry
      if (present(domain_symmetry_factor)) reaction_rate = reaction_rate * domain_symmetry_factor

   end subroutine calculate_reaction_rate

   ! --------------------------------------------------------------------------

end module radiation_reaction_rate
