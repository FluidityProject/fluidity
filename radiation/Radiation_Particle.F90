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

module radiation_particle
   
   !!< Radiation particle type module
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module   
   
   use radiation_particle_data_type
   use radiation_materials
   use radiation_materials_interpolation
   
   implicit none
   
   private 

   public :: particle_type, &
             keff_type, &
             create, &
             destroy
        
   interface set
      module procedure set_all_keff
   end interface set

   interface destroy
      module procedure particles_destroy, &
                       particle_destroy
   end interface destroy

   interface create
      module procedure particles_create_from_options, &
                       particle_create_from_options
   end interface create
         
contains

   ! --------------------------------------------------------------------------

   subroutine particles_create_from_options(states_all, &
                                            particles) 
      
      !!< Create each particle type.
      
      type(state_type), dimension(:), intent(in) :: states_all
      type(particle_type), dimension(:), allocatable, intent(inout) :: particles
            
      ! local variable
      integer :: p 
      
      ewrite(1,*) 'Create each particle type'
                     
      allocate(particles(option_count('/embedded_models/radiation/particle_type')))
            
      particle_type_loop: do p = 1,option_count('/embedded_models/radiation/particle_type')
             
         ! set the particle option path
         particles(p)%option_path = '/embedded_models/radiation/particle_type['//int2str(p - 1)//']'
         
         ! set the particle name
         call get_option(trim(particles(p)%option_path)//'/name',particles(p)%name)
         
         ! create each particle type
         call create(states_all, &
                     particles(p))
                           
      end do particle_type_loop
            
      ewrite(1,*) 'Finished creating each particle'
      
   end subroutine particles_create_from_options

   ! --------------------------------------------------------------------------

   subroutine particle_create_from_options(states_all, &
                                           particle) 
      
      !!< Create a particle type.
      
      type(state_type), dimension(:), intent(in), target :: states_all
      type(particle_type), intent(inout) :: particle
                  
      ewrite(1,*) 'Create particle type ',trim(particle%name)
      
      ! set the state pointer for this particle
      ! currently this is set to states_all(1) until particles actually
      ! have their own state
      check_associated: if (associated(particle%state)) then
      
         nullify(particle%state)
      
      else check_associated

         allocate(particle%state)
      
      end if check_associated
            
      particle%state => states_all(1)
      
      ! create the radmats component   
      call create(particle%particle_radmat, &
                  trim(particle%option_path), &
                  trim(particle%name))
      
      ! read the radmats component   
      call particle_radmat_read(particle%particle_radmat)
      
      ! create the radmats_ii component   
      call create(particle%particle_radmat_ii, &
                  particle%particle_radmat, &
                  particle%state)
      
      ! set the keff components to 1.0
      call set(particle%keff, &
               value = 1.0)
                  
      ewrite(1,*) 'Finished creating particle type ',trim(particle%name)
      
   end subroutine particle_create_from_options

   ! --------------------------------------------------------------------------
   
   subroutine set_all_keff(keff, &
                           value)
      
      !!< Set all the keff components to value
      
      type(keff_type), intent(out) :: keff
      real, intent(in) :: value
      
      keff%keff_new          = value
      keff%keff_iter         = value
      keff%keff_iter_coupled = value
      keff%keff_iter_control = value
   
   end subroutine set_all_keff
   
   ! --------------------------------------------------------------------------

   subroutine particles_destroy(particles) 
   
      !!< Destroy the particles data type (deallocate, zero, unitialise) 

      type(particle_type), dimension(:), allocatable, intent(inout) :: particles
      
      ! local variable
      integer :: p 
      
      ewrite(1,*) 'Destroy each particle type'
      
      particle_type_loop: do p = 1,size(particles)
      
         call destroy(particles(p))
               
      end do particle_type_loop 
      
      if (allocated(particles)) deallocate(particles)

      ewrite(1,*) 'Finished destroying each particle type'
          
   end subroutine particles_destroy

   ! --------------------------------------------------------------------------

   subroutine particle_destroy(particle) 
   
      !!< Destroy a particle data type (deallocate, zero, unitialise) 

      type(particle_type), intent(inout) :: particle
            
      ewrite(1,*) 'Destroy particle type ',trim(particle%name)
      
      ! destroy the radmats component      
      call destroy(particle%particle_radmat)
      
      ! destroy the radmats_ii component   
      call destroy(particle%particle_radmat_ii)
      
      ! unitialise the option path and name
      particle%name = '/uninitialised_name/'
      
      particle%option_path = 'uninitialised_path'
      
      ! nullify the state pointer
      check_associated: if (associated(particle%state)) then
      
         nullify(particle%state)
      
      end if check_associated      
      
      ! set the created flag
      particle%created = .false.

      ewrite(1,*) 'Finished destroying particle type'
                      
   end subroutine particle_destroy

   ! --------------------------------------------------------------------------

end module radiation_particle
