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

module radiation_materials_interpolation_destroy
   
   use radiation_materials_interpolation_data_types

   implicit none
   
   private 

   public :: destroy
   
   interface destroy
      module procedure particle_radmat_ii_destroy, &
                       energy_group_set_all_destroy, &
                       energy_group_set_destroy, &
                       region_id_ii_all_destroy, &
                       region_id_ii_destroy, &
                       dataset_ii_destroy, &
                       ii_size_destroy
   end interface destroy

contains

   ! --------------------------------------------------------------------------

   subroutine particle_radmat_ii_destroy(particle_radmat_ii) 
   
      !!< Destroy the particle_radmat_ii type

      type(particle_radmat_ii_type), intent(inout) :: particle_radmat_ii
      
      energy_group_set_ii_deallocate: if (allocated(particle_radmat_ii%energy_group_set_ii)) then
            
         call destroy(particle_radmat_ii%energy_group_set_ii)
         
      end if energy_group_set_ii_deallocate

      ! set the flags      
      particle_radmat_ii%created = .false.
      
      particle_radmat_ii%formed = .false.
      
   end subroutine particle_radmat_ii_destroy

   ! --------------------------------------------------------------------------
   
   subroutine energy_group_set_all_destroy(energy_group_set_ii)
      
      !!< Destroy all the energy_group_set_ii type
      
      type(energy_group_set_ii_type), dimension(:), allocatable, intent(inout) :: energy_group_set_ii
      
      ! local variables
      integer :: g_set
      
      energy_group_set_loop: do g_set = 1,size(energy_group_set_ii)
         
         call destroy(energy_group_set_ii(g_set))
         
      end do energy_group_set_loop
      
      deallocate(energy_group_set_ii)
   
   end subroutine energy_group_set_all_destroy

   ! --------------------------------------------------------------------------
   
   subroutine energy_group_set_destroy(energy_group_set_ii)
      
      !!< Destroy the energy_group_set_ii type
      
      type(energy_group_set_ii_type), intent(inout) :: energy_group_set_ii
            
      region_id_ii_deallocate: if (allocated(energy_group_set_ii%region_id_ii)) then
         
         call destroy(energy_group_set_ii%region_id_ii)
                  
      end if region_id_ii_deallocate
      
      call destroy(energy_group_set_ii%ii_size)
   
   end subroutine energy_group_set_destroy

   ! --------------------------------------------------------------------------

   subroutine region_id_ii_all_destroy(region_id_ii)
      
      !!< Destroy all the region_id_ii type 
       
      type(region_id_ii_type), dimension(:), allocatable, intent(inout) :: region_id_ii
      
      ! local variables
      integer :: n
         
      node_loop: do n = 1,size(region_id_ii) 
                     
         call destroy(region_id_ii(n))
         
      end do node_loop
         
      deallocate(region_id_ii)
      
   end subroutine region_id_ii_all_destroy

   ! --------------------------------------------------------------------------

   subroutine region_id_ii_destroy(region_id_ii)

      !!< Destroy the region_id_ii type

      type(region_id_ii_type), intent(inout) :: region_id_ii
      
      call destroy(region_id_ii%dataset_ii)
      
   end subroutine region_id_ii_destroy

   ! --------------------------------------------------------------------------

   subroutine dataset_ii_destroy(dataset_ii)

      !!< Destroy the dataset_ii type 

      type(dataset_ii_type), intent(inout) :: dataset_ii
            
      fraction_deallocate: if (allocated(dataset_ii%physical_radmat_ii%fraction)) then
               
         deallocate(dataset_ii%physical_radmat_ii%fraction)
            
      end if fraction_deallocate

      coordinate_deallocate: if (allocated(dataset_ii%physical_radmat_ii%radmat_base_coordinate)) then
               
         deallocate(dataset_ii%physical_radmat_ii%radmat_base_coordinate)
            
      end if coordinate_deallocate
      
   end subroutine dataset_ii_destroy

   ! --------------------------------------------------------------------------

   subroutine ii_size_destroy(ii_size)
      
      !!< Destroy the ii_size type
      
      type(ii_size_type), intent(inout) :: ii_size
      
      region_id_ii_size_deallocate: if (allocated(ii_size%region_id_ii_size)) then
         
         deallocate(ii_size%region_id_ii_size)
               
      end if region_id_ii_size_deallocate
      
      ! set the flag
      ii_size%size_set = .false.
      
   end subroutine ii_size_destroy

   ! --------------------------------------------------------------------------

end module radiation_materials_interpolation_destroy
