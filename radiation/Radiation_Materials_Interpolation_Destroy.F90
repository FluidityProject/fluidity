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
      module procedure np_radmat_ii_destroy, &
                       region_id_vele_ii_all_destroy, &
                       region_id_vele_ii_destroy, &
                       dataset_vele_ii_destroy, &
                       np_radmat_ii_size_destroy
   end interface destroy

contains

   ! --------------------------------------------------------------------------

   subroutine np_radmat_ii_destroy(np_radmat_ii) 
   
      !!< Destroy the np_radmat_ii type

      type(np_radmat_ii_type), intent(inout) :: np_radmat_ii
            
      region_id_ii_deallocate: if (allocated(np_radmat_ii%region_id_vele_ii)) then
         
         call destroy(np_radmat_ii%region_id_vele_ii)
                  
      end if region_id_ii_deallocate
      
      call destroy(np_radmat_ii%np_radmat_ii_size)

      ! set the flags      
      np_radmat_ii%created = .false.
      
      np_radmat_ii%formed = .false.
      
   end subroutine np_radmat_ii_destroy

   ! --------------------------------------------------------------------------

   subroutine region_id_vele_ii_all_destroy(region_id_vele_ii)
      
      !!< Destroy the region_id_vele_ii type for all vele
       
      type(region_id_vele_ii_type), dimension(:), allocatable, intent(inout) :: region_id_vele_ii
      
      ! local variables
      integer :: vele
         
      vele_loop: do vele = 1,size(region_id_vele_ii) 
                     
         call destroy(region_id_vele_ii(vele))
         
      end do vele_loop
         
      deallocate(region_id_vele_ii)
      
   end subroutine region_id_vele_ii_all_destroy

   ! --------------------------------------------------------------------------

   subroutine region_id_vele_ii_destroy(region_id_vele_ii)

      !!< Destroy the region_id_vele_ii type for one vele

      type(region_id_vele_ii_type), intent(inout) :: region_id_vele_ii
      
      call destroy(region_id_vele_ii%dataset_vele_ii)
      
   end subroutine region_id_vele_ii_destroy

   ! --------------------------------------------------------------------------

   subroutine dataset_vele_ii_destroy(dataset_vele_ii)

      !!< Destroy the dataset_vele_ii type 

      type(dataset_vele_ii_type), intent(inout) :: dataset_vele_ii
            
      fraction_deallocate: if (allocated(dataset_vele_ii%physical_radmat_vele_ii%fraction)) then
               
         deallocate(dataset_vele_ii%physical_radmat_vele_ii%fraction)
            
      end if fraction_deallocate

      coordinate_deallocate: if (allocated(dataset_vele_ii%physical_radmat_vele_ii%radmat_base_coordinate)) then
               
         deallocate(dataset_vele_ii%physical_radmat_vele_ii%radmat_base_coordinate)
            
      end if coordinate_deallocate
      
   end subroutine dataset_vele_ii_destroy

   ! --------------------------------------------------------------------------

   subroutine np_radmat_ii_size_destroy(np_radmat_ii_size)
      
      !!< Destroy the np_radmat_ii_size type
      
      type(np_radmat_ii_size_type), intent(inout) :: np_radmat_ii_size
      
      region_id_ii_size_deallocate: if (allocated(np_radmat_ii_size%region_id_vele_ii_size)) then
         
         deallocate(np_radmat_ii_size%region_id_vele_ii_size)
               
      end if region_id_ii_size_deallocate
      
      ! set the flag
      np_radmat_ii_size%size_set = .false.
      
   end subroutine np_radmat_ii_size_destroy

   ! --------------------------------------------------------------------------

end module radiation_materials_interpolation_destroy
