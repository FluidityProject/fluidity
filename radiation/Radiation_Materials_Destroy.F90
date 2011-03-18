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

module radiation_materials_destroy

   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud

   use radiation_materials_data_types     
   
   implicit none
   
   private 

   public :: destroy

   interface destroy
      module procedure particle_radmat_destroy, &
                       particle_radmat_size_destroy, &
                       delayed_lambda_spectrum_destroy, &
                       dataset_radmat_destroy, &
                       physical_radmat_destroy, &
                       radmat_destroy
   end interface destroy
     
contains

   ! --------------------------------------------------------------------------
   
   subroutine particle_radmat_destroy(particle_radmat)
      
      !!< Destroy any allocated memory associated with this particle_radmat,
      !!< set the option path to uninitialised and the flags as needed
    
      type(particle_radmat_type), intent(inout) :: particle_radmat
      
      ! local variables
      integer :: dmat
               
      particle_radmat%option_path = "/uninitialised_path/"

      particle_radmat%name = "/uninitialised_name/"
                  
      deallacate_data_set: if (allocated(particle_radmat%dataset_radmats)) then

         Data_set_loop: do dmat = 1,size(particle_radmat%dataset_radmats)
         
            call destroy(particle_radmat%dataset_radmats(dmat))
         
         end do Data_set_loop 
      
         deallocate(particle_radmat%dataset_radmats)
      
      end if deallacate_data_set
      
      call destroy(particle_radmat%delayed_lambda_spectrum)
      
      call destroy(particle_radmat%particle_radmat_size)
      
      particle_radmat%max_number_of_scatter_moments = 0
      
      particle_radmat%created = .false.
 
      particle_radmat%readin = .false.
           
   end subroutine particle_radmat_destroy
   
   ! --------------------------------------------------------------------------

   subroutine particle_radmat_size_destroy(particle_radmat_size)
   
      !! Destroy any allocated memory associated with this particle_radmat_size and set sizes to zero
      
      type(particle_radmat_size_type), intent(inout) :: particle_radmat_size
      
      if (allocated(particle_radmat_size%number_of_physical_radmats)) deallocate(particle_radmat_size%number_of_physical_radmats)
      if (allocated(particle_radmat_size%number_of_radmats))          deallocate(particle_radmat_size%number_of_radmats)
      if (allocated(particle_radmat_size%number_of_scatter_moments))  deallocate(particle_radmat_size%number_of_scatter_moments)
      
      particle_radmat_size%number_of_energy_groups       = 0
      particle_radmat_size%number_of_delayed_groups      = 0
      particle_radmat_size%total_number_radmats          = 0
      particle_radmat_size%total_number_physical_radmats = 0      
      particle_radmat_size%total_number_dataset_radmats  = 0
      particle_radmat_size%size_set                      = .false.

   end subroutine particle_radmat_size_destroy

   ! --------------------------------------------------------------------------
   
   subroutine delayed_lambda_spectrum_destroy(delayed_lambda_spectrum)
   
      type(delayed_lambda_spectrum_type), intent(inout) :: delayed_lambda_spectrum

      deallocate_lambda: if (allocated(delayed_lambda_spectrum%lambda)) then
            
         deallocate(delayed_lambda_spectrum%lambda)
            
      end if deallocate_lambda
         
      delayed_lambda_spectrum%lambda_set = .false.
         
      deallocate_Dspectrum: if (allocated(delayed_lambda_spectrum%spectrum)) then
            
         deallocate(delayed_lambda_spectrum%spectrum)
            
      end if deallocate_Dspectrum 

      delayed_lambda_spectrum%spectrum_set = .false.
   
   end subroutine delayed_lambda_spectrum_destroy
   
   ! --------------------------------------------------------------------------

   subroutine dataset_radmat_destroy(dataset_radmat)
   
      type(dataset_radmat_type), intent(inout) :: dataset_radmat
      
      ! local variables
      integer :: pmat
               
      dataset_radmat%option_path = "/uninitialised_path/"

      dataset_radmat%name = "/uninitialised_name/"

      dataset_radmat%file_name = "/uninitialised_name/"
         
      deallacate_physical_mat: if (allocated(dataset_radmat%physical_radmats)) then

         physical_radmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
         
            call destroy(dataset_radmat%physical_radmats(pmat))

         end do physical_radmat_loop
      
         deallocate(dataset_radmat%physical_radmats)
      
      end if deallacate_physical_mat
      
   end subroutine dataset_radmat_destroy
   
   ! --------------------------------------------------------------------------
   
   subroutine physical_radmat_destroy(physical_radmat)
   
      type(physical_radmat_type), intent(inout) :: physical_radmat
      
      ! local variables
      integer :: rmat
            
      physical_radmat%option_path = "/uninitialised_path/"

      physical_radmat%name = "/uninitialised_name/"

      deallacate_rmat: if (allocated(physical_radmat%radmats)) then
         
         radmat_loop: do rmat = 1,size(physical_radmat%radmats)
         
            call destroy(physical_radmat%radmats(rmat))
      
         end do radmat_loop
      
         deallocate(physical_radmat%radmats)
      
      end if deallacate_rmat
   
   end subroutine physical_radmat_destroy  
   
   ! --------------------------------------------------------------------------

   subroutine radmat_destroy(radmat)
   
      type(radmat_type), intent(inout) :: radmat
      
      if (allocated(radmat%total))                         deallocate(radmat%total)
      if (allocated(radmat%absorption))                    deallocate(radmat%absorption)
      if (allocated(radmat%scatter))                       deallocate(radmat%scatter)
      if (allocated(radmat%removal))                       deallocate(radmat%removal)
      if (allocated(radmat%transport))                     deallocate(radmat%transport)
      if (allocated(radmat%diffusion))                     deallocate(radmat%diffusion)
      if (allocated(radmat%fission))                       deallocate(radmat%fission)
      if (allocated(radmat%production))                    deallocate(radmat%production)
      if (allocated(radmat%power))                         deallocate(radmat%power)
      if (allocated(radmat%energy_released_per_fission))   deallocate(radmat%energy_released_per_fission)
      if (allocated(radmat%particle_released_per_fission)) deallocate(radmat%particle_released_per_fission)
      if (allocated(radmat%prompt_spectrum))               deallocate(radmat%prompt_spectrum)
      if (allocated(radmat%velocity))                      deallocate(radmat%velocity)
      if (allocated(radmat%beta))                          deallocate(radmat%beta)
            
      radmat%total_set                         = .false.  
      radmat%absorption_set                    = .false.  
      radmat%scatter_set                       = .false.  
      radmat%removal_set                       = .false.  
      radmat%transport_set                     = .false.  
      radmat%diffusion_set                     = .false.  
      radmat%fission_set                       = .false.  
      radmat%production_set                    = .false.  
      radmat%power_set                         = .false.  
      radmat%energy_released_per_fission_set   = .false.  
      radmat%particle_released_per_fission_set = .false.  
      radmat%prompt_spectrum_set               = .false.  
      radmat%velocity_set                      = .false.  
      radmat%beta_set                          = .false.                   
   
   end subroutine radmat_destroy   
   
   ! --------------------------------------------------------------------------

end module radiation_materials_destroy   
