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
!    License as published by the Free Software Foundation; either
!    version 2.1 of the License, or (at your option) any later version.
!
!    This library is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

subroutine test_radiation_materials_destroy
 
   !!< Test the procedures contained within the module radiation_materials_destroy in Radiation_Materials_Destroy.F90

   use futils   
   use spud
   use unittest_tools 
    
   use radiation_materials_data_types
   use radiation_materials_destroy
   use radiation_materials_create   
   use create_unittest_input
       
   implicit none

   ! local variables
   logical :: has_failed,has_warned
   type(np_radmat_type) :: np_radmat_1
   character(len=10000) :: error_message    
     
   ! none of these tests use warnings
   has_warned = .false.

   ! create hard coded np_radmat_input1
   call create_np_radmat_input1(np_radmat_1)
   
   ! test the destruction of it
   call test_destroy(np_radmat_1)
   
   call report_test("[test_radiation_materials_destroy for np_radmat_input1", &
                    has_failed, &
                    has_warned, &
                    "failed with error message "//trim(error_message))   
    
contains

   ! --------------------------------------------------------------------------
   
   subroutine test_destroy(np_radmat)
      
      !!< Test the destruction of the np_radmat passed, as well as one of each component type within it
      
      type(np_radmat_type), intent(inout) :: np_radmat
            
      ! local variables
      type(dataset_radmat_type) :: dataset_radmat
      type(physical_radmat_type) :: physical_radmat
      type(radmat_type) :: radmat
      type(np_radmat_size_type) :: np_radmat_size
      type(delayed_lambda_spectrum_type) :: delayed_lambda_spectrum
      
      has_failed = .false.
      
      error_message = 'no error'
      
      ! create other data types using the size arrays of the above hard coded np_radmat
      call allocate(dataset_radmat, &
                    np_radmat%np_radmat_size, &
                    dmat = 1)

      call allocate(physical_radmat, &
                    np_radmat%np_radmat_size, &
                    dmat = 1, &
                    pmat = 1)

      call allocate(radmat, &
                    np_radmat%np_radmat_size%number_of_energy_groups, &
                    np_radmat%np_radmat_size%number_of_scatter_moments(1), &
                    np_radmat%np_radmat_size%number_of_delayed_groups, &
                    allocate_all = .true.)
   
      np_radmat_size%total_number_dataset_radmats  = np_radmat%np_radmat_size%total_number_dataset_radmats
      np_radmat_size%total_number_physical_radmats = np_radmat%np_radmat_size%total_number_physical_radmats
      call allocate(np_radmat_size)

      call allocate(delayed_lambda_spectrum,& 
                   np_radmat%np_radmat_size%number_of_delayed_groups, &
                   np_radmat%np_radmat_size%number_of_energy_groups)
   
      ! set the data types values from the hard coded np_radmat components
      dataset_radmat          = np_radmat%dataset_radmats(1)
      physical_radmat         = np_radmat%dataset_radmats(1)%physical_radmats(1)
      radmat                  = np_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)
      np_radmat_size          = np_radmat%np_radmat_size
      delayed_lambda_spectrum = np_radmat%delayed_lambda_spectrum
   
      ! test the desutruction of each data type
      ! NOTE there is a order such as to check component destruction type first  
      !      and this procedure will return on the first fail      
      call test_radmat_destroy(radmat)
      
      if (has_failed) return

      call test_physical_radmat_destroy(physical_radmat)

      if (has_failed) return

      call test_dataset_radmat_destroy(dataset_radmat)

      if (has_failed) return
      
      call test_np_radmat_size_destroy(np_radmat_size)

      if (has_failed) return
      
      call test_delayed_lambda_spectrum_destroy(delayed_lambda_spectrum)

      if (has_failed) return

      call test_np_radmat_destroy(np_radmat)

   end subroutine test_destroy
   
   ! --------------------------------------------------------------------------
   
   subroutine test_np_radmat_destroy(np_radmat)

      type(np_radmat_type), intent(inout) :: np_radmat
      
      call destroy(np_radmat)

      if (allocated(np_radmat%dataset_radmats)) then
         has_failed = .true.
         error_message = ' failed as allocated(np_radmat%dataset_radmats) true'
         return
      end if 

      if (np_radmat%option_path /= "/uninitialised_path/") then
         has_failed = .true.
         error_message = 'test_np_radmat_destroy failed as np_radmat%option_path /= "/uninitialised_path/"'
         return
      end if 

      if (np_radmat%created) then
         has_failed = .true.
         error_message = 'test_np_radmat_destroy failed as np_radmat%created is true'
         return
      end if 

      if (np_radmat%readin) then
         has_failed = .true.
         error_message = 'test_np_radmat_destroy failed as np_radmat%readin is true'
         return
      end if 
      
      call test_np_radmat_size_destroy(np_radmat%np_radmat_size)
      if (has_failed) then
         error_message = 'test_np_radmat_destroy failed as '//trim(error_message)
         return
      end if

      call test_delayed_lambda_spectrum_destroy(np_radmat%delayed_lambda_spectrum)
      if (has_failed) then
         error_message = 'test_np_radmat_destroy failed as '//trim(error_message)
         return
      end if
      
   end subroutine test_np_radmat_destroy
   
   ! --------------------------------------------------------------------------
   
   subroutine test_dataset_radmat_destroy(dataset_radmat)

      type(dataset_radmat_type), intent(inout) :: dataset_radmat
      
      call destroy(dataset_radmat)

      if (allocated(dataset_radmat%physical_radmats)) then
         has_failed = .true.
         error_message = 'test_dataset_radmat_destroy failed as allocated(dataset_radmat%physical_radmats) true'
         return
      end if 

      if (dataset_radmat%option_path /= "/uninitialised_path/") then
         has_failed = .true.
         error_message = 'test_dataset_radmat_destroy failed as dataset_radmat%option_path /= "/uninitialised_path/"'
         return
      end if 
      
   end subroutine test_dataset_radmat_destroy
   
   ! --------------------------------------------------------------------------
   
   subroutine test_physical_radmat_destroy(physical_radmat)

      type(physical_radmat_type), intent(inout) :: physical_radmat
      
      call destroy(physical_radmat)
      
      if (allocated(physical_radmat%radmats)) then
         has_failed = .true.
         error_message = 'test_physical_radmat_destroy failed as allocated(physical_radmat%radmats) true'
         return
      end if 

      if (physical_radmat%option_path /= "/uninitialised_path/") then
         has_failed = .true.
         error_message = 'test_physical_radmat_destroy failed as physical_radmat%option_path /= "/uninitialised_path/"'
         return
      end if 
      
   end subroutine test_physical_radmat_destroy
   
   ! --------------------------------------------------------------------------
   
   subroutine test_radmat_destroy(radmat)

      type(radmat_type), intent(inout) :: radmat
      
      call destroy(radmat)
      
      if (allocated(radmat%total)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%total) true'
         return
      end if 
      
      if (allocated(radmat%absorption)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%absorption) true'
         return
      end if 
      
      if (allocated(radmat%scatter)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%scatter) true'
         return
      end if 
      
      if (allocated(radmat%removal)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%removal) true'
         return
      end if 
      
      if (allocated(radmat%transport)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%transport) true'
         return
      end if 
      
      if (allocated(radmat%diffusion)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%diffusion) true'
         return
      end if 
      
      if (allocated(radmat%fission)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%fission) true'
         return
      end if 
      
      if (allocated(radmat%production)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%production) true'
         return
      end if 
      
      if (allocated(radmat%power)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%power) true'
         return
      end if 
      
      if (allocated(radmat%energy_released_per_fission)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%energy_released_per_fission) true'
         return
      end if 
      
      if (allocated(radmat%np_released_per_fission)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%np_released_per_fission) true'
         return
      end if 
      
      if (allocated(radmat%prompt_spectrum)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%prompt_spectrum) true'
         return
      end if 
      
      if (allocated(radmat%velocity)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%velocity) true'
         return
      end if 
      
      if (allocated(radmat%beta)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as allocated(radmat%beta) true'
         return
      end if 
      
      if (radmat%total_set) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%total_set true'
         return
      end if 
      
      if (radmat%absorption_set) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%absorption_set true'
         return
      end if 
      
      if (radmat%scatter_set) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%scatter_set true'
         return
      end if 
      
      if (radmat%removal_set) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%removal_set true'
         return
      end if 
      
      if (radmat%transport_set(1)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%transport_set(1) true'
         return
      end if 
      
      if (radmat%transport_set(2)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%transport_set(2) true'
         return
      end if 
      
      if (radmat%transport_set(3)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%transport_set(3) true'
         return
      end if 
      
      if (radmat%diffusion_set(1)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%diffusion_set(1) true'
         return
      end if 
      
      if (radmat%diffusion_set(2)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%diffusion_set(2) true'
         return
      end if 
      
      if (radmat%diffusion_set(3)) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%diffusion_set(3) true'
         return
      end if 
      
      if (radmat%fission_set) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%fission_set true'
         return
      end if 
      
      if (radmat%production_set) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%production_set true'
         return
      end if 
      
      if (radmat%power_set) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%power_set true'
         return
      end if 
      
      if (radmat%energy_released_per_fission_set) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%energy_released_per_fission_set true'
         return
      end if 
      
      if (radmat%np_released_per_fission_set) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%np_released_per_fission_set true'
         return
      end if 
      
      if (radmat%prompt_spectrum_set) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%prompt_spectrum_set true'
         return
      end if 
      
      if (radmat%velocity_set) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%velocity_set true'
         return
      end if 
      
      if (radmat%beta_set) then
         has_failed = .true.
         error_message = 'test_radmat_destroy failed as radmat%beta_set true'
         return
      end if       

   end subroutine test_radmat_destroy
   
   ! --------------------------------------------------------------------------
   
   subroutine test_np_radmat_size_destroy(np_radmat_size)

      type(np_radmat_size_type), intent(inout) :: np_radmat_size
      
      call destroy(np_radmat_size)
      
      if (np_radmat_size%size_set) then
         has_failed = .true.
         error_message = 'test_np_radmat_size_destroy failed as np_radmat_size%size_set true'
         return
      end if       

      if (allocated(np_radmat_size%number_of_physical_radmats)) then
         has_failed = .true.
         error_message = 'test_np_radmat_size_destroy failed as allocated(np_radmat_size%number_of_physical_radmats) true'
         return
      end if 

      if (allocated(np_radmat_size%number_of_radmats)) then
         has_failed = .true.
         error_message = 'test_np_radmat_size_destroy failed as allocated(np_radmat_size%number_of_radmats) true'
         return
      end if 

      if (allocated(np_radmat_size%number_of_scatter_moments)) then
         has_failed = .true.
         error_message = 'test_np_radmat_size_destroy failed as allocated(np_radmat_size%number_of_scatter_moments) true'
         return
      end if 

      if (np_radmat_size%number_of_energy_groups /= 0) then
         has_failed = .true.
         error_message = 'test_np_radmat_size_destroy failed as np_radmat_size%number_of_energy_groups /= 0'
         return
      end if 

      if (np_radmat_size%number_of_delayed_groups /= 0) then
         has_failed = .true.
         error_message = 'test_np_radmat_size_destroy failed as np_radmat_size%number_of_delayed_groups /= 0'
         return
      end if 

      if (np_radmat_size%total_number_radmats /= 0) then
         has_failed = .true.
         error_message = 'test_np_radmat_size_destroy failed as np_radmat_size%total_number_radmats /= 0'
         return
      end if 

      if (np_radmat_size%total_number_physical_radmats /= 0) then
         has_failed = .true.
         error_message = 'test_np_radmat_size_destroy failed as np_radmat_size%total_number_physical_radmats /= 0'
         return
      end if 

      if (np_radmat_size%total_number_dataset_radmats /= 0) then
         has_failed = .true.
         error_message = 'test_np_radmat_size_destroy failed as np_radmat_size%total_number_dataset_radmats /= 0'
         return
      end if 
      
   end subroutine test_np_radmat_size_destroy
   
   ! --------------------------------------------------------------------------
   
   subroutine test_delayed_lambda_spectrum_destroy(delayed_lambda_spectrum)

      type(delayed_lambda_spectrum_type) :: delayed_lambda_spectrum
      
      call destroy(delayed_lambda_spectrum)
      
      if (allocated(delayed_lambda_spectrum%lambda)) then
         has_failed = .true.
         error_message = 'test_delayed_lambda_spectrum_destroy failed as allocated(delayed_lambda_spectrum%lambda) true'
         return
      end if 
      
      if (allocated(delayed_lambda_spectrum%spectrum)) then
         has_failed = .true.
         error_message = 'test_delayed_lambda_spectrum_destroy failed as allocated(delayed_lambda_spectrum%spectrum) true'
         return
      end if 
      
      if (delayed_lambda_spectrum%lambda_set) then
         has_failed = .true.
         error_message = 'test_delayed_lambda_spectrum_destroy failed as delayed_lambda_spectrum%lambda_set true'
         return
      end if       
      
      if (delayed_lambda_spectrum%spectrum_set) then
         has_failed = .true.
         error_message = 'test_delayed_lambda_spectrum_destroy failed as delayed_lambda_spectrum%spectrum_set true'
         return
      end if       
      
   end subroutine test_delayed_lambda_spectrum_destroy
   
   ! --------------------------------------------------------------------------

end subroutine test_radiation_materials_destroy
