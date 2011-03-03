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

module radiation_materials_operations
   
   !!< module containing procedures that perform operations with the np_radmat data type
   
   use futils
   use unittest_tools

   use radiation_materials_data_types 
   
   implicit none
   
   private 
             
   public :: radmat_type_equals

   interface radmat_type_equals
      module procedure radmat_equals, &
                       physical_radmat_equals, &
                       dataset_radmat_equals, &
                       delayed_lambda_spectrum_equals, &
                       np_radmat_size_equals, &
                       np_radmat_equals
   end interface 
                  
contains

   ! --------------------------------------------------------------------------

   function np_radmat_equals(np_radmat_1, &
                             np_radmat_2, &
                             error_message_np_radmat_equals) result (equals)

      !!< Compare two np_radmat_type data types to see if they are equal
      !!< any real components are checked to a default tolerance  
      !!< Return an error message of what failed
      
      type(np_radmat_type), intent(in) :: np_radmat_1
      type(np_radmat_type), intent(in) :: np_radmat_2
      character(len=*), intent(inout) :: error_message_np_radmat_equals
      
      logical :: equals
      
      ! local variable
      integer :: dmat
      character(len=10000) :: error_message_np_radmat_size_equals
      character(len=10000) :: error_message_delayed_lambda_spectrum_equals
      character(len=10000) :: error_message_dataset_radmat_equals
            
      equals = .true.
      
      error_message_np_radmat_size_equals = 'no error'

      ! check the option path
      if (np_radmat_1%option_path /= np_radmat_2%option_path) then
         equals = .false.
         error_message_np_radmat_equals = 'np_radmat_1%option_path /= np_radmat_2%option_path'
         return
      end if 
      
      ! check the flags
      if (np_radmat_1%created .neqv. np_radmat_2%created) then
         equals = .false.
         error_message_np_radmat_equals = 'np_radmat_1%created .neqv. np_radmat_2%created'
         return 
      end if 
      if (np_radmat_1%readin .neqv. np_radmat_2%readin) then
         equals = .false.
         error_message_np_radmat_equals = 'np_radmat_1%readin .neqv. np_radmat_2%readin'
         return 
      end if 
      
      created_if: if (np_radmat_1%created) then 
      
         ! check the sizes
         equals = radmat_type_equals(np_radmat_1%np_radmat_size,np_radmat_2%np_radmat_size,error_message_np_radmat_size_equals)
         if (.not. equals) then
            error_message_np_radmat_equals = 'np_radmat_1%np_radmat_size /= np_radmat_2%np_radmat_size, with error message '&
                                           &//trim(error_message_np_radmat_size_equals)
            return
         end if
      
         ! check the delayed lambda and spectrum
         equals = radmat_type_equals(np_radmat_1%delayed_lambda_spectrum,np_radmat_2%delayed_lambda_spectrum,error_message_delayed_lambda_spectrum_equals)
         if (.not. equals) then 
            error_message_np_radmat_equals = 'np_radmat_1%delayed_lambda_spectrum /= np_radmat_2%delayed_lambda_spectrum, with error message '&
                                           &//trim(error_message_delayed_lambda_spectrum_equals)
            return
         end if
      
         ! check the allocated status of the datasets
         if (allocated(np_radmat_1%dataset_radmats) .neqv. allocated(np_radmat_2%dataset_radmats)) then
            equals = .false.
            error_message_np_radmat_equals = 'allocated(np_radmat_1%dataset_radmats) .neqv. allocated(np_radmat_2%dataset_radmats)'
            return
         end if 

         allocated_if: if (allocated(np_radmat_1%dataset_radmats)) then
      
            ! check have same size of dataset_radmats
            if (size(np_radmat_1%dataset_radmats) /= size(np_radmat_2%dataset_radmats)) then
               equals = .false.
               error_message_np_radmat_equals = 'size(np_radmat_1%dataset_radmats) /= size(np_radmat_2%dataset_radmats)'
               return
            end if 
      
            ! now check each dataset_radmat type associated with the np_radmat type
            dmat_loop: do dmat = 1,size(np_radmat_1%dataset_radmats)
      
               equals = radmat_type_equals(np_radmat_1%dataset_radmats(dmat),np_radmat_2%dataset_radmats(dmat),error_message_dataset_radmat_equals)
           
               if (.not. equals) then
                  error_message_np_radmat_equals = 'np_radmat_1%dataset_radmats('//int2str(dmat)//') /= np_radmat_2%dataset_radmats('//int2str(dmat)//'), with error message '&
                                                 &//trim(error_message_dataset_radmat_equals)
                  return
               end if
      
            end do dmat_loop
      
         end if allocated_if
      
      end if created_if
      
   end function np_radmat_equals

   ! --------------------------------------------------------------------------

   function np_radmat_size_equals(np_radmat_size_1, &
                                  np_radmat_size_2, &
                                  error_message_np_radmat_size_equals) result (equals)

      !!< Compare two np_radmat_size_type data types to see if they are equal
      !!< any real components are checked to a default tolerance  
      !!< Return an error message of what failed
      
      type(np_radmat_size_type), intent(in) :: np_radmat_size_1
      type(np_radmat_size_type), intent(in) :: np_radmat_size_2
      character(len=*), intent(inout) :: error_message_np_radmat_size_equals
      
      logical :: equals
      
      ! local variable
      integer :: i
      
      equals = .true.

      error_message_np_radmat_size_equals = 'no error'
      
      ! check the set flags
      if (np_radmat_size_1%size_set .neqv. np_radmat_size_2%size_set) then  
         equals = .false.
         error_message_np_radmat_size_equals = 'np_radmat_size_1%size_set .neqv. np_radmat_size_2%size_se'
         return 
      end if
      
      ! check the integers
      if (np_radmat_size_1%total_number_radmats /= np_radmat_size_2%total_number_radmats) then           
         equals = .false.
         error_message_np_radmat_size_equals = 'np_radmat_size_1%total_number_radmats /= np_radmat_size_2%total_number_radmats'
         return 
      end if
      if (np_radmat_size_1%total_number_physical_radmats /= np_radmat_size_2%total_number_physical_radmats) then  
         equals = .false.
         error_message_np_radmat_size_equals = 'np_radmat_size_1%total_number_physical_radmats /= np_radmat_size_2%total_number_physical_radmats'
         return 
      end if
      if (np_radmat_size_1%total_number_dataset_radmats /= np_radmat_size_2%total_number_dataset_radmats) then   
         equals = .false.
         error_message_np_radmat_size_equals = 'np_radmat_size_1%total_number_dataset_radmats /= np_radmat_size_2%total_number_dataset_radmats'
         return 
      end if
      if (np_radmat_size_1%number_of_energy_groups /= np_radmat_size_2%number_of_energy_groups) then        
         equals = .false.
         error_message_np_radmat_size_equals = 'np_radmat_size_1%number_of_energy_groups /= np_radmat_size_2%number_of_energy_groups'
         return 
      end if
      if (np_radmat_size_1%number_of_delayed_groups /= np_radmat_size_2%number_of_delayed_groups) then       
         equals = .false.
         error_message_np_radmat_size_equals = 'np_radmat_size_1%number_of_delayed_groups /= np_radmat_size_2%number_of_delayed_groups'
         return 
      end if
            
      ! check that have same allocated status 
      if (allocated(np_radmat_size_1%number_of_physical_radmats) .neqv. allocated(np_radmat_size_2%number_of_physical_radmats)) then  
         equals = .false.
         error_message_np_radmat_size_equals = 'allocated(np_radmat_size_1%number_of_physical_radmats) .neqv. allocated(np_radmat_size_2%number_of_physical_radmats'
         return
      end if
      if (allocated(np_radmat_size_1%number_of_radmats) .neqv. allocated(np_radmat_size_2%number_of_radmats)) then           
         equals = .false.
         error_message_np_radmat_size_equals = 'allocated(np_radmat_size_1%number_of_radmats) .neqv. allocated(np_radmat_size_2%number_of_radmats)'
         return
      end if
      if (allocated(np_radmat_size_1%number_of_radmats_base) .neqv. allocated(np_radmat_size_2%number_of_radmats_base)) then      
         equals = .false.
         error_message_np_radmat_size_equals = 'allocated(np_radmat_size_1%number_of_radmats_base) .neqv. allocated(np_radmat_size_2%number_of_radmats_base)'
         return
      end if
      if (allocated(np_radmat_size_1%number_of_scatter_moments)  .neqv. allocated(np_radmat_size_2%number_of_scatter_moments)) then   
         equals = .false.
         error_message_np_radmat_size_equals = 'allocated(np_radmat_size_1%number_of_scatter_moments)  .neqv. allocated(np_radmat_size_2%number_of_scatter_moments)'
         return
      end if
     
      ! check the size
      if (allocated(np_radmat_size_1%number_of_physical_radmats) .and. &
         &(size(np_radmat_size_1%number_of_physical_radmats) /= size(np_radmat_size_2%number_of_physical_radmats))) then  
         equals = .false.
         error_message_np_radmat_size_equals = 'size(np_radmat_size_1%number_of_physical_radmats) /= size(np_radmat_size_2%number_of_physical_radmats)'
         return 
      end if   

      if (allocated(np_radmat_size_1%number_of_radmats) .and. &
         &(size(np_radmat_size_1%number_of_radmats) /= size(np_radmat_size_2%number_of_radmats))) then  
         equals = .false.
         error_message_np_radmat_size_equals = 'size(np_radmat_size_1%number_of_radmats) /= size(np_radmat_size_2%number_of_radmats)'
         return 
      end if  

      if (allocated(np_radmat_size_1%number_of_radmats_base) .and. &
         &(size(np_radmat_size_1%number_of_radmats_base) /= size(np_radmat_size_2%number_of_radmats_base))) then  
         equals = .false.
         error_message_np_radmat_size_equals = 'size(np_radmat_size_1%number_of_radmats_base) /= size(np_radmat_size_2%number_of_radmats_base)'
         return 
      end if   

      if (allocated(np_radmat_size_1%number_of_scatter_moments)   .and. &
         &(size(np_radmat_size_1%number_of_scatter_moments) /= size(np_radmat_size_2%number_of_scatter_moments))) then  
         equals = .false.
         error_message_np_radmat_size_equals = 'size(np_radmat_size_1%number_of_scatter_moments) /= size(np_radmat_size_2%number_of_scatter_moments)'
         return 
      end if   
         
      ! check the values
      if (allocated(np_radmat_size_1%number_of_physical_radmats) .and. np_radmat_size_1%size_set) then
         loop1: do i = 1,size(np_radmat_size_1%number_of_physical_radmats)
            if (np_radmat_size_1%number_of_physical_radmats(i) /= np_radmat_size_2%number_of_physical_radmats(i))  then 
               equals = .false.
               error_message_np_radmat_size_equals = 'np_radmat_size_1%number_of_physical_radmats('//int2str(i)//') /= np_radmat_size_2%number_of_physical_radmats('//int2str(i)//')'
               return 
            end if
         end do loop1
      end if

      if (allocated(np_radmat_size_1%number_of_radmats) .and. np_radmat_size_1%size_set) then      
         loop2: do i = 1,size(np_radmat_size_1%number_of_radmats)      
            if (np_radmat_size_1%number_of_radmats(i) /= np_radmat_size_2%number_of_radmats(i)) then  
               equals = .false.
               error_message_np_radmat_size_equals = 'np_radmat_size_1%number_of_radmats('//int2str(i)//') /= np_radmat_size_2%number_of_radmats('//int2str(i)//')'
               return 
            end if
         end do loop2
      end if

      if (allocated(np_radmat_size_1%number_of_radmats_base) .and. np_radmat_size_1%size_set) then      
         loop3: do i = 1,size(np_radmat_size_1%number_of_radmats_base)
            if (np_radmat_size_1%number_of_radmats_base(i) /= np_radmat_size_2%number_of_radmats_base(i)) then  
               equals = .false.
               error_message_np_radmat_size_equals = 'np_radmat_size_1%number_of_radmats_base('//int2str(i)//') /= np_radmat_size_2%number_of_radmats_base('//int2str(i)//')'
               return 
            end if
         end do loop3
      end if
      
      if (allocated(np_radmat_size_1%number_of_scatter_moments) .and. np_radmat_size_1%size_set ) then
         loop4: do i = 1,size(np_radmat_size_1%number_of_scatter_moments)
            if (np_radmat_size_1%number_of_scatter_moments(i) /= np_radmat_size_2%number_of_scatter_moments(i)) then  
               equals = .false.
               error_message_np_radmat_size_equals = 'np_radmat_size_1%number_of_scatter_moments('//int2str(i)//') /= np_radmat_size_2%number_of_scatter_moments('//int2str(i)//')'
               return 
            end if
         end do loop4
      end if
                  
   end function np_radmat_size_equals

   ! --------------------------------------------------------------------------

   function delayed_lambda_spectrum_equals(delayed_lambda_spectrum_1, &
                                           delayed_lambda_spectrum_2, &
                                           error_message_delayed_lambda_spectrum_equals) result (equals)

      !!< Compare two delayed_lambda_spectrum_type data types to see if they are equal
      !!< any real components are checked to a default tolerance  
      !!< Return an error message of what failed
      
      type(delayed_lambda_spectrum_type), intent(in) :: delayed_lambda_spectrum_1
      type(delayed_lambda_spectrum_type), intent(in) :: delayed_lambda_spectrum_2
      character(len=*), intent(inout) :: error_message_delayed_lambda_spectrum_equals
      
      logical :: equals
      
      equals = .true.

      error_message_delayed_lambda_spectrum_equals = 'no error'
      
      ! check the set flags
      if (delayed_lambda_spectrum_1%lambda_set .neqv. delayed_lambda_spectrum_2%lambda_set) then    
         equals = .false.
         error_message_delayed_lambda_spectrum_equals = 'delayed_lambda_spectrum_1%lambda_set .neqv. delayed_lambda_spectrum_2%lambda_set'
         return 
      end if       
      if (delayed_lambda_spectrum_1%spectrum_set .neqv. delayed_lambda_spectrum_2%spectrum_set) then  
         equals = .false.
         error_message_delayed_lambda_spectrum_equals = 'delayed_lambda_spectrum_1%spectrum_set .neqv. delayed_lambda_spectrum_2%spectrum_set'
         return 
      end if       
      
      ! check that have same allocated status 
      if (allocated(delayed_lambda_spectrum_1%lambda) .neqv. allocated(delayed_lambda_spectrum_2%lambda)) then    
         equals = .false.
         error_message_delayed_lambda_spectrum_equals = 'allocated(delayed_lambda_spectrum_1%lambda) .neqv. allocated(delayed_lambda_spectrum_2%lambda)'
         return
      end if       
      if (allocated(delayed_lambda_spectrum_1%spectrum) .neqv. allocated(delayed_lambda_spectrum_2%spectrum)) then  
         equals = .false.
         error_message_delayed_lambda_spectrum_equals = 'allocated(delayed_lambda_spectrum_1%spectrum) .neqv. allocated(delayed_lambda_spectrum_2%spectrum)'
         return
      end if       
      
      ! check the size
      if (allocated(delayed_lambda_spectrum_1%lambda) .and. &
         &(size(delayed_lambda_spectrum_1%lambda) /= size(delayed_lambda_spectrum_2%lambda))) then  
         equals = .false.  
         error_message_delayed_lambda_spectrum_equals = 'size(delayed_lambda_spectrum_1%lambda) /= size(delayed_lambda_spectrum_2%lambda)'       
         return
      end if       
         
      if (allocated(delayed_lambda_spectrum_1%spectrum) .and. &
         &(size(delayed_lambda_spectrum_1%spectrum) /= size(delayed_lambda_spectrum_2%spectrum))) then  
         equals = .false.
         error_message_delayed_lambda_spectrum_equals = 'size(delayed_lambda_spectrum_1%spectrum) /= size(delayed_lambda_spectrum_2%spectrum)'
         return 
      end if       
      
      ! check the values
      if (allocated(delayed_lambda_spectrum_1%lambda) .and. delayed_lambda_spectrum_1%lambda_set .and. &
         &(.not. fequals(delayed_lambda_spectrum_1%lambda,delayed_lambda_spectrum_2%lambda))) then  
         equals = .false.
         error_message_delayed_lambda_spectrum_equals = 'delayed_lambda_spectrum_1%lambda /= delayed_lambda_spectrum_2%lambda'
         return 
      end if       
         
      if (allocated(delayed_lambda_spectrum_1%spectrum) .and. delayed_lambda_spectrum_1%spectrum_set .and. &
         &(.not. fequals(delayed_lambda_spectrum_1%spectrum,delayed_lambda_spectrum_2%spectrum))) then  
         equals = .false.
         error_message_delayed_lambda_spectrum_equals = 'delayed_lambda_spectrum_1%spectrum /= delayed_lambda_spectrum_2%spectrum'
         return
      end if 
            
   end function delayed_lambda_spectrum_equals

   ! --------------------------------------------------------------------------

   function dataset_radmat_equals(dataset_radmat_1, &
                                  dataset_radmat_2, &
                                  error_message_dataset_radmat_equals) result (equals)

      !!< Compare two dataset_radmat_type data types to see if they are equal
      !!< any real components are checked to a default tolerance  
      !!< Return an error message of what failed
      
      type(dataset_radmat_type), intent(in) :: dataset_radmat_1
      type(dataset_radmat_type), intent(in) :: dataset_radmat_2
      character(len=*), intent(inout) :: error_message_dataset_radmat_equals
      
      logical :: equals

      ! local variables
      integer :: pmat
      character(len=10000) :: error_message_physical_radmat_equals
      
      equals = .true.

      error_message_dataset_radmat_equals = 'no error'
            
      ! check the option path
      if (dataset_radmat_1%option_path /= dataset_radmat_2%option_path) then  
         equals = .false.
         error_message_dataset_radmat_equals = 'dataset_radmat_1%option_path /= dataset_radmat_2%option_path'
         return
      end if 
      
      ! check that have same allocated status of physical_radmats component
      if (allocated(dataset_radmat_1%physical_radmats) .neqv. allocated(dataset_radmat_2%physical_radmats)) then  
         equals = .false.
         error_message_dataset_radmat_equals = 'allocated(dataset_radmat_1%physical_radmats) .neqv. allocated(dataset_radmat_2%physical_radmats)'
         return
      end if 
      
      allocated_if: if (allocated(dataset_radmat_1%physical_radmats)) then
      
         ! check have same size of physical_radmats
         if (size(dataset_radmat_1%physical_radmats) /= size(dataset_radmat_2%physical_radmats)) then  
            equals = .false.
            error_message_dataset_radmat_equals = 'size(dataset_radmat_1%physical_radmats) /= size(dataset_radmat_2%physical_radmats)'
            return
         end if 
      
         ! now check each physical_radmat type associated with the dataset_radmat type
         pmat_loop: do pmat = 1,size(dataset_radmat_1%physical_radmats)
      
            equals = radmat_type_equals(dataset_radmat_1%physical_radmats(pmat),dataset_radmat_2%physical_radmats(pmat),error_message_physical_radmat_equals)
        
            if (.not. equals) then
               error_message_dataset_radmat_equals = 'dataset_radmat_1%physical_radmats('//int2str(pmat)//') /= dataset_radmat_2%physical_radmats('//int2str(pmat)//'), with error message '&
                                                   &//trim(error_message_physical_radmat_equals)
               return
            end if 
      
         end do pmat_loop
      
      end if allocated_if

   end function dataset_radmat_equals

   ! --------------------------------------------------------------------------

   function physical_radmat_equals(physical_radmat_1, &
                                   physical_radmat_2, &
                                   error_message_physical_radmat_equals) result (equals)

      !!< Compare two physical_radmat_type data types to see if they are equal
      !!< any real components are checked to a default tolerance  
      !!< Return an error message of what failed
      
      type(physical_radmat_type), intent(in) :: physical_radmat_1
      type(physical_radmat_type), intent(in) :: physical_radmat_2
      character(len=*), intent(inout) :: error_message_physical_radmat_equals
      
      logical :: equals

      ! local variables
      integer :: rmat
      character(len=10000) :: error_message_radmat_equals
      
      equals = .true.
      
      error_message_physical_radmat_equals = 'no error'
      
      ! check the option path
      if (physical_radmat_1%option_path /= physical_radmat_2%option_path) then 
         equals = .false.
         error_message_physical_radmat_equals = 'physical_radmat_1%option_path /= physical_radmat_2%option_path'
         return
      end if

      ! check that have same allocated status of physical_radmats component
      if (allocated(physical_radmat_1%radmats) .neqv. allocated(physical_radmat_2%radmats)) then 
         equals = .false.
         error_message_physical_radmat_equals = 'allocated(physical_radmat_1%radmats) .neqv. allocated(physical_radmat_2%radmats)'
         return
      end if 

      allocated_if: if (allocated(physical_radmat_1%radmats)) then
      
         ! check have same size of radmats
         if (size(physical_radmat_1%radmats) /= size(physical_radmat_2%radmats)) then
            equals = .false.
            error_message_physical_radmat_equals = 'size(physical_radmat_1%radmats) /= size(physical_radmat_2%radmats)'
            return
         end if 
      
         ! now check each radmat type associated with the physical_radmat type
         rmat_loop: do rmat = 1,size(physical_radmat_1%radmats)
      
            equals = radmat_type_equals(physical_radmat_1%radmats(rmat),physical_radmat_2%radmats(rmat),error_message_radmat_equals)

            if (.not. equals) then
               error_message_physical_radmat_equals = 'physical_radmat_1%radmats('//int2str(rmat)//'),physical_radmat_2%radmats('//int2str(rmat)//'), with error message '&
                                                    &//trim(error_message_radmat_equals)
               return
            end if 
      
         end do rmat_loop
      
      end if allocated_if

   end function physical_radmat_equals

   ! --------------------------------------------------------------------------
    
   function radmat_equals(radmat_1, &
                          radmat_2, &
                          error_message_radmat_equals) result (equals)
      
      !!< Compare two radmat_type data types to see if they are equal
      !!< any real components are checked to a default tolerance  
      !!< Return an error message of what failed
      
      type(radmat_type), intent(in) :: radmat_1
      type(radmat_type), intent(in) :: radmat_2
      character(len=*), intent(inout) :: error_message_radmat_equals
      
      logical :: equals
      
      ! local variables
      integer :: m
      
      equals = .true.
      
      error_message_radmat_equals = 'no error'

      ! check the set flags
      if (radmat_1%total_set .neqv. radmat_2%total_set) then
         equals = .false.
         error_message_radmat_equals = 'radmat_1%total_set .neqv. radmat_2%total_set'
         return
      end if       
      if (radmat_1%absorption_set .neqv. radmat_2%absorption_set) then      
         equals = .false.
         error_message_radmat_equals = 'radmat_1%absorption_set .neqv. radmat_2%absorption_set'
         return
      end if       
      if (radmat_1%scatter_set .neqv. radmat_2%scatter_set) then         
         equals = .false.
         error_message_radmat_equals = 'radmat_1%scatter_set .neqv. radmat_2%scatter_set'
         return
      end if       
      if (radmat_1%removal_set .neqv. radmat_2%removal_set) then         
         equals = .false.
         error_message_radmat_equals = 'radmat_1%removal_set .neqv. radmat_2%removal_set'
         return
      end if       
      if (radmat_1%transport_set(1) .neqv. radmat_2%transport_set(1)) then    
         equals = .false.
         error_message_radmat_equals = 'radmat_1%transport_set(1) .neqv. radmat_2%transport_set(1)'
         return
      end if       
      if (radmat_1%transport_set(2) .neqv. radmat_2%transport_set(2)) then    
         equals = .false.
         error_message_radmat_equals = 'radmat_1%transport_set(2) .neqv. radmat_2%transport_set(2)'
         return
      end if       
      if (radmat_1%transport_set(3) .neqv. radmat_2%transport_set(3)) then    
         equals = .false.
         error_message_radmat_equals = 'radmat_1%transport_set(3) .neqv. radmat_2%transport_set(3)'
         return
      end if       
      if (radmat_1%diffusion_set(1) .neqv. radmat_2%diffusion_set(1)) then    
         equals = .false.
         error_message_radmat_equals = 'radmat_1%diffusion_set(1) .neqv. radmat_2%diffusion_set(1)'
         return
      end if       
      if (radmat_1%diffusion_set(2) .neqv. radmat_2%diffusion_set(2)) then    
         equals = .false.
         error_message_radmat_equals = 'radmat_1%diffusion_set(2) .neqv. radmat_2%diffusion_set(2)'
         return
      end if       
      if (radmat_1%diffusion_set(3) .neqv. radmat_2%diffusion_set(3)) then    
         equals = .false.
         error_message_radmat_equals = 'radmat_1%diffusion_set(3) .neqv. radmat_2%diffusion_set(3)'
         return      
      end if       
      if (radmat_1%fission_set .neqv. radmat_2%fission_set) then         
         equals = .false.
         error_message_radmat_equals = 'radmat_1%fission_set .neqv. radmat_2%fission_set'
         return
      end if       
      if (radmat_1%production_set .neqv. radmat_2%production_set) then      
         equals = .false.
         error_message_radmat_equals = 'radmat_1%production_set .neqv. radmat_2%production_set'
         return
      end if       
      if (radmat_1%power_set .neqv. radmat_2%power_set) then           
         equals = .false.
         error_message_radmat_equals = 'radmat_1%power_set .neqv. radmat_2%power_set'
         return
      end if       
      if (radmat_1%energy_released_per_fission_set .neqv. radmat_2%energy_released_per_fission_set) then           
         equals = .false.
         error_message_radmat_equals = 'radmat_1%energy_released_per_fission_set .neqv. radmat_2%energy_released_per_fission_set'
         return
      end if       
      if (radmat_1% np_released_per_fission_set.neqv. radmat_2%np_released_per_fission_set) then           
         equals = .false.
         error_message_radmat_equals = 'radmat_1% np_released_per_fission_set.neqv. radmat_2%np_released_per_fission_set'
         return
      end if       
      if (radmat_1%prompt_spectrum_set .neqv. radmat_2%prompt_spectrum_set) then 
         equals = .false.
         error_message_radmat_equals = 'radmat_1%prompt_spectrum_set .neqv. radmat_2%prompt_spectrum_set'
         return
      end if       
      if (radmat_1%velocity_set .neqv. radmat_2%velocity_set) then        
         equals = .false.
         error_message_radmat_equals = 'radmat_1%velocity_set .neqv. radmat_2%velocity_set'
         return
      end if       
      if (radmat_1%beta_set .neqv. radmat_2%beta_set) then            
         equals = .false.
         error_message_radmat_equals = 'radmat_1%beta_set .neqv. radmat_2%beta_set'
         return
      end if 


            
      ! check the real arrays - allocated, size then values
      if (allocated(radmat_1%total) .neqv. allocated(radmat_2%total)) then           
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%total) .neqv. allocated(radmat_2%total)'
         return      
      end if       
      if (allocated(radmat_1%absorption) .neqv. allocated(radmat_2%absorption)) then      
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%absorption) .neqv. allocated(radmat_2%absorption)'
         return      
      end if       
      if (allocated(radmat_1%scatter) .neqv. allocated(radmat_2%scatter)) then         
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%scatter) .neqv. allocated(radmat_2%scatter)'
         return      
      end if       
      if (allocated(radmat_1%removal) .neqv. allocated(radmat_2%removal)) then         
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%removal) .neqv. allocated(radmat_2%removal)'
         return      
      end if       
      if (allocated(radmat_1%transport) .neqv. allocated(radmat_2%transport)) then       
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%transport) .neqv. allocated(radmat_2%transport)'
         return      
      end if       
      if (allocated(radmat_1%diffusion) .neqv. allocated(radmat_2%diffusion)) then       
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%diffusion) .neqv. allocated(radmat_2%diffusion)'
         return      
      end if       
      if (allocated(radmat_1%fission) .neqv. allocated(radmat_2%fission)) then         
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%fission) .neqv. allocated(radmat_2%fission)'
         return      
      end if       
      if (allocated(radmat_1%production) .neqv. allocated(radmat_2%production)) then      
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%production) .neqv. allocated(radmat_2%production)'
         return      
      end if       
      if (allocated(radmat_1%power) .neqv. allocated(radmat_2%power)) then           
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%power) .neqv. allocated(radmat_2%power)'
         return      
      end if       
      if (allocated(radmat_1%energy_released_per_fission) .neqv. allocated(radmat_2%energy_released_per_fission)) then           
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%energy_released_per_fission) .neqv. allocated(radmat_2%energy_released_per_fission)'
         return      
      end if       
      if (allocated(radmat_1%np_released_per_fission) .neqv. allocated(radmat_2%np_released_per_fission)) then           
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%np_released_per_fission) .neqv. allocated(radmat_2%np_released_per_fission)'
         return      
      end if       
      if (allocated(radmat_1%prompt_spectrum) .neqv. allocated(radmat_2%prompt_spectrum)) then 
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%prompt_spectrum) .neqv. allocated(radmat_2%prompt_spectrum)'
         return      
      end if       
      if (allocated(radmat_1%velocity) .neqv. allocated(radmat_2%velocity)) then        
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%velocity) .neqv. allocated(radmat_2%velocity)'
         return      
      end if 
      if (allocated(radmat_1%beta) .neqv. allocated(radmat_2%beta)) then        
         equals = .false.
         error_message_radmat_equals = 'allocated(radmat_1%beta) .neqv. allocated(radmat_2%beta)'
         return      
      end if 


 
      if (allocated(radmat_1%total) .and. (size(radmat_1%total) /= size(radmat_2%total))) then           
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%total) /= size(radmat_2%total)'
         return      
      end if       
      if (allocated(radmat_1%absorption) .and. (size(radmat_1%absorption) /= size(radmat_2%absorption))) then      
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%absorption) /= size(radmat_2%absorption)'
         return      
      end if       
      if (allocated(radmat_1%scatter) .and. (size(radmat_1%scatter) /= size(radmat_2%scatter))) then         
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%scatter) /= size(radmat_2%scatter)'
         return      
      end if       
      if (allocated(radmat_1%removal) .and. (size(radmat_1%removal) /= size(radmat_2%removal))) then         
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%removal) /= size(radmat_2%removal)'
         return      
      end if       
      if (allocated(radmat_1%transport) .and. (size(radmat_1%transport) /= size(radmat_2%transport))) then       
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%transport) /= size(radmat_2%transport)'
         return      
      end if       
      if (allocated(radmat_1%diffusion) .and. (size(radmat_1%diffusion) /= size(radmat_2%diffusion))) then       
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%diffusion) /= size(radmat_2%diffusion)'
         return      
      end if       
      if (allocated(radmat_1%fission) .and. (size(radmat_1%fission) /= size(radmat_2%fission))) then         
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%fission) /= size(radmat_2%fission)'
         return      
      end if       
      if (allocated(radmat_1%production) .and. (size(radmat_1%production) /= size(radmat_2%production))) then      
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%production) /= size(radmat_2%production)'
         return      
      end if       
      if (allocated(radmat_1%power) .and. (size(radmat_1%power) /= size(radmat_2%power))) then           
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%power) /= size(radmat_2%power)'
         return      
      end if       
      if (allocated(radmat_1%energy_released_per_fission) .and. (size(radmat_1%energy_released_per_fission) /= size(radmat_2%energy_released_per_fission))) then           
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%energy_released_per_fission) /= size(radmat_2%energy_released_per_fission)'
         return      
      end if       
      if (allocated(radmat_1%np_released_per_fission) .and. (size(radmat_1%np_released_per_fission) /= size(radmat_2%np_released_per_fission))) then           
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%np_released_per_fission) /= size(radmat_2%np_released_per_fission)'
         return      
      end if       
      if (allocated(radmat_1%prompt_spectrum) .and. (size(radmat_1%prompt_spectrum) /= size(radmat_2%prompt_spectrum))) then 
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%prompt_spectrum) /= size(radmat_2%prompt_spectrum)'
         return      
      end if       
      if (allocated(radmat_1%velocity) .and. (size(radmat_1%velocity) /= size(radmat_2%velocity))) then        
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%velocity) /= size(radmat_2%velocity)'
         return      
      end if             
      if (allocated(radmat_1%beta) .and. (size(radmat_1%beta) /= size(radmat_2%beta))) then        
         equals = .false.
         error_message_radmat_equals = 'size(radmat_1%beta) /= size(radmat_2%beta)'
         return      
      end if             



      if (allocated(radmat_1%total) .and. radmat_1%total_set .and. (.not. fequals(radmat_1%total,radmat_2%total))) then           
         equals = .false.
         error_message_radmat_equals = 'radmat_1%total /= radmat_2%tota'
         return
      end if       
      if (allocated(radmat_1%absorption) .and. radmat_1%absorption_set .and. (.not. fequals(radmat_1%absorption,radmat_2%absorption))) then      
         equals = .false.
         error_message_radmat_equals = 'radmat_1%absorption /= radmat_2%absorption'
         return      
      end if             
      mom_loop: do m = 1,size(radmat_1%scatter,3)
      if (allocated(radmat_1%scatter) .and. radmat_1%scatter_set .and. (.not. fequals(radmat_1%scatter(:,:,m), radmat_2%scatter(:,:,m)))) then   
         equals = .false.
         error_message_radmat_equals = 'radmat_1%scatter(:,:,'//int2str(m)//') /= radmat_2%scatter(:,:,'//int2str(m)//')'
         return            
      end if       
      end do mom_loop          
      if (allocated(radmat_1%removal) .and. radmat_1%removal_set .and. (.not. fequals(radmat_1%removal, radmat_2%removal))) then          
         equals = .false.
         error_message_radmat_equals = 'radmat_1%removal /= radmat_2%removal'
         return
      end if       
      if (allocated(radmat_1%transport) .and. radmat_1%transport_set(1) .and. (.not. fequals(radmat_1%transport(:,1), radmat_2%transport(:,1)))) then   
         equals = .false.
         error_message_radmat_equals = 'radmat_1%transport(:,1) /= radmat_2%transport(:,1)'
         return
      end if       
      if (allocated(radmat_1%transport) .and. radmat_1%transport_set(2) .and. (.not. fequals(radmat_1%transport(:,2), radmat_2%transport(:,2)))) then   
         equals = .false.
         error_message_radmat_equals = 'radmat_1%transport(:,2)/= radmat_2%transport(:,2)'
         return
      end if       
      if (allocated(radmat_1%transport) .and. radmat_1%transport_set(3) .and. (.not. fequals(radmat_1%transport(:,3), radmat_2%transport(:,3)))) then   
         equals = .false.
         error_message_radmat_equals = 'radmat_1%transport(:,3) /= radmat_2%transport(:,3)'
         return
      end if       
      if (allocated(radmat_1%diffusion) .and. radmat_1%diffusion_set(1) .and. (.not. fequals(radmat_1%diffusion(:,1), radmat_2%diffusion(:,1)))) then   
         equals = .false.
         error_message_radmat_equals = 'radmat_1%diffusion(:,1) /= radmat_2%diffusion(:,1)'
         return
      end if       
      if (allocated(radmat_1%diffusion) .and. radmat_1%diffusion_set(2) .and. (.not. fequals(radmat_1%diffusion(:,2), radmat_2%diffusion(:,2)))) then   
         equals = .false.
         error_message_radmat_equals = 'radmat_1%diffusion(:,2) /= radmat_2%diffusion(:,2)'
         return
      end if       
      if (allocated(radmat_1%diffusion) .and. radmat_1%diffusion_set(3) .and. (.not. fequals(radmat_1%diffusion(:,3), radmat_2%diffusion(:,3)))) then   
         equals = .false.
         error_message_radmat_equals = 'radmat_1%diffusion(:,3) /= radmat_2%diffusion(:,3)'
         return
      end if       
      if (allocated(radmat_1%fission) .and. radmat_1%fission_set .and. (.not. fequals(radmat_1%fission, radmat_2%fission))) then          
         equals = .false.
         error_message_radmat_equals = 'radmat_1%fission /= radmat_2%fission'
         return
      end if        
      if (allocated(radmat_1%production) .and. radmat_1%production_set .and. (.not. fequals(radmat_1%production, radmat_2%production))) then       
         equals = .false.
         error_message_radmat_equals = 'radmat_1%production /= radmat_2%production'
         return
      end if       
      if (allocated(radmat_1%power) .and. radmat_1%power_set .and. (.not. fequals(radmat_1%power, radmat_2%power))) then            
         equals = .false.
         error_message_radmat_equals = 'radmat_1%power /= radmat_2%power'
         return      
      end if       
      if (allocated(radmat_1%energy_released_per_fission) .and. radmat_1%energy_released_per_fission_set .and. (.not. fequals(radmat_1%energy_released_per_fission, radmat_2%energy_released_per_fission))) then            
         equals = .false.
         error_message_radmat_equals = 'admat_1%energy_released_per_fission /= radmat_2%energy_released_per_fission'
         return      
      end if       
      if (allocated(radmat_1%np_released_per_fission) .and. radmat_1%np_released_per_fission_set .and. (.not. fequals(radmat_1%np_released_per_fission, radmat_2%np_released_per_fission))) then            
         equals = .false.
         error_message_radmat_equals = 'radmat_1%np_released_per_fission /= radmat_2%np_released_per_fission'
         return      
      end if       
      if (allocated(radmat_1%prompt_spectrum) .and. radmat_1%prompt_spectrum_set .and. (.not. fequals(radmat_1%prompt_spectrum,radmat_2%prompt_spectrum))) then  
         equals = .false.
         error_message_radmat_equals = 'radmat_1%prompt_spectrum /= radmat_2%prompt_spectrum'
         return
      end if       
      if (allocated(radmat_1%velocity) .and. radmat_1%velocity_set .and. (.not. fequals(radmat_1%velocity, radmat_2%velocity))) then         
         equals = .false.
         error_message_radmat_equals = 'radmat_1%velocity /= radmat_2%velocity'
         return
      end if 
      if (allocated(radmat_1%beta) .and. radmat_1%beta_set .and. (.not. fequals(radmat_1%beta, radmat_2%beta))) then         
         equals = .false.
         error_message_radmat_equals = 'radmat_1%beta /= radmat_2%beta'
         return
      end if 

   end function radmat_equals
    
   ! --------------------------------------------------------------------------

end module radiation_materials_operations  
