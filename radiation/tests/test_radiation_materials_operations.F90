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

subroutine test_radiation_materials_operations
 
   !!< Test the procedures contained within the module radiation_materials_operations 
   !!< in Radiation_Materials_Operations.F90

   use futils   
   use unittest_tools 
    
   use radiation_materials_data_types
   use radiation_materials_create
   use radiation_materials_destroy
   use radiation_materials_operations
   use create_unittest_input
       
   implicit none

   ! local variables
   logical :: has_failed,has_warned
   type(particle_radmat_type) :: particle_radmat_1
   type(particle_radmat_type) :: particle_radmat_2

   ! none of these tests use warnings
   has_warned = .false.
   
   ! create two particle_radmat
   call create_particle_radmat_input1(particle_radmat_1)
   call create_particle_radmat_input1(particle_radmat_2)
   
   ! now compare
   call test_radmat_equals(particle_radmat_1%dataset_radmats(1)%physical_radmats(1)%radmats(1), &
                           particle_radmat_2%dataset_radmats(1)%physical_radmats(1)%radmats(1))

   call test_physical_radmat_equals(particle_radmat_1%dataset_radmats(1)%physical_radmats(1), &
                                    particle_radmat_2%dataset_radmats(1)%physical_radmats(1))

   call test_dataset_radmat_equals(particle_radmat_1%dataset_radmats(1), &
                                   particle_radmat_2%dataset_radmats(1))

   call test_delayed_lambda_spectrum_equals(particle_radmat_1%delayed_lambda_spectrum, &
                                            particle_radmat_2%delayed_lambda_spectrum)

   call test_particle_radmat_size_equals(particle_radmat_1%particle_radmat_size, &
                                   particle_radmat_2%particle_radmat_size)

   call test_particle_radmat_equals(particle_radmat_1, &
                                    particle_radmat_2)
   
   ! destroy the particle_radmat 
   call destroy(particle_radmat_1)   
   call destroy(particle_radmat_2)   

contains

   ! --------------------------------------------------------------------------
   
   subroutine test_radmat_equals(radmat_1, &
                                 radmat_2)
   
      !!< Test the procedure that checks that two radmat are equal
      
      type(radmat_type), intent(in) :: radmat_1
      type(radmat_type), intent(in) :: radmat_2
      
      ! local variable
      character(len=10000) :: error_message_radmat_equals
            
      has_failed = .not. radmat_type_equals(radmat_1, &
                                            radmat_2, &
                                            error_message_radmat_equals) 

      call report_test("[test_radmat_equals]", &
                       has_failed, &
                       has_warned, &
                       "failed as values not he same, with error message "//trim(error_message_radmat_equals))   
      
   end subroutine test_radmat_equals
   
   ! --------------------------------------------------------------------------
   
   subroutine test_physical_radmat_equals(physical_radmat_1, &
                                          physical_radmat_2)
   
      !!< Test the procedure that checks that two physical_radmat are equal
      
      type(physical_radmat_type), intent(in) :: physical_radmat_1
      type(physical_radmat_type), intent(in) :: physical_radmat_2

      ! local variable
      character(len=10000) :: error_message_physical_radmat_equals
            
      has_failed = .not. radmat_type_equals(physical_radmat_1, &
                                            physical_radmat_2, &
                                            error_message_physical_radmat_equals) 
      
      call report_test("[test_physical_radmat_equals]", &
                       has_failed, &
                       has_warned, &
                       "failed as values not he same, with error message "//trim(error_message_physical_radmat_equals))   
      
   end subroutine test_physical_radmat_equals
   
   ! --------------------------------------------------------------------------
   
   subroutine test_dataset_radmat_equals(dataset_radmat_1, &
                                         dataset_radmat_2)
   
      !!< Test the procedure that checks that two dataset_radmat are equal
      
      type(dataset_radmat_type), intent(in) :: dataset_radmat_1
      type(dataset_radmat_type), intent(in) :: dataset_radmat_2

      ! local variable
      character(len=10000) :: error_message_dataset_radmat_equals
            
      has_failed = .not. radmat_type_equals(dataset_radmat_1, &
                                            dataset_radmat_2, &
                                            error_message_dataset_radmat_equals) 
      
      call report_test("[test_dataset_radmat_equals]", &
                       has_failed, &
                       has_warned, &
                       "failed as values not he same, with error message "//trim(error_message_dataset_radmat_equals))   
      
   end subroutine test_dataset_radmat_equals
   
   ! --------------------------------------------------------------------------
   
   subroutine test_delayed_lambda_spectrum_equals(delayed_lambda_spectrum_1, &
                                                  delayed_lambda_spectrum_2)
   
      !!< Test the procedure that checks that two delayed_lambda_spectrum are equal
      
      type(delayed_lambda_spectrum_type), intent(in) :: delayed_lambda_spectrum_1
      type(delayed_lambda_spectrum_type), intent(in) :: delayed_lambda_spectrum_2

      ! local variable
      character(len=10000) :: error_message_delayed_lambda_spectrum_equals
            
      has_failed = .not. radmat_type_equals(delayed_lambda_spectrum_1, &
                                            delayed_lambda_spectrum_2, &
                                            error_message_delayed_lambda_spectrum_equals) 
      
      call report_test("[test_delayed_lambda_spectrum_equals]", &
                       has_failed, &
                       has_warned, &
                       "failed as values not he same, with error message "//trim(error_message_delayed_lambda_spectrum_equals))   
      
   end subroutine test_delayed_lambda_spectrum_equals 
   
   ! --------------------------------------------------------------------------
   
   subroutine test_particle_radmat_size_equals(particle_radmat_size_1, &
                                               particle_radmat_size_2)
   
      !!< Test the procedure that checks that two particle_radmat_size are equal
      
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size_1
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size_2

      ! local variable
      character(len=10000) :: error_message_particle_radmat_size_equals
            
      has_failed = .not. radmat_type_equals(particle_radmat_size_1, &
                                            particle_radmat_size_2, &
                                            error_message_particle_radmat_size_equals) 
      
      call report_test("[test_particle_radmat_size_equals]", &
                       has_failed, &
                       has_warned, &
                       "failed as values not he same, with error message "//trim(error_message_particle_radmat_size_equals))   
      
   end subroutine test_particle_radmat_size_equals
   
   ! --------------------------------------------------------------------------
   
   subroutine test_particle_radmat_equals(particle_radmat_1, &
                                          particle_radmat_2)
   
      !!< Test the procedure that checks that two particle_radmat are equal
      
      type(particle_radmat_type), intent(in) :: particle_radmat_1
      type(particle_radmat_type), intent(in) :: particle_radmat_2

      ! local variable
      character(len=10000) :: error_message_particle_radmat_equals
            
      has_failed = .not. radmat_type_equals(particle_radmat_1, &
                                            particle_radmat_2, &
                                            error_message_particle_radmat_equals) 
      
      call report_test("[test_particle_radmat_equals]", &
                       has_failed,has_warned, &
                       "failed as values not he same, with error message "//trim(error_message_particle_radmat_equals))   
      
   end subroutine test_particle_radmat_equals
   
   ! --------------------------------------------------------------------------

end subroutine test_radiation_materials_operations
