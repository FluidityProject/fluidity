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

subroutine test_radiation_materials_read
   
   !!< Test the procedures contained within the module radiation_materials_read in Radiation_Materials_Read.F90
   
   use futils
   use spud
   use unittest_tools  

   use radiation_materials_data_types
   use radiation_materials_create
   use radiation_materials_destroy
   use radiation_materials_read
   use radiation_materials_operations
   use create_unittest_input 
         
   implicit none
   
   ! local variables
   logical :: has_failed,has_warned
   character(len=*), parameter :: rad_input_test_dir = 'data/'
   type(np_radmat_type) :: np_radmat_1
      
   ! none of these tests use warnings
   has_warned = .false.

   ! create an np_radmat for comparison
   call create_np_radmat_input1(np_radmat_1)
   
   call test_np_radmat_read(filename           = rad_input_test_dir//'radiation_materials_read_input1.flml', &
                            np_radmat_expected = np_radmat_1)

   call test_read_in_dataset_radmat(filename           = rad_input_test_dir//'radiation_materials_read_input1.flml', &
                                    np_radmat_expected = np_radmat_1)
         
   ! get rid of the comparison np_radmat
   call destroy(np_radmat_1)   
   
contains

   ! --------------------------------------------------------------------------

   subroutine test_np_radmat_read(filename, &
                                  np_radmat_expected)
                             
      !!< Test the procedure that reads the np_radmat using the flml options
      
      character(len=*), intent(in) :: filename      
      type(np_radmat_type), intent(in) :: np_radmat_expected      
       
      ! local variables
      type(np_radmat_type) :: np_radmat_readin
      character(len=10000) :: error_message_np_radmat_equals

      ! load the options from the correct input file
      call load_options(trim(filename))

      ! intiialise the read in data type from the options
      call create(np_radmat_readin, &
                  np_radmat_option_path = "/radiation/particle_type[0]")
                                                      
      call np_radmat_read(np_radmat_readin)         

      has_failed = .not. radmat_type_equals(np_radmat_readin, &
                                            np_radmat_expected, &
                                            error_message_np_radmat_equals) 

      call report_test("[test_np_radmat_read]", &
                       has_failed, &
                       has_warned, &
                       "failed as np_radmat values not the same, with error message "//trim(error_message_np_radmat_equals))   
            
      ! remove the options
      call clear_options
      
      call destroy(np_radmat_readin)
                                                   
   end subroutine test_np_radmat_read
   
   ! --------------------------------------------------------------------------

   subroutine test_read_in_dataset_radmat(filename, &
                                          np_radmat_expected)
                             
      !!< Test the procedure that reads the dataset_radmat using the flml options
      
      character(len=*), intent(in) :: filename      
      type(np_radmat_type), intent(in) :: np_radmat_expected      
       
      ! local variables
      type(np_radmat_type) :: np_radmat_readin
      character(len=10000) :: error_message_dataset_radmat_equals

      ! load the options from the correct input file
      call load_options(trim(filename))

      ! intiialise the read in data type from the options
      call create(np_radmat_readin, &
                  np_radmat_option_path = "/radiation/particle_type[0]")
                                                      
      call read_in_dataset_radmat(np_radmat_readin%dataset_radmats(1), &
                                  np_radmat_readin%delayed_lambda_spectrum, &
                                  read_delayed_lambda_spectrum = .true., &
                                  read_velocity_data           = .false., &
                                  read_power_data              = .true.)         

      has_failed = .not. radmat_type_equals(np_radmat_readin%dataset_radmats(1), &
                                            np_radmat_expected%dataset_radmats(1), &
                                            error_message_dataset_radmat_equals) 

      call report_test("[test_read_in_dataset_radmat]", &
                       has_failed, &
                       has_warned, &
                       "failed as np_radmat values not the same, with error message "//trim(error_message_dataset_radmat_equals))   
            
      ! remove the options
      call clear_options
      
      call destroy(np_radmat_readin)
                                                   
   end subroutine test_read_in_dataset_radmat
   
   ! --------------------------------------------------------------------------

end subroutine test_radiation_materials_read
