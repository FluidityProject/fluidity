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
   type(particle_radmat_type) :: particle_radmat_1
      
   ! none of these tests use warnings
   has_warned = .false.

   ! create an particle_radmat for comparison
   call create_particle_radmat_input1(particle_radmat_1)
   
   call test_particle_radmat_read(filename                 = rad_input_test_dir//'radiation_materials_read_input1.flml', &
                                  particle_radmat_expected = particle_radmat_1)

   call test_read_in_dataset_radmat(filename                 = rad_input_test_dir//'radiation_materials_read_input1.flml', &
                                    particle_radmat_expected = particle_radmat_1)
         
   ! get rid of the comparison particle_radmat
   call destroy(particle_radmat_1)   
   
contains

   ! --------------------------------------------------------------------------

   subroutine test_particle_radmat_read(filename, &
                                        particle_radmat_expected)
                             
      !!< Test the procedure that reads the particle_radmat using the flml options
      
      character(len=*), intent(in) :: filename      
      type(particle_radmat_type), intent(in) :: particle_radmat_expected      
       
      ! local variables
      type(particle_radmat_type) :: particle_radmat_readin
      character(len=10000) :: error_message_particle_radmat_equals

      ! load the options from the correct input file
      call load_options(trim(filename))

      ! intiialise the read in data type from the options
      call create(particle_radmat_readin, &
                  particle_option_path = "/embedded_models/radiation/particle_type[0]", &
                  particle_name = "neutron")
                                                      
      call particle_radmat_read(particle_radmat_readin)         

      has_failed = .not. radmat_type_equals(particle_radmat_readin, &
                                            particle_radmat_expected, &
                                            error_message_particle_radmat_equals) 

      call report_test("[test_particle_radmat_read]", &
                       has_failed, &
                       has_warned, &
                       "failed as particle_radmat values not the same, with error message "//trim(error_message_particle_radmat_equals))   
            
      ! remove the options
      call clear_options
      
      call destroy(particle_radmat_readin)
                                                   
   end subroutine test_particle_radmat_read
   
   ! --------------------------------------------------------------------------

   subroutine test_read_in_dataset_radmat(filename, &
                                          particle_radmat_expected)
                             
      !!< Test the procedure that reads the dataset_radmat using the flml options
      
      character(len=*), intent(in) :: filename      
      type(particle_radmat_type), intent(in) :: particle_radmat_expected      
       
      ! local variables
      type(particle_radmat_type) :: particle_radmat_readin
      character(len=10000) :: error_message_dataset_radmat_equals

      ! load the options from the correct input file
      call load_options(trim(filename))

      ! intiialise the read in data type from the options
      call create(particle_radmat_readin, &
                  particle_option_path = "/embedded_models/radiation/particle_type[0]", &
                  particle_name = "neutron")
                                                      
      call read_in_dataset_radmat(particle_radmat_readin%dataset_radmats(1), &
                                  particle_radmat_readin%delayed_lambda_spectrum, &
                                  read_delayed_lambda_spectrum = .true., &
                                  read_velocity_data           = .false., &
                                  read_power_data              = .true.)         

      has_failed = .not. radmat_type_equals(particle_radmat_readin%dataset_radmats(1), &
                                            particle_radmat_expected%dataset_radmats(1), &
                                            error_message_dataset_radmat_equals) 

      call report_test("[test_read_in_dataset_radmat]", &
                       has_failed, &
                       has_warned, &
                       "failed as particle_radmat values not the same, with error message "//trim(error_message_dataset_radmat_equals))   
            
      ! remove the options
      call clear_options
      
      call destroy(particle_radmat_readin)
                                                   
   end subroutine test_read_in_dataset_radmat
   
   ! --------------------------------------------------------------------------

end subroutine test_radiation_materials_read
