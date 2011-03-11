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

subroutine test_radiation_materials_read_format_radmats
   
   !!< Test the procedures contained within the module radiation_materials_read_format_radmats 
   !!< in Radiation_Materials_Read_Format_Radmats.F90
   
   use futils
   use unittest_tools  

   use radiation_materials_read_format_radmats
   use radiation_materials_data_types
   use radiation_materials_create
   use radiation_materials_destroy
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
   
   call test_read_format_radmats(filename                         = rad_input_test_dir//'radiation_materials_read_input1.radmats', &
                                 read_delayed_lambda_spectrum     = .true., &
                                 read_velocity_data               = .false., &
                                 read_power_data                  = .true., &
                                 problem_dimension                = 3, &
                                 record_len                       = 132, &
                                 particle_radmat_size             = particle_radmat_1%particle_radmat_size, &
                                 dataset_radmat_expected          = particle_radmat_1%dataset_radmats(1), &
                                 delayed_lambda_spectrum_expected = particle_radmat_1%delayed_lambda_spectrum)
         
   ! get rid of the comparison particle_radmat 
   call destroy(particle_radmat_1)   
   
contains

   ! --------------------------------------------------------------------------

   subroutine test_read_format_radmats(filename, &
                                       read_delayed_lambda_spectrum, &
                                       read_velocity_data, &
                                       read_power_data, &
                                       problem_dimension, &
                                       record_len, &
                                       particle_radmat_size, &
                                       dataset_radmat_expected, &
                                       delayed_lambda_spectrum_expected)
                             
      !!< Test the procedure that reads the dataset_radmat in a format_radmats file via comparing to a reference 
      !!< This may also include a read in of the delayed lambda and spectrum
      
      character(len=*), intent(in) :: filename      
      logical, intent(in) :: read_delayed_lambda_spectrum
      logical, intent(in) :: read_velocity_data
      logical, intent(in) :: read_power_data
      integer, intent(in) :: problem_dimension
      integer, intent(in) :: record_len
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size      
      type(dataset_radmat_type), intent(in) :: dataset_radmat_expected
      type(delayed_lambda_spectrum_type), intent(in) :: delayed_lambda_spectrum_expected  
       
      ! local variables
      type(dataset_radmat_type) :: dataset_radmat_format_radmats
      type(delayed_lambda_spectrum_type) :: delayed_lambda_spectrum_format_radmats  
      character(len=10000) :: error_message_delayed_lambda_spectrum_equals
      character(len=10000) :: error_message_dataset_radmat_equals

      ! intiialise the read in data type
      call allocate(dataset_radmat_format_radmats, &
                    particle_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_format_radmats)
      
      ! set the option paths 
      dataset_radmat_format_radmats%option_path = '/radiation/particle_type[0]/radiation_material_data_set_from_file[0]'
      dataset_radmat_format_radmats%physical_radmats(1)%option_path =  '/radiation/particle_type[0]/radiation_material_data_set_from_file[0]/physical_material[0]'
      
      call allocate(delayed_lambda_spectrum_format_radmats, &
                    particle_radmat_size%number_of_delayed_groups, &
                    particle_radmat_size%number_of_energy_groups)

      call zero(delayed_lambda_spectrum_format_radmats)
                                                
      call read_format_radmats(trim(filename), &
                               dataset_radmat_format_radmats, &
                               delayed_lambda_spectrum_format_radmats, &
                               read_delayed_lambda_spectrum, &
                               read_velocity_data, &
                               read_power_data, &
                               problem_dimension, &
                               record_len)         

      has_failed = .not. radmat_type_equals(delayed_lambda_spectrum_format_radmats, &
                                            delayed_lambda_spectrum_expected, &
                                            error_message_delayed_lambda_spectrum_equals) 

      call report_test("[test_read_format_radmats check delayed_lambda_spectrum values]", &
                       has_failed, &
                       has_warned, &
                       "failed as delayed_lambda_spectrum values not the same, with error message "//trim(error_message_delayed_lambda_spectrum_equals))   
      
      has_failed = .not. radmat_type_equals(dataset_radmat_format_radmats, &
                                            dataset_radmat_expected, &
                                            error_message_dataset_radmat_equals) 

      call report_test("[test_read_format_radmats check dataset_radmat values]", &
                       has_failed, &
                       has_warned, &
                       "failed as dataset_radmat values not the same, with error message "//trim(error_message_dataset_radmat_equals))   
      
      call destroy(dataset_radmat_format_radmats)
      call destroy(delayed_lambda_spectrum_format_radmats)
                                                   
   end subroutine test_read_format_radmats
   
   ! --------------------------------------------------------------------------

end subroutine test_radiation_materials_read_format_radmats
