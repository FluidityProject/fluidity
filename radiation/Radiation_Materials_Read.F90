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

module radiation_materials_read

   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud

   use radiation_materials_data_types 
   use radiation_materials_read_wims9plus    
   
   implicit none
   
   private 

   public :: np_radmat_read, &
             read_in_dataset_radmat 
        
contains

   ! --------------------------------------------------------------------------

   subroutine np_radmat_read(np_radmat)
      
      !!< Read and set the data into the neutral particle radmat object
      !!< Initially checks that this object is created - if not abort
      
      type(np_radmat_type), intent(inout) :: np_radmat
      
      ! local variables
      integer :: dmat 
      logical :: read_delayed_lambda_spectrum
      logical :: read_velocity_data
      logical :: read_power_data  
      character(len=OPTION_PATH_LEN) :: dataset_file_name
      character(len=OPTION_PATH_LEN) :: dataset_file_name_to_read_delayed_lambda
      
      ewrite(1,*) 'Read in radiation np_radmat values for',trim(np_radmat%option_path)
      
      created: if (.not. np_radmat%created) then 
         
         FLAbort('Cannot read in a np_radmat type if it is not created yet')
         
      end if created
                  
      ! see if there is a need for the velocity data for this object (ie. time run) -set the module variable
      read_velocity_data = .false.      
      if (have_option(trim(np_radmat%option_path)//'/time_run')) read_velocity_data = .true.
         
      ! see if there is a need for the power data for this object (ie. RadPowerTemplate exists and/or an eigenvalue run 
      ! is normalised to a prescribed power and/or RadPowerTemplate has been invoked in the input - set the module variable
      check_need_power: if (have_option(trim(np_radmat%option_path)//'/eigenvalue_run/flux_normalisation/total_power') .or. &
                            have_option(trim(np_radmat%option_path)//'/scalar_field::RadPowerTemplate') .or. &
                            have_option(trim(np_radmat%option_path)//'/region_id_material_mapping/scalar_field::RadPowerTemplate') .or. &
                            have_option(trim(np_radmat%option_path)//'/link_with_multimaterial/scalar_field::RadPowerTemplate') .or. &
                            have_option(trim(np_radmat%option_path)//'/link_with_porous_media/scalar_field::RadPowerTemplate')) then
            
         read_power_data = .true.
            
      else check_need_power
         
         read_power_data = .false.
         
      end if check_need_power         
      
      ewrite(1,*) 'Read in velocity: ',read_velocity_data
      ewrite(1,*) 'Read in power: ',read_power_data 
      
      ! check which data set we read the delayed lambda and spectrum from
      dataset_file_name_to_read_delayed_lambda = ''
      have_delayed: if (have_option(trim(np_radmat%option_path)//'/delayed_neutron_precursor')) then
      
         call get_option(trim(np_radmat%option_path)//'/delayed_neutron_precursor/read_delayed_lambda_spectrum_from_data_set/file_name',&
                         &dataset_file_name_to_read_delayed_lambda)
      
      end if have_delayed
 
      Data_set_loop: do dmat = 1,size(np_radmat%dataset_radmats)
         
         ! read in the materials for this data set      
         read_delayed_lambda_spectrum = .false.          
         call get_option(trim(np_radmat%dataset_radmats(dmat)%option_path)//'/file_name',dataset_file_name)
         
         ewrite(1,*) 'Read in dataset ',trim(dataset_file_name)
                           
         if (trim(dataset_file_name) == trim(dataset_file_name_to_read_delayed_lambda)) read_delayed_lambda_spectrum = .true.
       
         ewrite(1,*) 'Read in delayed lambda/spectrum: ',read_delayed_lambda_spectrum
                                        
         call read_in_dataset_radmat(np_radmat%dataset_radmats(dmat), &
                                     np_radmat%delayed_lambda_spectrum, &
                                     read_delayed_lambda_spectrum, &
                                     read_velocity_data, &
                                     read_power_data )
      
      end do Data_set_loop
      
      np_radmat%readin = .true.
      
      ewrite(1,*) 'Finished read in radiation np_radmat values for',trim(np_radmat%option_path)
      
   end subroutine np_radmat_read

   ! --------------------------------------------------------------------------
            
   subroutine read_in_dataset_radmat(dataset_radmat, &
                                     delayed_lambda_spectrum, &
                                     read_delayed_lambda_spectrum, &
                                     read_velocity_data, & 
                                     read_power_data)
      
      !!< Read in a specific dataset_radmat format.
      
      type(dataset_radmat_type), intent(inout) :: dataset_radmat
      type(delayed_lambda_spectrum_type), intent(inout) :: delayed_lambda_spectrum
      logical, intent(in) :: read_delayed_lambda_spectrum
      logical, intent(in) :: read_velocity_data
      logical, intent(in) :: read_power_data
      
      ! lcocal variables
      character(len=OPTION_PATH_LEN+15) :: filename
      integer :: problem_dimension
      integer :: record_len 
      
      ! get the problem dimension from the options
      call get_option("/geometry/dimension",problem_dimension) 
    
      read_data_set_if: if (have_option(trim(dataset_radmat%option_path)//'/format_wims9plus')) then

         ! find the wim9 file name from the options
         call get_option(trim(dataset_radmat%option_path)//'/file_name',filename)
         filename = trim(filename) // ".wims9plus_data"

         ! set the wims9plus record length used for reading from the options
         call get_option(trim(dataset_radmat%option_path)//"/format_wims9plus/maximum_record_length",record_len) 
                  
         call read_wims9plus(trim(filename), &
                             dataset_radmat, &
                             delayed_lambda_spectrum, &
                             read_delayed_lambda_spectrum, &
                             read_velocity_data, &
                             read_power_data, &
                             problem_dimension, &
                             record_len)         
      
      else read_data_set_if 
      
         FLAbort('Unknown radiation material data read in file format')
                       
      end if read_data_set_if     
        
   end subroutine read_in_dataset_radmat
   
   ! --------------------------------------------------------------------------

end module radiation_materials_read   
