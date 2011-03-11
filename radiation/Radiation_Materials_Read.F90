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
   use radiation_materials_read_format_radmats    
   
   implicit none
   
   private 

   public :: particle_radmat_read, &
             read_in_dataset_radmat 
        
contains

   ! --------------------------------------------------------------------------

   subroutine particle_radmat_read(particle_radmat)
      
      !!< Read and set the data into the particle radmat type
      !!< Initially checks that this type is created - if not abort
      
      type(particle_radmat_type), intent(inout) :: particle_radmat
      
      ! local variables
      integer :: dmat 
      logical :: read_delayed_lambda_spectrum
      logical :: read_velocity_data
      logical :: read_power_data  
      character(len=OPTION_PATH_LEN) :: dataset_file_name
      character(len=OPTION_PATH_LEN) :: dataset_file_name_to_read_delayed_lambda
      
      ewrite(1,*) 'Read in radiation particle_radmat values for',trim(particle_radmat%name)
      
      created: if (.not. particle_radmat%created) then 
         
         FLAbort('Cannot read in a particle_radmat type if it is not created yet')
         
      end if created
                  
      ! see if there is a need for the velocity data for this object (ie. time run) 
      read_velocity_data = .false.      
      if (have_option(trim(particle_radmat%option_path)//'/time_run')) read_velocity_data = .true.
         
      ! see if there is a need for the power data for this object 
      check_need_power: if (have_option(trim(particle_radmat%option_path)//'/eigenvalue_run/flux_normalisation/total_power')) then
            
         read_power_data = .true.
            
      else check_need_power
         
         read_power_data = .false.
         
      end if check_need_power         
      
      ewrite(1,*) 'Read in velocity: ',read_velocity_data
      ewrite(1,*) 'Read in power: ',read_power_data 
      
      ! check which data set we read the delayed lambda and spectrum from
      dataset_file_name_to_read_delayed_lambda = ''
      have_delayed: if (have_option(trim(particle_radmat%option_path)//'/delayed_neutron_precursor')) then
      
         call get_option(trim(particle_radmat%option_path)//'/delayed_neutron_precursor/read_delayed_lambda_spectrum_from_data_set/file_name',&
                         &dataset_file_name_to_read_delayed_lambda)
      
      end if have_delayed
 
      data_set_loop: do dmat = 1,size(particle_radmat%dataset_radmats)
         
         ! read in the materials for this data set      
         read_delayed_lambda_spectrum = .false.          
         call get_option(trim(particle_radmat%dataset_radmats(dmat)%option_path)//'/file_name',dataset_file_name)
         
         ewrite(1,*) 'Read in dataset ',trim(dataset_file_name)
                           
         if (trim(dataset_file_name) == trim(dataset_file_name_to_read_delayed_lambda)) read_delayed_lambda_spectrum = .true.
       
         ewrite(1,*) 'Read in delayed lambda/spectrum: ',read_delayed_lambda_spectrum
                                        
         call read_in_dataset_radmat(particle_radmat%dataset_radmats(dmat), &
                                     particle_radmat%delayed_lambda_spectrum, &
                                     read_delayed_lambda_spectrum, &
                                     read_velocity_data, &
                                     read_power_data )
      
      end do data_set_loop
      
      particle_radmat%readin = .true.
      
      ewrite(1,*) 'Finished read in radiation particle_radmat values for ',trim(particle_radmat%name)
      
   end subroutine particle_radmat_read

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
    
      read_data_set_if: if (have_option(trim(dataset_radmat%option_path)//'/format_radmats')) then

         ! find the format radmats file name from the options
         call get_option(trim(dataset_radmat%option_path)//'/file_name',filename)

         ! set the radmats file record length used for reading from the options
         call get_option(trim(dataset_radmat%option_path)//'/format_radmats/maximum_record_length',record_len) 
                  
         call read_format_radmats(trim(filename), &
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
