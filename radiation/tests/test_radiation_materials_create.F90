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

subroutine test_radiation_materials_create
 
   !!< Test the procedures contained within the module radiation_materials_create in Radiation_Materials_Create.F90

   use futils   
   use spud
   use unittest_tools 
    
   use radiation_materials_data_types
   use radiation_materials_create
   use radiation_materials_destroy
       
   implicit none

   ! local variables
   logical :: has_failed,has_warned
   character(len=*), parameter :: rad_input_test_dir = 'data/'
   
   character(len=10000) :: filename,error_message      
   integer :: number_of_energy_groups_expected
   integer :: number_of_delayed_groups_expected
   integer :: total_number_dataset_radmats_expected
   integer :: total_number_physical_radmats_expected
   integer :: total_number_radmats_expected
   integer, dimension(:), allocatable :: number_of_physical_radmats_expected 
   integer, dimension(:), allocatable :: number_of_radmats_expected 
   integer, dimension(:), allocatable :: number_of_radmats_base_expected 
   integer, dimension(:), allocatable :: number_of_scatter_moments_expected     
   character(len=10000) np_radmat_option_path_expected    
   character(len=10000), dimension(:), allocatable :: dataset_radmats_option_path_expected    
   character(len=10000), dimension(:), allocatable :: physical_radmats_option_path_expected    

   ! none of these tests use warnings
   has_warned = .false.
   
   ! test all the procedures in one call to the main procedure that will induce a call to all the others   
   filename = rad_input_test_dir//'radiation_materials_read_input1.flml'
   number_of_energy_groups_expected       = 2
   number_of_delayed_groups_expected      = 6
   total_number_dataset_radmats_expected  = 1
   total_number_physical_radmats_expected = 1
   total_number_radmats_expected          = 10
   allocate(number_of_physical_radmats_expected(total_number_dataset_radmats_expected))
   allocate(number_of_radmats_expected(total_number_physical_radmats_expected))
   allocate(number_of_radmats_base_expected(total_number_dataset_radmats_expected))
   allocate(number_of_scatter_moments_expected(total_number_dataset_radmats_expected))   
   number_of_physical_radmats_expected = (/1/)
   number_of_radmats_expected          = (/10/)
   number_of_radmats_base_expected     = (/1,11/)
   number_of_scatter_moments_expected  = (/2/)
   np_radmat_option_path_expected      = '/radiation/particle_type[0]'
   allocate(dataset_radmats_option_path_expected(total_number_dataset_radmats_expected))
   allocate(physical_radmats_option_path_expected(total_number_physical_radmats_expected))
   dataset_radmats_option_path_expected  = (/'/radiation/particle_type[0]/radiation_material_data_set_from_file[0]'/)
   physical_radmats_option_path_expected = (/'/radiation/particle_type[0]/radiation_material_data_set_from_file[0]/physical_material[0]'/)
   call test_np_radmat_create_from_options()
   call report_test("[test_np_radmat_create_from_options for filename "//trim(filename)//"]", &
                    has_failed, &
                    has_warned, &
                    "failed with error message "//trim(error_message))   
   deallocate(number_of_physical_radmats_expected)
   deallocate(number_of_radmats_expected)
   deallocate(number_of_radmats_base_expected)
   deallocate(number_of_scatter_moments_expected)   
   deallocate(dataset_radmats_option_path_expected)
   deallocate(physical_radmats_option_path_expected)
            
contains

   ! --------------------------------------------------------------------------
   
   subroutine test_np_radmat_create_from_options()
                             
      !!< Test the procedure that creates the np_radmat using the flml options
      !!< Create involves setting the size, allocating, zeroing and set option path
                   
      ! local variables
      type(np_radmat_type) :: np_radmat_create
      integer :: dmat,pmat,rmat,m
      
      ! load the options from the correct input file
      call load_options(trim(filename))

      ! create the data type from the options
      call create(np_radmat_create, &
                  np_radmat_option_path = "/radiation/particle_type[0]")
      
      error_message = 'no error'
                                                      
      has_failed = .false.

      if (.not. np_radmat_create%np_radmat_size%size_set) then
         has_failed = .true.
         error_message = '.not. np_radmat_create%np_radmat_size%size_set'
         return
      end if
      
      if (number_of_energy_groups_expected /= np_radmat_create%np_radmat_size%number_of_energy_groups) then
         has_failed = .true.
         error_message = 'number_of_energy_groups_expected /= np_radmat_create%np_radmat_size%number_of_energy_groups'
         return
      end if
      
      if (number_of_delayed_groups_expected /= np_radmat_create%np_radmat_size%number_of_delayed_groups) then
         has_failed = .true.
         error_message = 'number_of_delayed_groups_expected /= np_radmat_create%np_radmat_size%number_of_delayed_groups'
         return
      end if
      
      if (total_number_dataset_radmats_expected /= np_radmat_create%np_radmat_size%total_number_dataset_radmats) then
         has_failed = .true.
         error_message = 'total_number_dataset_radmats_expected /= np_radmat_create%np_radmat_size%total_number_dataset_radmats'
         return
      end if
      
      if (total_number_physical_radmats_expected /= np_radmat_create%np_radmat_size%total_number_physical_radmats) then
         has_failed = .true.
         error_message = 'total_number_physical_radmats_expected /= np_radmat_create%np_radmat_size%total_number_physical_radmats'
         return
      end if
      
      if (total_number_radmats_expected /= np_radmat_create%np_radmat_size%total_number_radmats) then
         has_failed = .true.
         error_message = 'total_number_radmats_expected /= np_radmat_create%np_radmat_size%total_number_radmats'
         return
      end if
      
      do dmat = 1,total_number_dataset_radmats_expected  
         if (number_of_physical_radmats_expected(dmat) /= np_radmat_create%np_radmat_size%number_of_physical_radmats(dmat)) then
            has_failed = .true.
            error_message = 'number_of_physical_radmats_expected('//int2str(dmat)//') /= np_radmat_create%np_radmat_size%number_of_physical_radmats('//int2str(dmat)//')'
            return
         end if
      end do

      do pmat = 1,total_number_physical_radmats_expected       
         if (number_of_radmats_expected(pmat) /= np_radmat_create%np_radmat_size%number_of_radmats(pmat)) then
            has_failed = .true.
            error_message = 'number_of_radmats_expected('//int2str(pmat)//') /= np_radmat_create%np_radmat_size%number_of_radmats('//int2str(pmat)//')'
            return
         end if
      end do

      do dmat = 1,total_number_dataset_radmats_expected        
         if (number_of_radmats_base_expected(dmat) /= np_radmat_create%np_radmat_size%number_of_radmats_base(dmat)) then
            has_failed = .true.
            error_message = 'number_of_radmats_base_expected('//int2str(dmat)//') /= np_radmat_create%np_radmat_size%number_of_radmats_base('//int2str(dmat)//')'
            return
         end if
      end do      

      do dmat = 1,total_number_dataset_radmats_expected           
         if (number_of_scatter_moments_expected(dmat) /= np_radmat_create%np_radmat_size%number_of_scatter_moments(dmat)) then
            has_failed = .true.
            error_message = 'number_of_scatter_moments_expected('//int2str(dmat)//') /= np_radmat_create%np_radmat_size%number_of_scatter_moments('//int2str(dmat)//')'
            return
         end if
      end do
                  
      if (np_radmat_option_path_expected /= np_radmat_create%option_path) then
         has_failed = .true.
         error_message = 'np_radmat_option_path_expected /= np_radmat_create%option_path'
         return
      end if
      
      do dmat = 1,total_number_dataset_radmats_expected
         if (trim(dataset_radmats_option_path_expected(dmat)) /= trim(np_radmat_create%dataset_radmats(dmat)%option_path)) then
            has_failed = .true.
            error_message = 'trim(dataset_radmats_option_path_expected('//int2str(dmat)//')) /= trim(np_radmat_create%dataset_radmats('//int2str(dmat)//')%option_path)'
            return
         end if
      end do

      do dmat = 1,total_number_dataset_radmats_expected    
         do pmat = 1,number_of_physical_radmats_expected(dmat)       
            if (trim(physical_radmats_option_path_expected(number_of_radmats_base_expected(dmat) - 1 + pmat)) /=  &
                trim(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%option_path)) then
               has_failed = .true.
               error_message = 'trim(physical_radmats_option_path_expected('//int2str(number_of_radmats_base_expected(dmat) - 1 + pmat)//')) /= '//&
                               &'trim(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')%physical_radmats_option_path)'
               return
            end if
         end do
      end do
      
      if (size(np_radmat_create%dataset_radmats) /= total_number_dataset_radmats_expected) then
         has_failed = .true.
         error_message = 'size(np_radmat_create%dataset_radmats) /= total_number_dataset_radmats_expected'
         return
      end if
      
      do dmat = 1,total_number_dataset_radmats_expected       
         if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats) /= number_of_physical_radmats_expected(dmat)) then
            has_failed = .true.
            error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats) /= number_of_physical_radmats_expected('//int2str(dmat)//')'
            return
         end if      
      end do 
      
      do dmat = 1,total_number_dataset_radmats_expected       
         do pmat = 1,number_of_physical_radmats_expected(dmat)      
            if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats) /= number_of_radmats_expected(number_of_radmats_base_expected(dmat) - 1 + pmat)) then
               has_failed = .true.
               error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')%radmats) /='//&
                               &'number_of_radmats_expected('//int2str(number_of_radmats_base_expected(dmat) - 1 + pmat)//')'
               return
            end if
         end do      
      end do 
      
      do dmat = 1,total_number_dataset_radmats_expected       
         do pmat = 1,number_of_physical_radmats_expected(dmat)      
            do rmat = 1,number_of_radmats_expected(number_of_radmats_base_expected(dmat) - 1 + pmat)
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%total) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%total) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%absorption) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%absorption) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%scatter,1) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%scatter,1) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%scatter,2) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%scatter,2) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%scatter,3) /= number_of_scatter_moments_expected(dmat)) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%scatter,3) /= number_of_scatter_moments_expected('//int2str(dmat)//')'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%removal,1) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%removal,1) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%removal,2) /= number_of_scatter_moments_expected(dmat)) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%removal,2) /= number_of_scatter_moments_expected('//int2str(dmat)//')'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%transport,1) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%transport,1) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%transport,2) /= 3) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%transport,2) /= 3)'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%diffusion,1) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%diffusion,1) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%diffusion,2) /= 3) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%diffusion,2) /= 3)'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%fission) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%fission) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%production) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%production) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%power) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%power) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%energy_released_per_fission) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%np_released_per_fission) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%np_released_per_fission) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%prompt_spectrum) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%prompt_spectrum) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%velocity) /= number_of_energy_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%velocity) /= number_of_energy_groups_expected'
                  return
               end if
               if (size(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%beta) /= number_of_delayed_groups_expected) then
                  has_failed = .true.
                  error_message = 'size(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%beta) /= number_of_delayed_groups_expected'
                  return
               end if


               if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%total,0.0)) then
                  has_failed = .true.
                  error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%total,0.0)'
                  return
               end if
               if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%absorption,0.0)) then
                  has_failed = .true.
                  error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%absorption,0.0)'
                  return
               end if
               do m = 1,number_of_scatter_moments_expected(dmat)
                  if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%scatter(:,:,m),0.0)) then
                     has_failed = .true.
                     error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                     &'%radmats('//int2str(rmat)//')%scatter(:,:,'//int2str(m)//'),0.0)'
                     return
                  end if
               end do 
               if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%removal,0.0)) then
                  has_failed = .true.
                  error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%removal,0.0)'
                  return
               end if
               if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%transport,0.0)) then
                  has_failed = .true.
                  error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%transport,0.0)'
                  return
               end if
               if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%diffusion,0.0)) then
                  has_failed = .true.
                  error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%diffusion,0.0)'
                  return
               end if
               if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%fission,0.0)) then
                  has_failed = .true.
                  error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%fission,0.0)'
                  return
               end if
               if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%production,0.0)) then
                  has_failed = .true.
                  error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%production,0.0)'
                  return
               end if
               if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%power,0.0)) then
                  has_failed = .true.
                  error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%power,0.0)'
                  return
               end if
               if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission,0.0)) then
                  has_failed = .true.
                  error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%energy_released_per_fission,0.0)'
                  return
               end if
               if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%np_released_per_fission,0.0)) then
                  has_failed = .true.
                  error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%np_released_per_fission,0.0)'
                  return
               end if
               if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%prompt_spectrum,0.0)) then
                  has_failed = .true.
                  error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%prompt_spectrum,0.0)'
                  return
               end if
               if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%velocity,0.0)) then
                  has_failed = .true.
                  error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%velocity,0.0)'
                  return
               end if
               if (fnequals(np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%beta,0.0)) then
                  has_failed = .true.
                  error_message = 'fnequals(np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%beta,0.0)'
                  return
               end if


               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%total_set) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%total_set'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%absorption_set) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%absorption_set'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%scatter_set) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%scatter_set'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%removal_set) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%removal_set'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%transport_set(1)) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%transport_set(1)'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%transport_set(2)) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%transport_set(2)'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%transport_set(3)) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%transport_set(3)'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%diffusion_set(1)) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%diffusion_set(1)'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%diffusion_set(3)) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%diffusion_set(3)'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%diffusion_set(3)) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%diffusion_set(3)'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%fission_set) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%fission_set'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%production_set) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%production_set'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%power_set) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%power_set'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission_set) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%energy_released_per_fission_set'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%np_released_per_fission_set) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%np_released_per_fission_set'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%prompt_spectrum_set) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%prompt_spectrum_set'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%velocity_set) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%velocity_set'
                  return
               end if
               if (np_radmat_create%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%beta_set) then
                  has_failed = .true.
                  error_message = 'np_radmat_create%dataset_radmats('//int2str(dmat)//')%physical_radmats('//int2str(pmat)//')'//& 
                                  &'%radmats('//int2str(rmat)//')%beta_set'
                  return
               end if            
            end do       
         end do      
      end do 
      
      if (size(np_radmat_create%delayed_lambda_spectrum%lambda) /= number_of_delayed_groups_expected) then
         has_failed = .true.
         error_message = 'size(np_radmat_create%delayed_lambda_spectrum%lambda) /= number_of_delayed_groups_expected'
         return
      end if

      if (size(np_radmat_create%delayed_lambda_spectrum%spectrum,1) /= number_of_energy_groups_expected) then
         has_failed = .true.
         error_message = 'size(np_radmat_create%delayed_lambda_spectrum%spectrum,1) /= number_of_energy_groups_expected'
         return
      end if

      if (size(np_radmat_create%delayed_lambda_spectrum%spectrum,2) /= number_of_delayed_groups_expected) then
         has_failed = .true.
         error_message = 'size(np_radmat_create%delayed_lambda_spectrum%spectrum,2) /= number_of_delayed_groups_expected'
         return
      end if

      if (fnequals(np_radmat_create%delayed_lambda_spectrum%lambda,0.0)) then
         has_failed = .true.
         error_message = 'fnequals(np_radmat_create%delayed_lambda_spectrum%lambda,0.0)'
         return
      end if

      if (fnequals(np_radmat_create%delayed_lambda_spectrum%spectrum,0.0)) then
         has_failed = .true.
         error_message = 'fnequals(np_radmat_create%delayed_lambda_spectrum%spectrum,0.0)'
         return
      end if

      if (np_radmat_create%delayed_lambda_spectrum%lambda_set) then
         has_failed = .true.
         error_message = 'np_radmat_create%delayed_lambda_spectrum%lambda_set'
         return
      end if

      if (np_radmat_create%delayed_lambda_spectrum%spectrum_set) then
         has_failed = .true.
         error_message = 'np_radmat_create%delayed_lambda_spectrum%lambda_set'
         return
      end if

      if (np_radmat_create%readin) then
         has_failed = .true.
         error_message = 'np_radmat_create%readin'
         return
      end if
            
      if (.not. np_radmat_create%created) then
         has_failed = .true.
         error_message = '.not. np_radmat_create%created'
         return
      end if
            
      ! remove the options
      call clear_options
      
      call destroy(np_radmat_create)
          
   end subroutine test_np_radmat_create_from_options
   
   ! --------------------------------------------------------------------------

end subroutine test_radiation_materials_create
