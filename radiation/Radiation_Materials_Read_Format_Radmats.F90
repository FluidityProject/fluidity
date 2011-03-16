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

module radiation_materials_read_format_radmats
   
   use futils
   
   use radiation_materials_data_types
   use radiation_materials_read_format_radmats_base
   
   implicit none
   
   private
   
   public :: read_format_radmats

contains

   ! --------------------------------------------------------------------------

   subroutine read_format_radmats(format_radmats_filename, &
                                  dataset_radmat, &
                                  delayed_lambda_spectrum, &
                                  read_delayed_lambda_spectrum, &
                                  read_velocity_data, &
                                  read_power_data, &
                                  problem_dimension, &
                                  record_len)
      
      !!< Read a format_radmats material data set file, which also may include the data objects DELAYED, POWER and VELOCITY
      !!< Due to the formatting of the format_radmats file associated with the scatter data this will only work with
      !!< number_of_energy_groups < 10,000,000. Note that we may read in either SIGTRAN or DIFFUSIONX (Y,Z) numbers within a material
      !!< but not both. We may also read in either SIGTOT or SIGABS within a material but not both. We may also read in either NU or PRODUCTION 
      !!< within a material but not both. Note that the problem dimension is taken into account so that for 1d or 2d only X and X,Y SIGTRAN or DIFFUSION
      !!< terms are searched for within the file. If no SIGTRAN or DIFFUSION is found then the total is used to form SIGTRAN and DIFFUSION.

      character(len=*), intent(in) :: format_radmats_filename      
      type(dataset_radmat_type), intent(inout) :: dataset_radmat
      type(delayed_lambda_spectrum_type), intent(inout) :: delayed_lambda_spectrum  
      logical, intent(in) :: read_delayed_lambda_spectrum
      logical, intent(in) :: read_velocity_data
      logical, intent(in) :: read_power_data
      integer, intent(in) :: problem_dimension
      integer, intent(in) :: record_len
          
      ! local variables 
      integer :: pmat,rmat,idim,iostat
      integer :: format_radmats_file_unit
      integer :: number_of_energy_groups
      integer :: number_of_scatter_moments
      integer :: number_of_delayed_groups      
      integer :: number_of_energy_groups_format_radmats
      integer :: number_of_delayed_groups_format_radmats
      integer :: number_of_scatter_moments_format_radmats
      integer :: number_non_thermal_groups_format_radmats
      integer :: total_number_of_radmats_format_radmats  
      integer :: total_number_of_radmats
      integer :: line_number
      integer :: line_number_macro
      integer :: first_keyword_found
      logical :: format_radmats_file_exists,exit_if_eof
      integer, dimension(:), allocatable :: keyword_find ! the specific keyword integer ids to find
      character(len=record_len) :: line_string
      ! the keyword list to use to search through the format_radmats input file, the 21 is due to the scatter moment data
      character(len=21), dimension(:), allocatable :: keyword_list

      ewrite(1,*) 'Read in format_radmats data with max keyword record length ',record_len
                  
      ! deduce the number of energy groups from the size of a totals array
      number_of_energy_groups = size(dataset_radmat%physical_radmats(1)%radmats(1)%total)
      
      ! deduce the number of scatter moments from the size of a scatter array - note that this is dataset dependent
      number_of_scatter_moments = size(dataset_radmat%physical_radmats(1)%radmats(1)%scatter,3)
      
      ! deduce the number of delayed groups from the size of a beta array
      number_of_delayed_groups = size(dataset_radmat%physical_radmats(1)%radmats(1)%beta)

      ! deduce the total number of radmats associated with this dataset from the options
      total_number_of_radmats = 0

      total_material_loop: do pmat = 1,size(dataset_radmat%physical_radmats)            
  
         total_number_of_radmats = total_number_of_radmats + size(dataset_radmat%physical_radmats(pmat)%radmats)
      
      end do total_material_loop
      
      ! initiliase the base module via creating the format_radmats keyword list
      call keyword_list_initialise_format_radmats(keyword_list)
            
      ! get a free io unit for format_radmats file
      format_radmats_file_unit = free_unit()
      
      inquire(file=trim(format_radmats_filename),exist=format_radmats_file_exists)
      
      no_format_radmats_file: if (.not. format_radmats_file_exists) then
         
         ewrite(-1,*) "Input error: format_radmats file ",trim(format_radmats_filename)
         FLExit("was not found")
         
      end if no_format_radmats_file
      
      ! open a sequential access file
      open(unit   = format_radmats_file_unit, &
           file   = trim(format_radmats_filename), &
           status = 'old', &
           form   = 'formatted', &
           action = 'read', &
           access = 'sequential', &
           iostat = iostat)
      
      check_open_fine: if (iostat /= 0) then
      
         ewrite(-1,*) "Error opening format_radmats file ",trim(format_radmats_filename)
         FLAbort("IOSTAT returned non zero")
      
      end if check_open_fine
      
      ! calculate the number of radmats, energy groups, number of scatter moments, number of delayed groups, number of format_radmats non thermal groups 
      call deduce_number_of_items_format_radmats(format_radmats_file_unit, &
                                                 total_number_of_radmats_format_radmats, &
                                                 number_of_energy_groups_format_radmats, &
                                                 number_non_thermal_groups_format_radmats, &
                                                 number_of_scatter_moments_format_radmats, &
                                                 number_of_delayed_groups_format_radmats, &
                                                 record_len, &
                                                 keyword_list)

      check_material_consistency: if (total_number_of_radmats /= total_number_of_radmats_format_radmats) then
      
         ewrite(-1,*) "Input error: format_radmats file ",trim(format_radmats_filename)," has ",total_number_of_radmats_format_radmats,&
                     &" total materials, whereas the simulation data base has ",total_number_of_radmats
         FLExit("Total number of materials inconsistency")
      
      end if check_material_consistency

      check_energy_group_consistency: if (number_of_energy_groups /= number_of_energy_groups_format_radmats) then
      
         ewrite(-1,*) "Input error: format_radmats file ",trim(format_radmats_filename)," has ",number_of_energy_groups_format_radmats,&
                     &" energy groups, whereas the simulation has ",number_of_energy_groups
         FLExit("Energy group inconsistency")
      
      end if check_energy_group_consistency

      check_scatter_moment_consistency: if (number_of_scatter_moments /= number_of_scatter_moments_format_radmats) then
      
         ewrite(-1,*) "Input error: format_radmats file ",trim(format_radmats_filename)," has ",number_of_scatter_moments_format_radmats,&
                     &" scatter moments, whereas the simulation has ",number_of_scatter_moments
         FLExit("Scatter moment inconsistency")
      
      end if check_scatter_moment_consistency
      
      check_delayed_group_consistency: if (number_of_delayed_groups /= number_of_delayed_groups_format_radmats) then
      
         ewrite(-1,*) "Input error: format_radmats file ",trim(format_radmats_filename)," has ",number_of_delayed_groups_format_radmats,&
                     &" energy groups, whereas the simulation has ",number_of_delayed_groups
         FLExit("Delayed group inconsistency")
      
      end if check_delayed_group_consistency
                  
      ! read in each radmat within each physical material - not the delayed, velocity or power yet, they are done after

      ! go back to start of file
      call rewind_file_seq(format_radmats_file_unit, &
                           line_number)

      physical_material_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
      
         radmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)
                        
            ! find the next MACRO keyword line
            allocate(keyword_find(1))
            keyword_find(1) = 1
            exit_if_eof = .true.
            call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                                    line_number, &
                                                    first_keyword_found, &
                                                    exit_if_eof, &
                                                    keyword_find, &
                                                    keyword_list, &
                                                    line_string)

            deallocate(keyword_find)
                                    
            ! keep the line number associated with this MACRO line
            line_number_macro = line_number
                                    
            ! read the scatter momemnts of this macro - then put the file back to the start of this macro
            call read_scatter_mom_format_radmats(format_radmats_file_unit, &
                                                 line_number, &
                                                 dataset_radmat%physical_radmats(pmat)%radmats(rmat), &
                                                 record_len, &
                                                 keyword_list)
                        
            ! read in either the absorption or total and form the other using the scatter data already found - then put the file back to the start of this macro 
            call read_absorption_and_total_format_radmats(format_radmats_file_unit, &
                                                          line_number, &
                                                          dataset_radmat%physical_radmats(pmat)%radmats(rmat), &
                                                          record_len, &
                                                          keyword_list)
                                                                 
            ! calc the removal moments from the total and scatter data - then put the file back to the start of this macro
            call calculate_removal_moments_format_radmats(format_radmats_file_unit, &
                                                          line_number, &
                                                          dataset_radmat%physical_radmats(pmat)%radmats(rmat), &
                                                          record_len, &
                                                          keyword_list)
            
            ! read or associate the sigtran and diffusion cross sections for each geometric dimension - then put the file back to the start of this macro
            geom_dim_loop: do idim = 1,problem_dimension 
            
               call read_or_associate_transport_and_diffusion_format_radmats(format_radmats_file_unit, &
                                                                             line_number, &
                                                                             dataset_radmat%physical_radmats(pmat)%radmats(rmat), &
                                                                             idim, &
                                                                             record_len, &
                                                                             keyword_list)
                                                                                    
            end do geom_dim_loop
                                                
            ! read SIGFISS, FISSNU, PRODUCTION and FISSCHI  - then put the file back to the start of this macro               
            call read_fission_data_format_radmats(format_radmats_file_unit, &
                                                  line_number, &
                                                  dataset_radmat%physical_radmats(pmat)%radmats(rmat), &
                                                  number_non_thermal_groups_format_radmats, &
                                                  record_len, &
                                                  keyword_list)

            ! go back to the line number of the current macro
            call go_to_line_seq(format_radmats_file_unit, &
                                line_number, &
                                go_to_line_number = line_number_macro, &
                                exit_if_eor       = .true.,&
                                exit_if_eof       = .true.) 
                                                               
         end do radmat_loop
         
      end do physical_material_loop
      
      ! Read in the VELOCITY, DELAYED and POWER/ENERGY_RELEASED data if needed
            
      ! Read in the VELOCITY which may not be material dependent if needed
      need_velocity: if (read_velocity_data) then

         allocate(keyword_find(1))
         keyword_find(1) = 14
         call read_velocity_style_data_format_radmats(format_radmats_file_unit, &
                                                      dataset_radmat, &
                                                      number_of_energy_groups, &
                                                      keyword_find, &
                                                      record_len, &
                                                      keyword_list)

         deallocate(keyword_find)
                  
      end if need_velocity
      
      ! read in the power (or energy released per fission) data needed for this material data set if needed
      need_power: if (read_power_data) then
      
         call read_power_data_format_radmats(format_radmats_file_unit, &
                                             dataset_radmat, &
                                             record_len, &
                                             keyword_list)
              
      end if need_power
      
      ! Read in the Delayed data needed from this material data set if needed if needed
      need_delayed_lambda_spectrum: if (read_delayed_lambda_spectrum) then
               
         call read_delayed_data_format_radmats(format_radmats_file_unit, &
                                               dataset_radmat, &
                                               record_len, &
                                               keyword_list, &
                                               delayed_lambda_spectrum, &
                                               read_delayed_lambda_spectrum)
         
      end if need_delayed_lambda_spectrum
      
      ! cleanup the base via deallocating the keyword list
      call keyword_list_cleanup(keyword_list)
                                 
      close(format_radmats_file_unit)
            
   end subroutine read_format_radmats

   ! --------------------------------------------------------------------------

end module radiation_materials_read_format_radmats

