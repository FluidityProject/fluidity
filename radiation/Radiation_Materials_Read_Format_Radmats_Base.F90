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

module radiation_materials_read_format_radmats_base
   
   !!< Base procedures associated with reading in a format_radmats dataset file into stored data type arrays
   
   use futils
   
   use radiation_materials_data_types
   
   implicit none
   
   private
   
   public :: keyword_list_initialise_format_radmats, &
             keyword_list_cleanup, &
             deduce_number_of_items_format_radmats, &
             read_scatter_mom_format_radmats, &
             read_absorption_and_total_format_radmats, &
             calculate_removal_moments_format_radmats, &
             read_or_form_trans_and_diff_format_radmats, &
             read_fission_data_format_radmats, &
             read_velocity_style_data_format_radmats, &
             read_delayed_data_format_radmats, &
             read_power_data_format_radmats, &
             calculate_power_format_radmats, &
             read_format_radmats_siga_style_xsection, &
             read_format_radmats_sigs_style_xsection, &
             count_number_lines_with_keyword, &
             find_line_with_any_keyword, &
             find_line_with_any_desired_keyword, &
             all_upper_case, &
             make_character_upper_case, &
             read_words_from_string, &
             string_word_count, &
             number_substrings_within_string, &
             read_next_line_seq, &
             read_previous_line_seq, &
             rewind_file_seq, &
             go_to_line_seq, &
             desired_keyword_not_first_found

   integer :: iostat,status
   
contains

   ! --------------------------------------------------------------------------
   
   subroutine keyword_list_initialise_format_radmats(keyword_list)
      
      !!< Form the keyword list used when searching the format_radmats file

      character(len=*), dimension(:), allocatable, intent(inout) :: keyword_list     
      
      check_allocated: if (allocated(keyword_list)) then
      
         FLAbort("Cannot allocate already allocated keyword_list in module radiation_materials_read_format_radmats_base")
      
      else check_allocated
         
         allocate(keyword_list(22),STAT=status)
               
         if (status /= 0) FLAbort("Issue allocating memory for keyword_list in module radiation_materials_read_format_radmats_base")
         
         keyword_list(1)  = 'MACRO'
         keyword_list(2)  = 'SIGABS'
         keyword_list(3)  = 'SIGTOT'
         keyword_list(4)  = 'SIGTRAN ' ! the space at the end is important
         keyword_list(5)  = 'SIGTRANY'
         keyword_list(6)  = 'SIGTRANZ'
         keyword_list(7)  = 'DIFFUSIONX'
         keyword_list(8)  = 'DIFFUSIONY' 
         keyword_list(9)  = 'DIFFUSIONZ' 
         keyword_list(10) = 'SIGFISS' 
         keyword_list(11) = 'FISSNU' 
         keyword_list(12) = 'PRODUCTION' 
         keyword_list(13) = 'FISSCHI' 
         keyword_list(14) = 'VELOCITY' 
         keyword_list(15) = 'POWER' 
         keyword_list(16) = 'ENERGY_RELEASED' 
         keyword_list(17) = 'BETA' 
         keyword_list(18) = 'LAMBDA' 
         keyword_list(19) = 'DELAYED_SPECTRUM' 
         keyword_list(20) = 'GROUP' 
         keyword_list(21) = 'DELAYED_GROUPS' 
         keyword_list(22) = 'SCATTER'
                     
      end if check_allocated
       
   end subroutine keyword_list_initialise_format_radmats
   
   ! --------------------------------------------------------------------------
   
   subroutine keyword_list_cleanup(keyword_list)
   
      !!< Cleanup the keyword list that was used when searching the format_radmats file

      character(len=*), dimension(:), allocatable, intent(inout)  :: keyword_list  
            
      if (allocated(keyword_list)) deallocate(keyword_list)
   
   end subroutine keyword_list_cleanup
   
   ! --------------------------------------------------------------------------
   
   subroutine deduce_number_of_items_format_radmats(format_radmats_file_unit, &
                                                    total_number_of_radmats_format_radmats, &
                                                    number_of_energy_groups_format_radmats, &
                                                    number_non_thermal_groups_format_radmats, &
                                                    number_of_scatter_moments_format_radmats, &
                                                    number_of_delayed_groups_format_radmats, &
                                                    record_len, &
                                                    keyword_list)
                                                  
      !!< Deduce the number of MACROS then from the first MACRO block, the number of GROUP statements associated with the first SCATTER how many groups there are in 
      !!< format_radmats file with io unit format_radmats_file_unit. The number of scatter moments is deduced via counting of keyword SCATTER in first MACRO.
      !!< Then following the first MACRO find, there is a search for the first FISSCHI (which may not be there or in the first MACRO) 
      !!< to deduce the number_non_thermal_groups_format_radmats - which counts lines between it and the next keyword
      !!< A search for DELAYED_GROUPS is then used to determine the number of delayed groups in the format_radmats file 
      
      integer, intent(in) ::format_radmats_file_unit 
      integer, intent(out) :: total_number_of_radmats_format_radmats
      integer, intent(out) :: number_of_energy_groups_format_radmats
      integer, intent(out) :: number_non_thermal_groups_format_radmats
      integer, intent(out) :: number_of_scatter_moments_format_radmats
      integer, intent(out) :: number_of_delayed_groups_format_radmats
      integer, intent(in) :: record_len
      character(len=*), dimension(:), intent(in) :: keyword_list
            
      ! local variables
      logical :: exit_if_eof
      integer :: line
      integer :: number_of_words
      integer :: first_keyword_found
      integer :: line_number,line_number_first_macro,keyword_to_count,keyword_to_stop
      integer :: line_number_difference,line_number_first_fisschi         
      integer, dimension(:), allocatable :: keyword_find ! the specific keyword integer ids to find
      character(len=record_len) :: line_string
      logical :: end_of_file
           
      ! read from start of file - initialise line number 
      line_number                              = 0
      line_number_first_fisschi                = 0
      line_number_difference                   = 0
      total_number_of_radmats_format_radmats   = 0
      number_of_energy_groups_format_radmats   = 0
      number_non_thermal_groups_format_radmats = 0
      number_of_scatter_moments_format_radmats = 0
      number_of_delayed_groups_format_radmats  = 0

      ! go back to start of file
      call rewind_file_seq(format_radmats_file_unit, &
                           line_number)
      
      ! first count the number of macros - being keyword 1 in the keyword_list
      keyword_to_count = 1
      keyword_to_stop = 0
      call count_number_lines_with_keyword(format_radmats_file_unit, &
                                           line_number, &
                                           total_number_of_radmats_format_radmats, &
                                           keyword_to_count, &
                                           record_len, &
                                           keyword_list)

      ! go back to start of file
      call rewind_file_seq(format_radmats_file_unit, &
                           line_number)
      
      ! find the first MACRO keyword line
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
      
      line_number_first_macro = line_number
                        
      ! find the first 'SCATTER' keyword line 
      allocate(keyword_find(1))
      keyword_find(1) = 22
      exit_if_eof = .true.
      call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                              line_number, &
                                              first_keyword_found, &
                                              exit_if_eof, &
                                              keyword_find, &
                                              keyword_list, &
                                              line_string)

      deallocate(keyword_find)
            
      ! now count the number of lines found with the keyword GROUP before either any other keyword or EoF
      keyword_to_count = 20
      keyword_to_stop = 0
      call count_number_lines_with_keyword(format_radmats_file_unit, &
                                           line_number, &
                                           number_of_energy_groups_format_radmats, &
                                           keyword_to_count, &
                                           record_len, &
                                           keyword_list, &
                                           keyword_to_stop=keyword_to_stop)                                                  
                       
      ! go back to the line number of the first macro
      call go_to_line_seq(format_radmats_file_unit, &
                          line_number, &
                          go_to_line_number=line_number_first_macro, &
                          exit_if_eor=.true.,&
                          exit_if_eof=.true.) 
            
      ! Deduce the number of scatter moment via now count number of lines with 
      ! SCATTER in first MACRO (stop if find next MACRO)
      keyword_to_count = 22
      keyword_to_stop = 1
      call count_number_lines_with_keyword(format_radmats_file_unit, &
                                           line_number, &
                                           number_of_scatter_moments_format_radmats, &
                                           keyword_to_count, &
                                           record_len, &
                                           keyword_list, &
                                           keyword_to_stop=keyword_to_stop)                                                  
                       
      ! go back to the line number of the first macro
      call go_to_line_seq(format_radmats_file_unit, &
                          line_number, &
                          go_to_line_number=line_number_first_macro, &
                          exit_if_eor=.true.,&
                          exit_if_eof=.true.) 
      
      ! now search for FISSCHI
      allocate(keyword_find(1))
      keyword_find(1) = 13
      exit_if_eof = .false.
      call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                              line_number, &
                                              first_keyword_found, &
                                              exit_if_eof, &
                                              keyword_find, &
                                              keyword_list, &
                                              line_string)

      deallocate(keyword_find)
      
      have_fisschi: if (first_keyword_found == 1) then
         
         ! save the FISSCHI line
         line_number_first_fisschi = line_number
         
         ! find the next keyword or end of file
         exit_if_eof = .false.
         call find_line_with_any_keyword(format_radmats_file_unit, &
                                         line_number, &
                                         first_keyword_found, &
                                         exit_if_eof, &
                                         keyword_list, &
                                         line_string) 
                           
         line_number_difference = line_number - line_number_first_fisschi
         
         ! check that the line number difference makes sense 
         check_line_number_diff_again: if (line_number_difference < 2) then
      
            FLExit("Error reading format_radmats formatted file as line difference between a FISSCHI and next keyword line is less than 2")
      
         end if check_line_number_diff_again
      
         ! put the file back to the FISSCHI line
         call go_to_line_seq(format_radmats_file_unit, &
                             line_number, &
                             go_to_line_number = line_number_first_fisschi, &
                             exit_if_eor       = .true.,&
                             exit_if_eof       = .true.) 

         ! now read in each line with FISSCHI data on and count the number of words (reals) 
         
         ! initialise the count
         number_non_thermal_groups_format_radmats = 0
         
         fisschi_line_loop: do line = 1,line_number_difference - 1
         
            ! read the next line
            call read_next_line_seq(format_radmats_file_unit, &
                                    line_string, &
                                    line_number, &
                                    end_of_file, &
                                    exit_if_eor=.true., &
                                    exit_if_eof=.true.)
            
            ! find the number of words within the line_string
            call string_word_count(line_string, &
                                   number_of_words)
                        
            number_non_thermal_groups_format_radmats = &
            number_non_thermal_groups_format_radmats + &
            number_of_words
            
         end do fisschi_line_loop
         
      end if have_fisschi
      
      ! go back to start of file
      call rewind_file_seq(format_radmats_file_unit, &
                           line_number)
      
      ! now find the number of delayed groups by searching for DELAYED_GROUPS then reading a value of the next line
      allocate(keyword_find(1))
      keyword_find(1) = 21
      exit_if_eof = .false.
      call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                              line_number, &
                                              first_keyword_found, &
                                              exit_if_eof, &
                                              keyword_find, &
                                              keyword_list, &
                                              line_string)
                                              
      deallocate(keyword_find)
      
      have_delayed_spectrum: if (first_keyword_found == 1) then
         
         ! read the next line
         call read_next_line_seq(format_radmats_file_unit, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor=.true., &
                                 exit_if_eof=.true.)
         
         ! read the value of the line_string
         read(line_string,*) number_of_delayed_groups_format_radmats
      
      end if have_delayed_spectrum

      ! go back to start of file
      call rewind_file_seq(format_radmats_file_unit, &
                           line_number)
         
   end subroutine deduce_number_of_items_format_radmats

   ! --------------------------------------------------------------------------

   subroutine read_scatter_mom_format_radmats(format_radmats_file_unit, &
                                              line_number, &
                                              radmat, &
                                              record_len, &
                                              keyword_list)
            
      !!< read in the format_radmats style scatter moments for this macro tagged from line_number
      
      integer, intent(in) :: format_radmats_file_unit
      integer, intent(inout) :: line_number
      type(radmat_type), intent(inout) :: radmat 
      integer, intent(in) :: record_len
      character(len=*), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      logical :: exit_if_eof
      logical :: found_correct_scatter
      integer :: m,m_check
      integer :: number_of_scatter_moments
      integer :: line_number_macro
      integer :: first_keyword_found
      integer :: desired_keyword
      integer :: m_read 
      integer, dimension(:), allocatable :: keyword_find ! the specific keyword integer ids to find      
      character(len=record_len) :: line_string   
      character(len=record_len) :: cdummy,c_m_read   

      ! find the number of scatter moments from the size of scatter
      number_of_scatter_moments = size(radmat%scatter,3)
      
      ! save the line number of this macro so as to return to it
      line_number_macro = line_number
      
      allocate(keyword_find(2))    
      exit_if_eof = .true.      
      
      ! search for MACRO as well so as to not search to far each time
      keyword_find(1) = 1

      read_scatter_mom_loop: do m = 1,number_of_scatter_moments
                                 
         keyword_find(2) = 22
         
         ! initialise a flag of whether or not the correct SCATTER is found
         found_correct_scatter = .false.
         
         ! check each line sequentially that has the keyword SCATTER to if this is the moment we want
         check_scatter_loop: do m_check = 1,number_of_scatter_moments
         
            call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                                    line_number, &
                                                    first_keyword_found, &
                                                    exit_if_eof, &
                                                    keyword_find, &
                                                    keyword_list, &
                                                    line_string)          
          
            ! exit if not found SCATTER as this is always needed
            desired_keyword = 2
            desired_if: if (first_keyword_found /= desired_keyword) then 
            
               call desired_keyword_not_first_found(keyword_find, &
                                                    desired_keyword, &
                                                    first_keyword_found, &
                                                    keyword_list)
         
            end if desired_if
            
            ! find which moment in the file this is
            read(line_string,*) cdummy,cdummy,c_m_read
            
            read(c_m_read,*) m_read
       
            ! exit if this is the correct moment to read
            found: if (m == m_read) then
               
               found_correct_scatter = .true.
               
               exit check_scatter_loop
            
            end if found
            
         end do check_scatter_loop
         
         ! exit if not found the correct SCATTER
         not_found: if (.not. found_correct_scatter) then
            
            ewrite(1,*) 'Error searching for SCATTER MOMENT block ',m
            FLExit('Error reading radiation radmats file as SCATTER MOMENT block not found')
         
         end if not_found
            
         ! now read in the SCATTER MOMENT 
         call read_format_radmats_sigs_style_xsection(radmat%scatter(:,:,m), &
                                                      line_number, &
                                                      format_radmats_file_unit, &
                                                      record_len, &
                                                      keyword_list)
         
         ! put the file back to the start of this macro
         call go_to_line_seq(format_radmats_file_unit, &
                             line_number, &
                             go_to_line_number = line_number_macro, &
                             exit_if_eor       = .true.,&
                             exit_if_eof       = .true.) 
              
      end do read_scatter_mom_loop
      
      deallocate(keyword_find)
      
      ! set tag to say scatter is now set in for this radmat
      radmat%scatter_set = .true.
     
   end subroutine read_scatter_mom_format_radmats
   
   ! --------------------------------------------------------------------------

   subroutine read_absorption_and_total_format_radmats(format_radmats_file_unit, &
                                                       line_number, &
                                                       radmat, &
                                                       record_len, &
                                                       keyword_list)
            
      !!< read in either the absorption or total and form the other using the scatter data already found from a format_radmats file for
      !!< radmat. This relies on the scatter data being set, if it isnt then read the scatter data.
      
      integer, intent(in) :: format_radmats_file_unit
      integer, intent(inout) :: line_number
      type(radmat_type), intent(inout) :: radmat 
      integer, intent(in) :: record_len
      character(len=*), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      logical :: exit_if_eof
      integer :: g,number_of_energy_groups,line_number_macro,first_keyword_found
      integer, dimension(:), allocatable :: keyword_find ! the specific keyword integer ids to find      
      character(len=record_len) :: line_string 
                   
      ! find the number_of_energy_groups from the size of total
      number_of_energy_groups = size(radmat%total)

      ! save the line number of this macro so as to return to it
      line_number_macro = line_number
      
      scatter_not_set: if (.not. radmat%scatter_set) then
   
         ! read the scatter momemnts and set the line number back to the MACRO
         call read_scatter_mom_format_radmats(format_radmats_file_unit, &
                                              line_number, &
                                              radmat, &
                                              record_len, &
                                              keyword_list)
      
      end if scatter_not_set
      
      allocate(keyword_find(2))    
      exit_if_eof = .false.      
      
      ! search for MACRO as well so as to not search to far each time
      keyword_find(1) = 1
      
      ! find sigabs if there
      keyword_find(2) = 2

      call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                              line_number, &
                                              first_keyword_found, &
                                              exit_if_eof, &
                                              keyword_find, &
                                              keyword_list, &
                                              line_string)  

      ! if not found SIGABS then search for SIGTOT, if neither (or both) then exit
      not_found_siga: if (first_keyword_found /= 2) then 
               
         ! search for SIGTOT
         keyword_find(2) = 3
      
         ! put the file back to the start of this macro
         call go_to_line_seq(format_radmats_file_unit, &
                             line_number, &
                             go_to_line_number = line_number_macro, &
                             exit_if_eor       = .true.,&
                             exit_if_eof       = .true.) 
                     
         call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                                 line_number, &
                                                 first_keyword_found, &
                                                 exit_if_eof, &
                                                 keyword_find, &
                                                 keyword_list, &
                                                 line_string)    
               
         ! exit if neither SIGABS or SIGTOT found else read in sigtot and calc sigabs
         not_found_sigtot_or_siga: if (first_keyword_found /= 2) then
               
            ewrite(-1,*) "Error reading format_radmats file"
            FLExit("Neither the SIGABS or SIGTOT cross sections found")
               
         else not_found_sigtot_or_siga
               
            ! read in the SIGTOT 
            call read_format_radmats_siga_style_xsection(radmat%total, &
                                                         line_number, &
                                                         format_radmats_file_unit, &
                                                         record_len)
                  
            ! calc the sigabs from the sigtot and sigs moment 1
            siga_group_loop: do g = 1,number_of_energy_groups 
                  
               radmat%absorption(g) = radmat%total(g) - sum(radmat%scatter(g,:,1))
                  
            end do siga_group_loop
                  
         end if not_found_sigtot_or_siga
            
      else not_found_siga
                        
         ! read in the SIGABS 
         call read_format_radmats_siga_style_xsection(radmat%absorption, &
                                                      line_number, &
                                                      format_radmats_file_unit, &
                                                      record_len)
               
         ! check if there is also sigtot and if so exit as cannot have both, else form the sigtot
         keyword_find(2) = 3

         ! put the file back to the start of this macro
         call go_to_line_seq(format_radmats_file_unit, &
                             line_number, &
                             go_to_line_number = line_number_macro, &
                             exit_if_eor       = .true.,&
                             exit_if_eof       = .true.)
                              
         call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                                 line_number, &
                                                 first_keyword_found, &
                                                 exit_if_eof, &
                                                 keyword_find, &
                                                 keyword_list, &
                                                 line_string)    
                              
         sigtot_and_siga: if (first_keyword_found == 2) then
               
            ewrite(-1,*) "Error reading format_radmats file"
            FLExit("Both the SIGABS and SIGTOT cross sections found")
                              
         else sigtot_and_siga
               
            ! calc the sigtot from sigabs and sigs moment 1
            sigtot_group_loop: do g = 1,number_of_energy_groups 
                  
               radmat%total(g) = radmat%absorption(g) + sum(radmat%scatter(g,:,1))
                  
            end do sigtot_group_loop
                                 
         end if sigtot_and_siga
               
      end if not_found_siga

      deallocate(keyword_find)
            
      ! set tag to say absorption and total are now set in for this radmat
      radmat%absorption_set = .true.
      radmat%total_set = .true.

      ! put the file back to the start of this macro
      call go_to_line_seq(format_radmats_file_unit, &
                          line_number, &
                          go_to_line_number = line_number_macro, &
                          exit_if_eor       = .true.,&
                          exit_if_eof       = .true.) 
      
   end subroutine read_absorption_and_total_format_radmats

   ! --------------------------------------------------------------------------

   subroutine calculate_removal_moments_format_radmats(format_radmats_file_unit, &
                                                       line_number, &
                                                       radmat, &
                                                       record_len, &
                                                       keyword_list)
            
      !!< Calculate the removal moments from the total and scatter data for this radmat. If the total is not set then 
      !!< make sure it is (which in itself make sures that the scatter data is set)
      
      integer, intent(in) :: format_radmats_file_unit
      integer, intent(inout) :: line_number
      type(radmat_type), intent(inout) :: radmat 
      integer, intent(in) :: record_len
      character(len=*), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      integer :: g,m,number_of_energy_groups,number_of_scatter_moments
      
      ! find the number of groups and scatter moments from size of scatter
      number_of_energy_groups = size(radmat%scatter,1)
      number_of_scatter_moments = size(radmat%scatter,3)
      
      total_not_set: if (.not. radmat%total_set) then

         ! read in either the absorption or total and form the other using the scatter data already found 
         call read_absorption_and_total_format_radmats(format_radmats_file_unit, &
                                                       line_number, &
                                                       radmat, &
                                                       record_len, &
                                                       keyword_list)
      
      end if total_not_set

      ! now form the removal cross section for each moment
      removal_form_g_loop: do g = 1,number_of_energy_groups
               
         removal_form_m_loop: do m = 1,number_of_scatter_moments
               
            radmat%removal(g,m) = radmat%total(g) - radmat%scatter(g,g,m)
               
         end do removal_form_m_loop
            
      end do removal_form_g_loop

      ! set tag to say removal is now set in for this radmat
      radmat%removal_set = .true.
      
   end subroutine calculate_removal_moments_format_radmats

   ! --------------------------------------------------------------------------

   subroutine read_or_form_trans_and_diff_format_radmats(format_radmats_file_unit, &
                                                         line_number, &
                                                         radmat, &
                                                         idim, &
                                                         record_len, &
                                                         keyword_list)
            
      !!< Read from a format_radmats file or form the sigtran and diffusion cross sections for this geometric dimension idim. 
      !!< This could depend on the total already being set so check that it is and if not set it. For the second and third dimension if 
      !!< neither transport or diffusion found then set to the first dimension values.
      
      integer, intent(in) :: format_radmats_file_unit
      integer, intent(inout) :: line_number
      type(radmat_type), intent(inout) :: radmat 
      integer, intent(in) :: idim
      integer, intent(in) :: record_len
      character(len=*), dimension(:), intent(in) :: keyword_list
      
      ! local variables 
      logical :: exit_if_eof
      integer :: line_number_macro,first_keyword_found
      integer, dimension(:), allocatable :: keyword_find ! the specific keyword integer ids to find       
      character(len=record_len) :: line_string  

      ! save the line number of this macro so as to return to it
      line_number_macro = line_number
             
      total_not_set: if (.not. radmat%total_set) then
      
         ! read in either the absorption or total and form the other using the scatter data already found 
         call read_absorption_and_total_format_radmats(format_radmats_file_unit, &
                                                       line_number, &
                                                       radmat, &
                                                       record_len, &
                                                       keyword_list)
      
      end if total_not_set

      allocate(keyword_find(2))    
      exit_if_eof = .false.      
      
      ! search for MACRO as well so as to not search to far each time
      keyword_find(1) = 1

      ! read SIGTRAN if there 
      which_dim_tran: if (idim == 1) then 
         
         keyword_find(2) = 4
      
      else if (idim == 2) then which_dim_tran

         keyword_find(2) = 5
      
      else if (idim == 3) then which_dim_tran
      
         keyword_find(2) = 6
      
      end if which_dim_tran
      
      call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                              line_number, &
                                              first_keyword_found, &
                                              exit_if_eof, &
                                              keyword_find, &
                                              keyword_list, &
                                              line_string)
            
      ! if not found check for diffusion, if neither found set from previous dim or total. 
      ! If found sigtran check if also diffusion as cannot have both
      not_found_sigtran_if: if (first_keyword_found /= 2) then 
       
         ! read DIFFUSION if there - else set to = 1/3sigtranx or vice versa if there
         which_dim_diff: if (idim == 1) then 
         
            keyword_find(2) = 7
      
         else if (idim == 2) then which_dim_diff

            keyword_find(2) = 8
      
         else if (idim == 3) then which_dim_diff
      
            keyword_find(2) = 9
      
         end if which_dim_diff         
         
         ! put the file back to the start of this macro
         call go_to_line_seq(format_radmats_file_unit, &
                             line_number, &
                             go_to_line_number=line_number_macro, &
                             exit_if_eor=.true.,&
                             exit_if_eof=.true.) 
                             
         call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                                 line_number, &
                                                 first_keyword_found, &
                                                 exit_if_eof, &
                                                 keyword_find, &
                                                 keyword_list, &
                                                 line_string)            
               
         ! if not found sigtran or diffusion then set from sigtot for idim 1 or to idim 1 values for the other two idim 
         ! else read diffusion and set sigtran
         no_sigtran_or_diffusion: if (first_keyword_found /= 2) then
                       
            is_this_idim1: if (idim ==1) then
                  
               ! set the transport cross section to the total
               radmat%transport(:,idim) = radmat%total(:)
               
               ! set the diffusion coeff from the sigtran coeff
               radmat%diffusion(:,idim) = 1.0/(3.0*radmat%transport(:,idim))
            
            else is_this_idim1

               ! set the transport cross section to that of idim 1
               radmat%transport(:,idim) = radmat%transport(:,1)
               
               ! set the diffusion coeff from the sigtran coeff
               radmat%diffusion(:,idim) = 1.0/(3.0*radmat%transport(:,idim))
                        
            end if is_this_idim1
                                 
         else no_sigtran_or_diffusion

            ! now read in the DIFFUSION 
            call read_format_radmats_siga_style_xsection(radmat%diffusion(:,idim), &
                                                         line_number, &
                                                         format_radmats_file_unit, &
                                                         record_len)

            ! set the sigtran coeff from the diffusion coeff
            radmat%transport(:,idim) = 1.0/(3.0*radmat%diffusion(:,idim))
                                    
         end if no_sigtran_or_diffusion
               
      else not_found_sigtran_if
                        
         ! now read in the SIGTRAN 
         call read_format_radmats_siga_style_xsection(radmat%transport(:,idim), &
                                                      line_number, &
                                                      format_radmats_file_unit, &
                                                      record_len)
            
        ! check if DIFFUSION there also as cannot have both, else form diffusionx from sigtran
         which_dim_diff_again: if (idim == 1) then 
         
            keyword_find(2) = 7
      
         else if (idim == 2) then which_dim_diff_again

            keyword_find(2) = 8
      
         else if (idim == 3) then which_dim_diff_again
      
            keyword_find(2) = 9
      
         end if which_dim_diff_again    
         
         ! put the file back to the start of this macro
         call go_to_line_seq(format_radmats_file_unit, &
                             line_number, &
                             go_to_line_number=line_number_macro, &
                             exit_if_eor=.true.,&
                             exit_if_eof=.true.) 

         call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                                 line_number, &
                                                 first_keyword_found, &
                                                 exit_if_eof, &
                                                 keyword_find, &
                                                 keyword_list, &
                                                 line_string)             
                                             
         not_found_diffusion_if: if (first_keyword_found /= 2) then 
                                             
            ! set the diffusion coeff from the sigtran coeff
            radmat%diffusion(:,idim) = 1.0/(3.0*radmat%transport(:,idim))
                           
         else not_found_diffusion_if
                        
            ewrite(-1,*) "Error reading format_radmats file"
            FLExit("Both the SIGTRAN and DIFFUSION cross sections found")

         end if not_found_diffusion_if             
            
      end if not_found_sigtran_if 

      deallocate(keyword_find)    
                        
      ! set tag to say transport and diffusion are now set in for this radmat for this dimension idim
      radmat%transport_set(idim) = .true.
      radmat%diffusion_set(idim) = .true.

      ! put the file back to the start of this macro
      call go_to_line_seq(format_radmats_file_unit, &
                          line_number, &
                          go_to_line_number=line_number_macro, &
                          exit_if_eor=.true.,&
                          exit_if_eof=.true.) 
      
   end subroutine read_or_form_trans_and_diff_format_radmats

   ! --------------------------------------------------------------------------

   subroutine read_fission_data_format_radmats(format_radmats_file_unit, &
                                               line_number, &
                                               radmat, &
                                               number_non_thermal_groups_format_radmats, &
                                               record_len, &
                                               keyword_list)
            
      !!< Read in the fission, nu, production and prompt spectrum data for this radmat from a format_radmats file
      !!< If nu is read then the production is formed. Both the nu and production cannot both be present.  
      
      integer, intent(in) :: format_radmats_file_unit
      integer, intent(inout) :: line_number     
      type(radmat_type), intent(inout) :: radmat 
      integer, intent(in) :: number_non_thermal_groups_format_radmats
      integer, intent(in) :: record_len
      character(len=*), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      logical :: exit_if_eof      
      integer :: line_number_macro,first_keyword_found
      real :: sum_prompt_spectrum
      integer, dimension(:), allocatable :: keyword_find ! the specific keyword integer ids to find       
      character(len=record_len) :: line_string    

      ! save the line number of this macro so as to return to it
      line_number_macro = line_number
                  
      ! read the fission 
      allocate(keyword_find(2))    
      exit_if_eof = .false.      
      
      ! search for MACRO as well so as to not search to far each time
      keyword_find(1) = 1
      keyword_find(2) = 10

      call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                              line_number, &
                                              first_keyword_found, &
                                              exit_if_eof, &
                                              keyword_find, &
                                              keyword_list, &
                                              line_string)
            
      read_fission: if (first_keyword_found == 2) then 

         call read_format_radmats_siga_style_xsection(radmat%fission, &
                                                      line_number, &
                                                      format_radmats_file_unit, &
                                                      record_len)
         
         ! change the set flag
         radmat%fission_set = .true.
         
      end if read_fission
                                      
      ! put the file back to the start of this macro
      call go_to_line_seq(format_radmats_file_unit, &
                          line_number, &
                          go_to_line_number=line_number_macro, &
                          exit_if_eor=.true.,&
                          exit_if_eof=.true.) 

      ! read FISSNU or production. First search for FISSNU
      keyword_find(2) = 11

      call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                              line_number, &
                                              first_keyword_found, &
                                              exit_if_eof, &
                                              keyword_find, &
                                              keyword_list, &
                                              line_string)
            
      found_fissnu: if (first_keyword_found == 2) then 
         
         ! if the fission is not set then reading the FISSNU in doesnt make sense so exit
         fissnu_no_sigfiss: if (.not. radmat%fission_set) then 

            ewrite(-1,*) "Error reading format_radmats file"            
            FLExit("Found FISSNU without finding SIGFISS, this makes no sense")
         
         end if fissnu_no_sigfiss
         
         ! read FISSNU 
         call read_format_radmats_siga_style_xsection(radmat%particle_released_per_fission, &
                                                      line_number, &
                                                      format_radmats_file_unit, &
                                                      record_len)
         ! change the set flag
         radmat%particle_released_per_fission_set = .true.
         
         ! times fissnu by the fission to give the production variable         
         radmat%production = radmat%particle_released_per_fission * radmat%fission 

         ! change the set flag
         radmat%production_set = .true.
         
         ! check if the PRODUCTION keyword is also present for this macro - if so exit as cannot have both
      
         ! put the file back to the start of this macro
         call go_to_line_seq(format_radmats_file_unit, &
                             line_number, &
                             go_to_line_number=line_number_macro, &
                             exit_if_eor=.true.,&
                             exit_if_eof=.true.) 

         ! set to PRODUCTION keyword
         keyword_find(2) = 12
         
         ! search for PRODUCTION
         call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                                 line_number, &
                                                 first_keyword_found, &
                                                 exit_if_eof, &
                                                 keyword_find, &
                                                 keyword_list, &
                                                 line_string)
         
         found_production_also: if (first_keyword_found == 2) then
         
            ewrite(-1,*) "Error reading format_radmats file"
            FLExit("Both the FISSNU and PRODUCTION cross sections found")
         
         end if found_production_also
      
      else found_fissnu
      
         ! if not found the FISSNU then search for PRODUCTION

         ! put the file back to the start of this macro
         call go_to_line_seq(format_radmats_file_unit, &
                             line_number, &
                             go_to_line_number=line_number_macro, &
                             exit_if_eor=.true.,&
                             exit_if_eof=.true.) 

         ! set to PRODUCTION keyword
         keyword_find(2) = 12
         
         ! search for PRODUCTION
         call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                                 line_number, &
                                                 first_keyword_found, &
                                                 exit_if_eof, &
                                                 keyword_find, &
                                                 keyword_list, &
                                                 line_string)
      
         found_production: if (first_keyword_found == 2) then
            
            ! read the PRODUCTION variable
            call read_format_radmats_siga_style_xsection(radmat%production, &
                                                         line_number, &
                                                         format_radmats_file_unit, &
                                                         record_len)

            ! change the set flag
            radmat%production_set = .true.
                        
         else found_production
            
            ! if fission was found then we should expect FISSNU or PRODUCTION also
            sigfiss_no_fissnu_or_production: if (radmat%fission_set) then 
            
               ewrite(-1,*) "Error reading format_radmats file"
               FLExit("Neither the FISSNU and PRODUCTION cross sections found when SIGFISS was")
            
            end if sigfiss_no_fissnu_or_production
            
         end if found_production
      
      end if found_fissnu
               
      ! read the fission prompt spectrum and renormilse to account for computational io error
      ! also the fission prompt spectrum is only printed for non thermal groups in format_radmats, we have the whole array which is 
      ! already initialised to 0.0 so need to only read in for non thermal groups, the number of which has already been determined
      keyword_find(2) = 13
      
      ! put the file back to the start of this macro
      call go_to_line_seq(format_radmats_file_unit, &
                          line_number, &
                          go_to_line_number=line_number_macro, &
                          exit_if_eor=.true.,&
                          exit_if_eof=.true.) 

      call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                              line_number, &
                                              first_keyword_found, &
                                              exit_if_eof, &
                                              keyword_find, &
                                              keyword_list, &
                                              line_string)
            
      read_fisschi: if (first_keyword_found == 2) then

         ! if the production is not set then reading the prompt spectrum in doesnt make sense so exit
         fisschi_no_production: if (.not. radmat%production_set) then 

            ewrite(-1,*) "Error reading format_radmats file"            
            FLExit("Found FISSCHI without setting PRODUCTION, this makes no sense")
         
         end if fisschi_no_production
                        
         call read_format_radmats_siga_style_xsection(radmat%prompt_spectrum, &
                                                      line_number, &
                                                      format_radmats_file_unit, &
                                                      record_len, &
                                                      number_groups_to_read_in=number_non_thermal_groups_format_radmats)
                  
         !renormalise spectrum
         sum_prompt_spectrum = sum(radmat%prompt_spectrum)
                  
         ! check sum of spectrum is close to 1.0 - this also avoids a divide by zero by accident
         check_sum_spectrum: if (abs(sum_prompt_spectrum - 1.0) < 1.0e-03) then
                     
            radmat%prompt_spectrum = radmat%prompt_spectrum / sum_prompt_spectrum
                  
         else check_sum_spectrum
                  
            ewrite(-1,*) "Read in error for format_radmats prompt spectrum"
            FLExit("Poor spectrum for material")
                  
         end if check_sum_spectrum

         ! change the set flag
         radmat%prompt_spectrum_set = .true.
      
      else read_fisschi
 
         ! if PRODUCTION was found then we should expect FISSCHI also
         production_not_fisschi: if (radmat%production_set) then 
            
            ewrite(-1,*) "Error reading format_radmats file"
            FLExit("FISSCHI not found when PRODUCTION was set")
            
         end if production_not_fisschi
               
      end if read_fisschi
      
      deallocate(keyword_find)

      ! put the file back to the start of this macro
      call go_to_line_seq(format_radmats_file_unit, &
                          line_number, &
                          go_to_line_number=line_number_macro, &
                          exit_if_eor=.true.,&
                          exit_if_eof=.true.) 
      
   end subroutine read_fission_data_format_radmats

   ! --------------------------------------------------------------------------
   
   subroutine read_velocity_style_data_format_radmats(format_radmats_file_unit, &
                                                      dataset_radmat, &
                                                      number_of_groups, &
                                                      keyword_find, &
                                                      record_len, &
                                                      keyword_list)
      
      !!< Read in the VELOCITY/POWER/ENERGY_RELEASED/BETA style data which may not be material dependent from a format_radmats format file
      !!< What is read in is determined by the input keyword character id. The number_groups could refer to energy or delayed.
      !!< The expected structure of the data after a keyword is rigid in that there should be no blank lines, no leading spaces, 
      !!< all data for each group on one line and ascending sequential macro data. IF ALL keyword found then all materials have 
      !!< same data and the input file is only expected to have this (any further numbers are not considered).
      
      integer, intent(in) :: format_radmats_file_unit
      type(dataset_radmat_type), intent(inout) :: dataset_radmat
      integer, intent(in) :: number_of_groups
      integer, dimension(1), intent(in) :: keyword_find
      integer, intent(in) :: record_len
      character(len=*), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      integer :: g,pmat,rmat,line_number,first_keyword_found
      real, dimension(:), allocatable :: values
      logical :: all,exit_if_eof
      character(len=record_len) :: line_string
      logical :: end_of_file

      ! go back to start of file
      call rewind_file_seq(format_radmats_file_unit, &
                           line_number)
                
      ! first find the line with the keyword  - which could be anywhere in the file   
      exit_if_eof = .true.      
      call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                              line_number, &
                                              first_keyword_found, &
                                              exit_if_eof, &
                                              keyword_find, &
                                              keyword_list, &
                                              line_string)
            
      ! now inspect the next line for either the keyword ALL at the start
      call read_next_line_seq(format_radmats_file_unit, &
                              line_string, &
                              line_number, &
                              end_of_file, &
                              exit_if_eor=.true., &
                              exit_if_eof=.true.)
            
      inspect_string: if (line_string(1:3) == 'ALL') then 
         
         all = .true.
      
      else 
         
         all = .false.
      
      end if inspect_string

      ! initialise the values array
      allocate(values(number_of_groups))
         
      values = 0.0
      
      read_all_values: if (all) then
                        
         ! read the next line
         call read_next_line_seq(format_radmats_file_unit, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor=.true., &
                                 exit_if_eof=.true.)
         
         ! read the values of the line_string
         read(line_string,*) ( values(g), g = 1,number_of_groups )
               
      else read_all_values
         
         ! go back one line if not allmacro
         call go_to_line_seq(format_radmats_file_unit, &
                             line_number, &
                             go_to_line_number=line_number-1, &
                             exit_if_eor=.true., &
                             exit_if_eof=.true.)
      
      end if read_all_values
         
      ! set the values for each radmat into the relevant arrays depending on the keyword
      physical_material_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
      
         radmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)

            read_values: if (.not. all) then
                                       
               ! read the next line
               call read_next_line_seq(format_radmats_file_unit, &
                                       line_string, &
                                       line_number, &
                                       end_of_file, &
                                       exit_if_eor=.true., &
                                       exit_if_eof=.true.)
               
               ! read the values from the line_string
               read(line_string,*) ( values(g), g = 1,number_of_groups )
                    
            end if read_values
               
            set_data: if (keyword_find(1) == 14) then
      
               dataset_radmat%physical_radmats(pmat)%radmats(rmat)%velocity = values

               dataset_radmat%physical_radmats(pmat)%radmats(rmat)%velocity_set = .true.
      
            else if (keyword_find(1) == 15) then set_data
      
               dataset_radmat%physical_radmats(pmat)%radmats(rmat)%power = values
               
               dataset_radmat%physical_radmats(pmat)%radmats(rmat)%power_set = .true.

            else if (keyword_find(1) == 16) then set_data
      
               dataset_radmat%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission = values
               
               dataset_radmat%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission_set = .true.
      
            else if (keyword_find(1) == 17) then set_data
      
               dataset_radmat%physical_radmats(pmat)%radmats(rmat)%beta = values

               dataset_radmat%physical_radmats(pmat)%radmats(rmat)%beta_set = .true.
      
            else set_data
                  
               FLAbort("Error trying to read velocity style data from a format_radmats format file with unrecognised first keyword") 
                  
            end if set_data 
               
         end do radmat_loop
         
      end do physical_material_loop
               
      deallocate(values)

      ! go back to start of file
      call rewind_file_seq(format_radmats_file_unit, &
                           line_number)
      
   end subroutine read_velocity_style_data_format_radmats
   
   ! --------------------------------------------------------------------------
   
   subroutine read_delayed_data_format_radmats(format_radmats_file_unit, &
                                               dataset_radmat, &
                                               record_len, &
                                               keyword_list, &
                                               delayed_lambda_spectrum, &
                                               read_delayed_lambda_spectrum)
   
      !!< Read in the delayed data for this radiation material data set from a format_radmats formatted file
      !!< The read in of the delayed lambda and spectrum into this data depends on an inout logical
      
      integer, intent(in) :: format_radmats_file_unit
      type(dataset_radmat_type), intent(inout) :: dataset_radmat
      integer, intent(in) :: record_len 
      character(len=*), dimension(:), intent(in) :: keyword_list
      type(delayed_lambda_spectrum_type), intent(inout) :: delayed_lambda_spectrum 
      logical, intent(in) :: read_delayed_lambda_spectrum
      
      ! local variables
      logical :: exit_if_eof
      integer :: d,g,line_number,first_keyword_found
      integer :: number_of_energy_groups,number_of_delayed_groups
      integer, dimension(:), allocatable :: keyword_find
      character(len=record_len) :: line_string  
      logical :: end_of_file   
       
      allocate(keyword_find(1)) 
      
      ! find the number of groups and delayed groups from the size of the delayed spectrum
      number_of_energy_groups = size(delayed_lambda_spectrum%spectrum,1)
      number_of_delayed_groups = size(delayed_lambda_spectrum%spectrum,2)
           
      ! read the delayed BETA
      keyword_find(1) = 17     
      call read_velocity_style_data_format_radmats(format_radmats_file_unit, &
                                                   dataset_radmat, &
                                                   number_of_delayed_groups, &
                                                   keyword_find, &
                                                   record_len, &
                                                   keyword_list)
      
      ! if needed from this data set read in the delayed lambda and spectrum
      read_lambda_spectrum: if (read_delayed_lambda_spectrum) then
      
         keyword_find(1) = 18
         exit_if_eof = .true.

         ! go back to start of file
         call rewind_file_seq(format_radmats_file_unit, &
                              line_number)

         call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                                 line_number, &
                                                 first_keyword_found, &
                                                 exit_if_eof, &
                                                 keyword_find, &
                                                 keyword_list, &
                                                 line_string)
         
         ! read the next line 
         call read_next_line_seq(format_radmats_file_unit, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor=.true., &
                                 exit_if_eof=.true.)
         
         ! read the values from line_string
         read(line_string,*) ( delayed_lambda_spectrum%lambda(d), d = 1,number_of_delayed_groups )
                  
         delayed_lambda_spectrum%lambda_set = .true.
         
         keyword_find(1) = 19
         exit_if_eof = .true.
 
         ! go back to start of file
         call rewind_file_seq(format_radmats_file_unit, &
                              line_number)

         call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                                 line_number, &
                                                 first_keyword_found, &
                                                 exit_if_eof, &
                                                 keyword_find, &
                                                 keyword_list, &
                                                 line_string)
         
         ! read spectrum for each delayed group
         delayed_group_loop: do d = 1,number_of_delayed_groups
         
            ! read the next line 
            call read_next_line_seq(format_radmats_file_unit, &
                                    line_string, &
                                    line_number, &
                                    end_of_file, &
                                    exit_if_eor=.true., &
                                    exit_if_eof=.true.)
            
            ! read the values from the line_string
            read(line_string,*) ( delayed_lambda_spectrum%spectrum(g,d), g = 1,number_of_energy_groups )
                  
         end do delayed_group_loop
         
         delayed_lambda_spectrum%spectrum_set = .true.
                  
      end if read_lambda_spectrum
            
      deallocate(keyword_find)        

      ! go back to start of file
      call rewind_file_seq(format_radmats_file_unit, &
                           line_number)
      
   end subroutine read_delayed_data_format_radmats

   ! --------------------------------------------------------------------------

   subroutine read_power_data_format_radmats(format_radmats_file_unit, &
                                             dataset_radmat, &
                                             record_len, &
                                             keyword_list)
   
      !!< Read in the power or energy released per fission data (then form power for each material if fission present)
      !!< cannot have both and will exit if both found
      
      integer, intent(in) :: format_radmats_file_unit
      type(dataset_radmat_type), intent(inout) :: dataset_radmat
      integer, intent(in) :: record_len 
      character(len=*), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      logical :: exit_if_eof
      integer :: line_number,first_keyword_found
      integer :: number_of_energy_groups
      integer, dimension(:), allocatable :: keyword_find
      character(len=record_len) :: line_string  
       
      ! go back to start of file
      call rewind_file_seq(format_radmats_file_unit, &
                           line_number)
      
      ! determine the number of groups
      number_of_energy_groups = size(dataset_radmat%physical_radmats(1)%radmats(1)%power)
      
      allocate(keyword_find(1)) 
      keyword_find(1) = 15
      exit_if_eof = .false.
      
      ! search for POWER first
      call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                              line_number, &
                                              first_keyword_found, &
                                              exit_if_eof, &
                                              keyword_find, &
                                              keyword_list, &
                                              line_string)
      
      found_power: if (first_keyword_found == 1) then
         
         ! read POWER 
         call read_velocity_style_data_format_radmats(format_radmats_file_unit, &
                                                      dataset_radmat, &
                                                      number_of_energy_groups, &
                                                      keyword_find, &
                                                      record_len, &
                                                      keyword_list)                           
         
         ! check if the ENERGY_RELEASED keyword is also present - if so exit as cannot have both
      
         ! go back to start of file
         call rewind_file_seq(format_radmats_file_unit, &
                              line_number)

         ! set to ENERGY_RELEASED keyword
         keyword_find(1) = 16
         
         ! search for ENERGY_RELEASED
         call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                                 line_number, &
                                                 first_keyword_found, &
                                                 exit_if_eof, &
                                                 keyword_find, &
                                                 keyword_list, &
                                                 line_string)
         
         found_energy_released_also: if (first_keyword_found == 1) then
         
            ewrite(-1,*) "Error reading format_radmats file"
            FLExit("Both the POWER and ENERGY_RELEASED data found when the flml options suggest it is needed")
         
         end if found_energy_released_also
      
      else found_power
      
         ! if not found the POWER then search for ENERGY_RELEASED

         ! go back to start of file
         call rewind_file_seq(format_radmats_file_unit, &
                              line_number)

         ! set to ENERGY_RELEASED keyword
         keyword_find(1) = 16
         
         ! search for ENERGY_RELEASED
         call find_line_with_any_desired_keyword(format_radmats_file_unit, &
                                                 line_number, &
                                                 first_keyword_found, &
                                                 exit_if_eof, &
                                                 keyword_find, &
                                                 keyword_list, &
                                                 line_string)
      
         found_energy_released: if (first_keyword_found == 1) then
            
            ! read the ENERGY_RELEASED variable
            call read_velocity_style_data_format_radmats(format_radmats_file_unit, &
                                                         dataset_radmat, &
                                                         number_of_energy_groups, &
                                                         keyword_find, &
                                                         record_len, &
                                                         keyword_list)         
            
           ! calculate the power variable from the fission and energy released per fission for this dataset_radmat
           call calculate_power_format_radmats(dataset_radmat)
                        
         else found_energy_released
         
            ewrite(-1,*) "Error reading format_radmats file"
            FLExit("Neither the POWER and ENERGY_RELEASED data found when the flml options suggest it is needed")
         
         end if found_energy_released
      
      end if found_power
      
      deallocate(keyword_find)
      
      ! go back to start of file
      call rewind_file_seq(format_radmats_file_unit, &
                           line_number)
      
   end subroutine read_power_data_format_radmats

   ! --------------------------------------------------------------------------

   subroutine calculate_power_format_radmats(dataset_radmat)
   
      !!< Calculate the power variable from the fission and energy released per fission for this dataset_radmat
      !!< This is called directly after the energy_released_per_fission has been read in.
      !!< This only happens for materials that have already had the fission set.
            
      type(dataset_radmat_type), intent(inout) :: dataset_radmat
      
      ! local variables
      integer pmat,rmat
      
      physical_material_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
      
         radmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)
            
            fission_set: if (dataset_radmat%physical_radmats(pmat)%radmats(rmat)%fission_set) then
                                            
               dataset_radmat%physical_radmats(pmat)%radmats(rmat)%power = &
               dataset_radmat%physical_radmats(pmat)%radmats(rmat)%fission * &
               dataset_radmat%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission
                              
               dataset_radmat%physical_radmats(pmat)%radmats(rmat)%power_set = .true.
                           
            end if fission_set
                                                         
         end do radmat_loop
         
      end do physical_material_loop
      
   end subroutine calculate_power_format_radmats

   ! --------------------------------------------------------------------------

   subroutine read_format_radmats_siga_style_xsection(siga_style_xsection, &
                                                      line_number, &
                                                      format_radmats_file_unit, &
                                                      record_len, &
                                                      number_groups_to_read_in)
      
      !!< Actually read in a absorption style cross section that is off 
      !!< length number_of_energy_groups from a format_radmats format file 
      
      real, dimension(:), intent(inout) :: siga_style_xsection
      integer, intent(inout) :: line_number
      integer, intent(in) :: format_radmats_file_unit
      integer, intent(in) :: record_len
      integer, optional, intent(in) :: number_groups_to_read_in
      
      ! local variables
      integer :: w
      integer :: group_count
      integer :: g_read_in
      integer :: number_of_energy_groups
      logical :: end_of_file
      character(len=record_len) :: line_string   
      character(len=record_len), dimension(:), allocatable :: words  
       
      ! find the number of groups from the size of sigabs
      number_of_energy_groups = size(siga_style_xsection)
      
      ! the line number is assumed at the line with the KEYWORD SIGABS (or other etc.)
      
      ! decude the number of expected groups to read in
      diff_groups: if (present(number_groups_to_read_in)) then
      
         g_read_in = number_groups_to_read_in

      else diff_groups
      
         g_read_in = number_of_energy_groups
      
      end if diff_groups
      
      ! initialise the group count
      group_count = 0
      
      ! keep reading lines until have all necessary data or error            
      read_sigabs_loop: do
                                    
         ! read the next line
         call read_next_line_seq(format_radmats_file_unit, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor=.true., &
                                 exit_if_eof=.true.)
         
         ! read the words from the line_string - words are substrings seperated 
         ! by any number of blank spaces
         call read_words_from_string(line_string, &
                                     words)
         
         ! read the real data from the words array         
         word_loop: do w = 1,size(words) 
            
            group_count = group_count + 1
                        
            ! read the value from the word substring, if this isnt a real
            ! it will fail here
            read(words(w),*) siga_style_xsection(group_count)

            ! exit if read all necessary group data
            if (group_count == g_read_in) exit read_sigabs_loop
            
         end do word_loop
                                 
      end do read_sigabs_loop  
            
      if (allocated(words)) deallocate(words)
      
   end subroutine read_format_radmats_siga_style_xsection
   
   ! --------------------------------------------------------------------------

   subroutine read_format_radmats_sigs_style_xsection(scatter, &
                                                      line_number, &
                                                      format_radmats_file_unit, &
                                                      record_len, &
                                                      keyword_list)

      !!< Actually read in a scatter moment style cross section that is off length 
      !!< number_of_energy_groups*number_of_energy_groups from a format_radmats format file 
      !!< This assumes there is only ever one first-last pair for each group
      
      real, dimension(:,:), intent(inout) :: scatter
      integer, intent(inout) :: line_number
      integer, intent(in) :: format_radmats_file_unit
      integer, intent(in) :: record_len
      character(len=*), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      logical :: exit_if_eof
      logical :: found_correct_g
      integer :: g,gstart,gend,g_to_read,g_check,g_found
      integer :: line_number_scatter
      integer :: number_of_energy_groups,first_keyword_found
      character(len=record_len) :: cdummy,c_gstart,c_gend,c_g_found
      character(len=record_len) :: line_string
               
      ! find the number of groups from the size of the scatter
      number_of_energy_groups = size(scatter,1)
      
      ! save the liner number with the keyword SCATTER 
      line_number_scatter = line_number
                  
      group_loop: do g = 1,number_of_energy_groups
         
         exit_if_eof = .true.
         
         ! have a flag to say whether or not the correct group number is found
         found_correct_g = .false.
         
         group_check_loop: do g_check = 1,number_of_energy_groups
         
            ! find the next line with GROUP on it - which should be the next keyword         
            call find_line_with_any_keyword(format_radmats_file_unit, &
                                            line_number, &
                                            first_keyword_found, &
                                            exit_if_eof, &
                                            keyword_list, &
                                            line_string) 
                  
            ! check that we found a GROUP, rather than anything else
            check_found: if (first_keyword_found /= 20) then
            
               ewrite(-1,*) "Did not find GROUP keyword for group",g_check,"when reading scatter block"
               FLExit("Error reading format_radmats style scatter block")
         
            end if check_found 
            
            ! now read the GROUP line to find which GROUP this is
            read(line_string,*) cdummy,c_g_found
            
            read(c_g_found,*) g_found
            
            ! if this is the correct GROUP then exit loop
            found: if (g == g_found) then
               
               found_correct_g = .true.
               
               exit group_check_loop 
            
            end if found
         
         end do group_check_loop
         
         ! exit if not found the correct g number
         not_found: if (.not. found_correct_g) then
            
            ewrite(1,*) 'Error searching for GROUP block ',g
            FLExit('Error reading radiation radmats file as GROUP block not found')
         
         end if not_found
            
         ! now read the GROUP line (in the string) to find gstart and gend
         
         ! read the sub strings of the string
         read(line_string,*) cdummy,cdummy,cdummy,c_gstart,cdummy,c_gend
         
         ! read the integers of the substrings
         read(c_gstart,*) gstart
         read(c_gend,*) gend
       
         ! calc the number of groups to read in
         g_to_read = gend - gstart + 1
         
         ! read this group g to other groups as needed (gstart:gend) scatter moments 
         call read_format_radmats_siga_style_xsection(scatter(g,gstart:gend), &
                                                      line_number, &
                                                      format_radmats_file_unit, &
                                                      record_len, &
                                                      number_groups_to_read_in=g_to_read)

         ! put the file back to the start of this SCATTER block
         call go_to_line_seq(format_radmats_file_unit, &
                             line_number, &
                             go_to_line_number = line_number_scatter, &
                             exit_if_eor       = .true.,&
                             exit_if_eof       = .true.) 
      
      end do group_loop
                  
   end subroutine read_format_radmats_sigs_style_xsection

   ! --------------------------------------------------------------------------
   
   subroutine count_number_lines_with_keyword(file_unit, &
                                             line_number, &
                                             keyword_count, &
                                             keyword_to_count, &
                                             record_len, &
                                             keyword_list, &
                                             keyword_to_stop)
      
      !!< Count the number of lines found with specific keyword before either finding keyword_to_stop (if present) or EoF
      !!< If keyword_to_stop is set to 0 then the counting stops if ANY other keyword is found

      integer, intent(in) :: file_unit  
      integer, intent(inout) :: line_number
      integer, intent(inout) :: keyword_count
      integer, intent(in) :: keyword_to_count          
      integer, intent(in) :: record_len
      character(len=*), dimension(:), intent(in) :: keyword_list   
      integer, intent(in), optional :: keyword_to_stop
            
      ! local variables 
      logical :: exit_if_eof
      integer :: first_keyword_found
      character(len=record_len) :: line_string
      
      ! some initial checks
      initial_check: if (present(keyword_to_stop)) then
      
         ! exit if keyword_to_stop is < 0
         if (keyword_to_stop < 0) FLAbort("In subroutine count_number_keyword and keyword_to_stop cannot be negative")
      
         ! exit if keyword_to_stop is to big
         if (keyword_to_stop > size(keyword_list)) FLAbort("In subroutine count_number_keyword and keyword_to_stop set to number to big")
      
      end if initial_check
      
      ! initialise the keyword_count
      keyword_count = 0
      
      ! no need to exit if eof
      exit_if_eof = .false.
      
      count_loop: do 
      
         call find_line_with_any_keyword(file_unit, &
                                         line_number, &
                                         first_keyword_found, &
                                         exit_if_eof, &
                                         keyword_list, &
                                         line_string) 
         
         found: if (first_keyword_found == keyword_to_count) then
         
            keyword_count = keyword_count + 1
         
         else if (first_keyword_found == 0) then found
            
            ! end of file
            exit count_loop
         
         else found
            
            stop_present: if (present(keyword_to_stop)) then
            
               one_stop: if (keyword_to_stop > 0) then
            
                  check: if (first_keyword_found == keyword_to_stop) then
               
                     exit count_loop
               
                  else check
               
                     cycle count_loop
               
                  end if check
            
               else one_stop
            
                  exit count_loop
            
               end if one_stop

            end if stop_present
            
         end if found
         
      end do count_loop  
          
   end subroutine count_number_lines_with_keyword
   
   ! --------------------------------------------------------------------------

   subroutine find_line_with_any_keyword(file_unit, &
                                         line_number, &
                                         first_keyword_found, &
                                         exit_if_eof, &
                                         keyword_list, &
                                         line_string) 
      
      !!< Find the first line that contains any of the keywords to set the variable 
      
      integer, intent(in) :: file_unit
      integer, intent(inout) :: line_number
      integer, intent(out) :: first_keyword_found
      logical, intent(in) :: exit_if_eof
      character(len=*), dimension(:), intent(in) :: keyword_list
      character(len=*), intent(inout) :: line_string
      
      ! local variables
      integer :: i                  
      integer, dimension(:), allocatable :: keyword_find ! the specific keyword integer ids to find
      
      allocate(keyword_find(size(keyword_list)))
      
      set_keyword_find: do i = 1,size(keyword_list)
      
         keyword_find(i) = i
         
      end do set_keyword_find
            
      call find_line_with_any_desired_keyword(file_unit, &
                                              line_number, &
                                              first_keyword_found, &
                                              exit_if_eof, &
                                              keyword_find, &
                                              keyword_list, &
                                              line_string)
      
      deallocate(keyword_find)
      
   end subroutine find_line_with_any_keyword
    
   ! --------------------------------------------------------------------------
      
   subroutine find_line_with_any_desired_keyword(file_unit, &
                                                 line_number, &
                                                 first_keyword_found, &
                                                 exit_if_eof, &
                                                 keyword_find, &
                                                 keyword_list, &
                                                 line_string)  
      
      !!< Sequentially read a file with the input unit given, from the current line number line tagged line_number
      !!< until finiding a line which contains a desired keyword from the list of keywords.
      !!< The last line number read, the last line read as a string, which (if any) keyword was found and an
      !!< end of file logical are returned
      
      ! if exit_if_eof is true then flexit occurs if the end of file is reached
      
      integer, intent(in) :: file_unit
      integer, intent(inout) :: line_number
      integer, intent(out) :: first_keyword_found
      logical, intent(in) :: exit_if_eof
      integer, dimension(:), intent(in) :: keyword_find
      character(len=*), dimension(:), intent(in) :: keyword_list
      character(len=*), intent(inout) :: line_string   
         
      ! local variables
      integer :: k
      logical :: end_of_file
      character(len=len(line_string)) :: first_word
                  
      first_keyword_found = 0
                
      read_line: do
       
         call read_next_line_seq(file_unit, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor=.true., &
                                 exit_if_eof=exit_if_eof)
      
         if (end_of_file) exit read_line
         
         ! cycle blank lines
         if (trim(line_string) == '') cycle read_line
         
         ! initialise the first_word
         first_word = ''
          
         ! read the first word of the line_string
         read(line_string,*) first_word
         
         ! make sure the first_word is all upper case as all the keywords are upper case
         call all_upper_case(first_word)
                                         
         ! check if the first word is any of the keywords - exit when found a match
         each_keyword_loop: do k = 1,size(keyword_find)

            found_keyword: if (trim(first_word) == trim(keyword_list(keyword_find(k)))) then 
                 
               first_keyword_found = k
                  
               exit read_line
               
            end if found_keyword
               
         end do each_keyword_loop 
         
      end do  read_line  
            
   end subroutine find_line_with_any_desired_keyword

   ! --------------------------------------------------------------------------
   
   subroutine all_upper_case(string)
      
      !!< Make the input string all upper case
      
      character(len=*), intent(inout) :: string
      
      ! local variables
      integer :: i
      character(len=1) :: c
            
      char_loop: do i = 1,len(string)   
      
         c = string(i:i)
         
         call make_character_upper_case(c)
         
         string(i:i) = c
      
      end do char_loop

   end subroutine all_upper_case
   
   ! --------------------------------------------------------------------------
    
   subroutine make_character_upper_case(c)
      
      !!< Make the character if needed upper case using the ASCII sequence 
      !!< function assuming decimal value 
      
      character(len=*), intent(inout) :: c
      
      ! local variables
      integer :: i
      
      ! first assert that the length of c is 1
      assert(len(c) == 1)
      
      ! find the integer corresponding to this character in the ASCII sequence
      i = iachar(c)
      
      ! the ASCII decimal alphabet lower case characters are from 97 to 122
      ! with the corresponding upper case being -32
      lower_found: if ((97 <= i) .and. (i <= 122)) then
      
         c = achar(i - 32)
      
      end if lower_found      
   
   end subroutine make_character_upper_case
   
   ! --------------------------------------------------------------------------
   
   subroutine read_words_from_string(string, &
                                     words)
   
      !!< Read all the words from a input string. Words are any sub string  
      !!< seperated by any number of spaces
      
      character(len=*), intent(in) :: string
      character(len=*), dimension(:), allocatable, intent(inout) :: words
      
      ! local variables
      integer :: w
      integer :: number_of_words
      
      ! deallocate words if it is already allocated
      if (allocated(words)) deallocate(words)
      
      ! find the number of words within the string
      call string_word_count(string, &
                             number_of_words)
      
      ! initialise the words
      allocate(words(number_of_words))
      
      words = ''
      
      ! read the words from the string      
      read(string,*) (words(w), w = 1,number_of_words)
                  
   end subroutine read_words_from_string 
   
   ! --------------------------------------------------------------------------
   
   subroutine string_word_count(string, &
                                number_of_words)
      
      !!< Count the number of words within a string. Words are any sub string
      !!< seperated by any number of spaces

      character(len=*), intent(in) :: string
      integer, intent(out) :: number_of_words
      
      ! local variables
      integer :: c
      logical :: found_space
      
      ! initialise the word count
      number_of_words = 0
      
      ! if the string is empty then return
      if (len_trim(string) == 0) return
      
      ! sweep each character of the string inspecting for spaces ' '
      found_space = .true.
      char_loop: do c = 1,len_trim(string)
         
         ! check if the character is a space
         check_for_space: if (string(c:c) == ' ') then
            
            found_space = .true.
                     
         else if (found_space) then
               
             number_of_words = number_of_words + 1
             
             found_space = .false.
                           
         end if check_for_space

      end do char_loop
      
   end subroutine string_word_count
   
   ! --------------------------------------------------------------------------
   
   function number_substrings_within_string(main_string,sub_string) 
      
      !!< count the number of whole non overlapping sub_string's within the main_string
      
      character(len=*), intent(in) :: main_string
      character(len=*), intent(in) :: sub_string
      
      integer :: number_substrings_within_string
      
      ! local variables
      integer :: substring_lower_bound
      integer :: substring_pos
      
      number_substrings_within_string = 0
      
      substring_lower_bound = 1
      
      substring_pos = 0

      count_sub_loop: do 

         substring_pos = index(trim(main_string(substring_lower_bound:)),trim(sub_string)) 
      
         found_substring_pos: if (substring_pos > 0) then
              
            number_substrings_within_string = number_substrings_within_string + 1
              
            substring_lower_bound = substring_lower_bound + substring_pos + len(sub_string) - 1
       
         else found_substring_pos
           
            exit count_sub_loop
           
         end if found_substring_pos

      end do count_sub_loop
   
   end function number_substrings_within_string

   ! --------------------------------------------------------------------------
   
   subroutine read_next_line_seq(file_unit, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor, &
                                 exit_if_eof)
   
      !!< Read the next line in a sequential file and increase the line number tag by 1
      !!< FLExit if end of record if allowed, FLExit if end of file if allowed and FLAbort if other error
      
      integer, intent(in) :: file_unit
      character(len=*), intent(inout) :: line_string
      integer, intent(inout) :: line_number
      logical, intent(inout) :: end_of_file
      logical, intent(in) :: exit_if_eor
      logical, intent(in) :: exit_if_eof
      
      end_of_file = .false.
      
      line_number = line_number + 1    
      
      read(file_unit,'(a)',iostat=iostat) line_string
      
      iostat_if: if (iostat == -2) then
       
         if (exit_if_eor) FLExit("Problem reading file using read_next_line_seq, iostat = -2") 
      
      else if (iostat == -1) then iostat_if
         
         end_of_file = .true.
         
         if (exit_if_eof) FLExit("Problem reading file using read_next_line_seq, iostat = -1") 
       
      else if (iostat > 0) then iostat_if 
      
         FLAbort("Problem reading file using read_next_line_seq, iostat > 0") 
      
      end if iostat_if
                        
   end subroutine read_next_line_seq

   ! --------------------------------------------------------------------------
   
   subroutine read_previous_line_seq(file_unit, &
                                     line_string, &
                                     line_number, &
                                     exit_if_eor, &
                                     exit_if_eof)
   
      !!< Read the previous line of sequential file and decrease the line number tag by 1
      !!< FLExit if end of record if allowed, FLExit if end of file if allowed and FLAbort if other error
      
      integer, intent(in) :: file_unit
      character(len=*), intent(inout) :: line_string
      integer, intent(inout) :: line_number
      logical, intent(in) :: exit_if_eor
      logical, intent(in) :: exit_if_eof
      
      line_number = line_number - 1    
      
      ! two back spaces are required         
      backspace(unit=file_unit,iostat=iostat)

      iostat_if1: if (iostat == -2) then
       
         if (exit_if_eor) FLExit("Problem backspaceing file using read_previous_line_seq, iostat = -2") 
      
      else if (iostat == -1) then iostat_if1

         if (exit_if_eof) FLExit("Problem backspaceing file using read_previous_line_seq, iostat = -1") 
       
      else if (iostat > 0) then iostat_if1 
      
         FLAbort("Problem backspaceing file using read_previous_line_seq, iostat > 0") 
      
      end if iostat_if1
           
      backspace(unit=file_unit,iostat=iostat)

      iostat_if2: if (iostat == -2) then
       
         if (exit_if_eor) FLExit("Problem backspaceing file using read_previous_line_seq, iostat = -2") 
      
      else if (iostat == -1) then iostat_if2

         if (exit_if_eof) FLExit("Problem backspaceing file using read_previous_line_seq, iostat = -1") 
       
      else if (iostat > 0) then iostat_if2 
      
         FLAbort("Problem backspaceing file using read_previous_line_seq, iostat > 0") 
      
      end if iostat_if2

      read(file_unit,'(a)',iostat=iostat) line_string

      iostat_if3: if (iostat == -2) then
       
         if (exit_if_eor) FLExit("Problem reading file using read_previous_line_seq, iostat = -2") 
      
      else if (iostat == -1) then iostat_if3

         if (exit_if_eof) FLExit("Problem reading file using read_previous_line_seq, iostat = -1") 
       
      else if (iostat > 0) then iostat_if3 
      
         FLAbort("Problem reading file using read_previous_line_seq, iostat > 0") 
      
      end if iostat_if3
         
   end subroutine read_previous_line_seq

   ! --------------------------------------------------------------------------
   
   subroutine rewind_file_seq(file_unit, &
                              line_number)
   
      !!< Return to the first line of file using rewind and set line number tag to 1
      
      integer, intent(in) :: file_unit
      integer, intent(out) :: line_number

      line_number = 0
      
      rewind(unit=file_unit,iostat=iostat)
      
      if (.not. iostat == 0) FLAbort("Problem with rewind of file using rewind_file_seq")
      
   end subroutine rewind_file_seq

   ! --------------------------------------------------------------------------
   
   subroutine go_to_line_seq(file_unit, &
                             line_number, &
                             go_to_line_number, &
                             exit_if_eor, &
                             exit_if_eof)
   
      !!< From the current line go to the desired line of a sequential file
      !!< FLExit if end of record if allowed, FLExit if end of file if allowed and FLAbort if other error
      
      integer, intent(in) :: file_unit
      integer, intent(inout) :: line_number
      integer, intent(in) :: go_to_line_number
      logical, intent(in) :: exit_if_eor
      logical, intent(in) :: exit_if_eof
      
      ! local variables
      integer :: l,difference_lines

      direction: if (go_to_line_number == line_number) then
      
         ! nothing to do
      
      else if (go_to_line_number > line_number) then
      
         difference_lines = go_to_line_number - line_number
         
         line_loop1: do l = 1,difference_lines
         
            line_number = line_number + 1
            
            read(file_unit,*,iostat=iostat)
            
            iostat_if1: if (iostat == -2) then
       
               if (exit_if_eor) FLExit("Problem reading file using go_to_line_seq, iostat = -2") 
      
            else if (iostat == -1) then iostat_if1

               if (exit_if_eof) FLExit("Problem reading file using go_to_line_seq, iostat = -1") 
       
            else if (iostat > 0) then iostat_if1 
      
               FLAbort("Problem reading file using go_to_line_seq, iostat > 0") 
      
            end if iostat_if1  
                           
         end do line_loop1
         
      else direction
      
         difference_lines = line_number - go_to_line_number
         
         line_loop2: do l = 1,difference_lines
         
            line_number = line_number - 1
            
            backspace(unit=file_unit,iostat=iostat)
            
            iostat_if2: if (iostat == -2) then
       
               if (exit_if_eor) FLExit("Problem backspaceing file using go_to_line_seq, iostat = -2") 
      
            else if (iostat == -1) then iostat_if2

               if (exit_if_eof) FLExit("Problem backspaceing file using go_to_line_seq, iostat = -1") 
       
            else if (iostat > 0) then iostat_if2 
      
               FLAbort("Problem backspaceing file using go_to_line_seq, iostat > 0") 
      
            end if iostat_if2          
              
         end do line_loop2
      
      end if direction   
      
   end subroutine go_to_line_seq

   ! --------------------------------------------------------------------------
   
   subroutine desired_keyword_not_first_found(keyword_find, &
                                              desired_keyword, &
                                              first_keyword_found, &
                                              keyword_list)
   
      !!< exit as the desired keyword was not the first found
      
      integer, dimension(:), intent(in) :: keyword_find
      integer, intent(in) :: desired_keyword
      integer, intent(in) :: first_keyword_found
      character(len=*), dimension(:), intent(in) :: keyword_list
      
      ewrite(-1,*) "Desired keyword ",trim(keyword_list(keyword_find(desired_keyword)))," not found first"
      ewrite(-1,*) "Found keyword ",trim(keyword_list(keyword_find(first_keyword_found)))," first instead"
      FLExit("Keyword find error in read file")
      
   end subroutine desired_keyword_not_first_found

   ! --------------------------------------------------------------------------

end module radiation_materials_read_format_radmats_base
