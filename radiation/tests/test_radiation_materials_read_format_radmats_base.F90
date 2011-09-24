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

subroutine test_radiation_materials_read_format_radmats_base
   
   !!< Test the procedures contained within the module 
   !!< radiation_materials_read_format_radmats_base 
   !!< in Radiation_Materials_Read_Format_Radmats_Base.F90
   
   use futils
   use unittest_tools  

   use radiation_materials_read_format_radmats_base 
   use radiation_materials_data_types
   use radiation_materials_create
   use radiation_materials_destroy
   use create_unittest_input 
         
   implicit none
   
   ! local variables
   logical :: has_failed,has_warned
   character(len=*), parameter :: rad_input_test_dir = 'data/'
   character(len=21), dimension(:), allocatable :: keyword_list     
   type(particle_radmat_type) :: particle_radmat_1
   integer :: file_unit_1,file_unit_2,file_unit_3,file_unit_4
      
   ! none of these tests use warnings
   has_warned = .false.
   
   ! open the input files that will be used
   call open_format_radmats_input(file_unit = file_unit_1, &
                                  file_name = rad_input_test_dir//'radiation_materials_read_input1.radmats')

   call open_format_radmats_input(file_unit = file_unit_2, &
                                  file_name = rad_input_test_dir//'radiation_materials_read_input2.radmats')

   call open_format_radmats_input(file_unit = file_unit_3, &
                                  file_name = rad_input_test_dir//'radiation_materials_read_input3.radmats')

   call open_format_radmats_input(file_unit = file_unit_4, &
                                  file_name = rad_input_test_dir//'radiation_materials_read_input4.radmats')
   
   call test_keyword_list_initialise_format_radmats(keyword_list) 
   
   call test_deduce_number_of_items_format_radmats(file_unit                 = file_unit_1, &
                                                   total_number_of_radmats   = 10,&
                                                   number_of_scatter_moments = 2, &
                                                   number_of_energy_groups   = 2, &
                                                   number_non_thermal_groups = 1, &
                                                   number_of_delayed_groups  = 6, &
                                                   record_len                = 132, &
                                                   keyword_list              = keyword_list)
   
   call test_deduce_number_of_items_format_radmats(file_unit                 = file_unit_2, &
                                                   total_number_of_radmats   = 10 ,&
                                                   number_of_scatter_moments = 2, &
                                                   number_of_energy_groups   = 6, &
                                                   number_non_thermal_groups = 3, &
                                                   number_of_delayed_groups  = 8, &
                                                   record_len                = 132, &
                                                   keyword_list              = keyword_list)
   
   call test_deduce_number_of_items_format_radmats(file_unit                 = file_unit_3, &
                                                   total_number_of_radmats   = 10 ,&
                                                   number_of_scatter_moments = 2, &
                                                   number_of_energy_groups   = 11, &
                                                   number_non_thermal_groups = 5, &
                                                   number_of_delayed_groups  = 0, &
                                                   record_len                = 132, &
                                                   keyword_list              = keyword_list)
   
   call test_deduce_number_of_items_format_radmats(file_unit                 = file_unit_4, &
                                                   total_number_of_radmats   = 10 ,&
                                                   number_of_scatter_moments = 2, &
                                                   number_of_energy_groups   = 172, &
                                                   number_non_thermal_groups = 92, &
                                                   number_of_delayed_groups  = 0, &
                                                   record_len                = 132, &
                                                   keyword_list              = keyword_list)

   ! create the particle_radmat's associated with the input files opened above
   call create_particle_radmat_input1(particle_radmat_1)
   
   call test_read_scatter_mom_format_radmats(file_unit            = file_unit_1, &
                                             record_len           = 132, &
                                             dataset_radmat       = particle_radmat_1%dataset_radmats(1), &
                                             particle_radmat_size = particle_radmat_1%particle_radmat_size, &
                                             keyword_list         = keyword_list) 
                                       
   call test_read_absorption_and_total_format_radmats(file_unit            = file_unit_1, &
                                                      record_len           = 132, &
                                                      dataset_radmat       = particle_radmat_1%dataset_radmats(1), &
                                                      particle_radmat_size = particle_radmat_1%particle_radmat_size, &
                                                      keyword_list         = keyword_list) 
   
   call test_calculate_removal_moments_format_radmats(file_unit            = file_unit_1, &
                                                      record_len           = 132, &
                                                      dataset_radmat       = particle_radmat_1%dataset_radmats(1), &
                                                      particle_radmat_size = particle_radmat_1%particle_radmat_size, &
                                                      keyword_list   = keyword_list) 
      
   call test_read_or_form_trans_and_diff_format_radmats(file_unit            = file_unit_1, &
                                                        record_len           = 132, &
                                                        dataset_radmat       = particle_radmat_1%dataset_radmats(1), &
                                                        particle_radmat_size = particle_radmat_1%particle_radmat_size, &
                                                        keyword_list         = keyword_list, &
                                                        problem_dim          = 3)

   call test_read_fission_data_format_radmats(file_unit                 = file_unit_1, &
                                              record_len                = 132, &
                                              dataset_radmat            = particle_radmat_1%dataset_radmats(1), &
                                              particle_radmat_size      = particle_radmat_1%particle_radmat_size, &
                                              keyword_list              = keyword_list, &
                                              number_non_thermal_groups = 1 ) 
         
   call test_read_velocity_style_data_format_radmats(file_unit           = file_unit_1, &
                                                     keyword_find        = (/16/), &
                                                     record_len          = 132, &
                                                     dataset_radmat      = particle_radmat_1%dataset_radmats(1), &
                                                     particle_radmat_size= particle_radmat_1%particle_radmat_size, &
                                                     keyword_list        = keyword_list)
   
   call test_read_velocity_style_data_format_radmats(file_unit           = file_unit_1, &
                                                     keyword_find        = (/17/), &
                                                     record_len          = 132, &
                                                     dataset_radmat      = particle_radmat_1%dataset_radmats(1), &
                                                     particle_radmat_size= particle_radmat_1%particle_radmat_size, &
                                                     keyword_list        = keyword_list)
   
   ! read the delayed data - beta, lambda and spectrum
   call test_read_delayed_data_format_radmats(file_unit               = file_unit_1, &
                                              record_len              = 132, &
                                              dataset_radmat          = particle_radmat_1%dataset_radmats(1), &
                                              particle_radmat_size    = particle_radmat_1%particle_radmat_size, &
                                              delayed_lambda_spectrum = particle_radmat_1%delayed_lambda_spectrum, &
                                              keyword_list            = keyword_list) 
   
   ! read the energy released through the power rotuine
   call test_read_power_data_format_radmats(file_unit            = file_unit_1, &
                                            record_len           = 132, &
                                            dataset_radmat       = particle_radmat_1%dataset_radmats(1), &
                                            particle_radmat_size = particle_radmat_1%particle_radmat_size, &
                                            test_power           = .false., &
                                            keyword_list         = keyword_list)

   call test_calculate_power_format_radmats(dataset_radmat       = particle_radmat_1%dataset_radmats(1), &
                                            particle_radmat_size = particle_radmat_1%particle_radmat_size) 
   
   ! get rid of the comparison particle_radmat's 
   call destroy(particle_radmat_1)   

   call test_read_format_radmats_siga_style_xsection(file_unit           = file_unit_1, &
                                                     keyword_line_number = 110, &
                                                     number_of_groups    = 2, &
                                                     record_len          = 132, &
                                                     values_expected     = (/8.110421E-06,&
                                                                           &1.643668E-04/), &
                                                     number_groups_to_read_in = 2)

   call test_read_format_radmats_siga_style_xsection(file_unit           = file_unit_2, &
                                                     keyword_line_number = 301, &
                                                     number_of_groups    = 6, &
                                                     record_len          = 132, &
                                                     values_expected     = (/1.319194E-01,&
                                                                           &2.411816E-01, &
                                                                           &3.394288E-01, &
                                                                           &3.457756E-01, &
                                                                           &0.0, &
                                                                           &0.0/), &
                                                     number_groups_to_read_in = 4) ! this is intentially 4, less than number of groups

   call test_read_format_radmats_siga_style_xsection(file_unit           = file_unit_3, &
                                                     keyword_line_number = 210, &
                                                     number_of_groups    = 11, &
                                                     record_len          = 132, &
                                                     values_expected     = (/2.738271E+00, &
                                                                            &2.499367E+00, &
                                                                            &2.442142E+00, &
                                                                            &2.437577E+00, &
                                                                            &2.432429E+00, &
                                                                            &2.432092E+00, &
                                                                            &2.436949E+00, &
                                                                            &2.440309E+00, &
                                                                            &2.440377E+00, &
                                                                            &2.438720E+00, &
                                                                            &0.0/), &
                                                     number_groups_to_read_in = 10) ! this is intentially 10, less than number of groups

   call test_read_format_radmats_siga_style_xsection(file_unit           = file_unit_4, &
                                                     keyword_line_number = 2, &
                                                     number_of_groups    = 172, &
                                                     record_len          = 132, &
                                                     values_expected     = (/ &
 & 8.067299E-03,  6.038252E-03,  5.409664E-03,  5.968075E-03,  7.794893E-03,  1.007881E-02,  3.085696E-03,  2.508856E-04, &
 & 1.246804E-04,  1.170927E-04,  1.101318E-04,  8.987960E-05,  8.768449E-05,  8.835270E-05,  8.829044E-05,  8.277333E-05, &
 & 6.385334E-05,  3.538791E-05,  3.413362E-05,  3.358710E-05,  3.340005E-05,  3.180470E-05,  3.018787E-05,  2.948743E-05, &
 & 2.919953E-05,  2.932990E-05,  2.985447E-05,  3.122614E-05,  3.146842E-05,  3.239765E-05,  3.448074E-05,  3.824649E-05, &
 & 4.158788E-05,  4.501959E-05,  5.346744E-05,  5.797070E-05,  6.766545E-05,  7.434823E-05,  8.332580E-05,  8.529446E-05, &
 & 8.688278E-05,  9.702778E-05,  1.058806E-04,  1.166753E-04,  1.229638E-04,  1.241919E-04,  1.515886E-04,  1.604488E-04, &
 & 1.735822E-04,  1.989856E-04,  2.330108E-04,  2.261882E-04,  2.872929E-04,  2.276664E-04,  3.208860E-04,  4.138338E-04, &
 & 4.372362E-04,  4.386278E-04,  4.670197E-04,  5.642580E-04,  5.043866E-04,  4.964821E-04,  1.004975E-03,  8.220389E-04, &
 & 5.121720E-04,  1.590378E-03,  1.040307E-03,  6.230862E-04,  2.889554E-03,  8.162664E-04,  1.125000E-03,  6.190391E-04, &
 & 7.096766E-04,  1.457212E-03,  1.163117E-02,  1.115339E-03,  4.095611E-04,  6.373236E-04,  1.285890E-03,  1.238334E-02, &
 & 1.469368E-03,  6.760540E-04,  1.679034E-03,  3.935889E-04,  9.094165E-04,  3.415513E-03,  3.908598E-04,  1.971280E-02, &
 & 9.017391E-04,  1.196401E-03,  5.350327E-04,  1.695218E-04,  6.623784E-04,  3.787110E-04,  5.305748E-04,  2.315603E-04, &
 & 2.075522E-04,  2.152663E-04,  2.265663E-04,  2.606361E-04,  3.793263E-04,  8.698702E-04,  5.660166E-04,  3.079374E-04, &
 & 2.928140E-04,  2.969397E-04,  3.061724E-04,  3.200321E-04,  3.338119E-04,  3.432857E-04,  3.665170E-04,  4.040914E-04, &
 & 4.503087E-04,  5.965551E-04,  1.095904E-03,  1.630850E-03,  1.811831E-03,  1.819887E-03,  1.752076E-03,  1.598265E-03, &
 & 1.387934E-03,  1.263386E-03,  1.192312E-03,  1.103408E-03,  1.042265E-03,  1.007934E-03,  9.667788E-04,  9.309679E-04, &
 & 9.063189E-04,  8.805951E-04,  8.701653E-04,  8.695509E-04,  8.781026E-04,  8.998969E-04,  9.680308E-04,  1.088646E-03, &
 & 1.223964E-03,  1.304494E-03,  1.437788E-03,  1.668983E-03,  1.826163E-03,  2.091639E-03,  2.571622E-03,  2.833519E-03, &
 & 2.959912E-03,  3.095748E-03,  3.067261E-03,  2.897151E-03,  2.822073E-03,  2.856051E-03,  2.935983E-03,  3.100944E-03, &
 & 3.239625E-03,  3.421451E-03,  3.722558E-03,  3.938990E-03,  4.221414E-03,  4.486343E-03,  4.736976E-03,  5.152412E-03, &
 & 5.621122E-03,  6.179927E-03,  6.846588E-03,  7.523764E-03,  8.249716E-03,  9.194903E-03,  1.049205E-02,  1.243735E-02, &
 & 1.504008E-02,  1.775471E-02,  2.132482E-02,  3.129871E-02/), &
                                                     number_groups_to_read_in = 172)

   call test_read_format_radmats_sigs_style_xsection(file_unit           = file_unit_2, &
                                                     keyword_line_number = 303, &
                                                     number_of_groups    = 6, &
                                                     record_len          = 132, &
                                                     values_expected     = reshape((/ &
!!!GROUP      1 FIRST      1 LAST      3
!!!  1.340728E-01  2.518163E-02  1.290233E-05
!!!GROUP      2 FIRST      2 LAST      3
!!!  2.400110E-01  3.143471E-02
!!!GROUP      3 FIRST      3 LAST      4
!!!  3.544767E-01  5.875821E-03
!!!GROUP      4 FIRST      3 LAST      6
!!!  3.001369E-05  3.363762E-01  2.930392E-02  2.042756E-05
!!!GROUP      5 FIRST      4 LAST      6
!!!  5.222967E-04  3.323453E-01  3.644267E-02
!!!GROUP      6 FIRST      4 LAST      6
!!! -1.882949E-13  1.958612E-03  3.808241E-01
& 1.340728E-01,0.000000E-00,0.000000E-00,0.000000E-00,0.000000E-00,0.000000E-00, &
& 2.518163E-02,2.400110E-01,0.000000E-00,0.000000E-00,0.000000E-00,0.000000E-00, &
& 1.290233E-05,3.143471E-02,3.544767E-01,3.001369E-05,0.000000E-00,0.000000E-00, &
& 0.000000E-00,0.000000E-00,5.875821E-03,3.363762E-01,5.222967E-04,-1.882949E-13, &
& 0.000000E-00,0.000000E-00,0.000000E-00,2.930392E-02,3.323453E-01,1.958612E-03, &
& 0.000000E-00,0.000000E-00,0.000000E-00,2.042756E-05,3.644267E-02,3.808241E-01 &
                                                                         & /),(/6,6/)), &
                                                      keyword_list       = keyword_list)

   call test_read_format_radmats_sigs_style_xsection(file_unit           = file_unit_3, &
                                                     keyword_line_number = 119, &
                                                     number_of_groups    = 11, &
                                                     record_len          = 132, &
                                                     values_expected     = reshape((/ &
!!!GROUP      1 FIRST      1 LAST      8
!!!  1.916469E-02 -6.137192E-03 -8.297635E-07 -3.124828E-08 -3.743999E-11  1.630492E-13  4.069419E-15  7.149059E-15
!!!GROUP      2 FIRST      2 LAST     10
!!!  2.799424E-02 -5.242567E-03 -7.767272E-10  7.369533E-11  4.031214E-13  1.055027E-14  2.202174E-14  2.569627E-15
!!!  4.603376E-16
!!!GROUP      3 FIRST      3 LAST      6
!!!  2.395448E-02 -4.710529E-03 -2.046192E-11 -6.321095E-16
!!!GROUP      4 FIRST      4 LAST      5
!!!  2.235033E-02 -4.565034E-03
!!!GROUP      5 FIRST      5 LAST      7
!!!  2.052691E-02 -2.820397E-03 -1.093021E-13
!!!GROUP      6 FIRST      5 LAST     11
!!!  2.866047E-05  2.828996E-02 -4.865639E-03 -5.883996E-03 -1.009163E-05 -1.606018E-07 -1.281481E-08
!!!GROUP      7 FIRST      6 LAST     11
!!!  6.954149E-04  5.089897E-02 -3.442162E-02 -3.814381E-04 -1.049360E-05 -3.559836E-08
!!!GROUP      8 FIRST      6 LAST     11
!!! -7.848687E-08  1.129392E-04  2.964792E-02 -1.178535E-02 -1.532755E-03 -4.955521E-05
!!!GROUP      9 FIRST      6 LAST     11
!!!  1.046262E-16 -8.379096E-14 -4.034751E-04  3.047878E-02 -1.343326E-02 -9.057373E-04
!!!GROUP     10 FIRST      8 LAST     11
!!! -2.080068E-07 -5.521337E-04  1.696245E-02 -5.232302E-03
!!!GROUP     11 FIRST      8 LAST     11
!!! -5.258472E-10 -3.853524E-06 -1.894755E-03 -5.618228E-03
&  1.916469E-02, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00,  &
& -6.137192E-03, 2.799424E-02, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00,  &
& -8.297635E-07,-5.242567E-03, 2.395448E-02, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00,  &
& -3.124828E-08,-7.767272E-10,-4.710529E-03, 2.235033E-02, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00,  &
& -3.743999E-11, 7.369533E-11,-2.046192E-11,-4.565034E-03, 2.052691E-02, 2.866047E-05, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00,  &
&  1.630492E-13, 4.031214E-13,-6.321095E-16, 0.000000E-00,-2.820397E-03, 2.828996E-02, 6.954149E-04,-7.848687E-08, 1.046262E-16, 0.000000E-00, 0.000000E-00,  &
&  4.069419E-15, 1.055027E-14, 0.000000E-00, 0.000000E-00,-1.093021E-13,-4.865639E-03, 5.089897E-02, 1.129392E-04,-8.379096E-14, 0.000000E-00, 0.000000E-00,  &
&  7.149059E-15, 2.202174E-14, 0.000000E-00, 0.000000E-00, 0.000000E-00,-5.883996E-03,-3.442162E-02, 2.964792E-02,-4.034751E-04,-2.080068E-07,-5.258472E-10,  &
&  0.000000E-00, 2.569627E-15, 0.000000E-00, 0.000000E-00, 0.000000E-00,-1.009163E-05,-3.814381E-04,-1.178535E-02, 3.047878E-02,-5.521337E-04,-3.853524E-06,  &
&  0.000000E-00, 4.603376E-16, 0.000000E-00, 0.000000E-00, 0.000000E-00,-1.606018E-07,-1.049360E-05,-1.532755E-03,-1.343326E-02, 1.696245E-02,-1.894755E-03,  &
&  0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00, 0.000000E-00,-1.281481E-08,-3.559836E-08,-4.955521E-05,-9.057373E-04,-5.232302E-03,-5.618228E-03  &
                                                                         & /),(/11,11/)), &
                                                     keyword_list       = keyword_list)

   call test_count_number_lines_with_keyword(file_unit              = file_unit_1, &
                                             record_len             = 132, &
                                             keyword_list           = keyword_list, &
                                             keyword_to_count       = 2, &
                                             keyword_to_stop        = 19, &
                                             keyword_count_expected = 10, &
                                             line_number_expected   = 199) 

   call test_find_line_with_any_keyword(file_unit                    = file_unit_1, &
                                        record_len                   = 132, &
                                        keyword_list                 = keyword_list, &
                                        first_keyword_found_expected = 1, &
                                        line_number_expected         = 1) 
   
   ! this will use file_unit_1
   call test_find_line_with_any_desired_keyword(record_len   = 132, &
                                                keyword_list = keyword_list) 
   
   call test_all_upper_case(test_string                    = 'Fluidity', &
                            test_string_all_upper_case_ref = 'FLUIDITY')
   
   call test_all_upper_case(test_string                    = 'Fluidity rAdiatiOn', &
                            test_string_all_upper_case_ref = 'FLUIDITY RADIATION')
   
   call test_make_character_upper_case(test_character                = 'a', &
                                       test_character_upper_case_ref = 'A')
   
   call test_make_character_upper_case(test_character                = 'b', &
                                       test_character_upper_case_ref = 'B')
  
   call test_make_character_upper_case(test_character                = 'c', &
                                       test_character_upper_case_ref = 'C')
   
   call test_make_character_upper_case(test_character                = 'd', &
                                       test_character_upper_case_ref = 'D')
   
   call test_make_character_upper_case(test_character                = 'e', &
                                       test_character_upper_case_ref = 'E')
   
   call test_make_character_upper_case(test_character                = 'f', &
                                       test_character_upper_case_ref = 'F')
   
   call test_make_character_upper_case(test_character                = 'g', &
                                       test_character_upper_case_ref = 'G')
   
   call test_make_character_upper_case(test_character                = 'h', &
                                       test_character_upper_case_ref = 'H')
   
   call test_make_character_upper_case(test_character                = 'i', &
                                       test_character_upper_case_ref = 'I')
   
   call test_make_character_upper_case(test_character                = 'j', &
                                       test_character_upper_case_ref = 'J')
   
   call test_make_character_upper_case(test_character                = 'k', &
                                       test_character_upper_case_ref = 'K')
   
   call test_make_character_upper_case(test_character                = 'l', &
                                       test_character_upper_case_ref = 'L')
   
   call test_make_character_upper_case(test_character                = 'm', &
                                       test_character_upper_case_ref = 'M')
   
   call test_make_character_upper_case(test_character                = 'n', &
                                       test_character_upper_case_ref = 'N')
   
   call test_make_character_upper_case(test_character                = 'o', &
                                       test_character_upper_case_ref = 'O')
   
   call test_make_character_upper_case(test_character                = 'p', &
                                       test_character_upper_case_ref = 'P')
   
   call test_make_character_upper_case(test_character                = 'q', &
                                       test_character_upper_case_ref = 'Q')
   
   call test_make_character_upper_case(test_character                = 'r', &
                                       test_character_upper_case_ref = 'R')
   
   call test_make_character_upper_case(test_character                = 's', &
                                       test_character_upper_case_ref = 'S')
   
   call test_make_character_upper_case(test_character                = 't', &
                                       test_character_upper_case_ref = 'T')
   
   call test_make_character_upper_case(test_character                = 'u', &
                                       test_character_upper_case_ref = 'U')
   
   call test_make_character_upper_case(test_character                = 'v', &
                                       test_character_upper_case_ref = 'V')
   
   call test_make_character_upper_case(test_character                = 'w', &
                                       test_character_upper_case_ref = 'W')
   
   call test_make_character_upper_case(test_character                = 'x', &
                                       test_character_upper_case_ref = 'X')
   
   call test_make_character_upper_case(test_character                = 'y', &
                                       test_character_upper_case_ref = 'Y')
   
   call test_make_character_upper_case(test_character                = 'z', &
                                       test_character_upper_case_ref = 'Z')
   
   ! the extra appended spaces in the character array constructor are necessary to compile - not for the actual test   
   call test_read_words_from_string(test_string    = 'A   random     sentence illustrating   a lack of imagination but   with too many spaces  3.14', &
                                    expected_words = (/'A           ', &
                                                       'random      ', &
                                                       'sentence    ', &
                                                       'illustrating', &
                                                       'a           ', &
                                                       'lack        ', &
                                                       'of          ', &
                                                       'imagination ', &
                                                       'but         ', &
                                                       'with        ', &
                                                       'too         ', &
                                                       'many        ', &
                                                       'spaces      ', &
                                                       '3.14        '/))

   call test_string_word_count(test_string              = 'A   random     sentence illustrating   a lack of imagination but   with too many spaces  3.14', &
                               expected_number_of_words = 14)
   
   call test_number_substrings_within_string(test_string       = 'A long time ago in a galaxy far, far away....', &
                                             test_substring    = 'far', &
                                             number_substrings = 2) 

   call test_number_substrings_within_string(test_string       = ' E+E-E- E  E+ E-EE E E E-E+E-', &
                                             test_substring    = 'E', &
                                             number_substrings = 13) 

   call test_number_substrings_within_string(test_string       = 'bish bash bosh', &
                                             test_substring    = 'sh', &
                                             number_substrings = 3) 
   
   ! this will use file_unit_1
   call test_read_next_line_seq(record_len = 132) 

   ! this will use file_unit_1
   call test_read_previous_line_seq(record_len = 132) 

   ! this will use file_unit_1
   call test_rewind_file_seq(record_len = 132) 

   ! this will use file_unit_1
   call test_go_to_line_seq(record_len = 132) 

   call test_keyword_list_cleanup(keyword_list) 
   
   close(file_unit_1)
   close(file_unit_2)
   close(file_unit_3)
   close(file_unit_4)

   !call test_desired_keyword_not_first_found() ! cannot be tested

contains

   ! --------------------------------------------------------------------------

   subroutine test_keyword_list_initialise_format_radmats(keyword_list)
      
      !!< Test the procedure that initialises the keyword list
      
      character(len=21), dimension(:), allocatable, intent(inout) :: keyword_list     
      
      ! local variables
      integer :: k
      character(len=21), dimension(22) :: keyword_list_check     
                   
      call keyword_list_initialise_format_radmats(keyword_list = keyword_list)
      
      ! check that the size of the keyword_list is correct
      first_check: if (size(keyword_list) /= 22) then
         
        has_failed = .true.
         
      else first_check

        has_failed = .false.
      
      end if first_check

      call report_test("[test_keyword_list_initialise_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed keyword_list wrong size")   
      
      ! check that each keyword is correct within the list
      keyword_list_check(1)  = 'MACRO'
      keyword_list_check(2)  = 'SIGABS'
      keyword_list_check(3)  = 'SIGTOT'
      keyword_list_check(4)  = 'SIGTRAN ' ! the space at the end is important
      keyword_list_check(5)  = 'SIGTRANY'
      keyword_list_check(6)  = 'SIGTRANZ'
      keyword_list_check(7)  = 'DIFFUSIONX'
      keyword_list_check(8)  = 'DIFFUSIONY' 
      keyword_list_check(9)  = 'DIFFUSIONZ' 
      keyword_list_check(10) = 'SIGFISS' 
      keyword_list_check(11) = 'FISSNU' 
      keyword_list_check(12) = 'PRODUCTION'       
      keyword_list_check(13) = 'FISSCHI' 
      keyword_list_check(14) = 'VELOCITY' 
      keyword_list_check(15) = 'POWER' 
      keyword_list_check(16) = 'ENERGY_RELEASED'       
      keyword_list_check(17) = 'BETA' 
      keyword_list_check(18) = 'LAMBDA' 
      keyword_list_check(19) = 'DELAYED_SPECTRUM' 
      keyword_list_check(20) = 'GROUP'  
      keyword_list_check(21) = 'DELAYED_GROUPS'  
      keyword_list_check(22) = 'SCATTER'
     
      keyword_loop: do k = 1,size(keyword_list_check)
      
         second_check: if (keyword_list(k) /= keyword_list_check(k)) then
         
           has_failed = .true.
         
         else second_check

           has_failed = .false.
      
         end if second_check

         call report_test("[test_keyword_list_initialise_format_radmats]", &
                         has_failed, &
                         has_warned, &
                         "failed keyword_list has incorrect element "//int2str(k))   
      
      end do keyword_loop
            
   end subroutine test_keyword_list_initialise_format_radmats
   
   ! --------------------------------------------------------------------------

   subroutine test_keyword_list_cleanup(keyword_list)
      
      !!< Test the procedure that deallocates the keyword list
      
      character(len=21), dimension(:), allocatable, intent(inout) :: keyword_list   
            
      call keyword_list_cleanup(keyword_list)    
      
      check: if (allocated(keyword_list)) then
      
         has_failed = .true.
         
         deallocate(keyword_list)
      
      else check
      
         has_failed = .false.
      
      end if check
      
      call report_test("[test_keyword_list_cleanup]", &
                       has_failed, &
                       has_warned, &
                       "failed")   
      
   end subroutine test_keyword_list_cleanup
   
   ! --------------------------------------------------------------------------

   subroutine test_deduce_number_of_items_format_radmats(file_unit, &
                                                         total_number_of_radmats, &
                                                         number_of_scatter_moments, &
                                                         number_of_energy_groups, &
                                                         number_non_thermal_groups, &
                                                         number_of_delayed_groups, &
                                                         record_len, &
                                                         keyword_list)
      
      !!< Test the procedure that deduces from a format_radmats input file the number of rad mats, the number of scatter moments, 
      !!< the number of energy groups, the number of non thermal energy groups and the number of delayed groups
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: total_number_of_radmats
      integer, intent(in) :: number_of_scatter_moments
      integer, intent(in) :: number_of_energy_groups
      integer, intent(in) :: number_non_thermal_groups
      integer, intent(in) :: number_of_delayed_groups
      integer, intent(in) :: record_len
      character(len=21), dimension(:), intent(in) :: keyword_list   
       
      ! local variables
      integer :: total_number_of_radmats_format_radmats
      integer :: number_of_energy_groups_format_radmats
      integer :: number_non_thermal_groups_format_radmats
      integer :: number_of_scatter_moments_format_radmats
      integer :: number_of_delayed_groups_format_radmats
                                                
      call deduce_number_of_items_format_radmats(file_unit, &
                                                 total_number_of_radmats_format_radmats, &
                                                 number_of_energy_groups_format_radmats, &
                                                 number_non_thermal_groups_format_radmats, &
                                                 number_of_scatter_moments_format_radmats, &
                                                 number_of_delayed_groups_format_radmats, &
                                                 record_len, &
                                                 keyword_list)
     
      ! first check
      check_number_of_energy_groups: if (number_of_energy_groups_format_radmats == number_of_energy_groups) then
            
         has_failed = .false.
                                            
      else check_number_of_energy_groups

         has_failed = .true.
                            
      end if check_number_of_energy_groups
            
      call report_test("[test_deduce_number_of_items_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong number of energy groups")   
            
      ! second check
      check_number_of_energy_thermal_groups: if (number_non_thermal_groups_format_radmats == number_non_thermal_groups) then
            
         has_failed = .false.
                                            
      else check_number_of_energy_thermal_groups

         has_failed = .true.
                            
      end if check_number_of_energy_thermal_groups
            
      call report_test("[test_deduce_number_of_items_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong number of thermal energy groups")   
            
      ! third check
      check_number_of_delayed_groups: if (number_of_delayed_groups_format_radmats == number_of_delayed_groups) then
            
         has_failed = .false.
                                            
      else check_number_of_delayed_groups

         has_failed = .true.
                            
      end if check_number_of_delayed_groups
            
      call report_test("[test_deduce_number_of_items_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong number of delayed groups")   

      ! fourth check
      check_number_of_scatter_moments: if (number_of_scatter_moments_format_radmats == number_of_scatter_moments) then
            
         has_failed = .false.
                                            
      else check_number_of_scatter_moments

         has_failed = .true.
                            
      end if check_number_of_scatter_moments
            
      call report_test("[test_deduce_number_of_items_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong number of scatter moments")   

      ! fourth check
      check_number_of_macro: if (total_number_of_radmats_format_radmats == total_number_of_radmats) then
            
         has_failed = .false.
                                            
      else check_number_of_macro

         has_failed = .true.
                            
      end if check_number_of_macro
            
      call report_test("[test_deduce_number_of_items_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong number of radiation materials")   
                                          
   end subroutine test_deduce_number_of_items_format_radmats
   
   ! --------------------------------------------------------------------------

   subroutine test_read_scatter_mom_format_radmats(file_unit, &
                                                   record_len, &
                                                   dataset_radmat, &
                                                   particle_radmat_size, &
                                                   keyword_list)
      
      !!< Test the procedure that reads in the scatter data numbers
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size      
      character(len=21), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      integer :: m
      integer :: pmat
      integer :: rmat
      integer :: line_number
      integer :: line_number_macro
      integer :: first_keyword_found
      logical :: exit_if_eof
      integer, dimension(:), allocatable :: keyword_find ! the specific keyword integer ids to find
      character(len=record_len) :: line_string
      type(dataset_radmat_type) :: dataset_radmat_format_radmats
                              
      ! intiialise the read in data type
      call allocate(dataset_radmat_format_radmats, &
                    particle_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_format_radmats)
      
      line_number = 0
      rewind(unit=file_unit)                  

      ! now check the read in values with the expected values to a default tolerance
      pmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)

            ! find the next MACRO keyword line
            allocate(keyword_find(1))
            keyword_find(1) = 1
            exit_if_eof = .false.
            call find_line_with_any_desired_keyword(file_unit, &
                                                    line_number, &
                                                    first_keyword_found, &
                                                    exit_if_eof, &
                                                    keyword_find, &
                                                    keyword_list, &
                                                    line_string)
            
            not_found: if (first_keyword_found /= keyword_find(1)) then
               
               has_failed = .true.
               
               call report_test("[test_read_scatter_mom_format_radmats]", &
                                has_failed, &
                                has_warned, &
                                "failed to find keyword MACRO as needed")   
                        
            end if not_found
            
            deallocate(keyword_find)
                                    
            ! keep the line number associated with this MACRO line
            line_number_macro = line_number
                          
            ! now read in the required values
            call read_scatter_mom_format_radmats(file_unit, &
                                                 line_number, &
                                                 dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat), &
                                                 record_len, &
                                                 keyword_list)
            
            ! check each moment matrix
            mom_loop: do m = 1,size(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%scatter,3)

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%scatter(:,:,m), &
                                      dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%scatter(:,:,m))

               if (has_failed) exit pmat_loop
            
            end do mom_loop

            ! go back to the line number of the current macro
            call go_to_line_seq(file_unit, &
                                line_number, &
                                go_to_line_number = line_number_macro, &
                                exit_if_eor       = .true.,&
                                exit_if_eof       = .true.) 
                                    
         end do rmat_loop
            
      end do pmat_loop
            
      call report_test("[test_read_scatter_mom_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_format_radmats)
      
   end subroutine test_read_scatter_mom_format_radmats
   
   ! --------------------------------------------------------------------------

   subroutine test_read_absorption_and_total_format_radmats(file_unit, &
                                                            record_len, &
                                                            dataset_radmat, &
                                                            particle_radmat_size, &
                                                            keyword_list)
      
      !!< Test the procedure that reads in the absorption or total data numbers. This will also implicitly read in the scatter data.
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size      
      character(len=21), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      integer :: pmat
      integer :: rmat
      integer :: line_number
      integer :: line_number_macro
      integer :: first_keyword_found
      logical :: exit_if_eof
      integer, dimension(:), allocatable :: keyword_find ! the specific keyword integer ids to find
      character(len=record_len) :: line_string
      type(dataset_radmat_type) :: dataset_radmat_format_radmats
                              
      ! intiialise the read in data type
      call allocate(dataset_radmat_format_radmats, &
                    particle_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_format_radmats)
      
      line_number = 0
      rewind(unit=file_unit)                  

      ! now check the read in values with the expected values to a default tolerance
      pmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)

            ! find the next MACRO keyword line
            allocate(keyword_find(1))
            keyword_find(1) = 1
            exit_if_eof = .false.
            call find_line_with_any_desired_keyword(file_unit, &
                                                    line_number, &
                                                    first_keyword_found, &
                                                    exit_if_eof, &
                                                    keyword_find, &
                                                    keyword_list, &
                                                    line_string)
    
            not_found: if (first_keyword_found /= keyword_find(1)) then
               
               has_failed = .true.
               
               call report_test("[test_read_absorption_and_total_format_radmats]", &
                                has_failed, &
                                has_warned, &
                                "failed to find keyword MACRO as needed")   
                        
            end if not_found
            
            deallocate(keyword_find)
                                    
            ! keep the line number associated with this MACRO line
            line_number_macro = line_number
                          
            ! now read in the required values
            call read_absorption_and_total_format_radmats(file_unit, &
                                                          line_number, &
                                                          dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat), &
                                                          record_len, &
                                                          keyword_list)
            
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%total, &
                                   dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%total)

            if (has_failed) exit pmat_loop
            
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%absorption, &
                                   dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%absorption)

            if (has_failed) exit pmat_loop
            
            ! go back to the line number of the current macro
            call go_to_line_seq(file_unit, &
                                line_number, &
                                go_to_line_number = line_number_macro, &
                                exit_if_eor       = .true.,&
                                exit_if_eof       = .true.) 
                                    
         end do rmat_loop
            
      end do pmat_loop
            
      call report_test("[test_read_absorption_and_total_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_format_radmats)
            
   end subroutine test_read_absorption_and_total_format_radmats
   
   ! --------------------------------------------------------------------------

   subroutine test_calculate_removal_moments_format_radmats(file_unit, &
                                                            record_len, &
                                                            dataset_radmat, &
                                                            particle_radmat_size, &
                                                            keyword_list)
      
      !!< Test the procedure that calculates the removal moments. This will also read in the total and absorption (which itself will 
      !!< read in the scatter data)
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size      
      character(len=21), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      integer :: m
      integer :: pmat
      integer :: rmat
      integer :: line_number
      integer :: line_number_macro
      integer :: first_keyword_found
      logical :: exit_if_eof
      integer, dimension(:), allocatable :: keyword_find ! the specific keyword integer ids to find
      character(len=record_len) :: line_string
      type(dataset_radmat_type) :: dataset_radmat_format_radmats
                              
      ! intiialise the read in data type
      call allocate(dataset_radmat_format_radmats, &
                    particle_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_format_radmats)
      
      line_number = 0
      rewind(unit=file_unit)                  

      ! now check the read in values with the expected values to a default tolerance
      pmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)

            ! find the next MACRO keyword line
            allocate(keyword_find(1))
            keyword_find(1) = 1
            exit_if_eof = .false.
            call find_line_with_any_desired_keyword(file_unit, &
                                                    line_number, &
                                                    first_keyword_found, &
                                                    exit_if_eof, &
                                                    keyword_find, &
                                                    keyword_list, &
                                                    line_string)
    
            not_found: if (first_keyword_found /= keyword_find(1)) then
               
               has_failed = .true.
               
               call report_test("[test_calculate_removal_moments_format_radmats]", &
                                has_failed, &
                                has_warned, &
                                "failed to find keyword MACRO as needed")   
                        
            end if not_found
            
            deallocate(keyword_find)
                                    
            ! keep the line number associated with this MACRO line
            line_number_macro = line_number
                          
            ! now read in the required values
            call calculate_removal_moments_format_radmats(file_unit, &
                                                          line_number, &
                                                          dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat), &
                                                          record_len, &
                                                          keyword_list)
            
            ! check each moment matrix
            mom_loop: do m = 1,size(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%removal,2)

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%removal(:,m), &
                                      dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%removal(:,m))

               if (has_failed) exit pmat_loop
            
            end do mom_loop

            if (has_failed) exit pmat_loop
                        
            ! go back to the line number of the current macro
            call go_to_line_seq(file_unit, &
                                line_number, &
                                go_to_line_number = line_number_macro, &
                                exit_if_eor       = .true.,&
                                exit_if_eof       = .true.) 
                                    
         end do rmat_loop
            
      end do pmat_loop
            
      call report_test("[test_calculate_removal_moments_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_format_radmats)
                  
   end subroutine test_calculate_removal_moments_format_radmats
   
   ! --------------------------------------------------------------------------

   subroutine test_read_or_form_trans_and_diff_format_radmats(file_unit, &
                                                              record_len, &
                                                              dataset_radmat, &
                                                              particle_radmat_size, &
                                                              keyword_list, &
                                                              problem_dim)
      
      !!< Test the procedure that reads the transport or diffusion data. This will also read in the total and absorption (which itself will 
      !!< read in the scatter data)
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size      
      character(len=21), dimension(:), intent(in) :: keyword_list
      integer, intent(in) :: problem_dim  
      
      ! local variables
      integer :: idim
      integer :: pmat
      integer :: rmat
      integer :: line_number
      integer :: line_number_macro
      integer :: first_keyword_found
      logical :: exit_if_eof
      integer, dimension(:), allocatable :: keyword_find ! the specific keyword integer ids to find
      character(len=record_len) :: line_string
      type(dataset_radmat_type) :: dataset_radmat_format_radmats
                              
      ! intiialise the read in data type
      call allocate(dataset_radmat_format_radmats, &
                    particle_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_format_radmats)
      
      line_number = 0
      rewind(unit=file_unit)                  

      ! now check the read in values with the expected values to a default tolerance
      pmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)

            ! find the next MACRO keyword line
            allocate(keyword_find(1))
            keyword_find(1) = 1
            exit_if_eof = .false.
            call find_line_with_any_desired_keyword(file_unit, &
                                                    line_number, &
                                                    first_keyword_found, &
                                                    exit_if_eof, &
                                                    keyword_find, &
                                                    keyword_list, &
                                                    line_string)
    
            not_found: if (first_keyword_found /= keyword_find(1)) then
               
               has_failed = .true.
               
               call report_test("[test_read_or_form_trans_and_diff_format_radmats]", &
                                has_failed, &
                                has_warned, &
                                "failed to find keyword MACRO as needed")   
                        
            end if not_found
            
            deallocate(keyword_find)
                                    
            ! keep the line number associated with this MACRO line
            line_number_macro = line_number
            
            ! now read in the required values for each dimension
            dim_loop: do idim = 1,problem_dim              

               call read_or_form_trans_and_diff_format_radmats(file_unit, &
                                                               line_number, &
                                                               dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat), &
                                                               idim, &
                                                               record_len, &
                                                               keyword_list)
            
               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%transport(:,idim), &
                                      dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%transport(:,idim))

               if (has_failed) exit pmat_loop

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%diffusion(:,idim), &
                                      dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%diffusion(:,idim))

               if (has_failed) exit pmat_loop

               ! go back to the line number of the current macro
               call go_to_line_seq(file_unit, &
                                   line_number, &
                                   go_to_line_number = line_number_macro, &
                                   exit_if_eor       = .true.,&
                                   exit_if_eof       = .true.) 
            
            end do dim_loop
                                                            
         end do rmat_loop
            
      end do pmat_loop
            
      call report_test("[test_read_or_form_trans_and_diff_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_format_radmats)
      
   end subroutine test_read_or_form_trans_and_diff_format_radmats
   
   ! --------------------------------------------------------------------------

   subroutine test_read_fission_data_format_radmats(file_unit, &
                                                    record_len, &
                                                    dataset_radmat, &
                                                    particle_radmat_size, &
                                                    keyword_list, &
                                                    number_non_thermal_groups)
      
      !!< Test the procedure that reads the fission data. 
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size      
      character(len=21), dimension(:), intent(in) :: keyword_list
      integer, intent(in) :: number_non_thermal_groups  
      
      ! local variables
      integer :: pmat
      integer :: rmat
      integer :: line_number
      integer :: line_number_macro
      integer :: first_keyword_found
      logical :: exit_if_eof
      integer, dimension(:), allocatable :: keyword_find ! the specific keyword integer ids to find
      character(len=record_len) :: line_string
      type(dataset_radmat_type) :: dataset_radmat_format_radmats
                              
      ! intiialise the read in data type
      call allocate(dataset_radmat_format_radmats, &
                    particle_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_format_radmats)
      
      line_number = 0
      rewind(unit=file_unit)                  

      ! now check the read in values with the expected values to a default tolerance
      pmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)

            ! find the next MACRO keyword line
            allocate(keyword_find(1))
            keyword_find(1) = 1
            exit_if_eof = .false.
            call find_line_with_any_desired_keyword(file_unit, &
                                                    line_number, &
                                                    first_keyword_found, &
                                                    exit_if_eof, &
                                                    keyword_find, &
                                                    keyword_list, &
                                                    line_string)
    
            not_found: if (first_keyword_found /= keyword_find(1)) then
               
               has_failed = .true.
               
               call report_test("[test_read_fission_data_format_radmats]", &
                                has_failed, &
                                has_warned, &
                                "failed to find keyword MACRO as needed")   
                        
            end if not_found
            
            deallocate(keyword_find)
                                    
            ! keep the line number associated with this MACRO line
            line_number_macro = line_number
                          
            ! now read in the required values
            call read_fission_data_format_radmats(file_unit, &
                                                  line_number, &
                                                  dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat), &
                                                  number_non_thermal_groups, &
                                                  record_len, &
                                                  keyword_list)
            
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%fission, &
                                   dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%fission)

            if (has_failed) exit pmat_loop
            
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%production, &
                                   dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%production)

            if (has_failed) exit pmat_loop
            
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%particle_released_per_fission, &
                                   dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%particle_released_per_fission)

            if (has_failed) exit pmat_loop
            
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%prompt_spectrum, &
                                   dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%prompt_spectrum)

            if (has_failed) exit pmat_loop
                                    
            ! go back to the line number of the current macro
            call go_to_line_seq(file_unit, &
                                line_number, &
                                go_to_line_number = line_number_macro, &
                                exit_if_eor       = .true.,&
                                exit_if_eof       = .true.) 
                                    
         end do rmat_loop
            
      end do pmat_loop
            
      call report_test("[test_read_fission_data_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_format_radmats)
      
   end subroutine test_read_fission_data_format_radmats
   
   ! --------------------------------------------------------------------------

   subroutine test_read_velocity_style_data_format_radmats(file_unit, &
                                                           keyword_find, &
                                                           record_len, &
                                                           dataset_radmat, &
                                                           particle_radmat_size, &
                                                           keyword_list)
      
      !!< Test the procedure that reads in a velocity style data numbers
      
      integer, intent(in) :: file_unit
      integer, dimension(1), intent(in) :: keyword_find
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size
      character(len=21), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      integer :: rmat
      integer :: pmat
      integer :: number_of_groups
      type(dataset_radmat_type) :: dataset_radmat_format_radmats
          
      ! intiialise the read in data type
      call allocate(dataset_radmat_format_radmats, &
                    particle_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_format_radmats)
            
      number_of_groups_if: if (keyword_find(1) == 17) then
            
         number_of_groups = particle_radmat_size%number_of_delayed_groups
            
      else number_of_groups_if
               
         number_of_groups = particle_radmat_size%number_of_energy_groups
            
      end if number_of_groups_if
            
      ! now read in the required values
      call read_velocity_style_data_format_radmats(file_unit, &
                                                   dataset_radmat_format_radmats, &
                                                   number_of_groups, &
                                                   keyword_find, &
                                                   record_len, &
                                                   keyword_list)

      ! now check the read in values with the expected values to a default tolerance
      pmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)
              
            check_values: if (keyword_find(1) == 14) then
            
               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%velocity, &
                                      dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%velocity)
                                                              
            else if (keyword_find(1) == 15) then 

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%power, &
                                      dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%power)

            else if (keyword_find(1) == 16) then 

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission, &
                                      dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission)
                                        
            else if (keyword_find(1) == 17) then 

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%beta, &
                                      dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%beta)

            end if check_values

            if (has_failed) exit pmat_loop
            
         end do rmat_loop
            
      end do pmat_loop
            
      call report_test("[test_read_velocity_style_data_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_format_radmats)
            
   end subroutine test_read_velocity_style_data_format_radmats
   
   ! --------------------------------------------------------------------------

   subroutine test_read_delayed_data_format_radmats(file_unit, &
                                                    record_len, &
                                                    dataset_radmat, &
                                                    particle_radmat_size, &
                                                    delayed_lambda_spectrum, &
                                                    keyword_list)
      
      !!< Test the procedure that reads in the delayed data numbers
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size      
      type(delayed_lambda_spectrum_type), intent(in) :: delayed_lambda_spectrum      
      character(len=21), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      integer :: pmat
      integer :: rmat
      type(dataset_radmat_type) :: dataset_radmat_format_radmats
      type(delayed_lambda_spectrum_type) :: delayed_lambda_spectrum_format_radmats
      
      ! intiialise the read in data types
      call allocate(dataset_radmat_format_radmats, &
                    particle_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_format_radmats)
 
      call allocate(delayed_lambda_spectrum_format_radmats, &
                    particle_radmat_size%number_of_delayed_groups, &
                    particle_radmat_size%number_of_energy_groups)

      call zero(delayed_lambda_spectrum_format_radmats)
                        
      ! now read in the required values
      call read_delayed_data_format_radmats(file_unit, &
                                            dataset_radmat_format_radmats, &
                                            record_len, &
                                            keyword_list, &
                                            delayed_lambda_spectrum_format_radmats, &
                                            read_delayed_lambda_spectrum = .true.)

      ! now check the read in values with the expected values to a default tolerance
            
      pmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)
                          
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%beta, &
                                   dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%beta)

            if (has_failed) exit pmat_loop
                                    
         end do rmat_loop
            
      end do pmat_loop
            
      call report_test("[test_read_delayed_data_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "beta failed values read in not what expected,")   
            
      has_failed =  fnequals(delayed_lambda_spectrum%lambda, &
                             delayed_lambda_spectrum_format_radmats%lambda)

      call report_test("[test_read_delayed_data_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "lambda failed values read in not what expected")   
            
      has_failed =  fnequals(delayed_lambda_spectrum%spectrum, &
                             delayed_lambda_spectrum_format_radmats%spectrum)
            
      call report_test("[test_read_delayed_data_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "spectrum failed values read in not what expected")   
                        
      call destroy(dataset_radmat_format_radmats)
      call destroy(delayed_lambda_spectrum_format_radmats)
                  
   end subroutine test_read_delayed_data_format_radmats
   
   ! --------------------------------------------------------------------------

   subroutine test_read_power_data_format_radmats(file_unit, &
                                                  record_len, &
                                                  dataset_radmat, &
                                                  particle_radmat_size, &
                                                  test_power, &
                                                  keyword_list)
      
      !!< Test the procedure that reads in the power data numbers
      !!< The power variable may not be read in but formed from the fission and energy released
      !!< but as this routine will not read the fission we may not always want to test power (only energy released)
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size      
      logical, intent(in) :: test_power
      character(len=21), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      integer :: pmat
      integer :: rmat
      type(dataset_radmat_type) :: dataset_radmat_format_radmats
                              
      ! intiialise the read in data type
      call allocate(dataset_radmat_format_radmats,particle_radmat_size,dmat=1)
      call zero(dataset_radmat_format_radmats)
                        
      ! now read in the required values
      call read_power_data_format_radmats(file_unit, &
                                          dataset_radmat_format_radmats, &
                                          record_len, &
                                          keyword_list)

      ! now check the read in values with the expected values to a default tolerance
      pmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)
                          
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission, &
                                   dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission)

            if (has_failed) exit pmat_loop
               
            check_power: if (test_power) then                                            

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%power, &
                                      dataset_radmat_format_radmats%physical_radmats(pmat)%radmats(rmat)%power)

               if (has_failed) exit pmat_loop
                  
            end if check_power
                  
         end do rmat_loop
            
      end do pmat_loop
            
      call report_test("[test_read_power_data_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_format_radmats)
      
   end subroutine test_read_power_data_format_radmats
   
   ! --------------------------------------------------------------------------

   subroutine test_calculate_power_format_radmats(dataset_radmat, &
                                                  particle_radmat_size)
      
      !!< Test the procedure that calculates the power data from the fission and energy released
      
      type(dataset_radmat_type), intent(inout) :: dataset_radmat
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size      
      
      ! local variables
      type(dataset_radmat_type) :: dataset_radmat_tmp
      integer :: pmat,rmat
      
      ! intiialise the tmp data type
      call allocate(dataset_radmat_tmp,particle_radmat_size,dmat=1)
      call zero(dataset_radmat_tmp)
      
      ! set tmp
      dataset_radmat_tmp = dataset_radmat
      
      ! zero the power            
      pmat_loop_zero: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop_zero: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)
               
            dataset_radmat%physical_radmats(pmat)%radmats(rmat)%power = 0.0
                         
         end do rmat_loop_zero
            
      end do pmat_loop_zero
      
      ! calc the power
      call calculate_power_format_radmats(dataset_radmat)
      
      ! check the power
      pmat_loop_check: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop_check: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)
               
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%power, &
                                   dataset_radmat_tmp%physical_radmats(pmat)%radmats(rmat)%power)

            if (has_failed) exit pmat_loop_check
                         
         end do rmat_loop_check
            
      end do pmat_loop_check

      call report_test("[test_calculate_power_format_radmats]", &
                       has_failed, &
                       has_warned, &
                       "failed values calculated not what expected for pmat "//int2str(pmat)//" for rmat "//int2str(rmat))   
            
      call destroy(dataset_radmat_tmp)
      
   end subroutine test_calculate_power_format_radmats
   
   ! --------------------------------------------------------------------------

   subroutine test_read_format_radmats_siga_style_xsection(file_unit, &
                                                           keyword_line_number, &
                                                           number_of_groups, &
                                                           record_len, &
                                                           values_expected, &
                                                           number_groups_to_read_in)
      
      !!< Test the procedure that reads in the absorption (or SIGA) style cross section values
      !!< This assumes that the keyword is at the keyword_line_number which is read to first
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: keyword_line_number 
      integer, intent(in) :: number_of_groups
      integer, intent(in) :: record_len  
      real, dimension(number_of_groups) :: values_expected   
      integer, intent(in), optional :: number_groups_to_read_in
      
      ! local variables
      integer :: line_number
      real, dimension(:), allocatable :: values_format_radmats
                  
      line_number = 0
      rewind(unit=file_unit)
            
      ! go to the desired keyword line
      call go_to_line_seq(file_unit, &
                          line_number, &
                          go_to_line_number=keyword_line_number, &
                          exit_if_eor=.false., &
                          exit_if_eof=.false.)
            
      ! intiialise the values read in array
      allocate(values_format_radmats(number_of_groups))
      values_format_radmats = 0.0
            
      ! now read in the required values
      call read_format_radmats_siga_style_xsection(values_format_radmats, &
                                                   line_number, &
                                                   file_unit, &
                                                   record_len, &
                                                   number_groups_to_read_in=number_groups_to_read_in)

      ! now check the read in values with the expected values to a default tolerance
      check_values: if (fequals(values_expected,values_format_radmats)) then
            
         has_failed = .false.
                                            
      else check_values

         has_failed = .true.
                            
      end if check_values
            
      call report_test("[test_read_format_radmats_siga_style_xsection]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      deallocate(values_format_radmats)                                  
      
   end subroutine test_read_format_radmats_siga_style_xsection
   
   ! --------------------------------------------------------------------------

   subroutine test_read_format_radmats_sigs_style_xsection(file_unit, &
                                                           keyword_line_number, &
                                                           number_of_groups, &
                                                           record_len, &
                                                           values_expected, &
                                                           keyword_list)
      
      !!< Test the procedure that reads in the scatter moment style cross section values
      !!< This read in assumes that the current line position is at the line with the keyword (eg. SCATTER MOMENT      1)
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: keyword_line_number 
      integer, intent(in) :: number_of_groups
      integer, intent(in) :: record_len  
      real, dimension(number_of_groups,number_of_groups) :: values_expected   
      character(len=21), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      integer :: line_number
      real, dimension(:,:), allocatable :: values_format_radmats
                        
      line_number = 0
      rewind(unit=file_unit)
            
      ! go to the desired keyword line
      call go_to_line_seq(file_unit, &
                          line_number, &
                          go_to_line_number=keyword_line_number, &
                          exit_if_eor=.false., &
                          exit_if_eof=.false.)
            
      ! intiialise the values read in array
      allocate(values_format_radmats(number_of_groups,number_of_groups))
      values_format_radmats = 0.0
            
      ! now read in the required values
      call read_format_radmats_sigs_style_xsection(values_format_radmats, &
                                                   line_number, &
                                                   file_unit, &
                                                   record_len, &
                                                   keyword_list)

      ! now check the read in values with the expected values to a default tolerance
      check_values: if (fequals(values_expected,values_format_radmats)) then
            
         has_failed = .false.
                                           
      else check_values

         has_failed = .true.
                            
      end if check_values
            
      call report_test("[test_read_format_radmats_sigs_style_xsection]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      deallocate(values_format_radmats)
            
   end subroutine test_read_format_radmats_sigs_style_xsection
   
   ! --------------------------------------------------------------------------

   subroutine test_count_number_lines_with_keyword(file_unit, &
                                                   record_len, &
                                                   keyword_list, &
                                                   keyword_to_count, &
                                                   keyword_to_stop, &
                                                   keyword_count_expected, &
                                                   line_number_expected)
      
      !!< Test the procedure that counts the number of lines read of a sequential file that contain a keyword
      !!< and stops counting if eof or find a line with a stopping keyword
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len
      character(len=21), dimension(:), intent(in) :: keyword_list   
      integer, intent(in) :: keyword_to_count
      integer, intent(in) :: keyword_to_stop
      integer, intent(in) :: keyword_count_expected
      integer, intent(in) :: line_number_expected
      
      ! local variables
      integer :: keyword_count
      integer :: line_number
                  
      line_number = 0
      rewind(unit=file_unit)
                              
      call count_number_lines_with_keyword(file_unit, &
                                           line_number, &
                                           keyword_count, &
                                           keyword_to_count, &
                                           record_len, &
                                           keyword_list, &
                                           keyword_to_stop=keyword_to_stop)                                                  

      check_keyword_count: if (keyword_count == keyword_count_expected) then
            
         has_failed = .false.
                                            
      else check_keyword_count

         has_failed = .true.
                            
      end if check_keyword_count
            
      call report_test("[test_count_number_keyword]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong number of keyword count")   
                        
      check_correct_line_number: if (line_number == line_number_expected) then
            
         has_failed = .false.
            
      else check_correct_line_number
            
         has_failed = .true.
            
      end if check_correct_line_number

      call report_test("[test_count_number_keyword]", &
                       has_failed, &
                       has_warned, &
                       "failed as keyword count finished on wrong line number")                                                                     
           
   end subroutine test_count_number_lines_with_keyword
   
   ! --------------------------------------------------------------------------

   subroutine test_find_line_with_any_keyword(file_unit, &
                                              record_len, &
                                              keyword_list, &
                                              first_keyword_found_expected, &
                                              line_number_expected)
            
      !!< Test the procedure that reads a sequential file until finding a line with any keyword
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len
      character(len=21), dimension(:), intent(in) :: keyword_list   
      integer, intent(in) :: first_keyword_found_expected
      integer, intent(in) :: line_number_expected
      
      ! local variables
      integer :: line_number
      integer :: first_keyword_found
      logical :: exit_if_eof
      character(len=record_len) :: line_string
            
      line_number = 0
      rewind(unit=file_unit)
            
      ! dont exit if end of file
      exit_if_eof = .false.
            
      ! first search for word one - being MACRO in the format_radmats list            
      call find_line_with_any_keyword(file_unit, &
                                      line_number, &
                                      first_keyword_found, &
                                      exit_if_eof, &
                                      keyword_list, &
                                      line_string)            

      check_first_keyword_found: if (first_keyword_found == first_keyword_found_expected) then
            
         has_failed = .false.
                                            
      else check_first_keyword_found

         has_failed = .true.
                            
      end if check_first_keyword_found
            
      call report_test("[test_find_line_with_any_keyword]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong first keyword found first or end of file")   
                        
      check_correct_line_number: if (line_number == line_number_expected) then
            
         has_failed = .false.
            
      else check_correct_line_number
            
         has_failed = .true.
            
      end if check_correct_line_number

      call report_test("[test_find_line_with_any_keyword]", &
                       has_failed, &
                       has_warned, &
                       "failed as first keyword found with wrong line number")                                                            
      
   end subroutine test_find_line_with_any_keyword
   
   ! --------------------------------------------------------------------------

   subroutine test_find_line_with_any_desired_keyword(record_len, &
                                                      keyword_list)
      
      !!< Test the procedure that reads a sequential access file until finding a line with a desired keyword in file_unit_1

      integer, intent(in) :: record_len
      character(len=21), dimension(:), intent(in) :: keyword_list   
      
      ! local variables
      integer :: line_number
      integer :: first_keyword_found
      logical :: exit_if_eof
      integer, dimension(:), allocatable  :: keyword_find
      character(len=record_len) :: line_string
                  
      line_number = 0
      rewind(unit=file_unit_1)
            
      ! dont exit if end of file
      exit_if_eof = .false.
            
      ! first search for word one - being MACRO in the format_radmats list
      allocate(keyword_find(1))
      keyword_find(1) = 1  
            
      call find_line_with_any_desired_keyword(file_unit_1, &
                                              line_number, &
                                              first_keyword_found, &
                                              exit_if_eof, &
                                              keyword_find, &
                                              keyword_list, &
                                              line_string)            

      check_first_keyword_found: if (first_keyword_found == 1) then 
            
         has_failed = .false.
                            
      else check_first_keyword_found

         has_failed = .true.
                            
      end if check_first_keyword_found
            
      call report_test("[test_find_line_with_any_desired_keyword]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong first keyword found first or end of file")   
                        
      check_correct_line_number: if (line_number == 1) then
            
         has_failed = .false.
            
      else check_correct_line_number
            
         has_failed = .true.
            
      end if check_correct_line_number

      call report_test("[test_find_line_with_any_desired_keyword]", &
                       has_failed, &
                       has_warned, &
                       "failed as first keyword found with wrong line number")   
            
      deallocate(keyword_find)  
            
      ! now search for three keywords
      allocate(keyword_find(3))            
      keyword_find(1) = 1  
      keyword_find(2) = 22  
      keyword_find(3) = 10  
            
      call find_line_with_any_desired_keyword(file_unit_1, &
                                              line_number, &
                                              first_keyword_found, &
                                              exit_if_eof, &
                                              keyword_find, &
                                              keyword_list, &
                                              line_string)            

      check_second_keyword_found: if (first_keyword_found == 2) then 
            
         has_failed = .false.
                            
      else check_second_keyword_found

         has_failed = .true.
                            
      end if check_second_keyword_found
            
      call report_test("[test_find_line_with_any_desired_keyword]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong keyword found first or end of file")   
                        
      check_correct_line_number_again: if (line_number == 6) then
            
         has_failed = .false.
            
      else check_correct_line_number_again
            
         has_failed = .true.
          
      end if check_correct_line_number_again

      call report_test("[test_find_line_with_any_desired_keyword]", &
                       has_failed, &
                       has_warned, &
                       "failed as keyword found with wrong line number")   
            
      deallocate(keyword_find)  
                                            
   end subroutine test_find_line_with_any_desired_keyword
   
   ! --------------------------------------------------------------------------

   subroutine test_all_upper_case(test_string, &
                                  test_string_all_upper_case_ref)
      
      !!< Test the procedure that returns the all upper case version of a string
      
      character(len=*), intent(in) :: test_string
      character(len=*), intent(in) :: test_string_all_upper_case_ref
      
      ! local variables
      character(len=len(test_string)) :: test_string_all_upper_case 
            
      ! form the upper case
      test_string_all_upper_case = test_string
      call all_upper_case(test_string_all_upper_case)

      first_check: if (trim(test_string_all_upper_case) == trim(test_string_all_upper_case_ref)) then
                       
         has_failed = .false.
      
      else first_check
      
         has_failed = .true.
      
      end if first_check
   
      call report_test("[test_all_upper_case]", &
                      has_failed, &
                      has_warned, &
                      "failed for test_string "//trim(test_string))   
      
   end subroutine test_all_upper_case

   ! --------------------------------------------------------------------------

   subroutine test_make_character_upper_case(test_character, &
                                             test_character_upper_case_ref)
      
      !!< Test the procedure that returns the all upper case version of a character
      
      character(len=1), intent(in) :: test_character
      character(len=1), intent(in) :: test_character_upper_case_ref
      
      ! local variables
      character(len=1) ::  test_character_upper_case 
      
      ! form the upper case
      test_character_upper_case = test_character
      call make_character_upper_case( test_character_upper_case)

      first_check: if (trim( test_character_upper_case) == trim(test_character_upper_case_ref)) then
                       
         has_failed = .false.
      
      else first_check
      
         has_failed = .true.
      
      end if first_check
   
      call report_test("[make_character_upper_case]", &
                      has_failed, &
                      has_warned, &
                      "failed for test_character "//trim(test_character))   
      
   end subroutine test_make_character_upper_case

   ! --------------------------------------------------------------------------

   subroutine test_read_words_from_string(test_string, &
                                          expected_words)
      
      !!< Test the procedure that reads the words from a string
      
      character(len=*), intent(in) :: test_string
      character(len=*), dimension(:), intent(in) :: expected_words
      
      ! local variables
      integer :: w
      character(len=len(test_string)), dimension(:), allocatable :: words 
            
      ! read the words
      call read_words_from_string(test_string, &
                                  words)
      
      ! check that there is the right number of words
      first_check: if (size(expected_words) == size(words)) then
                       
         has_failed = .false.
      
      else first_check
      
         has_failed = .true.
      
      end if first_check
   
      call report_test("[test_read_words_from_string]", &
                       has_failed, &
                       has_warned, &
                       "failed as incorrect number of words for test_string "//trim(test_string))   
      
      ! check each word with the expected
      word_loop: do w = 1,size(words)
      
         second_check: if (trim(expected_words(w)) == trim(words(w))) then
                       
            has_failed = .false.
      
         else second_check
      
            has_failed = .true.
      
         end if second_check
   
         call report_test("[test_read_words_from_string]", &
                          has_failed, &
                          has_warned, &
                          "failed for test_string "//trim(test_string)//" for word number "//int2str(w))   
      
      end do word_loop
         
   end subroutine test_read_words_from_string

   ! --------------------------------------------------------------------------

   subroutine test_string_word_count(test_string, &
                                     expected_number_of_words)
      
      !!< Test the procedure that counts the number of words within a string
      
      character(len=*), intent(in) :: test_string
      integer, intent(in) ::expected_number_of_words 
            
      ! local variables
      integer :: number_of_words 
            
      ! count the words
      call string_word_count(test_string, &
                             number_of_words)
      
      ! check that there is the right number of words
      first_check: if (expected_number_of_words == number_of_words) then
                       
         has_failed = .false.
      
      else first_check
      
         has_failed = .true.
      
      end if first_check
   
      call report_test("[test_string_word_count]", &
                      has_failed, &
                      has_warned, &
                      "failed as incorrect number of words for test_string "//trim(test_string))   
               
   end subroutine test_string_word_count

   ! --------------------------------------------------------------------------

   subroutine test_number_substrings_within_string(test_string,  &
                                                   test_substring, &
                                                   number_substrings)
      
      !!< Test the procedure that returns the number of substrings within a string 
      
      character(len=*), intent(in) :: test_string,test_substring  
      integer, intent(in) :: number_substrings
                   
      first_check: if (number_substrings_within_string(trim(test_string),trim(test_substring)) == number_substrings) then
                       
         has_failed = .false.
      
      else first_check
      
         has_failed = .true.
      
      end if first_check
   
      call report_test("[test_number_substrings_within_string]", &
                      has_failed, &
                      has_warned, &
                      "failed for test_string "//trim(test_string)//" and subtring "//trim(test_substring))   
      
   end subroutine test_number_substrings_within_string   

   ! --------------------------------------------------------------------------

   subroutine test_read_next_line_seq(record_len)
      
      !!< Test the procedure that reads the next line of a sequential file and increases the line number
      !!< count by one and checks if eof and eor and should exit if this happens. This will use file_unit_1.

      integer, intent(in) :: record_len
      
      ! local variables
      integer :: line
      integer :: line_number
      logical :: end_of_file
      character(len=record_len) :: line_string
                  
      line_number = 0
      rewind(unit=file_unit_1)
            
      ! read some lines
      line_loop1: do line = 1,146
               
         call read_next_line_seq(file_unit_1, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor = .false., &
                                 exit_if_eof = .false.)
            
      end do line_loop1
            
      ! test that we have the correct line            
      check_correct_line1: if (trim(line_string) == '  3.608233E-01  4.655511E-03') then
               
         has_failed = .false.
            
      else check_correct_line1
            
         has_failed = .true.
            
      end if check_correct_line1
            
      call report_test("[test_read_next_line_seq]", &
                       has_failed, &
                       has_warned, &
                       "failed as incorrect first line_string read")   
            
      ! read some more lines
      line_loop2: do line = 1,5
            
         call read_next_line_seq(file_unit_1, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor = .false., &
                                 exit_if_eof = .false.)
            
      end do line_loop2
            
      ! test that we have the correct line            
      check_correct_line2: if (trim(line_string) == '  2.453666E-02 -1.374331E-03') then
               
         has_failed = .false.
            
      else check_correct_line2
            
         has_failed = .true.
            
      end if check_correct_line2
            
      call report_test("[test_read_next_line_seq]", &
                       has_failed, &
                       has_warned, &
                       "failed as incorrect second line_string read")   

      ! read some more lines
      line_loop3: do line = 1,55
            
         call read_next_line_seq(file_unit_1, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor = .false., &
                                 exit_if_eof = .false.)
           
      end do line_loop3
            
      ! test we have the end of file
      check_eof: if (end_of_file) then
               
         has_failed = .false.
            
      else check_eof
            
         has_failed = .true.
            
      end if check_eof
            
      call report_test("[test_read_next_line_seq]", &
                       has_failed, &
                       has_warned, &
                       "failed as not end of file")   
                                       
   end subroutine test_read_next_line_seq
   
   ! --------------------------------------------------------------------------

   subroutine test_read_previous_line_seq(record_len)
      
      !!< Test the procedure that reads the previous line of a sequential file and decrease the line number
      !!< count by one and checks if eof and eor and should exit if this happens
      !!< such as to read backwards we must read forwards a few times to begin with. This will use file_unit_1.
      
      integer, intent(in) :: record_len
      
      ! local variables
      integer :: line
      integer :: line_number
      logical :: end_of_file
      character(len=record_len) :: line_string
                    
      line_number = 0
      rewind(unit=file_unit_1)
            
      ! read some lines forward
      line_loop1: do line = 1,80
               
         call read_next_line_seq(file_unit_1, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor = .false., &
                                 exit_if_eof = .false.)
            
      end do line_loop1
            
      ! read some lines backward
      line_loop2: do line = 1,4
               
         call read_previous_line_seq(file_unit_1, &
                                     line_string, &
                                     line_number, &
                                     exit_if_eor = .false., &
                                     exit_if_eof = .false.)
            
      end do line_loop2
          
      ! test that we have the correct line            
      check_correct_line1: if (trim(line_string) == '  2.527262E-02 -1.623783E-03') then
               
         has_failed = .false.
            
      else check_correct_line1
            
         has_failed = .true.
            
      end if check_correct_line1
            
      call report_test("[test_read_previous_line_seq]", &
                       has_failed, &
                       has_warned, &
                       "failed as incorrect first line_string read")   
                                          
   end subroutine test_read_previous_line_seq
   
   ! --------------------------------------------------------------------------

   subroutine test_rewind_file_seq(record_len)
      
      !!< Test the procedure that rewinds a sequential file and changes the line number
      !!< via first opening then reading a few lines. This will use file_unit_1.

      integer, intent(in) :: record_len
            
      ! local variables
      integer :: line
      integer :: line_number
      logical :: end_of_file
      character(len=record_len) :: line_string
                  
      line_number = 0
      rewind(unit=file_unit_1)
            
      ! read some lines forward
      line_loop1: do line = 1,122
               
         call read_next_line_seq(file_unit_1, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor = .false., &
                                 exit_if_eof = .false.)
            
      end do line_loop1
            
      ! rewind file
      call rewind_file_seq(file_unit_1, &
                           line_number) 
            
      ! read next line
      call read_next_line_seq(file_unit_1, &
                              line_string, &
                              line_number, &
                              end_of_file, &
                              exit_if_eor = .false., &
                              exit_if_eof = .false.)
        
      ! test that we have the correct line            
      check_correct_line1: if (trim(line_string) == 'MACRO      1') then
               
         has_failed = .false.
            
      else check_correct_line1
            
         has_failed = .true.
            
      end if check_correct_line1
            
      call report_test("[test_rewind_file_seq]", &
                       has_failed, &
                       has_warned, &
                       "failed as incorrect first line_string read")   
                                                         
   end subroutine test_rewind_file_seq
   
   ! --------------------------------------------------------------------------

   subroutine test_go_to_line_seq(record_len)
      
      !!< Test the procedure that goes to a line number in a sequential file and changes the line number
      !!< via also reading a line and checking the string. This will use file_unit_1.

      integer, intent(in) :: record_len
            
      ! local variables
      integer :: line_number
      logical :: end_of_file
      character(len=record_len) :: line_string
                  
      line_number = 0
      rewind(unit=file_unit_1)
            
      ! first test
      call go_to_line_seq(file_unit_1, &
                          line_number, &
                          go_to_line_number = 9, &
                          exit_if_eor       = .false., &
                          exit_if_eof       = .false.)

      ! read next line
      call read_next_line_seq(file_unit_1, &
                              line_string, &
                              line_number, &
                              end_of_file, &
                              exit_if_eor = .false., &
                              exit_if_eof = .false.)
                    
      ! test that we have the correct line            
      check_correct_line1: if (trim(line_string) == '  4.883432E-06  3.137488E-01') then
               
         has_failed = .false.
            
      else check_correct_line1
            
         has_failed = .true.
            
      end if check_correct_line1
            
      call report_test("[test_go_to_line_seq]", &
                       has_failed, &
                       has_warned, &
                       "failed as incorrect first line_string read")   
            
      ! second test
      call go_to_line_seq(file_unit_1, &
                          line_number, &
                          go_to_line_number = 191, &
                          exit_if_eor       = .false., &
                          exit_if_eof       = .false.)

      ! read next line
      call read_next_line_seq(file_unit_1, &
                              line_string, &
                              line_number, &
                              end_of_file, &
                              exit_if_eor = .false., &
                              exit_if_eof = .false.)
                    
      ! test that we have the correct line            
      check_correct_line2: if (trim(line_string) == 'DELAYED_GROUPS') then
               
         has_failed = .false.
            
      else check_correct_line2
            
         has_failed = .true.
            
      end if check_correct_line2
            
      call report_test("[test_go_to_line_seq]", &
                       has_failed, &
                       has_warned, &
                       "failed as incorrect second line_string read")   

      ! third test backwards
      call go_to_line_seq(file_unit_1, &
                          line_number, &
                          go_to_line_number = 54, &
                          exit_if_eor       = .false., &
                          exit_if_eof       = .false.)

      ! read next line
      call read_next_line_seq(file_unit_1, &
                              line_string, &
                              line_number, &
                              end_of_file, &
                              exit_if_eor = .false., &
                              exit_if_eof = .false.)
                  
      ! test that we have the correct line            
      check_correct_line3: if (trim(line_string) == '  1.907523E-02 -7.976619E-04') then
               
         has_failed = .false.
            
      else check_correct_line3
            
         has_failed = .true.
            
      end if check_correct_line3
           
      call report_test("[test_go_to_line_seq]", &
                       has_failed, &
                       has_warned, &
                       "failed as incorrect third line_string read")                                
            
   end subroutine test_go_to_line_seq
   
   ! --------------------------------------------------------------------------

   !subroutine test_desired_keyword_not_first_found()
            
      ! cannot think how to test this procedure as it just ewrites and flexits
            
   !end subroutine test_desired_keyword_not_first_found
   
   ! --------------------------------------------------------------------------

   subroutine open_format_radmats_input(file_unit, &
                                        file_name)

      !!< Open the input files that are used for the unittesting
      
      integer, intent(inout) :: file_unit
      character(len=*), intent(in) :: file_name
      
      ! local variables
      logical :: file_exists
      integer :: iostat
      
      ! get a free io unit
      file_unit = free_unit()
            
      inquire(file=trim(file_name),exist=file_exists)
      
      no_input: if (.not. file_exists) then
      
         has_failed = .true.
         
         call report_test("[open_format_radmats_input]", &
                          has_failed, &
                          has_warned, &
                          "failed as input file not found, for input file "//trim(file_name))  
         
      else no_input
                  
         open(unit   = file_unit, &
              file   = trim(file_name), &
              status = 'old', &
              form   = 'formatted', &
              action = 'read', &
              access = 'sequential', &
              iostat = iostat)        
        
         check_open_fine: if (iostat /= 0) then
            
            has_failed = .true.
         
            call report_test("[open_format_radmats_input]", &
                             has_failed, &
                             has_warned, &
                             "failed as input file open iostat /= 0, for input file "//trim(file_name))   
                        
         end if check_open_fine
            
      end if no_input
   
   end subroutine open_format_radmats_input
   
   ! --------------------------------------------------------------------------

end subroutine test_radiation_materials_read_format_radmats_base
