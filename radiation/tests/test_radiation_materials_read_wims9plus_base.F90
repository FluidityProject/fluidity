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

subroutine test_radiation_materials_read_wims9plus_base
   
   !!< Test the procedures contained within the module radiation_materials_read_wims9plus_base in Radiation_Materials_Read_Wims9Plus_Base.F90
   
   use futils
   use unittest_tools  

   use radiation_materials_read_wims9plus_base 
   use radiation_materials_data_types
   use radiation_materials_create
   use radiation_materials_destroy
   use create_unittest_input 
         
   implicit none
   
   ! local variables
   logical :: has_failed,has_warned
   character(len=*), parameter :: rad_input_test_dir = 'data/'
   character(len=21), dimension(:), allocatable :: keyword_list_two_moments     
   type(np_radmat_type) :: np_radmat_1
   integer :: file_unit_1,file_unit_2,file_unit_3,file_unit_4
      
   ! none of these tests use warnings
   has_warned = .false.
   
   ! open the input files that will be used
   call open_wims9plus_input(file_unit = file_unit_1, &
                             file_name = rad_input_test_dir//'radiation_materials_read_input1.wims9plus_data')

   call open_wims9plus_input(file_unit = file_unit_2, &
                             file_name = rad_input_test_dir//'radiation_materials_read_input2.wims9plus_data')

   call open_wims9plus_input(file_unit = file_unit_3, &
                             file_name = rad_input_test_dir//'radiation_materials_read_input3.wims9plus_data')

   call open_wims9plus_input(file_unit = file_unit_4, &
                             file_name = rad_input_test_dir//'radiation_materials_read_input4.wims9plus_data')
   
   call test_keyword_list_initialise_wims9plus_two_scatter_moments(keyword_list_two_moments) 
   
   call test_deduce_number_of_items_wims9plus(file_unit                 = file_unit_1, &
                                              total_number_of_radmats   = 10,&
                                              number_of_scatter_moments = 2, &
                                              number_of_energy_groups   = 2, &
                                              number_non_thermal_groups = 1, &
                                              number_of_delayed_groups  = 6, &
                                              record_len                = 132, &
                                              keyword_list              = keyword_list_two_moments)
   
   call test_deduce_number_of_items_wims9plus(file_unit                 = file_unit_2, &
                                              total_number_of_radmats   = 10 ,&
                                              number_of_scatter_moments = 2, &
                                              number_of_energy_groups   = 6, &
                                              number_non_thermal_groups = 3, &
                                              number_of_delayed_groups  = 8, &
                                              record_len                = 132, &
                                              keyword_list              = keyword_list_two_moments)
   
   call test_deduce_number_of_items_wims9plus(file_unit                 = file_unit_3, &
                                              total_number_of_radmats   = 10 ,&
                                              number_of_scatter_moments = 2, &
                                              number_of_energy_groups   = 11, &
                                              number_non_thermal_groups = 5, &
                                              number_of_delayed_groups  = 0, &
                                              record_len                = 132, &
                                              keyword_list              = keyword_list_two_moments)
   
   call test_deduce_number_of_items_wims9plus(file_unit                 = file_unit_4, &
                                              total_number_of_radmats   = 10 ,&
                                              number_of_scatter_moments = 2, &
                                              number_of_energy_groups   = 172, &
                                              number_non_thermal_groups = 92, &
                                              number_of_delayed_groups  = 0, &
                                              record_len                = 132, &
                                              keyword_list              = keyword_list_two_moments)

   ! create the np_radmat's associated with the input files opened above
   call create_np_radmat_input1(np_radmat_1)
   
   call test_read_scatter_mom_wims9plus(file_unit      = file_unit_1, &
                                        record_len     = 132, &
                                        dataset_radmat = np_radmat_1%dataset_radmats(1), &
                                        np_radmat_size = np_radmat_1%np_radmat_size, &
                                        keyword_list   = keyword_list_two_moments) 
                                       
   call test_read_absorption_and_total_wims9plus(file_unit      = file_unit_1, &
                                                 record_len     = 132, &
                                                 dataset_radmat = np_radmat_1%dataset_radmats(1), &
                                                 np_radmat_size = np_radmat_1%np_radmat_size, &
                                                 keyword_list   = keyword_list_two_moments) 
   
   call test_calculate_removal_moments_wims9plus(file_unit      = file_unit_1, &
                                                 record_len     = 132, &
                                                 dataset_radmat = np_radmat_1%dataset_radmats(1), &
                                                 np_radmat_size = np_radmat_1%np_radmat_size, &
                                                 keyword_list   = keyword_list_two_moments) 
      
   call test_read_or_associate_transport_and_diffusion_wims9plus(file_unit      = file_unit_1, &
                                                                 record_len     = 132, &
                                                                 dataset_radmat = np_radmat_1%dataset_radmats(1), &
                                                                 np_radmat_size = np_radmat_1%np_radmat_size, &
                                                                 keyword_list   = keyword_list_two_moments, &
                                                                 problem_dim    = 3)

   call test_read_fission_data_wims9plus(file_unit                 = file_unit_1, &
                                         record_len                = 132, &
                                         dataset_radmat            = np_radmat_1%dataset_radmats(1), &
                                         np_radmat_size            = np_radmat_1%np_radmat_size, &
                                         keyword_list              = keyword_list_two_moments, &
                                         number_non_thermal_groups = 1 ) 
      
   call test_read_velocity_style_data_wims9plus(file_unit           = file_unit_1, &
                                                keyword_line_number = 1055, &
                                                keyword_find        = (/14/), &
                                                record_len          = 132, &
                                                dataset_radmat      = np_radmat_1%dataset_radmats(1), &
                                                np_radmat_size      = np_radmat_1%np_radmat_size, &
                                                keyword_list        = keyword_list_two_moments)
   
   call test_read_velocity_style_data_wims9plus(file_unit           = file_unit_1, &
                                                keyword_line_number = 1044, &
                                                keyword_find        = (/16/), &
                                                record_len          = 132, &
                                                dataset_radmat      = np_radmat_1%dataset_radmats(1), &
                                                np_radmat_size      = np_radmat_1%np_radmat_size, &
                                                keyword_list        = keyword_list_two_moments)
   
   call test_read_velocity_style_data_wims9plus(file_unit           = file_unit_1, &
                                                keyword_line_number = 1060, &
                                                keyword_find        = (/17/), &
                                                record_len          = 132, &
                                                dataset_radmat      = np_radmat_1%dataset_radmats(1), &
                                                np_radmat_size      = np_radmat_1%np_radmat_size, &
                                                keyword_list        = keyword_list_two_moments)
   
   ! read the delayed data - beta, lambda and spectrum
   call test_read_delayed_data_wims9plus(file_unit               = file_unit_1, &
                                         record_len              = 132, &
                                         dataset_radmat          = np_radmat_1%dataset_radmats(1), &
                                         np_radmat_size          = np_radmat_1%np_radmat_size, &
                                         delayed_lambda_spectrum = np_radmat_1%delayed_lambda_spectrum, &
                                         keyword_list            = keyword_list_two_moments) 
   
   ! read the energy released through the power rotuine
   call test_read_power_data_wims9plus(file_unit      = file_unit_1, &
                                       record_len     = 132, &
                                       dataset_radmat = np_radmat_1%dataset_radmats(1), &
                                       np_radmat_size = np_radmat_1%np_radmat_size, &
                                       test_power     = .false., &
                                       keyword_list   = keyword_list_two_moments)

   call test_calculate_power_wims9plus(dataset_radmat = np_radmat_1%dataset_radmats(1), &
                                       np_radmat_size = np_radmat_1%np_radmat_size) 
   
   ! get rid of the comparison np_radmat's 
   call destroy(np_radmat_1)   

   call test_read_wims9plus_siga_style_xsection(file_unit           = file_unit_1, &
                                                keyword_line_number = 976, &
                                                number_of_groups    = 2, &
                                                record_len          = 132, &
                                                values_expected     = (/8.110421E-06,&
                                                                      &1.643668E-04/), &
                                                number_groups_to_read_in = 2)

   call test_read_wims9plus_siga_style_xsection(file_unit           = file_unit_2, &
                                                keyword_line_number = 1969, &
                                                number_of_groups    = 6, &
                                                record_len          = 132, &
                                                values_expected     = (/1.319194E-01,&
                                                                      &2.411816E-01, &
                                                                      &3.394288E-01, &
                                                                      &3.457756E-01, &
                                                                      &0.0, &
                                                                      &0.0/), &
                                                number_groups_to_read_in = 4) ! this is intentially 4, less than number of groups

   call test_read_wims9plus_siga_style_xsection(file_unit           = file_unit_3, &
                                                keyword_line_number = 3001, &
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

   call test_read_wims9plus_siga_style_xsection(file_unit           = file_unit_4, &
                                                keyword_line_number = 82, &
                                                number_of_groups    = 172, &
                                                record_len          = 132, &
                                                values_expected     = (/ &
 & 2.089550E-01,  2.485903E-01,  2.644520E-01,  2.984862E-01,  3.185046E-01,  3.539481E-01,  3.858338E-01,  4.903702E-01, &
 & 5.163930E-01,  5.248332E-01,  4.775521E-01,  4.099437E-01,  4.847418E-01,  4.518094E-01,  5.865291E-01,  7.013625E-01, &
 & 4.122334E-01,  2.872880E-01,  2.669537E-01,  2.683544E-01,  2.947675E-01,  3.370069E-01,  4.822639E-01,  7.277249E-01, &
 & 8.941768E-01,  9.970101E-01,  1.003968E+00,  1.098100E+00,  1.285010E+00,  1.397355E+00,  1.613572E+00,  1.941810E+00, &
 & 2.153307E+00,  2.320652E+00,  2.554434E+00,  2.769150E+00,  3.076663E+00,  3.359049E+00,  3.620728E+00,  3.874037E+00, &
 & 4.026507E+00,  4.518819E+00,  5.114379E+00,  5.616333E+00,  6.321034E+00,  6.956129E+00,  7.844631E+00,  8.632944E+00, &
 & 9.619372E+00,  1.058024E+01,  1.178970E+01,  1.327812E+01,  1.461586E+01,  1.587692E+01,  1.665795E+01,  1.808908E+01, &
 & 1.940919E+01,  2.083848E+01,  2.235409E+01,  2.511900E+01,  2.886036E+01,  3.164473E+01,  3.626163E+01,  4.261498E+01, &
 & 4.654334E+01,  5.173311E+01,  5.880960E+01,  6.267722E+01,  6.682759E+01,  7.075751E+01,  7.284831E+01,  7.476070E+01, &
 & 7.758511E+01,  8.074845E+01,  8.351167E+01,  8.692976E+01,  9.038523E+01,  9.393508E+01,  9.761094E+01,  1.020737E+02, &
 & 1.088284E+02,  1.158879E+02,  1.232408E+02,  1.302733E+02,  1.347374E+02,  1.386949E+02,  1.432659E+02,  1.495224E+02, &
 & 1.572276E+02,  1.617291E+02,  1.677197E+02,  1.731993E+02,  1.775074E+02,  1.822171E+02,  1.869414E+02,  1.913817E+02, &
 & 1.928876E+02,  1.944672E+02,  1.967196E+02,  2.008826E+02,  2.035590E+02,  2.047278E+02,  2.066640E+02,  2.087655E+02, &
 & 2.108638E+02,  2.130535E+02,  2.151775E+02,  2.175135E+02,  2.190806E+02,  2.199592E+02,  2.214925E+02,  2.230495E+02, &
 & 2.241339E+02,  2.257568E+02,  2.278691E+02,  2.292958E+02,  2.300670E+02,  2.307459E+02,  2.312166E+02,  2.319174E+02, &
 & 2.328348E+02,  2.335032E+02,  2.339801E+02,  2.347147E+02,  2.353331E+02,  2.357657E+02,  2.364720E+02,  2.373150E+02, &
 & 2.380760E+02,  2.394885E+02,  2.407173E+02,  2.421575E+02,  2.436422E+02,  2.454617E+02,  2.489747E+02,  2.529494E+02, &
 & 2.560681E+02,  2.574273E+02,  2.591718E+02,  2.613451E+02,  2.623516E+02,  2.636213E+02,  2.652172E+02,  2.659079E+02, &
 & 2.662866E+02,  2.668592E+02,  2.676411E+02,  2.680734E+02,  2.679538E+02,  2.673873E+02,  2.669008E+02,  2.659517E+02, &
 & 2.652504E+02,  2.646380E+02,  2.634840E+02,  2.632521E+02,  2.632120E+02,  2.622048E+02,  2.613729E+02,  2.598454E+02, &
 & 2.610915E+02,  2.593094E+02,  2.618567E+02,  2.575697E+02,  2.589347E+02,  2.511682E+02,  2.502860E+02,  2.678501E+02, &
 & 2.569172E+02,  2.541621E+02,  2.636447E+02,  2.454069E+02/), &
                                                number_groups_to_read_in = 172)

   call test_read_wims9plus_sigs_style_xsection(file_unit           = file_unit_2, &
                                                keyword_line_number = 1971, &
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
                                                 keyword_list       = keyword_list_two_moments)

   call test_read_wims9plus_sigs_style_xsection(file_unit           = file_unit_3, &
                                                keyword_line_number = 2910, &
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
                                                keyword_list       = keyword_list_two_moments)

   call test_count_number_lines_with_keyword(file_unit              = file_unit_1, &
                                             record_len             = 132, &
                                             keyword_list           = keyword_list_two_moments, &
                                             keyword_to_count       = 2, &
                                             keyword_to_stop        = 1, &
                                             keyword_count_expected = 49, &
                                             line_number_expected   = 867) 

   call test_find_line_with_any_keyword(file_unit                    = file_unit_1, &
                                        record_len                   = 132, &
                                        keyword_list                 = keyword_list_two_moments, &
                                        first_keyword_found_expected = 2, &
                                        line_number_expected         = 61) 
   
   ! this will use file_unit_1
   call test_find_line_with_any_desired_keyword(record_len   = 132, &
                                                 keyword_list = keyword_list_two_moments) 
   
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

   call test_keyword_list_cleanup(keyword_list_two_moments) 
   
   close(file_unit_1)
   close(file_unit_2)
   close(file_unit_3)
   close(file_unit_4)

   !call test_desired_keyword_not_first_found() ! cannot be tested

contains

   ! --------------------------------------------------------------------------

   subroutine test_keyword_list_initialise_wims9plus_two_scatter_moments(keyword_list)
      
      !!< Test the procedure that initialises the keyword list assuming two scatter moments
      
      character(len=21), dimension(:), allocatable, intent(inout) :: keyword_list     
      
      ! local variables
      integer :: k
      character(len=21), dimension(23) :: keyword_list_check     
                   
      call keyword_list_initialise_wims9plus(number_of_scatter_moments = 2, &
                                             keyword_list = keyword_list)
      
      ! check that the size of the keyword_list is correct
      first_check: if (size(keyword_list) /= 23) then
         
        has_failed = .true.
         
      else first_check

        has_failed = .false.
      
      end if first_check

      call report_test("[test_keyword_list_initialise_wims9plus_two_scatter_moments]", &
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
      keyword_list_check(9)  = 'DIFFUSIONX' 
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
      keyword_list_check(22) = 'SCATTER MOMENT      1'
      keyword_list_check(23) = 'SCATTER MOMENT      2'     
     
      keyword_loop: do k = 1,size(keyword_list_check)
      
         second_check: if (keyword_list(k) /= keyword_list_check(k)) then
         
           has_failed = .true.
         
         else second_check

           has_failed = .false.
      
         end if second_check

         call report_test("[test_keyword_list_initialise_wims9plus_two_scatter_moments]", &
                         has_failed, &
                         has_warned, &
                         "failed keyword_list has incorrect element "//int2str(k))   
      
      end do keyword_loop
            
   end subroutine test_keyword_list_initialise_wims9plus_two_scatter_moments
   
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

   subroutine test_deduce_number_of_items_wims9plus(file_unit, &
                                                    total_number_of_radmats, &
                                                    number_of_scatter_moments, &
                                                    number_of_energy_groups, &
                                                    number_non_thermal_groups, &
                                                    number_of_delayed_groups, &
                                                    record_len, &
                                                    keyword_list)
      
      !!< Test the procedure that deduces from a wims9plus input file the number of rad mats, the number of scatter moments, 
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
      integer :: total_number_of_radmats_wims9plus
      integer :: number_of_energy_groups_wims9plus
      integer :: number_non_thermal_groups_wims9plus
      integer :: number_of_scatter_moments_wims9plus
      integer :: number_of_delayed_groups_wims9plus
                                                
      call deduce_number_of_items_wims9plus(file_unit, &
                                            total_number_of_radmats_wims9plus, &
                                            number_of_energy_groups_wims9plus, &
                                            number_non_thermal_groups_wims9plus, &
                                            number_of_scatter_moments_wims9plus, &
                                            number_of_delayed_groups_wims9plus, &
                                            record_len, &
                                            keyword_list)
     
      ! first check
      check_number_of_energy_groups: if (number_of_energy_groups_wims9plus == number_of_energy_groups) then
            
         has_failed = .false.
                                            
      else check_number_of_energy_groups

         has_failed = .true.
                            
      end if check_number_of_energy_groups
            
      call report_test("[test_deduce_number_of_items_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong number of energy groups")   
            
      ! second check
      check_number_of_energy_thermal_groups: if (number_non_thermal_groups_wims9plus == number_non_thermal_groups) then
            
         has_failed = .false.
                                            
      else check_number_of_energy_thermal_groups

         has_failed = .true.
                            
      end if check_number_of_energy_thermal_groups
            
      call report_test("[test_deduce_number_of_items_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong number of thermal energy groups")   
            
      ! third check
      check_number_of_delayed_groups: if (number_of_delayed_groups_wims9plus == number_of_delayed_groups) then
            
         has_failed = .false.
                                            
      else check_number_of_delayed_groups

         has_failed = .true.
                            
      end if check_number_of_delayed_groups
            
      call report_test("[test_deduce_number_of_items_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong number of delayed groups")   

      ! fourth check
      check_number_of_scatter_moments: if (number_of_scatter_moments_wims9plus == number_of_scatter_moments) then
            
         has_failed = .false.
                                            
      else check_number_of_scatter_moments

         has_failed = .true.
                            
      end if check_number_of_scatter_moments
            
      call report_test("[test_deduce_number_of_items_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong number of scatter moments")   

      ! fourth check
      check_number_of_macro: if (total_number_of_radmats_wims9plus == total_number_of_radmats) then
            
         has_failed = .false.
                                            
      else check_number_of_macro

         has_failed = .true.
                            
      end if check_number_of_macro
            
      call report_test("[test_deduce_number_of_items_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "failed as wrong number of radiation materials")   
                                          
   end subroutine test_deduce_number_of_items_wims9plus
   
   ! --------------------------------------------------------------------------

   subroutine test_read_scatter_mom_wims9plus(file_unit, &
                                             record_len, &
                                             dataset_radmat, &
                                             np_radmat_size, &
                                             keyword_list)
      
      !!< Test the procedure that reads in the scatter data numbers
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(np_radmat_size_type), intent(in) :: np_radmat_size      
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
      type(dataset_radmat_type) :: dataset_radmat_wims9plus
                              
      ! intiialise the read in data type
      call allocate(dataset_radmat_wims9plus, &
                    np_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_wims9plus)
      
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
               
               call report_test("[test_read_scatter_mom_wims9plus]", &
                                has_failed, &
                                has_warned, &
                                "failed to find keyword MACRO as needed")   
                        
            end if not_found
            
            deallocate(keyword_find)
                                    
            ! keep the line number associated with this MACRO line
            line_number_macro = line_number
                          
            ! now read in the required values
            call read_scatter_mom_wims9plus(file_unit, &
                                            line_number, &
                                            dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat), &
                                            record_len, &
                                            keyword_list)
            
            ! check each moment matrix
            mom_loop: do m = 1,size(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%scatter,3)

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%scatter(:,:,m), &
                                      dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%scatter(:,:,m))

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
            
      call report_test("[test_read_scatter_mom_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_wims9plus)
      
   end subroutine test_read_scatter_mom_wims9plus
   
   ! --------------------------------------------------------------------------

   subroutine test_read_absorption_and_total_wims9plus(file_unit, &
                                                       record_len, &
                                                       dataset_radmat, &
                                                       np_radmat_size, &
                                                       keyword_list)
      
      !!< Test the procedure that reads in the absorption or total data numbers. This will also implicitly read in the scatter data.
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(np_radmat_size_type), intent(in) :: np_radmat_size      
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
      type(dataset_radmat_type) :: dataset_radmat_wims9plus
                              
      ! intiialise the read in data type
      call allocate(dataset_radmat_wims9plus, &
                    np_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_wims9plus)
      
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
               
               call report_test("[test_read_absorption_and_total_wims9plus]", &
                                has_failed, &
                                has_warned, &
                                "failed to find keyword MACRO as needed")   
                        
            end if not_found
            
            deallocate(keyword_find)
                                    
            ! keep the line number associated with this MACRO line
            line_number_macro = line_number
                          
            ! now read in the required values
            call read_absorption_and_total_wims9plus(file_unit, &
                                                     line_number, &
                                                     dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat), &
                                                     record_len, &
                                                     keyword_list)
            
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%total, &
                                   dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%total)

            if (has_failed) exit pmat_loop
            
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%absorption, &
                                   dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%absorption)

            if (has_failed) exit pmat_loop
            
            ! go back to the line number of the current macro
            call go_to_line_seq(file_unit, &
                                line_number, &
                                go_to_line_number = line_number_macro, &
                                exit_if_eor       = .true.,&
                                exit_if_eof       = .true.) 
                                    
         end do rmat_loop
            
      end do pmat_loop
            
      call report_test("[test_read_absorption_and_total_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_wims9plus)
            
   end subroutine test_read_absorption_and_total_wims9plus
   
   ! --------------------------------------------------------------------------

   subroutine test_calculate_removal_moments_wims9plus(file_unit, &
                                                       record_len, &
                                                       dataset_radmat, &
                                                       np_radmat_size, &
                                                       keyword_list)
      
      !!< Test the procedure that calculates the removal moments. This will also read in the total and absorption (which itself will 
      !!< read in the scatter data)
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(np_radmat_size_type), intent(in) :: np_radmat_size      
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
      type(dataset_radmat_type) :: dataset_radmat_wims9plus
                              
      ! intiialise the read in data type
      call allocate(dataset_radmat_wims9plus, &
                    np_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_wims9plus)
      
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
               
               call report_test("[test_calculate_removal_moments_wims9plus]", &
                                has_failed, &
                                has_warned, &
                                "failed to find keyword MACRO as needed")   
                        
            end if not_found
            
            deallocate(keyword_find)
                                    
            ! keep the line number associated with this MACRO line
            line_number_macro = line_number
                          
            ! now read in the required values
            call calculate_removal_moments_wims9plus(file_unit, &
                                                     line_number, &
                                                     dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat), &
                                                     record_len, &
                                                     keyword_list)
            
            ! check each moment matrix
            mom_loop: do m = 1,size(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%removal,2)

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%removal(:,m), &
                                      dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%removal(:,m))

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
            
      call report_test("[test_calculate_removal_moments_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_wims9plus)
                  
   end subroutine test_calculate_removal_moments_wims9plus
   
   ! --------------------------------------------------------------------------

   subroutine test_read_or_associate_transport_and_diffusion_wims9plus(file_unit, &
                                                                       record_len, &
                                                                       dataset_radmat, &
                                                                       np_radmat_size, &
                                                                       keyword_list, &
                                                                       problem_dim)
      
      !!< Test the procedure that reads the transport or diffusion data. This will also read in the total and absorption (which itself will 
      !!< read in the scatter data)
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(np_radmat_size_type), intent(in) :: np_radmat_size      
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
      type(dataset_radmat_type) :: dataset_radmat_wims9plus
                              
      ! intiialise the read in data type
      call allocate(dataset_radmat_wims9plus, &
                    np_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_wims9plus)
      
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
               
               call report_test("[test_read_or_associate_transport_and_diffusion_wims9plus]", &
                                has_failed, &
                                has_warned, &
                                "failed to find keyword MACRO as needed")   
                        
            end if not_found
            
            deallocate(keyword_find)
                                    
            ! keep the line number associated with this MACRO line
            line_number_macro = line_number
            
            ! now read in the required values for each dimension
            dim_loop: do idim = 1,problem_dim              

               call read_or_associate_transport_and_diffusion_wims9plus(file_unit, &
                                                                        line_number, &
                                                                        dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat), &
                                                                        idim, &
                                                                        record_len, &
                                                                        keyword_list)
            
               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%transport(:,idim), &
                                      dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%transport(:,idim))

               if (has_failed) exit pmat_loop

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%diffusion(:,idim), &
                                      dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%diffusion(:,idim))

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
            
      call report_test("[test_read_or_associate_transport_and_diffusion_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_wims9plus)
      
   end subroutine test_read_or_associate_transport_and_diffusion_wims9plus
   
   ! --------------------------------------------------------------------------

   subroutine test_read_fission_data_wims9plus(file_unit, &
                                               record_len, &
                                               dataset_radmat, &
                                               np_radmat_size, &
                                               keyword_list, &
                                               number_non_thermal_groups)
      
      !!< Test the procedure that reads the fission data. 
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(np_radmat_size_type), intent(in) :: np_radmat_size      
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
      type(dataset_radmat_type) :: dataset_radmat_wims9plus
                              
      ! intiialise the read in data type
      call allocate(dataset_radmat_wims9plus, &
                    np_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_wims9plus)
      
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
               
               call report_test("[test_read_fission_data_wims9plus]", &
                                has_failed, &
                                has_warned, &
                                "failed to find keyword MACRO as needed")   
                        
            end if not_found
            
            deallocate(keyword_find)
                                    
            ! keep the line number associated with this MACRO line
            line_number_macro = line_number
                          
            ! now read in the required values
            call read_fission_data_wims9plus(file_unit, &
                                             line_number, &
                                             dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat), &
                                             number_non_thermal_groups, &
                                             record_len, &
                                             keyword_list)
            
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%fission, &
                                   dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%fission)

            if (has_failed) exit pmat_loop
            
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%production, &
                                   dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%production)

            if (has_failed) exit pmat_loop
            
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%np_released_per_fission, &
                                   dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%np_released_per_fission)

            if (has_failed) exit pmat_loop
            
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%prompt_spectrum, &
                                   dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%prompt_spectrum)

            if (has_failed) exit pmat_loop
                                    
            ! go back to the line number of the current macro
            call go_to_line_seq(file_unit, &
                                line_number, &
                                go_to_line_number = line_number_macro, &
                                exit_if_eor       = .true.,&
                                exit_if_eof       = .true.) 
                                    
         end do rmat_loop
            
      end do pmat_loop
            
      call report_test("[test_read_fission_data_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_wims9plus)
      
   end subroutine test_read_fission_data_wims9plus
   
   ! --------------------------------------------------------------------------

   subroutine test_read_velocity_style_data_wims9plus(file_unit, &
                                                      keyword_line_number, &
                                                      keyword_find, &
                                                      record_len, &
                                                      dataset_radmat, &
                                                      np_radmat_size, &
                                                      keyword_list)
      
      !!< Test the procedure that reads in a velocity style data numbers
      !!< This assumes that the keyword is at the keyword_line_number which is read to first
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: keyword_line_number 
      integer, dimension(1), intent(in) :: keyword_find
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(np_radmat_size_type), intent(in) :: np_radmat_size
      character(len=21), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      integer :: line_number
      integer :: rmat
      integer :: pmat
      integer :: number_of_groups
      type(dataset_radmat_type) :: dataset_radmat_wims9plus
                  
      line_number = 0
      rewind(unit=file_unit)
            
      ! go to the desired keyword line
      call go_to_line_seq(file_unit, &
                          line_number, &
                          go_to_line_number=keyword_line_number, &
                          exit_if_eor=.false., &
                          exit_if_eof=.false.)
          
      ! intiialise the read in data type
      call allocate(dataset_radmat_wims9plus, &
                    np_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_wims9plus)
            
      number_of_groups_if: if (keyword_find(1) == 17) then
            
         number_of_groups = np_radmat_size%number_of_delayed_groups
            
      else number_of_groups_if
               
         number_of_groups = np_radmat_size%number_of_energy_groups
            
      end if number_of_groups_if
            
      ! now read in the required values
      call read_velocity_style_data_wims9plus(file_unit, &
                                              dataset_radmat_wims9plus, &
                                              number_of_groups, &
                                              keyword_find, &
                                              record_len, &
                                              keyword_list)

      ! now check the read in values with the expected values to a default tolerance
      pmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)
              
            check_values: if (keyword_find(1) == 14) then
            
               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%velocity, &
                                      dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%velocity)
                                                              
            else if (keyword_find(1) == 15) then 

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%power, &
                                      dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%power)

            else if (keyword_find(1) == 16) then 

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission, &
                                      dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission)
                                        
            else if (keyword_find(1) == 17) then 

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%beta, &
                                      dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%beta)

            end if check_values

            if (has_failed) exit pmat_loop
            
         end do rmat_loop
            
      end do pmat_loop
            
      call report_test("[test_read_velocity_style_data_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_wims9plus)
            
   end subroutine test_read_velocity_style_data_wims9plus
   
   ! --------------------------------------------------------------------------

   subroutine test_read_delayed_data_wims9plus(file_unit, &
                                               record_len, &
                                               dataset_radmat, &
                                               np_radmat_size, &
                                               delayed_lambda_spectrum, &
                                               keyword_list)
      
      !!< Test the procedure that reads in the delayed data numbers
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(np_radmat_size_type), intent(in) :: np_radmat_size      
      type(delayed_lambda_spectrum_type), intent(in) :: delayed_lambda_spectrum      
      character(len=21), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      integer :: pmat
      integer :: rmat
      type(dataset_radmat_type) :: dataset_radmat_wims9plus
      type(delayed_lambda_spectrum_type) :: delayed_lambda_spectrum_wims9plus
      
      ! intiialise the read in data types
      call allocate(dataset_radmat_wims9plus, &
                    np_radmat_size, &
                    dmat=1)

      call zero(dataset_radmat_wims9plus)
 
      call allocate(delayed_lambda_spectrum_wims9plus, &
                    np_radmat_size%number_of_delayed_groups, &
                    np_radmat_size%number_of_energy_groups)

      call zero(delayed_lambda_spectrum_wims9plus)
                        
      ! now read in the required values
      call read_delayed_data_wims9plus(file_unit, &
                                       dataset_radmat_wims9plus, &
                                       record_len, &
                                       keyword_list, &
                                       delayed_lambda_spectrum_wims9plus, &
                                       read_delayed_lambda_spectrum = .true.)

      ! now check the read in values with the expected values to a default tolerance
            
      pmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)
                          
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%beta, &
                                   dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%beta)

            if (has_failed) exit pmat_loop
                                    
         end do rmat_loop
            
      end do pmat_loop
            
      call report_test("[test_read_delayed_data_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "beta failed values read in not what expected,")   
            
      has_failed =  fnequals(delayed_lambda_spectrum%lambda, &
                             delayed_lambda_spectrum_wims9plus%lambda)

      call report_test("[test_read_delayed_data_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "lambda failed values read in not what expected")   
            
      has_failed =  fnequals(delayed_lambda_spectrum%spectrum, &
                             delayed_lambda_spectrum_wims9plus%spectrum)
            
      call report_test("[test_read_delayed_data_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "spectrum failed values read in not what expected")   
                        
      call destroy(dataset_radmat_wims9plus)
      call destroy(delayed_lambda_spectrum_wims9plus)
                  
   end subroutine test_read_delayed_data_wims9plus
   
   ! --------------------------------------------------------------------------

   subroutine test_read_power_data_wims9plus(file_unit, &
                                             record_len, &
                                             dataset_radmat, &
                                             np_radmat_size, &
                                             test_power, &
                                             keyword_list)
      
      !!< Test the procedure that reads in the power data numbers
      !!< The power variable may not be read in but formed from the fission and energy released
      !!< but as this routine will not read the fission we may not always want to test power (only energy released)
      
      integer, intent(in) :: file_unit
      integer, intent(in) :: record_len  
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      type(np_radmat_size_type), intent(in) :: np_radmat_size      
      logical, intent(in) :: test_power
      character(len=21), dimension(:), intent(in) :: keyword_list
      
      ! local variables
      integer :: pmat
      integer :: rmat
      type(dataset_radmat_type) :: dataset_radmat_wims9plus
                              
      ! intiialise the read in data type
      call allocate(dataset_radmat_wims9plus,np_radmat_size,dmat=1)
      call zero(dataset_radmat_wims9plus)
                        
      ! now read in the required values
      call read_power_data_wims9plus(file_unit, &
                                      dataset_radmat_wims9plus, &
                                      record_len, &
                                      keyword_list)

      ! now check the read in values with the expected values to a default tolerance
      pmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)
                          
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission, &
                                   dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%energy_released_per_fission)

            if (has_failed) exit pmat_loop
               
            check_power: if (test_power) then                                            

               has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%power, &
                                      dataset_radmat_wims9plus%physical_radmats(pmat)%radmats(rmat)%power)

               if (has_failed) exit pmat_loop
                  
            end if check_power
                  
         end do rmat_loop
            
      end do pmat_loop
            
      call report_test("[test_read_power_data_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      call destroy(dataset_radmat_wims9plus)
      
   end subroutine test_read_power_data_wims9plus
   
   ! --------------------------------------------------------------------------

   subroutine test_calculate_power_wims9plus(dataset_radmat, &
                                             np_radmat_size)
      
      !!< Test the procedure that calculates the power data from the fission and energy released
      
      type(dataset_radmat_type), intent(inout) :: dataset_radmat
      type(np_radmat_size_type), intent(in) :: np_radmat_size      
      
      ! local variables
      type(dataset_radmat_type) :: dataset_radmat_tmp
      integer :: pmat,rmat
      
      ! intiialise the tmp data type
      call allocate(dataset_radmat_tmp,np_radmat_size,dmat=1)
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
      call calculate_power_wims9plus(dataset_radmat)
      
      ! check the power
      pmat_loop_check: do pmat = 1,size(dataset_radmat%physical_radmats)
            
         rmat_loop_check: do rmat = 1,size(dataset_radmat%physical_radmats(pmat)%radmats)
               
            has_failed =  fnequals(dataset_radmat%physical_radmats(pmat)%radmats(rmat)%power, &
                                   dataset_radmat_tmp%physical_radmats(pmat)%radmats(rmat)%power)

            if (has_failed) exit pmat_loop_check
                         
         end do rmat_loop_check
            
      end do pmat_loop_check

      call report_test("[test_calculate_power_wims9plus]", &
                       has_failed, &
                       has_warned, &
                       "failed values calculated not what expected for pmat "//int2str(pmat)//" for rmat "//int2str(rmat))   
            
      call destroy(dataset_radmat_tmp)
      
   end subroutine test_calculate_power_wims9plus
   
   ! --------------------------------------------------------------------------

   subroutine test_read_wims9plus_siga_style_xsection(file_unit, &
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
      real, dimension(:), allocatable :: values_wims9plus
                  
      line_number = 0
      rewind(unit=file_unit)
            
      ! go to the desired keyword line
      call go_to_line_seq(file_unit, &
                          line_number, &
                          go_to_line_number=keyword_line_number, &
                          exit_if_eor=.false., &
                          exit_if_eof=.false.)
            
      ! intiialise the values read in array
      allocate(values_wims9plus(number_of_groups))
      values_wims9plus = 0.0
            
      ! now read in the required values
      call read_wims9plus_siga_style_xsection(values_wims9plus, &
                                              line_number, &
                                              file_unit, &
                                              record_len, &
                                              number_groups_to_read_in=number_groups_to_read_in)

      ! now check the read in values with the expected values to a default tolerance
      check_values: if (fequals(values_expected,values_wims9plus)) then
            
         has_failed = .false.
                                            
      else check_values

         has_failed = .true.
                            
      end if check_values
            
      call report_test("[test_read_wims9plus_siga_style_xsection]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      deallocate(values_wims9plus)                                  
      
   end subroutine test_read_wims9plus_siga_style_xsection
   
   ! --------------------------------------------------------------------------

   subroutine test_read_wims9plus_sigs_style_xsection(file_unit, &
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
      real, dimension(:,:), allocatable :: values_wims9plus
                        
      line_number = 0
      rewind(unit=file_unit)
            
      ! go to the desired keyword line
      call go_to_line_seq(file_unit, &
                          line_number, &
                          go_to_line_number=keyword_line_number, &
                          exit_if_eor=.false., &
                          exit_if_eof=.false.)
            
      ! intiialise the values read in array
      allocate(values_wims9plus(number_of_groups,number_of_groups))
      values_wims9plus = 0.0
            
      ! now read in the required values
      call read_wims9plus_sigs_style_xsection(values_wims9plus, &
                                              line_number, &
                                              file_unit, &
                                              record_len, &
                                              keyword_list)

      ! now check the read in values with the expected values to a default tolerance
      check_values: if (fequals(values_expected,values_wims9plus)) then
            
         has_failed = .false.
                                           
      else check_values

         has_failed = .true.
                            
      end if check_values
            
      call report_test("[test_read_wims9plus_sigs_style_xsection]", &
                       has_failed, &
                       has_warned, &
                       "failed values read in not what expected")   
                        
      deallocate(values_wims9plus)
            
   end subroutine test_read_wims9plus_sigs_style_xsection
   
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
            
      ! first search for word one - being MACRO in the wims9plus list            
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
            
      ! first search for word one - being MACRO in the wims9plus list
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
                        
      check_correct_line_number: if (line_number == 867) then
            
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
                        
      check_correct_line_number_again: if (line_number == 872) then
            
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
      line_loop1: do line = 1,122
               
         call read_next_line_seq(file_unit_1, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor = .false., &
                                 exit_if_eof = .false.)
            
      end do line_loop1
            
      ! test that we have the correct line            
      check_correct_line1: if (trim(line_string) == '  1.288531E-04  1.792035E-03') then
               
         has_failed = .false.
            
      else check_correct_line1
            
         has_failed = .true.
            
      end if check_correct_line1
            
      call report_test("[test_read_next_line_seq]", &
                       has_failed, &
                       has_warned, &
                       "failed as incorrect first line_string read")   
            
      ! read some more lines
      line_loop2: do line = 1,921
            
         call read_next_line_seq(file_unit_1, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor = .false., &
                                 exit_if_eof = .false.)
            
      end do line_loop2
            
      ! test that we have the correct line            
      check_correct_line2: if (trim(line_string) == 'END') then
               
         has_failed = .false.
            
      else check_correct_line2
            
         has_failed = .true.
            
      end if check_correct_line2
            
      call report_test("[test_read_next_line_seq]", &
                       has_failed, &
                       has_warned, &
                       "failed as incorrect second line_string read")   

      ! read some more lines
      line_loop3: do line = 1,29
            
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
      line_loop1: do line = 1,122
               
         call read_next_line_seq(file_unit_1, &
                                 line_string, &
                                 line_number, &
                                 end_of_file, &
                                 exit_if_eor = .false., &
                                 exit_if_eof = .false.)
            
      end do line_loop1
            
      ! read some lines backward
      line_loop2: do line = 1,13
               
         call read_previous_line_seq(file_unit_1, &
                                     line_string, &
                                     line_number, &
                                     exit_if_eor = .false., &
                                     exit_if_eof = .false.)
            
      end do line_loop2
          
      ! test that we have the correct line            
      check_correct_line1: if (trim(line_string) == '  2.011938E+00  5.258210E-02') then
               
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
      check_correct_line1: if (trim(line_string) == 'MESH_MATERIALS ') then
               
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
      check_correct_line1: if (trim(line_string) == 'EIGENVALUE ') then
               
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
                          go_to_line_number = 595, &
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
      check_correct_line2: if (trim(line_string) == 'GROUP      1 FIRST      1 LAST      2') then
               
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
                          go_to_line_number = 55, &
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
      check_correct_line3: if (trim(line_string) == '  B10          2   2.206270E-04  B11          4   8.924460E-04  C           13   8.148154E-02  He3         23   3.718941E-11') then
               
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

   subroutine open_wims9plus_input(file_unit, &
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
         
         call report_test("[open_wims9plus_input]", &
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
         
            call report_test("[open_wims9plus_input]", &
                             has_failed, &
                             has_warned, &
                             "failed as input file open iostat /= 0, for input file "//trim(file_name))   
                        
         end if check_open_fine
            
      end if no_input
   
   end subroutine open_wims9plus_input
   
   ! --------------------------------------------------------------------------

end subroutine test_radiation_materials_read_wims9plus_base
