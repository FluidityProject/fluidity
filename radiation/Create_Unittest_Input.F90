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

module create_unittest_input

   !!< This module is only called from the unittests to load hard wired materials for comparing against

   use radiation_materials_data_types
   use radiation_materials_create, only: allocate,zero 
   
   implicit none
   
   private
   
   public :: create_particle_radmat_input1
   
contains

   ! --------------------------------------------------------------------------

   subroutine create_particle_radmat_input1(particle_radmat)

      !!< Create the particle_radmat associated with the input file rad_input_test_dir//'test_radiation_materials_format_radmats_input1'
      !!< this effectively hard codes in the values so that they can be compared against the read in value.
      !!< When used this implicity is also testing the allocate, zero and destroy procedures.
      
      type(particle_radmat_type), intent(inout) :: particle_radmat
      
      ! local variables
      integer :: dmat,pmat,rmat,igrp
      
      particle_radmat%particle_radmat_size%total_number_dataset_radmats  = 1
      particle_radmat%particle_radmat_size%total_number_physical_radmats = 1
      call allocate(particle_radmat%particle_radmat_size)

      particle_radmat%particle_radmat_size%total_number_radmats       = 10  
      particle_radmat%particle_radmat_size%number_of_energy_groups    = 2
      particle_radmat%particle_radmat_size%number_of_delayed_groups   = 6
      particle_radmat%particle_radmat_size%number_of_scatter_moments  = (/2/)
      particle_radmat%particle_radmat_size%number_of_physical_radmats = (/1/)
      particle_radmat%particle_radmat_size%number_of_radmats          = (/10/)
      particle_radmat%particle_radmat_size%number_of_radmats_base     = (/1/)   
      particle_radmat%particle_radmat_size%size_set                   = .true.

      call allocate(particle_radmat)

      call zero(particle_radmat)
      
      ! set the option paths
      particle_radmat%option_path = '/embedded_models/radiation/particle_type[0]'
      particle_radmat%dataset_radmats(1)%option_path = '/embedded_models/radiation/particle_type[0]/material_data_set[0]'
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%option_path =  '/embedded_models/radiation/particle_type[0]/material_data_set[0]/from_file/physical_material[0]'

      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%velocity = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%velocity = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%velocity = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%velocity = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%velocity = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%velocity = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%velocity = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%velocity = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%velocity = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%velocity = (/0.0,0.0/)    

      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%beta = &
                    & (/2.211e-4,1.4673e-3,1.3132e-3,2.6465e-3,7.705e-4,2.814e-4/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%beta = &
                    & (/2.211e-4,1.4673e-3,1.3132e-3,2.6465e-3,7.705e-4,2.814e-4/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%beta = &
                    & (/2.211e-4,1.4673e-3,1.3132e-3,2.6465e-3,7.705e-4,2.814e-4/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%beta = &
                    & (/2.211e-4,1.4673e-3,1.3132e-3,2.6465e-3,7.705e-4,2.814e-4/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%beta = &
                    & (/2.211e-4,1.4673e-3,1.3132e-3,2.6465e-3,7.705e-4,2.814e-4/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%beta = &
                    & (/2.211e-4,1.4673e-3,1.3132e-3,2.6465e-3,7.705e-4,2.814e-4/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%beta = &
                    & (/2.211e-4,1.4673e-3,1.3132e-3,2.6465e-3,7.705e-4,2.814e-4/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%beta = &
                    & (/2.211e-4,1.4673e-3,1.3132e-3,2.6465e-3,7.705e-4,2.814e-4/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%beta = &
                    & (/2.211e-4,1.4673e-3,1.3132e-3,2.6465e-3,7.705e-4,2.814e-4/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%beta = &
                    & (/2.211e-4,1.4673e-3,1.3132e-3,2.6465e-3,7.705e-4,2.814e-4/) 

      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%energy_released_per_fission = (/1.0,2.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%energy_released_per_fission = (/10.0,20.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%energy_released_per_fission = (/100.0,200.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%energy_released_per_fission = (/1000.0,2000.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%energy_released_per_fission = (/10000.0,20000.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%energy_released_per_fission = (/10000.0,20000.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%energy_released_per_fission = (/1000.0,2000.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%energy_released_per_fission = (/100.0,200.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%energy_released_per_fission = (/10.0,20.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%energy_released_per_fission = (/1.0,2.0/) 

      particle_radmat%delayed_lambda_spectrum%lambda = (/0.0124,0.0305,0.111,0.301,1.14,3.01/)  
      particle_radmat%delayed_lambda_spectrum%spectrum(1,:) = (/1.0,1.0,1.0,1.0,1.0,1.0/)  
      particle_radmat%delayed_lambda_spectrum%spectrum(2,:) = (/0.0,0.0,0.0,0.0,0.0,0.0/)  

      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%power = (/1.0*1.548450E-04,2.0*3.373658E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%power = (/10.0*1.624533E-04,20.0*3.539991E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%power = (/100.0*1.570189E-04,200.0*3.421179E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%power = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%power = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%power = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%power = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%power = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%power = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%power = (/0.0,0.0/) 

      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%absorption = (/7.970915E-04,4.262647E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%absorption = (/8.360552E-04,4.469902E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%absorption = (/8.082244E-04,4.321858E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%absorption = (/8.074012E-06,2.305797E-04/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%absorption = (/6.172237E-04,4.237469E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%absorption = (/7.957230E-06,1.608090E-04/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%absorption = (/8.110421E-06,1.643668E-04/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%absorption = (/8.001987E-06,1.618206E-04/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%absorption = (/1.739165E-03,3.923710E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%absorption = (/7.418196E-06,2.110243E-04/) 

      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%fission = (/1.548450E-04,3.373658E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%fission = (/1.624533E-04,3.539991E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%fission = (/1.570189E-04,3.421179E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%fission = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%fission = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%fission = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%fission = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%fission = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%fission = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%fission = (/0.0,0.0/) 

      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%particle_released_per_fission = (/2.447095E+00,2.438252E+00/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%particle_released_per_fission = (/2.447095E+00,2.438252E+00/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%particle_released_per_fission = (/2.447095E+00,2.438252E+00/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%particle_released_per_fission = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%particle_released_per_fission = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%particle_released_per_fission = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%particle_released_per_fission = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%particle_released_per_fission = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%particle_released_per_fission = (/0.0,0.0/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%particle_released_per_fission = (/0.0,0.0/) 

      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%prompt_spectrum = (/1.0,0.0/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%prompt_spectrum = (/1.0,0.0/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%prompt_spectrum = (/1.0,0.0/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%prompt_spectrum = (/0.0,0.0/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%prompt_spectrum = (/0.0,0.0/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%prompt_spectrum = (/0.0,0.0/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%prompt_spectrum = (/0.0,0.0/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%prompt_spectrum = (/0.0,0.0/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%prompt_spectrum = (/0.0,0.0/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%prompt_spectrum = (/0.0,0.0/)   

      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%scatter(1,:,1) = (/2.726214E-01,2.687237E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%scatter(2,:,1) = (/4.883432E-06,3.137488E-01/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%scatter(1,:,2) = (/1.894029E-02,-7.920593E-04/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%scatter(2,:,2) = (/3.370167E-06,5.277442E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%scatter(1,:,1) = (/2.794389E-01,2.753700E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%scatter(2,:,1) = (/5.010780E-06,3.216049E-01/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%scatter(1,:,2) = (/1.941256E-02,-8.116679E-04/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%scatter(2,:,2) = (/3.455521E-06,5.408857E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%scatter(1,:,1) = (/2.745693E-01,2.706226E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%scatter(2,:,1) = (/4.919814E-06,3.159931E-01/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%scatter(1,:,2) = (/1.907523E-02,-7.976619E-04/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%scatter(2,:,2) = (/3.394552E-06,5.314986E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%scatter(1,:,1) = (/3.713450E-01,5.507865E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%scatter(2,:,1) = (/1.398190E-06,4.142296E-01/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%scatter(1,:,2) = (/2.527262E-02,-1.623783E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%scatter(2,:,2) = (/9.958659E-07,-5.366928E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%scatter(1,:,1) = (/3.710299E-01,5.457175E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%scatter(2,:,1) = (/3.527893E-06,4.133705E-01/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%scatter(1,:,2) = (/2.525239E-02,-1.610079E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%scatter(2,:,2) = (/2.513011E-06,-1.202694E-04/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%scatter(1,:,1) = (/2.824581E-01,3.235187E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%scatter(2,:,1) = (/2.112539E-06,3.189540E-01/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%scatter(1,:,2) = (/1.960485E-02,-9.547625E-04/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%scatter(2,:,2) = (/1.504733E-06,-1.349342E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%scatter(1,:,1) = (/2.869382E-01,3.265550E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%scatter(2,:,1) = (/2.215178E-06,3.276513E-01/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%scatter(1,:,2) = (/1.992113E-02,-9.637701E-04/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%scatter(2,:,2) = (/1.577834E-06,-1.267117E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%scatter(1,:,1) = (/2.837670E-01,3.244058E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%scatter(2,:,1) = (/2.141722E-06,3.214269E-01/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%scatter(1,:,2) = (/1.969726E-02,-9.573941E-04/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%scatter(2,:,2) = (/1.525517E-06,-1.325964E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%scatter(1,:,1) = (/3.608233E-01,4.655511E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%scatter(2,:,1) = (/3.665012E-06,4.076633E-01/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%scatter(1,:,2) = (/2.453666E-02,-1.374331E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%scatter(2,:,2) = (/2.610771E-06,4.349558E-04/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%scatter(1,:,1) = (/3.295907E-01,4.783548E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%scatter(2,:,1) = (/1.359111E-06,3.809175E-01/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%scatter(1,:,2) = (/2.248608E-02,-1.410293E-03/) 
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%scatter(2,:,2) = (/9.680496E-07,-4.735995E-03/) 

      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%transport(:,1) = (/2.320508E-01,3.107263E-01/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%transport(:,1) = (/2.378406E-01,3.186084E-01/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%transport(:,1) = (/2.337054E-01,3.129783E-01/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%transport(:,1) = (/3.289176E-01,4.142774E-01/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%transport(:,1) = (/3.283896E-01,4.133193E-01/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%transport(:,1) = (/2.421057E-01,3.169462E-01/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%transport(:,1) = (/2.455837E-01,3.254330E-01/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%transport(:,1) = (/2.431245E-01,3.193613E-01/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%transport(:,1) = (/3.180434E-01,4.069460E-01/)   
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%transport(:,1) = (/2.891242E-01,3.806539E-01/)   
      
      ! form the total, removal, diffusion, transport, production
      dmat_loop: do dmat = 1,size(particle_radmat%dataset_radmats)
            
         pmat_loop: do pmat = 1,size(particle_radmat%dataset_radmats(dmat)%physical_radmats)
            
            rmat_loop: do rmat = 1,size(particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats)
               
               group_loop: do igrp = 1,size(particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%total)
                  
                  particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%total(igrp) = &
                  particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%absorption(igrp) + &
                  sum(particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%scatter(igrp,:,1))
                  
                  particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%removal(igrp,:) = &
                  particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%total(igrp) - &
                  particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%scatter(igrp,igrp,:)
                                    
               end do group_loop
               
               particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%transport(:,2) = &
               particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%transport(:,1) 

               particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%transport(:,3) = &
               particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%transport(:,1) 

               particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%diffusion = &
               1.0/(3.0*particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%transport) 

               particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%production = &
               particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%fission * &
               particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%radmats(rmat)%particle_released_per_fission 
                         
            end do rmat_loop
            
         end do pmat_loop
      
      end do dmat_loop
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%total_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%total_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%total_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%total_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%total_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%total_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%total_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%total_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%total_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%total_set = .true.
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%absorption_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%absorption_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%absorption_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%absorption_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%absorption_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%absorption_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%absorption_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%absorption_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%absorption_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%absorption_set = .true.
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%scatter_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%scatter_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%scatter_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%scatter_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%scatter_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%scatter_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%scatter_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%scatter_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%scatter_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%scatter_set = .true.
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%removal_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%removal_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%removal_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%removal_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%removal_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%removal_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%removal_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%removal_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%removal_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%removal_set = .true.
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%transport_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%transport_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%transport_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%transport_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%transport_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%transport_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%transport_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%transport_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%transport_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%transport_set = .true.
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%diffusion_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%diffusion_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%diffusion_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%diffusion_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%diffusion_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%diffusion_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%diffusion_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%diffusion_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%diffusion_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%diffusion_set = .true.
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%fission_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%fission_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%fission_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%fission_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%fission_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%fission_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%fission_set = .false.
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%production_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%production_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%production_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%production_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%production_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%production_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%production_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%production_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%production_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%production_set = .false.
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%power_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%power_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%power_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%power_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%power_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%power_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%power_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%power_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%power_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%power_set = .false.
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%energy_released_per_fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%energy_released_per_fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%energy_released_per_fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%energy_released_per_fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%energy_released_per_fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%energy_released_per_fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%energy_released_per_fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%energy_released_per_fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%energy_released_per_fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%energy_released_per_fission_set = .true.
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%particle_released_per_fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%particle_released_per_fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%particle_released_per_fission_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%particle_released_per_fission_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%particle_released_per_fission_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%particle_released_per_fission_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%particle_released_per_fission_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%particle_released_per_fission_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%particle_released_per_fission_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%particle_released_per_fission_set = .false.
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%prompt_spectrum_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%prompt_spectrum_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%prompt_spectrum_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%prompt_spectrum_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%prompt_spectrum_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%prompt_spectrum_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%prompt_spectrum_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%prompt_spectrum_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%prompt_spectrum_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%prompt_spectrum_set = .false.
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%velocity_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%velocity_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%velocity_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%velocity_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%velocity_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%velocity_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%velocity_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%velocity_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%velocity_set = .false.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%velocity_set = .false.
       
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(1)%beta_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(2)%beta_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(3)%beta_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(4)%beta_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(5)%beta_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(6)%beta_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(7)%beta_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(8)%beta_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(9)%beta_set = .true.
      particle_radmat%dataset_radmats(1)%physical_radmats(1)%radmats(10)%beta_set = .true.

      particle_radmat%delayed_lambda_spectrum%lambda_set = .true.
      particle_radmat%delayed_lambda_spectrum%spectrum_set = .true.

      particle_radmat%created = .true.
      particle_radmat%readin = .true.      
       
   end subroutine create_particle_radmat_input1 
    
   ! --------------------------------------------------------------------------

end module create_unittest_input
