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

module radiation_materials_create
   
   !!< Module with procedures for the creation of particle_radmat types. This includes 
   !!< deducing the size of the arrays from the spud options dictionary, allocating, 
   !!< setting the option path as needed and zeroing
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud

   use radiation_materials_data_types 
   
   implicit none
   
   private 
             
   public :: allocate, &
             zero, &
             set_option_path_and_name, &
             create

   interface allocate
      module procedure particle_radmat_allocate, &
                       particle_radmat_size_from_options_allocate, &
                       particle_radmat_size_allocate, &
                       delayed_lambda_spectrum_allocate, &
                       dataset_radmat_allocate, &
                       physical_radmat_allocate, &
                       radmat_allocate
   end interface allocate

   interface zero
      module procedure particle_radmat_zero, &
                       delayed_lambda_spectrum_zero, &
                       dataset_radmat_zero, &
                       physical_radmat_zero, &
                       radmat_zero
   end interface zero

   interface set_option_path_and_name
      module procedure particle_radmat_set_option_path_and_name, &
                       dataset_radmat_set_option_path_and_name, &
                       physical_radmat_set_option_path_and_name
   end interface set_option_path_and_name

   interface create
      module procedure particle_radmat_create_from_options, &
                       particle_radmat_size_create_from_options
   end interface create
                      
contains

   ! --------------------------------------------------------------------------

   subroutine particle_radmat_create_from_options(particle_radmat, &
                                                  particle_option_path, &
                                                  particle_name)
      
      !!< Create from the spud options path the particle_radmat type. 
      !!< This involves determining the size arrays, allocating, 
      !!< setting the options path and zeroing
      
      type(particle_radmat_type), intent(inout) :: particle_radmat
      character(len=*), intent(in) :: particle_option_path
      character(len=*), intent(in) :: particle_name
      
      ewrite(1,*) 'Create radiation particle_radmat for ',trim(particle_name)
      
      call create(particle_radmat%particle_radmat_size, &
                  trim(particle_option_path), &
                  trim(particle_name))     
         
      call allocate(particle_radmat)
      
      call set_option_path_and_name(particle_radmat, &
                                    trim(particle_option_path))
         
      call zero(particle_radmat)
      
      ! set the maximum number of scatter moments
      particle_radmat%max_number_of_scatter_moments = maxval(particle_radmat%particle_radmat_size%number_of_scatter_moments)

      particle_radmat%created = .true.
      
      ewrite(1,*) 'Finished creating radiation particle_radmat for ',trim(particle_option_path) 
      
   end subroutine particle_radmat_create_from_options
   
   ! --------------------------------------------------------------------------

   subroutine particle_radmat_size_create_from_options(particle_radmat_size, &
                                                       particle_option_path, &
                                                       particle_name)
      
      !!< Create from the options dictionary the particle_radmat_size type --> allocate, zero and set
      
      type(particle_radmat_size_type), intent(inout) :: particle_radmat_size
      character(len=*), intent(in) :: particle_option_path
      character(len=*), intent(in) :: particle_name
      
      ! local variables
      integer :: dmat,pmat,i,status,g_set
      integer :: number_of_energy_group_set
      integer :: number_of_energy_groups_g_set
      integer :: interpolation_dimensions
      integer :: interpolation_values_shape(2)
      integer, dimension(:), allocatable :: number_interpolation_values
      character(len=OPTION_PATH_LEN) :: dataset_radmat_option_path
      character(len=OPTION_PATH_LEN) :: physical_radmat_option_path
      character(len=OPTION_PATH_LEN) :: interpolation_dim_path
      character(len=OPTION_PATH_LEN) :: energy_group_set_path

      ewrite(1,*) 'Create radiation particle_radmat_size data type for ',trim(particle_name)      

      ! allocate the sizes array from the options given by the option path
      call allocate(particle_radmat_size, &
                    particle_option_path)
                    
      ! now set the particle_radmat_size data type values  
      particle_radmat_size%number_of_scatter_moments  = 0  
      particle_radmat_size%number_of_radmats          = 0                            
      particle_radmat_size%number_of_radmats_base     = 0 
      particle_radmat_size%number_of_physical_radmats = 0
               
      particle_radmat_size%number_of_radmats_base(1)  = 1
      
      dataset_loop: do dmat = 1,particle_radmat_size%total_number_dataset_radmats 

         dataset_radmat_option_path = &
         &trim(particle_option_path)//"/material_data_set["//int2str(dmat - 1)//"]"

         ! the number of scatter moments is data set dependent   
         call get_option(trim(dataset_radmat_option_path)//"/from_file/number_of_scatter_moments",particle_radmat_size%number_of_scatter_moments(dmat))  
                 
         ! deduce the number of physical radmats within this data set
         particle_radmat_size%number_of_physical_radmats(dmat) = option_count(trim(dataset_radmat_option_path)//"/from_file/physical_material")
                                   
         physical_radmat_loop: do pmat = 1,particle_radmat_size%number_of_physical_radmats(dmat)  
      
            physical_radmat_option_path = &
            &trim(dataset_radmat_option_path)//"/from_file/physical_material["//int2str(pmat - 1)//"]"
    
            ! deduce the number of physical radmats interpolation dimensions within this physical radmat
            interpolation_dimensions = option_count(trim(physical_radmat_option_path)//"/interpolation_dimension")

            allocate(number_interpolation_values(interpolation_dimensions),STAT=status)
            
            if (status /= 0) FLAbort("Issue allocating memory for number_interpolation_values in set_particle_radmat_size_from_options")
            
            ! loop each interpolation dimension and get the number of interpolation values per dimension
            interpolation_dim_loop: do i = 1,interpolation_dimensions
                   
              ! get the interpolation dimension path
              interpolation_dim_path = trim(physical_radmat_option_path)//"/interpolation_dimension["//int2str(i - 1)//"]"
              
              ! get the shape of the float array of interpolation values
              interpolation_values_shape = option_shape(trim(interpolation_dim_path)//"/interpolation_values")
              
              ! set the number of interpolation values from the shape array
              number_interpolation_values(i) = interpolation_values_shape(1)

            end do interpolation_dim_loop
                              
            particle_radmat_size%number_of_radmats(sum(particle_radmat_size%number_of_physical_radmats(1:dmat-1)) + pmat) = product(number_interpolation_values)
                        
            deallocate(number_interpolation_values)
                  
         end do physical_radmat_loop
                       
      end do dataset_loop

      ! get the number of energy groups via summing the number within each energy group set

      ! deduce the number of energy group sets
      number_of_energy_group_set = option_count(trim(particle_option_path)//'/energy_group_set')
      
      particle_radmat_size%number_of_energy_groups = 0
      
      energy_group_set_loop: do g_set = 1,number_of_energy_group_set
            
         ! set the energy_group_set path
         energy_group_set_path = trim(particle_option_path)//'/energy_group_set['//int2str(g_set - 1)//']'
            
         ! get the number_energy_groups within this set
         call get_option(trim(energy_group_set_path)//'/number_of_energy_groups',number_of_energy_groups_g_set)         
         
         particle_radmat_size%number_of_energy_groups = particle_radmat_size%number_of_energy_groups +  &
                                                        number_of_energy_groups_g_set
         
      end do energy_group_set_loop
             
      ! get the number of delayed groups for this particle type
      delayed_if: if (have_option(trim(particle_option_path)//"/delayed_neutron_precursor")) then
                                                                         
         ! get the number of delayed groups and energy groups
         call get_option(trim(particle_option_path)//&
                        &"/delayed_neutron_precursor/number_delayed_neutron_precursor_groups",particle_radmat_size%number_of_delayed_groups)
         
      else
         
         particle_radmat_size%number_of_delayed_groups = 0
         
      end if delayed_if  
      
      particle_radmat_size%size_set = .true.
                                    
   end subroutine particle_radmat_size_create_from_options

   ! --------------------------------------------------------------------------
   
   subroutine particle_radmat_size_from_options_allocate(particle_radmat_size, &
                                                         particle_option_path)
   
      !!< Allocate the particle_radmat_size type from the options given by the path
            
      type(particle_radmat_size_type), intent(inout) :: particle_radmat_size
      character(len=*), intent(in) :: particle_option_path
      
      ! local variables
      integer :: dmat,pmat,i,status
      integer :: interpolation_dimensions
      integer :: interpolation_values_shape(2)
      integer, dimension(:), allocatable :: number_interpolation_values
      character(len=OPTION_PATH_LEN) :: dataset_radmat_option_path
      character(len=OPTION_PATH_LEN) :: physical_radmat_option_path
      character(len=OPTION_PATH_LEN) :: interpolation_dim_path
      
      ewrite(1,*) 'Allocate radiation particle_radmat_size from options for ',trim(particle_option_path)
      
      ! deduce the sizes of arrays to allocate
      particle_radmat_size%total_number_dataset_radmats  = option_count(trim(particle_option_path)//"/material_data_set")                             
      particle_radmat_size%total_number_radmats          = 0
      particle_radmat_size%total_number_physical_radmats = 0

      dataset_loop: do dmat = 1,particle_radmat_size%total_number_dataset_radmats 

         dataset_radmat_option_path = &
         &trim(particle_option_path)//"/material_data_set["//int2str(dmat - 1)//"]"
    
         particle_radmat_size%total_number_physical_radmats = particle_radmat_size%total_number_physical_radmats + &  
                                     option_count(trim(dataset_radmat_option_path)//"/from_file/physical_material")
                      
         physical_radmat_loop: do pmat = 1,option_count(trim(dataset_radmat_option_path)//"/from_file/physical_material")
      
            physical_radmat_option_path = &
            &trim(dataset_radmat_option_path)//"/from_file/physical_material["//int2str(pmat - 1)//"]"

            ! deduce the number of physical radmats interpolation dimensions within this physical radmat
            interpolation_dimensions = option_count(trim(physical_radmat_option_path)//"/interpolation_dimension")

            allocate(number_interpolation_values(interpolation_dimensions),STAT=status)
            
            if (status /= 0) FLAbort("Issue allocating memory for number_interpolation_values in set_particle_radmat_size_from_options")
            
            ! loop each interpolation dimension and get the number of interpolation values per dimension
            interpolation_dim_loop: do i = 1,interpolation_dimensions
                   
              ! get the interpolation dimension path
              interpolation_dim_path = trim(physical_radmat_option_path)//"/interpolation_dimension["//int2str(i - 1)//"]"

              ! get the shape of the float array of interpolation values
              interpolation_values_shape = option_shape(trim(interpolation_dim_path)//"/interpolation_values")
              
              ! set the number of interpolation values from the shape array
              number_interpolation_values(i) = interpolation_values_shape(1)
     
            end do interpolation_dim_loop
                  
            particle_radmat_size%total_number_radmats = particle_radmat_size%total_number_radmats + product(number_interpolation_values)
                                
            deallocate(number_interpolation_values)
                  
         end do physical_radmat_loop
                       
      end do dataset_loop

      ! allocate the arrays sizes now knowing the sizes of them
      call allocate(particle_radmat_size)  
       
   end subroutine particle_radmat_size_from_options_allocate
   
   ! --------------------------------------------------------------------------
   
   subroutine particle_radmat_size_allocate(particle_radmat_size)
   
      !!< Allocate the arrays that will hold info about the size of the particle_radmat type
      
      type(particle_radmat_size_type), intent(inout) :: particle_radmat_size
      
      ! local variable
      integer :: status
      
      ewrite(1,*) 'Allocate radiation particle_radmat_size'
            
      allocate(particle_radmat_size%number_of_physical_radmats(particle_radmat_size%total_number_dataset_radmats),STAT=status)
      if (status /= 0) FLAbort("Issue allocating memory for particle_radmat_size%number_of_physical_radmats in particle_radmat_size_allocate")
      
      allocate(particle_radmat_size%number_of_radmats(particle_radmat_size%total_number_physical_radmats),STAT=status)
      if (status /= 0) FLAbort("Issue allocating memory for particle_radmat_size%number_of_radmats in particle_radmat_size_allocate")
      
      allocate(particle_radmat_size%number_of_radmats_base(particle_radmat_size%total_number_dataset_radmats),STAT=status)
      if (status /= 0) FLAbort("Issue allocating memory for particle_radmat_size%number_of_radmats_base in particle_radmat_size_allocate")
      
      ! the number of scatter moments is data set dependent  
      allocate(particle_radmat_size%number_of_scatter_moments(particle_radmat_size%total_number_dataset_radmats),STAT=status)
      if (status /= 0) FLAbort("Issue allocating memory for particle_radmat_size%number_of_scatter_moments in particle_radmat_size_allocate")
      
   end subroutine particle_radmat_size_allocate
   
   ! --------------------------------------------------------------------------

   subroutine particle_radmat_allocate(particle_radmat)
      
      !!< Allocate a neutral particle radmat type, particle_radmat

      type(particle_radmat_type), intent(inout) :: particle_radmat
      
      ! local variables
      integer :: status,dmat,total_number_dataset_radmats
      
      ewrite(1,*) 'Allocate particle_radmat'
      
      ! first check that the particle_radmat%particle_radmat_size arrays have been set
      not_set: if (.not. particle_radmat%particle_radmat_size%size_set) then
      
         FLAbort('Cannot allocate a particle_radmat type if the size type particle_radmat%particle_radmat_size is not already set')
      
      end if not_set
      
      total_number_dataset_radmats = particle_radmat%particle_radmat_size%total_number_dataset_radmats
                                    
      allocate(particle_radmat%dataset_radmats(total_number_dataset_radmats),STAT=status)

      if (status /= 0) FLAbort("Issue allocating memory for particle_radmat%dataset_radmats")

      dataset_loop: do dmat = 1,total_number_dataset_radmats 
                 
         call allocate(particle_radmat%dataset_radmats(dmat), &
                       particle_radmat%particle_radmat_size, &
                       dmat)
      
      end do dataset_loop
      
      call allocate(particle_radmat%delayed_lambda_spectrum, &
                    particle_radmat%particle_radmat_size%number_of_delayed_groups, &
                    particle_radmat%particle_radmat_size%number_of_energy_groups)
                        
   end subroutine particle_radmat_allocate

   ! --------------------------------------------------------------------------
   
   subroutine delayed_lambda_spectrum_allocate(delayed_lambda_spectrum, &
                                               number_of_delayed_groups, &
                                               number_of_energy_groups)
      
      !!< Allocate the delayed lambda and spectrum type
      
      type(delayed_lambda_spectrum_type), intent(inout) :: delayed_lambda_spectrum
      integer, intent(in) :: number_of_delayed_groups
      integer, intent(in) :: number_of_energy_groups
      
      ! local variable
      integer :: status
      
      allocate(delayed_lambda_spectrum%lambda(number_of_delayed_groups),STAT=status)
         
      if (status /= 0) FLAbort("Issue allocating memory for delayed_lambda_spectrum%lambda")
         
      allocate(delayed_lambda_spectrum%spectrum(number_of_energy_groups,number_of_delayed_groups),STAT=status)
         
      if (status /= 0) FLAbort("Issue allocating memory for delayed_lambda_spectrum%spectrum")
                                
   end subroutine delayed_lambda_spectrum_allocate
   
   ! --------------------------------------------------------------------------

   subroutine dataset_radmat_allocate(dataset_radmat, &
                                      particle_radmat_size, &
                                      dmat)

      type(dataset_radmat_type), intent(inout) :: dataset_radmat      
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size
      integer, intent(in) :: dmat
      
      ! local variables
      integer :: pmat,status,number_of_physical_radmats
                                    
      number_of_physical_radmats = particle_radmat_size%number_of_physical_radmats(dmat)

      allocate(dataset_radmat%physical_radmats(number_of_physical_radmats),STAT=status)
      
      if (status /= 0) FLAbort("Issue allocating memory for dataset_radmat%physical_radmats")
      
      physical_radmat_loop: do pmat = 1,number_of_physical_radmats  
                     
         call allocate(dataset_radmat%physical_radmats(pmat), &
                       particle_radmat_size, &
                       dmat, &
                       pmat)
      
      end do physical_radmat_loop
                                          
   end subroutine dataset_radmat_allocate  
   
   ! --------------------------------------------------------------------------
   
   subroutine physical_radmat_allocate(physical_radmat, &
                                       particle_radmat_size, &
                                       dmat, &
                                       pmat)
      
      !!< Allocate a physical_radmat type
      
      type(physical_radmat_type), intent(inout) :: physical_radmat
      type(particle_radmat_size_type), intent(in) :: particle_radmat_size
      integer, intent(in) :: dmat      
      integer, intent(in) :: pmat
      
      ! local variables
      integer :: rmat,status,number_of_radmats
      
      number_of_radmats = particle_radmat_size%number_of_radmats(particle_radmat_size%number_of_radmats_base(dmat) - 1 + pmat)
            
      allocate(physical_radmat%radmats(number_of_radmats),STAT=status)
      
      if (status /= 0) FLAbort("Issue allocating memory for physical_radmat%radmats")
      
      radmat_loop: do rmat = 1,number_of_radmats
         
         call allocate(physical_radmat%radmats(rmat), &
                       particle_radmat_size%number_of_energy_groups, &
                       particle_radmat_size%number_of_scatter_moments(dmat), &
                       particle_radmat_size%number_of_delayed_groups, &
                       allocate_all = .true.)
      
      end do radmat_loop
        
   end subroutine physical_radmat_allocate
   
   ! --------------------------------------------------------------------------

   subroutine radmat_allocate(radmat, &
                              number_of_energy_groups, &
                              number_of_scatter_moments, &
                              number_of_delayed_groups, &
                              allocate_all, &
                              allocate_total, &
                              allocate_absorption, &
                              allocate_scatter, &
                              allocate_removal, &
                              allocate_transport, &
                              allocate_diffusion, &
                              allocate_fission, &
                              allocate_production, &
                              allocate_power, &
                              allocate_energy_released_per_fission, &
                              allocate_particle_released_per_fission, &
                              allocate_prompt_spectrum, &
                              allocate_velocity, &
                              allocate_beta)
      
      !!< Allocate the arrays within a radmat type as needed
      
      type(radmat_type), intent(inout) :: radmat
      integer, intent(in) :: number_of_energy_groups
      integer, intent(in) :: number_of_scatter_moments      
      integer, intent(in) :: number_of_delayed_groups
      logical, intent(in), optional :: allocate_all
      logical, intent(in), optional :: allocate_total
      logical, intent(in), optional :: allocate_absorption
      logical, intent(in), optional :: allocate_scatter
      logical, intent(in), optional :: allocate_removal
      logical, intent(in), optional :: allocate_transport
      logical, intent(in), optional :: allocate_diffusion
      logical, intent(in), optional :: allocate_fission
      logical, intent(in), optional :: allocate_production
      logical, intent(in), optional :: allocate_power
      logical, intent(in), optional :: allocate_energy_released_per_fission
      logical, intent(in), optional :: allocate_particle_released_per_fission
      logical, intent(in), optional :: allocate_prompt_spectrum
      logical, intent(in), optional :: allocate_velocity
      logical, intent(in), optional :: allocate_beta
      
      ! local variable
      integer :: status
      
      total: if (allocate_all .or. allocate_total) then
         allocate(radmat%total(number_of_energy_groups),STAT=status)
         if (status /= 0) FLAbort("Issue allocating memory for radmat%total")
      end if total
      
      absorption: if (allocate_all .or. allocate_absorption) then
         allocate(radmat%absorption(number_of_energy_groups),STAT=status)
         if (status /= 0) FLAbort("Issue allocating memory for radmat%absorption")
      end if absorption

      scatter: if (allocate_all .or. allocate_scatter) then      
         allocate(radmat%scatter(number_of_energy_groups,number_of_energy_groups,number_of_scatter_moments),STAT=status)
         if (status /= 0) FLAbort("Issue allocating memory for radmat%scatter")
      end if scatter

      removal: if (allocate_all .or. allocate_removal) then      
         allocate(radmat%removal(number_of_energy_groups,number_of_scatter_moments),STAT=status)
         if (status /= 0) FLAbort("Issue allocating memory for radmat%removal")
      end if removal

      transport: if (allocate_all .or. allocate_transport) then      
         allocate(radmat%transport(number_of_energy_groups,3),STAT=status)
         if (status /= 0) FLAbort("Issue allocating memory for radmat%transport")
      end if transport

      diffusion: if (allocate_all .or. allocate_diffusion) then      
         allocate(radmat%diffusion(number_of_energy_groups,3),STAT=status)
         if (status /= 0) FLAbort("Issue allocating memory for radmat%diffusion")
      end if diffusion 

      fission: if (allocate_all .or. allocate_fission) then      
         allocate(radmat%fission(number_of_energy_groups),STAT=status)
         if (status /= 0) FLAbort("Issue allocating memory for (radmat%fission")
      end if fission

      production: if (allocate_all .or. allocate_production) then      
         allocate(radmat%production(number_of_energy_groups),STAT=status)
         if (status /= 0) FLAbort("Issue allocating memory for radmat%production")
      end if production

      power: if (allocate_all .or. allocate_power) then      
         allocate(radmat%power(number_of_energy_groups),STAT=status)
         if (status /= 0) FLAbort("Issue allocating memory for radmat%power")
      end if power

      energy_released_per_fission: if (allocate_all .or. allocate_energy_released_per_fission) then      
         allocate(radmat%energy_released_per_fission(number_of_energy_groups),STAT=status)
         if (status /= 0) FLAbort("Issue allocating memory for radmat%energy_released_per_fission")
      end if energy_released_per_fission

      particle_released_per_fission: if (allocate_all .or. allocate_particle_released_per_fission) then      
         allocate(radmat%particle_released_per_fission(number_of_energy_groups),STAT=status)
         if (status /= 0) FLAbort("Issue allocating memory for radmat%particle_released_per_fission")
      end if particle_released_per_fission

      prompt_spectrum: if (allocate_all .or. allocate_prompt_spectrum) then      
         allocate(radmat%prompt_spectrum(number_of_energy_groups),STAT=status)
         if (status /= 0) FLAbort("Issue allocating memory for radmat%prompt_spectrum")
      end if prompt_spectrum

      velocity: if (allocate_all .or. allocate_velocity) then      
         allocate(radmat%velocity(number_of_energy_groups),STAT=status)
         if (status /= 0) FLAbort("Issue allocating memory for radmat%velocity")
      end if velocity

      beta: if (allocate_all .or. allocate_beta) then                                                                                                    
         allocate(radmat%beta(number_of_delayed_groups),STAT=status)         
         if (status /= 0) FLAbort("Issue allocating memory for radmat%beta")
      end if beta
                                             
   end subroutine radmat_allocate

   ! --------------------------------------------------------------------------

   subroutine particle_radmat_set_option_path_and_name(particle_radmat, &
                                                       particle_option_path)
      
      !!< Set the option path and name for a neutral particle radmat type, particle_radmat

      type(particle_radmat_type), intent(inout) :: particle_radmat
      character(len=*), intent(in) :: particle_option_path
      
      ! local variables
      integer :: dmat
      character(len=OPTION_PATH_LEN) :: dataset_radmat_option_path
      
      ewrite(1,*) 'Set radiation particle_radmat option path for ',trim(particle_option_path)
            
      particle_radmat%option_path = particle_option_path 

      call get_option(trim(particle_radmat%option_path)//'/name',particle_radmat%name)
      
      dataset_loop: do dmat = 1,size(particle_radmat%dataset_radmats) 
         
         ! - 1 needed as options count from 0
         dataset_radmat_option_path = &
         &trim(particle_radmat%option_path)//"/material_data_set["//int2str(dmat - 1)//"]"
                 
         call set_option_path_and_name(particle_radmat%dataset_radmats(dmat), &
                                       trim(dataset_radmat_option_path))
      
      end do dataset_loop
                              
   end subroutine particle_radmat_set_option_path_and_name

   ! --------------------------------------------------------------------------

   subroutine dataset_radmat_set_option_path_and_name(dataset_radmat, &
                                                      dataset_radmat_option_path)

      !!< Set the option path and name for a dataset_radmat type

      type(dataset_radmat_type), intent(inout) :: dataset_radmat      
      character(len=*), intent(in) :: dataset_radmat_option_path
         
      ! local variables
      integer :: pmat
      character(len=OPTION_PATH_LEN) :: physical_radmat_option_path
                                    
      dataset_radmat%option_path = dataset_radmat_option_path 

      call get_option(trim(dataset_radmat%option_path)//'/from_file/file_name',dataset_radmat%file_name)

      call get_option(trim(dataset_radmat%option_path)//'/name',dataset_radmat%name)
      
      physical_radmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)  
         
         ! - 1 needed as options count from 0  
         physical_radmat_option_path = &
         &trim(dataset_radmat%option_path)//"/from_file/physical_material["//int2str(pmat - 1)//"]"
                     
         call set_option_path_and_name(dataset_radmat%physical_radmats(pmat), &
                                       trim(physical_radmat_option_path))
      
      end do physical_radmat_loop
                                          
   end subroutine dataset_radmat_set_option_path_and_name  
   
   ! --------------------------------------------------------------------------
   
   subroutine physical_radmat_set_option_path_and_name(physical_radmat, &
                                                       physical_radmat_option_path)

      !!< Set the option path and name for a physical_radmat type
      
      type(physical_radmat_type), intent(inout) :: physical_radmat
      character(len=*), intent(in) :: physical_radmat_option_path
      
      physical_radmat%option_path = physical_radmat_option_path 
      
      call get_option(trim(physical_radmat%option_path)//'/name',physical_radmat%name)
        
   end subroutine physical_radmat_set_option_path_and_name
   
   ! --------------------------------------------------------------------------
   
   subroutine particle_radmat_zero(particle_radmat)
      
      !!< Zero the material data arrays associated with this particle_radmat type
      
      type(particle_radmat_type), intent(inout) :: particle_radmat
      
      ! local variable
      integer :: dmat
      
      ewrite(1,*) 'Zero radiation particle_radmat'
      
      dataset_loop: do dmat = 1,size(particle_radmat%dataset_radmats)
                 
         call zero(particle_radmat%dataset_radmats(dmat))
      
      end do dataset_loop 
      
      call zero(particle_radmat%delayed_lambda_spectrum)
        
   end subroutine particle_radmat_zero

   ! --------------------------------------------------------------------------
   
   subroutine delayed_lambda_spectrum_zero(delayed_lambda_spectrum)
   
      !!< Zero the lambda and spectrum arrays of this delayed_lambda_spectrum type  
      !!< only if they are allocated
      
      type(delayed_lambda_spectrum_type), intent(inout) :: delayed_lambda_spectrum
      
      if (allocated(delayed_lambda_spectrum%lambda))   delayed_lambda_spectrum%lambda   = 0.0
      if (allocated(delayed_lambda_spectrum%spectrum)) delayed_lambda_spectrum%spectrum = 0.0
   
   end subroutine delayed_lambda_spectrum_zero
      
   ! --------------------------------------------------------------------------
   
   subroutine dataset_radmat_zero(dataset_radmat)
   
      !!< Zero the material data associated with this dataset_radmat type
      
      type(dataset_radmat_type), intent(inout) :: dataset_radmat
      
      ! local variable
      integer :: pmat
      
      physical_radmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)  
                     
         call zero(dataset_radmat%physical_radmats(pmat))
      
      end do physical_radmat_loop      
            
   end subroutine dataset_radmat_zero
      
   ! --------------------------------------------------------------------------
   
   subroutine physical_radmat_zero(physical_radmat)
   
      !!< Zero the material data associated with this physical_radmat type
      
      type(physical_radmat_type), intent(inout) :: physical_radmat
      
      ! local variable
      integer :: rmat
      
      radmat_loop: do rmat = 1,size(physical_radmat%radmats)  
                     
         call zero(physical_radmat%radmats(rmat))
      
      end do radmat_loop      
            
   end subroutine physical_radmat_zero
      
   ! --------------------------------------------------------------------------
   
   subroutine radmat_zero(radmat)
   
      !!< Zero the material data associated with this radmat type
      !!< only if they are allocated
      
      type(radmat_type), intent(inout) :: radmat
            
      if (allocated(radmat%total))                         radmat%total                         = 0.0
      if (allocated(radmat%absorption))                    radmat%absorption                    = 0.0
      if (allocated(radmat%scatter))                       radmat%scatter                       = 0.0 
      if (allocated(radmat%removal))                       radmat%removal                       = 0.0 
      if (allocated(radmat%transport))                     radmat%transport                     = 0.0
      if (allocated(radmat%diffusion))                     radmat%diffusion                     = 0.0   
      if (allocated(radmat%fission))                       radmat%fission                       = 0.0
      if (allocated(radmat%production))                    radmat%production                    = 0.0
      if (allocated(radmat%power))                         radmat%power                         = 0.0
      if (allocated(radmat%energy_released_per_fission))   radmat%energy_released_per_fission   = 0.0
      if (allocated(radmat%particle_released_per_fission)) radmat%particle_released_per_fission = 0.0
      if (allocated(radmat%prompt_spectrum))               radmat%prompt_spectrum               = 0.0
      if (allocated(radmat%velocity))                      radmat%velocity                      = 0.0     
      if (allocated(radmat%beta))                          radmat%beta                          = 0.0
                  
   end subroutine radmat_zero
      
   ! --------------------------------------------------------------------------
   
end module radiation_materials_create   
   
