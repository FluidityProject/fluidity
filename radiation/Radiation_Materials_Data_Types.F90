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

module radiation_materials_data_types
   
   !!< define the radiation materials data type used for each particle type
   
   ! this currently uses dynamic allocatable embedded data types, but may be more flexible if
   ! using dynamic pointer embedded data types with reference counting
   
   use global_parameters, only : OPTION_PATH_LEN
   
   implicit none
   
   private 

   public :: particle_radmat_type, &
             delayed_lambda_spectrum_type, &
             dataset_radmat_type, &
             physical_radmat_type, &
             radmat_type, &
             particle_radmat_size_type

   ! the basic radiation material type (eg. uranium oxide at a certain temperature)
   type radmat_type         
      !! Cross-sections
      real, dimension(:),     allocatable :: total 
      real, dimension(:),     allocatable :: absorption
      real, dimension(:,:,:), allocatable :: scatter 
      real, dimension(:,:),   allocatable :: removal 
      real, dimension(:,:),   allocatable :: transport
      real, dimension(:,:),   allocatable :: diffusion   
      real, dimension(:),     allocatable :: fission
      real, dimension(:),     allocatable :: production
      real, dimension(:),     allocatable :: power
      real, dimension(:),     allocatable :: energy_released_per_fission
      real, dimension(:),     allocatable :: particle_released_per_fission      
      real, dimension(:),     allocatable :: prompt_spectrum
      real, dimension(:),     allocatable :: velocity
      
      !! the delayed beta are radmat dependent as it could change with burnup 
      real, dimension(:),   allocatable :: beta
      
      !! logical flags for whether data has been set in or not 
      !! set means it has been allocated, initialised and then values read in or calculated
      logical :: total_set                         = .false.
      logical :: absorption_set                    = .false.
      logical :: scatter_set                       = .false.
      logical :: removal_set                       = .false.
      logical, dimension(3) :: transport_set       = .false.
      logical, dimension(3) :: diffusion_set       = .false.
      logical :: fission_set                       = .false.
      logical :: production_set                    = .false. 
      logical :: power_set                         = .false.
      logical :: energy_released_per_fission_set   = .false.
      logical :: particle_released_per_fission_set = .false.      
      logical :: prompt_spectrum_set               = .false.
      logical :: velocity_set                      = .false.
      logical:: beta_set                           = .false.            
   end type radmat_type

   
   ! a physical radiation material (eg. uranium oxide tabulated at a range of temperatures)
   type physical_radmat_type                      
      !! the physical_radmat_type name
      character(len=OPTION_PATH_LEN) :: name="/uninitialised_name/"
      !! path to options in the options tree
      character(len=OPTION_PATH_LEN) :: option_path="/uninitialised_path/"
      !! the individual radiation materials within this physical material
      type(radmat_type), dimension(:), allocatable :: radmats
   end type physical_radmat_type

   
   ! a data set of radiation physical materials (eg. uranium oxide, steel and water)
   type dataset_radmat_type                              
      !! the dataset_radmat_type FILE name
      character(len=OPTION_PATH_LEN) :: file_name="/uninitialised_file_name/"
      !! path to options in the options tree
      character(len=OPTION_PATH_LEN) :: option_path="/uninitialised_path/"
      !! the physical radmats within this dataset
      type(physical_radmat_type), dimension(:), allocatable :: physical_radmats
   end type dataset_radmat_type

   
   ! the delayed precursor lambda and spectrum data that is considered the same regardless of the physical material
   type delayed_lambda_spectrum_type
      !! the delayed lambda and spectrum
      real, dimension(:),   allocatable :: lambda
      real, dimension(:,:), allocatable :: spectrum  

      !! logical flags for whether data has been set in or not 
      !! set means it has been allocated, initialised and then values read in or calculated
      logical :: lambda_set = .false.
      logical :: spectrum_set = .false.         
   end type delayed_lambda_spectrum_type


   ! the size of the particle radmat data type needed for initialisation from either options or other
   type particle_radmat_size_type
      integer :: total_number_radmats
      integer :: total_number_physical_radmats
      integer :: total_number_dataset_radmats      
      integer, dimension(:), allocatable :: number_of_physical_radmats ! data set dependent
      integer, dimension(:), allocatable :: number_of_radmats ! data set and physical mat dependent
      integer, dimension(:), allocatable :: number_of_radmats_base ! base pointer of above for start of each dataset list of physical mat     

      integer, dimension(:), allocatable :: number_of_scatter_moments ! data set dependent
      integer :: number_of_energy_groups
      integer :: number_of_delayed_groups      

      logical :: size_set = .false.
   end type particle_radmat_size_type

   
   ! the collection of radiation data sets associated with a particle type
   type particle_radmat_type            
      !! the particle_radmat_type name
      character(len=OPTION_PATH_LEN) :: name="/uninitialised_name/"
      !! path to options in the options tree
      character(len=OPTION_PATH_LEN) :: option_path="/uninitialised_path/"
      !! the size tags of the particle_radmat
      type(particle_radmat_size_type) ::  particle_radmat_size  
      !! the datasets associated with this neutral particle object
      type(dataset_radmat_type), dimension(:), allocatable :: dataset_radmats
      !! the delayed lambda and spectrum data for this neutral particle object
      type(delayed_lambda_spectrum_type) :: delayed_lambda_spectrum
      !! The maximum number of scatter moments within the datasets 
      integer :: max_number_of_scatter_moments
      !! a flag to say that this type has been created
      logical :: created = .false.
      !! a flag to say that this type has been read in
      logical :: readin = .false.      
   end type particle_radmat_type

end module radiation_materials_data_types
