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

module radiation_materials_interpolation_data_types

   use global_parameters, only : OPTION_PATH_LEN
   
   implicit none
   
   private 

   public :: particle_radmat_ii_type, &
             particle_radmat_ii_size_type, &
             region_id_vele_ii_size_type, &
             region_id_vele_ii_type, &
             dataset_vele_ii_type, &
             physical_radmat_vele_ii_type

   
   ! the interpolation instructions for a particular physical material that is associated with an vele 
   type physical_radmat_vele_ii_type
      ! the physical radmat number 
      integer :: physical_radmat_number
      ! the ii fraction for each dimension of this physical material
      real, dimension(:), allocatable :: fraction
      ! the ii radmat base coordinate
      integer, dimension(:), allocatable :: radmat_base_coordinate
   end type physical_radmat_vele_ii_type

   
   ! the interpolation instructions associated with a particular particle dataset_radmat that is associated with an vele
   type dataset_vele_ii_type
      ! the dataset radmat number
      integer :: dataset_radmat_number 
      ! the interpolation instructions for the physical material that is associated with an vele 
      type(physical_radmat_vele_ii_type) :: physical_radmat_vele_ii
   end type dataset_vele_ii_type

   
   ! the interpolation instructions associated with region id options for a particular vele
   type region_id_vele_ii_type
      ! the interpolation instructions of the dataset mapped via region ids to this vele
      type(dataset_vele_ii_type) :: dataset_vele_ii
   end type region_id_vele_ii_type

   
   ! the interpolation instruction sizes associated with region id mapping
   type region_id_vele_ii_size_type
      ! the fraction size for the ii (also the radmat_base_coordinate size)
      integer :: fraction_size
   end type region_id_vele_ii_size_type


   ! the interpolation instructions sizes 
   type particle_radmat_ii_size_type
      ! the interpolation instruction size for each vele associated with region id mapping
      type(region_id_vele_ii_size_type), dimension(:), allocatable :: region_id_vele_ii_size  
     ! a logical flag to say if these interpolation instruction size have been set
      logical :: size_set=.false.      
   end type particle_radmat_ii_size_type

   
   ! the particle radiation material (radmat) interpolation instructions (ii) type
   type particle_radmat_ii_type
      ! the interpolation instructions for each vele associated with region id mapping
      type(region_id_vele_ii_type), dimension(:), allocatable :: region_id_vele_ii  
      ! the interpolation instructions sizes 
      type(particle_radmat_ii_size_type) :: particle_radmat_ii_size 
      ! a logical flag to say if these ii have been created
      logical :: created=.false.
      ! a logical flag to say if these ii have been formed
      logical :: formed=.false.
   end type particle_radmat_ii_type

   
end module radiation_materials_interpolation_data_types
