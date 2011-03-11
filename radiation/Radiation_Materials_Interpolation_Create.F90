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

module radiation_materials_interpolation_create

   use state_module  
   use fields
   use spud
   use global_parameters, only : OPTION_PATH_LEN
      
   use radiation_materials_interpolation_data_types
   use radiation_materials
   use radiation_extract_flux_field

   implicit none
   
   private 

   public :: allocate, &
             zero, &
             set, &
             create
   
   interface allocate
      module procedure particle_radmat_ii_allocate, &
                       region_id_vele_ii_allocate_all, &
                       region_id_vele_ii_allocate, &
                       dataset_vele_ii_allocate, &
                       particle_radmat_ii_size_allocate_from_options, &
                       particle_radmat_ii_size_allocate
   end interface allocate
   
   interface zero
      module procedure particle_radmat_ii_zero, &
                       region_id_vele_ii_zero, &
                       particle_radmat_ii_size_zero, &
                       region_id_vele_ii_size_zero
   end interface zero
   
   interface set
      module procedure set_particle_radmat_ii_size, &
                       set_region_id_vele_ii_size
   end interface set
   
   interface create
      module procedure particle_radmat_ii_create_from_options, &
                       particle_radmat_ii_size_create_from_options
   end interface create

contains

   ! --------------------------------------------------------------------------

   subroutine particle_radmat_ii_create_from_options(particle_radmat_ii, &
                                                     particle_radmat, &
                                                     state) 
   
      !!< Create from the options the particle_radmat_ii data type for a 
      !!< particular particle determined by the particle_radmat
      
      type(particle_radmat_ii_type), intent(inout) :: particle_radmat_ii
      type(particle_radmat_type), intent(in) :: particle_radmat      
      type(state_type), intent(in) :: state
      
      ewrite(1,*) 'Create radiation particle_radmat_ii from options for ', trim(particle_radmat%name)
      
      call create(particle_radmat_ii%particle_radmat_ii_size, &
                  particle_radmat, &
                  state)
      
      call allocate(particle_radmat_ii)
      
      call zero(particle_radmat_ii)
      
      particle_radmat_ii%created = .true.
      
      ewrite(1,*) 'Finished creating radiation particle_radmat_ii from options for ', trim(particle_radmat%name)
      
   end subroutine particle_radmat_ii_create_from_options

   ! --------------------------------------------------------------------------

   subroutine particle_radmat_ii_size_create_from_options(particle_radmat_ii_size, &
                                                          particle_radmat, &
                                                          state)
      
      !!< Create the particle_radmat_ii_size from the options (including the mesh size)

      type(particle_radmat_ii_size_type), intent(inout) :: particle_radmat_ii_size
      type(particle_radmat_type), intent(in) :: particle_radmat
      type(state_type), intent(in) :: state

      ewrite(1,*) 'Create radiation particle_radmat_ii_size from options for ', trim(particle_radmat%name)
                 
      call allocate(particle_radmat_ii_size, &
                    particle_radmat, &
                    state)
                      
      call zero(particle_radmat_ii_size)
            
      call set(particle_radmat_ii_size, &
               particle_radmat, &
               state)
                                   
      particle_radmat_ii_size%size_set = .true.
     
      ewrite(1,*) 'Finished creating radiation particle_radmat_ii_size from options for ', trim(particle_radmat%name)
      
   end subroutine particle_radmat_ii_size_create_from_options

   ! --------------------------------------------------------------------------

   subroutine particle_radmat_ii_size_allocate_from_options(particle_radmat_ii_size, &
                                                            particle_radmat, &
                                                            state)
      
      !!< Allocate the particle_radmat_ii_size from options

      type(particle_radmat_ii_size_type), intent(inout) :: particle_radmat_ii_size
      type(particle_radmat_type), intent(in) :: particle_radmat
      type(state_type), intent(in) :: state

      ! local variables
      integer :: vele_np
      logical :: region_id_size_allocate
      type(scalar_field), pointer :: particle_flux
      character(len=OPTION_PATH_LEN) :: region_id_path
       
      ewrite(1,*) 'Allocate radiation particle_radmat_ii_size from options'
      
      ! extract the g 1 np flux field to find the number of volume elements    
      call extract_flux_group_g(state, &
                                trim(particle_radmat%name), &
                                g = 1, &  
                                particle_flux = particle_flux)
      
      ! determine the number of volume elements of the radiation problem
      vele_np = ele_count(particle_flux)

      region_id_path = trim(particle_radmat%option_path)//'/region_id_material_mapping'
          
      link_region_id: if (have_option(trim(region_id_path))) then
      
         region_id_size_allocate = .true.
         
      else link_region_id   
      
         region_id_size_allocate = .false.
     
      end if link_region_id

      call allocate(particle_radmat_ii_size, &
                    vele_np, &
                    region_id_size_allocate)
       
   end subroutine particle_radmat_ii_size_allocate_from_options

   ! --------------------------------------------------------------------------

   subroutine particle_radmat_ii_size_allocate(particle_radmat_ii_size, &
                                               vele_np, &
                                               region_id_size_allocate)
      
      !!< Allocate the particle_radmat_ii_size
         
      type(particle_radmat_ii_size_type), intent(inout) :: particle_radmat_ii_size
      integer, intent(in) :: vele_np
      logical, intent(in) :: region_id_size_allocate
      
      ! local variables
      integer :: status
      
      ewrite(1,*) 'Allocate radiation particle_radmat_ii_size'
      
      ! allocate the component sizes as needed
      region_id_size_allocate_if: if (region_id_size_allocate) then
      
         allocate(particle_radmat_ii_size%region_id_vele_ii_size(vele_np),STAT=status)
              
         if (status /= 0) FLAbort('Issue allocating particle_radmat_ii_size%region_id_vele_ii_size in particle_radmat_ii_size_allocate')
      
      end if region_id_size_allocate_if      
   
   end subroutine particle_radmat_ii_size_allocate

   ! --------------------------------------------------------------------------

   subroutine particle_radmat_ii_size_zero(particle_radmat_ii_size)
      
      !!< Zero the relevant components of the particle_radmat_ii_size type

      type(particle_radmat_ii_size_type), intent(inout) :: particle_radmat_ii_size
      
      ewrite(1,*) 'Zero radiation particle_radmat_ii_size'
      
      zero_region_id_vele_ii_size: if (allocated(particle_radmat_ii_size%region_id_vele_ii_size)) then
      
         call zero(particle_radmat_ii_size%region_id_vele_ii_size)
      
      end if zero_region_id_vele_ii_size
   
   end subroutine particle_radmat_ii_size_zero

   ! --------------------------------------------------------------------------

   subroutine region_id_vele_ii_size_zero(region_id_vele_ii_size)
      
      !!< Zero the relevant components of the region_id_vele_ii_size type
      
      type(region_id_vele_ii_size_type), dimension(:), allocatable, intent(inout) :: region_id_vele_ii_size
      
      ! local variables
      integer :: vele
      
      allocated_if: if (allocated(region_id_vele_ii_size)) then
      
         vele_loop: do vele = 1,size(region_id_vele_ii_size)
            
            region_id_vele_ii_size(vele)%fraction_size = 0
            
         end do vele_loop
      
      end if allocated_if
         
   end subroutine region_id_vele_ii_size_zero

   ! --------------------------------------------------------------------------

   subroutine particle_radmat_ii_allocate(particle_radmat_ii)
      
      !!< Allocate the particle_radmat_ii from the particle_radmat_ii%particle_radmat_ii_size component

      type(particle_radmat_ii_type), intent(inout) :: particle_radmat_ii
      
      ! local variables
      integer :: status
      
      ewrite(1,*) 'Allocate radiation particle_radmat_ii'
      
      ! if the size data is not set then abort 
      no_size: if (.not. particle_radmat_ii%particle_radmat_ii_size%size_set) then
      
         FLAbort('Cannot allocate a particle_radmat_ii type if the size type particle_radmat_ii%particle_radmat_ii_size is not already set')
      
      end if no_size

      ! region id mapping allocate
      region_id_allocate: if (allocated(particle_radmat_ii%particle_radmat_ii_size%region_id_vele_ii_size)) then
                  
         allocate(particle_radmat_ii%region_id_vele_ii(size(particle_radmat_ii%particle_radmat_ii_size%region_id_vele_ii_size)),STAT=status)
         
         if (status /= 0) FLAbort('Issue allocating memory for particle_radmat_ii%region_id_vele_ii in particle_radmat_ii_allocate')

         call allocate(particle_radmat_ii%region_id_vele_ii, &
                       particle_radmat_ii%particle_radmat_ii_size%region_id_vele_ii_size)

      end if region_id_allocate
            
   end subroutine particle_radmat_ii_allocate

   ! --------------------------------------------------------------------------

   subroutine region_id_vele_ii_allocate_all(region_id_vele_ii, &
                                             region_id_vele_ii_size)
   
      !!< Allocate the region_id_vele_ii using the region_id_vele_ii_size for all vele 
      
      type(region_id_vele_ii_type), dimension(:), intent(inout) :: region_id_vele_ii
      type(region_id_vele_ii_size_type), dimension(:), intent(in) :: region_id_vele_ii_size
      
      ! local variables
      integer :: vele
                          
      ! allocate the arrays associated with each volume element vele
      vele_loop: do vele = 1,size(region_id_vele_ii)
         
         call allocate(region_id_vele_ii(vele), &
                       region_id_vele_ii_size(vele))
                     
      end do vele_loop
            
   end subroutine region_id_vele_ii_allocate_all

   ! --------------------------------------------------------------------------

   subroutine region_id_vele_ii_allocate(region_id_vele_ii, &
                                         region_id_vele_ii_size)
   
      !!< Allocate the region_id_vele_ii using the region_id_vele_ii_size for one vele 
      
      type(region_id_vele_ii_type), intent(inout) :: region_id_vele_ii
      type(region_id_vele_ii_size_type), intent(in) :: region_id_vele_ii_size
      
      ! local variables
      integer :: fraction_size
                          
      fraction_size = region_id_vele_ii_size%fraction_size
      
      call allocate(region_id_vele_ii%dataset_vele_ii, &
                    fraction_size)
                                                      
   end subroutine region_id_vele_ii_allocate

   ! --------------------------------------------------------------------------

   subroutine dataset_vele_ii_allocate(dataset_vele_ii, &
                                       fraction_size)
   
      !!< Allocate the dataset_vele_ii type
      
      type(dataset_vele_ii_type), intent(inout) :: dataset_vele_ii
      integer, intent(in) :: fraction_size
      
      ! local variables
      integer :: status

      allocate(dataset_vele_ii%physical_radmat_vele_ii%fraction(fraction_size),STAT=status)

      if (status /= 0) FLAbort('Issue allocating memory for dataset_vele_ii%physical_radmat_vele_ii%fraction in dataset_vele_ii_allocate')
                     
      allocate(dataset_vele_ii%physical_radmat_vele_ii%radmat_base_coordinate(fraction_size),STAT=status)

      if (status /= 0) FLAbort('Issue allocating memory for dataset_vele_ii%physical_radmat_vele_ii%radmat_base_coordinate in dataset_vele_ii_allocate')
                                 
   end subroutine dataset_vele_ii_allocate

   ! --------------------------------------------------------------------------
   
   subroutine particle_radmat_ii_zero(particle_radmat_ii)
      
      !!< Zero the particle_radmat_ii data type
   
      type(particle_radmat_ii_type), intent(inout) :: particle_radmat_ii
      
      ewrite(1,*) 'Zero radiation particle_radmat_ii'
      
      zero_region_id_ii: if (allocated(particle_radmat_ii%region_id_vele_ii)) then
         
         call zero(particle_radmat_ii%region_id_vele_ii)
      
      end if zero_region_id_ii
   
   end subroutine particle_radmat_ii_zero
   
   ! --------------------------------------------------------------------------
   
   subroutine region_id_vele_ii_zero(region_id_vele_ii)
      
      !!< Zero the region_id_vele_ii data type for all vele
   
      type(region_id_vele_ii_type), dimension(:), allocatable, intent(inout) :: region_id_vele_ii
      
      ! local variable
      integer :: vele
      
      allocated_if: if (allocated(region_id_vele_ii)) then
      
         vele_loop: do vele = 1,size(region_id_vele_ii) 
            
            region_id_vele_ii(vele)%dataset_vele_ii%dataset_radmat_number = 0
             
            region_id_vele_ii(vele)%dataset_vele_ii%physical_radmat_vele_ii%physical_radmat_number = 0

            region_id_vele_ii(vele)%dataset_vele_ii%physical_radmat_vele_ii%radmat_base_coordinate = 0
            
            region_id_vele_ii(vele)%dataset_vele_ii%physical_radmat_vele_ii%fraction = 0.0
         
         end do vele_loop
      
      end if allocated_if
         
   end subroutine region_id_vele_ii_zero
   
   ! --------------------------------------------------------------------------

   subroutine set_particle_radmat_ii_size(particle_radmat_ii_size, &
                                          particle_radmat, &
                                          state)
      
      !!< Set the components of the particle_radmat_ii_size type

      type(particle_radmat_ii_size_type), intent(inout) :: particle_radmat_ii_size
      type(particle_radmat_type), intent(in) :: particle_radmat
      type(state_type), intent(in) :: state
      
      ewrite(1,*) 'Set radiation particle_radmat_ii_size'
      
      ! set the region id size components
      set_region_id_size: if (allocated(particle_radmat_ii_size%region_id_vele_ii_size)) then
      
         call set(particle_radmat_ii_size%region_id_vele_ii_size, &
                  particle_radmat, &
                  state)
      
      end if set_region_id_size
      
      particle_radmat_ii_size%size_set = .true.

   end subroutine set_particle_radmat_ii_size
   
   ! --------------------------------------------------------------------------

   subroutine set_region_id_vele_ii_size(region_id_vele_ii_size, &
                                         particle_radmat, &
                                         state)
      
      !!< Set the region id vele ii size
      
      type(region_id_vele_ii_size_type), dimension(:), allocatable, intent(inout) :: region_id_vele_ii_size
      type(particle_radmat_type), intent(in) :: particle_radmat
      type(state_type), intent(in) :: state
      
      ! local variables
      integer :: id
      integer :: dmat
      integer :: pmat
      integer :: vele
      integer :: imap
      integer :: region_id_mapping_shape(2)
      integer, dimension(:), allocatable :: region_id_mapping
      character(len=OPTION_PATH_LEN) :: region_id_mapping_path 
      character(len=OPTION_PATH_LEN) :: data_set_filename 
      character(len=OPTION_PATH_LEN) :: physical_material_name 
      type(scalar_field), pointer :: particle_flux
      
      ! first check that region_id_vele_ii_size is allocated
      check_allocated: if (.not. allocated(region_id_vele_ii_size)) then
         
         FLAbort('Cannot set region_id_vele_ii_size if it is not already allocated')
         
      end if check_allocated
      
      ! extract the g 1 np flux field to find the volume elements region id    
      call extract_flux_group_g(state, &
                                trim(particle_radmat%name), &
                                g = 1, &  
                                particle_flux = particle_flux)
       
      region_id_mapping_path = trim(particle_radmat%option_path)//'/region_id_material_mapping/region_to_physical_radiation_material_map'
      
      ! form each volume element ii
      vele_loop: do vele = 1,size(region_id_vele_ii_size)
 
         ! check each region to physical radiation material mapping within the options             
         region_mapping_loop: do imap = 1,option_count(trim(region_id_mapping_path)) 
                
            ! get the length of the region ids list
            region_id_mapping_shape = option_shape(trim(region_id_mapping_path)//'['//int2str(imap-1)//']/region_id')
                
            allocate(region_id_mapping(region_id_mapping_shape(1)))
                
            ! get the actual region ids associated with this mapping
            call get_option(trim(region_id_mapping_path)//'['//int2str(imap-1)//']/region_id',region_id_mapping)
                
            ! check each region id against the region id of this vele
            id_loop: do id = 1,size(region_id_mapping)
                   
               id_match: if (ele_region_id(particle_flux,vele) == region_id_mapping(id)) then
                      
                  ! get the data set and physical material name for this region id mapping
                  call get_option(trim(region_id_mapping_path)//'['//int2str(imap-1)//']/radiation_data_set_name/file_name',data_set_filename)
                  call get_option(trim(region_id_mapping_path)//'['//int2str(imap-1)//']/radiation_physical_material_name/name',physical_material_name)
                      
                  ! find the matching read in data set and physical material via names comparisons
                  dmat_loop: do dmat = 1,size(particle_radmat%dataset_radmats)
                     
                     dmat_match: if (trim(particle_radmat%dataset_radmats(dmat)%file_name) == trim(data_set_filename)) then
                            
                        pmat_loop: do pmat = 1,size(particle_radmat%dataset_radmats(dmat)%physical_radmats)
                               
                           pmat_match: if (trim(particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%name) == trim(physical_material_name)) then
                                 
                              ! set the fraction size via counting the number of interpolation dimensions
                              region_id_vele_ii_size(vele)%fraction_size = &
                              option_count(trim(particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%option_path//'/interpolation_dimension'))
                                                                    
                              ! found a fraction size for this vele for region id mapping so exit loop
                              exit region_mapping_loop
                                
                           end if pmat_match
                              
                        end do pmat_loop
                        
                    end if dmat_match
                       
                 end do dmat_loop
                           
               end if id_match
                   
            end do id_loop
                
            deallocate(region_id_mapping)
                
         end do region_mapping_loop
             
         if (allocated(region_id_mapping)) deallocate(region_id_mapping)
             
         no_region_id_map_found: if (region_id_vele_ii_size(vele)%fraction_size == 0) then
                
            ewrite(-1,*) 'Error for vele',vele
            FLAbort('Error in set_region_id_vele_ii_size as no region id mapping found for the vele')
             
         end if no_region_id_map_found
   
      end do vele_loop
         
   end subroutine set_region_id_vele_ii_size
   
   ! --------------------------------------------------------------------------

end module radiation_materials_interpolation_create
