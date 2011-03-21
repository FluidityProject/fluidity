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
      module procedure energy_group_set_ii_allocate, &
                       region_id_ii_allocate_all, &
                       region_id_ii_allocate, &
                       dataset_ii_allocate, &
                       ii_size_allocate_from_options, &
                       ii_size_allocate
   end interface allocate
   
   interface zero
      module procedure energy_group_set_ii_zero, &
                       region_id_ii_zero, &
                       ii_size_zero, &
                       region_id_ii_size_zero
   end interface zero
   
   interface set
      module procedure set_ii_size, &
                       set_region_id_ii_size
   end interface set
   
   interface create
      module procedure particle_radmat_ii_create_from_options, &
                       ii_size_create_from_options
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
      
      ! local variables
      integer :: g_set
      integer :: number_of_energy_group_set
      
      ewrite(1,*) 'Create radiation particle_radmat_ii from options for ', trim(particle_radmat%name)
                                                                                 
      ! deduce the number of energy group sets
      number_of_energy_group_set = option_count(trim(particle_radmat%option_path)//'/energy_discretisation/energy_group_set')
            
      allocate(particle_radmat_ii%energy_group_set_ii(number_of_energy_group_set))
      
      ! create each energy_group_set_ii as needed
      energy_group_set_loop: do g_set = 1,number_of_energy_group_set
         
         ! create the size type
         call create(particle_radmat_ii%energy_group_set_ii(g_set)%ii_size, &
                     particle_radmat, &
                     state, &
                     g_set)
      
         call allocate(particle_radmat_ii%energy_group_set_ii(g_set))
      
         call zero(particle_radmat_ii%energy_group_set_ii(g_set))
      
      end do energy_group_set_loop
      
      particle_radmat_ii%created = .true.
      
      ewrite(1,*) 'Finished creating radiation particle_radmat_ii from options for ', trim(particle_radmat%name)
      
   end subroutine particle_radmat_ii_create_from_options

   ! --------------------------------------------------------------------------

   subroutine ii_size_create_from_options(ii_size, &
                                          particle_radmat, &
                                          state, &
                                          g_set)
      
      !!< Create the ii_size from the options (including the material mesh size)

      type(ii_size_type), intent(inout) :: ii_size
      type(particle_radmat_type), intent(in) :: particle_radmat
      type(state_type), intent(in) :: state
      integer, intent(in) :: g_set

      ewrite(1,*) 'Create radiation ii_size from options for ', trim(particle_radmat%name),' for g_set ',g_set
                 
      call allocate(ii_size, &
                    particle_radmat, &
                    state, &
                    g_set)
                      
      call zero(ii_size)
            
      call set(ii_size, &
               particle_radmat, &
               state, &
               g_set)
                                   
      ii_size%size_set = .true.
     
      ewrite(1,*) 'Finished creating radiation ii_size from options for ', trim(particle_radmat%name),' for g_set ',g_set
      
   end subroutine ii_size_create_from_options

   ! --------------------------------------------------------------------------

   subroutine ii_size_allocate_from_options(ii_size, &
                                            particle_radmat, &
                                            state, &
                                            g_set)
      
      !!< Allocate the ii_size from options

      type(ii_size_type), intent(inout) :: ii_size
      type(particle_radmat_type), intent(in) :: particle_radmat
      type(state_type), intent(in) :: state
      integer, intent(in) :: g_set

      ! local variables
      integer :: fn_space_dof
      logical :: region_id_size_allocate
      type(mesh_type), pointer :: material_fn_space
      character(len=OPTION_PATH_LEN) :: region_id_path
      character(len=OPTION_PATH_LEN) :: material_fn_space_name
       
      ewrite(1,*) 'Allocate radiation particle_radmat_ii_size from options for ', trim(particle_radmat%name),' for g_set ',g_set
      
      ! get the material fn space name
      call get_option(trim(particle_radmat%option_path)//'/energy_discretisation/energy_group_set['//int2str(g_set-1)//']/angular_discretisation/method/parity/angular_moment_set[0]/mesh/name',material_fn_space_name)
      
      ! extract the material fn_space of this energy group set of this particle type 
      material_fn_space => extract_mesh(state, trim(material_fn_space_name))
            
      ! determine the number of material fn_space dof 
      fn_space_dof = node_count(material_fn_space)

      region_id_path = trim(particle_radmat%option_path)//'/region_id_material_mapping'
          
      link_region_id: if (have_option(trim(region_id_path))) then
      
         region_id_size_allocate = .true.
         
      else link_region_id   
      
         region_id_size_allocate = .false.
     
      end if link_region_id

      call allocate(ii_size, &
                    fn_space_dof, &
                    region_id_size_allocate)
       
   end subroutine ii_size_allocate_from_options

   ! --------------------------------------------------------------------------

   subroutine ii_size_allocate(ii_size, &
                               fn_space_dof, &
                               region_id_size_allocate)
      
      !!< Allocate the ii_size
         
      type(ii_size_type), intent(inout) :: ii_size
      integer, intent(in) :: fn_space_dof
      logical, intent(in) :: region_id_size_allocate
      
      ! local variables
      integer :: status
      
      ewrite(1,*) 'Allocate radiation ii_size'
      
      ! allocate the component sizes as needed
      region_id_size_allocate_if: if (region_id_size_allocate) then
      
         allocate(ii_size%region_id_ii_size(fn_space_dof),STAT=status)
              
         if (status /= 0) FLAbort('Issue allocating ii_size%region_id_ii_size in particle_radmat_ii_size_allocate')
      
      end if region_id_size_allocate_if      
   
   end subroutine ii_size_allocate

   ! --------------------------------------------------------------------------

   subroutine ii_size_zero(ii_size)
      
      !!< Zero the relevant components of the ii_size type

      type(ii_size_type), intent(inout) :: ii_size
      
      ewrite(1,*) 'Zero radiation ii_size'
      
      zero_region_id_ii_size: if (allocated(ii_size%region_id_ii_size)) then
      
         call zero(ii_size%region_id_ii_size)
      
      end if zero_region_id_ii_size
   
   end subroutine ii_size_zero

   ! --------------------------------------------------------------------------

   subroutine region_id_ii_size_zero(region_id_ii_size)
      
      !!< Zero the relevant components of the region_id_ii_size type
      
      type(region_id_ii_size_type), dimension(:), allocatable, intent(inout) :: region_id_ii_size
      
      ! local variables
      integer :: n
      
      allocated_if: if (allocated(region_id_ii_size)) then
      
         node_loop: do n = 1,size(region_id_ii_size)
            
            region_id_ii_size(n)%fraction_size = 0
            
         end do node_loop
      
      end if allocated_if
         
   end subroutine region_id_ii_size_zero

   ! --------------------------------------------------------------------------

   subroutine energy_group_set_ii_allocate(energy_group_set_ii)
      
      !!< Allocate the energy_group_set_ii from the energy_group_set_ii%ii_size component

      type(energy_group_set_ii_type), intent(inout) :: energy_group_set_ii
      
      ! local variables
      integer :: status
      
      ewrite(1,*) 'Allocate radiation energy_group_set_ii'
      
      ! if the size data is not set then abort 
      no_size: if (.not. energy_group_set_ii%ii_size%size_set) then
      
         FLAbort('Cannot allocate a energy_group_set_ii type if the size type energy_group_set_ii%ii_size is not already set')
      
      end if no_size

      ! region id mapping allocate
      region_id_allocate: if (allocated(energy_group_set_ii%ii_size%region_id_ii_size)) then
                  
         allocate(energy_group_set_ii%region_id_ii(size(energy_group_set_ii%ii_size%region_id_ii_size)),STAT=status)
         
         if (status /= 0) FLAbort('Issue allocating memory for energy_group_set_ii%region_id_ii in energy_group_set_ii_allocate')

         call allocate(energy_group_set_ii%region_id_ii, &
                       energy_group_set_ii%ii_size%region_id_ii_size)

      end if region_id_allocate
            
   end subroutine energy_group_set_ii_allocate

   ! --------------------------------------------------------------------------

   subroutine region_id_ii_allocate_all(region_id_ii, &
                                        region_id_ii_size)
   
      !!< Allocate the region_id_ii using the region_id_ii_size
      
      type(region_id_ii_type), dimension(:), intent(inout) :: region_id_ii
      type(region_id_ii_size_type), dimension(:), intent(in) :: region_id_ii_size
      
      ! local variables
      integer :: n
                          
      node_loop: do n = 1,size(region_id_ii)
         
         call allocate(region_id_ii(n), &
                       region_id_ii_size(n))
                     
      end do node_loop
            
   end subroutine region_id_ii_allocate_all

   ! --------------------------------------------------------------------------

   subroutine region_id_ii_allocate(region_id_ii, &
                                    region_id_ii_size)
   
      !!< Allocate the region_id_ii using the region_id_ii_size 
      
      type(region_id_ii_type), intent(inout) :: region_id_ii
      type(region_id_ii_size_type), intent(in) :: region_id_ii_size
      
      ! local variables
      integer :: fraction_size
                          
      fraction_size = region_id_ii_size%fraction_size
      
      call allocate(region_id_ii%dataset_ii, &
                    fraction_size)
                                                      
   end subroutine region_id_ii_allocate

   ! --------------------------------------------------------------------------

   subroutine dataset_ii_allocate(dataset_ii, &
                                  fraction_size)
   
      !!< Allocate the dataset_ii type
      
      type(dataset_ii_type), intent(inout) :: dataset_ii
      integer, intent(in) :: fraction_size
      
      ! local variables
      integer :: status

      allocate(dataset_ii%physical_radmat_ii%fraction(fraction_size),STAT=status)

      if (status /= 0) FLAbort('Issue allocating memory for dataset_ii%physical_radmat_ii%fraction in dataset_ii_allocate')
                     
      allocate(dataset_ii%physical_radmat_ii%radmat_base_coordinate(fraction_size),STAT=status)

      if (status /= 0) FLAbort('Issue allocating memory for dataset_ii%physical_radmat_ii%radmat_base_coordinate in dataset_ii_allocate')
                                 
   end subroutine dataset_ii_allocate

   ! --------------------------------------------------------------------------
   
   subroutine energy_group_set_ii_zero(energy_group_set_ii)
      
      !!< Zero the energy_group_set_ii data type
   
      type(energy_group_set_ii_type), intent(inout) :: energy_group_set_ii
      
      ewrite(1,*) 'Zero radiation energy_group_set_ii'
      
      zero_region_id_ii: if (allocated(energy_group_set_ii%region_id_ii)) then
         
         call zero(energy_group_set_ii%region_id_ii)
      
      end if zero_region_id_ii
   
   end subroutine energy_group_set_ii_zero
   
   ! --------------------------------------------------------------------------
   
   subroutine region_id_ii_zero(region_id_ii)
      
      !!< Zero the region_id_ii data type for all vele
   
      type(region_id_ii_type), dimension(:), allocatable, intent(inout) :: region_id_ii
      
      ! local variable
      integer :: vele
      
      allocated_if: if (allocated(region_id_ii)) then
      
         vele_loop: do vele = 1,size(region_id_ii) 
            
            region_id_ii(vele)%dataset_ii%dataset_radmat_number = 0
             
            region_id_ii(vele)%dataset_ii%physical_radmat_ii%physical_radmat_number = 0

            region_id_ii(vele)%dataset_ii%physical_radmat_ii%radmat_base_coordinate = 0
            
            region_id_ii(vele)%dataset_ii%physical_radmat_ii%fraction = 0.0
         
         end do vele_loop
      
      end if allocated_if
         
   end subroutine region_id_ii_zero
   
   ! --------------------------------------------------------------------------

   subroutine set_ii_size(ii_size, &
                          particle_radmat, &
                          state, &
                          g_set)
      
      !!< Set the components of the ii_size type associated with a energy group set for a particle type

      type(ii_size_type), intent(inout) :: ii_size
      type(particle_radmat_type), intent(in) :: particle_radmat
      type(state_type), intent(in) :: state
      integer, intent(in) :: g_set
      
      ewrite(1,*) 'Set radiation ii_size'
      
      ! set the region id size components
      set_region_id_size: if (allocated(ii_size%region_id_ii_size)) then
      
         call set(ii_size%region_id_ii_size, &
                  particle_radmat, &
                  state, &
                  g_set)
      
      end if set_region_id_size
      
      ii_size%size_set = .true.

   end subroutine set_ii_size
   
   ! --------------------------------------------------------------------------

   subroutine set_region_id_ii_size(region_id_ii_size, &
                                    particle_radmat, &
                                    state, &
                                    g_set)
      
      !!< Set the region id vele ii size associated with a energy group set for a particle type 
      !!< This assumes that the material fn space is discontinuous which should have already 
      !!< been asserted in options checking
      
      type(region_id_ii_size_type), dimension(:), allocatable, intent(inout) :: region_id_ii_size
      type(particle_radmat_type), intent(in) :: particle_radmat
      type(state_type), intent(in) :: state
      integer, intent(in) :: g_set
      
      ! local variables
      integer :: iloc
      integer :: id
      integer :: dmat
      integer :: pmat
      integer :: vele
      integer :: imap
      integer :: region_id_mapping_shape(2)
      integer :: fraction_size_set
      integer, dimension(:), allocatable :: region_id_mapping
      integer, dimension(:), pointer :: element_nodes   
      character(len=OPTION_PATH_LEN) :: region_id_mapping_path 
      character(len=OPTION_PATH_LEN) :: data_set_name 
      character(len=OPTION_PATH_LEN) :: physical_material_name 
      character(len=OPTION_PATH_LEN) :: material_fn_space_name
      type(mesh_type), pointer :: material_fn_space
      type(scalar_field), pointer :: particle_flux
         
      ! first check that region_id_ii_size is allocated
      check_allocated: if (.not. allocated(region_id_ii_size)) then
         
         FLAbort('Cannot set region_id_ii_size if it is not already allocated')
         
      end if check_allocated

      ! extract the first group within this group set to find the volume elements region id    
      call extract_flux_group_from_group_set(state, &
                                             trim(particle_radmat%option_path), &
                                             g = 1, &  
                                             g_set = g_set, & 
                                             particle_flux = particle_flux)
             
      ! get the material fn space name
      call get_option(trim(particle_radmat%option_path)//'/energy_discretisation/energy_group_set['//int2str(g_set-1)//']/angular_discretisation/method/parity/angular_moment_set[0]/mesh/name',material_fn_space_name)
      
      ! extract the material fn_space of this energy group set of this particle type 
      material_fn_space => extract_mesh(state, trim(material_fn_space_name))
       
      region_id_mapping_path = trim(particle_radmat%option_path)//'/region_id_material_mapping/region_to_physical_radiation_material_map'
                        
      ! form each volume element ii
      vele_loop: do vele = 1,size(region_id_ii_size)
 
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
                  call get_option(trim(region_id_mapping_path)//'['//int2str(imap-1)//']/data_set/name',data_set_name)
                  call get_option(trim(region_id_mapping_path)//'['//int2str(imap-1)//']/physical_material/name',physical_material_name)
                      
                  ! find the matching read in data set and physical material via names comparisons
                  dmat_loop: do dmat = 1,size(particle_radmat%dataset_radmats)
                     
                     dmat_match: if (trim(particle_radmat%dataset_radmats(dmat)%name) == trim(data_set_name)) then
                            
                        pmat_loop: do pmat = 1,size(particle_radmat%dataset_radmats(dmat)%physical_radmats)
                               
                           pmat_match: if (trim(particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%name) == trim(physical_material_name)) then
                                 
                              ! set the fraction size via counting the number of interpolation dimensions
                              fraction_size_set = &
                              option_count(trim(particle_radmat%dataset_radmats(dmat)%physical_radmats(pmat)%option_path//'/interpolation_dimension'))
                              
                              ! set the size for each local discontinuous node of the material fn space
                                 
                              element_nodes => ele_nodes(material_fn_space, vele)
                              
                              local_node_loop: do iloc = 1, size(element_nodes)
                                                               
                                 region_id_ii_size(element_nodes(iloc))%fraction_size = fraction_size_set
                              
                              end do local_node_loop
                                                                    
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
             
         no_region_id_map_found: if (region_id_ii_size(vele)%fraction_size == 0) then
                
            ewrite(-1,*) 'Error for vele',vele
            FLAbort('Error in set_region_id_ii_size as no region id mapping found for the vele')
             
         end if no_region_id_map_found
   
      end do vele_loop
         
   end subroutine set_region_id_ii_size
   
   ! --------------------------------------------------------------------------

end module radiation_materials_interpolation_create
