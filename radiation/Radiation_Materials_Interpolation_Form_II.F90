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

module radiation_materials_interpolation_form_ii
   
   use futils
   use global_parameters, only : OPTION_PATH_LEN
   use spud
   use state_module  
   use fields
   
   use radiation_materials_data_types
   use radiation_materials_interpolation_data_types
   use radiation_materials_interpolation_create, only : zero
   use radiation_extract_flux_field

   implicit none
   
   private 

   public :: form
   
   interface form
      module procedure form_particle_radmat_ii, &
                       form_energy_group_set_ii, &
                       form_region_id_ii_all, &
                       form_dataset_ii, &
                       form_physical_radmat_ii, &
                       form_fraction_and_radmat_base_coordinate
   end interface form
   
contains

   ! --------------------------------------------------------------------------

   subroutine form_particle_radmat_ii(particle_radmat_ii, &
                                      particle_radmat, &
                                      state) 
   
      !!< Form the material data interpolation/mixing instructions for each DoF of the material mesh
      !!< of each group set for this particle type

      type(particle_radmat_ii_type), intent(inout) :: particle_radmat_ii
      type(particle_radmat_type), intent(in) :: particle_radmat
      type(state_type), intent(in) :: state
            
      ! local variables
      integer :: g_set
     
      ! first check that the particle_radmat_ii data type has been created
      check_created: if (.not. particle_radmat_ii%created) then
         
         FLAbort('Cannot form particle_radmat_ii if it is not created')
         
      end if check_created
      
      energy_group_set_loop: do g_set = 1,size(particle_radmat_ii%energy_group_set_ii)
         
         call form(particle_radmat_ii%energy_group_set_ii(g_set), &
                   particle_radmat, &
                   state, &
                   g_set)
      
      end do energy_group_set_loop
               
      ! set the flag for formed
      particle_radmat_ii%formed = .true.
      
   end subroutine form_particle_radmat_ii

   ! --------------------------------------------------------------------------

   subroutine form_energy_group_set_ii(energy_group_set_ii, &
                                       particle_radmat, &
                                       state, &
                                       g_set)
      
      type(energy_group_set_ii_type), intent(inout) :: energy_group_set_ii
      type(particle_radmat_type), intent(in) :: particle_radmat
      type(state_type), intent(in) :: state
      integer, intent(in) :: g_set

      ! local variables
      character(len=OPTION_PATH_LEN) :: region_id_path
      
      ! region id ii form
      region_id_ii_form: if (allocated(energy_group_set_ii%region_id_ii)) then

         region_id_path = trim(particle_radmat%option_path)//'/region_id_material_mapping'
         
         ! zero the region_id_ii components
         call zero(energy_group_set_ii%region_id_ii)
            
         ! form the values for region id ii
         call form(particle_radmat, &
                   state, &
                   energy_group_set_ii%region_id_ii, &
                   g_set)
                                  
      end if region_id_ii_form
      
   end subroutine form_energy_group_set_ii

   ! --------------------------------------------------------------------------

   subroutine form_region_id_ii_all(particle_radmat, &
                                    state, &
                                    region_id_ii, &
                                    g_set)
      
      !!< Form region_id_vele_ii values for all the material mesh dof
      
      type(particle_radmat_type), intent(in) :: particle_radmat
      type(state_type), intent(in) :: state
      type(region_id_ii_type), dimension(:), allocatable, intent(inout) :: region_id_ii
      integer, intent(in) :: g_set
      
      ! local variables
      logical :: found_mapping
      integer :: id
      integer :: vele
      integer :: imap
      integer :: iloc
      integer :: region_id_mapping_shape(2)
      integer, dimension(:), allocatable :: region_id_mapping
      integer, dimension(:), pointer :: element_nodes      
      character(len=OPTION_PATH_LEN) :: region_id_mapping_path 
      character(len=OPTION_PATH_LEN) :: data_set_filename 
      character(len=OPTION_PATH_LEN) :: physical_material_name 
      character(len=OPTION_PATH_LEN) :: material_fn_space_name
      type(mesh_type), pointer :: material_fn_space
      type(scalar_field), pointer :: particle_flux
      
      ! first check that region_id_ii is allocated
      check_allocated: if (.not. allocated(region_id_ii)) then
         
         FLAbort('Cannot form region_id_ii if it is not already allocated')
         
      end if check_allocated
            
      ! extract the first group within this group set to find the volume elements region id    
      call extract_flux_group_from_group_set(state, &
                                             trim(particle_radmat%option_path), &
                                             g = 1, &  
                                             g_set = g_set, & 
                                             particle_flux = particle_flux)
             
      region_id_mapping_path = trim(particle_radmat%option_path)//'/region_id_material_mapping/region_to_physical_radiation_material_map'

      ! get the material fn space name
      call get_option(trim(particle_radmat%option_path)//'/energy_group_set['//int2str(g_set-1)//']/angular_discretisation/mesh/name',material_fn_space_name)
      
      ! extract the material fn_space of this energy group set of this particle type 
      material_fn_space => extract_mesh(state, trim(material_fn_space_name))
      
      ! form each volume element ii
      vele_loop: do vele = 1,size(region_id_ii)
             
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
                  call get_option(trim(region_id_mapping_path)//'['//int2str(imap-1)//&
                                 &']/radiation_data_set_name/file_name',data_set_filename)

                  call get_option(trim(region_id_mapping_path)//'['//int2str(imap-1)//&
                                 &']/radiation_physical_material_name/name',physical_material_name)
                  
                  ! form the ii for each discontinuous dof
                  element_nodes => ele_nodes(material_fn_space, vele)
                              
                  local_node_loop: do iloc = 1, size(element_nodes)
                  
                     call form(region_id_ii(element_nodes(iloc))%dataset_ii, &
                               particle_radmat, &
                               trim(data_set_filename), &
                               trim(physical_material_name), &
                               found_mapping, &
                               element_nodes(iloc), &
                               state)
                  
                  end do local_node_loop
                  
                  ! if found a mapping for this vele exit the region mapping loop
                  if (found_mapping) exit region_mapping_loop
                                                 
               end if id_match
                   
            end do id_loop
                
            deallocate(region_id_mapping)
                
         end do region_mapping_loop
             
         if (allocated(region_id_mapping)) deallocate(region_id_mapping)
                
      end do vele_loop
      
   end subroutine form_region_id_ii_all

   ! --------------------------------------------------------------------------

   subroutine form_dataset_ii(dataset_ii, &
                              particle_radmat, &
                              data_set_filename, &
                              physical_material_name, &
                              found_mapping, &
                              inode, &
                              state)
   
      !!< Form the dataset_ii type associated with this dof 
      
      type(dataset_ii_type), intent(inout) :: dataset_ii
      type(particle_radmat_type), intent(in) :: particle_radmat
      character(len=*), intent(in) :: data_set_filename
      character(len=*), intent(in) :: physical_material_name
      logical, intent(out) :: found_mapping
      integer, intent(in) :: inode
      type(state_type), intent(in) :: state
      
      ! local variables
      integer :: dmat
      
      found_mapping = .false.
      
      ! find the matching read in data set and physical material via names comparisons
      dmat_loop: do dmat = 1,size(particle_radmat%dataset_radmats)
         
         dmat_match: if (trim(particle_radmat%dataset_radmats(dmat)%file_name) == trim(data_set_filename)) then

            dataset_ii%dataset_radmat_number = dmat
            
            call form(dataset_ii%physical_radmat_ii, &
                      particle_radmat%dataset_radmats(dmat), &
                      physical_material_name, &
                      found_mapping, &
                      inode, &
                      state)
            
            ! if found a mapping exit the dmat loop
            if (found_mapping) exit dmat_loop
            
        end if dmat_match
           
     end do dmat_loop
      
   end subroutine form_dataset_ii

   ! --------------------------------------------------------------------------

   subroutine form_physical_radmat_ii(physical_radmat_ii, &
                                      dataset_radmat, &
                                      physical_material_name, &
                                      found_mapping, &
                                      inode, &
                                      state)
   
      !!< Form the physical_radmat_vele_ii type associated with this dof
      
      type(physical_radmat_ii_type), intent(inout) :: physical_radmat_ii
      type(dataset_radmat_type), intent(in) :: dataset_radmat
      character(len=*), intent(in) :: physical_material_name
      logical, intent(inout) :: found_mapping
      integer, intent(in) :: inode
      type(state_type), intent(in) :: state

      ! local variables
      integer :: pmat
      
      pmat_loop: do pmat = 1,size(dataset_radmat%physical_radmats)
       
         pmat_match: if (trim(dataset_radmat%physical_radmats(pmat)%name) == trim(physical_material_name)) then
         
            physical_radmat_ii%physical_radmat_number = pmat
            
            call form(physical_radmat_ii%fraction, &
                      physical_radmat_ii%radmat_base_coordinate, &
                      dataset_radmat%physical_radmats(pmat), &
                      inode, &
                      state)
                    
            ! found a mapping so set flag and exit
            found_mapping = .true.
                  
            exit pmat_loop
        
         end if pmat_match
      
      end do pmat_loop
      
   end subroutine form_physical_radmat_ii

   ! --------------------------------------------------------------------------

   subroutine form_fraction_and_radmat_base_coordinate(fraction, &
                                                       radmat_base_coordinate, &
                                                       physical_radmat, &
                                                       inode, &
                                                       state)
   
      !!< Form the ii fraction and radmat_base_coordinate associated with this dof
      
      real, dimension(:), intent(inout) :: fraction
      integer, dimension(:), intent(inout) :: radmat_base_coordinate
      type(physical_radmat_type), intent(in) :: physical_radmat
      integer, intent(in) :: inode
      type(state_type), intent(in) :: state

      ! local variables
      integer :: status
      integer :: v
      integer :: f
      integer :: interpolation_values_shape(2)
      real, dimension(:), allocatable :: interpolation_values
      real :: actual_value
      character(len=OPTION_PATH_LEN) :: interpolation_dimension_path      
      character(len=OPTION_PATH_LEN) :: scalar_field_name      
      type(scalar_field), pointer :: sfield     
      
      ! form the fraction and base coordinate for each interpolation dimension
      fraction_loop: do f = 1,size(fraction)
         
         interpolation_dimension_path = trim(physical_radmat%option_path)//'/interpolation_dimension['//int2str(f-1)//']'
          
         attribute_none_if: if (have_option(trim(interpolation_dimension_path)//'/interpolation_attribute_none')) then
         
            fraction(f) = 0.0

            radmat_base_coordinate(f) = 1
         
         else attribute_none_if
            
            ! get the actual values
            attribute_if: if (have_option(trim(interpolation_dimension_path)//'/interpolation_attribute_prescribed')) then 
         
               call get_option(trim(interpolation_dimension_path)//'/interpolation_attribute_prescribed',actual_value)
            
            else if (have_option(trim(interpolation_dimension_path)//'/interpolation_attribute_scalar_field')) then attribute_if
            
              ! assume the scalar field is element wise (ie. only one value per element)
              
              call get_option(trim(interpolation_dimension_path)//'/interpolation_attribute_scalar_field/name',scalar_field_name)
              
              sfield => extract_scalar_field(state, &
                                             trim(scalar_field_name), &
                                             stat=status)
                            
              actual_value = node_val(sfield, inode)
              
            end if attribute_if
            
            ! find the radmat_base_coordinate and fraction from inspecting the interpolation values
            interpolation_values_shape = option_shape(trim(interpolation_dimension_path)//'/interpolation_values')

            allocate(interpolation_values(interpolation_values_shape(1)))

            interpolation_values = 0.0

            call get_option(trim(interpolation_dimension_path)//'/interpolation_values',interpolation_values)
                        
            ! check each interpolation value with the actual value bar the last to bound it 
            ! - if more than one interpolation values
            more_than_one: if (size(interpolation_values) > 1) then

               ! set the base coordinate initially as 1 to bound it 
               radmat_base_coordinate(f) = 1
            
               value_loop: do v = 1,size(interpolation_values) - 1
   
                  greater_if: if (actual_value > interpolation_values(v)) then
   
                     radmat_base_coordinate(f) = v
   
                  else greater_if
   
                     exit value_loop
   
                  end if greater_if
   
               end do value_loop
                        
               ! form the ii fraction knowing the base coordinate
               fraction(f) = max(0.0,min(1.0,(actual_value - interpolation_values(radmat_base_coordinate(f)) / &
               interpolation_values(radmat_base_coordinate(f)+1) - interpolation_values(radmat_base_coordinate(f)))))
            
            else more_than_one
               
               ! the size of interpolation values == 1 case

               radmat_base_coordinate(f) = 1
               
               fraction(f) = 0.0
            
            end if more_than_one

            deallocate(interpolation_values)
                     
         end if attribute_none_if
       
      end do fraction_loop
      
   end subroutine form_fraction_and_radmat_base_coordinate

   ! --------------------------------------------------------------------------

end module radiation_materials_interpolation_form_ii
