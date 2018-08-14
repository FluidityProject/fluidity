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
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!    Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA

#include "fdebug.h"

module particle_diagnostics

  use fldebug
  use global_parameters, only : OPTION_PATH_LEN
  use particles, only : get_particles, get_particle_arrays
  use spud
  use parallel_tools, only: getnprocs
  use state_module
  use halos
  use fields
  use detector_data_types
  use field_options
  
  implicit none

  private

  public :: calculate_diagnostics_from_particles, calculate_ratio_from_particles, calculate_numbers_from_particles

  contains

  subroutine calculate_diagnostics_from_particles(states, state_index, s_field)

    !!! Subroutine to determine which method is beind used (currently only ratio)

    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field

    character(len= OPTION_PATH_LEN) :: lmethod
    character(len= OPTION_PATH_LEN) :: name
    integer :: i

    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/name", lmethod, default = "ratio")

    select case(trim(lmethod))
       
     case("ratio")
        call calculate_ratio_from_particles(states, state_index, s_field)
       
     case("num_particles")
        call calculate_numbers_from_particles(states, state_index, s_field)
        
    end select

  end subroutine calculate_diagnostics_from_particles

  subroutine calculate_ratio_from_particles(states, state_index, s_field)

    !!Calculate s_field using the ratio method from particles
    !!First determine which particle groups/subgroups/attributes are being used
    !!Then determine the closest node for each particle and store attribute values
    !!Finally use the ratio method to calculate field values and place on
    !!Diagnostic scalar field

    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field
    type(halo_type), pointer :: halo

    character(len=OPTION_PATH_LEN) :: lgroup, lattribute
    type(detector_linked_list), allocatable, dimension(:) :: p_array
    type(detector_type), pointer :: particle
    integer :: attribute_number, i, p_allocated
    real, allocatable, dimension(:) :: node_values
    real, allocatable, dimension(:) :: node_part_count ! real instead of integer, so we can use halo_accumulate
    integer :: element, node_number, dim
    real, allocatable, dimension(:) :: local_crds
    integer, dimension(:), pointer :: nodes
    integer :: nprocs
    real :: att_value
    real :: ratio_val

    integer :: group_attribute
    integer, allocatable, dimension(:) :: group_arrays
    integer :: min_thresh
    integer :: max_thresh

    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/group/name", lgroup)
    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/group/attribute/name", lattribute)
    ewrite(2,*) "Calculate diagnostic field from particle group: ", trim(lgroup), ", attribute: ", trim(lattribute)

    !Initialize field as 0
    s_field%val(:) = 0

    !Set minimum and maximum particle thresholds per control volume
    min_thresh = 20 !likely change to be set in schema file
    max_thresh = 40 !!!

    !Call subroutine to check if particles are allocated
    call get_particles(p_array, p_allocated)

    ! contribution of local particles to non-owned nodes are summed
    ! into the owner in the halo_accumulate calls below
    ! we might be safe to assume we only need to add into halo 1 nodes (as these
    ! are the only ones that make up owned elements), but let's include halo 2 to be sure
    ! Only run this if nprocs > 1
    nprocs = getnprocs()
    if (nprocs>1) then
       halo => s_field%mesh%halos(2)
    end if

    if (p_allocated.eq.0) then
       !Particles not yet setup, Initial field value will be 0
       
    else
       !Particles are initialized, call subroutine to get relevant particle arrays and attributes
       call get_particle_arrays(lgroup, group_arrays, group_attribute, lattribute=lattribute)

       !Allocate arrays to store summed attribute values and particle counts at nodes
       allocate(node_values(node_count(s_field)))
       allocate(node_part_count(node_count(s_field)))
       node_values(:) = 0
       node_part_count(:) = 0

       !Loop over particle arrays
       do i = 1,size(group_arrays)
          particle => p_array(group_arrays(i))%first
          if (associated(particle)) then !Only work on arrays if local_particles exist on this processor
             allocate(local_crds(size(particle%local_coords)))
          end if
          do while(associated(particle))
             !Get element, local_crds and attribute value of each particle
             element = particle%element
             local_crds = particle%local_coords
             att_value = particle%attributes(group_attribute)

             !Find nodes for specified element
             nodes => ele_nodes(s_field, element)

             ! work out nearest node based on max. local coordinate
             node_number = nodes(maxloc(local_crds, dim=1))
             
             !Store particle attribute value on closest node
             node_values(node_number) = node_values(node_number) + att_value
             !Increase particle count for this node by 1
             node_part_count(node_number) = node_part_count(node_number) + 1.0
             particle => particle%next
             
          end do
          if (allocated(local_crds)) then
             deallocate(local_crds)
          end if
       end do

       if (nprocs>1) then
          call halo_accumulate(halo, node_part_count)
          call halo_accumulate(halo, node_values)
       end if

       !Loop over all nodes
       do i = 1,node_count(s_field)
!!$          !Count number of particles per node and ensure thresholds are not broken
!!$          if (node_part_count(i)<min_thresh) then
!!$             !make sure this is only called for nodes on this proc
!!$             call spawn_particles(node_part_count(i), node_values(i), i, p_array, group_arrays, group_attribute, min_thresh)
!!$          end if
!!$          if (node_part_count(i)>max_thresh) then
!!$             call delete_particles(node_part_count(i), node_values(i), p_array, group_arrays, group_attribute, max_thresh)
!!$          end if
          
          !Determine field value from ratio method
          ratio_val = node_values(i)/node_part_count(i)
          !Store value on field (if node has at least one particle)
          if (node_part_count(i).eq.0) then
             s_field%val(i) = 0
          else
             s_field%val(i) = ratio_val
          end if
       end do
       deallocate(node_values)
       deallocate(node_part_count)

       ! all values in owned nodes should now be correct
       ! now we need to make sure the halo is updated accordingly
       if (nprocs>1) then
          call halo_update(s_field)
       end if
    end if

  end subroutine calculate_ratio_from_particles

  subroutine calculate_numbers_from_particles(states, state_index, s_field)

    !!Calculate s_field using the ratio method from particles
    !!First determine which particle groups/subgroups/attributes are being used
    !!Then determine the closest node for each particle and store attribute values
    !!Finally use the ratio method to calculate field values and place on
    !!Diagnostic scalar field

    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field
    type(halo_type), pointer :: halo

    character(len=OPTION_PATH_LEN) :: lgroup
    type(detector_linked_list), allocatable, dimension(:) :: p_array
    type(detector_type), pointer :: particle
    integer :: i, p_allocated
    real, allocatable, dimension(:) :: node_part_count ! real instead of integer, so we can use halo_accumulate
    integer :: element, node_number, dim
    real, allocatable, dimension(:) :: local_crds
    integer, dimension(:), pointer :: nodes
    integer :: nprocs

    integer :: group_attribute
    integer, allocatable, dimension(:) :: group_arrays


    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/group/name", lgroup)
    ewrite(2,*) "Calculate diagnostic field from particle group: ", trim(lgroup)

    !Initialize field as 0
    s_field%val(:) = 0

    !Call subroutine to check if particles are allocated
    call get_particles(p_array, p_allocated)

    ! contribution of local particles to non-owned nodes are summed
    ! into the owner in the halo_accumulate calls below
    ! we might be safe to assume we only need to add into halo 1 nodes (as these
    ! are the only ones that make up owned elements), but let's include halo 2 to be sure
    ! Only run this if nprocs > 1
    nprocs = getnprocs()
    if (nprocs>1) then
       halo => s_field%mesh%halos(2)
    end if

    if (p_allocated.eq.0) then
       !Particles not yet setup, Initial field value will be 0
       
    else
       !Particles are initialized, call subroutine to get relevant particle arrays and attributes
       call get_particle_arrays(lgroup, group_arrays, group_attribute)

       !Allocate arrays to store summed attribute values and particle counts at nodes
       allocate(node_part_count(node_count(s_field)))
       node_part_count(:) = 0

       !Loop over particle arrays
       do i = 1,size(group_arrays)
          particle => p_array(group_arrays(i))%first
          do while(associated(particle))
             !Get element, local_crds of each particle
             element = particle%element
             local_crds = particle%local_coords
             !Find nodes for specified element
             nodes => ele_nodes(s_field, element)
             ! work out nearest node based on max. local coordinate
             node_number = nodes(maxloc(local_crds, dim=1))
             !Increase particle count for this node by 1
             node_part_count(node_number) = node_part_count(node_number) + 1.0
             particle => particle%next
             
          end do
       end do

       if (nprocs>1) then
          call halo_accumulate(halo, node_part_count)
       end if

       !Loop over all nodes
       do i = 1,node_count(s_field)
          !Store value on field 
          s_field%val(i) = node_part_count(i)
       end do
       deallocate(node_part_count)

       ! all values in owned nodes should now be correct
       ! now we need to make sure the halo is updated accordingly
       if (nprocs>1) then
          call halo_update(s_field)
       end if
    end if

  end subroutine calculate_numbers_from_particles

!!$  subroutine spawn_particles(node_part_count, node_values, i, p_array, group_arrays, group_attribute, min_thresh)
!!$    !Subroutine to calculate the number of particles in each control volume, and spawn additional
!!$    !particles within a control volume if a minimum number of particles is not met
!!$
!!$    type(detector_linked_list), intent(inout), dimension(:) :: p_array
!!$    real, intent(inout) :: node_values
!!$    real, intent(inout) :: node_part_count
!!$    integer, intent(in), dimension(:) :: group_arrays
!!$    integer, intent(in) :: group_attribute
!!$    integer, intent(in) :: min_thresh
!!$    integer, intent(in) :: i !Node number
!!$    
!!$    type(detector_type), pointer :: particle
!!$    type(detector_type), pointer :: temp_part
!!$
!!$    real :: ratio_val
!!$    real :: rand_num
!!$
!!$    real :: val
!!$    character(len=OPTION_PATH_LEN) :: name
!!$    integer :: name_len
!!$    character(len=OPTION_PATH_LEN) :: dum_name
!!$    character(len=OPTION_PATH_LEN) :: particle_name
!!$    integer :: id
!!$    integer :: k
!!$
!!$    ratio_val = node_values/node_part_count
!!$
!!$    do k = 1, min_thresh-node_part_count
!!$       !Need to calculate random location within this control volume determine which element it is in
!!$       !then determine local coords for location and finally global coords/position
!!$       !Also need to calculate random attribute value based off current ratio value (weighted by ratio)
!!$       
!!$       rand_num = (rand()+ratio_val)/2.0
!!$       if (rand_num>=0.5) then
!!$          val = 1
!!$       else
!!$          val = 0
!!$       end if
!!$       
!!$       temp_part => p_array(group_arrays(val+1))%last
!!$       id = temp_part%id_number
!!$       name = trim(temp_part%name)
!!$       
!!$       name_len = len_trim(name)
!!$       write(dum_name, fmt) id
!!$       name_cut = 1+len_trim(dum_name)
!!$       
!!$       write(particle_name, fmt) name(1:name_len-name_cut)//"_", id+1
!!$       !!make new linked list of particles within each cv? (
!!$       allocate(particle)
!!$       allocate(particle%position(size(temp_part%position)))
!!$       allocate(particle%local_coords(size(temp_part%local_coords)))
!!$       particle%name=trim(particle_name)
!!$       particle%id_number = id+1
!!$       particle%position = !(from element and lcoords)
!!$       particle%element = !(specific element)
!!$       particle%local_coords = !(random-closest to node)
!!$       particle%type = LAGRANGIAN_DETECTOR
!!$
!!$       !get elements node is contained in from function  node_neigh_mesh(mesh, node_number) result (node_neigh)
!!$       !or function node_neigh_scalar(field, node_number) result (node_neigh) (in Fields_Base.F90)
!!$       
!!$       !p_array(group_arrays())%last%next => particle
!!$       !p_array(group_arrays())%last => particle
!!$       !p_array(group_arrays())%last%next => null()
!!$       !p_array(group_arrays())%length = p_array(group_arrays())%length+1
!!$       !p_array(group_arrays())%detector_names(id+1) = particle%name
!!$
!!$       !do all this in new routine?!?
!!$       
!!$       allocate(particle%attributes(size(temp_part%attributes)))
!!$       allocate(particle%old_attributes(size(temp_part%old_attributes)))
!!$       allocate(particle%old_fields(size(temp_part%old_fields)))
!!$       
!!$       !initialise as 0? or call calculate routines for individual particle
!!$       particle%attributes(:) = 0
!!$       particle%old_attributes(:) = 0
!!$       particle%old_fields(:) = 0
!!$       particle%attribute(group_attribute) = val
!!$       
!!$       !update node_part_count and node_values with new particle (update node_part_count out of loop)
!!$       node_values = node_values + val
!!$       deallocate(particle) !!?? will this kill the previous particle?
!!$    end do
!!$    
!!$
!!$    
!!$
!!$  end subroutine spawn_particles
!!$
!!$  subroutine delete_particles(node_part_count, node_values, p_array, group_arrays, group_attribute, max_thresh)
!!$    !Subroutine to calculate the number of particles in each control volume, and delete
!!$    !particles within a control volume if a maximum number of particles is exceeded
!!$    type(detector_linked_list), intent(inout), dimension(:) :: p_array
!!$    real, intent(inout) :: node_values
!!$    real, intent(inout) :: node_part_count
!!$    integer, intent(in), dimension(:) :: group_arrays
!!$    integer, intent(in) :: group_attribute
!!$    integer, intent(in) :: max_thresh
!!$    
!!$    type(detector_type), pointer :: particle
!!$    type(detector_type), pointer :: temp_part
!!$
!!$  end subroutine delete_particles
!!$
end module particle_diagnostics
