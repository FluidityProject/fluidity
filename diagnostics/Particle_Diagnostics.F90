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
  use state_module
  use fields
  use fields_base
  use detector_data_types
  use field_options
  
  implicit none

  private

  public :: calculate_diagnostics_from_particles, calculate_particles_from_ratio

  contains

  subroutine calculate_diagnostics_from_particles(states, state_index, s_field)

    !!! in here need to determine which method is being used (currently only ratio)
    !!! then need to determine the particle group, subgroup and attribute used for calculation
    !!! interpolate these attribute values onto nodes for field

    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field

    character(len= OPTION_PATH_LEN) :: lmethod
    character(len= OPTION_PATH_LEN) :: name
    integer :: i

    !ewrite(2,*) "index: ", state_index, "size: ", size(states)
    !ewrite(2,*) "halos", size(states(state_index)%halos) !!!??? Need to import halos?
    !do i = 1,size(states(state_index)%halos)
    !   name = trim(states(state_index)%halo_names(i))
    !   ewrite(2,*) "halo name: ", name
    !end do
    
    !states(state_index)

    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/name", lmethod, default = "ratio")

    select case(trim(lmethod))
       
      case("ratio")
        call calculate_particles_from_ratio(states, state_index, s_field)

!!$      case default
!!$        if(present(stat)) then
!!$          stat = 1
!!$          return
!!$        end if
!!$        FLAbort("Diagnostic algorithm method " // trim(lmethod) // " not found")
    end select

  end subroutine calculate_diagnostics_from_particles

  subroutine calculate_particles_from_ratio(states, state_index, s_field)

    !!Calculate s_field using the ratio method from particles
    !!First determine which particle groups/subgroups/attributes are being used

    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field

    character(len=OPTION_PATH_LEN) :: lgroup, lattribute
    type(detector_linked_list), allocatable, dimension(:) :: p_array
    type(detector_type), pointer :: particle
    integer :: attribute_number, i, p_allocated
    integer, allocatable, dimension(:) :: field_nodes, field_nodes_count
    integer :: element, node_number, dim
    real, allocatable, dimension(:) :: local_crds
    integer, dimension(:), pointer :: nodes
    real :: att_value
    real :: ratio_val

    integer :: group_attribute
    integer, allocatable, dimension(:) :: group_arrays

    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/group/name", lgroup)
    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/group/attribute/name", lattribute)
    ewrite(2,*) "Calculate diagnostic field from particle group: ", trim(lgroup), ", attribute: ", trim(lattribute)
    call get_particles(p_array, p_allocated)

    if (p_allocated.eq.0) then
       !Particles not yet setup, get initial field values from options tree
       
    else
       call get_option("/geometry/dimension",dim)
       call get_particle_arrays(lgroup, lattribute, group_arrays, group_attribute)

       allocate(field_nodes(node_count(s_field)))
       allocate(field_nodes_count(node_count(s_field)))
       field_nodes(:) = 0
       field_nodes_count(:) = 0
       !initialise field as 0
       s_field%val(:) = 0
       !ewrite(2,*) "147258"
       !ewrite(2,*) "Node count", node_count(s_field)!!!
       !ewrite(2,*) "Element count", element_count(s_field)

      ! do i = 1,element_count(s_field)
       !   ewrite
      ! end do

       !ewrite(2,*) "258147"
       do i = 1,size(group_arrays)
          particle => p_array(group_arrays(i))%first
          if (associated(particle)) then
             allocate(local_crds(size(particle%local_coords)))
          end if
          do while(associated(particle))
             element = particle%element
             local_crds = particle%local_coords
             nodes => ele_nodes(s_field, element)
             att_value = particle%attributes(group_attribute)

             !ewrite(2,*) "element number", element

             call get_particle_node(nodes, local_crds, dim, node_number)
             field_nodes(node_number) = field_nodes(node_number) + att_value
             field_nodes_count(node_number) = field_nodes_count(node_number) + 1
             particle => particle%next
             
          end do
          if (allocated(local_crds)) then
             deallocate(local_crds)
          end if
       end do
       do i = 1,node_count(s_field)
          !ewrite(2,*) i, field_nodes_count(i)

          ratio_val = field_nodes(i)*1.0/field_nodes_count(i)
          if (field_nodes_count(i).eq.0) then
             s_field%val(i) = 0
          else
             s_field%val(i) = ratio_val
          end if
       end do
       deallocate(field_nodes)
       deallocate(field_nodes_count)
    end if

  end subroutine calculate_particles_from_ratio

  subroutine get_particle_node(nodes, local_crds, dim, node_number)

    !!Given local_nodes, local coords and dimension
    !!returns node number of closest node

    integer, dimension(:), pointer, intent(in) :: nodes
    real, dimension(:), intent(in) :: local_crds
    integer, intent(in) :: dim
    integer, intent(out) :: node_number

    real :: maxval

    select case(dim)
       
    case(1)
       maxval = local_crds(1)
       if (local_crds(2)>maxval) then
          node_number = nodes(2)
       else
          node_number = nodes(1)
       end if    
    case(2)
       maxval = local_crds(1)
       node_number = nodes(1)
       if (local_crds(2)>maxval) then
          maxval = local_crds(2)
          node_number = nodes(2)
          if (local_crds(3)>maxval) then
             node_number = nodes(3)
          end if
       end if
    case(3)
       maxval = local_crds(1)
       node_number = nodes(1)
       if (local_crds(2)>maxval) then
          maxval = local_crds(2)
          node_number = nodes(2)
          if (local_crds(3)>maxval) then
             maxval = local_crds(3)
             node_number = nodes(3)
             if (local_crds(4)>maxval) then
                node_number = nodes(4)
             end if
          end if
       end if
    end select

  end subroutine get_particle_node

end module particle_diagnostics
