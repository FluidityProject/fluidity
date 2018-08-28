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
  use futils, only: int2str
  use particles, only : get_particles, get_particle_arrays, update_list_lengths
  use spud
  use parallel_tools, only: getnprocs, allsum
  use state_module
  use halos
  use elements, only: eval_shape
  use fields_base, only: ele_val, ele_loc
  use fields_calculations, only: dot_product
  use fields
  use pickers
  use detector_data_types
  use detector_tools, only: temp_insert, insert, allocate, deallocate
  use field_options
  
  implicit none

  private

  public :: calculate_diagnostics_from_particles, calculate_ratio_from_particles, calculate_numbers_from_particles

  contains

  subroutine calculate_diagnostics_from_particles(states, state_index, s_field)

    !!! Subroutine to determine which method is beind used (currently only ratio)
    use particles, only: particle_lists
    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field

    character(len= OPTION_PATH_LEN) :: lmethod

    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/name", lmethod, default = "ratio")

    select case(trim(lmethod))
       
     case("ratio")
        call calculate_ratio_from_particles(states, state_index, s_field)
       
     case("num_particles")
        call calculate_numbers_from_particles(s_field)
        
    end select

  end subroutine calculate_diagnostics_from_particles

  subroutine calculate_ratio_from_particles(states, state_index, s_field)

    !!Calculate s_field using the ratio method from particles
    !!First determine which particle groups/subgroups/attributes are being used
    !!Then determine the closest node for each particle and store attribute values
    !!Finally use the ratio method to calculate field values and place on
    !!Diagnostic scalar field

    use particles, only: particle_lists
    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field
    type(halo_type), pointer :: halo

    character(len=OPTION_PATH_LEN) :: lgroup, lattribute
    type(vector_field), pointer :: xfield
    type(detector_linked_list), allocatable, dimension(:,:) :: node_particles
    type(detector_linked_list), allocatable, dimension(:) :: p_array
    type(detector_type), pointer :: particle
    integer :: i, p_allocated, j
    real, allocatable, dimension(:) :: node_values
    real, allocatable, dimension(:) :: node_part_count ! real instead of integer, so we can use halo_accumulate
    integer :: element, node_number
    real, allocatable, dimension(:) :: local_crds
    integer, dimension(:), pointer :: nodes
    integer :: nprocs
    real :: att_value
    real :: ratio_val

    integer, allocatable, dimension(:) :: summed_particles, add_particles, remove_particles
    integer :: group_attribute
    integer, allocatable, dimension(:) :: group_arrays
    integer :: min_thresh
    integer :: max_thresh

!!!
    integer :: dum_counter

    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/group/name", lgroup)
    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/group/attribute/name", lattribute)
    ewrite(2,*) "Calculate diagnostic field from particle group: ", trim(lgroup), ", attribute: ", trim(lattribute)

    !Initialize field as 0
    s_field%val(:) = 0

    !Set minimum and maximum particle thresholds per control volume
    if (have_option("/particles/particle_group::"//trim(lgroup)//"/particle_spawning")) then
       call get_option("/particles/particle_group::"//trim(lgroup)//"/particle_spawning/min_cv_threshhold", min_thresh)
       call get_option("/particles/particle_group::"//trim(lgroup)//"/particle_spawning/max_cv_threshhold", max_thresh)
    end if
    !min_thresh = 20 !likely change to be set in schema file
    !max_thresh = 40 !!!

    xfield=>extract_vector_field(states(1), "Coordinate")

    !Call subroutine to check if particles are allocated

    ! contribution of local particles to non-owned nodes are summed
    ! into the owner in the halo_accumulate calls below
    ! we might be safe to assume we only need to add into halo 1 nodes (as these
    ! are the only ones that make up owned elements), but let's include halo 2 to be sure
    ! Only run this if nprocs > 1
    nprocs = getnprocs()
    if (nprocs>1) then
       halo => s_field%mesh%halos(2)
    end if

    if (allocated(particle_lists)) then
       !Particles are initialized, call subroutine to get relevant particle arrays and attributes
       call get_particle_arrays(lgroup, group_arrays, group_attribute, lattribute=lattribute)

       !Allocate arrays to store summed attribute values and particle counts at nodes
       allocate(node_particles(size(group_arrays),node_count(s_field)))
       allocate(node_values(node_count(s_field)))
       allocate(node_part_count(node_count(s_field)))
       allocate(summed_particles(size(group_arrays)))
       allocate(add_particles(size(group_arrays)))
       allocate(remove_particles(size(group_arrays)))
       node_values(:) = 0
       node_part_count(:) = 0
       summed_particles(:) = 0
       
       !Loop over particle arrays
       do i = 1,size(group_arrays)
          particle => particle_lists(group_arrays(i))%first
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
             call temp_insert(particle,node_particles(i,node_number))
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
       if (have_option("/particles/particle_group::"//trim(lgroup)//"/particle_spawning")) then
          do i = 1,node_count(s_field)
             !Count number of particles per node and ensure thresholds are not broken
             if (node_part_count(i)<min_thresh) then
                !make sure this is only called for nodes on this proc
                do j = 1,size(group_arrays)
                   if (node_part_count(i)>0.and.node_particles(j,i)%length>0) then
                      do while(node_part_count(i)<min_thresh)
                         call spawn_particles(node_part_count(i), node_values(i), node_particles(:,i), group_arrays, group_attribute, xfield, add_particles)
                         summed_particles=summed_particles+add_particles
                      end do
                   end if
                end do
             end if
             if (node_part_count(i)>max_thresh) then
                do while(node_part_count(i)>max_thresh)
                   call delete_particles(node_part_count(i), node_values(i), node_particles(:,i), group_arrays, group_attribute,remove_particles)
                   summed_particles=summed_particles-remove_particles
                end do
             end if

             !Determine field value from ratio method
             ratio_val = node_values(i)/node_part_count(i)
             !Store value on field (if node has at least one particle)
             if (node_part_count(i).eq.0) then
                s_field%val(i) = 0
             else
                s_field%val(i) = ratio_val
             end if
          end do

          do j=1,size(group_arrays)
             call allsum(summed_particles(j))
             particle_lists(group_arrays(j))%total_num_det=particle_lists(group_arrays(j))%total_num_det+summed_particles(j)
          end do
       else
          do i = 1,node_count(s_field)
             !Determine field value from ratio method
             ratio_val = node_values(i)/node_part_count(i)
             !Store value on field (if node has at least one particle)
             if (node_part_count(i).eq.0) then
                s_field%val(i) = 0
             else
                s_field%val(i) = ratio_val
             end if
          end do
       end if

       
       deallocate(node_values)
       deallocate(node_part_count)
       deallocate(node_particles)
       deallocate(add_particles)
       deallocate(remove_particles)
       deallocate(summed_particles)

       ! all values in owned nodes should now be correct
       ! now we need to make sure the halo is updated accordingly
       if (nprocs>1) then
          call halo_update(s_field)
       end if
    else
       !Particles not yet setup, Initial field value will be 0
    end if
    

  end subroutine calculate_ratio_from_particles

  subroutine calculate_numbers_from_particles(s_field)

    !!Calculate s_field using the ratio method from particles
    !!First determine which particle groups/subgroups/attributes are being used
    !!Then determine the closest node for each particle and store attribute values
    !!Finally use the ratio method to calculate field values and place on
    !!Diagnostic scalar field

    use particles, only: particle_lists
    type(scalar_field), intent(inout) :: s_field
    type(halo_type), pointer :: halo

    character(len=OPTION_PATH_LEN) :: lgroup
    type(detector_linked_list), allocatable, dimension(:) :: p_array
    type(detector_type), pointer :: particle
    integer :: i, p_allocated
    real, allocatable, dimension(:) :: node_part_count ! real instead of integer, so we can use halo_accumulate
    integer :: element, node_number
    real, allocatable, dimension(:) :: local_crds
    integer, dimension(:), pointer :: nodes
    integer :: nprocs

    integer :: group_attribute
    integer, allocatable, dimension(:) :: group_arrays


    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/group/name", lgroup)
    ewrite(2,*) "Calculate diagnostic field from particle group: ", trim(lgroup)

    !Initialize field as 0
    s_field%val(:) = 0

    ! contribution of local particles to non-owned nodes are summed
    ! into the owner in the halo_accumulate calls below
    ! we might be safe to assume we only need to add into halo 1 nodes (as these
    ! are the only ones that make up owned elements), but let's include halo 2 to be sure
    ! Only run this if nprocs > 1
    nprocs = getnprocs()
    if (nprocs>1) then
       halo => s_field%mesh%halos(2)
    end if

    if (allocated(particle_lists)) then
       !Particles are initialized, call subroutine to get relevant particle arrays and attributes
       call get_particle_arrays(lgroup, group_arrays, group_attribute)

       !Allocate arrays to store summed attribute values and particle counts at nodes
       allocate(node_part_count(node_count(s_field)))
       node_part_count(:) = 0

       !Loop over particle arrays
       do i = 1,size(group_arrays)
          particle => particle_lists(group_arrays(i))%first
          if (associated(particle)) then !Only work on arrays if local_particles exist on this processor
             allocate(local_crds(size(particle%local_coords)))
          end if
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
          if (allocated(local_crds)) then
             deallocate(local_crds)
          end if
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
    else
       !Particles not yet setup, Initial field value will be 0
    end if

  end subroutine calculate_numbers_from_particles

  subroutine spawn_particles(node_part_count, node_values, node_particles, group_arrays, group_attribute, xfield, add_particles)
    !Subroutine to calculate the number of particles in each control volume, and spawn additional
    !particles within a control volume if a minimum number of particles is not met

    use particles, only: particle_lists
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    real, intent(inout) :: node_values
    real, intent(inout) :: node_part_count
    integer, intent(in), dimension(:) :: group_arrays
    integer, intent(in) :: group_attribute
    type(vector_field), pointer, intent(in) :: xfield
    integer, dimension(:), intent(inout) :: add_particles
    
    type(detector_type), pointer :: particle
    type(detector_type), pointer :: temp_part

    real, dimension(:), allocatable :: node_coord, node_vector
    real :: pertubation

    character(len=OPTION_PATH_LEN) :: name
    character(len=OPTION_PATH_LEN) :: particle_name
    integer :: id, name_len, tot_len
    integer :: j
    integer, dimension(3) :: attribute_size

    add_particles(:) = 0
    
    do j = 1,size(group_arrays)
       temp_part => particle_lists(group_arrays(j))%last
       if (associated(temp_part)) then
          id = temp_part%id_number
          name_len = len(int2str(id+1))+1 !id length + '_'
          tot_len = len(trim(temp_part%name))
          name = trim(temp_part%name(1:tot_len-name_len))
          attribute_size(1) = size(temp_part%attributes)
          attribute_size(2) = size(temp_part%old_attributes)
          attribute_size(3) = size(temp_part%old_fields)
          allocate(node_coord(size(temp_part%local_coords)))
          allocate(node_vector(size(temp_part%local_coords)))
       end if

       particle => node_particles(j)%first
       do while(associated(particle))
          !duplicate particles, pertubating coords towards cv
          !get new lcoords from global coords
          !
          particle_name = trim(name)//int2str(id+1)
          temp_part => null()
          call allocate(temp_part, size(particle%position), size(particle%local_coords), attribute_size)

          temp_part%name = trim(particle_name)
          temp_part%id_number = id+1
          temp_part%list_id = particle%list_id

          !pertubate previous position towards cv by 1/100th of element size
          !global coords of node, node_coords-part_coords and normalize to get vector * 1/100th of mesh width and add to part_coords to get new_coords
          call set_local_node(particle%local_coords, node_coord)
          node_vector=node_coord-particle%local_coords
          pertubation = rand()/2 !currently will move particle randomly between it's current position and the closest node. Move set distance instead??
          !can also spawn randomly within control volume (from node position, move randomly in any direction up to distance of 1/2mesh size etc)
          temp_part%local_coords = particle%local_coords + node_vector*pertubation
          temp_part%element = particle%element
          call local_to_global(xfield, temp_part%local_coords, temp_part%element, temp_part%position)
          temp_part%type = LAGRANGIAN_DETECTOR
          temp_part%attributes(:) = 0
          temp_part%old_attributes(:) = 0
          temp_part%old_fields(:) = 0
          temp_part%attributes(group_attribute) = particle%attributes(group_attribute)

          node_values = node_values + temp_part%attributes(group_attribute)
          node_part_count = node_part_count + 1

          temp_part%previous => particle_lists(group_arrays(j))%last
          particle_lists(group_arrays(j))%last%next => temp_part
          particle_lists(group_arrays(j))%last => temp_part
          particle_lists(group_arrays(j))%last%next => null()
          particle_lists(group_arrays(j))%length = particle_lists(group_arrays(j))%length + 1
          add_particles(j) = add_particles(j) + 1
          id = id + 1

          particle => particle%temp_next

       end do
       if (allocated(node_coord)) then
          deallocate(node_coord)
          deallocate(node_vector)
       end if

    end do

  end subroutine spawn_particles

  subroutine set_local_node(local_coords, node_coord)

    real, dimension(:), intent(inout) :: node_coord
    real, dimension(:) :: local_coords

    node_coord(:) = 0
    select case(size(local_coords))
    case(2)
       node_coord(:) = 0
       node_coord(1) = 1.0
       if (local_coords(2)>local_coords(1)) then
          node_coord(:) = 0
          node_coord(2) = 1.0
       end if
    case(3)
       node_coord(:) = 0
       node_coord(1) = 1.0
       if (local_coords(2)>local_coords(1).or.local_coords(3)>local_coords(1)) then
          node_coord(:) = 0
          node_coord(2) = 1.0
          if (local_coords(3)>local_coords(2)) then
             node_coord(:) = 0
             node_coord(3) = 1.0
          end if
       end if
    case(4)
       node_coord(:) = 0
       node_coord(1) = 1.0
       if (local_coords(2)>local_coords(1).or.local_coords(3)>local_coords(1).or. &
            local_coords(4)>local_coords(1)) then
          node_coord(:) = 0
          node_coord(2) = 1.0
          if (local_coords(3)>local_coords(2).or.local_coords(4)>local_coords(2)) then
             node_coord(:) = 0
             node_coord(3) = 1.0
             if (local_coords(4)>local_coords(3)) then
                node_coord(:) = 0
                node_coord(4) = 1.0
             end if
          end if
       end if
    end select

  end subroutine set_local_node

  subroutine local_to_global(coordinates, l_coords, ele_number, global_coord)
    
    type(vector_field), pointer, intent(in) :: coordinates
    real, dimension(:), intent(in) :: l_coords
    integer, intent(in) :: ele_number
    real, dimension(:), intent(inout) :: global_coord

    integer :: i
    
    real, dimension(coordinates%dim,coordinates%dim+1) :: coord_ele
    real, dimension(coordinates%mesh%shape%loc) :: l_shape
    real, dimension(coordinates%mesh%shape%quadrature%vertices) :: dum_lcoords
    
    dum_lcoords = l_coords/sum(l_coords) ! You will already have this setup. Local coordinate array
    l_shape = eval_shape(coordinates%mesh%shape, dum_lcoords) ! Evaluate shape function at local_coordinates
    coord_ele=ele_val(coordinates, ele_number) ! global coordinate values of field at nodes of ele_number
    do i=1,size(global_coord)
       global_coord(i) = dot_product(l_shape,coord_ele(i,:)) ! should give global coordinates
    enddo


  end subroutine local_to_global

  subroutine delete_particles(node_part_count, node_values, node_particles, group_arrays, group_attribute, remove_particles)
    !Subroutine to calculate the number of particles in each control volume, and delete
    !particles within a control volume if a maximum number of particles is exceeded

    use particles, only: particle_lists
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    real, intent(inout) :: node_values
    real, intent(inout) :: node_part_count
    integer, intent(in), dimension(:) :: group_arrays
    integer, intent(in) :: group_attribute
    integer, dimension(:), intent(inout) :: remove_particles
    
    type(detector_type), pointer :: particle
    type(detector_type), pointer :: temp_part

    integer :: j
    real :: rand1, rand2, rand3

    remove_particles(:) = 0
    
    do j = 1,size(group_arrays)
       particle =>node_particles(j)%first
       do while(associated(particle))
          rand1=rand()
          if (rand1>0.5) then
             rand2=rand()
             if (rand2>0.5) then
                if (associated(particle%previous)) then
                   particle%previous%next => particle%next
                else
                   particle_lists(group_arrays(j))%first => particle%next
                end if
                if (associated(particle%next)) then
                   particle%next%previous => particle%previous
                else
                   particle_lists(group_arrays(j))%last => particle%previous
                end if
                temp_part =>particle%temp_next
                remove_particles(j)=remove_particles(j)+1
                particle_lists(group_arrays(j))%length = particle_lists(group_arrays(j))%length - 1
                node_values = node_values - particle%attributes(group_attribute)
                node_part_count = node_part_count - 1
                call deallocate(particle)
             else
                rand3=rand()
                if (rand3>0.5) then
                   if (associated(particle%previous)) then
                      particle%previous%next => particle%next
                   else
                      particle_lists(group_arrays(j))%first => particle%next
                   end if
                   if (associated(particle%next)) then
                      particle%next%previous => particle%previous
                   else
                      particle_lists(group_arrays(j))%last => particle%previous
                   end if
                   temp_part =>particle%temp_next
                   remove_particles(j)=remove_particles(j)+1
                   particle_lists(group_arrays(j))%length = particle_lists(group_arrays(j))%length - 1
                   node_values = node_values - particle%attributes(group_attribute)
                   node_part_count = node_part_count - 1
                   call deallocate(particle)
                end if
             end if
          else
             rand2=rand()
             if (rand2>0.5) then
                rand3=rand()
                if (rand3>0.5) then
                   if (associated(particle%previous)) then
                      particle%previous%next => particle%next
                   else
                      particle_lists(group_arrays(j))%first => particle%next
                   end if
                   if (associated(particle%next)) then
                      particle%next%previous => particle%previous
                   else
                      particle_lists(group_arrays(j))%last => particle%previous
                   end if
                   temp_part =>particle%temp_next
                   remove_particles(j)=remove_particles(j)+1
                   particle_lists(group_arrays(j))%length = particle_lists(group_arrays(j))%length - 1
                   node_values = node_values - particle%attributes(group_attribute)
                   node_part_count = node_part_count - 1
                   call deallocate(particle)
                end if
             end if
          end if
          if (associated(particle)) then
             particle => particle%temp_next
          else
             particle =>temp_part
             temp_part => null()
          end if
       end do
    end do
    
    
  end subroutine delete_particles

end module particle_diagnostics
