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
  use particles, only : get_particles, get_particle_arrays, update_list_lengths, &
       & set_particle_attributes, initialise_constant_particle_attributes
  use spud
  use parallel_tools, only: getnprocs, allsum
  use state_module
  use halos
  use elements, only: eval_shape
  use fields_base, only: ele_val, ele_loc
  use fields_calculations, only: dot_product
  use parallel_fields, only: node_owned
  use fields
  use pickers
  use detector_data_types
  use detector_tools, only: temp_insert, insert, allocate, deallocate
  use field_options
  use multimaterial_module, only: calculate_sum_material_volume_fractions
  use diagnostic_fields_new, only: calculate_dependencies, calculate_diagnostic_variable 
  
  implicit none

  private

  public :: initialise_particle_diagnostics, update_particle_diagnostics, &
       & calculate_diagnostics_from_particles, calculate_ratio_from_particles, &
       & calculate_numbers_from_particles, initialise_constant_particle_diagnostics, &
       & initialise_particle_diagnostic_fields, initialise_particle_material_fields

  contains

  subroutine initialise_constant_particle_diagnostics(state)
    !subroutine to initialise constant particle attributes and
    !multimaterial fields if 'from_particles'

    use particles, only: particle_lists
    type(state_type), dimension(:), intent(inout) :: state

    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path
    type(scalar_field), pointer :: s_field
    integer :: i, k
    integer :: particle_groups, list_counter, particle_materials
    integer, dimension(:), allocatable :: particle_arrays

    type(detector_type), pointer :: particle

    !Check if there are particles
    particle_groups = option_count("/particles/particle_group")

    if (particle_groups==0) return

    !Set up particle_lists
    allocate(particle_arrays(particle_groups))
    particle_arrays(:) = 0
    do i = 1,particle_groups
       particle_arrays(i) = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
    end do

    !Initialise constant particle attributes
    list_counter=1
    do i = 1,particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1,particle_arrays(i)
          particle => particle_lists(list_counter)%first
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          if (option_count(trim(subgroup_path) // "/attributes/attribute/constant").gt.0) then
             call initialise_constant_particle_attributes(state, subgroup_path, particle_lists(list_counter))
          end if
          list_counter = list_counter + 1
       end do
    end do

    !Check if MVF field is generated from particles
    particle_materials = option_count("material_phase/scalar_field::MaterialVolumeFraction/diagnostic/algorithm::from_particles")
    if (particle_materials.gt.0) then
       !Initialise MultialVolumeFraction fields dependent on particles
       do i = 1,size(state)
          s_field => extract_scalar_field(state(i), "MaterialVolumeFraction")     
          if (have_option(trim(s_field%option_path)//"/diagnostic/algorithm::from_particles")) then
             call calculate_diagnostics_from_particles(state, i, s_field)
          end if
       end do
       k = size(state)
       
       !Initialise internal MaterialVolumeFraction field 
       s_field => extract_scalar_field(state(k), "MaterialVolumeFraction")
       call calculate_sum_material_volume_fractions(state, s_field)
       call scale(s_field, -1.0)
       call addto(s_field, 1.0)
    end if
    
  end subroutine initialise_constant_particle_diagnostics

  subroutine initialise_particle_diagnostics(state)
    !subroutine to initialise particle attributes, diagnostic fields
    !dependent on particles, and diagnostic fields to be set after
    !particles are initialised
    use particles, only: particle_lists
    type(state_type), dimension(:), intent(inout) :: state
    type(state_type), dimension(size(state)) :: calculated_states

    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path, name
    type(vector_field), pointer :: xfield
    type(scalar_field), pointer :: s_field
    real :: current_time
    integer :: i, k
    integer :: dim, particle_groups, list_counter
    integer, dimension(:), allocatable :: particle_arrays

    type(detector_type), pointer :: particle
    
    !Check if there are particles
    particle_groups = option_count("/particles/particle_group")

    if (particle_groups==0) return

    !Set up particle_lists
    allocate(particle_arrays(particle_groups))
    particle_arrays(:) = 0
    do i = 1,particle_groups
       particle_arrays(i) = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
    end do

    !Allocate parameters
    xfield=>extract_vector_field(state(1), "Coordinate")
    call get_option("/geometry/dimension",dim)
    call get_option("/timestepping/current_time", current_time)

    !Initialise particle attributes
    list_counter=1
    do i = 1,particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1,particle_arrays(i)
          particle => particle_lists(list_counter)%first
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          if (option_count(trim(subgroup_path) // "/attributes/attribute/python").gt.0 .or. &
              & option_count(trim(subgroup_path) // "/attributes/attribute/python_fields").gt.0) then
             call set_particle_attributes(state, dim, current_time, subgroup_path, particle_lists(list_counter))
          end if
          list_counter = list_counter + 1
       end do
    end do

    !Initialise diagnostic fields generated from particles
    do i = 1,size(state)
       do k = 1,scalar_field_count(state(i))
          s_field => extract_scalar_field(state(i),k)     
          if (have_option(trim(s_field%option_path)//"/diagnostic/algorithm::from_particles")) then
             call get_option(trim(s_field%option_path)//"/name", name)
             if (name=="MaterialVolumeFraction") cycle
             call calculate_diagnostics_from_particles(state, i, s_field)
          end if
       end do
    end do
    k = size(state)

    !Initialise diagnostic fields with init_after_particles
    do i = 1,size(state)
       do k = 1,scalar_field_count(state(i))
          s_field => extract_scalar_field(state(i),k)     
          if (have_option(trim(s_field%option_path)//"/diagnostic/init_after_particles")) then
             ! Calculate dependencies
             call calculate_dependencies(state, i, s_field, &
          & dep_states_mask = calculated_states, exclude_nonrecalculated = .false.)
             ! Calculate the diagnostic
             ewrite(2, *) "Calculating diagnostic field: "//trim(state(i)%name)//"::"//trim(s_field%name)
             call calculate_diagnostic_variable(state, i, s_field)
             ! Mark the field as calculated
             call insert(calculated_states(i), s_field, s_field%name)
          end if
       end do
    end do

  end subroutine initialise_particle_diagnostics

  subroutine update_particle_diagnostics(state, time)
    !!Routine to loop over particle arrays and update particle attributes
    !!and diagnostic fields which depend on particles
    use particles, only: particle_lists
    type(state_type), dimension(:), intent(inout) :: state
    real, intent(in) :: time
    type(vector_field), pointer :: xfield
    type(scalar_field), pointer :: s_field
    type(detector_type), pointer :: particle
    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path

    integer :: i, k
    integer :: dim, particle_groups, list_counter, particle_materials
    integer, dimension(:), allocatable :: particle_arrays

    !Check whether there are any particles.
    particle_groups = option_count("/particles/particle_group")
    
    if (particle_groups==0) return

    !Set up particle_lists
    allocate(particle_arrays(particle_groups))
    particle_arrays(:) = 0
    do i = 1,particle_groups
       particle_arrays(i) = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
    end do

    ewrite(2,*), "In update_particle_diagnostics"

    !Allocate parameters
    xfield=>extract_vector_field(state(1), "Coordinate")
    call get_option("/geometry/dimension",dim)
    list_counter = 1

    !Update particle attributes by array
    do i = 1,particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1, particle_arrays(i)
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          if (particle_lists(list_counter)%length.eq.0) then
             list_counter = list_counter + 1
             cycle
          end if
          particle => particle_lists(list_counter)%first
          if (size(particle%attributes).ne.0) then
             call set_particle_attributes(state, dim, time, subgroup_path, particle_lists(list_counter))
          end if
          list_counter = list_counter + 1
       end do
    end do

    !Update diagnostic fields with algorithm "from_particles"
    do i = 1,size(state)
       do k = 1,scalar_field_count(state(i))
          s_field => extract_scalar_field(state(i),k)     
          if (have_option(trim(s_field%option_path)//"/diagnostic/algorithm::from_particles")) then
             call calculate_diagnostics_from_particles(state, i, s_field)
          end if
       end do
    end do
    k = size(state)

    !Update internal MaterialVolumeFraction field
    particle_materials = option_count("material_phase/scalar_field::MaterialVolumeFraction/diagnostic/algorithm::from_particles")
    if (particle_materials.gt.0) then
       s_field => extract_scalar_field(state(k), "MaterialVolumeFraction")
       call calculate_sum_material_volume_fractions(state, s_field)
       call scale(s_field, -1.0)
       call addto(s_field, 1.0)
    end if
    
  end subroutine update_particle_diagnostics

  subroutine initialise_particle_material_fields(state)
    !subroutine to initialise multimaterial fields from
    !particles after mesh adapt  

    type(state_type), dimension(:), intent(inout) :: state

    type(scalar_field), pointer :: s_field
    integer :: i, k
    integer :: particle_groups, particle_materials

    !Check if there are particles
    particle_groups = option_count("/particles/particle_group")

    if (particle_groups==0) return

    !Check if MVF field is generated from particles
    particle_materials = option_count("material_phase/scalar_field::MaterialVolumeFraction/diagnostic/algorithm::from_particles")
    if (particle_materials.gt.0) then
       !Initialise MultialVolumeFraction fields dependent on particles
       do i = 1,size(state)
          s_field => extract_scalar_field(state(i), "MaterialVolumeFraction")     
          if (have_option(trim(s_field%option_path)//"/diagnostic/algorithm::from_particles")) then
             call calculate_diagnostics_from_particles(state, i, s_field)
          end if
       end do
       k = size(state)
       
       !Initialise internal MaterialVolumeFraction field 
       s_field => extract_scalar_field(state(k), "MaterialVolumeFraction")
       call calculate_sum_material_volume_fractions(state, s_field)
       call scale(s_field, -1.0)
       call addto(s_field, 1.0)
    end if
    
  end subroutine initialise_particle_material_fields

  subroutine initialise_particle_diagnostic_fields(state)
    !subroutine to initialise diagnostic fields dependent on 
    !particles after mesh adapt
    type(state_type), dimension(:), intent(inout) :: state
    type(state_type), dimension(size(state)) :: calculated_states

    character(len = OPTION_PATH_LEN) :: name
    type(scalar_field), pointer :: s_field
    integer :: i, k, particle_groups
    
    !Check if there are particles
    particle_groups = option_count("/particles/particle_group")

    if (particle_groups==0) return

    !Initialise diagnostic fields generated from particles
    do i = 1,size(state)
       do k = 1,scalar_field_count(state(i))
          s_field => extract_scalar_field(state(i),k)     
          if (have_option(trim(s_field%option_path)//"/diagnostic/algorithm::from_particles")) then
             call get_option(trim(s_field%option_path)//"/name", name)
             if (name=="MaterialVolumeFraction") cycle
             call calculate_diagnostics_from_particles(state, i, s_field)
          end if
       end do
    end do
    k = size(state)

    !Initialise diagnostic fields with init_after_particles
    do i = 1,size(state)
       do k = 1,scalar_field_count(state(i))
          s_field => extract_scalar_field(state(i),k)     
          if (have_option(trim(s_field%option_path)//"/diagnostic/init_after_particles")) then
             ! Calculate dependencies
             call calculate_dependencies(state, i, s_field, &
          & dep_states_mask = calculated_states, exclude_nonrecalculated = .false.)
             ! Calculate the diagnostic
             ewrite(2, *) "Calculating diagnostic field: "//trim(state(i)%name)//"::"//trim(s_field%name)
             call calculate_diagnostic_variable(state, i, s_field)
             ! Mark the field as calculated
             call insert(calculated_states(i), s_field, s_field%name)
          end if
       end do
    end do

  end subroutine initialise_particle_diagnostic_fields
    
  subroutine calculate_diagnostics_from_particles(states, state_index, s_field)

    !!! Subroutine to determine which method is being used
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
    type(scalar_field), intent(inout) :: s_field
    integer, intent(in) :: state_index
    type(halo_type), pointer :: halo

    character(len=OPTION_PATH_LEN) :: lgroup, lattribute
    type(vector_field), pointer :: xfield
    type(detector_linked_list), allocatable, dimension(:,:) :: node_particles
    type(detector_type), pointer :: particle
    integer :: i
    real, allocatable, dimension(:) :: node_values
    real, allocatable, dimension(:) :: node_part_count ! real instead of integer, so we can use halo_accumulate
    integer :: element, node_number
    real, allocatable, dimension(:) :: local_crds
    integer, dimension(:), pointer :: nodes
    integer :: nprocs
    real :: att_value
    real :: ratio_val

    integer :: group_attribute
    integer, allocatable, dimension(:) :: group_arrays
    
    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/group/name", lgroup)
    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/group/attribute/name", lattribute)
    ewrite(2,*) "Calculate diagnostic field from particle group: ", trim(lgroup), ", attribute: ", trim(lattribute)

    !Initialize field as 0
    s_field%val(:) = 0

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
       node_values(:) = 0
       node_part_count(:) = 0
       
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

       !Spawn and Delete particles if enabled
       if (state_index==1.and.have_option("/particles/particle_group::"//trim(lgroup)//"/particle_spawning")) then
          call particle_cv_check(states, s_field, node_values, node_part_count, node_particles, group_arrays, lgroup, group_attribute, halo, nprocs)
       end if
       
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
       
       deallocate(node_values)
       deallocate(node_part_count)
       deallocate(node_particles)

       ! all values in owned nodes should now be correct
       ! now we need to make sure the halo is updated accordingly
       if (nprocs>1) then
          call halo_update(s_field)
       end if
    else
       !Particles not yet setup, Initial field value will be 0
    end if

  end subroutine calculate_ratio_from_particles

  subroutine particle_cv_check(states, s_field, node_values, node_part_count, node_particles, group_arrays, lgroup, group_attribute, halo, nprocs)
    !Routine to check particle numbers fall within CV thresholds
    !Spawns or deletes particles if numbers exceed or fall below CV thresholds

    use particles, only: particle_lists
    type(state_type), dimension(:), target, intent(in) :: states
    type(scalar_field), intent(in) :: s_field
    type(halo_type), pointer, intent(inout) :: halo
    character(len=OPTION_PATH_LEN), intent(in) :: lgroup
    type(detector_linked_list), dimension(:,:), intent(inout) :: node_particles
    real, dimension(:), intent(inout) :: node_values, node_part_count
    integer, dimension(:), intent(inout) :: group_arrays
    integer, intent(in) :: group_attribute
    integer, intent(in) :: nprocs
    type(vector_field), pointer :: xfield
    integer :: i, j

    integer, allocatable, dimension(:) :: summed_particles, add_particles, remove_particles
    real, dimension(:), allocatable :: temp_values, temp_part_count
    integer :: min_thresh
    integer :: max_thresh
    
    xfield=>extract_vector_field(states(1), "Coordinate")

    !Set minimum and maximum particle thresholds per control volume
    if (have_option("/particles/particle_group::"//trim(lgroup)//"/particle_spawning")) then
       call get_option("/particles/particle_group::"//trim(lgroup)//"/particle_spawning/min_cv_threshhold", min_thresh)
       call get_option("/particles/particle_group::"//trim(lgroup)//"/particle_spawning/max_cv_threshhold", max_thresh)
    end if

    allocate(summed_particles(size(group_arrays)))
    allocate(add_particles(size(group_arrays)))
    allocate(remove_particles(size(group_arrays)))
    allocate(temp_part_count(size(node_part_count)))
    allocate(temp_values(size(node_values)))

    summed_particles(:) = 0
    temp_part_count(:) = 0
    temp_values(:) = 0

    !Loop over all nodes
    do i = 1,node_count(s_field)
       !Count number of particles per node and ensure thresholds are not broken
       if (node_part_count(i)<min_thresh) then
          if (node_part_count(i)>0) then
             if (node_part_count(i)<min_thresh/2) then
                call multi_spawn_particles(temp_part_count(i), temp_values(i), node_particles(:,i), group_arrays, group_attribute, xfield, add_particles, summed_particles, i, 3)
             else 
                call spawn_particles(temp_part_count(i), temp_values(i), node_particles(:,i), group_arrays, group_attribute, xfield, add_particles, i)
                summed_particles=summed_particles+add_particles
             end if
          else if (node_part_count(i)==0) then!need to add something here for parallel (only if node owned?) check in more detail...
             if (node_owned(s_field, i)) then
                call spawn_zero_particles(temp_part_count(:), temp_values(:), node_particles(:,:), group_arrays, group_attribute, xfield, add_particles, i, min_thresh)
                summed_particles=summed_particles+add_particles
             end if
          end if
       end if
    end do
    if (nprocs>1) then
       call halo_accumulate(halo, temp_part_count)
       call halo_accumulate(halo, temp_values)
    end if

    node_part_count(:) = node_part_count(:) + temp_part_count(:)
    node_values(:) = node_values(:) + temp_values(:)
    temp_part_count(:) = 0
    temp_values(:) = 0

    do i = 1,node_count(s_field)
       if (node_part_count(i)>max_thresh) then
          call delete_particles(temp_part_count(i), temp_values(i), node_particles(:,i), group_arrays, group_attribute, remove_particles)
          summed_particles=summed_particles-remove_particles
       end if
    end do
    if (nprocs>1) then
       call halo_accumulate(halo, temp_part_count)
       call halo_accumulate(halo, temp_values)
    end if
    node_part_count(:) = node_part_count(:) + temp_part_count(:)
    node_values(:) = node_values(:) + temp_values(:)
    
    do j=1,size(group_arrays)
       call allsum(summed_particles(j))
       particle_lists(group_arrays(j))%total_num_det=particle_lists(group_arrays(j))%total_num_det+summed_particles(j)
    end do
    
    deallocate(add_particles)
    deallocate(remove_particles)
    deallocate(summed_particles)
    deallocate(temp_part_count)
    deallocate(temp_values)
    
  end subroutine particle_cv_check

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
    type(detector_type), pointer :: particle
    integer :: i
    real, allocatable, dimension(:) :: node_part_count ! real instead of integer, so we can use halo_accumulate
    integer :: element, node_number
    real, allocatable, dimension(:) :: local_crds
    integer, dimension(:), pointer :: nodes
    integer :: nprocs

    integer :: group_attribute
    integer, allocatable, dimension(:) :: group_arrays


    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/method/group/name", lgroup)
    ewrite(2,*) "Calculate number of particles from particle group: ", trim(lgroup)

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

  subroutine multi_spawn_particles(node_part_count, node_values, node_particles, group_arrays, group_attribute, xfield, add_particles, summed_particles, node_num, mult)
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    real, intent(inout) :: node_values
    real, intent(inout) :: node_part_count
    integer, intent(in), dimension(:) :: group_arrays
    integer, intent(in) :: group_attribute
    integer, intent(in) :: node_num
    type(vector_field), pointer, intent(in) :: xfield
    integer, dimension(:), intent(inout) :: add_particles
    integer, dimension(:), intent(inout) :: summed_particles
    integer, intent(in) :: mult

    integer :: i

    do i = 1,mult
       call spawn_particles(node_part_count, node_values, node_particles, group_arrays, group_attribute, xfield, add_particles, node_num)
       summed_particles=summed_particles+add_particles
    end do

  end subroutine multi_spawn_particles

  subroutine spawn_particles(node_part_count, node_values, node_particles, group_arrays, group_attribute, xfield, add_particles, node_num)
    !Subroutine to calculate the number of particles in each control volume, and spawn additional
    !particles within a control volume if a minimum number of particles is not met

    use particles, only: particle_lists
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    real, intent(inout) :: node_values
    real, intent(inout) :: node_part_count
    integer, intent(in), dimension(:) :: group_arrays
    integer, intent(in) :: group_attribute
    integer, intent(in) :: node_num
    type(vector_field), pointer, intent(in) :: xfield
    integer, dimension(:), intent(inout) :: add_particles
    
    type(detector_type), pointer :: particle
    type(detector_type), pointer :: temp_part

    real, dimension(:), allocatable :: node_coord
    real :: rand_lcoord

    character(len=OPTION_PATH_LEN) :: name
    character(len=OPTION_PATH_LEN) :: particle_name
    integer, dimension(:), pointer :: ele_nums
    integer, allocatable, dimension(:) :: node_numbers
    integer :: ele_spawn
    integer :: id, name_len, tot_len
    integer :: j, i
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
          allocate(node_numbers(xfield%dim+1))
       end if

       particle => node_particles(j)%first
       do while(associated(particle))
          !duplicate particles, spawn in random element adjacent
          !to CV ensuring particle falls within CV
          particle_name = trim(name)//int2str(id+1)
          temp_part => null()
          call allocate(temp_part, size(particle%position), size(particle%local_coords), attribute_size)

          temp_part%name = trim(particle_name)
          temp_part%id_number = id+1
          temp_part%list_id = particle%list_id

          !Get ele numbers of adjacent elements
          ele_nums => node_neigh(xfield, node_num)
          
          !randomly select adjacent element to node
          ele_spawn=floor(rand()*size(ele_nums)+1)
          temp_part%element = ele_nums(ele_spawn)
          node_numbers(:) = ele_nodes(xfield, ele_nums(ele_spawn))
          !randomly select local coords within the element, ensuring coords are within cv
          do i = 1,xfield%dim+1
             if (node_num==node_numbers(i)) then
                rand_lcoord=rand()/(1/0.45)+0.55!lcoords for cv range from 0.55<x<1
                node_coord(:)=(1-rand_lcoord)/(xfield%dim)
                node_coord(i)=rand_lcoord            
                temp_part%local_coords=node_coord
             end if
          end do
          
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
          deallocate(node_numbers)
       end if

    end do

  end subroutine spawn_particles

  subroutine spawn_zero_particles(node_part_count, node_values, node_particles, group_arrays, group_attribute, xfield, add_particles, node_num, min_thresh)
    !Subroutine to spawn particles in a control volume with 0 'parent' particles in CV

    use particles, only: particle_lists
    type(detector_linked_list), intent(inout), dimension(:,:) :: node_particles
    real, intent(inout), dimension(:) :: node_values
    real, intent(inout), dimension(:) :: node_part_count
    integer, intent(in), dimension(:) :: group_arrays
    integer, intent(in) :: group_attribute
    type(vector_field), pointer, intent(in) :: xfield
    integer, dimension(:), intent(inout) :: add_particles
    integer, intent(in) :: node_num
    integer, intent(in) :: min_thresh
    integer, dimension(:), pointer :: ele_nums
    integer, allocatable, dimension(:,:) :: node_numbers
    integer, allocatable, dimension(:) :: part_counter
    real :: part_total
    real :: rng_spawn
    logical :: spawned
    integer :: prev_parts, ele_spawn
    
    type(detector_type), pointer :: particle
    type(detector_type), pointer :: temp_part

    real, dimension(:), allocatable :: node_coord
    real :: rand_lcoord

    character(len=OPTION_PATH_LEN) :: name
    character(len=OPTION_PATH_LEN) :: particle_name
    integer :: id, name_len, tot_len
    integer :: i, j, k
    integer, dimension(3) :: attribute_size
    
    add_particles(:) = 0

    !Get ele numbers of adjacent elements
    ele_nums => node_neigh(xfield, node_num)

    allocate(node_numbers(size(ele_nums),xfield%dim+1))
    allocate(part_counter(size(group_arrays)))
    part_counter(:) = 0

    !Get node numbers for each adjacent element
    do j = 1,size(ele_nums)
       node_numbers(j,:) = ele_nodes(xfield,ele_nums(j))
    end do

    !Count the total number of particles in each surrounding CV
    do i = 1,size(group_arrays)
       do j = 1,size(ele_nums)
          do k=1,xfield%dim+1
             part_counter(i) = part_counter(i) + node_particles(i,node_numbers(j,k))%length
          end do
       end do
    end do

    part_total=sum(part_counter)
    if (part_total==0) then
       ewrite(2,*) "There are no particles present in adjacent CV's"
       return
    end if

    !If adjacent CV's contain particles, spawn particles in current CV
    !with weighted values depending on the number and group of
    !particles in adjacent CV's
    do i = 1,min_thresh !spawn particles until min_thresh is met
       spawned=.false.
       k=1
       do while (spawned.eqv..false.)
          rng_spawn=rand()
          prev_parts=0
          !Adjust the weight of the spawn chance dependent on the
          !group being spawned
          do j = 1,k-1
             prev_parts=prev_parts+part_counter(j)
          end do
          !Check if random spawn conditions are met
          if (rng_spawn<=part_counter(k)/(part_total-prev_parts)) then
             !Spawn particle from this group
             particle => particle_lists(group_arrays(k))%last
             id = particle%id_number
             name_len = len(int2str(id+1))+1 !id length + '_'
             tot_len = len(trim(particle%name))
             name = trim(particle%name(1:tot_len-name_len))
             attribute_size(1) = size(particle%attributes)
             attribute_size(2) = size(particle%old_attributes)
             attribute_size(3) = size(particle%old_fields)
             allocate(node_coord(size(particle%local_coords)))

             particle_name = trim(name)//int2str(id+1)
             temp_part => null()
             call allocate(temp_part, size(particle%position), size(particle%local_coords), attribute_size)
             temp_part%name = trim(particle_name)
             temp_part%id_number = id+1
             temp_part%list_id = particle%list_id
             
             !randomly select adjacent element to node
             ele_spawn=floor(rand()*size(ele_nums)+1)
             temp_part%element = ele_nums(ele_spawn)
             !randomly select local coords within the element, ensuring coords are within cv
             do j = 1,xfield%dim+1
                if (node_num==node_numbers(ele_spawn,j)) then
                   rand_lcoord=rand()/(1/0.45)+0.55!lcoords for cv range from 0.55<x<1
                   node_coord(:)=(1-rand_lcoord)/(xfield%dim)
                   node_coord(j)=rand_lcoord            
                   temp_part%local_coords= node_coord
                end if
             end do
 
             call local_to_global(xfield, temp_part%local_coords, temp_part%element, temp_part%position)
             temp_part%type = LAGRANGIAN_DETECTOR
             temp_part%attributes(:) = 0
             temp_part%old_attributes(:) = 0
             temp_part%old_fields(:) = 0
             temp_part%attributes(group_attribute) = particle%attributes(group_attribute)
             node_values(node_num) = node_values(node_num) + temp_part%attributes(group_attribute)
             node_part_count(node_num) = node_part_count(node_num) + 1
             
             temp_part%previous => particle_lists(group_arrays(k))%last
             particle_lists(group_arrays(k))%last%next => temp_part
             particle_lists(group_arrays(k))%last => temp_part
             particle_lists(group_arrays(k))%last%next => null()
             particle_lists(group_arrays(k))%length = particle_lists(group_arrays(k))%length + 1
             add_particles(k) = add_particles(k) + 1
             particle => null()

             if (allocated(node_coord)) then
                deallocate(node_coord)
             end if
             
             spawned=.true.
          else
             k=k+1
          end if
       end do
    end do

  end subroutine spawn_zero_particles

  subroutine local_to_global(coordinates, l_coords, ele_number, global_coord)
    !Subroutine to convert local coordinates within a specified element
    !to global coordinates within the model domain
    
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
                !Remove from particle list
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
                !Remove from temp list
                if (associated(particle%temp_previous)) then
                   particle%temp_previous%temp_next => particle%temp_next
                else
                   node_particles(j)%first => particle%temp_next
                end if
                if (associated(particle%temp_next)) then
                   particle%temp_next%temp_previous => particle%temp_previous
                else
                   node_particles(j)%last => particle%temp_previous
                end if
                remove_particles(j)=remove_particles(j)+1
                node_particles(j)%length = node_particles(j)%length -1
                particle_lists(group_arrays(j))%length = particle_lists(group_arrays(j))%length - 1
                node_values = node_values - particle%attributes(group_attribute)
                node_part_count = node_part_count - 1
                call deallocate(particle)
             else
                rand3=rand()
                if (rand3>0.5) then
                   !Remove from particle list
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
                   !Remove from temp list
                   if (associated(particle%temp_previous)) then
                      particle%temp_previous%temp_next => particle%temp_next
                   else
                      node_particles(j)%first => particle%temp_next
                   end if
                   if (associated(particle%temp_next)) then
                      particle%temp_next%temp_previous => particle%temp_previous
                   else
                      node_particles(j)%last => particle%temp_previous
                   end if
                   remove_particles(j)=remove_particles(j)+1
                   node_particles(j)%length = node_particles(j)%length -1
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
                   !Remove from particle list
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
                   !Remove from temp list
                   if (associated(particle%temp_previous)) then
                      particle%temp_previous%temp_next => particle%temp_next
                   else
                      node_particles(j)%first => particle%temp_next
                   end if
                   if (associated(particle%temp_next)) then
                      particle%temp_next%temp_previous => particle%temp_previous
                   else
                      node_particles(j)%last => particle%temp_previous
                   end if
                   remove_particles(j)=remove_particles(j)+1
                   node_particles(j)%length = node_particles(j)%length -1
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
