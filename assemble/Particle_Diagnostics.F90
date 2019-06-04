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
  use global_parameters, only : OPTION_PATH_LEN, FIELD_NAME_LEN
  use futils, only: int2str
  use spud
  use parallel_tools, only: getnprocs, allsum
  use elements, only: eval_shape
  use fields_base, only: ele_val, ele_loc
  use parallel_fields, only: node_owned
  use fields_calculations, only: dot_product
  use fields
  use state_module
  use halos
  use detector_data_types
  use pickers
  use field_options
  use detector_tools, only: temp_insert, insert, allocate, deallocate, temp_deallocate
  use particles, only : get_particle_arrays, update_list_lengths, &
& update_particle_subgroup_attributes_and_fields, initialise_constant_particle_attributes
  use multimaterial_module, only: calculate_sum_material_volume_fractions
  
  implicit none

  private

  public :: initialise_particle_diagnostics, update_particle_diagnostics, &
       & calculate_diagnostics_from_particles, calculate_ratio_from_particles, &
       & calculate_numbers_from_particles, initialise_constant_particle_diagnostics, &
       & initialise_particle_diagnostic_fields_post_adapt, particle_cv_check

  contains

  subroutine initialise_constant_particle_diagnostics(state)
    !subroutine to initialise constant particle attributes and
    !MVF fields if 'from_particles'

    use particles, only: particle_lists
    type(state_type), dimension(:), intent(inout) :: state

    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path
    type(scalar_field), pointer :: s_field
    integer :: i, k
    integer :: particle_groups, list_counter, particle_materials
    integer, dimension(:), allocatable :: particle_arrays

    type(detector_type), pointer :: particle

    ewrite(2,*) "In initialise_constant_particle_diagnostics"

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

    !Initialise fields based on number of particles if present
    do i = 1,size(state)
       do k = 1,scalar_field_count(state(i))
          s_field => extract_scalar_field(state(i),k)     
          if (have_option(trim(s_field%option_path)//"/diagnostic/algorithm::number_of_particles")) then
             call calculate_diagnostics_from_particles(state, i, s_field)
          end if
       end do
    end do
    k = size(state)

    deallocate(particle_arrays)
    
  end subroutine initialise_constant_particle_diagnostics

  subroutine initialise_particle_diagnostics(state)
    !subroutine to initialise particle attributes, diagnostic fields
    !dependent on particles, and diagnostic fields to be set after
    !particles are initialised
    use particles, only: particle_lists
    type(state_type), dimension(:), intent(inout) :: state
    
    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path, name
    type(vector_field), pointer :: xfield
    type(scalar_field), pointer :: s_field
    real :: current_time, dt
    integer :: i, k
    integer :: dim, particle_groups, list_counter
    integer, dimension(:), allocatable :: particle_arrays
    type(detector_type), pointer :: particle
    logical :: from_file

    ewrite(2,*) "In initialise_particle_diagnostics"
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
    call get_option("/timestepping/timestep", dt)

    !Initialise particle attributes
    list_counter=1
    do i = 1,particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do k = 1,particle_arrays(i)
          particle => particle_lists(list_counter)%first
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
          !Only initialise particle attributes if not loaded from a file
          from_file=have_option(trim(subgroup_path)//"/initial_position/from_file")
          if (from_file) then
             cycle
          end if
          if (option_count(trim(subgroup_path) // "/attributes/attribute/python").gt.0 .or. &
              & option_count(trim(subgroup_path) // "/attributes/attribute/python_fields").gt.0) then
             call update_particle_subgroup_attributes_and_fields(state, current_time, dt, subgroup_path, particle_lists(list_counter))
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

    deallocate(particle_arrays)

  end subroutine initialise_particle_diagnostics

  subroutine update_particle_diagnostics(state, time, dt)
    !!Routine to loop over particle arrays and update particle attributes
    !!and diagnostic fields which depend on particles
    use particles, only: particle_lists
    type(state_type), dimension(:), intent(inout) :: state
    real, intent(in) :: time
    real, intent(in) :: dt
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
             call update_particle_subgroup_attributes_and_fields(state, time, dt, subgroup_path, particle_lists(list_counter))
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

    !Update diagnostic fields with algorithm "number_of_particles"
    do i = 1,size(state)
       do k = 1,scalar_field_count(state(i))
          s_field => extract_scalar_field(state(i),k)     
          if (have_option(trim(s_field%option_path)//"/diagnostic/algorithm::number_of_particles")) then
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

    deallocate(particle_arrays)
    
  end subroutine update_particle_diagnostics

  subroutine initialise_particle_diagnostic_fields_post_adapt(state)
    !subroutine to initialise diagnostic fields dependent on 
    !particles after mesh adapt
    type(state_type), dimension(:), intent(inout) :: state

    character(len = OPTION_PATH_LEN) :: name
    type(scalar_field), pointer :: s_field
    integer :: i, k, particle_groups, particle_materials
    
    !Check if there are particles
    particle_groups = option_count("/particles/particle_group")

    if (particle_groups==0) return

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

    !Initialise diagnostic fields based on the number of particles
    do i = 1,size(state)
       do k = 1,scalar_field_count(state(i))
          s_field => extract_scalar_field(state(i),k)     
          if (have_option(trim(s_field%option_path)//"/diagnostic/algorithm::number_of_particles")) then
             call calculate_diagnostics_from_particles(state, i, s_field)
          end if
       end do
    end do
    k = size(state)

  end subroutine initialise_particle_diagnostic_fields_post_adapt
    
  subroutine calculate_diagnostics_from_particles(states, state_index, s_field)

    !!! Subroutine to determine which method is being used
    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field

    character(len= OPTION_PATH_LEN) :: lmethod

    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/name", lmethod, default = "from_particles")

    select case(trim(lmethod))
       
    case("from_particles")
        call calculate_ratio_from_particles(states, state_index, s_field)
       
    case("number_of_particles")
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
    type(detector_type), pointer :: particle
    integer :: i, j
    real, allocatable, dimension(:) :: node_values
    real, allocatable, dimension(:) :: node_part_count
    integer :: element, node_number
    real, allocatable, dimension(:) :: local_crds
    integer, dimension(:), pointer :: nodes
    integer :: nprocs
    real :: att_value, ratio_val

    integer :: group_attribute
    integer, allocatable, dimension(:) :: group_arrays
    
    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/particle_group/name", lgroup)
    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/particle_group/particle_attribute/name", lattribute)
    ewrite(2,*) "Calculate diagnostic field from particle group: ", trim(lgroup), ", attribute: ", trim(lattribute)

    !Initialize field as 0
    s_field%val(:) = 0

    xfield=>extract_vector_field(states(1), "Coordinate")

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
       allocate(node_values(node_count(s_field)))
       allocate(node_part_count(node_count(s_field)))
       node_values(:) = 0
       node_part_count(:) = 0

       if (have_option(trim(complete_field_path(s_field%option_path))// "/algorithm/interpolation/weighted_distance")) then

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

                !Distribute particle values across nodes of element by weighted distance function
                do j = 1, size(nodes)
                   node_number = nodes(j)
                   node_values(node_number) = node_values(node_number) + att_value*local_crds(j)
                   node_part_count(node_number) = node_part_count(node_number) + 1.0*local_crds(j)
                end do
                particle => particle%next
                
             end do
             if (allocated(local_crds)) then
                deallocate(local_crds)
             end if
          end do

       else if (have_option(trim(complete_field_path(s_field%option_path))// "/algorithm/interpolation/nearest_neighbour")) then
       
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
                particle => particle%next
                
             end do
             if (allocated(local_crds)) then
                deallocate(local_crds)
             end if
          end do
       end if
          
       if (nprocs>1) then
          call halo_accumulate(halo, node_part_count)
          call halo_accumulate(halo, node_values)
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
       deallocate(group_arrays)

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
    type(detector_type), pointer :: particle
    integer :: i
    real, allocatable, dimension(:) :: node_part_count ! real instead of integer, so we can use halo_accumulate
    integer :: element, node_number
    real, allocatable, dimension(:) :: local_crds
    integer, dimension(:), pointer :: nodes
    integer :: nprocs

    integer, allocatable, dimension(:) :: group_arrays

    logical :: have_subgroup, have_attribute

    have_subgroup=.false.
    have_attribute=.false.

    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/particle_group/name", lgroup)

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
       call get_particle_arrays(lgroup, group_arrays)

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
       deallocate(group_arrays)
    else
       !Particles not yet setup, Initial field value will be 0
    end if

  end subroutine calculate_numbers_from_particles

  subroutine particle_cv_check(states)
    !Routine to check particle numbers fall within CV thresholds
    !Spawns or deletes particles if numbers exceed or fall below CV thresholds
    
    type(state_type), dimension(:), target, intent(in) :: states
    
    type(mesh_type), pointer :: mesh
    integer, allocatable, dimension(:) :: group_arrays
    integer :: particle_groups, k
    character(len=OPTION_PATH_LEN) :: mesh_name
    character(len=OPTION_PATH_LEN) :: particle_group
    
    logical :: have_attribute_caps, have_radius
    integer :: min_thresh, max_thresh
    real :: radius, cap_percent

    ewrite(2,*) "In particle_cv_check"

    !Check if there are particles
    particle_groups = option_count('/particles/particle_group')

    if (particle_groups==0) return
    do k = 1, particle_groups
       if (have_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning")) then
          !Get the mesh particles will be spawned to
          call get_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/mesh/name", mesh_name)
          mesh => extract_mesh(states(1), trim(mesh_name))
          
          !Set minimum and maximum particle thresholds per control volume
          call get_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/min_cv_threshhold", min_thresh)
          call get_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/max_cv_threshhold", max_thresh)

          !Check option for where particles should be spawned within a control volume.
          if (have_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/spawn_location/radius")) then
             have_radius=.true.
             call get_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/spawn_location/radius", radius)
          else
             have_radius=.false.
          end if

          !Check option on whether certain particle subgroup spawning should be capped if one subgroup is dominant
          have_attribute_caps= .false.
          if (have_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/subgroup_spawning_caps")) then
             have_attribute_caps=.true.
             call get_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/subgroup_spawning_caps/percentage", cap_percent)
          end if

          !Need to get particle arrays
          call get_option("/particles/particle_group["//int2str(k-1)//"]/name", particle_group)
          call get_particle_arrays(particle_group, group_arrays)
          if (have_attribute_caps.eqv..true.) then
             call spawn_delete_particles(states, mesh, group_arrays, max_thresh, min_thresh, have_radius, radius, cap_percent)
          else
             call spawn_delete_particles(states, mesh, group_arrays, max_thresh, min_thresh, have_radius, radius)
          end if
          deallocate(group_arrays)
       end if
    end do
    
  end subroutine particle_cv_check

  subroutine spawn_delete_particles(states, mesh, group_arrays, max_thresh, min_thresh, have_radius, radius, cap_percent)
    use particles, only: particle_lists
    !Spawn and delete particles based on the parameters assigned in particle_cv_check
    type(state_type), dimension(:), target, intent(in) :: states
    type (mesh_type), intent(in) :: mesh
    integer, dimension(:), intent(in) :: group_arrays
    integer, intent(in) :: min_thresh, max_thresh
    logical, intent(in) :: have_radius
    real, intent(in) :: radius
    real, optional, intent(in) :: cap_percent
    
    type(halo_type), pointer :: halo
    type(detector_linked_list), allocatable, target, dimension(:,:) :: node_particles
    type(detector_linked_list), pointer :: del_node_particles
    real, allocatable, dimension(:) :: node_part_count
    real, allocatable, dimension(:) :: local_crds
    type(detector_type), pointer :: particle
    integer :: element, node_number
    integer, dimension(:), pointer :: nodes
    integer :: nprocs
    type(vector_field), pointer :: xfield
    integer :: i, j, mult

    integer, allocatable, dimension(:) :: summed_particles, add_particles, remove_particles
    real, dimension(:), allocatable :: temp_part_count
    integer :: total_particles

    xfield=>extract_vector_field(states(1), "Coordinate")
    nprocs = getnprocs()
    if (nprocs>1) then
       halo => mesh%halos(2)
    end if
    allocate(node_particles(size(group_arrays),node_count(mesh)))
    allocate(node_part_count(node_count(mesh)))
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
          
          !Find nodes for specified element
          nodes => ele_nodes(mesh, element)
          
          ! work out nearest node based on max. local coordinate
          node_number = nodes(maxloc(local_crds, dim=1))
          
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
    end if

    allocate(summed_particles(size(group_arrays)))
    allocate(add_particles(size(group_arrays)))
    allocate(remove_particles(size(group_arrays)))
    allocate(temp_part_count(size(node_part_count)))

    summed_particles(:) = 0
    temp_part_count(:) = 0


    !Loop over all nodes
    do i = 1,node_count(mesh)
       !Count number of particles per node and ensure thresholds are not broken
       if (node_part_count(i)<min_thresh) then
          if (node_part_count(i)>0) then
             if (node_part_count(i)<min_thresh/2) then
                mult = nint(min_thresh/node_part_count(i)*1.0)
                call multi_spawn_particles(temp_part_count(i), node_particles(:,i), group_arrays, xfield, summed_particles, i, mult, have_radius, radius, cap_percent)
             else
                call spawn_particles(temp_part_count(i), node_particles(:,i), group_arrays, xfield, add_particles, i, have_radius, radius, cap_percent)
                summed_particles=summed_particles+add_particles
             end if
          else if (node_part_count(i)==0) then
             if (node_owned(mesh, i)) then
                call spawn_zero_particles(temp_part_count(:), node_particles(:,:), group_arrays, xfield, add_particles, i, cap_percent, max_thresh)
                summed_particles=summed_particles+add_particles
             end if
          end if
       end if
    end do
    if (nprocs>1) then
       call halo_accumulate(halo, temp_part_count)
    end if

    node_part_count(:) = node_part_count(:) + temp_part_count(:)
    temp_part_count(:) = 0

    do i = 1,node_count(mesh)
       if (node_part_count(i)>max_thresh) then
          if (node_part_count(i)>2*max_thresh) then
             mult = nint(node_part_count(i)*1.0/max_thresh)
             call multi_delete_particles(mult, temp_part_count(i), node_particles(:,i), group_arrays, summed_particles)
          else
             call delete_particles(temp_part_count(i), node_particles(:,i), group_arrays, remove_particles)
             summed_particles=summed_particles-remove_particles
          end if
       end if
    end do
    if (nprocs>1) then
       call halo_accumulate(halo, temp_part_count)
    end if
    node_part_count(:) = node_part_count(:) + temp_part_count(:)
    
    do j=1,size(group_arrays)
       call allsum(summed_particles(j))
       particle_lists(group_arrays(j))%total_num_det=particle_lists(group_arrays(j))%total_num_det+summed_particles(j)
    end do

    !Sanity check
    do j=1,size(group_arrays)
       total_particles = particle_lists(group_arrays(j))%length
       call allsum(total_particles)
       assert(total_particles==particle_lists(group_arrays(j))%total_num_det)
    end do

    deallocate(node_part_count)
    deallocate(add_particles)
    deallocate(remove_particles)
    deallocate(summed_particles)
    deallocate(temp_part_count)


    do i=1,node_count(mesh)
       do j=1,size(group_arrays)
          del_node_particles => node_particles(j,i)
          call temp_deallocate(del_node_particles)
       end do
    end do
    deallocate(node_particles)


  end subroutine spawn_delete_particles

  subroutine multi_spawn_particles(node_part_count, node_particles, group_arrays, xfield, summed_particles, node_num, mult, have_radius, radius, cap_percent)
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    real, intent(inout) :: node_part_count
    integer, intent(in), dimension(:) :: group_arrays
    integer, intent(in) :: node_num
    type(vector_field), pointer, intent(in) :: xfield
    integer, dimension(:), intent(inout) :: summed_particles
    integer, intent(in) :: mult
    logical, intent(in) :: have_radius
    real, intent(in) :: radius
    real, optional, intent(in) :: cap_percent

    integer :: i, power, j
    integer, allocatable, dimension(:) :: add_particles
    logical :: power_set

    allocate(add_particles(size(group_arrays)))

    power_set=.false.
    j=0
    do while (power_set.eqv..false.)
       if (mult>=2**j) then
          j=j+1
       else
          power=j-1
          power_set=.true.
       end if
    end do

    do i = 1,power
       call spawn_particles(node_part_count, node_particles, group_arrays, xfield, add_particles, node_num, have_radius, radius, cap_percent)
       summed_particles=summed_particles+add_particles
    end do

    deallocate(add_particles)

  end subroutine multi_spawn_particles

  subroutine spawn_particles(node_part_count, node_particles, group_arrays, xfield, add_particles, node_num, have_radius, radius, cap_percent)
    !Subroutine to calculate the number of particles in each control volume, and spawn additional
    !particles within a control volume if a minimum number of particles is not met

    use particles, only: particle_lists
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    real, intent(inout) :: node_part_count
    integer, intent(in), dimension(:) :: group_arrays
    integer, intent(in) :: node_num
    type(vector_field), pointer, intent(in) :: xfield
    integer, dimension(:), intent(inout) :: add_particles
    logical, intent(in) :: have_radius
    real, intent(in) :: radius
    real, optional, intent(in) :: cap_percent
    
    type(detector_type), pointer :: particle
    type(detector_type), pointer :: temp_part

    real, dimension(:), allocatable :: node_coord

    character(len=OPTION_PATH_LEN) :: name
    character(len=OPTION_PATH_LEN) :: particle_name
    integer, dimension(:), allocatable :: node_numbers
    integer, dimension(:), pointer :: ele_nums
    integer :: id, name_len, tot_len, group_spawn, ele_spawn
    integer :: j, i, k, l
    integer, dimension(3) :: attribute_size
    logical :: spawn_group, coords_set
    real :: max_lcoord, rand_lcoord, sum_coords
    real, dimension(:), allocatable :: rand_lcoords

    !Check if particle node values of above or below a threshold, only spawn from
    !group containing those node values if they are

    spawn_group = .false.
    if (present(cap_percent)) then
       do i=1,size(group_arrays)
          if ((((node_particles(i)%length*1.0)/sum(node_particles(:)%length)*1.0)*100)>cap_percent) then
             spawn_group = .true.
             group_spawn = i
          end if
       end do
    end if


    add_particles(:) = 0
    !Loop over particle groups for spawning
    do j = 1,size(group_arrays)
       temp_part => particle_lists(group_arrays(j))%last
       if (associated(temp_part)) then
          if (spawn_group.eqv..true.) then
             if (j/=group_spawn) cycle
          end if
             
          id = temp_part%id_number
          name_len = len(int2str(id+1))+1 !id length + '_'
          tot_len = len(trim(temp_part%name))
          name = trim(temp_part%name(1:tot_len-name_len))
          attribute_size(1) = size(temp_part%attributes)
          attribute_size(2) = size(temp_part%old_attributes)
          attribute_size(3) = size(temp_part%old_fields)
          allocate(node_coord(size(temp_part%local_coords)))
          allocate(node_numbers(size(temp_part%local_coords)))
       end if

       !Duplicate parent particles
       particle => node_particles(j)%first
       do while(associated(particle))
          particle_name = trim(name)//int2str(id+1)
          temp_part => null()
          call allocate(temp_part, size(particle%position), size(particle%local_coords), attribute_size)

          temp_part%name = trim(particle_name)
          temp_part%id_number = id+1
          temp_part%list_id = particle%list_id

          !Check if particles are spawning within a radius around their parent, or randomly within the CV
          if (have_radius.eqv..true.) then
             !max_lcoord=maxval(particle%local_coords)
             temp_part%element=particle%element
             node_numbers(:) = ele_nodes(xfield, temp_part%element)
             !randomly select local coords within radius around parent particle
             !ensuring new particle falls within cv
             allocate(rand_lcoords(size(node_coord)))
             coords_set=.false.
             !temporary dummy counter l for troubleshooting
             l = 0
             do while (coords_set.eqv..false.)
                rand_lcoord=0
                do i = 1,size(node_coord)
                   if (node_numbers(i)==node_num) then
                      k=i
                      cycle
                   end if
                   rand_lcoords(i) = (rand()-0.5)*2*(radius/(size(particle%local_coords)-1))
                   rand_lcoord = rand_lcoord + rand_lcoords(i)
                end do
                rand_lcoords(k) = -rand_lcoord
                node_coord(:) = particle%local_coords(:) + rand_lcoords(:)
                if (sum(node_coord)/=1) then
                   sum_coords=0
                   do i = 1,size(node_coord)
                      if (i==k) cycle
                      sum_coords=sum_coords+node_coord(i)
                   end do
                   node_coord(k)=1-sum_coords
                end if
                coords_set=.true.
                l = l + 1
                if ((maxloc(node_coord, 1)/=k).or.node_coord(k)>1.0) then !Not within CV or element, try again
                   coords_set=.false.
                end if
                if (minval(node_coord)<0.0) then !not within element, try again
                   coords_set=.false.
                end if
                if (l>100) then !Catch to stop infinite looping
                   ewrite(2,*) "loop limit reached, setting spawned coords to parent coords"
                   node_coord(:) = particle%local_coords(:)
                   coords_set=.true.
                end if
             end do
             temp_part%local_coords = node_coord
             deallocate(rand_lcoords)
          else !Particles will spawn randomly in the CV
             !randomly select adjacent element to node
             !Get ele numbers of adjacent elements
             ele_nums => node_neigh(xfield, node_num)
             ele_spawn=floor(rand()*size(ele_nums)+1)
             temp_part%element = ele_nums(ele_spawn)
             node_numbers(:) = ele_nodes(xfield, ele_nums(ele_spawn))
             max_lcoord=rand()/(1/0.45)+0.51!lcoords for cv range from 0.51<x<1
             call set_spawned_lcoords(max_lcoord, node_coord, node_num, node_numbers)
             temp_part%local_coords=node_coord
          end if
          
          !Convert local particle coordinates to global coordinates
          call local_to_global(xfield, temp_part%local_coords, temp_part%element, temp_part%position)
          temp_part%type = LAGRANGIAN_DETECTOR
          !Copy parent particle attributes
          temp_part%attributes(:) = particle%attributes(:)
          temp_part%old_attributes(:) = particle%old_attributes(:)
          temp_part%old_fields(:) = particle%old_fields(:)
          node_part_count = node_part_count + 1

          !add particle to relevant particle list
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

  subroutine spawn_zero_particles(node_part_count, node_particles, group_arrays, xfield, add_particles, node_num, cap_percent, max_thresh)
    !Subroutine to spawn particles in a control volume with 0 'parent' particles in CV
    use particles, only: particle_lists
    type(detector_linked_list), intent(inout), dimension(:,:) :: node_particles
    real, intent(inout), dimension(:) :: node_part_count
    integer, intent(in), dimension(:) :: group_arrays
    type(vector_field), pointer, intent(in) :: xfield
    integer, dimension(:), intent(inout) :: add_particles
    integer, intent(in) :: node_num
    real, optional, intent(in) :: cap_percent
    integer, intent(in) :: max_thresh
    
    integer, dimension(:), pointer :: ele_nums
    integer, allocatable, dimension(:,:) :: node_numbers
    integer, allocatable, dimension(:) :: part_counter
    real :: part_total
    logical :: spawn_group
    
    type(detector_type), pointer :: particle
    type(detector_type), pointer :: temp_part

    real, dimension(:), allocatable :: node_coord
    real :: max_lcoord

    character(len=OPTION_PATH_LEN) :: name
    character(len=OPTION_PATH_LEN) :: particle_name
    integer :: id, name_len, tot_len
    integer :: i, j, k, group_spawn
    integer, dimension(3) :: attribute_size
    integer, dimension(:), allocatable :: remove_particles
    
    add_particles(:) = 0

    !Get ele numbers of adjacent elements
    ele_nums => node_neigh(xfield, node_num)

    allocate(node_numbers(size(ele_nums),xfield%dim+1))
    allocate(node_coord(xfield%dim+1))
    allocate(part_counter(size(group_arrays)))
    part_counter(:) = 0

    !Get node numbers for each adjacent element
    do j = 1,size(ele_nums)
       node_numbers(j,:) = ele_nodes(xfield,ele_nums(j))
    end do

    !Count the total number of particles in each surrounding CV
    do j = 1,size(ele_nums)
       do i = 1,size(group_arrays)
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

    !Check if one group of particles makes up cap_percent% of particles in surrounding CV's
    !only spawn particles from this group if true
    spawn_group = .false.
    if (present(cap_percent)) then
       do i=1,size(group_arrays)
          if ((part_counter(i)*1.0)/(sum(part_counter(:))*1.0)>cap_percent) then
             spawn_group=.true.
             group_spawn=i
          end if
       end do
    end if

    !If adjacent CV's contain particles, clone particles from adjacent
    !CV's into this CV

    do j = 1,size(group_arrays)
       temp_part => particle_lists(group_arrays(j))%last
       if (associated(temp_part)) then
          if (spawn_group.eqv..true.) then
             if (j/=group_spawn) cycle
          end if

          id = temp_part%id_number
          name_len = len(int2str(id+1))+1 !id length + '_'
          tot_len = len(trim(temp_part%name))
          name = trim(temp_part%name(1:tot_len-name_len))
          attribute_size(1) = size(temp_part%attributes)
          attribute_size(2) = size(temp_part%old_attributes)
          attribute_size(3) = size(temp_part%old_fields)
       end if
       !Duplicate parent particles in surrounding CV's
       do i = 1,size(ele_nums)
          do k = 1,xfield%dim+1
             particle => node_particles(j, node_numbers(i,k))%first
             if (node_numbers(i,k)==node_num) cycle !prevent from duplicating particles spawned in this routine
             do while(associated(particle))
                particle_name = trim(name)//int2str(id+1)
                temp_part => null()
                call allocate(temp_part, size(particle%position), size(particle%local_coords), attribute_size)
                temp_part%name = trim(particle_name)
                temp_part%id_number = id+1
                temp_part%list_id = particle%list_id
                temp_part%element=ele_nums(i)
                !randomly select local coords within the element, ensuring coords are within cv
                max_lcoord=rand()/(1/0.45)+0.51!lcoords for cv range from 0.51<x<1
                call set_spawned_lcoords(max_lcoord, node_coord, node_num, node_numbers(i,:))
                temp_part%local_coords = node_coord

                !Convert newly spawned particle's local coords to global coords
                call local_to_global(xfield, temp_part%local_coords, temp_part%element, temp_part%position)
                temp_part%type = LAGRANGIAN_DETECTOR
                !Copy parent particle attributes
                temp_part%attributes(:) = particle%attributes(:)
                temp_part%old_attributes(:) = particle%old_attributes(:)
                temp_part%old_fields(:) = particle%old_fields(:)

                node_part_count(node_num) = node_part_count(node_num) + 1

                !add particle to relevant particle list
                temp_part%previous => particle_lists(group_arrays(j))%last
                particle_lists(group_arrays(j))%last%next => temp_part
                particle_lists(group_arrays(j))%last => temp_part
                particle_lists(group_arrays(j))%last%next => null()
                particle_lists(group_arrays(j))%length = particle_lists(group_arrays(j))%length + 1

                !add particle to temp particle list
                call temp_insert(temp_part,node_particles(j,node_num))
                
                add_particles(j) = add_particles(j) + 1
                id = id + 1
                particle => particle%temp_next
             end do
          end do
       end do
    end do

    allocate(remove_particles(size(group_arrays)))
    remove_particles(:) = 0
    
    do while(node_part_count(node_num)>max_thresh)
       call delete_particles(node_part_count(node_num), node_particles(:,node_num), group_arrays, remove_particles)
       add_particles(:) = add_particles(:) - remove_particles(:)
    end do
    
    deallocate(remove_particles)
    deallocate(node_coord)
    deallocate(node_numbers)
    deallocate(part_counter)

  end subroutine spawn_zero_particles

  subroutine set_spawned_lcoords(max_lcoord, node_coord, node_num, node_numbers)
    !Subroutine to randomly set spawned particle local coordinates based off
    !the maximum local coordinate given

    real, intent(inout) :: max_lcoord
    real, dimension(:), intent(inout) :: node_coord
    integer, intent(in) :: node_num
    integer, dimension(:), intent(in) :: node_numbers

    real :: rand_lcoord_2, rand_lcoord_3
    integer :: i
    
    if (max_lcoord<0.51) then
       max_lcoord=0.51 !Ensures minimum lcoord is 0.51
    end if
    if (max_lcoord>0.999) then
       max_lcoord=0.999 !Ensures maximum lcoord is 0.999
    end if
    do i = 1,size(node_coord)
       if (node_num==node_numbers(i)) then
          select case(size(node_coord))
          case(2)
             node_coord(:)=1-max_lcoord
             node_coord(i)=max_lcoord
          case(3)
             rand_lcoord_2=(1-max_lcoord)*rand()
             select case(i)
             case(1)
                node_coord(1)=max_lcoord !Will fall in this nodes CV
                node_coord(2)=rand_lcoord_2
                node_coord(3)=1-max_lcoord-rand_lcoord_2
             case(2)
                node_coord(1)=rand_lcoord_2
                node_coord(2)=max_lcoord!Will fall in this nodes CV
                node_coord(3)=1-max_lcoord-rand_lcoord_2
             case(3)
                node_coord(1)=rand_lcoord_2
                node_coord(2)=1-max_lcoord-rand_lcoord_2
                node_coord(3)=max_lcoord!Will fall in this nodes CV
             end select
          case(4)
             rand_lcoord_2=(1-max_lcoord)*rand()
             rand_lcoord_3=(1-max_lcoord-rand_lcoord_2)*rand()
             select case(i)
             case(1)
                node_coord(1)=max_lcoord !Will fall in this nodes CV
                node_coord(2)=rand_lcoord_2
                node_coord(3)=rand_lcoord_3
                node_coord(4)=1-max_lcoord-rand_lcoord_2-rand_lcoord_3
             case(2)
                node_coord(1)=rand_lcoord_2
                node_coord(2)=max_lcoord !Will fall in this nodes CV
                node_coord(3)=rand_lcoord_3
                node_coord(4)=1-max_lcoord-rand_lcoord_2-rand_lcoord_3
             case(3)
                node_coord(1)=rand_lcoord_2
                node_coord(2)=rand_lcoord_3
                node_coord(3)=max_lcoord !Will fall in this nodes CV
                node_coord(4)=1-max_lcoord-rand_lcoord_2-rand_lcoord_3
             case(4)
                node_coord(1)=rand_lcoord_2
                node_coord(2)=rand_lcoord_3
                node_coord(3)=1-max_lcoord-rand_lcoord_2-rand_lcoord_3
                node_coord(4)=max_lcoord !Will fall in this nodes CV
             end select
          end select
       end if
    end do

  end subroutine set_spawned_lcoords

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

  subroutine multi_delete_particles(mult, node_part_count, node_particles, group_arrays, summed_particles)
    integer, intent(in) :: mult
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    real, intent(inout) :: node_part_count
    integer, intent(in), dimension(:) :: group_arrays
    integer, dimension(:), intent(inout) :: summed_particles

    integer :: i, power, j
    integer, allocatable, dimension(:) :: remove_particles
    logical :: power_set

    allocate(remove_particles(size(group_arrays)))

    power_set=.false.
    j=0
    do while (power_set.eqv..false.)
       if (mult>=2**j) then
          j=j+1
       else
          power=j-1
          power_set=.true.
       end if
    end do
    
    do i = 1,power
       call delete_particles(node_part_count, node_particles, group_arrays, remove_particles)
       summed_particles=summed_particles-remove_particles
    end do

    deallocate(remove_particles)

  end subroutine multi_delete_particles

  subroutine delete_particles(node_part_count, node_particles, group_arrays, remove_particles)
    !Subroutine to calculate the number of particles in each control volume, and delete
    !particles within a control volume if a maximum number of particles is exceeded

    use particles, only: particle_lists
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    real, intent(inout) :: node_part_count
    integer, intent(in), dimension(:) :: group_arrays
    integer, dimension(:), intent(inout) :: remove_particles
    
    type(detector_type), pointer :: particle
    type(detector_type), pointer :: temp_part

    integer :: j
    real :: random

    remove_particles(:) = 0

    !loop over all particles in each particle group, flip a coin, if coin is heads delete the particle
    do j = 1,size(group_arrays)
       particle =>node_particles(j)%first
       do while(associated(particle))
          random=rand()
          if (random>0.5) then !Delete the particle
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
             node_part_count = node_part_count - 1
             call deallocate(particle)
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
