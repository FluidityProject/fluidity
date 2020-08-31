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
  use parallel_tools, only: getnprocs, allsum, getprocno
  use elements, only: eval_shape
  use fields_base, only: ele_val, ele_loc, node_val
  use parallel_fields, only: node_owned
  use fields_calculations, only: dot_product
  use fields
  use profiler
  use state_module
  use halos
  use detector_data_types
  use pickers
  use field_options
  use detector_tools, only: temp_insert, insert, allocate, deallocate, temp_deallocate
  use particles, only : get_particle_arrays, initialise_constant_particle_attributes, &
       & update_particle_attributes_and_fields, get_particle_arrays
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
          if (option_count(trim(subgroup_path) // "/attributes/scalar_attribute/constant") + &
               & option_count(trim(subgroup_path) // "/attributes/vector_attribute/constant") + &
               & option_count(trim(subgroup_path) // "/attributes/tensor_attribute/constant").gt.0) then
             call initialise_constant_particle_attributes(state, subgroup_path, particle_lists(list_counter))
          end if
          list_counter = list_counter + 1
       end do
    end do

    !Check if MVF field is generated from particles
    particle_materials = option_count("material_phase/scalar_field::MaterialVolumeFraction/diagnostic/algorithm::from_particles")
    if (particle_materials.gt.0) then
       !Initialise MaterialVolumeFraction fields dependent on particles
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
    !subroutine to initialise particle attributes and diagnostic fields
    !dependent on particles
    type(state_type), dimension(:), intent(inout) :: state

    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path, name
    type(scalar_field), pointer :: s_field
    real :: current_time, dt
    integer :: i, k
    integer :: dim, particle_groups, list_counter
    type(detector_type), pointer :: particle
    logical :: from_file

    ewrite(2,*) "In initialise_particle_diagnostics"
    !Check if there are particles
    particle_groups = option_count("/particles/particle_group")

    if (particle_groups==0) return

    !Allocate parameters
    call get_option("/timestepping/current_time", current_time)
    call get_option("/timestepping/timestep", dt)

    !Initialise particle attributes
    call update_particle_attributes_and_fields(state, current_time, dt)

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

  end subroutine initialise_particle_diagnostics

  subroutine update_particle_diagnostics(state, time, dt)
    !!Routine to loop over particle arrays and update particle attributes
    !!and diagnostic fields which depend on particles
    use particles, only: particle_lists
    type(state_type), dimension(:), intent(inout) :: state
    real, intent(in) :: time
    real, intent(in) :: dt
    type(scalar_field), pointer :: s_field
    type(detector_type), pointer :: particle
    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path

    integer :: i, k
    integer :: dim, particle_groups, list_counter, particle_materials

    !Check whether there are any particles.
    particle_groups = option_count("/particles/particle_group")

    if (particle_groups==0) return

    !Update particle attributes
    call update_particle_attributes_and_fields(state, time, dt)

    call profiler_tic("update_particle_diagnostic_fields")
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

    call profiler_toc("update_particle_diagnostic_fields")

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

    case("from_particles")!Case to calculate diagnostic field from particle ratio method
        call calculate_ratio_from_particles(states, state_index, s_field)

    case("number_of_particles")!Case to calculate diagnostic field from the number of particles present
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

    !!Calculate s_field using the number of particles present in a control volume

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

    logical :: have_attribute_caps, have_radius, copy_parents
    integer :: min_thresh, max_thresh
    real :: radius, cap_percent

    ewrite(2,*) "In particle_cv_check"

    call profiler_tic("particles_spawn_delete")
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

          !Check option for what algorithm particle_attributes should be calculated with
          if (have_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/spawn_attributes/copy_parent_attributes")) then
             copy_parents=.true.
          else
             copy_parents=.false.
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
          !call routine with optional cap parameters when relevant
          if (have_attribute_caps.eqv..true.) then
             call spawn_delete_particles(states, mesh, group_arrays, max_thresh, min_thresh, have_radius, radius, copy_parents, cap_percent)
          else
             call spawn_delete_particles(states, mesh, group_arrays, max_thresh, min_thresh, have_radius, radius, copy_parents)
          end if
          deallocate(group_arrays)
       end if
    end do
    call profiler_toc("particles_spawn_delete")

  end subroutine particle_cv_check

  subroutine spawn_delete_particles(states, mesh, group_arrays, max_thresh, min_thresh, have_radius, radius, copy_parents, cap_percent)
    use particles, only: particle_lists
    !Subroutine to calculate the number of particles in each control volume, and call spawning
    !or deleting routines if threshold limits are broken
    !> Model state structure
    type(state_type), dimension(:), target, intent(in) :: states
    !> Model mesh
    type (mesh_type), intent(in) :: mesh
    !> Number of particle arrays present
    integer, dimension(:), intent(in) :: group_arrays
    !> Minimum and maximum control volume particle thresholds
    integer, intent(in) :: min_thresh, max_thresh
    !> Parameter to determine spawning scheme used
    logical, intent(in) :: have_radius
    !> Parameter to determine particle spawn location
    real, intent(in) :: radius
    !> Parameter to determine spawned particle attribute values
    logical, intent(in) :: copy_parents
    !> Parameter to determine if spawning/deleting will be capped per group 
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

    ewrite(2,*) "In spawn_delete_particles"

    !Get initial fields and number of processors
    xfield=>extract_vector_field(states(1), "Coordinate")
    nprocs = getnprocs()
    if (nprocs>1) then
       halo => mesh%halos(2)
    end if
    !Get number of nodes on the mesh
    allocate(node_particles(size(group_arrays),node_count(mesh)))
    allocate(node_part_count(node_count(mesh)))
    node_part_count(:) = 0

    !Loop over particle arrays and sum particles per node
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

    !Loop over all nodes and ensure thresholds are not broken. Check minimum threshold and spawn particles, then check maximum threshold and delete particles
    
    !Loop over all nodes
    do i = 1,node_count(mesh)
       !Count number of particles per node and ensure minimum threshold is not broken
       if (node_part_count(i)<min_thresh) then
          !Spawn particles if threshold is broken (but only for nodes with particles in CV)
          if (node_part_count(i)>0) then
             if (node_part_count(i)<min_thresh/2) then
                mult = nint(min_thresh/node_part_count(i)*1.0)
                call multi_spawn_particles(temp_part_count(i), node_particles(:,i), group_arrays, xfield, summed_particles, i, mult, have_radius, radius, copy_parents, cap_percent)
             else
                call spawn_particles(temp_part_count(i), node_particles(:,i), group_arrays, xfield, add_particles, i, have_radius, radius, copy_parents, cap_percent)
                summed_particles=summed_particles+add_particles
             end if
          end if
       end if
    end do
    !Spawn on nodes with 0 particles in CV
    do i = 1,node_count(mesh)
       if (node_part_count(i)<min_thresh) then
          if (node_part_count(i)==0) then
             if (node_owned(mesh, i)) then
                call spawn_zero_particles(temp_part_count(:), node_particles(:,:), group_arrays, xfield, add_particles, i, max_thresh, copy_parents)
                summed_particles=summed_particles+add_particles
             end if
          end if
       end if
    end do
    
    !Update halos
    if (nprocs>1) then
       call halo_accumulate(halo, temp_part_count)
    end if
    !Update node particle counter and reset dummy counter
    node_part_count(:) = node_part_count(:) + temp_part_count(:)
    temp_part_count(:) = 0

    !Loop over all nodes
    do i = 1,node_count(mesh)
       !Count number of particles per node and ensure maximum threshold is not broken
       if (node_part_count(i)>max_thresh) then
          !Delete particles if threshold is broken
          if (node_part_count(i)>2*max_thresh) then
             mult = nint(node_part_count(i)*1.0/max_thresh)
             call multi_delete_particles(mult, temp_part_count(i), node_particles(:,i), group_arrays, summed_particles, cap_percent)
          else
             call delete_particles(temp_part_count(i), node_particles(:,i), group_arrays, remove_particles, cap_percent)
             summed_particles=summed_particles-remove_particles
          end if
       end if
    end do

    !Update halos
    if (nprocs>1) then
       call halo_accumulate(halo, temp_part_count)
    end if
    !Update node particle counter
    node_part_count(:) = node_part_count(:) + temp_part_count(:)

    !Update particle_list parameters
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

  subroutine multi_spawn_particles(node_part_count, node_particles, group_arrays, xfield, summed_particles, node_num, mult, have_radius, radius, copy_parents, cap_percent)
    !Subroutine to call spawn particles multiple times based on the mult factor
    !calculated from the current number of particles and the minimum threshold
    !> Linked list of particles which exist on this node
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    !> Array counting particles on all nodes
    real, intent(inout) :: node_part_count
    !> Number of particle arrays present
    integer, intent(in), dimension(:) :: group_arrays
    !> Current node number we are working on
    integer, intent(in) :: node_num
    !> Input position field
    type(vector_field), pointer, intent(in) :: xfield
    !> Array to sum spawned particles
    integer, dimension(:), intent(inout) :: summed_particles
    !> Factor to determine number of spawn_particle calls
    integer, intent(in) :: mult
    !> Parameter to determine spawning scheme used
    logical, intent(in) :: have_radius
    !> Parameter to determine particle spawn location
    real, intent(in) :: radius
    !> Parameter to determine spawned particle attribute values
    logical, intent(in) :: copy_parents
    !> Parameter to determine if spawning/deleting will be capped per group 
    real, optional, intent(in) :: cap_percent

    integer :: i, power, j
    !> Array to count number of particles spawned per particle group
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
       call spawn_particles(node_part_count, node_particles, group_arrays, xfield, add_particles, node_num, have_radius, radius, copy_parents, cap_percent)
       summed_particles=summed_particles+add_particles
    end do

    deallocate(add_particles)

  end subroutine multi_spawn_particles

  subroutine spawn_particles(node_part_count, node_particles, group_arrays, xfield, add_particles, node_num, have_radius, radius, copy_parents, cap_percent)
    !Subroutine to spawn particles in a control volume based off the parent particles present

    use particles, only: particle_lists
    !> Linked list of particles which exist on this node
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    !> Array counting particles on all nodes
    real, intent(inout) :: node_part_count
    !> Number of particle arrays present
    integer, intent(in), dimension(:) :: group_arrays
    !> Current node number we are working on
    integer, intent(in) :: node_num
    !> Input position field
    type(vector_field), pointer, intent(in) :: xfield
    !> Array to count number of particles spawned per particle group
    integer, dimension(:), intent(inout) :: add_particles
    !> Parameter to determine spawning scheme used
    logical, intent(in) :: have_radius
    !> Parameter to determine particle spawn location
    real, intent(in) :: radius
    !> Parameter to determine spawned particle attribute values
    logical, intent(in) :: copy_parents
    !> Parameter to determine if spawning/deleting will be capped per group 
    real, optional, intent(in) :: cap_percent

    !> Dummy particles
    type(detector_type), pointer :: particle
    type(detector_type), pointer :: temp_part

    real, dimension(:), allocatable :: node_coord
    real, dimension(:), allocatable :: weighted_attributes, weighted_old_attributes, weighted_old_fields

    character(len=OPTION_PATH_LEN) :: name
    character(len=OPTION_PATH_LEN) :: particle_name
    integer, dimension(:), allocatable :: node_numbers, ele_num_part
    integer, dimension(:), pointer :: ele_nums
    integer :: id, name_len, tot_len, group_spawn, ele_spawn, proc_num
    integer :: j, i, k, l, m
    integer, dimension(3) :: attribute_size
    logical :: spawn_group, coords_set, rand_set
    real :: max_lcoord, rand_lcoord, sum_coords
    real, dimension(:), allocatable :: rand_lcoords

    proc_num = getprocno()

    !Check the ratio of particles present per group and determine which groups spawn
    spawn_group = .true.
    if (present(cap_percent)) then
       spawn_group = .false.
       do i=1,size(group_arrays)
          if ((((node_particles(i)%length*1.0)/sum(node_particles(:)%length)*1.0)*100)>=cap_percent) then
             spawn_group = .true.
             group_spawn = i
          end if
       end do
    end if

    add_particles(:) = 0
    !Return if no groups spawn
    if (spawn_group.eqv..false.) return

    !Count the number of particles per surrounding element for weighted spawning if
    !not spawning based on radius from parent
    if (have_radius.eqv..false.) then
       ele_nums => node_neigh(xfield, node_num)
       allocate(ele_num_part(size(ele_nums)))
       ele_num_part(:) = 0
       do j = 1,size(group_arrays)
          particle => node_particles(j)%first
          do while(associated(particle))
             do i = 1,size(ele_nums)
                if (ele_nums(i) == particle%element) then
                   ele_num_part(i) = ele_num_part(i) + 1
                end if
             end do
             particle => particle%temp_next
          end do
       end do
    end if

    !Loop over particle groups for spawning
    do j = 1,size(group_arrays)
       if (present(cap_percent)) then
          if (j/=group_spawn) cycle
       end if
       !Calculated weighted attribute parameters based on existing particles if
       !not copying parent attributes when spawning
       if (copy_parents.eqv..false.) then
          particle => node_particles(j)%first
          if (associated(particle)) then
             allocate(weighted_attributes(size(particle%attributes)))
             allocate(weighted_old_attributes(size(particle%old_attributes)))
             allocate(weighted_old_fields(size(particle%old_fields)))
             weighted_attributes(:) = 0
             weighted_old_attributes(:) = 0
             weighted_old_fields(:) = 0
          end if
          do while(associated(particle))
             weighted_attributes = weighted_attributes + particle%attributes
             weighted_old_attributes = weighted_old_attributes + particle%old_attributes
             weighted_old_fields = weighted_old_fields + particle%old_fields
             particle => particle%temp_next
          end do
          if (allocated(weighted_attributes)) then
             weighted_attributes = weighted_attributes/node_particles(j)%length
             weighted_old_attributes = weighted_old_attributes/node_particles(j)%length
             weighted_old_fields = weighted_old_fields/node_particles(j)%length
          end if
       end if

       !Determine id, size and shape of spawned particle parameters
       temp_part => particle_lists(group_arrays(j))%last
       if (associated(temp_part)) then

          id = particle_lists(group_arrays(j))%proc_part_count
          name_len = len(int2str(id+1))+1 !id length + '_'
          tot_len = len(trim(temp_part%name))
          name = trim(temp_part%name(1:tot_len-name_len))
          attribute_size(1) = size(temp_part%attributes)
          attribute_size(2) = size(temp_part%old_attributes)
          attribute_size(3) = size(temp_part%old_fields)
          allocate(node_coord(size(temp_part%local_coords)))
          allocate(node_numbers(size(temp_part%local_coords)))
       end if

       !Spawn particles
       particle => node_particles(j)%first
       do while(associated(particle))
          particle_name = trim(name)//int2str(id+1)
          temp_part => null()
          call allocate(temp_part, size(particle%position), size(particle%local_coords), attribute_size)

          temp_part%name = trim(particle_name)
          temp_part%id_number = id+1
          temp_part%list_id = particle%list_id
          temp_part%proc_id = proc_num

          !Check if particles are spawning within a radius around their parent, or randomly within the CV
          if (have_radius.eqv..true.) then
             !Spawn particles within radius of parent particle
             temp_part%element=particle%element
             node_numbers(:) = ele_nodes(xfield, temp_part%element)
             allocate(rand_lcoords(size(node_coord)))
             coords_set=.false.
             !dummy counter l to prevent infinite looping
             l = 0
             !randomly select local coords within input radius around parent particle
             !ensuring new particle falls within cv
             do while (coords_set.eqv..false.)
                rand_lcoord=0
                rand_set=.false.
                do while (rand_set.eqv..false.)
                   k = 1 + FLOOR(size(node_coord)*rand()) !randomly select a local_coord to be sum of other coords
                   do i = 1,size(node_coord)
                      if (node_numbers(i)==node_num) m=i
                      if (i==k) cycle
                      rand_lcoords(i) = (rand()-0.5)*2*(radius/(size(node_coord)-1)) !perturb each free local_coord by +/- radius/dim
                      rand_lcoord = rand_lcoord + rand_lcoords(i)!sum local coord perturbations
                   end do
                   rand_lcoords(k) = -rand_lcoord !calculate final coord
                   if (MAXVAL(rand_lcoords)>((radius/(size(node_coord)-1))*0.8)) then !Enforce the new particle to be at least rad/dim * 80% away
                      rand_set=.true.
                   end if
                end do
                node_coord(:) = particle%local_coords(:) + rand_lcoords(:) !add perturbation to spawned local_coordinates
                if (sum(node_coord)/=1) then !ensure sum of local_coordinates = 1
                   sum_coords=0
                   do i = 1,size(node_coord)
                      if (i==k) cycle
                      sum_coords=sum_coords+node_coord(i)
                   end do
                   node_coord(k)=1-sum_coords
                end if
                coords_set=.true.
                l = l + 1
                if ((maxloc(node_coord, 1)/=m).or.node_coord(m)>1.0) then !Not within CV or element, try again
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
          else !Particles will spawn randomly in the CV, prioritizing
               !elements with fewer particles
               !randomly select adjacent element to node
               !Get ele numbers of adjacent elements
             ele_spawn = MINLOC(ele_num_part, DIM=1)
             !ele_spawn=floor(rand()*size(ele_nums)+1)
             temp_part%element = ele_nums(ele_spawn)
             node_numbers(:) = ele_nodes(xfield, ele_nums(ele_spawn))
             max_lcoord = rand()/(1/0.49)+0.51!lcoords for cv range from 0.51<x<1
             call set_spawned_lcoords(max_lcoord, node_coord, node_num, node_numbers)
             temp_part%local_coords = node_coord
             ele_num_part(ele_spawn) = ele_num_part(ele_spawn) + 1
          end if

          !Convert local particle coordinates to global coordinates
          call local_to_global(xfield, temp_part%local_coords, temp_part%element, temp_part%position)
          !Copy parent particle attributes
          if (copy_parents.eqv..true.) then
             !Copy parent particle attributes
             temp_part%attributes(:) = particle%attributes(:)
             temp_part%old_attributes(:) = particle%old_attributes(:)
             temp_part%old_fields(:) = particle%old_fields(:)
          else
             temp_part%attributes(:) = weighted_attributes(:)
             temp_part%old_attributes(:) = weighted_old_attributes(:)
             temp_part%old_fields(:) = weighted_old_fields(:)
          end if
          node_part_count = node_part_count + 1

          !add particle to relevant particle list
          temp_part%previous => particle_lists(group_arrays(j))%last
          particle_lists(group_arrays(j))%last%next => temp_part
          particle_lists(group_arrays(j))%last => temp_part
          particle_lists(group_arrays(j))%last%next => null()
          particle_lists(group_arrays(j))%length = particle_lists(group_arrays(j))%length + 1
          add_particles(j) = add_particles(j) + 1
          id = id + 1
          particle_lists(group_arrays(j))%proc_part_count = particle_lists(group_arrays(j))%proc_part_count + 1

          particle => particle%temp_next

       end do
       if (allocated(node_coord)) then
          deallocate(node_coord)
          deallocate(node_numbers)
       end if
       if (allocated(weighted_attributes)) then
          deallocate(weighted_attributes)
          deallocate(weighted_old_attributes)
          deallocate(weighted_old_fields)
       end if
    end do
    if (allocated(ele_num_part)) deallocate(ele_num_part)

  end subroutine spawn_particles

  subroutine spawn_zero_particles(node_part_count, node_particles, group_arrays, xfield, add_particles, node_num, max_thresh, copy_parents)
    !Subroutine to spawn particles in a control volume with 0 'parent' particles in the CV

    use particles, only: particle_lists
    !> Linked list of particles which exist on this node
    type(detector_linked_list), intent(inout), dimension(:,:) :: node_particles
    !> Array counting particles on all nodes
    real, intent(inout), dimension(:) :: node_part_count
    !> Number of particle arrays present
    integer, intent(in), dimension(:) :: group_arrays
    !> Input position field
    type(vector_field), pointer, intent(in) :: xfield
    !> Array to count number of particles spawned per particle group
    integer, dimension(:), intent(inout) :: add_particles
    !> Current node number we are working on
    integer, intent(in) :: node_num
    !> Maximum number of particles allowed per control volume
    integer, intent(in) :: max_thresh
    !> Parameter to determine spawned particle attribute values
    logical, intent(in) :: copy_parents

    !> Array for element numbers associated with a node
    integer, dimension(:), pointer :: ele_nums
    !> Array for node numbers associated with an element
    integer, allocatable, dimension(:,:) :: node_numbers
    integer :: part_total

    !> Dummy particles
    type(detector_type), pointer :: particle
    type(detector_type), pointer :: temp_part

    !> Arrays for coordinates of nodes and distance from nodes
    real, dimension(:), allocatable :: node_coord, node_loc, distance
    real, dimension(:,:,:), allocatable :: ele_val
    real :: max_lcoord, ratio

    character(len=OPTION_PATH_LEN) :: name
    character(len=OPTION_PATH_LEN) :: particle_name
    integer :: id, name_len, tot_len
    integer :: i, j, k, l, m, proc_num
    integer, dimension(3) :: attribute_size
    integer, dimension(:), allocatable :: remove_particles
    integer, dimension(:), allocatable :: length_group

    !Variable arrays to calculated weighted attribute values when not
    !copying attributes from parents
    type VarSizedArray
       real, allocatable :: col(:)
    end type VarSizedArray
    type(VarSizedArray), allocatable :: weighted_attributes(:)
    type(VarSizedArray), allocatable :: weighted_old_attributes(:)
    type(VarSizedArray), allocatable :: weighted_old_fields(:)

    ewrite(2,*) "In spawn_zero_particles"

    proc_num = getprocno()
    add_particles(:) = 0

    !Get ele numbers of adjacent elements
    ele_nums => node_neigh(xfield, node_num)

    allocate(node_numbers(size(ele_nums),xfield%dim+1))
    allocate(node_coord(xfield%dim+1))
    allocate(node_loc(xfield%dim))
    allocate(ele_val(size(group_arrays),size(ele_nums),xfield%dim+1))

    !Get node numbers for each adjacent element
    do j = 1,size(ele_nums)
       node_numbers(j,:) = ele_nodes(xfield,ele_nums(j))
    end do

    !Initialise weighted attribute arrays based on particle attributes in
    !surrounding elements if not copying parent attributes
    if (copy_parents.eqv..false.) then
       allocate(weighted_attributes(size(group_arrays)))
       allocate(weighted_old_attributes(size(group_arrays)))
       allocate(weighted_old_fields(size(group_arrays)))
       allocate(length_group(size(group_arrays)))
       length_group(:) = 0
       do i = 1,size(group_arrays)
          temp_part => particle_lists(group_arrays(i))%last
          if (associated(temp_part)) then
             allocate(weighted_attributes(i)%col(size(temp_part%attributes)))
             allocate(weighted_old_attributes(i)%col(size(temp_part%old_attributes)))
             allocate(weighted_old_fields(i)%col(size(temp_part%old_fields)))
             weighted_attributes(i)%col(:)=0
             weighted_old_attributes(i)%col(:)=0
             weighted_old_fields(i)%col(:)=0
          end if
       end do
    end if
    node_loc = node_val(xfield,node_num)
    ele_val(:,:,:) = 0
    part_total = 0
    !Loop over elements adjacent to control volume
    do j = 1,size(ele_nums)
       do k = 1,xfield%dim+1 !loop over each node of the element
          if (node_numbers(j,k)==node_num) cycle !cycle if node is the node from our control volume
          do i = 1,size(group_arrays) !loop over particle grouos
             temp_part => node_particles(i,node_numbers(j,k))%first
             allocate(distance(node_particles(i,node_numbers(j,k))%length))
             distance(:)=0
             !loop over particles in this group and determine distance from the CV node
             do l = 1,node_particles(i,node_numbers(j,k))%length
                do m = 1,xfield%dim
                   distance(l) = distance(l) + abs(node_loc(m)-temp_part%position(m))**2
                end do
                distance(l) = SQRT(distance(l))
                ele_val(i,j,k) = ele_val(i,j,k) + 1/distance(l)**2!store distance values for weighting
                if (copy_parents.eqv..false.) then !copy attributes for weighting
                   weighted_attributes(i)%col(:) = weighted_attributes(i)%col(:) + temp_part%attributes(:)
                   weighted_old_attributes(i)%col(:) = weighted_old_attributes(i)%col(:) + temp_part%old_attributes(:)
                   weighted_old_fields(i)%col(:) = weighted_old_fields(i)%col(:) + temp_part%old_fields(:)
                end if
                temp_part => temp_part%temp_next
             end do
             part_total = part_total + node_particles(i,node_numbers(j,k))%length
             if (copy_parents.eqv..false.) then !determine the number of particles in each group being weighted
                length_group(i) = length_group(i) + node_particles(i,node_numbers(j,k))%length
             end if
             deallocate(distance)
          end do
       end do
    end do

    !Return if no particles in surrounding control volumes
    if (part_total==0) then
       ewrite(2,*) "There are no particles present in adjacent CV's"
       return
    end if

    !If adjacent CV's contain particles, clone weighted particles
    !from adjacent CV's into this CV
    do i = 1,size(group_arrays)!Loop over particle groups
       temp_part => particle_lists(group_arrays(i))%last
       if (associated(temp_part)) then

          id = particle_lists(group_arrays(i))%proc_part_count
          name_len = len(int2str(id+1))+1 !id length + '_'
          tot_len = len(trim(temp_part%name))
          name = trim(temp_part%name(1:tot_len-name_len))
          attribute_size(1) = size(temp_part%attributes)
          attribute_size(2) = size(temp_part%old_attributes)
          attribute_size(3) = size(temp_part%old_fields)
       end if
       if (copy_parents.eqv..false.) then !Weight surrounding attributes if not copying from parent
          weighted_attributes(i)%col(:) = weighted_attributes(i)%col(:)/length_group(i)
          weighted_old_attributes(i)%col(:) = weighted_old_attributes(i)%col(:)/length_group(i)
          weighted_old_fields(i)%col(:) = weighted_old_fields(i)%col(:)/length_group(i)
       end if
       !Duplicate parent particles infrom surrounding CV's weighting based on distance
       do j = 1,size(ele_nums)!loop over surrounding elements
          do k = 1,xfield%dim+1!loop over nodes attached to element
             particle => node_particles(i, node_numbers(j,k))%first
             if (node_numbers(j,k)==node_num) cycle !prevent from duplicating particles spawned in this routine
             if (.not. associated(particle)) cycle
             ratio = ele_val(i,j,k)/SUM(ele_val(:,j,k))!determine ratio of particles in given CV to all surrounding CV's
             if (ISNAN(ratio)) ratio=0
             !Spawn a number of particle from this CV based on given parameters:
             !maximum particle threshold/4 * 1/number of surrounding elements * ratio
             !of particles in given CV to all surrounding CV's
             do l = 1,NINT((max_thresh/4.0)*(1.0/size(ele_nums))*ratio)
                particle_name = trim(name)//int2str(id+1)
                temp_part => null()
                call allocate(temp_part, size(particle%position), size(particle%local_coords), attribute_size)
                temp_part%name = trim(particle_name)
                temp_part%id_number = id+1
                temp_part%list_id = particle%list_id
                temp_part%proc_id = proc_num
                temp_part%element=ele_nums(j)
                !randomly select local coords within the element, ensuring coords are within cv
                max_lcoord=rand()/(1/0.49)+0.51!lcoords for cv range from 0.51<x<1
                node_coord(:)=1
                node_coord(k)=0
                do while (MINLOC(node_coord,DIM=1)==k)
                   call set_spawned_lcoords(max_lcoord, node_coord, node_num, node_numbers(j,:))
                end do
                temp_part%local_coords = node_coord

                !Convert newly spawned particle's local coords to global coords
                call local_to_global(xfield, temp_part%local_coords, temp_part%element, temp_part%position)
                if (copy_parents.eqv..true.) then
                   !Copy parent particle attributes
                   temp_part%attributes(:) = particle%attributes(:)
                   temp_part%old_attributes(:) = particle%old_attributes(:)
                   temp_part%old_fields(:) = particle%old_fields(:)
                else
                   temp_part%attributes(:) = weighted_attributes(i)%col(:)
                   temp_part%old_attributes(:) = weighted_old_attributes(i)%col(:)
                   temp_part%old_fields(:) = weighted_old_fields(i)%col(:)
                end if

                node_part_count(node_num) = node_part_count(node_num) + 1

                !add particle to relevant particle list
                temp_part%previous => particle_lists(group_arrays(i))%last
                particle_lists(group_arrays(i))%last%next => temp_part
                particle_lists(group_arrays(i))%last => temp_part
                particle_lists(group_arrays(i))%last%next => null()
                particle_lists(group_arrays(i))%length = particle_lists(group_arrays(i))%length + 1

                !add particle to temp particle list
                call temp_insert(temp_part,node_particles(i,node_num))

                add_particles(i) = add_particles(i) + 1
                id = id + 1
                particle_lists(group_arrays(i))%proc_part_count = particle_lists(group_arrays(i))%proc_part_count + 1
             end do
          end do
       end do
    end do

    allocate(remove_particles(size(group_arrays)))
    remove_particles(:) = 0

    !Call delete_particles if number of particles now exceeds the maximum threshold
    do while(node_part_count(node_num)>max_thresh)
       call delete_particles(node_part_count(node_num), node_particles(:,node_num), group_arrays, remove_particles)
       add_particles(:) = add_particles(:) - remove_particles(:)
    end do

    deallocate(remove_particles)
    deallocate(node_coord)
    deallocate(node_numbers)
    deallocate(node_loc)
    deallocate(ele_val)
    if (allocated(weighted_attributes)) then
       deallocate(weighted_attributes)
       deallocate(weighted_old_attributes)
       deallocate(weighted_old_fields)
       deallocate(length_group)
    end if

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

  subroutine multi_delete_particles(mult, node_part_count, node_particles, group_arrays, summed_particles, cap_percent)
    !Subroutine to call delete particles multiple times based on the mult factor
    !calculated from the current number of particles and the maximum threshold
    !> Factor to determine number of delete_particle calls
    integer, intent(in) :: mult
    !> Linked list of particles which exist on this node
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    !> Array counting particles on all nodes
    real, intent(inout) :: node_part_count
    !> Number of particle arrays present
    integer, intent(in), dimension(:) :: group_arrays
    !> Array to sum deleted particles
    integer, dimension(:), intent(inout) :: summed_particles
    !> Parameter to determine if spawning/deleting will be capped per group 
    real, optional, intent(in) :: cap_percent

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
       call delete_particles(node_part_count, node_particles, group_arrays, remove_particles, cap_percent)
       summed_particles=summed_particles-remove_particles
    end do

    deallocate(remove_particles)

  end subroutine multi_delete_particles

  subroutine delete_particles(node_part_count, node_particles, group_arrays, remove_particles, cap_percent)
    !Subroutine to delete particles in a control volume based off the number of particles present

    use particles, only: particle_lists
    !> Linked list of particles which exist on this node
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    !> Array counting particles on all nodes
    real, intent(inout) :: node_part_count
    !> Number of particle arrays present
    integer, intent(in), dimension(:) :: group_arrays
    !> Array to sum number of deleted particles per particle group
    integer, dimension(:), intent(inout) :: remove_particles
    !> Parameter to determine if spawning/deleting will be capped per group
    real, optional, intent(in) :: cap_percent

    type(detector_type), pointer :: particle
    type(detector_type), pointer :: temp_part

    logical :: delete_group
    integer :: j
    real :: random

    !Check ratio of each particle group present and which groups will be deleted
    delete_group = .true.
    if (present(cap_percent)) then
       delete_group = .false.
       do j=1,size(group_arrays)
          if ((((node_particles(j)%length*1.0)/sum(node_particles(:)%length)*1.0)*100)>=cap_percent) then
             delete_group = .true.
          end if
       end do
    end if

    remove_particles(:) = 0

    !return if no group is being deleted
    if (delete_group.eqv..false.) return

    !loop over all particles in each particle group, flip a coin, if coin is heads (r>0.5) delete the particle
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
