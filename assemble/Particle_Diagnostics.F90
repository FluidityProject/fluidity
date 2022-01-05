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
  use fields_base, only: ele_val, ele_loc, node_val, eval_field
  use parallel_fields, only: node_owned
  use fields_calculations, only: dot_product
  use fields
  use profiler
  use state_module
  use halos
  use detector_data_types
  use pickers
  use field_options
  use detector_tools, only: temp_list_insert, insert, allocate, deallocate, temp_list_deallocate, &
       & remove, temp_list_remove
  use particles, only : get_particle_arrays, initialise_constant_particle_attributes, &
       & particle_lists
  use multimaterial_module, only: calculate_sum_material_volume_fractions

  implicit none

  private

  public :: calculate_particle_material_fields, calculate_diagnostic_fields_from_particles, &
       & initialise_constant_particle_diagnostics, particle_cv_check

  contains

  subroutine initialise_constant_particle_diagnostics(state)
    !subroutine to initialise constant particle attributes and
    !MVF fields if 'from_particles'

    type(state_type), dimension(:), intent(inout) :: state

    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path
    type(scalar_field), pointer :: s_field
    integer :: i, j
    integer :: particle_groups, list_counter, particle_materials
    integer, dimension(:), allocatable :: particle_arrays

    type(detector_type), pointer :: particle

    ewrite(1,*) "In initialise_constant_particle_diagnostics"

    !Check if there are particles
    particle_groups = option_count("/particles/particle_group")

    if (particle_groups==0) return

    !Set up particle_lists
    allocate(particle_arrays(particle_groups))
    particle_arrays = 0
    do i = 1,particle_groups
       particle_arrays(i) = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
    end do

    !Initialise constant particle attributes
    list_counter=1
    do i = 1,particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"
       do j = 1,particle_arrays(i)
          particle => particle_lists(list_counter)%first
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(j-1)//"]"
          if (option_count(trim(subgroup_path) // "/attributes/scalar_attribute/constant") + &
               & option_count(trim(subgroup_path) // "/attributes/vector_attribute/constant") + &
               & option_count(trim(subgroup_path) // "/attributes/tensor_attribute/constant")>0) then
             call initialise_constant_particle_attributes(state, subgroup_path, particle_lists(list_counter))
          end if
          list_counter = list_counter + 1
       end do
    end do

    deallocate(particle_arrays)

  end subroutine initialise_constant_particle_diagnostics

  subroutine calculate_diagnostic_fields_from_particles(state)
    !subroutine to calculate diagnostic fields which are dependent on particles
    !MVF fields are not calculated here
    type(state_type), dimension(:), intent(inout) :: state

    character(len = OPTION_PATH_LEN) :: name
    type(scalar_field), pointer :: s_field
    integer :: i, k
    integer :: particle_groups

    ewrite(1,*) "In calculate_diagnostic_fields_from_particles"
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
          call calculate_field_from_particles(state, i, s_field)
        end if
      end do
    end do

    !Initialise fields based on number of particles if present
    do i = 1,size(state)
      do k = 1,scalar_field_count(state(i))
        s_field => extract_scalar_field(state(i),k)
        if (have_option(trim(s_field%option_path)//"/diagnostic/algorithm::number_of_particles")) then
          call calculate_field_from_particles(state, i, s_field)
        end if
      end do
    end do

  end subroutine calculate_diagnostic_fields_from_particles

  subroutine calculate_particle_material_fields(state)
    !subroutine to initialise MVF fields from particles

    type(state_type), dimension(:), intent(inout) :: state

    type(scalar_field), pointer :: s_field
    integer :: i, k, particle_materials

    !Check if MVF field is generated from particles
    particle_materials = option_count("/material_phase/scalar_field::MaterialVolumeFraction/diagnostic/algorithm::from_particles")
    if (particle_materials==0) return

    k = 0

    !Initialise MaterialVolumeFraction fields dependent on particles
    do i = 1,size(state)
       s_field => extract_scalar_field(state(i), "MaterialVolumeFraction")
       if (have_option(trim(s_field%option_path)//"/diagnostic/algorithm::from_particles")) then
          call calculate_field_from_particles(state, i, s_field)
       else if(have_option(trim(s_field%option_path)//"/diagnostic/algorithm::Internal")) then
          k = i
       end if
    end do

    if (k==0) FLAbort("No diagnostic internal algorithm found.")

    !Initialise internal MaterialVolumeFraction field
    s_field => extract_scalar_field(state(k), "MaterialVolumeFraction")
    call calculate_sum_material_volume_fractions(state, s_field)
    call scale(s_field, -1.0)
    call addto(s_field, 1.0)

  end subroutine calculate_particle_material_fields

  subroutine calculate_field_from_particles(states, state_index, s_field)

    !!Calculate s_field using the ratio method, or from the number of particles present
    !!If using the ratio method first determine which particle groups/subgroups/attributes
    !!are being used, then determine the closest node for each particle and store attribute
    !!values. Finally use the ratio method to calculate field values and place on
    !!Diagnostic scalar field

    type(state_type), dimension(:), target, intent(inout) :: states
    integer, intent(in) :: state_index
    type(scalar_field), intent(inout) :: s_field
    type(halo_type), pointer :: halo

    character(len=OPTION_PATH_LEN) :: lgroup, lattribute
    type(vector_field), pointer :: xfield
    type(detector_type), pointer :: particle
    integer :: i, j
    real, allocatable, dimension(:) :: node_values
    real, allocatable, dimension(:) :: node_part_count ! real instead of integer, so we can use halo_accumulate
    integer :: element, node_number
    real, allocatable, dimension(:) :: local_crds
    integer, dimension(:), pointer :: nodes
    integer :: nprocs, att_n
    real :: att_value, ratio_val
    character(len= OPTION_PATH_LEN) :: lmethod
    logical :: from_particles

    integer :: group_attribute
    integer, allocatable, dimension(:) :: group_arrays

    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/name", lmethod, default = "from_particles")

    from_particles = trim(lmethod)=="from_particles"

    call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/particle_group/name", lgroup)
    if (from_particles) then
       if (have_option(trim(complete_field_path(s_field%option_path))// "/algorithm/particle_group/particle_attribute")) then
          call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/particle_group/particle_attribute/name", lattribute)
          att_n=0
          ewrite(2,*) "Calculate diagnostic field from particle group: ", trim(lgroup), ", attribute: ", trim(lattribute)
       else if (have_option(trim(complete_field_path(s_field%option_path))// "/algorithm/particle_group/particle_attribute_array")) then
          call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/particle_group/particle_attribute_array/name", lattribute)
          call get_option(trim(complete_field_path(s_field%option_path))// "/algorithm/particle_group/particle_attribute_array/attribute_index", att_n)
          ewrite(2,*) "Calculate diagnostic field from particle group: ", trim(lgroup), ", attribute_array: ", trim(lattribute), ", attribute_index: ", att_n
       end if
    else
       ewrite(2,*) "Calculate diagnostic field from number of particles on particle group: ", trim(lgroup)
    end if

    xfield=>extract_vector_field(states(1), "Coordinate")

    nodes => ele_nodes(s_field, 1)
    if (xfield%dim+1/=size(nodes)) then
       FLAbort("Can only generate particle diagnostic fields for a P1CV mesh")
    end if

    ! contribution of local particles to non-owned nodes are summed
    ! into the owner in the halo_accumulate calls below
    ! we might be safe to assume we only need to add into halo 1 nodes (as these
    ! are the only ones that make up owned elements), but let's include halo 2 to be sure
    ! Only run this if nprocs > 1
    nprocs = getnprocs()
    if (nprocs>1) then
       halo => s_field%mesh%halos(2)
    end if

    if (.not. allocated(particle_lists)) return !Particles not yet setup, Initial field value will be 0
    !Particles are initialized, call subroutine to get relevant particle arrays and attributes
    if (from_particles) then
       call get_particle_arrays(lgroup, group_arrays, group_attribute, att_n=att_n, lattribute=lattribute)
    else
       call get_particle_arrays(lgroup, group_arrays)
    end if


    !Allocate arrays to store summed attribute values and particle counts at nodes
    if (from_particles) then
       allocate(node_values(node_count(s_field)))
       node_values = 0
    end if
    allocate(node_part_count(node_count(s_field)))
    node_part_count = 0

    if (from_particles .and. have_option(trim(complete_field_path(s_field%option_path))// "/algorithm/interpolation/weighted_distance")) then
       !Calculate node values from attributes with the weighted_distance interpolation algorithm
       !Loop over particle arrays
       do i = 1,size(group_arrays)
          particle => particle_lists(group_arrays(i))%first
          if (.not. associated(particle)) cycle !Only work on arrays if local_particles exist on this processor
          allocate(local_crds(size(particle%local_coords)))
          do while(associated(particle))
             !Get element, local_crds and attribute value of each particle
             element = particle%element
             local_crds = particle%local_coords
             att_value = particle%attributes(group_attribute)

             !Find nodes for specified element
             nodes => ele_nodes(s_field, element)

             !Distribute particle values across nodes of element by weighted distance function
             node_values(nodes) = node_values(nodes) + att_value*local_crds
             node_part_count(nodes) = node_part_count(nodes) + local_crds
             particle => particle%next

          end do
          if (allocated(local_crds)) then
             deallocate(local_crds)
          end if
       end do

    else if (.not. from_particles .or. have_option(trim(complete_field_path(s_field%option_path))// "/algorithm/interpolation/nearest_neighbour")) then
       !Calculate node values from attributes with the nearest neighbour interpolation algorithm
       !or for the case where we are simply counting the number of particles (if not from_particles)

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
             if (from_particles) then
                att_value = particle%attributes(group_attribute)
             end if
             !Find nodes for specified element
             nodes => ele_nodes(s_field, element)
             ! work out nearest node based on max. local coordinate
             node_number = nodes(maxloc(local_crds, dim=1))
             if (from_particles) then
                !Store particle attribute value on closest node
                node_values(node_number) = node_values(node_number) + att_value
             end if
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
       if (from_particles) then
          call halo_accumulate(halo, node_values)
       end if
    end if

    if (from_particles) then
       !Store value on field (if node has at least one particle)
       where (node_part_count/=0)
          !Determine field value from ratio method
          s_field%val = node_values/node_part_count
       elsewhere
          s_field%val = 0
       end where
    else
       !Store number of particles on field
       s_field%val = node_part_count
    end if

    if (from_particles) then
       deallocate(node_values)
    end if
    deallocate(node_part_count)
    deallocate(group_arrays)

    ! all values in owned nodes should now be correct
    ! now we need to make sure the halo is updated accordingly
    if (nprocs>1) then
       call halo_update(s_field)
    end if

  end subroutine calculate_field_from_particles

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
    integer, allocatable :: seed(:)
    integer :: n
    real :: radius, cap_percent

    ewrite(1,*) "In particle_cv_check"

    call profiler_tic("particles_spawn_delete")
    !Check if there are particles
    particle_groups = option_count('/particles/particle_group')

    if (particle_groups==0) return
    call random_seed(size = n)
    allocate(seed(n))
    do k = 1,n
       seed(k) = k*n
    end do
    call random_seed(PUT=seed)
    do k = 1, particle_groups
       if (have_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning")) then
          !Get the mesh particles will be spawned to
          call get_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/mesh/name", mesh_name)
          mesh => extract_mesh(states(1), trim(mesh_name))

          !Set minimum and maximum particle thresholds per control volume
          call get_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/min_cv_threshhold", min_thresh)
          call get_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/max_cv_threshhold", max_thresh)

          !Check option for where particles should be spawned within a control volume.
          have_radius = have_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/spawn_location/radius")
          if (have_radius) call get_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/spawn_location/radius", radius)

          !Check option for what algorithm particle_attributes should be calculated with
          copy_parents = have_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/spawn_attributes/copy_parent_attributes")

          !Check option on whether certain particle subgroup spawning should be capped if one subgroup is dominant
          have_attribute_caps = have_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/subgroup_spawning_caps")
          if (have_attribute_caps) then
             call get_option("/particles/particle_group["//int2str(k-1)//"]/particle_spawning/subgroup_spawning_caps/percentage", cap_percent)
          end if

          !Need to get particle arrays
          call get_option("/particles/particle_group["//int2str(k-1)//"]/name", particle_group)
          call get_particle_arrays(particle_group, group_arrays)
          !call routine with optional cap parameters when relevant
          if (have_attribute_caps) then
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

    !Subroutine to calculate the number of particles in each control volume, and call spawning
    !or deleting routines if threshold limits are broken
    !> Model state structure
    type(state_type), dimension(:), target, intent(in) :: states
    !> Model mesh
    type (mesh_type), intent(in) :: mesh
    !> Indicies in particle_lists present
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

    !Counters to determine the new total number of particles in the relevant linked_list per processor
    integer, allocatable, dimension(:) :: summed_particles, add_particles, remove_particles
    !Temporary particle counter to determine how many particles have been spawned per group per control volume
    real, dimension(:), allocatable :: temp_part_count
    integer :: total_particles

    ewrite(1,*) "In spawn_delete_particles"

    !Get initial fields and number of processors
    xfield=>extract_vector_field(states(1), "Coordinate")
    nprocs = getnprocs()
    if (nprocs>1) then
       halo => mesh%halos(2)
    end if
    !Get number of nodes on the mesh
    !NOTE: this has the potential to be very memory intensive, and we should consider improving this data structure
    allocate(node_particles(size(group_arrays),node_count(mesh)))
    allocate(node_part_count(node_count(mesh)))
    node_part_count = 0

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

          !work out nearest node based on max. local coordinate
          node_number = nodes(maxloc(local_crds, dim=1))

          !Increase particle count for this node by 1
          node_part_count(node_number) = node_part_count(node_number) + 1.0
          call temp_list_insert(particle,node_particles(i,node_number))
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

    summed_particles = 0
    temp_part_count = 0

    !Loop over all nodes and ensure thresholds are not broken. Check minimum threshold and spawn particles, then check maximum threshold and delete particles

    !Loop over all nodes
    do i = 1,node_count(mesh)
       !Count number of particles per node and ensure minimum threshold is not broken
       !Spawn particles if threshold is broken (but only for nodes with particles in CV)
       if (node_part_count(i)>0) then
          if (node_part_count(i)<(min_thresh/2)) then
             mult = nint(min_thresh/node_part_count(i)*1.0)
             call multi_spawn_particles(temp_part_count(i), node_particles(:,i), group_arrays, xfield, summed_particles, i, mult, have_radius, radius, copy_parents, cap_percent)
          else if (node_part_count(i)<min_thresh) then
             call spawn_particles(temp_part_count(i), node_particles(:,i), group_arrays, xfield, add_particles, i, have_radius, radius, copy_parents, cap_percent)
             summed_particles=summed_particles+add_particles
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
    !Add number of spawned particles on each CV to particle count
    node_part_count = node_part_count(:) + temp_part_count(:)

    !Loop over all nodes
    do i = 1,node_count(mesh)
       !Count number of particles per node and ensure maximum threshold is not broken
       !Delete particles if threshold is broken
       if (node_part_count(i)>(2*max_thresh)) then
          mult = nint(node_part_count(i)*1.0/max_thresh)
          call multi_delete_particles(mult, node_particles(:,i), group_arrays, summed_particles, cap_percent)
       else if (node_part_count(i)>max_thresh) then
          call delete_particles(node_particles(:,i), group_arrays, remove_particles, cap_percent=cap_percent)
          summed_particles=summed_particles-remove_particles
       end if
    end do

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
          call temp_list_deallocate(del_node_particles)
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
    !> Indicies in particle_lists present
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
    !Loop to find the largest power of 2 less than or equal to mult factor
    do while (.not. power_set)
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

    !> Linked list of particles which exist on this node
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    !> Array counting particles on all nodes
    real, intent(inout) :: node_part_count
    !> Indicies in particle_lists present
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
    real, dimension(:), allocatable :: average_attributes, average_old_attributes, average_old_fields

    integer, dimension(:), allocatable :: node_numbers, ele_num_part
    integer, dimension(:), pointer :: ele_nums
    integer :: id, group_spawn, ele_spawn, proc_num
    integer :: j, i, k, l, m, dim
    logical :: spawn_group, coords_set, rand_set
    real :: max_lcoord, rand_lcoord, sum_coords, rand_val
    real, dimension(:), allocatable :: rand_lcoords

    proc_num = getprocno()

    !Check the ratio of particles present per group and determine which groups spawn
    spawn_group = .true.
    if (present(cap_percent)) then
       spawn_group = .false.
       do i=1,size(group_arrays)
          !Check if a particle group makes up >= cap_percent of the total particles in the current control volume
          assert(cap_percent>50)
          if ((((node_particles(i)%length*1.0)/sum(node_particles(:)%length)*1.0)*100)>=cap_percent) then
             spawn_group = .true.
             group_spawn = i
          end if
       end do
    end if

    add_particles = 0
    !Return if no groups spawn
    if (.not. spawn_group) return

    !Count the number of particles per surrounding element for weighted spawning if
    !not spawning based on radius from parent
    if (.not. have_radius) then
       ele_nums => node_neigh(xfield, node_num)
       allocate(ele_num_part(size(ele_nums)))
       ele_num_part = 0
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
       !Calculated average attribute parameters based on existing particles if
       !not copying parent attributes when spawning
       if (.not. copy_parents) then
          particle => node_particles(j)%first
          if (associated(particle)) then
             allocate(average_attributes(size(particle%attributes)))
             allocate(average_old_attributes(size(particle%old_attributes)))
             allocate(average_old_fields(size(particle%old_fields)))
             average_attributes = 0
             average_old_attributes = 0
             average_old_fields = 0
          end if
          do while(associated(particle))
             average_attributes = average_attributes + particle%attributes
             average_old_attributes = average_old_attributes + particle%old_attributes
             average_old_fields = average_old_fields + particle%old_fields
             particle => particle%temp_next
          end do
          if (allocated(average_attributes)) then
             average_attributes = average_attributes/node_particles(j)%length
             average_old_attributes = average_old_attributes/node_particles(j)%length
             average_old_fields = average_old_fields/node_particles(j)%length
          end if
       end if

       !Determine id, size and shape of spawned particle parameters
       temp_part => particle_lists(group_arrays(j))%last
       if (associated(temp_part)) then
          id = particle_lists(group_arrays(j))%proc_part_count
          allocate(node_coord(size(temp_part%local_coords)))
          allocate(node_numbers(size(temp_part%local_coords)))
       end if

       dim = size(node_coord)-1

       !Spawn particles
       particle => node_particles(j)%first
       do while(associated(particle))
          temp_part => null()
          call allocate(temp_part, size(particle%position), size(particle%local_coords), particle_lists(group_arrays(j))%total_attributes)

          temp_part%id_number = id+1
          temp_part%list_id = particle%list_id
          temp_part%proc_id = proc_num

          !Check if particles are spawning within a radius around their parent, or randomly within the CV
          if (have_radius) then
             !Spawn particles within radius of parent particle
             temp_part%element=particle%element
             node_numbers(:) = ele_nodes(xfield, temp_part%element)
             allocate(rand_lcoords(size(node_coord)))
             coords_set=.false.
             !dummy counter l to prevent infinite looping
             l = 0
             !randomly select local coords within input radius around parent particle
             !ensuring new particle falls within cv
             do while (.not. coords_set)
                rand_set=.false.
                do while (.not. rand_set)
                   rand_lcoord=0
                   call random_number(rand_val)
                   k = 1 + floor(size(node_coord)*rand_val) !randomly select a local_coord to be sum of other coords
                   do i = 1,size(node_coord)
                      if (node_numbers(i)==node_num) m=i
                      if (i==k) cycle
                      call random_number(rand_val)
                      rand_lcoords(i) = (rand_val-0.5)*2*(radius/dim) !perturb each free local_coord by +/- radius/dim
                      rand_lcoord = rand_lcoord + rand_lcoords(i)!sum local coord perturbations
                   end do
                   rand_lcoords(k) = -rand_lcoord !calculate final coord
                   if (maxval(rand_lcoords)>((radius/dim)*0.8)) then !Enforce the new particle to be at least rad/dim * 80% away
                      rand_set=.true.
                   end if
                end do
                node_coord(:) = particle%local_coords(:) + rand_lcoords(:) !add perturbation to spawned local_coordinates
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
             ele_spawn = minloc(ele_num_part, DIM=1)
             temp_part%element = ele_nums(ele_spawn)
             node_numbers(:) = ele_nodes(xfield, ele_nums(ele_spawn))
             call random_number(rand_val)
             max_lcoord = rand_val/(1/0.49)+0.51!lcoords for cv range from 0.51<x<1
             call set_spawned_lcoords(max_lcoord, node_coord, node_num, node_numbers)
             temp_part%local_coords = node_coord
             ele_num_part(ele_spawn) = ele_num_part(ele_spawn) + 1
          end if

          !Convert local particle coordinates to global coordinates
          temp_part%position = eval_field(temp_part%element, xfield, temp_part%local_coords)
          !Copy parent particle attributes
          if (copy_parents) then
             !Copy parent particle attributes
             temp_part%attributes(:) = particle%attributes(:)
             temp_part%old_attributes(:) = particle%old_attributes(:)
             temp_part%old_fields(:) = particle%old_fields(:)
          else
             temp_part%attributes(:) = average_attributes(:)
             temp_part%old_attributes(:) = average_old_attributes(:)
             temp_part%old_fields(:) = average_old_fields(:)
          end if
          node_part_count = node_part_count + 1

          !add particle to relevant particle list
          call insert(temp_part, particle_lists(group_arrays(j)))
          add_particles(j) = add_particles(j) + 1
          id = id + 1
          particle_lists(group_arrays(j))%proc_part_count = particle_lists(group_arrays(j))%proc_part_count + 1

          particle => particle%temp_next

       end do
       if (allocated(node_coord)) then
          deallocate(node_coord)
          deallocate(node_numbers)
       end if
       if (allocated(average_attributes)) then
          deallocate(average_attributes)
          deallocate(average_old_attributes)
          deallocate(average_old_fields)
       end if
    end do
    if (allocated(ele_num_part)) deallocate(ele_num_part)

  end subroutine spawn_particles

  subroutine spawn_zero_particles(node_part_count, node_particles, group_arrays, xfield, add_particles, node_num, max_thresh, copy_parents)
    !Subroutine to spawn particles in a control volume with 0 'parent' particles in the CV

    !> Linked list of particles which exist on this node
    type(detector_linked_list), intent(inout), dimension(:,:) :: node_particles
    !> Array counting particles on all nodes
    real, intent(inout), dimension(:) :: node_part_count
    !> Indicies in particle_lists present
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
    real :: max_lcoord, ratio, rand_val

    integer :: id
    integer :: i, j, k, l, m, proc_num
    integer, dimension(:), allocatable :: remove_particles
    integer, dimension(:), allocatable :: length_group

    !Variable arrays to calculated average attribute values when not
    !copying attributes from parents
    type VarSizedArray
       real, allocatable :: col(:)
    end type VarSizedArray
    type(VarSizedArray), allocatable :: average_attributes(:)
    type(VarSizedArray), allocatable :: average_old_attributes(:)
    type(VarSizedArray), allocatable :: average_old_fields(:)

    ewrite(1,*) "In spawn_zero_particles"

    proc_num = getprocno()
    add_particles = 0

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

    !Initialise average attribute arrays based on particle attributes in
    !surrounding elements if not copying parent attributes
    if (.not. copy_parents) then
       allocate(average_attributes(size(group_arrays)))
       allocate(average_old_attributes(size(group_arrays)))
       allocate(average_old_fields(size(group_arrays)))
       allocate(length_group(size(group_arrays)))
       length_group = 0
       do i = 1,size(group_arrays)
          temp_part => particle_lists(group_arrays(i))%last
          if (associated(temp_part)) then
             allocate(average_attributes(i)%col(size(temp_part%attributes)))
             allocate(average_old_attributes(i)%col(size(temp_part%old_attributes)))
             allocate(average_old_fields(i)%col(size(temp_part%old_fields)))
             average_attributes(i)%col(:)=0
             average_old_attributes(i)%col(:)=0
             average_old_fields(i)%col(:)=0
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
                if (.not. copy_parents) then !copy attributes to take average
                   average_attributes(i)%col(:) = average_attributes(i)%col(:) + temp_part%attributes(:)
                   average_old_attributes(i)%col(:) = average_old_attributes(i)%col(:) + temp_part%old_attributes(:)
                   average_old_fields(i)%col(:) = average_old_fields(i)%col(:) + temp_part%old_fields(:)
                end if
                temp_part => temp_part%temp_next
             end do
             part_total = part_total + node_particles(i,node_numbers(j,k))%length
             if (.not. copy_parents) then !determine the number of particles in each group being weighted
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
       id = particle_lists(group_arrays(i))%proc_part_count
       if (.not. copy_parents) then !take average of surrounding attributes if not copying from parent
          average_attributes(i)%col(:) = average_attributes(i)%col(:)/length_group(i)
          average_old_attributes(i)%col(:) = average_old_attributes(i)%col(:)/length_group(i)
          average_old_fields(i)%col(:) = average_old_fields(i)%col(:)/length_group(i)
       end if
       !Duplicate parent particles in from surrounding CV's weighting based on distance
       do j = 1,size(ele_nums)!loop over surrounding elements
          do k = 1,xfield%dim+1!loop over nodes attached to element
             particle => node_particles(i, node_numbers(j,k))%first
             if (node_numbers(j,k)==node_num) cycle !prevent from duplicating particles spawned in this routine
             if (.not. associated(particle)) cycle
             if (sum(ele_val(:,j,k))==0) then
                ratio=0
             else
                ratio = ele_val(i,j,k)/sum(ele_val(:,j,k))!determine ratio of particles in given CV to all surrounding CV's
             end if
             !Spawn a number of particle from this CV based on given parameters:
             !maximum particle threshold/4 * 1/number of surrounding elements * ratio
             !of particles in given CV to all surrounding CV's
             do l = 1,nint((max_thresh/4.0)*(1.0/size(ele_nums))*ratio)
                temp_part => null()
                call allocate(temp_part, size(particle%position), size(particle%local_coords), particle_lists(group_arrays(i))%total_attributes)
                temp_part%id_number = id+1
                temp_part%list_id = particle%list_id
                temp_part%proc_id = proc_num
                temp_part%element=ele_nums(j)
                !randomly select local coords within the element, ensuring coords are within cv
                call random_number(rand_val)
                max_lcoord=rand_val/(1/0.49)+0.51!lcoords for cv range from 0.51<x<1
                node_coord(:)=1
                node_coord(k)=0
                do while (minloc(node_coord,DIM=1)==k)
                   call set_spawned_lcoords(max_lcoord, node_coord, node_num, node_numbers(j,:))
                end do
                temp_part%local_coords = node_coord

                !Convert newly spawned particle's local coords to global coords
                temp_part%position = eval_field(temp_part%element, xfield, temp_part%local_coords)
                if (copy_parents) then
                   !Copy parent particle attributes
                   temp_part%attributes(:) = particle%attributes(:)
                   temp_part%old_attributes(:) = particle%old_attributes(:)
                   temp_part%old_fields(:) = particle%old_fields(:)
                else
                   temp_part%attributes(:) = average_attributes(i)%col(:)
                   temp_part%old_attributes(:) = average_old_attributes(i)%col(:)
                   temp_part%old_fields(:) = average_old_fields(i)%col(:)
                end if

                node_part_count(node_num) = node_part_count(node_num) + 1

                !add particle to relevant particle list
                call insert(temp_part, particle_lists(group_arrays(i)))

                !add particle to temp particle list
                call temp_list_insert(temp_part,node_particles(i,node_num))

                add_particles(i) = add_particles(i) + 1
                id = id + 1
                particle_lists(group_arrays(i))%proc_part_count = particle_lists(group_arrays(i))%proc_part_count + 1
             end do
          end do
       end do
    end do

    allocate(remove_particles(size(group_arrays)))
    remove_particles = 0

    !Call delete_particles if number of particles now exceeds the maximum threshold
    do while(node_part_count(node_num)>max_thresh)
       call delete_particles(node_particles(:,node_num), group_arrays, remove_particles, node_part_count=node_part_count(node_num))
       add_particles(:) = add_particles(:) - remove_particles(:)
    end do

    deallocate(remove_particles)
    deallocate(node_coord)
    deallocate(node_numbers)
    deallocate(node_loc)
    deallocate(ele_val)
    if (allocated(average_attributes)) then
       deallocate(average_attributes)
       deallocate(average_old_attributes)
       deallocate(average_old_fields)
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

    integer, dimension(4,4) :: permutation = reshape([1,2,3,4, 2,1,3,4, 2,3,1,4, 2,3,4,1], [4,4])
    real, dimension(4) :: work
    real :: tmp_res, rand_val
    integer :: i, j

    max_lcoord = max(0.51, min(max_lcoord, 0.999))
    ! set up the node coordinates to be permuted
    ! looks like: [x, (1 - x) * rand(), (1 - x - y) * rand(), 1 - x - y - z]
    ! depending on the number of coordinates
    work(1) = max_lcoord
    tmp_res = 1 - max_lcoord
    do j = 2, size(node_coord) - 1
       call random_number(rand_val)
       work(j) = tmp_res * rand_val
       tmp_res = tmp_res - work(j)
    end do
    work(size(node_coord)) = tmp_res

    do i = 1,size(node_coord)
       ! determine the node index corresponding to the
       ! target node number
       if (node_num == node_numbers(i)) exit
    end do
    assert(i<=size(node_coord))
    ! set the coordinates according to the permutation for this index
    do j = 1, size(node_coord)
       node_coord(j) = work(permutation(j,i))
    end do

  end subroutine set_spawned_lcoords

  subroutine multi_delete_particles(mult, node_particles, group_arrays, summed_particles, cap_percent)
    !Subroutine to call delete particles multiple times based on the mult factor
    !calculated from the current number of particles and the maximum threshold
    !> Factor to determine number of delete_particle calls
    integer, intent(in) :: mult
    !> Linked list of particles which exist on this node
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    !> Indicies in particle_lists present
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
    !Loop to find the largest power of 2 less than or equal to mult factor
    do while (.not. power_set)
       if (mult>=2**j) then
          j=j+1
       else
          power=j-1
          power_set=.true.
       end if
    end do

    do i = 1,power
       call delete_particles(node_particles, group_arrays, remove_particles, cap_percent=cap_percent)
       summed_particles=summed_particles-remove_particles
    end do

    deallocate(remove_particles)

  end subroutine multi_delete_particles

  subroutine delete_particles(node_particles, group_arrays, remove_particles, cap_percent, node_part_count)
    !Subroutine to delete particles in a control volume based off the number of particles present

    !> Linked list of particles which exist on this node
    type(detector_linked_list), intent(inout), dimension(:) :: node_particles
    !> Indicies in particle_lists present
    integer, intent(in), dimension(:) :: group_arrays
    !> Array to sum number of deleted particles per particle group
    integer, dimension(:), intent(inout) :: remove_particles
    !> Parameter to determine if spawning/deleting will be capped per group
    real, optional, intent(in) :: cap_percent
    !> Array counting number of deleted particles
    real, optional, intent(inout) :: node_part_count

    type(detector_type), pointer :: particle
    type(detector_type), pointer :: temp_part

    logical :: delete_group
    integer :: j
    real :: rand_val

    !Check ratio of each particle group present and which groups will be deleted
    delete_group = .true.
    if (present(cap_percent)) then
       delete_group = .false.
       do j=1,size(group_arrays)
          !Check if a particle group makes up >= cap_percent of the total particles in the current control volume
          assert(cap_percent>50)
          if ((((node_particles(j)%length*1.0)/sum(node_particles(:)%length)*1.0)*100)>=cap_percent) then
             delete_group = .true.
          end if
       end do
    end if

    remove_particles = 0

    !return if no group is being deleted
    if (.not. delete_group) return

    !loop over all particles in each particle group, flip a coin, if coin is heads (r>0.5) delete the particle
    do j = 1,size(group_arrays)
       particle =>node_particles(j)%first
       do while(associated(particle))
          call random_number(rand_val)
          if (rand_val>0.5) then !Delete the particle
             !Remove from particle list
             call remove(particle,particle_lists(group_arrays(j)))
             temp_part =>particle%temp_next
             !Remove from temp list
             call temp_list_remove(particle,node_particles(j))
             remove_particles(j)=remove_particles(j)+1
             if (present(node_part_count)) node_part_count = node_part_count - 1
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
