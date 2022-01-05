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
!    You should have received a copy of the GNU Lesser General PublicS
!    License along with this library; if not, write to the Free Software
!    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
!    USA
#include "fdebug.h"

module particles
  use fldebug
  use iso_c_binding, only: C_NULL_CHAR, c_ptr, c_f_pointer
  use global_parameters, only:FIELD_NAME_LEN,OPTION_PATH_LEN, &
& PYTHON_FUNC_LEN, integer_size, real_size, is_active_process
  use futils, only: int2str, free_unit
  use elements
  use mpi_interfaces
  use parallel_tools
  use spud
  use embed_python, only: set_detectors_from_python, deallocate_c_array
  use parallel_fields
  use fields
  use profiler
  use state_module
  use field_options
  use detector_data_types
  use pickers
  use detector_tools
  use detector_parallel
  use detector_move_lagrangian
  use time_period
  use h5hut

  implicit none

  private

  public :: initialise_particles, move_particles, write_particles_loop, destroy_particles, &
            update_particle_attributes_and_fields, checkpoint_particles_loop, &
            get_particle_arrays, particle_lists, initialise_constant_particle_attributes, &
            initialise_particles_during_simulation


  ! One particle list for each subgroup
  type(detector_linked_list), allocatable, dimension(:), target, save :: particle_lists
  ! Timing info for group output
  type(time_period_type), allocatable, dimension(:), save :: output_CS

  !> Derived type to hold the number of scalar, vector and tensor attributes,
  !! old attributes and old fields for a particle subgroup
  type attr_counts_type
    integer, dimension(3) :: attrs, old_attrs, old_fields
  end type attr_counts_type

  !> Derived type to hold scalar, vector and tensor attributes
  type attr_vals_type
    real, dimension(:), allocatable :: s
    real, dimension(:,:), allocatable :: v
    real, dimension(:,:,:), allocatable :: t
  end type attr_vals_type

  interface allocate
    module procedure allocate_attr_vals
  end interface allocate

  interface deallocate
    module procedure deallocate_attr_vals
  end interface deallocate

contains
  !> Allocate an attr_vals_type structure, with the given number
  !! of scalar, vector and tensor attributes, and the geometric dimension
  !! of the problem.
  subroutine allocate_attr_vals(vals, dim, counts)
    !> Structure to allocate
    type(attr_vals_type), pointer :: vals
    !> Geometric dimension
    integer, intent(in) :: dim
    !> Counts of each rank of attribute
    integer, dimension(3), intent(in) :: counts

    allocate(vals)
    allocate(vals%s(counts(1)))
    allocate(vals%v(dim, counts(2)))
    allocate(vals%t(dim, dim, counts(3)))
  end subroutine allocate_attr_vals

  !> Deallocate an attr_vals_type structure
  subroutine deallocate_attr_vals(vals)
    type(attr_vals_type), pointer :: vals

    deallocate(vals%s)
    deallocate(vals%v)
    deallocate(vals%t)
    deallocate(vals)
  end subroutine deallocate_attr_vals

  !> Initialise particles and set up particle file headers (per particle array)
  subroutine initialise_particles(filename, state, global, setup_output, ignore_analytical, number_of_partitions)
    !> Experiment filename to prefix particle output files
    character(len=*), intent(in) :: filename
    !> Model state structure
    type(state_type), dimension(:), intent(in) :: state
    !> Use global/parallel picker queries to determine particle elements?
    logical, intent(in), optional :: global
    !> Whether to set up output files for particle lists
    logical, intent(in), optional :: setup_output
    !> Whether to ignore analytical particles (i.e. not from file)
    logical, intent(in), optional :: ignore_analytical
    !> Number of processes to use for reading particle data
    integer, intent(in), optional :: number_of_partitions

    character(len=FIELD_NAME_LEN) :: subname
    character(len=OPTION_PATH_LEN) :: group_path, subgroup_path

    type(vector_field), pointer :: xfield
    real :: current_time

    integer :: sub_particles
    integer :: i, j, k
    integer :: dim, particle_groups, total_arrays, list_counter
    integer, dimension(:), allocatable :: particle_arrays
    integer :: totaldet_global
    logical :: from_file, do_output, do_analytical, store_old_fields
    integer :: n_fields, n_oldfields, phase, f
    integer :: s_field, v_field, t_field ! field index variables
    integer :: s_oldfield, v_oldfield, t_oldfield
    integer, dimension(3) :: field_counts, old_field_counts
    type(attr_names_type) :: attr_names, old_attr_names, field_names, old_field_names
    type(attr_write_type) :: attr_write
    type(attr_counts_type) :: attr_counts
    type(field_phase_type) :: field_phases, old_field_phases

    ! field pointers to get their names, for old_field_names
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield

    character(len=*), dimension(3), parameter :: orders = ["scalar", "vector", "tensor"]
    character(len=*), dimension(3), parameter :: types = ["prescribed", "diagnostic", "prognostic"]

    ewrite(2,*) "In initialise_particles"

    do_output = .true.
    if (present(setup_output)) do_output = setup_output
    do_analytical = .true.
    if (present(ignore_analytical)) do_analytical = .not. ignore_analytical

    ! Check whether there are any particle groups to initialise
    particle_groups = option_count("/particles/particle_group")
    if (particle_groups == 0) return

    ! Set up particle lists
    allocate(particle_arrays(particle_groups))
    allocate(output_CS(particle_groups))

    total_arrays = 0
    do i = 1, particle_groups
      group_path = "/particles/particle_group["//int2str(i-1)//"]"

      call init_output_CS(output_CS(i), trim(group_path) // "/particle_io")

      ! count subgroups for this group
      particle_arrays(i) = option_count(trim(group_path) // "/particle_subgroup")
      total_arrays = total_arrays + particle_arrays(i)
    end do
    allocate(particle_lists(total_arrays))

    ! Allocate parameters from the coordinate field
    xfield => extract_vector_field(state(1), "Coordinate")
    call get_option("/geometry/dimension", dim)
    call get_option("/timestepping/current_time", current_time)

    ! calculate the number of fields and old fields that would have to be stored
    ! (each combination of field order and field type)
    field_counts(:) = 0
    old_field_counts(:) = 0
    do i = 1, 3
      do j = 1, 3
        field_counts(i) = field_counts(i) + &
             option_count("/material_phase/"//orders(i)//"_field/"//types(j)//"/particles/include_in_particles")
        old_field_counts(i) = old_field_counts(i) + &
             option_count("/material_phase/"//orders(i)//"_field/"//types(j)//"/particles/include_in_particles/store_old_field")
      end do
    end do
    n_fields = sum(field_counts)
    n_oldfields = sum(old_field_counts)

    ! allocate arrays to hold field names
    call allocate(field_names, field_counts)
    call allocate(old_field_names, old_field_counts)
    ! allocate arrays to hold field phases
    call allocate(field_phases, field_counts)
    call allocate(old_field_phases, old_field_counts)

    ! read the names of the fields if there are any
    ! this is both for fields that should be included in particles, and that
    ! should have their old values available to particles too
    s_field = 0
    v_field = 0
    t_field = 0
    s_oldfield = 0
    v_oldfield = 0
    t_oldfield = 0
    do i = 1, size(state)
      if (field_counts(1) > 0) then
        do j = 1, size(state(i)%scalar_names)
          sfield => extract_scalar_field(state(i), state(i)%scalar_names(j))
          if (sfield%option_path == "" .or. aliased(sfield)) then
            cycle
          else if (have_option(trim(complete_field_path(sfield%option_path)) // "/particles/include_in_particles")) then
            s_field = s_field + 1
            field_names%s(s_field) = state(i)%scalar_names(j)
            field_phases%s(s_field) = i

            if (have_option(trim(complete_field_path(sfield%option_path)) // "/particles/include_in_particles/store_old_field")) then
              s_oldfield = s_oldfield + 1
              old_field_names%s(s_oldfield) = state(i)%scalar_names(j)
              old_field_phases%s(s_oldfield) = i
            end if
          end if
        end do
      end if

      if (field_counts(2) > 0) then
        do j = 1, size(state(i)%vector_names)
          vfield => extract_vector_field(state(i), state(i)%vector_names(j))
          if (vfield%option_path == "" .or. aliased(vfield)) then
            cycle
          else if (have_option(trim(complete_field_path(vfield%option_path)) // "/particles/include_in_particles")) then
            v_field = v_field + 1
            field_names%v(v_field) = state(i)%vector_names(j)
            field_phases%v(v_field) = i

            if (have_option(trim(complete_field_path(vfield%option_path)) // "/particles/include_in_particles/store_old_field")) then
              v_oldfield = v_oldfield + 1
              old_field_names%v(v_oldfield) = state(i)%vector_names(j)
              old_field_phases%v(v_oldfield) = i
            end if
          end if
        end do
      end if

      if (field_counts(3) > 0) then
        do j = 1, size(state(i)%tensor_names)
          tfield => extract_tensor_field(state(i), state(i)%tensor_names(j))
          if (tfield%option_path == "" .or. aliased(tfield)) then
            cycle
          else if (have_option(trim(complete_field_path(tfield%option_path)) // "/particles/include_in_particles")) then
            t_field = t_field + 1
            field_names%t(t_field) = state(i)%tensor_names(j)
            field_phases%t(t_field) = i

            if (have_option(trim(complete_field_path(tfield%option_path)) // "/particles/include_in_particles/store_old_field")) then
              t_oldfield = t_oldfield + 1
              old_field_names%t(t_oldfield) = state(i)%tensor_names(j)
              old_field_phases%t(t_oldfield) = i
            end if
          end if
        end do
      end if
    end do
    assert(s_field + v_field + t_field == n_fields)
    assert(s_oldfield + v_oldfield + t_oldfield == n_oldfields)

    list_counter = 1
    do i = 1,particle_groups
       group_path = "/particles/particle_group["//int2str(i-1)//"]"

       do k = 1, particle_arrays(i)
          subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"

          ! If the option "from_file" exists, it means we are
          ! continuing the simulation after checkpointing and the
          ! reading of the particle positions must be done from a file
          from_file = have_option(trim(subgroup_path) // "/initial_position/from_file")

          ! But if we're flredecomping, we don't want to handle
          ! particles with analytically-specified positions (i.e. not
          ! from a file)
          if (.not. do_analytical .and. .not. from_file) cycle

          ! Set up the particle list structure
          call get_option(trim(subgroup_path) // "/name", subname)

          ! Register this I/O list with a global list of detectors/particles
          call register_detector_list(particle_lists(list_counter))

          !Set list id
          particle_lists(list_counter)%id = list_counter

          ! Find number of attributes, old attributes, and names of each
          call attr_names_and_count(trim(subgroup_path) // "/attributes/scalar_attribute", &
               attr_names%s, old_attr_names%s, attr_names%sn, old_attr_names%sn, attr_write%s, &
               attr_counts%attrs(1), attr_counts%old_attrs(1))
          call attr_names_and_count(trim(subgroup_path) // "/attributes/vector_attribute", &
               attr_names%v, old_attr_names%v, attr_names%vn, old_attr_names%vn, attr_write%v, &
               attr_counts%attrs(2), attr_counts%old_attrs(2))
          call attr_names_and_count(trim(subgroup_path) // "/attributes/tensor_attribute", &
               attr_names%t, old_attr_names%t, attr_names%tn, old_attr_names%tn, attr_write%t, &
               attr_counts%attrs(3), attr_counts%old_attrs(3))

          ! save names in the detector list -- this will allocate and assign values
          ! as expected
          particle_lists(list_counter)%attr_names = attr_names
          particle_lists(list_counter)%old_attr_names = old_attr_names
          particle_lists(list_counter)%attr_write = attr_write

          ! If any attributes are from fields, we'll need to store old fields too
          store_old_fields = .false.
          if (option_count(trim(subgroup_path) // "/attributes/scalar_attribute/value_on_advection/python_fields") > 0 .or. &
              option_count(trim(subgroup_path) // "/attributes/scalar_attribute_array/value_on_advection/python_fields") > 0 .or. &
              option_count(trim(subgroup_path) // "/attributes/vector_attribute/value_on_advection/python_fields") > 0 .or. &
              option_count(trim(subgroup_path) // "/attributes/vector_attribute_array/value_on_advection/python_fields") > 0 .or. &
              option_count(trim(subgroup_path) // "/attributes/tensor_attribute/value_on_advection/python_fields") > 0 .or. &
              option_count(trim(subgroup_path) // "/attributes/tensor_attribute_array/value_on_advection/python_fields") > 0) then
             store_old_fields = .true.
          end if

          if (store_old_fields) then
            attr_counts%old_fields(:) = old_field_counts(:)
            ! only copy old field names if they're required
            particle_lists(list_counter)%field_names = field_names
            particle_lists(list_counter)%old_field_names = old_field_names

            ! and the field phases so we can look them up later
            particle_lists(list_counter)%field_phases = field_phases
            particle_lists(list_counter)%old_field_phases = old_field_phases
          else
            attr_counts%old_fields(:) = 0

            ! allocate empty arrays for names and phases
            call allocate(particle_lists(list_counter)%field_names, [0, 0, 0])
            call allocate(particle_lists(list_counter)%old_field_names, [0, 0, 0])
            call allocate(particle_lists(list_counter)%field_phases, [0, 0, 0])
            call allocate(particle_lists(list_counter)%old_field_phases, [0, 0, 0])
          end if

          ! assign the total number of list slices for each kind of attribute
          ! this is used mostly for transferring detectors across processes
          particle_lists(list_counter)%total_attributes(1) = &
               total_attributes(attr_counts%attrs, dim)
          particle_lists(list_counter)%total_attributes(2) = &
               total_attributes(attr_counts%old_attrs, dim)
          particle_lists(list_counter)%total_attributes(3) = &
               total_attributes(attr_counts%old_fields, dim)

          ! Enable particles to drift with the mesh
          if (have_option("/particles/move_with_mesh")) then
            particle_lists(list_counter)%move_with_mesh = .true.
          end if

          if (is_active_process) then
            ! Read particles from options -- only if this process is currently active (as defined in flredecomp)
            if (from_file) then
              call get_option(trim(subgroup_path) // "/initial_position/from_file/number_of_particles", sub_particles)
              call read_particles_from_file(sub_particles, subname, subgroup_path, &
                   particle_lists(list_counter), xfield, dim, &
                   attr_counts, attr_names, old_attr_names, old_field_names, &
                   number_of_partitions)
            else
              call read_particles_from_python(subname, subgroup_path, &
                   particle_lists(list_counter), xfield, dim, &
                   current_time, state, attr_counts, sub_particles, global=global)
            end if
          end if

          particle_lists(list_counter)%total_num_det = sub_particles

          if (do_output) then
            ! Only set up output if we need to (i.e. actually running,
            ! not flredecomping)
            call set_particle_output_file(subname, filename, particle_lists(list_counter))
          end if

          ! Get options for lagrangian particle movement
          call read_detector_move_options(particle_lists(list_counter), "/particles")

          ! Make sure to deallocate attribute names before moving on
          call deallocate(attr_names)
          call deallocate(old_attr_names)

          list_counter = list_counter + 1
       end do
    end do

    ! And finally some sanity checks
    list_counter=1
    do i = 1,particle_groups
      group_path = "/particles/particle_group["//int2str(i-1)//"]"
      do k = 1,particle_arrays(i)
        subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
        call get_option(trim(subgroup_path)//"/name",subname)
        totaldet_global=particle_lists(list_counter)%length
        call allsum(totaldet_global)
        ewrite(2,*) "Found", particle_lists(list_counter)%length, "local and ", totaldet_global, "global particles for particle array ", trim(subname)

        assert(totaldet_global==particle_lists(list_counter)%total_num_det)
        list_counter = list_counter + 1
      end do
    end do

    deallocate(particle_arrays)
    call deallocate(field_names)
    call deallocate(old_field_names)
  end subroutine initialise_particles

  !> Initialise particles for times greater than 0
  subroutine initialise_particles_during_simulation(state, current_time, dt)
    !> Model state structure
    type(state_type), dimension(:), intent(in) :: state
    !> Current simulation time
    real, intent(in) :: current_time
    !> Current model timestep
    real, intent(in) :: dt

    integer :: i, k, j, dim, id_number
    integer :: particle_groups, particle_subgroups, list_counter, sub_particles
    integer :: nb_part_created
    integer, dimension(:), allocatable :: init_check

    type(vector_field), pointer :: xfield
    type(attr_counts_type) :: attr_counts
    type(attr_names_type) :: attr_names, old_attr_names, field_names, old_field_names
    type(attr_write_type) :: attr_write
    type(detector_type), pointer :: first_newly_init_part
    integer, dimension(3) :: field_counts, old_field_counts

    character(len=OPTION_PATH_LEN) :: group_path, subgroup_path
    character(len=FIELD_NAME_LEN) :: subname
    character(len=PYTHON_FUNC_LEN) :: script

    character(len=*), dimension(3), parameter :: orders = ["scalar", "vector", "tensor"]
    character(len=*), dimension(3), parameter :: types = ["prescribed", "diagnostic", "prognostic"]

    logical :: store_old_fields

    ! Check whether there are any particles.
    particle_groups = option_count("/particles/particle_group")
    if (particle_groups == 0) return

    ! calculate the number of fields and old fields that would have to be stored
    ! (each combination of field order and field type)
    old_field_counts(:) = 0
    do i = 1, 3
      do j = 1, 3
        old_field_counts(i) = old_field_counts(i) + &
             option_count("/material_phase/"//orders(i)//"_field/"//types(j)//"/particles/include_in_particles/store_old_field")
      end do
    end do

    ! Allocate parameters from the coordinate field
    xfield => extract_vector_field(state(1), "Coordinate")
    call get_option("/geometry/dimension", dim)

    list_counter = 1

    ewrite(2,*) "In initialise_particles_during_simulation"

    !Check if initialise_during_simulation is enabled
    do i = 1, particle_groups
      group_path = "/particles/particle_group["//int2str(i-1)//"]"
      particle_subgroups = option_count(trim(group_path) // "/particle_subgroup")
      do k = 1, particle_subgroups
        subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
        if (have_option(trim(subgroup_path)//"/initialise_during_simulation")) then
          ! Find number of attributes, old attributes, and names of each
          call attr_names_and_count(trim(subgroup_path) // "/attributes/scalar_attribute", &
            attr_names%s, old_attr_names%s, attr_names%sn, old_attr_names%sn, attr_write%s, &
            attr_counts%attrs(1), attr_counts%old_attrs(1))
          call attr_names_and_count(trim(subgroup_path) // "/attributes/vector_attribute", &
            attr_names%v, old_attr_names%v, attr_names%vn, old_attr_names%vn, attr_write%v, &
            attr_counts%attrs(2), attr_counts%old_attrs(2))
          call attr_names_and_count(trim(subgroup_path) // "/attributes/tensor_attribute", &
            attr_names%t, old_attr_names%t, attr_names%tn, old_attr_names%tn, attr_write%t, &
            attr_counts%attrs(3), attr_counts%old_attrs(3))

          ! If any attributes are from fields, we'll need to store old fields too
          store_old_fields = .false.
          if (option_count(trim(subgroup_path) // "/attributes/scalar_attribute/value_on_advection/python_fields") > 0 .or. &
              option_count(trim(subgroup_path) // "/attributes/scalar_attribute_array/value_on_advection/python_fields") > 0 .or. &
              option_count(trim(subgroup_path) // "/attributes/vector_attribute/value_on_advection/python_fields") > 0 .or. &
              option_count(trim(subgroup_path) // "/attributes/vector_attribute_array/value_on_advection/python_fields") > 0 .or. &
              option_count(trim(subgroup_path) // "/attributes/tensor_attribute/value_on_advection/python_fields") > 0 .or. &
              option_count(trim(subgroup_path) // "/attributes/tensor_attribute_array/value_on_advection/python_fields") > 0) then
            store_old_fields = .true.
          end if

          if (store_old_fields) then
            attr_counts%old_fields(:) = old_field_counts(:)
          else
            attr_counts%old_fields(:) = 0
          end if

          call get_option(trim(subgroup_path)//"/initialise_during_simulation/python", script)

          id_number = particle_lists(list_counter)%proc_part_count
          call get_option(trim(subgroup_path) // "/name", subname)
          first_newly_init_part => particle_lists(list_counter)%last
          call read_particles_from_python(subname, subgroup_path, particle_lists(list_counter), xfield, dim, &
            current_time, state, attr_counts, sub_particles, global=.true., &
            id_number=id_number, script=script, nb_part_created=nb_part_created)

          particle_lists(list_counter)%total_num_det = particle_lists(list_counter)%total_num_det + sub_particles

          if (nb_part_created > 0) then
            if (.not. associated(first_newly_init_part)) then
              first_newly_init_part => particle_lists(list_counter)%first
            else
              first_newly_init_part => first_newly_init_part%next
            end if
            call update_particle_subgroup_attributes_and_fields( &
              state, current_time, dt, subgroup_path, particle_lists(list_counter), &
              initial=.true., first_newly_init_part=first_newly_init_part, nb_part_created=nb_part_created)
          end if
        end if
        list_counter = list_counter + 1
      end do
    end do


  end subroutine initialise_particles_during_simulation

  !> Get the names and count of all attributes and old attributes for
  !! a given attribute rank for a particle subgroup
  subroutine attr_names_and_count(key, names, old_names, dims, old_dims, to_write, count, old_count)
    !> Prefix key to an attribute rank within a subgroup
    character(len=*), intent(in) :: key
    !> Output arrays for attribute names
    character(len=*), dimension(:), allocatable, intent(out) :: names, old_names
    !> Output arrays for attribute dimensions
    integer, dimension(:), allocatable, intent(out) :: dims, old_dims
    !> Output arrays for whether to write attributes
    logical, dimension(:), allocatable, intent(out) :: to_write
    !> Output attribute counts
    integer, intent(out) :: count, old_count

    integer :: i, old_i, single_count, array_count, single_old_count, array_old_count
    character(len=FIELD_NAME_LEN) :: array_key, subkey

    ! array-valued attribute name
    array_key = trim(key) // "_array"

    ! get option count so we can allocate the names array
    single_count = option_count(key)
    array_count = option_count(array_key)

    single_old_count = option_count(key//"/value_on_advection/python_fields/store_old_attribute")
    array_old_count = option_count(trim(array_key)//"/value_on_advection/python_fields/store_old_attribute")

    allocate(names(single_count + array_count))
    allocate(old_names(single_old_count + array_old_count))

    allocate(to_write(single_count + array_count))

    allocate(dims(single_count + array_count))
    allocate(old_dims(single_old_count + array_old_count))

    ! names for single-valued attributes
    old_i = 1
    do i = 1, single_count
      ! get the attribute's name
      write(subkey, "(a,'[',i0,']')") key, i-1
      call get_option(trim(subkey)//"/name", names(i))
      ! we set single-valued attributes to have a dimension of 0 to distinguish from
      ! a length 1 array attribute
      dims(i) = 0

      to_write(i) = .not. have_option(trim(subkey)//"/exclude_from_output")

      if (have_option(trim(subkey)//"/value_on_advection/python_fields/store_old_attribute")) then
        ! prefix with "old%" to distinguish from current attribute
        old_names(old_i) = "old%" // trim(names(i))
        old_dims(old_i) = 0
        old_i = old_i + 1
      end if
    end do

    ! names for array-valued attributes
    do i = 1, array_count
      write(subkey, "(a,'[',i0,']')") trim(array_key), i-1
      call get_option(trim(subkey)//"/name", names(i+single_count))
      call get_option(trim(subkey)//"/dimension", dims(i+single_count))

      to_write(i+single_count) = .not. have_option(trim(subkey)//"/exclude_from_output")

      if (have_option(trim(subkey)//"/value_on_advection/python_fields/store_old_attribute")) then
        old_names(old_i) = "old%" // trim(names(i+single_count))
        old_dims(old_i) = dims(i+single_count)
        old_i = old_i + 1
      end if
    end do

    ! compute the number of attribute arrays
    count = single_count + sum(dims)
    old_count = single_old_count + sum(old_dims)
  end subroutine attr_names_and_count

  !> Initialise particles which are defined by a Python function
  subroutine read_particles_from_python(subgroup_name, subgroup_path, &
      p_list, xfield, dim, current_time, state, attr_counts, n_particles, &
      global, id_number, script, nb_part_created)

    !> Name of the particles' subgroup
    character(len=FIELD_NAME_LEN), intent(in) :: subgroup_name
    !> Path prefix for the subgroup in options
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path
    !> Detector list to hold the particles
    type(detector_linked_list), intent(inout) :: p_list
    !> Coordinate vector field
    type(vector_field), pointer, intent(in) :: xfield
    !> Geometry dimension
    integer, intent(in) :: dim
    !> Current model time, for passing through to Python functions
    real, intent(in) :: current_time
    !> Model state structure
    type(state_type), dimension(:), intent(in) :: state
    !> Counts of attributes, old attributes, and old fields
    !! for each scalar, vector and tensor
    type(attr_counts_type), intent(in) :: attr_counts
    !> Number of particles being initialized
    integer, intent(out) :: n_particles
    !> Whether to consider this particle in a global element search
    logical, intent(in), optional :: global
    !> ID number of last particle currently in list
    integer, optional, intent(in) :: id_number
    !> Python script used by initialise_during_simulation
    character(len=PYTHON_FUNC_LEN), optional, intent(in) :: script
    !> Number of particles created as part of initialise_during_simulation
    integer, optional, intent(out) :: nb_part_created

    integer :: i, proc_num, stat, offset
    character(len=PYTHON_FUNC_LEN) :: func
    real, allocatable, dimension(:,:) :: coords ! all particle coordinates, from python
    ! if we don't know how many particles we're getting, we need a C pointer
    type(c_ptr) :: coord_ptr
    ! and a fortran pointer
    real, pointer, dimension(:,:) :: coord_array_ptr
    real :: dt

    proc_num = getprocno()

    ewrite(2,*) "Reading particles from options"

    if (present(script)) then
      func=script
      nb_part_created = 0
    else
      call get_option(trim(subgroup_path)//"/initial_position/python", func)
    end if
    call get_option("/timestepping/timestep", dt)

    call set_detectors_from_python(func, len(func), dim, current_time, coord_ptr, n_particles, stat)
    call c_f_pointer(coord_ptr, coord_array_ptr, [dim, n_particles])
    allocate(coords(dim, n_particles))
    if (n_particles==0) return
    coords = coord_array_ptr

    call deallocate_c_array(coord_ptr)
    offset = 0
    if (present(id_number)) offset = id_number
    do i = 1, n_particles
      call create_single_particle(p_list, xfield, coords(:,i), i+offset, proc_num, dim, attr_counts, &
          global=global, nb_part_created=nb_part_created)
    end do

    deallocate(coords)
  end subroutine read_particles_from_python

  !> Read attributes for all ranks from an H5Part file
  subroutine read_attrs(h5_id, dim, counts, names, vals, prefix)
    !> h5 file to read from
    !! it's assumed this has been set up to read from the right place!
    integer(kind=8), intent(in) :: h5_id
    !> spatial dimension
    integer, intent(in) :: dim
    !> counts of scalar/vector/tensor attributes
    integer, dimension(3), intent(in) :: counts
    !> attribute names to read from the file
    type(attr_names_type), intent(in) :: names
    !> SVT values to hold the output
    type(attr_vals_type), intent(inout) :: vals
    !> Optional prefix to attribute names
    character(len=*), intent(in), optional :: prefix

    integer :: i, j, k, ii, val_i
    integer(kind=8) :: h5_ierror
    character(len=FIELD_NAME_LEN) :: p

    p = ""
    if (present(prefix)) p = prefix

    val_i = 1
    scalar_attr_loop: do i = 1, size(names%s)
      if (names%sn(i) == 0) then
        ! single-valued attribute
        h5_ierror = h5pt_readdata_r8(h5_id, &
             trim(p)//trim(names%s(i)), vals%s(val_i))
        val_i = val_i + 1
      else
        do ii = 1, names%sn(i)
          ! inner loop for array-valued attribute
          h5_ierror = h5pt_readdata_r8(h5_id, &
               trim(p)//trim(names%s(i))//int2str(ii), vals%s(val_i))
          val_i = val_i + 1
        end do
      end if
    end do scalar_attr_loop

    val_i = 1
    vector_attr_loop: do i = 1, size(names%v)
      if (names%vn(i) == 0) then
        ! single-valued attribute
        do j = 1, dim
          h5_ierror = h5pt_readdata_r8(h5_id, &
               trim(p)//trim(names%v(i))//"_"//int2str(j-1), vals%v(j,val_i))
        end do
        val_i = val_i + 1
      else
        do ii = 1, names%vn(i)
          ! inner loop for array-valued attribute
          do j = 1, dim
            h5_ierror = h5pt_readdata_r8(h5_id, &
                 trim(p)//trim(names%v(i))//int2str(ii)//"_"//int2str(j-1), vals%v(j,val_i))
          end do
          val_i = val_i + 1
        end do
      end if
    end do vector_attr_loop

    val_i = 1
    tensor_attr_loop: do i = 1, size(names%t)
      if (names%tn(i) == 0) then
        ! single-valued attribute
        do j = 1, dim
          do k = 1, dim
            h5_ierror = h5pt_readdata_r8(h5_id, &
                 trim(p)//trim(names%t(i))//"_"//int2str((k-1)*dim + (j-1)), &
                 vals%t(j,k,val_i))
          end do
        end do
        val_i = val_i + 1
      else
        do ii = 1, names%tn(i)
          ! inner loop for array-valued attribute
          do j = 1, dim
            do k = 1, dim
              h5_ierror = h5pt_readdata_r8(h5_id, &
                   trim(p)//trim(names%t(i))//int2str(ii)//"_"//int2str((k-1)*dim + (j-1)), &
                   vals%t(j,k,val_i))
            end do
          end do
          val_i = val_i + 1
        end do
      end if
    end do tensor_attr_loop
  end subroutine read_attrs

  !> Read particles in the given subgroup from a checkpoint file
  subroutine read_particles_from_file(n_particles, subgroup_name, subgroup_path, &
       p_list, xfield, dim, &
       attr_counts, attr_names, old_attr_names, old_field_names, &
       n_partitions)
    !> Number of particles in this subgroup
    integer, intent(in) :: n_particles
    !> Name of the particles' subgroup
    character(len=FIELD_NAME_LEN), intent(in) :: subgroup_name
    !> Path prefix for the subgroup in options
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path
    !> Detector list to hold the particles
    type(detector_linked_list), intent(inout) :: p_list
    !> Coordinate vector field
    type(vector_field), pointer, intent(in) :: xfield
    !> Geometry dimension
    integer, intent(in) :: dim
    !> Counts of attributes, old attributes, and old fields
    !! for each scalar, vector and tensor
    type(attr_counts_type), intent(in) :: attr_counts
    !> Names of attributes to store on this subgroup
    type(attr_names_type), intent(in) :: attr_names
    !> The attributes for which an old value should be checkpointed
    type(attr_names_type), intent(in) :: old_attr_names
    !> Names of fields for which old values should be checkpointed
    type(attr_names_type), intent(in) :: old_field_names

    !> Optional parameter during recomposition to control the
    !! processes which are involved in reading from file
    integer, intent(in), optional :: n_partitions

    integer :: i
    integer :: ierr, commsize, rank, id(1), proc_id(1)
    integer :: input_comm, world_group, input_group ! opaque MPI pointers
    real, allocatable, dimension(:) :: positions ! particle coordinates
    character(len=OPTION_PATH_LEN) :: particles_cp_filename
    integer(kind=8) :: h5_ierror, h5_id, h5_prop, view_start, view_end ! h5hut state
    integer(kind=8), dimension(:), allocatable :: npoints, part_counter ! number of points for each rank to read
    type(attr_vals_type), pointer :: attr_vals, old_attr_vals, old_field_vals ! scalar/vector/tensor arrays

    ewrite(2,*) "Reading particles from file"

    ! create new mpi group for active particles only
    ! non-active processes are already not in this routine,
    ! so we don't have to worry about them
    if (present(n_partitions)) then
      call mpi_comm_group(MPI_COMM_FEMTOOLS, world_group, ierr)
      call mpi_group_incl(world_group, n_partitions, [(i, i=0, n_partitions-1)], input_group, ierr)
      call mpi_comm_create_group(MPI_COMM_FEMTOOLS, input_group, 0, input_comm, ierr)
    else
      input_comm = MPI_COMM_FEMTOOLS
    end if

    ! allocate arrays to hold positions and attributes for a single particle
    allocate(positions(dim))
    call allocate(attr_vals, dim, attr_counts%attrs)
    call allocate(old_attr_vals, dim, attr_counts%old_attrs)
    call allocate(old_field_vals, dim, attr_counts%old_fields)

    ! set up the checkpoint file for reading
    call get_option(trim(subgroup_path) // "/initial_position/from_file/file_name", particles_cp_filename)

    h5_prop = h5_createprop_file()
    ! because we're reading separate particle counts per core
    ! we can't use collective IO
    h5_ierror = h5_setprop_file_mpio_independent(h5_prop, input_comm)
    assert(h5_ierror == H5_SUCCESS)

    h5_id = h5_openfile(trim(particles_cp_filename), H5_O_RDONLY, h5_prop)
    h5_ierror = h5_closeprop(h5_prop)
    h5_ierror = h5_setstep(h5_id, int(1, 8))

    ! determine the number of particles we are to initiliase
    call mpi_comm_size(MPI_COMM_FEMTOOLS, commsize, ierr)
    call mpi_comm_rank(MPI_COMM_FEMTOOLS, rank, ierr)
    allocate(npoints(commsize))
    allocate(part_counter(commsize))
    h5_ierror = h5_readfileattrib_i8(h5_id, "npoints", npoints)
    h5_ierror = h5_readfileattrib_i8(h5_id, "part_counter", part_counter)
    h5_ierror = h5pt_setnpoints(h5_id, npoints(rank+1))

    ! figure out our local offset into the file
    h5_ierror = h5pt_getview(h5_id, view_start, view_end)

    do i = 1, npoints(rank+1)

      ! set view to read this particle
      h5_ierror = h5pt_setview(h5_id, int(view_start + i - 1, 8), int(view_start + i - 1, 8))

      if (dim >= 1) &
           h5_ierror = h5pt_readdata_r8(h5_id, "x", positions(1))
      if (dim >= 2) &
           h5_ierror = h5pt_readdata_r8(h5_id, "y", positions(2))
      if (dim >= 3) &
           h5_ierror = h5pt_readdata_r8(h5_id, "z", positions(3))

      h5_ierror = h5pt_readdata_i4(h5_id, "id", id(1))
      h5_ierror = h5pt_readdata_i4(h5_id, "proc_id", proc_id(1))

      ! batched reads of scalar, vector, tensor values of each kind of attribute
      call read_attrs(h5_id, dim, attr_counts%attrs, attr_names, attr_vals)
      call read_attrs(h5_id, dim, attr_counts%old_attrs, old_attr_names, old_attr_vals)
      call read_attrs(h5_id, dim, attr_counts%old_fields, old_field_names, old_field_vals, prefix="old%")

      ! don't use a global check for this particle
      call create_single_particle(p_list, xfield, &
           positions, id(1), proc_id(1), dim, &
           attr_counts, attr_vals, old_attr_vals, old_field_vals, global=.false.)
    end do

    ! reset proc_particle_count
    p_list%proc_part_count = part_counter(rank+1)

    h5_ierror = h5_closefile(h5_id)

    deallocate(positions)
    deallocate(npoints)
    deallocate(part_counter)

    deallocate(attr_vals)
    deallocate(old_attr_vals)
    deallocate(old_field_vals)

    if (present(n_partitions)) then
      call mpi_comm_free(input_comm, ierr)
      call mpi_group_free(input_group, ierr)
    end if
  end subroutine read_particles_from_file

  subroutine set_particle_output_file(subname, filename, p_list)
    !! Set up the particle output file for a single subgroup

    type(detector_linked_list), intent(inout) :: p_list
    character(len=*), intent(in) :: filename
    character(len=FIELD_NAME_LEN), intent(in) :: subname

    p_list%h5_id = h5_openfile(trim(filename) // '.particles.' // trim(subname) // '.h5part', H5_O_WRONLY, H5_PROP_DEFAULT)

    ! optionally set any file attributes here?
  end subroutine set_particle_output_file

  !> Allocate a single particle, populate and insert it into the given list
  !! In parallel, first check if the particle would be local and only allocate if it is
  subroutine create_single_particle(detector_list, xfield, position, id, proc_id, dim, &
      attr_counts, attr_vals, old_attr_vals, old_field_vals, global, nb_part_created)
    !> The detector list to hold the particle
    type(detector_linked_list), intent(inout) :: detector_list
    !> Coordinate vector field
    type(vector_field), pointer, intent(in) :: xfield
    !> Spatial position of the particle
    real, dimension(xfield%dim), intent(in) :: position
    !> Unique ID number for this particle
    integer, intent(in) :: id
    !> Procces ID on which this particle was created
    integer, intent(in) :: proc_id
    !> Geometry dimension
    integer, intent(in) :: dim
    !> Counts of scalar, vector and tensor attributes, old attributes
    !! and old fields to store on the particle
    type(attr_counts_type), intent(in) :: attr_counts
    !> If provided, initialise the particle's attributes directly
    type(attr_vals_type), intent(in), optional :: attr_vals, old_attr_vals, old_field_vals
    !> Whether to create this particle in a collective operation (true)
    !! or for the local processor only (false).
    !! This affects the inquiry of the element owning the particle
    logical, intent(in), optional :: global
    !> Number of particles created as part of initialise_during_simulation
    integer, intent(inout), optional :: nb_part_created

    type(detector_type), pointer :: detector
    type(element_type), pointer :: shape
    real, dimension(xfield%dim+1) :: lcoords
    integer :: element

    real :: dt
    logical :: picker_global = .true.

    if (present(global)) picker_global = global

    shape => ele_shape(xfield,1)
    assert(xfield%dim+1==local_coord_count(shape))

    ! Determine element and local_coords from position
    call picker_inquire(xfield, position, element, local_coord=lcoords, global=picker_global)
    call get_option("/timestepping/timestep", dt)
    ! If we're in parallel and don't own the element, skip this particle
    if (isparallel()) then
      if (element<0) return
      if (.not.element_owned(xfield,element)) return
    else
      ! In serial make sure the particle is in the domain
      ! unless we have the write_nan_outside override
      if (element<0 .and. .not.detector_list%write_nan_outside) then
        ewrite(-1,*) "Dealing with particle ", id, " proc_id:", proc_id
        FLExit("Trying to initialise particle outside of computational domain")
      end if
    end if

    ! Otherwise, allocate and insert particle
    allocate(detector)
    allocate(detector%position(xfield%dim))
    allocate(detector%local_coords(local_coord_count(shape)))
    call insert(detector, detector_list)

    if (present(nb_part_created)) then
      nb_part_created = nb_part_created + 1
    end if

    ! Populate particle
    detector%position = position
    detector%element = element
    detector%local_coords = lcoords
    detector%id_number = id
    detector%proc_id = proc_id
    detector%list_id = detector_list%id
    detector_list%proc_part_count = max(detector_list%proc_part_count,id)

    ! allocate space to store all attributes on the particle
    allocate(detector%attributes(total_attributes(attr_counts%attrs, dim)))
    allocate(detector%old_attributes(total_attributes(attr_counts%old_attrs, dim)))
    allocate(detector%old_fields(total_attributes(attr_counts%old_fields, dim)))

    ! copy attributes if they're present, otherwise initialise to zero
    call copy_attrs(detector%attributes, dim, attr_counts%attrs, attr_vals)
    call copy_attrs(detector%old_attributes, dim, attr_counts%old_attrs, old_attr_vals)
    call copy_attrs(detector%old_fields, dim, attr_counts%old_fields, old_field_vals)
  end subroutine create_single_particle

  !> Convert an array of scalar, vector and tensor attribute counts
  !! to the total number of attribute slices (i.e. 1 per scalar attribute,
  !! 'dim' per vector, and 'dim*dim' per tensor)
  function total_attributes(counts, dim)
    !> Counts of scalar, vector and tensor attributes
    integer, dimension(3), intent(in) :: counts
    !> Geometry dimension
    integer, intent(in) :: dim
    integer :: total_attributes

    total_attributes = counts(1) + dim*counts(2) + dim*dim*counts(3)
  end function total_attributes

  !> Copy from an attr_vals_type to attribute arrays
  subroutine copy_attrs(dest, dim, counts, vals)
    !! Destination attribute array
    real, dimension(:), intent(out) :: dest
    !! Geometric dimension
    integer, intent(in) :: dim
    !! Attribute counts for each rank
    integer, dimension(3), intent(in) :: counts
    !! The attr_vals to copy from, if present
    type(attr_vals_type), intent(in), optional :: vals

    integer :: cur, i, j, k

    if (present(vals)) then
      cur = 1
      scalar_copy_loop: do i = 1, counts(1)
        dest(cur) = vals%s(i)
        cur = cur + 1
      end do scalar_copy_loop

      vector_copy_loop: do i = 1, counts(2)
        do j = 1, dim
          dest(cur) = vals%v(j,i)
          cur = cur + 1
        end do
      end do vector_copy_loop

      tensor_copy_loop: do i = 1, counts(3)
        do j = 1, dim
          do k = 1, dim
            dest(cur) = vals%t(j,k,i)
            cur = cur + 1
          end do
        end do
      end do tensor_copy_loop
    else
      dest(:) = 0.
    end if
  end subroutine copy_attrs

  !> Call move_lagrangian_detectors on all tracked particle groups
  subroutine move_particles(state, dt)
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: dt

    integer :: i, particle_groups

    ewrite(2,*) "In move_particles"
    call profiler_tic("particle_advection")

    particle_groups = option_count("/particles/particle_group")
    if (particle_groups == 0) return

    do i = 1, size(particle_lists)
      call move_lagrangian_detectors(state, particle_lists(i), dt)
    end do

    call profiler_toc("particle_advection")
  end subroutine move_particles

  !> Initialise constant attribute values before diagnostic fields are set
  subroutine initialise_constant_particle_attributes(state, subgroup_path, p_list)
    !!Routine to initialise constant attributes for MVF field
    type(state_type), dimension(:), intent(in) :: state
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path
    type(detector_linked_list), intent(in) :: p_list

    type(detector_type), pointer :: particle
    character(len=OPTION_PATH_LEN) :: attr_key

    real, allocatable, dimension(:,:) :: attribute_array
    real :: constant
    real, allocatable, dimension(:) :: vconstant
    real, allocatable, dimension(:,:) :: tconstant
    integer :: i, j, nparticles, n, i_single, attr_idx, dim
    integer :: nscalar, nvector, ntensor

    !Check if this processor contains particles
    nparticles = p_list%length

    if (nparticles.eq.0) then
       return
    end if

    ! get attribute sizes from the detector list
    nscalar = size(p_list%attr_names%s)
    nvector = size(p_list%attr_names%v)
    ntensor = size(p_list%attr_names%t)

    call get_option("/geometry/dimension", dim)
    allocate(vconstant(dim))
    allocate(tconstant(dim, dim))

    particle => p_list%first
    allocate(attribute_array(size(particle%attributes),nparticles))
    attribute_array(:,:) = 0

    !Scalar constants
    i_single = 1
    attr_idx = 1
    do i = 1, nscalar
       n = p_list%attr_names%sn(i)
       if (n == 0) then
          ! single-valued attribute
          attr_key = trim(subgroup_path) // '/attributes/scalar_attribute['//int2str(i_single-1)//']'
          i_single = i_single + 1
          if (have_option(trim(attr_key)//'/constant')) then
             call get_option(trim(attr_key)//'/constant', constant)
             attribute_array(attr_idx:attr_idx,:) = constant
          end if
          attr_idx = attr_idx + 1
       end if
    end do

    !Vector constants
    i_single = 1
    do i=1, nvector
       n = p_list%attr_names%vn(i)
       if (n == 0) then
          ! single-valued attribute
          attr_key = trim(subgroup_path) // '/attributes/vector_attribute['//int2str(i_single-1)//']'
          i_single = i_single + 1
          if (have_option(trim(attr_key)//'/constant')) then
             call get_option(trim(attr_key)//'/constant', vconstant)
             ! broadcast vector constant out to all particles
             attribute_array(attr_idx:attr_idx+dim-1,:) = spread(vconstant, 2, nparticles)
          end if
          attr_idx = attr_idx + dim
       end if
    end do

    !Tensor constants
    i_single = 1
    do i=1, ntensor
       n = p_list%attr_names%tn(i)
       if (n == 0) then
          ! single-valued attribute
          attr_key = trim(subgroup_path) // '/attributes/tensor_attribute['//int2str(i_single-1)//']'
          i_single = i_single + 1
          if (have_option(trim(attr_key)//'/constant')) then
             call get_option(trim(attr_key)//'/constant', tconstant)
             ! flatten tensor, then broadcast out to all particles
             attribute_array(attr_idx:attr_idx+dim**2-1,:) = spread(reshape(tconstant, [dim**2]), 2, nparticles)
          end if
          attr_idx = attr_idx + dim**2
       end if
    end do

    !Set constant attribute values
    particle => p_list%first
    do j = 1,nparticles
       particle%attributes = attribute_array(:,j)
       particle => particle%next
    end do
    deallocate(vconstant)
    deallocate(tconstant)
    deallocate(attribute_array)

  end subroutine initialise_constant_particle_attributes

  !> Update attributes and fields for every subgroup of every particle group
  subroutine update_particle_attributes_and_fields(state, time, dt, initial)
    type(state_type), dimension(:), intent(in) :: state
    real, intent(in) :: time
    real, intent(in) :: dt
    logical, intent(in), optional :: initial
    character(len = OPTION_PATH_LEN) :: group_path, subgroup_path

    integer :: i, k
    integer :: particle_groups, particle_subgroups, list_counter

    ! Check whether there are any particles.
    particle_groups = option_count("/particles/particle_group")
    if (particle_groups == 0) return

    ewrite(2,*) "In update_particle_attributes_and_fields"

    list_counter = 1

    do i = 1, particle_groups
      group_path = "/particles/particle_group["//int2str(i-1)//"]"
      particle_subgroups = option_count(trim(group_path) // "/particle_subgroup")
      do k = 1, particle_subgroups
        subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"

        if (particle_lists(list_counter)%length==0) then
          list_counter = list_counter + 1
          cycle
        end if

        call update_particle_subgroup_attributes_and_fields( &
          state, time, dt, subgroup_path, particle_lists(list_counter), initial=initial)
        list_counter = list_counter + 1
      end do
    end do
  end subroutine update_particle_attributes_and_fields

  !> Copy a structure of attribute names by rank
  !! to a 2D character array for passing to C
  subroutine copy_names_to_array(attr_names, name_array, attr_counts, attr_dims, prefix)
    type(attr_names_type), intent(in) :: attr_names
    character, allocatable, dimension(:,:), intent(out) :: name_array
    integer, dimension(3), intent(out) :: attr_counts
    integer, allocatable, dimension(:), intent(out), optional :: attr_dims
    character(len=*), intent(in), optional :: prefix

    integer :: i, j, k, n, np
    character(len=FIELD_NAME_LEN) :: p

    p = ""
    if (present(prefix)) p = prefix

    np = len_trim(p)

    ! determine number of names for each rank
    ! so that we can allocate the names array
    attr_counts(1) = size(attr_names%s)
    attr_counts(2) = size(attr_names%v)
    attr_counts(3) = size(attr_names%t)
    allocate(name_array(FIELD_NAME_LEN, sum(attr_counts)))

    if (present(attr_dims)) then
      allocate(attr_dims(sum(attr_counts)))
    end if

    ! unfortunately, we have to use a character array for C
    ! interoperability, and we can't assign an array from
    ! a character scalar, so we have to do explicit lops
    j = 1
    do i = 1, attr_counts(1)
      ! copy prefix
      do k = 1, np
        name_array(k,j) = p(k:k)
      end do

      ! copy attribute name
      n = min(FIELD_NAME_LEN - np - 1, len_trim(attr_names%s(i)))
      do k = 1, n
        name_array(k+np,j) = attr_names%s(i)(k:k)
      end do

      ! null terminate
      name_array(np+n+1,j) = C_NULL_CHAR

      ! possibly copy in the dimension
      if (present(attr_dims)) then
        attr_dims(j) = attr_names%sn(i)
      end if
      j = j + 1
    end do
    do i = 1, attr_counts(2)
      ! copy prefix
      do k = 1, np
        name_array(k,j) = p(k:k)
      end do

      n = min(FIELD_NAME_LEN - np - 1, len_trim(attr_names%v(i)))
      do k = 1, n
        name_array(k+np,j) = attr_names%v(i)(k:k)
      end do
      name_array(np+n+1,j) = C_NULL_CHAR
      if (present(attr_dims)) then
        attr_dims(j) = attr_names%vn(i)
      end if
      j = j + 1
    end do
    do i = 1, attr_counts(3)
      ! copy prefix
      do k = 1, np
        name_array(k,j) = p(k:k)
      end do

      n = min(FIELD_NAME_LEN - np - 1, len_trim(attr_names%t(i)))
      do k = 1, n
        name_array(k+np,j) = attr_names%t(i)(k:k)
      end do
      name_array(np+n+1,j) = C_NULL_CHAR
      if (present(attr_dims)) then
        attr_dims(j) = attr_names%tn(i)
      end if
      j = j + 1
    end do
  end subroutine copy_names_to_array

  !> Set particle attributes for a single subgroup
  subroutine update_particle_subgroup_attributes_and_fields( &
        state, time, dt, subgroup_path, p_list, &
        initial, first_newly_init_part, nb_part_created)
    !> Model state structure
    type(state_type), dimension(:), intent(in) :: state
    !> Current model time
    real, intent(in) :: time
    !> Model timestep
    real, intent(in) :: dt
    !> Option path for the subgroup
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path
    !> Subgroup particle list
    type(detector_linked_list), intent(in) :: p_list
    !> Whether this is the first time attributes are having values filled
    logical, intent(in), optional :: initial
    !> Pointer to the first newly initialised particle
    type(detector_type), pointer, optional :: first_newly_init_part
    !> Number of particles to update
    integer, intent(in), optional :: nb_part_created

    character(len=PYTHON_FUNC_LEN) :: func
    type(detector_type), pointer :: particle
    character(len=OPTION_PATH_LEN) :: attr_key

    ! arrays into which all particle data for the
    ! subgroup is copied
    real, allocatable, dimension(:,:) :: positions
    real, allocatable, dimension(:,:) :: attribute_array
    real, allocatable, dimension(:,:) :: old_attributes
    real, allocatable, dimension(:,:) :: lcoords
    integer, allocatable, dimension(:) :: ele

    logical, allocatable, dimension(:) :: store_old_attr

    character, allocatable, dimension(:,:) :: old_attr_names, field_names, old_field_names
    integer, allocatable, dimension(:) :: old_attr_dims

    real :: constant
    real, allocatable, dimension(:) :: vconstant
    real, allocatable, dimension(:,:) :: tconstant
    integer :: i, j, nparticles, l, m, n, dim, attr_idx, i_single, i_array
    integer :: nscalar, nvector, ntensor
    integer, dimension(3) :: old_attr_counts, field_counts, old_field_counts
    logical :: is_array
    character(len=30) :: value_attr_str

    if (present(nb_part_created)) then
      nparticles = nb_part_created
    else
      nparticles = p_list%length
    end if
    ! return if no particles
    if (nparticles == 0) then
      return
    end if

    ! get attribute sizes from the detector list
    nscalar = size(p_list%attr_names%s)
    nvector = size(p_list%attr_names%v)
    ntensor = size(p_list%attr_names%t)

    ! return if no attributes
    if (nscalar+nvector+ntensor== 0) then
      return
    end if

    ! store all the old attribute names in a contiguous list
    ! for passing through to python functions
    ! we also allocate the array of old attribute dims here
    call copy_names_to_array(p_list%old_attr_names, old_attr_names, old_attr_counts, old_attr_dims)
    call copy_names_to_array(p_list%field_names, field_names, field_counts)
    call copy_names_to_array(p_list%old_field_names, old_field_names, old_field_counts, prefix="old%")

    call get_option("/geometry/dimension", dim)
    allocate(vconstant(dim))
    allocate(tconstant(dim, dim))

    ! allocate space to hold data for all particles in the group
    if (present(first_newly_init_part)) then
      particle => first_newly_init_part
    else
      particle => p_list%first
    end if

    allocate(positions(size(particle%position), nparticles))
    allocate(attribute_array(size(particle%attributes), nparticles))
    allocate(lcoords(size(particle%local_coords), nparticles))
    allocate(ele(nparticles))
    allocate(old_attributes(size(particle%old_attributes), nparticles))

    ! copy the data
    do i = 1, nparticles
      positions(:,i) = particle%position
      lcoords(:,i) = particle%local_coords
      ele(i) = particle%element
      ! copy current attributes in case we loaded from file
      attribute_array(:,i) = particle%attributes
      old_attributes(:,i) = particle%old_attributes
      particle => particle%next
    end do

    attr_idx = 1

    ! calculate new values for all attributes
    i_single = 1
    i_array = 1
    do i = 1, nscalar
      n = p_list%attr_names%sn(i)

      if (n == 0) then
        ! single-valued attribute
        attr_key = trim(subgroup_path) // '/attributes/scalar_attribute['//int2str(i_single-1)//']'
        i_single = i_single + 1
        is_array = .false.
        n = 1
      else
        ! array-valued attribute
        attr_key = trim(subgroup_path) // '/attributes/scalar_attribute_array['//int2str(i_array-1)//']'
        i_array = i_array + 1
        is_array = .true.
      end if

      if (present(initial)) then
        if (initial .and. have_option(trim(attr_key)//'/value_on_spawn')) then
          value_attr_str = "/value_on_spawn"
        end if
      else
        value_attr_str = "/value_on_advection"
      end if

      if (have_option(trim(attr_key)//trim(value_attr_str)//'/constant')) then
        call get_option(trim(attr_key)//trim(value_attr_str)//'/constant', constant)
        attribute_array(attr_idx:attr_idx+n-1,:) = constant
      else if (have_option(trim(attr_key)//trim(value_attr_str)//'/python')) then
        call get_option(trim(attr_key)//trim(value_attr_str)//'/python', func)
        call set_particle_scalar_attribute_from_python( &
             attribute_array(attr_idx:attr_idx+n-1,:), &
             positions(:,:), n, func, time, dt, is_array)
      else if (have_option(trim(attr_key)//trim(value_attr_str)//'/python_fields')) then
        call get_option(trim(attr_key)//trim(value_attr_str)//'/python_fields', func)
        call set_particle_scalar_attribute_from_python_fields( &
             p_list, state, positions(:,:), lcoords(:,:), ele(:), n, &
             attribute_array(attr_idx:attr_idx+n-1,:), &
             old_attr_names, old_attr_counts, old_attr_dims, old_attributes, &
             field_names, field_counts, old_field_names, old_field_counts, &
             func, time, dt, is_array, first_newly_init_part=first_newly_init_part)
      else if (have_option(trim(attr_key)//trim(value_attr_str)//'/from_checkpoint_file')) then
        ! don't do anything, the attribute was already loaded from file
      end if

      attr_idx = attr_idx + n
    end do

    i_single = 1
    i_array = 1
    do i = 1, nvector
      n = p_list%attr_names%vn(i)

      if (n == 0) then
        ! single-valued attribute
         attr_key = trim(subgroup_path) // '/attributes/vector_attribute['//int2str(i_single-1)//']'
         i_single = i_single + 1
         is_array = .false.
         n = 1
      else
        ! array-valued attribute
        attr_key = trim(subgroup_path) // '/attributes/vector_attribute_array['//int2str(i_array-1)//']'
        i_array = i_array + 1
        is_array = .true.
      end if

      if (have_option(trim(attr_key)//trim(value_attr_str)//'/constant')) then
        call get_option(trim(attr_key)//trim(value_attr_str)//'/constant', vconstant)
        ! broadcast vector constant out to all particles
        attribute_array(attr_idx:attr_idx+n*dim-1,:) = spread(vconstant, 2, nparticles)
      else if (have_option(trim(attr_key)//trim(value_attr_str)//'/python')) then
        call get_option(trim(attr_key)//trim(value_attr_str)//'/python', func)
        call set_particle_vector_attribute_from_python( &
             attribute_array(attr_idx:attr_idx+n*dim-1,:), &
             positions(:,:), n, func, time, dt, is_array)
      else if (have_option(trim(attr_key)//trim(value_attr_str)//'/python_fields')) then
        call get_option(trim(attr_key)//trim(value_attr_str)//'/python_fields', func)
        call set_particle_vector_attribute_from_python_fields( &
             p_list, state, positions(:,:), lcoords(:,:), ele(:), n, &
             attribute_array(attr_idx:attr_idx+n*dim-1,:), &
             old_attr_names, old_attr_counts, old_attr_dims, old_attributes, &
             field_names, field_counts, old_field_names, old_field_counts, &
             func, time, dt, is_array, first_newly_init_part=first_newly_init_part)
      else if (have_option(trim(attr_key)//trim(value_attr_str)//'/from_checkpoint_file')) then
        ! don't do anything, the attribute was already loaded from file
      end if

      attr_idx = attr_idx + n*dim
    end do

    i_single = 1
    i_array = 1
    do i = 1, ntensor
      n = p_list%attr_names%tn(i)

      if (n == 0) then
        ! single-valued attribute
        attr_key = trim(subgroup_path) // '/attributes/tensor_attribute['//int2str(i_single-1)//']'
        i_single = i_single + 1
        is_array = .false.
        n = 1
      else
        ! array-valued attribute
        attr_key = trim(subgroup_path) // '/attributes/tensor_attribute_array['//int2str(i_array-1)//']'
        i_array = i_array + 1
        is_array = .true.
      end if

      if (have_option(trim(attr_key)//trim(value_attr_str)//'/constant')) then
        call get_option(trim(attr_key)//trim(value_attr_str)//'/constant', tconstant)
        ! flatten tensor, then broadcast out to all particles
        attribute_array(attr_idx:attr_idx+n*dim**2-1,:) = spread(reshape(tconstant, [dim**2]), 2, nparticles)
      else if (have_option(trim(attr_key)//trim(value_attr_str)//'/python')) then
        call get_option(trim(attr_key)//trim(value_attr_str)//'/python', func)
        call set_particle_tensor_attribute_from_python( &
             attribute_array(attr_idx:attr_idx + n*dim**2 - 1,:), &
             positions(:,:), n, func, time, dt, is_array)
      else if (have_option(trim(attr_key)//trim(value_attr_str)//'/python_fields')) then
        call get_option(trim(attr_key)//trim(value_attr_str)//'/python_fields', func)
        call set_particle_tensor_attribute_from_python_fields( &
             p_list, state, positions(:,:), lcoords(:,:), ele(:), n, &
             attribute_array(attr_idx:attr_idx + n*dim**2 - 1,:), &
             old_attr_names, old_attr_counts, old_attr_dims, old_attributes, &
             field_names, field_counts, old_field_names, old_field_counts, &
             func, time, dt, is_array, first_newly_init_part=first_newly_init_part)
      else if (have_option(trim(attr_key)//trim(value_attr_str)//'/from_checkpoint_file')) then
        ! don't do anything, the attribute was already loaded from file
      end if

      attr_idx = attr_idx + n*dim**2
    end do

    ! Set attribute values and old_attribute values
    if (present(first_newly_init_part)) then
      particle => first_newly_init_part
    else
      particle => p_list%first
    end if

    if (size(particle%old_attributes) == 0) then
      ! no old attributes to store; only store current attributes
      do j = 1, nparticles
        particle%attributes = attribute_array(:,j)
        particle => particle%next
      end do
    else
      ! else, we have to figure out which attributes
      ! need to be stored as old attributes too
      allocate(store_old_attr(size(particle%attributes)))
      attr_idx = 1

      i_single = 1
      i_array = 1
      do i = 1, nscalar
        n = p_list%attr_names%sn(i)
        if (n == 0) then
          store_old_attr(attr_idx) = have_option(trim(subgroup_path) // &
               '/attributes/scalar_attribute['//int2str(i_single-1)//']/value_on_advection/python_fields/store_old_attribute')
          i_single = i_single + 1
          attr_idx = attr_idx + 1
        else
          store_old_attr(attr_idx:attr_idx+n-1) = have_option(trim(subgroup_path) // &
               '/attributes/scalar_attribute_array['//int2str(i_array-1)//']/value_on_advection/python_fields/store_old_attribute')
          i_array = i_array + 1
          attr_idx = attr_idx + n
        end if
      end do

      i_single = 1
      i_array = 1
      do i = 1, nvector
        n = p_list%attr_names%vn(i)
        if (n == 0) then
          store_old_attr(attr_idx:attr_idx+dim-1) = have_option(trim(subgroup_path) // &
               '/attributes/vector_attribute['//int2str(i_single-1)//']/value_on_advection/python_fields/store_old_attribute')
          i_single = i_single + 1
          attr_idx = attr_idx + dim
        else
          store_old_attr(attr_idx:attr_idx + n*dim-1) = have_option(trim(subgroup_path) // &
               '/attributes/vector_attribute_array['//int2str(i_array-1)//']/value_on_advection/python_fields/store_old_attribute')
          i_array = i_array + 1
          attr_idx = attr_idx + n*dim
        end if
      end do

      i_single = 1
      i_array = 1
      do i = 1, ntensor
        n = p_list%attr_names%tn(i)
        if (n == 0) then
          store_old_attr(attr_idx:attr_idx + dim**2 - 1) = have_option(trim(subgroup_path) // &
               '/attributes/tensor_attribute['//int2str(i_single-1)//']/value_on_advection/python_fields/store_old_attribute')
          i_single = i_single + 1
          attr_idx = attr_idx + dim**2
        else
          store_old_attr(attr_idx:attr_idx + n*dim**2 - 1) = have_option(trim(subgroup_path) // &
               '/attributes/tensor_attribute_array['//int2str(i_array-1)//']/value_on_advection/python_fields/store_old_attribute')
          i_array = i_array + 1
          attr_idx = attr_idx + n*dim**2
        end if
      end do

      do j = 1, nparticles
        ! store current attributes as usual
        particle%attributes = attribute_array(:,j)

        ! store the old attributes which are required
        m = 1
        do n = 1, size(particle%attributes)
          if (store_old_attr(n)) then
            particle%old_attributes(m) = particle%attributes(n)
            m = m + 1
          end if
        end do
        particle => particle%next
      end do

      deallocate(store_old_attr)
    end if

    ! update old field values on the particles
    call update_particle_subgroup_fields(state, ele, lcoords, p_list, old_field_counts, &
      first_newly_init_part=first_newly_init_part, nparticles=nparticles)

    deallocate(positions)
    deallocate(lcoords)
    deallocate(ele)
    deallocate(attribute_array)
    deallocate(old_attributes)
    deallocate(old_attr_names)
    deallocate(old_attr_dims)
    deallocate(vconstant)
    deallocate(tconstant)
  end subroutine update_particle_subgroup_attributes_and_fields

  !> Update old values of fields stored on particles
  subroutine update_particle_subgroup_fields( &
      state, ele, lcoords, p_list, counts, first_newly_init_part, nparticles)
    !! Model state structure
    type(state_type), dimension(:), intent(in) :: state
    !! Elements containing particles
    integer, dimension(:), intent(in) :: ele
    !! Local particle coordinates
    real, dimension(:,:), intent(in) :: lcoords
    !! Particle list
    type(detector_linked_list), intent(in) :: p_list
    !! Number of scalar/vector/tensor old fields
    integer, dimension(3), intent(in) :: counts
    !> Pointer to the first newly initialised particle
    type(detector_type), pointer, optional :: first_newly_init_part
    !> Number of particles to update
    integer, intent(in), optional :: nparticles

    integer :: i, dim, nparts
    real, allocatable, dimension(:,:) :: vals
    type(detector_type), pointer :: particle

    call get_option("/geometry/dimension", dim)

    if (present(nparticles)) then
      nparts = nparticles
    else
      nparts = p_list%length
    end if

    allocate(vals(counts(1) + dim*counts(2) + dim**2*counts(3), nparts))

    call evaluate_particle_fields(nparts, state, ele, lcoords, &
         p_list%old_field_names, p_list%old_field_phases, counts, vals, dim)

    ! assign back to particles
    if (present(first_newly_init_part)) then
      particle => first_newly_init_part
    else
      particle => p_list%first
    end if
    do i = 1, nparts
      particle%old_fields = vals(:,i)
      particle => particle%next
    end do

    deallocate(vals)
  end subroutine update_particle_subgroup_fields

  !! Write particle attributes for all groups that should output at the current time
  subroutine write_particles_loop(state, timestep, time)
    !> Model state structure
    type(state_type), dimension(:), intent(in) :: state
    !> Current timestep
    integer, intent(in) :: timestep
    !> Current model time
    real, intent(in) :: time

    integer :: i, k
    integer :: particle_groups, particle_subgroups, list_counter
    character(len=OPTION_PATH_LEN) :: group_path, subgroup_path
    logical :: output_group

    ewrite(1,*) "In write_particles_loop"

    particle_groups = option_count("/particles/particle_group")

    list_counter = 1
    do i = 1, particle_groups
      group_path = "/particles/particle_group["//int2str(i-1)//"]"
      particle_subgroups = option_count(trim(group_path) // "/particle_subgroup")

      output_group = should_output(output_CS(i), time, timestep, trim(group_path) // "/particle_io")
      if (output_group) then
        call update_output_CS(output_CS(i), time)
      else
        ! skip all subgroups
        list_counter = list_counter + particle_subgroups
        cycle
      end if

      do k = 1, particle_subgroups
        subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(k-1)//"]"
        call write_particles_subgroup(state, particle_lists(list_counter), timestep, time, trim(subgroup_path))
          list_counter = list_counter + 1
      end do
    end do
  end subroutine write_particles_loop

  !> Write particle attributes for a given subgroup
  subroutine write_particles_subgroup(state, detector_list, timestep, time, subgroup_path)
    !> Model state structure
    type(state_type), dimension(:), intent(in) :: state
    !> The particle subgroup data structure
    type(detector_linked_list), intent(inout) :: detector_list
    !> Current model timestep (to record in output file)
    integer, intent(in) :: timestep
    !> Current model time (to record in output file)
    real, intent(in) :: time
    !> Path prefix for the subgroup in options
    character(len=*), intent(in) :: subgroup_path

    integer :: dim, i, tot_atts
    integer(kind=8) :: h5_ierror
    real, dimension(:,:), allocatable :: positions, attrib_data
    integer, dimension(:), allocatable :: node_ids, proc_ids
    type(detector_type), pointer :: node

    ewrite(1,*) "In write_particles"

    ! create new step -- create them sequentially so they're easy to iterate
    h5_ierror = h5_setstep(detector_list%h5_id, h5_getnsteps(detector_list%h5_id) + 1)

    ! write time and timestep as step attributes
    h5_ierror = h5_writestepattrib_r8(detector_list%h5_id, "time", [time], int(1, 8))
    h5_ierror = h5_writestepattrib_i8(detector_list%h5_id, "timestep", [int(timestep, 8)], int(1, 8))

    ! set the number of particles this process is going to write
    h5_ierror = h5pt_setnpoints(detector_list%h5_id, int(detector_list%length, 8))

    ! set up arrays to hold all node data
    call get_option("/geometry/dimension", dim)
    tot_atts = detector_list%total_attributes(1)
    allocate(positions(detector_list%length, 3))
    allocate(attrib_data(detector_list%length, tot_atts))
    allocate(node_ids(detector_list%length))
    allocate(proc_ids(detector_list%length))

    node => detector_list%first
    position_loop: do i = 1, detector_list%length
      assert(size(node%position) == dim)
      assert(size(node%attributes) == tot_atts)

      positions(i,1:dim) = node%position(:)
      attrib_data(i,:) = node%attributes(:)
      node_ids(i) = node%id_number
      proc_ids(i) = node%proc_id

      node => node%next
    end do position_loop

    ! write out position
    if (dim >= 1) &
         h5_ierror = h5pt_writedata_r8(detector_list%h5_id, "x", positions(:,1))
    if (dim >= 2) &
         h5_ierror = h5pt_writedata_r8(detector_list%h5_id, "y", positions(:,2))
    if (dim >= 3) then
      h5_ierror = h5pt_writedata_r8(detector_list%h5_id, "z", positions(:,3))
    else
      positions(:,3) = 0.
      h5_ierror = h5pt_writedata_r8(detector_list%h5_id, "z", positions(:,3))
    end if

    h5_ierror = h5pt_writedata_i4(detector_list%h5_id, "id", node_ids(:))
    h5_ierror = h5pt_writedata_i4(detector_list%h5_id, "proc_id", proc_ids(:))

    call write_attrs(detector_list%h5_id, dim, detector_list%attr_names, attrib_data, to_write=detector_list%attr_write)

    h5_ierror = h5_flushstep(detector_list%h5_id)

    deallocate(proc_ids)
    deallocate(node_ids)
    deallocate(attrib_data)
    deallocate(positions)
  end subroutine write_particles_subgroup

  !> Write attributes with given names to an H5Part file
  subroutine write_attrs(h5_id, dim, names, vals, prefix, to_write)
    !> h5 file to write to
    integer(kind=8), intent(in) :: h5_id
    !> spatial dimension
    integer, intent(in) :: dim
    !> attribute names to write to the file
    type(attr_names_type), intent(in) :: names
    !> attribute values, ordered by position/rank as they are on particles
    real, dimension(:,:), intent(in) :: vals
    !> Optional prefix to attribute names
    character(len=*), intent(in), optional :: prefix
    !> Optional control of which attributes to write
    type(attr_write_type), intent(in), optional :: to_write

    integer :: i, j, k, att, ii
    integer(kind=8) :: h5_ierror
    character(len=FIELD_NAME_LEN) :: p
    logical :: write_attr

    p = ""
    if (present(prefix)) p = prefix

    write_attr = .true.

    ! write out attributes -- scalar, vector, tensor
    att = 1
    scalar_attr_loop: do i = 1, size(names%s)
      ! booleans aren't short-circuiting, so we have to stack here
      if (present(to_write)) then
        write_attr = to_write%s(i)
      end if

      if (names%sn(i) == 0) then
        ! single-valued attribute
        if (write_attr) &
             h5_ierror = h5pt_writedata_r8(h5_id, &
             trim(p)//trim(names%s(i)), vals(:,att))
        att = att + 1
      else
        do ii = 1, names%sn(i)
          ! inner loop for array-valued attribute
          if (write_attr) &
               h5_ierror = h5pt_writedata_r8(h5_id, &
               trim(p)//trim(names%s(i))//int2str(ii), vals(:,att))
          att = att + 1
        end do
      end if
    end do scalar_attr_loop

    vector_attr_loop: do i = 1, size(names%v)
      if (present(to_write)) then
        write_attr = to_write%v(i)
      end if

      if (names%vn(i) == 0) then
        do j = 1, dim
          if (write_attr) &
               h5_ierror = h5pt_writedata_r8(h5_id, &
               trim(p)//trim(names%v(i))//"_"//int2str(j-1), vals(:,att))
          att = att + 1
        end do
      else
        do ii = 1, names%vn(i)
          do j = 1, dim
            if (write_attr) &
                 h5_ierror = h5pt_writedata_r8(h5_id, &
                 trim(p)//trim(names%v(i))//int2str(ii)//"_"//int2str(j-1), vals(:,att))
            att = att + 1
          end do
        end do
      end if
    end do vector_attr_loop

    tensor_attr_loop: do i = 1, size(names%t)
      if (present(to_write)) then
        write_attr = to_write%t(i)
      end if

      if (names%tn(i) == 0) then
        do j = 1, dim
          do k = 1, dim
            if (write_attr) &
                 h5_ierror = h5pt_writedata_r8(h5_id, &
                 trim(p)//trim(names%t(i))//"_"//int2str((k-1)*dim + (j-1)), vals(:,att))
            att = att + 1
          end do
        end do
      else
        do ii = 1, names%tn(i)
          do j = 1, dim
            do k = 1, dim
              if (write_attr) &
                   h5_ierror = h5pt_writedata_r8(h5_id, &
                   trim(p)//trim(names%t(i))//int2str(ii)//"_"//int2str((k-1)*dim + (j-1)), vals(:,att))
              att = att + 1
            end do
          end do
        end do
      end if
    end do tensor_attr_loop
  end subroutine write_attrs

  !> Checkpoint all particles, by subgroup
  subroutine checkpoint_particles_loop(state, prefix, postfix, cp_no, number_of_partitions)
    !> Model state structure
    type(state_type), dimension(:), intent(in) :: state
    !> Checkpoint filename prefix
    character(len=*), intent(in) :: prefix
    !> Checkpoint filename postfix (e.g. flredecomp)
    character(len=*), intent(in) :: postfix
    !> Checkpoint number of the simulation
    integer, optional, intent(in) :: cp_no
    !> Only write data for this many processes (flredecomp to more partitions needs this)
    integer, optional, intent(in) :: number_of_partitions

    integer :: i, j, particle_groups, particle_subgroups, list_counter
    character(len=OPTION_PATH_LEN) :: group_path, subgroup_path, subgroup_path_name, name
    type(vector_field), pointer :: xfield

    integer :: output_comm, world_group, output_group, ierr

    ! create a new mpi group for active particles only
    ! otherwise the collectives (and especially file writing) will break

    xfield => extract_vector_field(state(1), "Coordinate")
    if (present(number_of_partitions)) then
      if (getprocno() > number_of_partitions) return

      call mpi_comm_group(MPI_COMM_FEMTOOLS, world_group, ierr)
      call mpi_group_incl(world_group, number_of_partitions, &
           [(i, i=0, number_of_partitions-1)], output_group, ierr)
      call mpi_comm_create_group(MPI_COMM_FEMTOOLS, output_group, 0, output_comm, ierr)
    else
      output_comm = MPI_COMM_FEMTOOLS
    end if

    particle_groups = option_count("/particles/particle_group")

    ewrite(1, *) "Checkpointing particles"

    assert(len_trim(prefix) > 0)

    list_counter = 1
    do i = 1, particle_groups
      group_path = "/particles/particle_group["//int2str(i-1)//"]"
      particle_subgroups = option_count(trim(group_path) // "/particle_subgroup")
      do j = 1, particle_subgroups
        ! set the path to this subgroup, and the path used in update_particle_options
        subgroup_path = trim(group_path) // "/particle_subgroup["//int2str(j-1)//"]"
        subgroup_path_name = trim(group_path) // "/particle_subgroup::"
        call get_option(trim(subgroup_path) // "/name", name)

        ! skip checkpointing this subgroup if we're coming from flredecomp
        ! and the particles weren't loaded from a file
        if (present(number_of_partitions) .and. &
             .not. have_option(trim(subgroup_path) // "/initial_position/from_file")) then
          cycle
        end if

        !Ensure all particles are local before checkpointing
        if (present(number_of_partitions)) then
           !Don't call distribute detectors if in flredecomp
        else
           call distribute_detectors(state(1),particle_lists(list_counter),positions = xfield)
        end if
        !Checkpoint particle group
        call checkpoint_particles_subgroup(state, prefix, postfix, cp_no, particle_lists(list_counter), &
             name, subgroup_path, subgroup_path_name, output_comm)
        list_counter = list_counter + 1
      end do
    end do

     ! clean up mpi structures
     if (present(number_of_partitions)) then
       call mpi_comm_free(output_comm, ierr)
       call mpi_group_free(output_group, ierr)
     end if
  end subroutine checkpoint_particles_loop

  !> Checkpoint a single particle subgroup
  subroutine checkpoint_particles_subgroup(state, prefix, postfix, cp_no, particle_list, &
       name, subgroup_path, subgroup_path_name, output_comm)
    !> Model state structure
    type(state_type), dimension(:), intent(in) :: state
    !> Checkpoint filename prefix
    character(len=*), intent(in) :: prefix
    !> Checkpoint filename postfix
    character(len=*), intent(in) :: postfix
    !> Checkpoint number of the simulation
    integer, optional, intent(in) :: cp_no
    !> Particle list for the subgroup
    type(detector_linked_list), intent(inout) :: particle_list
    !> Particle subgroup name
    character(len=*), intent(in) :: name
    !> Option subgroup path, and the prefix used for updating options
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path, subgroup_path_name
    !> MPI communicator to use for output/collectives
    integer, optional, intent(in) :: output_comm

    character(len=OPTION_PATH_LEN) :: particles_cp_filename
    integer :: i, dim, tot_attrs, tot_old_attrs, tot_old_fields
    real, dimension(:,:), allocatable :: positions, attr_data, old_attr_data, old_field_data
    integer(kind=8) :: h5_id, h5_prop, h5_ierror
    integer(kind=8), dimension(:), allocatable :: npoints
    integer, dimension(:), allocatable :: node_ids, proc_ids
    type(detector_type), pointer :: node

    integer :: comm, commsize, ierr
    comm = MPI_COMM_FEMTOOLS
    if (present(output_comm)) comm = output_comm

    ! we store the number of points per process, so gather them
    ! this is because h5part attributes must be agreed upon by all,
    ! so every process needs to know the size for every other process
    ! we also store the proc_part_count to ensure newly spawned
    ! particles after checkpointing remain unique
    call mpi_comm_size(comm, commsize, ierr)
    allocate(npoints(commsize*2))
    call mpi_allgather([int(particle_list%length, 8), int(particle_list%proc_part_count, 8)], &
         2, MPI_INTEGER8, &
         npoints, 2, MPI_INTEGER8, comm, ierr)

    ! construct a new particle checkpoint filename
    particles_cp_filename = trim(prefix)
    if(present(cp_no)) particles_cp_filename = trim(particles_cp_filename) // "_" // int2str(cp_no)
    particles_cp_filename = trim(particles_cp_filename) // "_" // trim(postfix)
    particles_cp_filename = trim(particles_cp_filename) // "_particles." // trim(name) // ".h5part"

    ! restrict h5 IO to the specified communicator
    h5_prop = h5_createprop_file()
    h5_ierror = h5_setprop_file_mpio_collective(h5_prop, comm)

    ! open output file
    h5_id = h5_openfile(trim(particles_cp_filename), H5_O_WRONLY, h5_prop)
    h5_ierror = h5_closeprop(h5_prop)
    ! write out number of points per process
    h5_ierror = h5_writefileattrib_i8(h5_id, "npoints", npoints(1::2), int(commsize, 8))
    ! write out per-processor particle count, for unique spawning
    h5_ierror = h5_writefileattrib_i8(h5_id, "part_counter", npoints(2::2), int(commsize, 8))
    ! write data in the first step
    h5_ierror = h5_setstep(h5_id, int(1, 8))
    ! the number of points this process is writing
    h5_ierror = h5pt_setnpoints(h5_id, int(particle_list%length, 8))

    ! get dimension of particle positions
    call get_option("/geometry/dimension", dim)

    tot_attrs = particle_list%total_attributes(1)
    tot_old_attrs = particle_list%total_attributes(2)
    tot_old_fields = particle_list%total_attributes(3)

    ! allocate arrays for node data
    allocate(positions(particle_list%length, dim))
    allocate(node_ids(particle_list%length))
    allocate(proc_ids(particle_list%length))
    allocate(attr_data(particle_list%length, tot_attrs))
    allocate(old_attr_data(particle_list%length, tot_old_attrs))
    allocate(old_field_data(particle_list%length, tot_old_fields))

    ! gather data off all particles
    node => particle_list%first
    positionloop_cp: do i = 1, particle_list%length
      ! collect positions
      assert(size(node%position) == dim)

      positions(i,:) = node%position(:)
      if (tot_attrs /= 0) &
           attr_data(i,:) = node%attributes(:)
      if (tot_old_attrs /= 0) &
           old_attr_data(i,:) = node%old_attributes(:)
      if (tot_old_fields /= 0) &
           old_field_data(i,:) = node%old_fields(:)

      ! collect node ids
      node_ids(i) = node%id_number
      proc_ids(i) = node%proc_id

      node => node%next
    end do positionloop_cp

    ! write out positions and ids
    if (dim >= 1) &
         h5_ierror = h5pt_writedata_r8(h5_id, "x", positions(:,1))
    if (dim >= 2) &
         h5_ierror = h5pt_writedata_r8(h5_id, "y", positions(:,2))
    if (dim >= 3) &
         h5_ierror = h5pt_writedata_r8(h5_id, "z", positions(:,3))
    h5_ierror = h5pt_writedata_i4(h5_id, "id", node_ids(:))
    h5_ierror = h5pt_writedata_i4(h5_id, "proc_id", proc_ids(:))

    call write_attrs(h5_id, dim, particle_list%attr_names, attr_data)
    call write_attrs(h5_id, dim, particle_list%old_attr_names, old_attr_data)
    call write_attrs(h5_id, dim, particle_list%old_field_names, old_field_data, prefix="old%")

    ! update schema file to read this subgroup from the checkpoint file
    call update_particle_subgroup_options(trim(particles_cp_filename), particle_list, name, &
         tot_attrs, subgroup_path_name)

    deallocate(old_field_data)
    deallocate(old_attr_data)
    deallocate(node_ids)
    deallocate(proc_ids)
    deallocate(attr_data)
    deallocate(positions)
    deallocate(npoints)

    h5_ierror = h5_closefile(h5_id)
  end subroutine checkpoint_particles_subgroup

  subroutine update_particle_subgroup_options(filename, particle_list, name, tot_atts, subgroup_path_name)
    !! Updates the initial options of particles in the schema file for reinitialization after checkpointing.
    !! Updates schema options for the initial number of particles and their initial positions.

    character(len = *), intent(in) :: filename
    character(len = *), intent(in) :: name

    type(detector_linked_list), intent(inout) :: particle_list
    integer, intent(in) :: tot_atts
    character(len=OPTION_PATH_LEN), intent(in) :: subgroup_path_name

    integer :: num_particles, j, stat
    logical :: particles_s, particles_v, particles_t

    character(len = 254) :: temp_string

    num_particles = particle_list%total_num_det

    temp_string=name

    ewrite(1,*) 'In update_particles_options'
    ewrite(1,*) temp_string

    call delete_option(trim(subgroup_path_name) // trim(temp_string) // "/initial_position")

    call set_option_attribute(trim(subgroup_path_name) // trim(temp_string) // "/initial_position/from_file/file_name", trim(filename), stat)
    call set_option(trim(subgroup_path_name) // trim(temp_string) // "/initial_position/from_file/number_of_particles", num_particles, stat = stat)

    if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
       FLAbort("Failed to set particles options filename when checkpointing particles with option path " // "/particles/particle_array::" // trim(temp_string))
    end if

    do j = 1, tot_atts
      particles_s = have_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/scalar_attribute["//int2str(j-1)//"]/constant")
      particles_v = have_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/vector_attribute["//int2str(j-1)//"]/constant")
      particles_t = have_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/tensor_attribute["//int2str(j-1)//"]/constant")
      if (particles_s) then
        call delete_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/scalar_attribute["//int2str(j-1)//"]/constant")
        call set_option_attribute(trim(subgroup_path_name) // trim(temp_string) // "/attributes/scalar_attribute["//int2str(j-1)// &
             "]/from_checkpoint_file/file_name", trim(filename) // "." // trim(temp_string), stat)
          if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
             FLAbort("Failed to set scalar field particles filename when checkpointing")
          end if
       else if (particles_v) then
          call delete_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/vector_attribute["//int2str(j-1)//"]/constant")
          call set_option_attribute(trim(subgroup_path_name) // trim(temp_string) // "/attributes/vector_attribute["//int2str(j-1)// &
               "]/from_checkpoint_file/file_name", trim(filename) // "." // trim(temp_string), stat)
          if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
             FLAbort("Failed to set vector field particles filename when checkpointing")
          end if
       else if (particles_t) then
          call delete_option(trim(subgroup_path_name) // trim(temp_string) // "/attributes/tensor_attribute["//int2str(j-1)//"]/constant")
          call set_option_attribute(trim(subgroup_path_name) // trim(temp_string) // "/attributes/tensor_attribute["//int2str(j-1)// &
               "]/from_checkpoint_file/file_name", trim(filename) // "." // trim(temp_string), stat)
          if(stat /= SPUD_NO_ERROR .and. stat /= SPUD_NEW_KEY_WARNING .and. stat /= SPUD_ATTR_SET_FAILED_WARNING) then
             FLAbort("Failed to set tensor field particles filename when checkpointing")
          end if
       end if
    end do

  end subroutine update_particle_subgroup_options

  subroutine get_particles(p_array, p_allocated)
    !Send particle arrays to another routine

    type(detector_linked_list), allocatable, dimension(:), intent(out) :: p_array
    integer, intent(out) :: p_allocated

    integer :: i, particle_groups

    particle_groups = option_count("/particles/particle_group")
    if (particle_groups==0) then
       FLAbort("No particle groups exist")
       return
    end if

    if (allocated(particle_lists)) then
       p_allocated = 1
       allocate(p_array(size(particle_lists)))
       do i = 1,size(particle_lists)
          p_array(i) = particle_lists(i)
       end do
    else
       p_allocated = 0
    end if

  end subroutine get_particles

  subroutine get_particle_arrays(lgroup, group_arrays, group_attribute, att_n, lattribute)
    !Read in a particle group and attribute name or particle subgroup, send back numbers of particle arrays and particle attribute

    character(len=OPTION_PATH_LEN), intent(in) :: lgroup
    character(len=OPTION_PATH_LEN), optional, intent(in) :: lattribute
    integer, allocatable, dimension(:), intent(out) :: group_arrays
    integer, optional, intent(out) :: group_attribute
    integer, optional, intent(in) :: att_n

    character(len=OPTION_PATH_LEN) :: group_name, attribute_name, subgroup_name
    integer :: particle_groups, array_counter, particle_subgroups, particle_attributes
    integer :: i, j, k, l

    logical :: found_attribute

    particle_groups = option_count("/particles/particle_group")

    found_attribute = .false.
    array_counter = 0
    do i = 1, particle_groups
       call get_option("/particles/particle_group["//int2str(i-1)//"]/name", group_name)
       particle_subgroups = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup")
       if (trim(group_name)==trim(lgroup)) then
          allocate(group_arrays(particle_subgroups))
          if (present(lattribute)) then
             if (att_n==0) then
                particle_attributes = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup["//int2str(0)//"]/attributes/scalar_attribute")
                do k = 1, particle_attributes
                   call get_option("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup["//int2str(0)// &
                        "]/attributes/scalar_attribute["//int2str(k-1)//"]/name", attribute_name)
                   if (trim(attribute_name)==trim(lattribute)) then
                      found_attribute = .true.
                      group_attribute = k
                   end if
                end do
                if (found_attribute.eqv..false.) then
                   FLExit("Could not find particle attribute "//trim(lattribute)//" in particle group "//trim(lgroup)//". Check attribute is a scalar.")
                end if
             else
                particle_attributes = option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup["//int2str(0)//"]/attributes/scalar_attribute_array")
                l = 0
                group_attribute = 0
                do k = 1, particle_attributes
                   call get_option("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup["//int2str(0)// &
                        "]/attributes/scalar_attribute_array["//int2str(k-1)//"]/name", attribute_name)
                   call get_option("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup["//int2str(0)// &
                        "]/attributes/scalar_attribute_array["//int2str(k-1)//"]/dimension", l)
                   if (trim(attribute_name)==trim(lattribute)) then
                      found_attribute = .true.
                      group_attribute = group_attribute + att_n
                   else
                      group_attribute = group_attribute + l
                   end if
                end do
                if (found_attribute.eqv..false.) then
                   FLExit("Could not find particle attribute "//trim(lattribute)//" in particle group "//trim(lgroup)//". Check attribute is a scalar.")
                end if
                group_attribute = group_attribute + option_count("/particles/particle_group["//int2str(i-1)//"]/particle_subgroup["//int2str(0)//"]/attributes/scalar_attribute")
             end if
          end if
          j=1
          do l = array_counter+1, array_counter+particle_subgroups
             group_arrays(j) = l
             j=j+1
          end do
          return
       end if
       array_counter = array_counter + particle_subgroups
    end do
    FLExit("Could not find particle group "//trim(lgroup))
  end subroutine get_particle_arrays

  subroutine update_list_lengths(list_num)
    integer, intent(in) :: list_num

    particle_lists(list_num)%total_num_det = particle_lists(list_num)%total_num_det + 1
    particle_lists(list_num)%length = particle_lists(list_num)%length + 1

  end subroutine update_list_lengths
  subroutine destroy_particles()
    type(detector_linked_list), pointer :: del_particle_lists
    integer :: i, particle_groups
    integer(kind=8) :: h5_ierror

    if (allocated(particle_lists)) then
      ! gracefully clean up output files and deallocate all particle arrays (detector lists)
      particle_groups = size(particle_lists)
      do i = 1, particle_groups
         if (particle_lists(i)%h5_id /= -1) h5_ierror = h5_closefile(particle_lists(i)%h5_id)
         del_particle_lists => particle_lists(i)
         call deallocate(del_particle_lists)
      enddo
    end if

    if (allocated(output_CS)) then
      deallocate(output_CS)
    end if
  end subroutine destroy_particles

end module particles
