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

module detector_tools
  use spud
  use fldebug
  use iso_c_binding, only: C_NULL_CHAR
  use global_parameters, only: OPTION_PATH_LEN, FIELD_NAME_LEN
  use futils, only: int2str
  use elements, only: local_coord_count
  use embed_python
  use integer_hash_table_module
  use parallel_tools
  use parallel_fields, only: element_owned
  use transform_elements
  use fields
  use state_module, only: state_type, extract_scalar_field, aliased, &
       & extract_vector_field, extract_tensor_field
  use field_options
  use detector_data_types
  use pickers

  implicit none

  private

  public :: insert, allocate, deallocate, copy, move, move_all, remove, &
            delete, delete_all, pack_detector, unpack_detector, &
            detector_value, set_detector_coords_from_python, &
            detector_buffer_size, set_particle_scalar_attribute_from_python, &
            set_particle_scalar_attribute_from_python_fields, &
            set_particle_vector_attribute_from_python, &
            set_particle_vector_attribute_from_python_fields, &
            set_particle_tensor_attribute_from_python, &
            set_particle_tensor_attribute_from_python_fields, &
            evaluate_particle_fields, temp_list_insert, &
            temp_list_deallocate, temp_list_remove

  interface insert
     module procedure insert_into_detector_list
  end interface

  ! Removes detector from a list without deallocating it
  interface remove
     module procedure remove_detector_from_list
  end interface

  ! Removes detector from a list and deallocates it
  interface delete
     module procedure delete_detector
  end interface

  ! Deletes all detectors from a given list
  interface delete_all
     module procedure delete_all_detectors
  end interface

  ! Move detector from one list to another
  interface move
     module procedure move_detector
  end interface

  ! Move all detectors in a list from one that list to another
  interface move_all
     module procedure move_all_detectors
  end interface

  interface copy
     module procedure detector_copy
  end interface

  interface allocate
     module procedure detector_allocate_from_params, detector_allocate_from_detector
  end interface

  interface deallocate
     module procedure detector_deallocate, detector_list_deallocate
  end interface

  ! Evaluate field at the location of the detector.
  interface detector_value
     module procedure detector_value_scalar, detector_value_vector
  end interface

contains

  subroutine detector_allocate_from_params(new_detector, ndims, local_coord_count, attribute_size)
    type(detector_type),  pointer, intent(out) :: new_detector
    integer, intent(in) :: ndims, local_coord_count
    integer, dimension(3), optional, intent(in) :: attribute_size !array to hold size of attributes

    assert(.not. associated(new_detector))

    ! allocate the memory for the new detector
    if (.not. associated(new_detector)) then
       allocate(new_detector)
    end if
    allocate(new_detector%position(ndims))
    allocate(new_detector%local_coords(local_coord_count))
    if (present(attribute_size)) then
       allocate(new_detector%attributes(attribute_size(1)))
       allocate(new_detector%old_attributes(attribute_size(2)))
       allocate(new_detector%old_fields(attribute_size(3)))
    else
       ! match the behaviour of create_single_detector, with empty attribute arrays
       allocate(new_detector%attributes(0))
       allocate(new_detector%old_attributes(0))
       allocate(new_detector%old_fields(0))
    end if

    assert(associated(new_detector))

  end subroutine detector_allocate_from_params

  subroutine detector_allocate_from_detector(new_detector, old_detector)
    type(detector_type), pointer, intent(in) :: old_detector
    type(detector_type),  pointer, intent(out) :: new_detector

    integer :: ndims, local_coord_count
    integer, dimension(3) :: attribute_size !array to hold size of attributes

    ndims = size(old_detector%position)
    local_coord_count = size(old_detector%local_coords)
    attribute_size(1) = size(old_detector%attributes)
    attribute_size(2) = size(old_detector%old_attributes)
    attribute_size(3) = size(old_detector%old_fields)

    ! allocate the memory for the new detector
    call detector_allocate_from_params(new_detector, ndims, local_coord_count, attribute_size)

  end subroutine detector_allocate_from_detector

  subroutine detector_deallocate(detector)
    type(detector_type), pointer :: detector

    if(associated(detector)) then
       if(allocated(detector%local_coords)) then
          deallocate(detector%local_coords)
       end if
       if(allocated(detector%position)) then
          deallocate(detector%position)
       end if
       if(allocated(detector%k)) then
          deallocate(detector%k)
       end if
       if(allocated(detector%update_vector)) then
          deallocate(detector%update_vector)
       end if
       if(allocated(detector%attributes)) then
          deallocate(detector%attributes)
       end if
       if(allocated(detector%old_attributes)) then
          deallocate(detector%old_attributes)
       end if
       if(allocated(detector%old_fields)) then
          deallocate(detector%old_fields)
       end if
       detector%next => null()
       detector%previous => null()
       detector%temp_next => null()
       detector%temp_previous => null()
       deallocate(detector)
    end if
    detector => null()

  end subroutine detector_deallocate

  subroutine detector_list_deallocate(detector_list)
    type(detector_linked_list), pointer :: detector_list

    type(rk_gs_parameters), pointer :: parameters

    ! Delete detectors
    if (detector_list%length > 0) then
       call delete_all(detector_list)
    end if

    ! Deallocate list information
    if (allocated(detector_list%sfield_list)) then
       deallocate(detector_list%sfield_list)
    end if
    if (allocated(detector_list%vfield_list)) then
       deallocate(detector_list%vfield_list)
    end if

    ! Deallocate move_parameters
    parameters => detector_list%move_parameters
    if (associated(parameters)) then
       if (allocated(parameters%timestep_weights)) then
          deallocate(parameters%timestep_weights)
       end if
       if (allocated(parameters%stage_matrix)) then
          deallocate(parameters%stage_matrix)
       end if
    end if

  end subroutine detector_list_deallocate

  subroutine temp_list_deallocate(detector_list)
    ! Removes all detectors from the temporary list

    type(detector_linked_list), pointer :: detector_list
    type(detector_type), pointer :: detector
    type(detector_type), pointer :: temp_detector

    if (detector_list%length==0) return

    detector => detector_list%first
    do while (associated(detector))
       temp_detector => detector
       detector => detector%temp_next
       call temp_list_remove(temp_detector, detector_list)
    end do

  end subroutine temp_list_deallocate

  subroutine detector_copy(new_detector, old_detector)
    ! Copies all the information from the old detector to
    ! the new detector
    type(detector_type), pointer, intent(in) :: old_detector
    type(detector_type),  pointer :: new_detector

    new_detector%position = old_detector%position
    new_detector%element = old_detector%element
    new_detector%id_number = old_detector%id_number
    new_detector%local_coords=old_detector%local_coords

  end subroutine detector_copy

  subroutine insert_into_detector_list(detector, current_list)
    ! Inserts detector at the end of a list
    type(detector_linked_list), intent(inout) :: current_list
    type(detector_type), pointer :: detector

    if (current_list%length == 0) then
       current_list%first => detector
       current_list%last => detector
       current_list%first%previous => null()
       current_list%last%next => null()
       current_list%length = 1
    else
       detector%previous => current_list%last
       current_list%last%next => detector
       current_list%last => detector
       current_list%last%next => null()
       current_list%length = current_list%length+1
    end if

  end subroutine insert_into_detector_list

  subroutine temp_list_insert(detector, current_list)
    ! Inserts detector at the end of a temporary detector list
    ! i.e. a list of detectors that is linked via the temp_next
    ! and temp_previous pointers
    type(detector_linked_list), intent(inout) :: current_list
    type(detector_type), pointer :: detector

    if (current_list%length == 0) then
       current_list%first => detector
       current_list%last => detector
       current_list%first%temp_previous => null()
       current_list%last%temp_next => null()
       current_list%length = 1
    else
       detector%temp_previous => current_list%last
       current_list%last%temp_next => detector
       current_list%last => detector
       current_list%last%temp_next => null()
       current_list%length = current_list%length+1
    end if

  end subroutine temp_list_insert

  subroutine remove_detector_from_list(detector, detector_list)
    !! Removes the detector from the list,
    !! but does not deallocated it
    type(detector_linked_list), intent(inout) :: detector_list
    type(detector_type), pointer :: detector

    if (detector_list%length==1) then
       detector_list%first => null()
       detector_list%last => null()
    else
       if (associated(detector%previous)) then
          detector%previous%next => detector%next
       else
          detector_list%first => detector%next
       end if

       if (associated(detector%next)) then
          detector%next%previous => detector%previous
       else
          detector_list%last => detector%previous
       end if
    end if

    detector_list%length = detector_list%length-1

  end subroutine remove_detector_from_list

  subroutine temp_list_remove(detector, detector_list)
    !! Removes the detector from the temporary detector list,
    !! but does not deallocated it
    type(detector_linked_list), intent(inout) :: detector_list
    type(detector_type), pointer :: detector

    if (detector_list%length==1) then
       detector_list%first => null()
       detector_list%last => null()
    else
       if (associated(detector%temp_previous)) then
          detector%temp_previous%temp_next => detector%temp_next
       else
          detector_list%first => detector%temp_next
       end if

       if (associated(detector%temp_next)) then
          detector%temp_next%temp_previous => detector%temp_previous
       else
          detector_list%last => detector%temp_previous
       end if
    end if
    detector_list%length = detector_list%length-1

  end subroutine temp_list_remove

  subroutine delete_detector(detector, detector_list)
    ! Removes and deallocates the given detector
    ! and outputs the next detector in the list as detector
    type(detector_type), pointer :: detector
    type(detector_linked_list), intent(inout), optional :: detector_list

    type(detector_type), pointer :: temp_detector

    if (present(detector_list)) then
       temp_detector => detector
       detector => detector%next
       call remove(temp_detector, detector_list)
       call deallocate(temp_detector)
    else
       call deallocate(detector)
    end if

  end subroutine delete_detector

  subroutine move_detector(detector, from_list, to_list)
    ! Move detector from one list to the other
    type(detector_linked_list), intent(inout) :: from_list
    type(detector_type), pointer :: detector
    type(detector_linked_list), intent(inout) :: to_list

    call remove(detector, from_list)
    call insert(detector, to_list)

  end subroutine move_detector

  subroutine move_all_detectors(from_list,to_list)
    ! Move all detectors from one list to the other
    type(detector_linked_list), intent(inout) :: from_list
    type(detector_linked_list), intent(inout) :: to_list
    type(detector_type), pointer :: detector

    do while (associated(from_list%first))
       detector => from_list%first
       call move(detector, from_list, to_list)
    end do

  end subroutine move_all_detectors

  subroutine delete_all_detectors(detector_list)
    ! Remove and deallocate all detectors in a list
    type(detector_linked_list), intent(inout) :: detector_list
    type(detector_type), pointer :: detector

    detector => detector_list%first
    do while (associated(detector))
       call delete(detector,detector_list)
    end do

  end subroutine delete_all_detectors

  function detector_buffer_size(ndims, have_update_vector, nstages, attribute_size)
    ! Returns the number of reals we need to pack a detector
    integer, intent(in) :: ndims
    logical, intent(in) :: have_update_vector
    integer, intent(in), optional :: nstages
    integer :: detector_buffer_size, det_params
    integer, dimension(3), optional, intent(in) :: attribute_size !array to hold size of attributes

    det_params = 3 !size of basic detector fields: detector element, id_number and proc_id

    ! common to everything is a position + basic fields
    detector_buffer_size = ndims + det_params

    if (present(attribute_size)) then
      detector_buffer_size = detector_buffer_size + sum(attribute_size)
    end if

    if (have_update_vector) then
      ! update vector adds ndims + nstages*ndims
      detector_buffer_size = detector_buffer_size + (nstages + 1)*ndims
    else
      ! otherwise, there's a list id
      detector_buffer_size = detector_buffer_size + 1
    end if

  end function detector_buffer_size

  subroutine pack_detector(detector,buff,ndims,nstages, attribute_size_in)
    ! Packs (serialises) detector into buff
    ! Basic fields are: element, position, id_number and type
    ! If nstages is given, the detector is still moving
    ! and we also pack update_vector and k
    type(detector_type), pointer, intent(in) :: detector
    real, dimension(:), intent(out) :: buff
    integer, intent(in) :: ndims
    integer, intent(in), optional :: nstages
    integer, dimension(3), optional, intent(in) :: attribute_size_in

    integer :: det_params, buf_pos
    integer, dimension(3) :: attribute_size

    assert(size(detector%position) == ndims)

    !Set size of basic detector fields: detector element, id_number, proc_id
    det_params = 3
    attribute_size(:) = 0

    if (present(attribute_size_in)) then
      attribute_size = attribute_size_in
    end if

    ! ensure buffer is big enough to receive this detector
    assert(size(buff) >= ndims + det_params + sum(attribute_size))

    ! Basic fields: ndims+det_params
    buff(1:ndims) = detector%position
    buff(ndims+1) = detector%element
    buff(ndims+2) = detector%id_number
    buff(ndims+3) = detector%proc_id

    buf_pos = ndims + det_params

    if (attribute_size(1) /= 0) then
       buff(buf_pos + 1:buf_pos + attribute_size(1)) = detector%attributes
       buf_pos = buf_pos + attribute_size(1)
    end if
    if (attribute_size(2) /= 0) then
       buff(buf_pos + 1:buf_pos + attribute_size(2)) = detector%old_attributes
       buf_pos = buf_pos + attribute_size(2)
    end if
    if (attribute_size(3) /= 0) then
       buff(buf_pos + 1:buf_pos + attribute_size(3)) = detector%old_fields
       buf_pos = buf_pos + attribute_size(3)
    end if

    ! Lagrangian advection fields: (nstages+1)*ndims
    if (present(nstages)) then
       assert(size(buff) == (nstages+2)*ndims + det_params + sum(attribute_size))
       assert(allocated(detector%update_vector))
       assert(allocated(detector%k))

       buff(buf_pos + 1:buf_pos + ndims) = detector%update_vector
       buf_pos = buf_pos + ndims

       buff(buf_pos + 1:buf_pos + nstages*ndims) = reshape(detector%k, (/nstages*ndims/))
    else
       assert(size(buff) == ndims + det_params + sum(attribute_size) + 1)
       buff(buf_pos + 1) = detector%list_id
    end if
  end subroutine pack_detector

  subroutine unpack_detector(detector,buff,ndims,global_to_local,coordinates,nstages,attribute_size_in)
    ! Unpacks the detector from buff and fills in the blanks
    type(detector_type), pointer :: detector
    real, dimension(:), intent(in) :: buff
    integer, intent(in) :: ndims
    type(integer_hash_table), intent(in), optional :: global_to_local
    type(vector_field), intent(in), optional :: coordinates
    integer, intent(in), optional :: nstages
    integer, dimension(3), optional, intent(in) :: attribute_size_in

    integer :: det_params, buf_pos
    integer, dimension(3) :: attribute_size
    !Set size of basic detector fields, being detector element, id_number
    det_params = 3
    attribute_size(:) = 0

    ! we default to assuming there are no attributes,
    ! but this can be overridden by the caller
    if (present(attribute_size_in)) then
      attribute_size = attribute_size_in
    end if

    ! allocate some arrays that we might not have
    ! set up beforehand
    if (.not. allocated(detector%position)) then
       allocate(detector%position(ndims))
    end if

    if (.not. allocated(detector%attributes)) then
       allocate(detector%attributes(attribute_size(1)))
       allocate(detector%old_attributes(attribute_size(2)))
       allocate(detector%old_fields(attribute_size(3)))
    end if

    ! Basic fields: ndims+4
    detector%position = buff(1:ndims)
    detector%element = buff(ndims+1)
    detector%id_number = buff(ndims+2)
    detector%proc_id = buff(ndims+3)

    buf_pos = ndims + det_params

    ! unpack attributes if necessary
    if (attribute_size(1) /= 0) then
       detector%attributes = buff(buf_pos + 1:buf_pos + attribute_size(1))
       buf_pos = buf_pos + attribute_size(1)
    end if
    if (attribute_size(2) /= 0) then
       detector%old_attributes = buff(buf_pos + 1:buf_pos + attribute_size(2))
       buf_pos = buf_pos + attribute_size(2)
    end if
    if (attribute_size(3) /= 0) then
       detector%old_fields = buff(buf_pos + 1:buf_pos + attribute_size(3))
       buf_pos = buf_pos + attribute_size(3)
    end if

    ! Reconstruct element number if global-to-local mapping is given
    if (present(global_to_local)) then
       assert(has_key(global_to_local, detector%element))
       detector%element=fetch(global_to_local,detector%element)

    ! Update local coordinates if coordinate field is given
       if (present(coordinates)) then
          if (.not. allocated(detector%local_coords)) then
             allocate(detector%local_coords(local_coord_count(ele_shape(coordinates,1))))
          end if
          detector%local_coords=local_coords(coordinates,detector%element,detector%position)
       end if
    end if

    ! Lagrangian advection fields: (nstages+1)*ndims
    if (present(nstages)) then
       assert(size(buff) == (nstages+2)*ndims + det_params + sum(attribute_size))

       ! update_vector, dimension(ndim)
       if (.not. allocated(detector%update_vector)) then
          allocate(detector%update_vector(ndims))
       end if
       detector%update_vector = buff(buf_pos + 1:buf_pos + ndims)
       buf_pos = buf_pos + ndims

       ! k, dimension(nstages:ndim)
       if (.not. allocated(detector%k)) then
          allocate(detector%k(nstages,ndims))
       end if
       detector%k = reshape(buff(buf_pos + 1:buf_pos + nstages*ndims), &
          (/nstages,ndims/))

       ! If update_vector still exists, we're not done moving
       detector%search_complete=.false.
    else
       assert(size(buff) == ndims + det_params + sum(attribute_size) + 1)

       detector%list_id = buff(buf_pos + 1)
       detector%search_complete = .true.
    end if

  end subroutine unpack_detector

  function detector_value_scalar(sfield, detector) result(value)
    !!< Evaluate field at the location of the detector.
    real :: value
    type(scalar_field), intent(in) :: sfield
    type(detector_type), intent(in) :: detector

    assert(detector%element>0)
    value = eval_field(detector%element, sfield, detector%local_coords)

  end function detector_value_scalar

  function detector_value_vector(vfield, detector) result(value)
    !!< Evaluate field at the location of the detector.
    type(vector_field), intent(in) :: vfield
    type(detector_type), intent(in) :: detector
    real, dimension(vfield%dim) :: value

    assert(detector%element>0)
    value = eval_field(detector%element, vfield, detector%local_coords)

  end function detector_value_vector

  subroutine set_detector_coords_from_python(values, ndete, func, time)
    !!< Given a list of positions and a time, evaluate the python function
    !!< specified in the string func at those points.
    real, dimension(:,:), target, intent(inout) :: values
    !! Func may contain any python at all but the following function must
    !! be defined:
    !!  def val(t)
    !! where t is the time. The result must be a float.
    character(len=*), intent(in) :: func
    real :: time

    real, dimension(:), pointer :: lvx,lvy,lvz
    real, dimension(0), target :: zero
    integer :: stat, dim, ndete

    call get_option("/geometry/dimension",dim)

    lvx=>values(1,:)
    lvy=>zero
    lvz=>zero
    if(dim>1) then
       lvy=>values(2,:)
       if(dim>2) then
          lvz => values(3,:)
       end if
    end if

    call set_detectors_from_python(func, len(func), dim, &
         ndete, time, dim,                               &
         lvx, lvy, lvz, stat)

    if (stat/=0) then
      ! FIX ME - Does not print
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1 , *) trim(func)
      FLExit("Dying")
    end if

  end subroutine set_detector_coords_from_python

  !> Evaluate a set of fields on particles
  subroutine evaluate_particle_fields(npart, state, ele, lcoords, names, phases, counts, vals, dim)
    !! Number of particles
    integer, intent(in) :: npart
    !! Model state structure
    type(state_type), dimension(:), intent(in) :: state
    !! Elements containing particles
    integer, dimension(:), intent(in) :: ele
    !! Local coordinates of particles in elements
    real, dimension(:,:), intent(in) :: lcoords
    !! Names of fields to be evaluated
    type(attr_names_type), intent(in) :: names
    !! Phases of each field
    type(field_phase_type), intent(in) :: phases
    !! Scalar/vector/tensor counts of each field
    integer, dimension(3), intent(in) :: counts
    !! Output array
    real, dimension(:,:), intent(out) :: vals
    !! Geometric dimension
    integer, intent(in) :: dim

    integer :: i, j
    integer :: field_idx, phase
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield

    field_idx = 1

    ! scalar fields
    do i = 1, counts(1)
      phase = phases%s(i)
      sfield => extract_scalar_field(state(phase), names%s(i))

      do j = 1, npart
        vals(field_idx,j) = eval_field(ele(j), sfield, lcoords(:,j))
      end do

      field_idx = field_idx + 1
    end do

    ! vector fields
    do i = 1, counts(2)
      phase = phases%v(i)
      vfield => extract_vector_field(state(phase), names%v(i))

      do j = 1, npart
        vals(field_idx:field_idx+dim-1,j) = eval_field(ele(j), vfield, lcoords(:,j))
      end do

      field_idx = field_idx + dim
    end do

    ! tensor fields
    do i = 1, counts(3)
      phase = phases%t(i)
      tfield => extract_tensor_field(state(phase), names%t(i))

      do j = 1, npart
        vals(field_idx:field_idx+dim**2-1,j) = reshape( &
             eval_field(ele(j), tfield, lcoords(:,j)), &
             [ dim**2 ])
      end do

      field_idx = field_idx + dim**2
    end do
  end subroutine evaluate_particle_fields

  !> Given the particle position and time, evaluate the specified python function for a group of particles
  !! specified in the string func at that location.
  subroutine set_particle_scalar_attribute_from_python(attributes, positions, natt, func, time, dt, is_array)
    !! (natt x nparts) array of calculated attribute values
    real, dimension(:,:), intent(out) :: attributes
    !! Current particle positions
    real, dimension(:,:), target, intent(in) :: positions
    !! Array dimension of this attribute
    integer, intent(in) :: natt
    !! Func may contain any python at all but the following function must
    !! be defined::
    !!  def val(X,t,dt)
    !! where X is position and t is the time. The result must be a float.
    character(len=*), intent(in) :: func
    !! Current simulation time
    real, intent(in) :: time
    !! Current simulation timestep
    real, intent(in) :: dt
    !! Whether this is an array-valued attribute
    logical, intent(in) :: is_array

    real, dimension(:), pointer :: lvx, lvy, lvz
    real, dimension(0), target :: zero
    integer :: stat, dim

    dim = size(positions, 1)
    select case(dim)
    case(1)
      lvx=>positions(1,:)
      lvy=>zero
      lvz=>zero
    case(2)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>zero
    case(3)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>positions(3,:)
    end select

    if (is_array) then
      ! call interface with additional array dimension argument
      call set_scalar_particles_from_python(func, len(func), dim, &
           size(attributes,2), natt, lvx, lvy, lvz, time, dt, attributes, stat)
    else
      call set_scalar_particles_from_python(func, len(func), dim, &
           size(attributes,2), lvx, lvy, lvz, time, dt, attributes, stat)
    end if
    if (stat/=0) then
      ! FIX ME - Does not print
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1 , *) trim(func)
      FLExit("Dying")
    end if
  end subroutine set_particle_scalar_attribute_from_python

  !> Given a particle position, time and field values, evaluate the python function
  !! specified in the string func at that location.
  subroutine set_particle_scalar_attribute_from_python_fields(particle_list, state, positions, lcoords, ele, natt, &
       attributes, old_attr_names, old_attr_counts, old_attr_Dims, old_attributes, field_names, field_counts, old_field_names, &
       old_field_counts, func, time, dt, is_array, first_newly_init_part)
    !! Particle list for which to evaluate the function
    type(detector_linked_list), intent(in) :: particle_list
    !! Model state structure
    type(state_type), dimension(:), intent(in) :: state
    !! Current particle positions
    real, dimension(:,:), target, intent(in) :: positions
    !! Local coordinates of particle positions
    real, dimension(:,:), intent(in) :: lcoords
    !! Elements containing particles
    integer, dimension(:), intent(in) :: ele
    !! Array dimension of this attribute
    integer, intent(in) :: natt
    !! Attribute values to set
    real, dimension(:,:), intent(out) :: attributes
    !! Names of attributes stored from the previous timestep
    character, dimension(:,:), intent(in) :: old_attr_names
    !! Number of each of scalar, vector, tensor old attributes
    integer, dimension(3), intent(in) :: old_attr_counts
    !! Array dimensions of old attributes
    integer, dimension(:), intent(in) :: old_attr_dims
    !! Attribute values from the previous timestep
    real, dimension(:,:), intent(in) :: old_attributes
    !! Names of fields that are to be passed to Python
    character, dimension(:,:), intent(in) :: field_names
    !! Number of each of scalar, vector, tensor fields
    integer, dimension(3), intent(in) :: field_counts
    !! Names of old fields that are to be passed to Python
    character, dimension(:,:), intent(in) :: old_field_names
    !! Number of each of salar, vector, tensor old fields
    integer, dimension(3), intent(in) :: old_field_counts
    !! Func may contain any python at all but the following function must
    !! be defined::
    !!  def val(X,t,dt,fields)
    !! where X is position t is the time, and fields is a dictionary where fields["FieldName"] and fields["OldFieldName" store
    !! the interpolated value of "FieldName" at the particle position at the current and previous time levels, and fields["OldAttributeName"] stores the attribute
    !! value at the previous time level. The result must be a float.
    character(len=*), intent(in) :: func
    !! Current model time
    real, intent(in) :: time
    !! Current model timestep
    real, intent(in) :: dt
    !! Whether this is an array-valued attribute
    logical, intent(in) :: is_array
    !> Pointer to the first newly initialised particle
    type(detector_type), pointer, optional :: first_newly_init_part

    ! locals
    integer :: i, j, field_idx
    integer :: dim, stat, phase
    integer :: nparts
    real, dimension(:), pointer :: lvx, lvy, lvz
    real, dimension(0), target :: zero
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    real, allocatable, dimension(:,:) :: field_vals, old_field_vals
    type(detector_type), pointer :: particle

    nparts = size(attributes, 2)
    dim = size(positions, 1)
    select case(dim)
    case(1)
      lvx=>positions(1,:)
      lvy=>zero
      lvz=>zero
    case(2)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>zero
    case(3)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>positions(3,:)
    end select

    ! allocate space to hold fields and old fields
    allocate(field_vals(field_counts(1) + dim*field_counts(2) + dim**2*field_counts(3), nparts))
    allocate(old_field_vals(old_field_counts(1) + dim*old_field_counts(2) + dim**2*old_field_counts(3), nparts))

    call evaluate_particle_fields(nparts, state, ele, lcoords, &
         particle_list%field_names, particle_list%field_phases, field_counts, &
         field_vals, dim)

    ! copy old fields off particles
    if (present(first_newly_init_part)) then
      particle => first_newly_init_part
    else
      particle => particle_list%first
    end if

    do i = 1, nparts
      old_field_vals(:,i) = particle%old_fields(:)
      particle => particle%next
    end do

    call set_scalar_particles_from_python_fields(func, len(func), dim, nparts, natt, &
         lvx, lvy, lvz, time, dt, field_counts, field_names, field_vals, old_field_counts, old_field_names, &
         old_field_vals, old_attr_counts, old_attr_names, old_attr_dims, old_attributes, is_array, attributes, stat)
    if (stat/=0) then
      ! FIX ME - Does not print
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1 , *) trim(func)
      FLExit("Dying")
    end if

    deallocate(field_vals)
    deallocate(old_field_vals)
  end subroutine set_particle_scalar_attribute_from_python_fields

  !> Given the particle position and time, evaluate the specified python function for a group of particles
  !! specified in the string func at that location.
  subroutine set_particle_vector_attribute_from_python(attributes, positions, natt, func, time, dt, is_array)
    !! Attribute values to set
    real, dimension(:,:), intent(out) :: attributes
    !! Current particle positions
    real, dimension(:,:), target, intent(in) :: positions
    !! Array dimension of this attribute
    integer, intent(in) :: natt
    !! Func may contain any python at all but the following function must
    !! be defined::
    !!  def val(X,t,dt)
    !! where X is position and t is the time. The result must be a float.
    character(len=*), intent(in) :: func
    !! Current simulation time
    real, intent(in) :: time
    !! Curretn simulation timestep
    real, intent(in) :: dt
    !! Whether this is an array-valued attribute
    logical, intent(in) :: is_array

    real, dimension(:), pointer :: lvx,lvy,lvz
    real, dimension(0), target :: zero
    integer :: stat, dim

    dim = size(positions, 1)
    select case(dim)
    case(1)
      lvx=>positions(1,:)
      lvy=>zero
      lvz=>zero
    case(2)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>zero
    case(3)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>positions(3,:)
    end select

    if (is_array) then
      ! call interface with additional array dimension argumetn
      call set_vector_particles_from_python(func, len(func), dim, &
            size(attributes,2), natt, lvx, lvy, lvz, time, dt, attributes, stat)
    else
      call set_vector_particles_from_python(func, len(func), dim, &
            size(attributes,2), lvx, lvy, lvz, time, dt, attributes, stat)
    end if
    if (stat/=0) then
      ! FIX ME - Does not print
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1 , *) trim(func)
      FLExit("Dying")
    end if
  end subroutine set_particle_vector_attribute_from_python

  !> Given a particle position, time and field values, evaluate the python function
  !! specified in the string func at that location.
  subroutine set_particle_vector_attribute_from_python_fields(particle_list, state, positions, lcoords, ele, natt, &
       attributes, old_attr_names, old_attr_counts, old_attr_dims, old_attributes, field_names, field_counts, old_field_names, &
       old_field_counts, func, time, dt, is_array, first_newly_init_part)
    !! Particle list for which to evaluate the function
    type(detector_linked_list), intent(in) :: particle_list
    !! Model state structure
    type(state_type), dimension(:), intent(in) :: state
    !! Current particle positions
    real, dimension(:,:), target, intent(in) :: positions
    !! Local coordinates of particle positions
    real, dimension(:,:), intent(in) :: lcoords
    !! Elements containing particles
    integer, dimension(:), intent(in) :: ele
    !! Array dimnesion of this attribute
    integer, intent(in) :: natt
    !! Attribute values to set
    real, dimension(:,:), intent(out) :: attributes
    !! Names of attributes stored from the previous timestep
    character, dimension(:,:), intent(in) :: old_attr_names
    !! Number of each of scalar, vector, tensor old attributes
    integer, dimension(3), intent(in) :: old_attr_counts
    !! Array dimensions of old attributes
    integer, dimension(:), intent(in) :: old_attr_dims
    !! Attribute values from the previous timestep
    real, dimension(:,:), intent(in) :: old_attributes
    !! Names of fields that are to be passed to Python
    character, dimension(:,:), intent(in) :: field_names
    !! Number of each of scalar, vector, tensor fields
    integer, dimension(3), intent(in) :: field_counts
    !! Names of old fields that are to be passed to Python
    character, dimension(:,:), intent(in) :: old_field_names
    !! Number of each of salar, vector, tensor old fields
    integer, dimension(3), intent(in) :: old_field_counts
    !! Func may contain any python at all but the following function must
    !! be defined::
    !!  def val(X,t,dt,fields)
    !! where X is position t is the time, and fields is a dictionary where fields["FieldName"] and fields["OldFieldName" store
    !! the interpolated value of "FieldName" at the particle position at the current and previous time levels, and fields["OldAttributeName"] stores the attribute
    !! value at the previous time level. The result must be a float.
    character(len=*), intent(in) :: func
    !! Current model time
    real, intent(in) :: time
    !! Current model timestep
    real, intent(in) :: dt
    !! Whether this is an array-valued attribute
    logical, intent(in) :: is_array
    !> Pointer to the first newly initialised particle
    type(detector_type), pointer, optional :: first_newly_init_part

    ! locals
    integer :: i, j, field_idx
    integer :: dim, stat, phase
    integer :: nparts
    real, dimension(:), pointer :: lvx, lvy, lvz
    real, dimension(0), target :: zero
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    real, allocatable, dimension(:,:) :: field_vals, old_field_vals
    type(detector_type), pointer :: particle

    nparts = size(attributes, 2)
    dim = size(positions, 1)
    select case(dim)
    case(1)
      lvx=>positions(1,:)
      lvy=>zero
      lvz=>zero
    case(2)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>zero
    case(3)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>positions(3,:)
    end select

    ! allocate space to hold fields and old fields
    allocate(field_vals(field_counts(1) + dim*field_counts(2) + dim**2*field_counts(3), nparts))
    allocate(old_field_vals(old_field_counts(1) + dim*old_field_counts(2) + dim**2*old_field_counts(3), nparts))

    call evaluate_particle_fields(nparts, state, ele, lcoords, &
         particle_list%field_names, particle_list%field_phases, field_counts, &
         field_vals, dim)

    ! copy old fields off particles
    if (present(first_newly_init_part)) then
      particle => first_newly_init_part
    else
      particle => particle_list%first
    end if

    do i = 1, nparts
      old_field_vals(:,i) = particle%old_fields(:)
      particle => particle%next
    end do

    call set_vector_particles_from_python_fields(func, len(func), dim, nparts, natt, &
         lvx, lvy, lvz, time, dt, field_counts, field_names, field_vals, old_field_counts, old_field_names, &
         old_field_vals, old_attr_counts, old_attr_names, old_attr_dims, old_attributes, is_array, attributes, stat)
    if (stat/=0) then
      ! FIX ME - Does not print
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1 , *) trim(func)
      FLExit("Dying")
    end if

    deallocate(field_vals)
    deallocate(old_field_vals)
  end subroutine set_particle_vector_attribute_from_python_fields

  !> Given the particle position and time, evaluate the specified python function for a group of particles
  !! specified in the string func at that location.
  subroutine set_particle_tensor_attribute_from_python(attributes, positions, natt, func, time, dt, is_array)
    !! Attribute values to set
    real, dimension(:,:), intent(inout) :: attributes
    !! Current particle positions
    real, dimension(:,:), target, intent(in) :: positions
    !! Array dimension of this attribute
    integer, intent(in) :: natt
    !! Func may contain any python at all but the following function must
    !! be defined::
    !!  def val(X,t,dt)
    !! where X is position and t is the time. The result must be a float.
    character(len=*), intent(in) :: func
    !! Current simulation time
    real, intent(in) :: time
    !! Curretn simulation timestep
    real, intent(in) :: dt
    !! Whether this is an array-valued attribute
    logical, intent(in) :: is_array

    real, dimension(:), pointer :: lvx,lvy,lvz
    real, dimension(0), target :: zero
    real, dimension(:,:,:,:), allocatable :: tensor_res
    integer :: i, j, k, dim, idx
    integer :: nparts
    integer :: stat

    nparts = size(attributes, 2)
    dim = size(positions, 1)
    select case(dim)
    case(1)
      lvx=>positions(1,:)
      lvy=>zero
      lvz=>zero
    case(2)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>zero
    case(3)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>positions(3,:)
    end select

    allocate(tensor_res(dim,dim,natt,nparts))
    if (is_array) then
      call set_tensor_particles_from_python(func, len(func), dim, &
            nparts, natt, lvx, lvy, lvz, time, dt, tensor_res, stat)
    else
      call set_tensor_particles_from_python(func, len(func), dim, &
            nparts, lvx, lvy, lvz, time, dt, tensor_res, stat)
    end if
    if (stat/=0) then
      ! FIX ME - Does not print
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1 , *) trim(func)
      FLExit("Dying")
    end if

    ! convert tensor to array of attributes
    attributes(:,:) = reshape(tensor_res, [natt * dim**2, nparts])
    deallocate(tensor_res)
  end subroutine set_particle_tensor_attribute_from_python
  !> Given a particle position, time and field values, evaluate the python function
  !! specified in the string func at that location.
  subroutine set_particle_tensor_attribute_from_python_fields(particle_list, state, positions, lcoords, ele, natt, &
       attributes, old_attr_names, old_attr_counts, old_attr_dims, old_attributes, field_names, field_counts, old_field_names, &
       old_field_counts, func, time, dt, is_array, first_newly_init_part)
    !! Particle list for which to evaluate the function
    type(detector_linked_list), intent(in) :: particle_list
    !! Model state structure
    type(state_type), dimension(:), intent(in) :: state
    !! Current particle positions
    real, dimension(:,:), target, intent(in) :: positions
    !! Local coordinates of particle positions
    real, dimension(:,:), intent(in) :: lcoords
    !! Elements containing particles
    integer, dimension(:), intent(in) :: ele
    !! Array dimension of this attribute
    integer, intent(in) :: natt
    !! Attribute values to set
    real, dimension(:,:), intent(out) :: attributes
    !! Names of attributes stored from the previous timestep
    character, dimension(:,:), intent(in) :: old_attr_names
    !! Number of each of scalar, vector, tensor old attributes
    integer, dimension(3), intent(in) :: old_attr_counts
    !! Array dimensions of old attributes
    integer, dimension(:), intent(in) :: old_attr_dims
    !! Attribute values from the previous timestep
    real, dimension(:,:), intent(in) :: old_attributes
    !! Names of fields that are to be passed to Python
    character, dimension(:,:), intent(in) :: field_names
    !! Number of each of scalar, vector, tensor fields
    integer, dimension(3), intent(in) :: field_counts
    !! Names of old fields that are to be passed to Python
    character, dimension(:,:), intent(in) :: old_field_names
    !! Number of each of salar, vector, tensor old fields
    integer, dimension(3), intent(in) :: old_field_counts
    !! Func may contain any python at all but the following function must
    !! be defined::
    !!  def val(X,t,dt,fields)
    !! where X is position t is the time, and fields is a dictionary where fields["FieldName"] and fields["OldFieldName" store
    !! the interpolated value of "FieldName" at the particle position at the current and previous time levels, and fields["OldAttributeName"] stores the attribute
    !! value at the previous time level. The result must be a float.
    character(len=*), intent(in) :: func
    !! Current model time
    real, intent(in) :: time
    !! Current model timestep
    real, intent(in) :: dt
    !! Whether this is an array-valued attribute
    logical, intent(in) :: is_array
    !> Pointer to the first newly initialised particle
    type(detector_type), pointer, optional :: first_newly_init_part

    ! locals
    integer :: i, j, field_idx
    integer :: dim, stat, phase
    integer :: nparts
    real, dimension(:), pointer :: lvx, lvy, lvz
    real, dimension(0), target :: zero
    type(scalar_field), pointer :: sfield
    type(vector_field), pointer :: vfield
    type(tensor_field), pointer :: tfield
    real, allocatable, dimension(:,:) :: field_vals, old_field_vals
    type(detector_type), pointer :: particle
    real, dimension(:,:,:,:), allocatable :: tensor_res

    nparts = size(attributes, 2)
    dim = size(positions, 1)
    select case(dim)
    case(1)
      lvx=>positions(1,:)
      lvy=>zero
      lvz=>zero
    case(2)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>zero
    case(3)
      lvx=>positions(1,:)
      lvy=>positions(2,:)
      lvz=>positions(3,:)
    end select

    ! allocate space to hold fields and old fields
    allocate(field_vals(field_counts(1) + dim*field_counts(2) + dim**2*field_counts(3), nparts))
    allocate(old_field_vals(old_field_counts(1) + dim*old_field_counts(2) + dim**2*old_field_counts(3), nparts))

    call evaluate_particle_fields(nparts, state, ele, lcoords, &
         particle_list%field_names, particle_list%field_phases, field_counts, &
         field_vals, dim)

    ! copy old fields off particles
    if (present(first_newly_init_part)) then
      particle => first_newly_init_part
    else
      particle => particle_list%first
    end if

    do i = 1, nparts
      old_field_vals(:,i) = particle%old_fields(:)
      particle => particle%next
    end do

    allocate(tensor_res(dim,dim,natt,nparts))

    call set_tensor_particles_from_python_fields(func, len(func), dim, nparts, natt, &
         lvx, lvy, lvz, time, dt, field_counts, field_names, field_vals, old_field_counts, old_field_names, &
         old_field_vals, old_attr_counts, old_attr_names, old_attr_dims, old_attributes, is_array, tensor_res, stat)
    if (stat/=0) then
      ! FIX ME - Does not print
       ewrite(-1, *) "Python error, Python string was:"
       ewrite(-1 , *) trim(func)
       FLExit("Dying")
    end if

    ! convert tensor to array of attributes
    attributes(:,:) = reshape(tensor_res, [natt * dim**2, nparts])

    deallocate(tensor_res)
    deallocate(field_vals)
    deallocate(old_field_vals)
  end subroutine set_particle_tensor_attribute_from_python_fields
end module detector_tools
