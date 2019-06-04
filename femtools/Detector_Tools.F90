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
  use embed_python, only: set_detectors_from_python, set_particles_from_python_fields, set_particles_from_python
  use integer_hash_table_module
  use parallel_tools
  use parallel_fields, only: element_owned
  use fields
  use state_module, only: state_type, extract_scalar_field, aliased
  use field_options
  use detector_data_types
  use pickers
  
  implicit none
  
  private

  public :: insert, allocate, deallocate, copy, move, move_all, remove, &
            delete, delete_all, pack_detector, unpack_detector, &
            detector_value, set_detector_coords_from_python, &
            detector_buffer_size, set_particle_attribute_from_python, &
            set_particle_attribute_from_python_fields, temp_insert, &
            temp_deallocate, temp_delete, temp_delete_all, temp_remove

  interface insert
     module procedure insert_into_detector_list
  end interface

  interface temp_insert
     module procedure insert_into_temp_detector_list
  end interface

  ! Removes detector from a list without deallocating it
  interface remove
     module procedure remove_detector_from_list
  end interface

  ! Removes detector from a temp list without deallocating it
  interface temp_remove
     module procedure remove_detector_from_temp_list
  end interface

  ! Removes detector from a list and deallocates it
  interface delete
     module procedure delete_detector
  end interface

  ! Removes detector from a temp list
  interface temp_delete
     module procedure delete_temp_detector
  end interface

  ! Deletes all detectors from a given list
  interface delete_all
     module procedure delete_all_detectors
  end interface

  ! Deletes all detectors from a given temp list
  interface temp_delete_all
     module procedure delete_all_temp_detectors
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

  interface temp_deallocate
     module procedure detector_list_temp_deallocate
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
    if (allocated(detector_list%detector_names)) then
       deallocate(detector_list%detector_names)
    end if
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

  subroutine detector_list_temp_deallocate(detector_list)
    type(detector_linked_list), pointer :: detector_list

    type(rk_gs_parameters), pointer :: parameters
    ! Delete detectors
    if (detector_list%length > 0) then
       call temp_delete_all(detector_list)
    end if
    ! Deallocate list information
    if (allocated(detector_list%detector_names)) then
       deallocate(detector_list%detector_names)
    end if
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

  end subroutine detector_list_temp_deallocate
    
  subroutine detector_copy(new_detector, old_detector)
    ! Copies all the information from the old detector to
    ! the new detector
    type(detector_type), pointer, intent(in) :: old_detector
    type(detector_type),  pointer :: new_detector
      
    new_detector%position = old_detector%position
    new_detector%element = old_detector%element
    new_detector%id_number = old_detector%id_number
    new_detector%type = old_detector%type
    new_detector%name = old_detector%name
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

  subroutine insert_into_temp_detector_list(detector, current_list)
    ! Inserts detector at the end of a list
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

  end subroutine insert_into_temp_detector_list

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

  subroutine remove_detector_from_temp_list(detector, detector_list)
    !! Removes the detector from the list,
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

       if (associated(detector%next)) then
          detector%temp_next%temp_previous => detector%temp_previous
       else
          detector_list%last => detector%temp_previous
       end if
    end if
    detector_list%length = detector_list%length-1

  end subroutine remove_detector_from_temp_list

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

  subroutine delete_temp_detector(detector, detector_list)
    ! Removes and deallocates the given detector
    ! and outputs the next detector in the list as detector
    type(detector_type), pointer :: detector
    type(detector_linked_list), intent(inout), optional :: detector_list
    
    type(detector_type), pointer :: temp_detector
    temp_detector => detector
    detector => detector%temp_next
    call temp_remove(temp_detector, detector_list)
    
  end subroutine delete_temp_detector

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

  subroutine delete_all_temp_detectors(detector_list)
    ! Remove and deallocate all detectors in a list
    type(detector_linked_list), intent(inout) :: detector_list
    type(detector_type), pointer :: detector

    detector => detector_list%first
    do while (associated(detector))
       call temp_delete(detector,detector_list)
    end do

  end subroutine delete_all_temp_detectors

  function detector_buffer_size(ndims, have_update_vector, nstages, attribute_size)
    ! Returns the number of reals we need to pack a detector
    integer, intent(in) :: ndims
    logical, intent(in) :: have_update_vector
    integer, intent(in), optional :: nstages
    integer :: detector_buffer_size
    integer, dimension(3), optional, intent(in) :: attribute_size !array to hold size of attributes
    if (present(attribute_size)) then
       if (have_update_vector) then
          assert(present(nstages))
          detector_buffer_size=(nstages+2)*ndims+3+sum(attribute_size)
       else
          detector_buffer_size=ndims+4+sum(attribute_size)
       end if
    else
       if (have_update_vector) then
          assert(present(nstages))
          detector_buffer_size=(nstages+2)*ndims+3
       else
          detector_buffer_size=ndims+4
       end if 
    end if

  end function detector_buffer_size

  subroutine pack_detector(detector,buff,ndims,nstages, attribute_size)
    ! Packs (serialises) detector into buff
    ! Basic fields are: element, position, id_number and type
    ! If nstages is given, the detector is still moving
    ! and we also pack update_vector and k
    type(detector_type), pointer, intent(in) :: detector
    real, dimension(:), intent(out) :: buff
    integer, intent(in) :: ndims
    integer, intent(in), optional :: nstages
    integer, dimension(3), optional, intent(in) :: attribute_size !array to hold size of attributes
    integer :: det_params
    assert(size(detector%position)==ndims)
    !Set size of basic detector fields, being detector element, id_number and type
    det_params = 3
    !Check if this is a particle carrying attributes
    if (present(attribute_size)) then
       assert(size(buff)>=ndims+det_params+sum(attribute_size))
       ! Basic fields: ndims+det_params
       buff(1:ndims) = detector%position
       buff(ndims+1) = detector%element
       buff(ndims+2) = detector%id_number
       buff(ndims+3) = detector%type
       if (attribute_size(1).ne.0) then
          buff(ndims+det_params+1:ndims+det_params+attribute_size(1)) = detector%attributes
       end if
       if (attribute_size(2).ne.0) then
          buff(ndims+det_params+1+attribute_size(1):ndims+det_params+attribute_size(1)+attribute_size(2)) &
               = detector%old_attributes
       end if
       if (attribute_size(3).ne.0) then
          buff(ndims+det_params+1+attribute_size(1)+attribute_size(2):ndims+det_params+sum(attribute_size)) &
               = detector%old_fields
       end if
       ! Lagrangian advection fields: (nstages+1)*ndims
       if (present(nstages)) then
          assert(size(buff)==(nstages+2)*ndims+det_params+sum(attribute_size))
          assert(allocated(detector%update_vector))
          assert(allocated(detector%k))
          
          buff(ndims+det_params+1+sum(attribute_size):2*ndims+det_params &
               +sum(attribute_size)) = detector%update_vector
          buff(2*ndims+det_params+1+sum(attribute_size):(nstages+2)*ndims &
               +det_params+sum(attribute_size)) = reshape(detector%k,(/nstages*ndims/))
       else
          assert(size(buff)==ndims+det_params+1+sum(attribute_size))
          buff(ndims+det_params+1+sum(attribute_size)) = detector%list_id
       end if
    else
       assert(size(buff)>=ndims+det_params)

       ! Basic fields: ndims+det_params
       buff(1:ndims) = detector%position
       buff(ndims+1) = detector%element
       buff(ndims+2) = detector%id_number
       buff(ndims+3) = detector%type

       ! Lagrangian advection fields: (nstages+1)*ndims
       if (present(nstages)) then
          assert(size(buff)==(nstages+2)*ndims+det_params)
          assert(allocated(detector%update_vector))
          assert(allocated(detector%k))
          
          buff(ndims+det_params+1:2*ndims+det_params) = detector%update_vector
          buff(2*ndims+det_params+1:(nstages+2)*ndims+det_params) = reshape(detector%k,(/nstages*ndims/))
       else
          assert(size(buff)==ndims+det_params+1)
          buff(ndims+det_params+1) = detector%list_id
       end if
    end if
    
  end subroutine pack_detector

  subroutine unpack_detector(detector,buff,ndims,global_to_local,coordinates,nstages,attribute_size)
    ! Unpacks the detector from buff and fills in the blanks
    type(detector_type), pointer :: detector
    real, dimension(:), intent(in) :: buff
    integer, intent(in) :: ndims
    type(integer_hash_table), intent(in), optional :: global_to_local
    type(vector_field), intent(in), optional :: coordinates
    integer, intent(in), optional :: nstages
    integer, dimension(3), optional, intent(in) :: attribute_size !array to hold size of attributes

    integer :: det_params
    !Set size of basic detector fields, being detector element, id_number and type
    det_params = 3
    
    !Check if this is a particle carrying attributes
    if (present(attribute_size)) then

       if (.not. allocated(detector%position)) then
          allocate(detector%position(ndims))
       end if
       if (.not. allocated(detector%attributes)) then
          allocate(detector%attributes(attribute_size(1)))
          allocate(detector%old_attributes(attribute_size(2)))
          allocate(detector%old_fields(attribute_size(3)))
       end if
       ! Basic fields: ndims+3
       detector%position = buff(1:ndims)
       detector%element = buff(ndims+1)
       detector%id_number = buff(ndims+2)
       detector%type = buff(ndims+3)
       if (attribute_size(1)/=0) then  
          detector%attributes = buff(ndims+det_params+1:ndims+det_params+attribute_size(1))
       end if
       if (attribute_size(2)/=0) then  
          detector%old_attributes = buff(ndims+det_params+1+attribute_size(1):ndims+det_params+attribute_size(1) &
               +attribute_size(2))
       end if
       if (attribute_size(3)/=0) then  
          detector%old_fields = buff(ndims+det_params+1+attribute_size(1)+attribute_size(2):ndims+det_params+ &
               sum(attribute_size))
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
          assert(size(buff)==(nstages+2)*ndims+det_params+sum(attribute_size))
          
          ! update_vector, dimension(ndim)
          if (.not. allocated(detector%update_vector)) then
             allocate(detector%update_vector(ndims))
          end if
          detector%update_vector = buff(ndims+det_params+1+sum(attribute_size) &
               :2*ndims+det_params+sum(attribute_size))
          
          ! k, dimension(nstages:ndim)
          if (.not. allocated(detector%k)) then
             allocate(detector%k(nstages,ndims))
          end if
          detector%k = reshape(buff(2*ndims+det_params+1+sum(attribute_size):(nstages+2) &
               *ndims+det_params+sum(attribute_size)),(/nstages,ndims/))
          
          ! If update_vector still exists, we're not done moving
          detector%search_complete=.false.
       else
          assert(size(buff)==ndims+det_params+1+sum(attribute_size))
          
          detector%list_id = buff(ndims+det_params+1+sum(attribute_size))
          detector%search_complete=.true.
       end if
    else

       if (.not. allocated(detector%position)) then
          allocate(detector%position(ndims))
       end if
       
       ! Basic fields: ndims+3
       detector%position = buff(1:ndims)
       detector%element = buff(ndims+1)
       detector%id_number = buff(ndims+2)
       detector%type = buff(ndims+3)
       
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
          assert(size(buff)==(nstages+2)*ndims+det_params)
          
          ! update_vector, dimension(ndim)
          if (.not. allocated(detector%update_vector)) then
             allocate(detector%update_vector(ndims))
          end if
          detector%update_vector = buff(ndims+det_params+1:2*ndims+det_params)
          
          ! k, dimension(nstages:ndim)
          if (.not. allocated(detector%k)) then
             allocate(detector%k(nstages,ndims))
          end if
          detector%k = reshape(buff(2*ndims+det_params+1:(nstages+2)*ndims+det_params),(/nstages,ndims/))
          
          ! If update_vector still exists, we're not done moving
          detector%search_complete=.false.
       else
          assert(size(buff)==ndims+det_params+1)
          
          detector%list_id = buff(ndims+det_params+1)
          detector%search_complete=.true.
       end if
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
    !! be defiled:
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
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1 , *) trim(func)
      FLExit("Dying")
    end if

  end subroutine set_detector_coords_from_python

  subroutine set_particle_attribute_from_python(attributes, positions, nparts, func, time, dt)
    !!< Given the particle position and time, evaluate the specified python function for a group of particles
    !!< specified in the string func at that location. 
    real, dimension(:), intent(inout) :: attributes
    !! Func may contain any python at all but the following function must
    !! be defined::
    !!  def val(X,t)
    !! where X is position and t is the time. The result must be a float. 
    character(len=*), intent(in) :: func !function for attributes to be set from
    integer, intent(in) :: nparts !number of particles
    real, intent(in) :: time !current simulation time
    real, intent(in) :: dt !current timestep
    real, dimension(:,:), target, intent(in) :: positions !current particle positions
    real, dimension(:), pointer :: lvx,lvy,lvz
    real, dimension(0), target :: zero
    integer :: stat

    select case(size(positions(:,1)))
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
    call set_particles_from_python(func, len(func), size(positions(:,1)), &
         nparts, lvx, lvy, lvz, time, dt, attributes, stat)
    if (stat/=0) then
      ewrite(-1, *) "Python error, Python string was:"
      ewrite(-1 , *) trim(func)
      FLExit("Dying")
    end if
  end subroutine set_particle_attribute_from_python

  subroutine set_particle_attribute_from_python_fields(particle_list, state, positions, lcoords, ele, nparts, attributes, old_att_names, old_attributes, func, time, dt)
    !!< Given a particle position, time and field values, evaluate the python function
    !!< specified in the string func at that location. 
    !! Func may contain any python at all but the following function must
    !! be defined::
    !!  def val(X,t, fields)
    !! where X is position t is the time, and fields is a dictionary where fields["FieldName"] and fields["OldFieldName" store
    !! the interpolated value of "FieldName" at the particle position at the current and previous time levels, and fields["OldAttributeName"] stores the attribute
    !! value at the previous time level. The result must be a float. 
    type(detector_linked_list), intent(in) :: particle_list
    type(state_type), dimension(:), intent(in) :: state
    real, dimension(:,:), intent(in) :: old_attributes !Attribute values from the previous timestep
    character, dimension(:,:), intent(in) :: old_att_names !Names of attributes stored from the previous timestep
    real, dimension(:), intent(inout) :: attributes !Attribute values to be set
    real, dimension(:,:), target, intent(in) :: positions !positions of particles
    real, dimension(:,:), intent(in) :: lcoords !local coordinates of the particle position
    integer, dimension(:), intent(in) :: ele !elements particles are found in
    character(len=*), intent(in) :: func !function for attributes to be set from
    integer, intent(in) :: nparts !number of particles
    real, intent(in) :: time
    real, intent(in) :: dt
    character, allocatable, dimension(:,:) :: field_names
    character, allocatable, dimension(:,:) :: old_field_names
    real, allocatable, dimension(:,:) :: field_vals
    real, allocatable, dimension(:,:) :: old_field_vals
    real :: value
    real, dimension(:), pointer :: lvx,lvy,lvz
    real, dimension(0), target :: zero
    character(len=FIELD_NAME_LEN) :: buffer
    type(detector_type), pointer :: particle

    type(scalar_field), pointer :: sfield
    character(len=FIELD_NAME_LEN) :: name
    integer :: phase, i, j, nfields
    integer :: p, f, stat, k, pfields, l, old_nfields, old_nattributes

    nfields=0
    do phase=1,size(state)
       do i=1, size(state(phase)%scalar_names)
          sfield => extract_scalar_field(state(phase),state(phase)%scalar_names(i))
          if (sfield%option_path=="".or.aliased(sfield)) then
             cycle
          else if (have_option(trim(complete_field_path(sfield%option_path)) // "/particles/include_in_particles")) then
             nfields=nfields+1
          end if
       end do
    end do

    select case(size(positions(:,1)))
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

    particle => particle_list%first
    
    old_nfields = size(particle%old_fields)
    old_nattributes = size(particle%old_attributes)

    allocate(field_names(FIELD_NAME_LEN,nfields))
    allocate(field_vals(nfields,nparts))
    allocate(old_field_names(FIELD_NAME_LEN,size(particle%old_fields)))
    allocate(old_field_vals(size(particle%old_fields),nparts))
    j=1
    l=1
    do phase=1,size(state)
       do f = 1, size(state(phase)%scalar_names)
          sfield => extract_scalar_field(state(phase),state(phase)%scalar_names(f))
          if (sfield%option_path=="".or.aliased(sfield)) then
             cycle
          else if (have_option(trim(complete_field_path(sfield%option_path)) // "/particles/include_in_particles")) then
             name=trim(state(phase)%scalar_names(f))
             do k = 1,len_trim(name)
                field_names(k,j)=name(k:k)
             end do
             field_names(k,j)= C_NULL_CHAR
             if (have_option(trim(complete_field_path(sfield%option_path)) // "/particles/include_in_particles/store_old_field")) then
                old_field_names(1,l) = 'O'
                old_field_names(2,l) = 'l'
                old_field_names(3,l) = 'd'
                do k = 4,len_trim(name)+3
                   old_field_names(k,l)=name(k-3:k-3)
                end do
                old_field_names(k,l)= C_NULL_CHAR
                do i = 1,nparts    
                   value = eval_field(ele(i), sfield, lcoords(:,i))
                   field_vals(j,i)=value
                   old_field_vals(l,i)=value
                end do
                l=l+1
             else
                do i = 1,nparts
                   value = eval_field(ele(i), sfield, lcoords(:,i))
                   field_vals(j,i)=value
                end do
             end if
             j=j+1
          end if
       end do
    end do
    particle => particle_list%first
    if (allocated(particle%old_fields)) then
       particle => particle_list%first
       do i = 1,nparts
          old_field_vals(:,i) = particle%old_fields
          particle=>particle%next
       end do
    end if
    
    call set_particles_from_python_fields(func, len(func), size(positions(:,1)), nparts, &
         lvx, lvy, lvz, time, dt, nfields, field_names, field_vals, old_nfields, old_field_names, &
         old_field_vals, old_nattributes, old_att_names, old_attributes, attributes, stat)
    if (stat/=0) then
       ewrite(-1, *) "Python error, Python string was:"
       ewrite(-1 , *) trim(func)
       FLExit("Dying")
    end if
    
    deallocate(field_names)
    deallocate(field_vals)
    deallocate(old_field_names)
    deallocate(old_field_vals)
    
  end subroutine set_particle_attribute_from_python_fields

end module detector_tools
