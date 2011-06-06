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
  use detector_data_types
  use fields
  use integer_hash_table_module
  
  implicit none
  
  private

  public :: insert, allocate, deallocate, copy, move, move_all, remove, &
            delete, delete_all, pack_detector, unpack_detector, &
            detector_value, set_detector_coords_from_python

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
     module procedure detector_deallocate
  end interface

  ! Evaluate field at the location of the detector.
  interface detector_value
     module procedure detector_value_scalar, detector_value_vector
  end interface

contains 

  subroutine detector_allocate_from_params(new_detector, ndims, local_coord_count)
    type(detector_type),  pointer, intent(out) :: new_detector
    integer, intent(in) :: ndims, local_coord_count
      
    assert(.not. associated(new_detector))
      
    ! allocate the memory for the new detector
    if (.not. associated(new_detector)) then
       allocate(new_detector)
    end if
    allocate(new_detector%position(ndims))
    allocate(new_detector%local_coords(local_coord_count))
      
    assert(associated(new_detector))
      
  end subroutine detector_allocate_from_params
    
  subroutine detector_allocate_from_detector(new_detector, old_detector)
    type(detector_type), pointer, intent(in) :: old_detector
    type(detector_type),  pointer, intent(out) :: new_detector
      
    integer :: ndims, local_coord_count
      
    ndims = size(old_detector%position)
    local_coord_count = size(old_detector%local_coords)
      
    ! allocate the memory for the new detector
    call detector_allocate_from_params(new_detector, ndims, local_coord_count)
      
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
       deallocate(detector)
    end if
    detector => null()

  end subroutine detector_deallocate
    
  subroutine detector_copy(new_detector, old_detector)
    ! Copies all the information from the old detector to
    ! the new detector
    type(detector_type), pointer, intent(in) :: old_detector
    type(detector_type),  pointer, intent(out) :: new_detector
      
    new_detector%position = old_detector%position
    new_detector%element = old_detector%element
    new_detector%id_number = old_detector%id_number
    new_detector%type = old_detector%type
    new_detector%name = old_detector%name
    new_detector%local_coords=old_detector%local_coords
      
  end subroutine detector_copy

  subroutine insert_into_detector_list(current_list,node)
    ! Inserts detector at the end of a list
    type(detector_linked_list), intent(inout) :: current_list
    type(detector_type), pointer :: node

    if (current_list%length == 0) then
       current_list%firstnode => node 
       current_list%lastnode => node 
       current_list%firstnode%previous => null()
       current_list%lastnode%next => null()
       current_list%length = 1
    else
       node%previous => current_list%lastnode
       current_list%lastnode%next => node
       current_list%lastnode => node
       current_list%lastnode%next => null()
       current_list%length = current_list%length+1
    end if

  end subroutine insert_into_detector_list

  subroutine remove_detector_from_list(detector_list,node)
    !! Removes the detector from the list, 
    !! but does not deallocated it
    type(detector_linked_list), intent(inout) :: detector_list
    type(detector_type), pointer, intent(inout) :: node

    if (detector_list%length==1) then
       detector_list%firstnode => null()
       detector_list%lastnode => null()
    else
       if (associated(node%previous)) then
          node%previous%next => node%next
       else
          detector_list%firstnode => node%next
       end if

       if (associated(node%next)) then
          node%next%previous => node%previous
       else
          detector_list%lastnode => node%previous
       end if
    end if

    detector_list%length = detector_list%length-1

  end subroutine remove_detector_from_list

  subroutine delete_detector(detector_list,node)
    ! Removes and deallocates the given detector 
    ! and outputs the next detector in the list as detector
    type(detector_linked_list), intent(inout) :: detector_list
    type(detector_type), pointer, intent(inout) :: node
    
    type(detector_type), pointer :: temp_node
      
    temp_node => node
    node => node%next
    call remove(detector_list,temp_node)
    call deallocate(temp_node)    
      
  end subroutine delete_detector

  subroutine move_detector(from_list,node,to_list)
    ! Move detector from one list to the other
    type(detector_linked_list), intent(inout) :: from_list
    type(detector_type), pointer, intent(inout) :: node
    type(detector_linked_list), intent(inout) :: to_list

    call remove(from_list,node)
    call insert(to_list,node)  

  end subroutine move_detector

  subroutine move_all_detectors(from_list,to_list)
    ! Move all detectors from one list to the other
    type(detector_linked_list), intent(inout) :: from_list
    type(detector_linked_list), intent(inout) :: to_list
    type(detector_type), pointer :: node

    do while (associated(from_list%firstnode))
       node => from_list%firstnode
       call move(from_list,node,to_list)   
    end do

  end subroutine move_all_detectors

  subroutine delete_all_detectors(detector_list)
    ! Remove and deallocate all detectors in a list   
    type(detector_linked_list), intent(inout) :: detector_list
    type(detector_type), pointer :: node

    node => detector_list%firstnode
    do while (associated(node))
       call delete(detector_list,node) 
    end do

  end subroutine delete_all_detectors

  subroutine pack_detector(detector,buff,ndims,nstages)
    ! Packs (serialises) detector into buff
    ! Basic fields are: element, position, id_number and type
    ! If nstages is given, the detector is still moving
    ! and we also pack update_vector and k
    type(detector_type), pointer, intent(in) :: detector
    real, dimension(:), intent(out) :: buff
    integer, intent(in) :: ndims
    integer, intent(in), optional :: nstages

    assert(detector%element>0)
    assert(size(buff)>=ndims+3)

    ! Basic fields: ndims+3
    buff(1) = detector%element
    buff(2:ndims+1) = detector%position
    buff(ndims+2) = detector%id_number
    buff(ndims+3) = detector%type

    ! Lagrangian advection fields: (nstages+1)*ndims
    if (present(nstages)) then
       assert(size(buff)==(nstages+2)*ndims+3)
       assert(allocated(detector%update_vector))
       assert(allocated(detector%k))

       buff(ndims+4:2*ndims+3) = detector%update_vector
       buff(2*ndims+4:(nstages+2)*ndims+3) = reshape(detector%k,(/nstages*ndims/))
    else
       assert(size(buff)==ndims+3)
    end if
    
  end subroutine pack_detector

  subroutine unpack_detector(detector,buff,ndims,global_to_local,coordinates,nstages)
    ! Unpacks the detector from buff and fills in the blanks
    type(detector_type), pointer, intent(inout) :: detector
    real, dimension(:), intent(in) :: buff
    integer, intent(in) :: ndims
    type(integer_hash_table), intent(in), optional :: global_to_local
    type(vector_field), intent(in), optional :: coordinates
    integer, intent(in), optional :: nstages    

    assert(size(buff)>=ndims+3)

    if (.not. allocated(detector%position)) then
       allocate(detector%position(ndims))
    end if

    ! Basic fields: ndims+3
    detector%element = buff(1)
    detector%position = reshape(buff(2:ndims+1),(/ndims/))
    detector%id_number = buff(ndims+2)
    detector%type = buff(ndims+3)

    assert(detector%element>0)

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
       assert(size(buff)==(nstages+2)*ndims+3)

       ! update_vector, dimension(ndim)
       if (.not. allocated(detector%update_vector)) then
          allocate(detector%update_vector(ndims))
       end if       
       detector%update_vector = reshape(buff(ndims+4:2*ndims+3),(/ndims/))

       ! k, dimension(nstages:ndim)
       if (.not. allocated(detector%k)) then
          allocate(detector%k(nstages,ndims))
       end if  
       detector%k = reshape(buff(2*ndims+4:(nstages+2)*ndims+3),(/nstages,ndims/))

       ! If update_vector still exists, we're not done moving
       detector%search_complete=.false.
    else
       assert(size(buff)==ndims+3)

       detector%search_complete=.true.
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

end module detector_tools
